#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>

#include <athread.h>

typedef struct {
    char *in_buf;
    char *out_buf;
    unsigned long size;
    unsigned long *out_size;
} SamMarkPara;

extern void slave_sam_markdup_cpe(void *arg);

static double now_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}

static int is_regular_file(const char *dir, const char *name) {
    char path[4096];
    struct stat st;

    if (strcmp(name, ".") == 0 || strcmp(name, "..") == 0)
        return 0;

    snprintf(path, sizeof(path), "%s/%s", dir, name);
    if (stat(path, &st) != 0)
        return 0;

    return S_ISREG(st.st_mode);
}

static char **build_file_list(const char *input_dir, int *count) {
    DIR *dp;
    struct dirent *entry;
    char **list = NULL;
    int capacity = 0;
    int n = 0;

    dp = opendir(input_dir);
    if (!dp) {
        fprintf(stderr, "Cannot open input directory %s: %s\n",
                input_dir, strerror(errno));
        return NULL;
    }

    while ((entry = readdir(dp)) != NULL) {
        if (!is_regular_file(input_dir, entry->d_name))
            continue;

        if (n >= capacity) {
            capacity = capacity ? capacity * 2 : 64;
            char **tmp = (char **)realloc(list, capacity * sizeof(char *));
            if (!tmp) {
                fprintf(stderr, "Out of memory\n");
                closedir(dp);
                for (int i = 0; i < n; ++i)
                    free(list[i]);
                free(list);
                return NULL;
            }
            list = tmp;
        }

        size_t len = strlen(input_dir) + 1 + strlen(entry->d_name) + 1;
        list[n] = (char *)malloc(len);
        if (!list[n]) {
            fprintf(stderr, "Out of memory\n");
            closedir(dp);
            for (int i = 0; i < n; ++i)
                free(list[i]);
            free(list);
            return NULL;
        }
        snprintf(list[n], len, "%s/%s", input_dir, entry->d_name);
        ++n;
    }
    closedir(dp);

    *count = n;
    return list;
}

static const char *basename_const(const char *path) {
    const char *p = strrchr(path, '/');
    return p ? p + 1 : path;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_dir> <output_dir>\n", argv[0]);
        return 1;
    }

    const char *input_dir = argv[1];
    const char *output_dir = argv[2];

    struct stat st;
    if (stat(output_dir, &st) != 0) {
        if (mkdir(output_dir, 0755) != 0) {
            fprintf(stderr, "Cannot create output directory %s: %s\n",
                    output_dir, strerror(errno));
            return 1;
        }
    } else if (!S_ISDIR(st.st_mode)) {
        fprintf(stderr, "%s exists and is not a directory\n", output_dir);
        return 1;
    }

    int file_count = 0;
    char **file_list = build_file_list(input_dir, &file_count);
    if (!file_list || file_count == 0) {
        fprintf(stderr, "No input files found in %s\n", input_dir);
        free(file_list);
        return 1;
    }

    printf("Found %d files in %s\n", file_count, input_dir);

    athread_init();

    unsigned long total_memory_used = 0;  /* Track MPE memory usage in bytes */
    int idx = 0;
    while (idx < file_count) {
        int batch_count = file_count - idx;
        if (batch_count > 64)
            batch_count = 64;

        char *in_bufs[64];
        char *out_bufs[64];
        unsigned long sizes[64];
        unsigned long out_sizes[64];
        SamMarkPara paras[64];

        memset(in_bufs, 0, sizeof(in_bufs));
        memset(out_bufs, 0, sizeof(out_bufs));
        memset(sizes, 0, sizeof(sizes));
        memset(out_sizes, 0, sizeof(out_sizes));

        for (int i = 0; i < 64; ++i) {
            paras[i].in_buf = NULL;
            paras[i].out_buf = NULL;
            paras[i].size = 0;
            paras[i].out_size = &out_sizes[i];
        }

        /* Open and read files in current batch into memory */
        for (int i = 0; i < batch_count; ++i) {
            const char *path = file_list[idx + i];
            FILE *fp = fopen(path, "rb");
            if (!fp) {
                fprintf(stderr, "Cannot open %s: %s\n", path, strerror(errno));
                continue;
            }

            if (fseeko(fp, 0, SEEK_END) != 0) {
                fprintf(stderr, "fseeko failed on %s\n", path);
                fclose(fp);
                continue;
            }
            off_t fsize = ftello(fp);
            if (fsize < 0) {
                fprintf(stderr, "ftello failed on %s\n", path);
                fclose(fp);
                continue;
            }
            if (fseeko(fp, 0, SEEK_SET) != 0) {
                fprintf(stderr, "fseeko failed on %s\n", path);
                fclose(fp);
                continue;
            }

            in_bufs[i] = (char *)malloc((size_t)fsize);
            out_bufs[i] = (char *)malloc((size_t)fsize);
            if (!in_bufs[i] || !out_bufs[i]) {
                fprintf(stderr, "Out of memory allocating buffers for %s\n", path);
                fclose(fp);
                free(in_bufs[i]);
                free(out_bufs[i]);
                in_bufs[i] = out_bufs[i] = NULL;
                continue;
            }
            
            /* Track memory allocation */
            total_memory_used += (size_t)fsize * 2;  /* in_buf + out_buf */
            printf("Allocated %.2f MB for file %s, total MPE memory: %.2f MB\n",
                   (size_t)fsize * 2 / (1024.0 * 1024.0),
                   path,
                   total_memory_used / (1024.0 * 1024.0));
            
            size_t nread = fread(in_bufs[i], 1, (size_t)fsize, fp);
            fclose(fp);
            if (nread != (size_t)fsize) {
                fprintf(stderr, "Short read on %s\n", path);
                free(in_bufs[i]);
                free(out_bufs[i]);
                in_bufs[i] = out_bufs[i] = NULL;
                continue;
            }

            sizes[i] = (unsigned long)fsize;
        }

        /* Set valid parameters for current batch */
        for (int i = 0; i < batch_count; ++i) {
            paras[i].in_buf = in_bufs[i];
            paras[i].out_buf = out_bufs[i];
            paras[i].size = sizes[i];
        }

        double t0 = now_ms();

        /* Launch 64 CPEs to process these buffers */
        __real_athread_spawn((void *)slave_sam_markdup_cpe, paras, 1);
        athread_join();

        double t1 = now_ms();
        printf("Processed batch %d (count=%d) in %.3f ms, MPE memory: %.2f MB\n",
               idx / 64, batch_count, t1 - t0, total_memory_used / (1024.0 * 1024.0));

        /* Write results to new files in output directory */
        for (int i = 0; i < batch_count; ++i) {
            if (!out_bufs[i])
                continue;

            const char *in_name = basename_const(file_list[idx + i]);
            char out_path[4096];
            snprintf(out_path, sizeof(out_path), "%s/%s", output_dir, in_name);

            FILE *fp = fopen(out_path, "wb");
            if (!fp) {
                fprintf(stderr, "Cannot open output %s: %s\n",
                        out_path, strerror(errno));
                continue;
            }

            unsigned long out_sz = out_sizes[i];
            if (out_sz == 0)
                out_sz = sizes[i];

            size_t nwritten = fwrite(out_bufs[i], 1, (size_t)out_sz, fp);
            if (nwritten != (size_t)out_sz) {
                fprintf(stderr, "Short write on %s\n", out_path);
            }
            fclose(fp);
        }

        /* Free buffers for current batch */
        for (int i = 0; i < batch_count; ++i) {
            if (in_bufs[i]) {
                total_memory_used -= sizes[i];
                free(in_bufs[i]);
            }
            if (out_bufs[i]) {
                total_memory_used -= sizes[i];
                free(out_bufs[i]);
            }
        }

        idx += batch_count;
    }

    for (int i = 0; i < file_count; ++i)
        free(file_list[i]);
    free(file_list);

    athread_halt();

    return 0;
}


