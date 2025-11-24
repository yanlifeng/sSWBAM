// main.c
// Sunway SAM 处理工具 - 主核代码
//
// 功能：
//   - 读取输入目录中的所有 SAM 文件
//   - 批量分配给 64 个 CPE（从核）并行处理
//   - 支持三种模式：排序、去重、排序+去重
//   - 将处理结果写回输出目录
//   - 统计各阶段耗时
//
// 编译：
//   make
//
// 运行：
//   ./sw_sam_process --all <input_dir> <output_dir>        # 排序 + 去重
//   ./sw_sam_process --sort <input_dir> <output_dir>       # 仅排序
//   ./sw_sam_process --markdup <input_dir> <output_dir>    # 仅去重
//
// 说明：
//   - input_dir: 输入 SAM 文件目录
//   - output_dir: 输出处理后的 SAM 文件目录
//   - --all: 先排序再去重（推荐用于完整流程）
//   - --sort: 仅按 RNAME（染色体）+ POS（位置）排序
//   - --markdup: 仅标记重复序列（需要输入已排序的文件）
//   - 单个文件大小限制：100MB（可调整 MAX_BUF_SIZE）

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <athread.h>
#include <pthread.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>
#include <unistd.h>
#include "../slave/sam_process_para.h"

extern void slave_sam_process_cpe(SamProcessPara paras[64]);

#define BATCH_SIZE      64
#define MAX_PATH_LEN    512
#define MAX_BASENAME    128
#define MAX_BUF_SIZE    (100UL * 1024UL * 1024UL)   // 100MB

static double now_ms()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

// 递归删除目录中的所有文件（不删除目录本身）
static int clear_directory(const char *path)
{
    DIR *dir = opendir(path);
    if (!dir) {
        return -1;
    }

    struct dirent *entry;
    char filepath[MAX_PATH_LEN];
    
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }
        
        snprintf(filepath, sizeof(filepath), "%s/%s", path, entry->d_name);
        
        struct stat st;
        if (stat(filepath, &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                // 递归删除子目录
                clear_directory(filepath);
                rmdir(filepath);
            } else {
                // 删除文件
                unlink(filepath);
            }
        }
    }
    
    closedir(dir);
    return 0;
}

// 准备输出目录：如果不存在则创建，如果存在则清空
static int prepare_output_directory(const char *path)
{
    struct stat st;
    
    if (stat(path, &st) == 0) {
        // 目录存在
        if (!S_ISDIR(st.st_mode)) {
            fprintf(stderr, "Error: %s exists but is not a directory\n", path);
            return -1;
        }
        
        // 清空目录
        printf("Output directory exists, clearing contents...\n");
        if (clear_directory(path) != 0) {
            fprintf(stderr, "Warning: Failed to clear directory %s\n", path);
        }
    } else {
        // 目录不存在，创建它
        printf("Creating output directory: %s\n", path);
        if (mkdir(path, 0755) != 0) {
            fprintf(stderr, "Error: Cannot create directory %s: %s\n", 
                    path, strerror(errno));
            return -1;
        }
    }
    
    return 0;
}

// 生成输出文件名
// 输入: input.sam, 模式: MODE_SORT_ONLY -> 输出: input.sorted.sam
// 输入: input.sam, 模式: MODE_MARKDUP_ONLY -> 输出: input.markdup.sam
// 输入: input.sam, 模式: MODE_ALL -> 输出: input.sorted.markdup.sam
static void generate_output_filename(const char *input_name, int mode, 
                                     char *output_name, size_t output_size)
{
    // 找到文件名中的 .sam 扩展名
    const char *ext = strstr(input_name, ".sam");
    size_t base_len;
    
    if (ext) {
        base_len = ext - input_name;
    } else {
        base_len = strlen(input_name);
    }
    
    // 复制基础文件名（不含扩展名）
    char base_name[MAX_BASENAME];
    if (base_len >= MAX_BASENAME) {
        base_len = MAX_BASENAME - 1;
    }
    strncpy(base_name, input_name, base_len);
    base_name[base_len] = '\0';
    
    // 根据模式生成输出文件名
    if (mode == MODE_SORT_ONLY) {
        snprintf(output_name, output_size, "%s.sorted.sam", base_name);
    } else if (mode == MODE_MARKDUP_ONLY) {
        snprintf(output_name, output_size, "%s.markdup.sam", base_name);
    } else if (mode == MODE_ALL) {
        snprintf(output_name, output_size, "%s.sorted.markdup.sam", base_name);
    } else {
        // 默认情况
        snprintf(output_name, output_size, "%s.processed.sam", base_name);
    }
}

// 处理一个批次（<=64 个文件）：
//  - 已经把每个文件读入 in_bufs[i]，大小 sizes[i]
//  - 对应的输出路径为 out_paths[i]
//  - 准备好 paras[i].in_buf/out_buf/size/mode，其他位置填 0。
//  - 调用从核处理（排序/去重/全流程），然后把 out_bufs 写回文件。
static void process_batch(int batch_count,
                          char in_paths[][MAX_PATH_LEN],
                          char out_paths[][MAX_PATH_LEN],
                          char *in_bufs[],
                          char *out_bufs[],
                          unsigned long sizes[],
                          int mode,
                          double *sort_ms_acc,
                          double *write_ms_acc)
{
    if (batch_count <= 0) return;

    SamProcessPara paras[64];
    unsigned long out_sizes[64];
    int i;

    for (i = 0; i < 64; ++i) {
        paras[i].in_buf  = 0;
        paras[i].out_buf = 0;
        paras[i].size    = 0;
        paras[i].out_buf_capacity = 0;
        paras[i].out_size = &(out_sizes[i]);
    }

    for (i = 0; i < batch_count; ++i) {
        paras[i].in_buf  = in_bufs[i];
        paras[i].out_buf = out_bufs[i];
        paras[i].size    = sizes[i];
        // 输出 buffer 容量是输入大小的 1.05 倍
        paras[i].out_buf_capacity = (unsigned long)((double)sizes[i] * 1.05);
        paras[i].mode    = mode;
    }

    double t0 = now_ms();

    // 并行处理（排序/去重/全流程）
    printf("  Spawning %d CPEs for parallel processing...\n", batch_count);
    __real_athread_spawn((void*)slave_sam_process_cpe, paras, 1);
    athread_join();

    double t1 = now_ms();
    *sort_ms_acc += (t1 - t0);
    printf("  CPE processing completed in %.3f ms\n", t1 - t0);
    
    // 立即释放输入 buffers 以节省内存（CPE 已经处理完毕，结果在 out_bufs 中）
    printf("  Freeing input buffers to save memory (batch_count=%d)...\n", batch_count);
    for (i = 0; i < batch_count; ++i) {
        printf("    Freeing in_buf[%d] at %p...\n", i, (void*)in_bufs[i]);
        if (in_bufs[i]) {
            free(in_bufs[i]);
            in_bufs[i] = NULL;
            printf("    in_buf[%d] freed\n", i);
        } else {
            printf("    in_buf[%d] is NULL, skipping\n", i);
        }
    }
    printf("  All input buffers freed\n");

    // 写回文件
    printf("  Writing results to output directory...\n");
    t0 = now_ms();
    int write_success = 0;
    int write_failed = 0;
    
    for (i = 0; i < batch_count; ++i) {
        printf("  [%d/%d] Processing file: %s\n", i+1, batch_count, out_paths[i]);
        
        unsigned long out_size = *(paras[i].out_size); // CPE 写回的实际输出长度
        printf("    Output size: %lu bytes (%.2f MB)\n", out_size, out_size / (1024.0 * 1024.0));
        
        // 检查 out_size 是否异常
        if (out_size > sizes[i] * 2) {
            fprintf(stderr, "    Error: Invalid out_size: %lu (input size: %lu)\n",
                    out_size, sizes[i]);
            fprintf(stderr, "    CPE processing may have failed. Skipping this file.\n");
            write_failed++;
            continue;
        }
        
        if (out_size == 0) {
            fprintf(stderr, "    Warning: CPE returned zero output size (input: %.2f MB)\n",
                    sizes[i] / (1024.0 * 1024.0));
            write_failed++;
            continue;
        }
        
        printf("    Opening output file...\n");
        FILE *fout = fopen(out_paths[i], "wb");
        if (!fout) {
            fprintf(stderr, "    Error: fopen failed: %s\n", strerror(errno));
            write_failed++;
            continue;
        }
        
        printf("    Writing %lu bytes...\n", out_size);
        size_t nwrite = fwrite(out_bufs[i], 1, (size_t)out_size, fout);
        printf("    fwrite returned: %zu bytes\n", nwrite);
        
        if (nwrite != out_size) {
            fprintf(stderr, "    Warning: fwrite incomplete: expect=%lu got=%zu\n",
                    out_size, nwrite);
            write_failed++;
        } else {
            printf("    Success: written %.2f MB\n", out_size / (1024.0 * 1024.0));
            write_success++;
        }
        
        fclose(fout);
        printf("    File closed\n");
        
        // 立即释放这个输出 buffer 以节省内存
        if (out_bufs[i]) {
            free(out_bufs[i]);
            out_bufs[i] = NULL;
            printf("    Output buffer freed\n");
        }
    }
    
    t1 = now_ms();
    *write_ms_acc += (t1 - t0);
    printf("  Write completed in %.3f ms (success: %d, failed: %d)\n", 
           t1 - t0, write_success, write_failed);
}

// 主函数：
//   argv[1] = 模式选项（--all, --sort, --markdup）
//   argv[2] = 输入目录
//   argv[3] = 输出目录
//
// 工作流程：
//   1. 解析命令行参数，确定处理模式
//   2. 遍历输入目录中的所有文件
//   3. 将文件读入内存缓冲区
//   4. 每凑满 64 个文件（或遍历完成），调用从核批量处理
//   5. 将处理结果写入输出目录
//   6. 输出统计信息（读取、处理、写入耗时）
int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <mode> <input_dir> <output_dir>\n"
                "Modes:\n"
                "  --all      : Sort + Mark duplicates (full pipeline)\n"
                "  --sort     : Sort only (by RNAME + POS)\n"
                "  --markdup  : Mark duplicates only (input must be sorted)\n"
                "\n"
                "Example:\n"
                "  %s --all /path/to/input /path/to/output\n"
                "  %s --sort /path/to/input /path/to/output\n"
                "  %s --markdup /path/to/sorted /path/to/marked\n",
                argv[0], argv[0], argv[0], argv[0]);
        return 1;
    }

    // 解析模式参数
    int mode = 0;
    const char *mode_str = argv[1];
    if (strcmp(mode_str, "--all") == 0) {
        mode = MODE_ALL;
    } else if (strcmp(mode_str, "--sort") == 0) {
        mode = MODE_SORT_ONLY;
    } else if (strcmp(mode_str, "--markdup") == 0) {
        mode = MODE_MARKDUP_ONLY;
    } else {
        fprintf(stderr, "Error: Invalid mode '%s'\n", mode_str);
        fprintf(stderr, "Valid modes: --all, --sort, --markdup\n");
        return 1;
    }

    const char *in_dir  = argv[2];
    const char *out_dir = argv[3];

    const char *mode_name = (mode == MODE_ALL) ? "Sort+Markdup" :
                            (mode == MODE_SORT_ONLY) ? "Sort" : "Markdup";
    printf("========================================\n");
    printf("Sunway SAM Processing Tool\n");
    printf("========================================\n");
    printf("Mode        : %s\n", mode_name);
    printf("Input dir   : %s\n", in_dir);
    printf("Output dir  : %s\n", out_dir);
    printf("========================================\n");

    // 准备输出目录
    if (prepare_output_directory(out_dir) != 0) {
        return 1;
    }

    DIR *dir = opendir(in_dir);
    if (!dir) {
        fprintf(stderr, "Error: Cannot open input directory %s: %s\n",
                in_dir, strerror(errno));
        return 1;
    }

    athread_init();
    printf("Athread initialized successfully\n\n");

    double read_ms  = 0.0;
    double sort_ms  = 0.0;
    double write_ms = 0.0;
    double total_start = now_ms();

    int total_files = 0;
    int total_batches = 0;

    char batch_basenames[BATCH_SIZE][MAX_BASENAME];
    char batch_inpaths[BATCH_SIZE][MAX_PATH_LEN];
    char batch_outpaths[BATCH_SIZE][MAX_PATH_LEN];

    char *in_bufs[BATCH_SIZE];
    char *out_bufs[BATCH_SIZE];
    unsigned long sizes[BATCH_SIZE];
    for (int i = 0; i < BATCH_SIZE; i++) sizes[i] = 0;

    int batch_count = 0;

    struct dirent *ent;
    while ((ent = readdir(dir)) != NULL) {
        if (ent->d_name[0] == '.') continue;  // skip . and ..

        // 路径/文件名
        snprintf(batch_basenames[batch_count], MAX_BASENAME,
                 "%s", ent->d_name);
        snprintf(batch_inpaths[batch_count], MAX_PATH_LEN,
                 "%s/%s", in_dir, ent->d_name);
        
        // 生成输出文件名
        char output_filename[MAX_BASENAME];
        generate_output_filename(ent->d_name, mode, output_filename, sizeof(output_filename));
        snprintf(batch_outpaths[batch_count], MAX_PATH_LEN,
                 "%s/%s", out_dir, output_filename);

        // 读文件到 in_buf
        struct stat st;
        if (stat(batch_inpaths[batch_count], &st) != 0) {
            fprintf(stderr, "stat failed: %s (%s)\n",
                    batch_inpaths[batch_count], strerror(errno));
            continue;
        }

        unsigned long fsize = (unsigned long)st.st_size;
        if (fsize > MAX_BUF_SIZE) {
            fprintf(stderr,
                    "Warning: File too large (> %lu MB): %s, size=%lu bytes - skipped\n",
                    (unsigned long)(MAX_BUF_SIZE / 1024 / 1024),
                    batch_inpaths[batch_count], fsize);
            continue;
        }
        
        printf("Reading file [%d]: %s (%.2f MB)\n", 
               total_files + 1, ent->d_name, fsize / (1024.0 * 1024.0));

        double t0 = now_ms();

        FILE *fin = fopen(batch_inpaths[batch_count], "rb");
        if (!fin) {
            fprintf(stderr, "fopen input failed: %s (%s)\n",
                    batch_inpaths[batch_count], strerror(errno));
            continue;
        }

        char *ibuf = 0;
        char *obuf = 0;
        if (fsize > 0) {
            // 输入和输出 buffer 都分配为 1.05 倍，因为：
            // 1. markdup 时 FLAG 字段可能变大（例如 "12" -> "1036"）
            // 2. MODE_ALL 模式下需要在 in_buf 和 out_buf 之间交换数据
            unsigned long buf_size = (unsigned long)((double)fsize * 1.05);
            ibuf = (char*)malloc((size_t)buf_size);
            obuf = (char*)malloc((size_t)buf_size);
            if (!ibuf || !obuf) {
                fprintf(stderr, "malloc buf failed for %s (size=%lu)\n",
                        batch_inpaths[batch_count], buf_size);
                if (ibuf) free(ibuf);
                if (obuf) free(obuf);
                fclose(fin);
                continue;
            }
            size_t nread = fread(ibuf, 1, (size_t)fsize, fin);
            if (nread != fsize) {
                fprintf(stderr,
                        "fread incomplete for %s: expect=%lu got=%zu\n",
                        batch_inpaths[batch_count], fsize, nread);
                free(ibuf);
                free(obuf);
                fclose(fin);
                continue;
            }
        }
        fclose(fin);

        double t1 = now_ms();
        read_ms += (t1 - t0);

        in_bufs[batch_count]  = ibuf;
        out_bufs[batch_count] = obuf;
        sizes[batch_count]    = fsize;

        batch_count++;
        total_files++;

        // 凑满 64 个就发一次从核
        if (batch_count == BATCH_SIZE) {
            total_batches++;
            printf("\n--- Processing Batch %d (%d files) ---\n", total_batches, batch_count);
            process_batch(batch_count,
                          batch_inpaths, batch_outpaths,
                          in_bufs, out_bufs, sizes,
                          mode,
                          &sort_ms, &write_ms);
            printf("Batch %d completed\n\n", total_batches);

            // 批处理后清理（buffers 已在 process_batch 中释放）
            int i;
            for (i = 0; i < batch_count; ++i) {
                in_bufs[i] = NULL;
                out_bufs[i] = NULL;
                sizes[i] = 0;
            }
            batch_count = 0;
        }
    }

    // 处理最后一批 (<64)
    if (batch_count > 0) {
        total_batches++;
        printf("\n--- Processing Final Batch %d (%d files) ---\n", total_batches, batch_count);
        process_batch(batch_count,
                      batch_inpaths, batch_outpaths,
                      in_bufs, out_bufs, sizes,
                      mode,
                      &sort_ms, &write_ms);
        printf("Final batch completed\n\n");
        int i;
        for (i = 0; i < batch_count; ++i) {
            // Buffers 已在 process_batch 中释放
            in_bufs[i] = NULL;
            out_bufs[i] = NULL;
            sizes[i] = 0;
        }
    }

    closedir(dir);

    double total_end = now_ms();
    double total_ms  = total_end - total_start;

    printf("\n========================================\n");
    printf("Processing Summary\n");
    printf("========================================\n");
    printf("Mode              : %s\n", mode_name);
    printf("Total batches     : %d\n", total_batches);
    printf("Files processed   : %d\n", total_files);
    printf("----------------------------------------\n");
    printf("Read time         : %.3f ms (%.2f%%)\n", read_ms, (read_ms / total_ms) * 100);
    printf("Process(CPE) time : %.3f ms (%.2f%%)\n", sort_ms, (sort_ms / total_ms) * 100);
    printf("Write time        : %.3f ms (%.2f%%)\n", write_ms, (write_ms / total_ms) * 100);
    printf("----------------------------------------\n");
    printf("Total time        : %.3f ms (%.2f s)\n", total_ms, total_ms / 1000.0);
    printf("========================================\n");

    return 0;
}

