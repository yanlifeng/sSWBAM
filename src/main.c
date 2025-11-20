// main.c
// Sunway SAM 排序工具 - 主核代码
//
// 功能：
//   - 读取输入目录中的所有 SAM 文件
//   - 批量分配给 64 个 CPE（从核）并行排序
//   - 将排序结果写回输出目录
//   - 统计各阶段耗时
//
// 编译：
//   make
//
// 运行：
//   ./sw_sam_sort <input_dir> <output_dir>
//
// 说明：
//   - input_dir: 输入 SAM 文件目录（通常是 out_regions_sam）
//   - output_dir: 输出排序后的 SAM 文件目录（通常是 out_regions_sorted_sam）
//   - 输出文件名格式：<原文件名>.sorted.sw.sam
//   - 每个 SAM 文件按 RNAME（染色体）+ POS（位置）排序
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
#include "../slave/sam_sort_para.h"

extern void slave_sam_sort_cpe(SamSortPara paras[64]);

static double now_ms()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

#define BATCH_SIZE      64
#define MAX_PATH_LEN    512
#define MAX_BASENAME    128
#define MAX_BUF_SIZE    (100UL * 1024UL * 1024UL)   // 100MB

// 处理一个批次（<=64 个文件）：
//  - 已经把每个文件读入 in_bufs[i]，大小 sizes[i]
//  - 对应的输出路径为 out_paths[i]
//  - 准备好 paras[i].in_buf/out_buf/size，其他位置填 0。
//  - 调用从核排序，然后把 out_bufs 写回文件。
static void process_batch(int batch_count,
                          char in_paths[][MAX_PATH_LEN],
                          char out_paths[][MAX_PATH_LEN],
                          char *in_bufs[],
                          char *out_bufs[],
                          unsigned long sizes[],
                          double *sort_ms_acc,
                          double *write_ms_acc)
{
    if (batch_count <= 0) return;

    SamSortPara paras[64];
    unsigned long out_sizes[64];
    int i;

    for (i = 0; i < 64; ++i) {
        paras[i].in_buf  = 0;
        paras[i].out_buf = 0;
        paras[i].size    = 0;
        paras[i].out_size = &(out_sizes[i]);
    }

    for (i = 0; i < batch_count; ++i) {
        paras[i].in_buf  = in_bufs[i];
        paras[i].out_buf = out_bufs[i];
        paras[i].size    = sizes[i];
    }

    double t0 = now_ms();

    // 并行排序
    __real_athread_spawn((void*)slave_sam_sort_cpe, paras, 1);
    athread_join();

    double t1 = now_ms();
    *sort_ms_acc += (t1 - t0);

    // 写回文件
    t0 = now_ms();
    for (i = 0; i < batch_count; ++i) {
        unsigned long out_size = *(paras[i].out_size); // CPE 写回的实际输出长度
        FILE *fout = fopen(out_paths[i], "wb");
        if (!fout) {
            fprintf(stderr, "fopen output failed: %s (%s)\n",
                    out_paths[i], strerror(errno));
            continue;
        }
        if (out_size > 0) {
            size_t nwrite = fwrite(out_bufs[i], 1, (size_t)out_size, fout);
            if (nwrite != out_size) {
                fprintf(stderr,
                        "fwrite incomplete for %s: expect=%lu got=%zu\n",
                        out_paths[i], out_size, nwrite);
            }
        }
        fclose(fout);
    }
    t1 = now_ms();
    *write_ms_acc += (t1 - t0);
}

// 主函数：
//   argv[1] = 输入目录（例如 out_regions_sam）
//   argv[2] = 输出目录（例如 out_regions_sorted_sam）
//
// 工作流程：
//   1. 遍历输入目录中的所有文件
//   2. 将文件读入内存缓冲区
//   3. 每凑满 64 个文件（或遍历完成），调用从核批量排序
//   4. 将排序结果写入输出目录
//   5. 输出统计信息（读取、排序、写入耗时）
int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s <input_dir> <output_dir>\n"
                "Example:\n"
                "  %s /path/to/out /path/to/out_sort\n",
                argv[0], argv[0]);
        return 1;
    }

    const char *in_dir  = argv[1];
    const char *out_dir = argv[2];

    DIR *dir = opendir(in_dir);
    if (!dir) {
        fprintf(stderr, "opendir failed: %s (%s)\n",
                in_dir, strerror(errno));
        return 1;
    }

    athread_init();

    double read_ms  = 0.0;
    double sort_ms  = 0.0;
    double write_ms = 0.0;
    double total_start = now_ms();

    int total_files = 0;

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
        snprintf(batch_outpaths[batch_count], MAX_PATH_LEN,
                 "%s/%s.sorted.sw.sam", out_dir, ent->d_name);

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
                    "File too large (> %lu MB): %s, size=%lu\n",
                    (unsigned long)(MAX_BUF_SIZE / 1024 / 1024),
                    batch_inpaths[batch_count], fsize);
            continue;
        }

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
            ibuf = (char*)malloc((size_t)fsize);
            obuf = (char*)malloc((size_t)fsize);
            if (!ibuf || !obuf) {
                fprintf(stderr, "malloc buf failed for %s size=%lu\n",
                        batch_inpaths[batch_count], fsize);
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
            process_batch(batch_count,
                          batch_inpaths, batch_outpaths,
                          in_bufs, out_bufs, sizes,
                          &sort_ms, &write_ms);

            // 批处理后释放缓冲
            int i;
            for (i = 0; i < batch_count; ++i) {
                if (in_bufs[i])  free(in_bufs[i]);
                if (out_bufs[i]) free(out_bufs[i]);
                sizes[i] = 0;
            }
            batch_count = 0;
        }
    }

    // 处理最后一批 (<64)
    if (batch_count > 0) {
        process_batch(batch_count,
                      batch_inpaths, batch_outpaths,
                      in_bufs, out_bufs, sizes,
                      &sort_ms, &write_ms);
        int i;
        for (i = 0; i < batch_count; ++i) {
            if (in_bufs[i])  free(in_bufs[i]);
            if (out_bufs[i]) free(out_bufs[i]);
            sizes[i] = 0;
        }
    }

    closedir(dir);

    double total_end = now_ms();
    double total_ms  = total_end - total_start;

    fprintf(stderr,
            "\n==== Summary ====\n"
            "Files processed : %d\n"
            "Read time       : %.3f ms\n"
            "Sort(CPE) time  : %.3f ms\n"
            "Write time      : %.3f ms\n"
            "Total time      : %.3f ms\n",
            total_files, read_ms, sort_ms, write_ms, total_ms);

    return 0;
}

