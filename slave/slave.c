// slave.c
// Sunway 从核代码：负责对自己那块 SAM buffer 排序（in_buf -> out_buf）
//
// 编译命令：
//   sw5cc -slave slave.c -o slave.o

#include <slave.h>
#include <stdlib.h>
#include <string.h>
#include "sam_sort_para.h"

typedef struct {
    unsigned long start;      // 在 in_buf 中的起始偏移
    unsigned long len;        // 这一行长度（包含 '\n' 如果有）
    const char   *rname;      // RNAME 指针（指向 in_buf 内部）
    int           rname_len;  // RNAME 长度
    long          pos;        // POS
    int           valid;      // 是否成功解析出 RNAME+POS
} LineInfo;

// 解析一行（不含 '\n'）的 RNAME + POS
static void parse_rname_pos(char *line_start, unsigned long text_len, LineInfo *info)
{
    info->rname     = 0;
    info->rname_len = 0;
    info->pos       = -1;
    info->valid     = 0;

    if (text_len == 0) return;

    char *p   = line_start;
    char *end = line_start + text_len;

    // SAM: 0:QNAME, 1:FLAG, 2:RNAME, 3:POS, ...
    int   field = 0;
    char *field_start = p;

    char *rname_start = 0;
    char *rname_end   = 0;
    char *pos_start   = 0;
    char *pos_end     = 0;

    while (p <= end) {
        if (p == end || *p == '\t') {
            if (field == 2) {          // RNAME
                rname_start = field_start;
                rname_end   = p;
            } else if (field == 3) {   // POS
                pos_start = field_start;
                pos_end   = p;
                break;                 // POS 后面不关心
            }
            field++;
            field_start = p + 1;
        }
        if (p == end) break;
        ++p;
    }

    if (!rname_start || !pos_start) {
        info->valid = 0;
        return;
    }

    info->rname     = rname_start;
    info->rname_len = (int)(rname_end - rname_start);

    // 解析 POS
    long value = 0;
    char *q    = pos_start;
    int   neg  = 0;
    if (q < pos_end && *q == '-') {
        neg = 1;
        ++q;
    }
    if (q == pos_end) {
        info->valid = 0;
        return;
    }
    while (q < pos_end && *q >= '0' && *q <= '9') {
        value = value * 10 + (*q - '0');
        ++q;
    }
    if (q != pos_end) {
        info->valid = 0;
        return;
    }
    if (neg) value = -value;
    info->pos   = value;
    info->valid = 1;
}

// 比较 RNAME（字典序）
static int cmp_rname(const LineInfo *a, const LineInfo *b)
{
    int len = (a->rname_len < b->rname_len) ? a->rname_len : b->rname_len;
    int r = 0;
    if (len > 0) r = memcmp(a->rname, b->rname, (unsigned long)len);
    if (r != 0) return r;
    if (a->rname_len < b->rname_len) return -1;
    if (a->rname_len > b->rname_len) return 1;
    return 0;
}

// 比较两个行：按 RNAME, 再按 POS
static int cmp_line(const LineInfo *a, const LineInfo *b)
{
    if (!a->valid && !b->valid) {
        if (a->start < b->start) return -1;
        if (a->start > b->start) return 1;
        return 0;
    } else if (!a->valid) {
        return -1;
    } else if (!b->valid) {
        return 1;
    }

    int cr = cmp_rname(a, b);
    if (cr < 0) return -1;
    if (cr > 0) return 1;
    if (a->pos < b->pos) return -1;
    if (a->pos > b->pos) return 1;

    // 完全相同，按原始顺序（start 小的在前）
    if (a->start < b->start) return -1;
    if (a->start > b->start) return 1;
    return 0;
}

// 快速排序（原地）
// arr: LineInfo 数组
// left, right: [left, right]
static void quicksort_lineinfo(LineInfo *arr, int left, int right)
{
    int i = left;
    int j = right;
    LineInfo pivot = arr[(left + right) / 2];

    while (i <= j) {
        while (cmp_line(&arr[i], &pivot) < 0) ++i;
        while (cmp_line(&arr[j], &pivot) > 0) --j;
        if (i <= j) {
            LineInfo tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            ++i;
            --j;
        }
    }
    if (left < j)  quicksort_lineinfo(arr, left, j);
    if (i < right) quicksort_lineinfo(arr, i, right);
}

// 解析一个 buffer 里的所有 SAM 行
static int parse_sam_lines(char *buf, unsigned long size, LineInfo **lines_out)
{
    *lines_out = 0;
    if (size == 0) return 0;

    // 1) 统计行数
    int  n_lines = 0;
    int  last_is_nl = 0;
    unsigned long i;
    for (i = 0; i < size; ++i) {
        if (buf[i] == '\n') {
            n_lines++;
            last_is_nl = 1;
        } else {
            last_is_nl = 0;
        }
    }
    if (!last_is_nl) n_lines++;  // 最后一行可能没有 '\n'

    if (n_lines <= 0) return 0;

    LineInfo *lines = (LineInfo*)malloc(sizeof(LineInfo) * (unsigned long)n_lines);
    if (!lines) {
        return 0;
    }

    // 2) 填写每行信息并解析 RNAME+POS
    int idx = 0;
    i = 0;
    while (i < size && idx < n_lines) {
        unsigned long start = i;
        while (i < size && buf[i] != '\n') ++i;
        unsigned long line_end = i;              // 不含 '\n'
        int has_nl = (i < size && buf[i] == '\n');
        if (has_nl) ++i;

        unsigned long len = has_nl ? (i - start) : (line_end - start);

        lines[idx].start = start;
        lines[idx].len   = len;
        lines[idx].rname = 0;
        lines[idx].rname_len = 0;
        lines[idx].pos   = -1;
        lines[idx].valid = 0;

        unsigned long text_len = line_end - start; // 不含 '\n'
        if (text_len > 0) {
            parse_rname_pos(buf + start, text_len, &lines[idx]);
        }
        idx++;
    }

    *lines_out = lines;
    return idx;
}

// 从核入口：每个 CPE 负责 paras[_PEN] 这一份
void slave_sam_sort_cpe(SamSortPara paras[64])
{
    int id = _PEN;
    SamSortPara *para = &paras[id];

    char         *in_buf  = para->in_buf;
    char         *out_buf = para->out_buf;
    unsigned long size    = para->size;

    if (!in_buf || !out_buf || size == 0) {
        return;
    }

    LineInfo *lines = 0;
    int n_lines = parse_sam_lines(in_buf, size, &lines);

    if (n_lines > 1 && lines) {
        quicksort_lineinfo(lines, 0, n_lines - 1);
    }

    // 按排序结果写入 out_buf
    unsigned long out_pos = 0;
    int i;
    for (i = 0; i < n_lines; ++i) {
        if (lines[i].len == 0) continue;
        // 简单边界保护：不应该超过原 size（正常情况 len 总和 == size）
        if (out_pos + lines[i].len > size) break;
        memcpy(out_buf + out_pos,
               in_buf  + lines[i].start,
               (unsigned long)lines[i].len);
        out_pos += lines[i].len;
    }

    *(para->out_size) = out_pos;
    if (lines) free(lines);
}

