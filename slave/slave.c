// slave.c
// Sunway 从核代码：负责对 SAM buffer 进行排序和/或去重
//
// 功能：
//   - MODE_SORT_ONLY: 仅排序（按 RNAME + POS）
//   - MODE_MARKDUP_ONLY: 仅标记重复（需要输入已排序）
//   - MODE_ALL: 先排序再标记重复（完整流程）
//
// 编译命令：
//   sw5cc -slave slave.c -o slave.o

#include <slave.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "sam_process_para.h"

// ==================== SAM 排序相关结构和函数 ====================

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

// ==================== SAM 去重相关结构和函数 ====================

/* SAM flags */
#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024
#define BAM_FSUPPLEMENTARY 2048

/* Reference name to ID mapping */
#define MAX_REFS 256

typedef struct {
    char names[MAX_REFS][256];
    int count;
} ref_map_t;

/* Compact SAM record for duplicate detection */
typedef struct {
    uint32_t line_offset;   /* start of this SAM line in input buffer */
    uint32_t flag_offset;   /* offset of FLAG field */
    int32_t  pos;           /* position */
    int32_t  mate_pos;      /* mate position */
    uint16_t line_len;      /* line length */
    uint16_t flag_len;      /* FLAG field length */
    int16_t  tid;           /* reference ID */
    int16_t  mate_tid;      /* mate reference ID */
    uint16_t flag;          /* original FLAG */
    uint16_t score;         /* quality score (truncated) */
    uint8_t  orientation;   /* orientation for duplicate detection */
    uint8_t  is_duplicate;  /* 0/1: marked as duplicate */
} sam_record_t;

typedef struct {
    sam_record_t *records;
    int count;
    int capacity;
    ref_map_t ref_map;
} record_list_t;

/* Get or create reference ID from name */
static int get_ref_id(ref_map_t *map, const char *rname, int len) {
    /* Special cases */
    if (len == 1 && rname[0] == '*')
        return -1;  /* Unmapped */
    
    /* Search existing */
    for (int i = 0; i < map->count; i++) {
        if (strncmp(map->names[i], rname, len) == 0 && map->names[i][len] == '\0')
            return i;
    }
    
    /* Add new */
    if (map->count >= MAX_REFS)
        return -1;
    
    if (len >= 256) len = 255;
    memcpy(map->names[map->count], rname, len);
    map->names[map->count][len] = '\0';
    return map->count++;
}

/* Initialize record list */
static record_list_t *record_list_init(void) {
    record_list_t *list = calloc(1, sizeof(record_list_t));
    if (!list) {
        return NULL;
    }
    
    list->capacity = 1024;
    list->records = malloc(list->capacity * sizeof(sam_record_t));
    if (!list->records) {
        free(list);
        return NULL;
    }
    list->count = 0;
    list->ref_map.count = 0;

    return list;
}

/* Free record list */
static void record_list_free(record_list_t *list) {
    if (!list) return;
    free(list->records);
    free(list);
}

/* Calculate quality score from quality string */
static int64_t calc_score(const char *qual, int len) {
    int64_t score = 0;
    for (int i = 0; i < len; i++) {
        int q = qual[i] - 33;  /* Phred+33 */
        if (q > 15) q = 15;    /* Cap at 15 as per markdup spec */
        if (q > 0) score += q;
    }
    return score;
}

/* Parse one line of SAM - extract only fields needed for duplicate detection */
static int parse_sam_line_markdup(const char *line, unsigned long line_offset, int len, 
                   sam_record_t *rec, ref_map_t *ref_map) {
    const char *p = line;
    const char *end = line + len;
    int field = 0;
    const char *field_start = p;
    
    memset(rec, 0, sizeof(sam_record_t));
    rec->line_offset = (uint32_t)line_offset;
    rec->line_len = (uint16_t)len;
    
    while (p < end && field < 11) {
        if (*p == '\t' || p == end - 1) {
            int field_len = (p == end - 1 && *p != '\t') ? (end - field_start) : (p - field_start);
            
            switch (field) {
                case 1: /* FLAG */
                    rec->flag_offset = (uint32_t)(line_offset + (field_start - line));
                    rec->flag_len = (uint16_t)field_len;
                    rec->flag = (uint16_t)atoi(field_start);
                    break;
                    
                case 2: /* RNAME */
                    rec->tid = (int16_t)get_ref_id(ref_map, field_start, field_len);
                    break;
                    
                case 3: /* POS */
                    rec->pos = (int32_t)atoi(field_start);
                    break;
                    
                case 6: /* RNEXT */
                    if (field_len == 1 && field_start[0] == '=') {
                        rec->mate_tid = rec->tid;
                    } else {
                        rec->mate_tid = (int16_t)get_ref_id(ref_map, field_start, field_len);
                    }
                    break;
                    
                case 7: /* PNEXT */
                    rec->mate_pos = (int32_t)atoi(field_start);
                    break;
                    
                case 10: /* QUAL */
                    {
                        int64_t s = calc_score(field_start, field_len);
                        if (s < 0) s = 0;
                        if (s > 65535) s = 65535;
                        rec->score = (uint16_t)s;
                    }
                    break;
            }
            
            field++;
            field_start = p + 1;
        }
        p++;
    }
    
    /* Calculate orientation */
    if (rec->flag & BAM_FPAIRED) {
        rec->orientation = (rec->flag & BAM_FREVERSE) ? 1 : 0;
        rec->orientation |= ((rec->flag & BAM_FMREVERSE) ? 1 : 0) << 1;
    }
    
    return (field >= 11) ? 0 : -1;
}

/* Comparison function for sorting records by position */
static int compare_records(const void *a, const void *b) {
    const sam_record_t *ra = (const sam_record_t *)a;
    const sam_record_t *rb = (const sam_record_t *)b;
    
    /* Sort by tid, pos, mate_tid, mate_pos, orientation */
    if (ra->tid != rb->tid) return ra->tid - rb->tid;
    if (ra->pos != rb->pos) return (ra->pos < rb->pos) ? -1 : 1;
    if (ra->mate_tid != rb->mate_tid) return ra->mate_tid - rb->mate_tid;
    if (ra->mate_pos != rb->mate_pos) return (ra->mate_pos < rb->mate_pos) ? -1 : 1;
    if (ra->orientation != rb->orientation) return ra->orientation - rb->orientation;
    return 0;
}

/* Mark duplicates by sorting instead of hashing - saves memory! */
static void mark_duplicates_sorted(record_list_t *list) {
    if (list->count <= 1) return;
    
    /* Sort records by position */
    qsort(list->records, list->count, sizeof(sam_record_t), compare_records);
    
    /* Scan sorted array to find duplicates */
    int i = 0;
    while (i < list->count) {
        /* Skip unmapped, secondary, supplementary */
        if (list->records[i].flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            i++;
            continue;
        }
        
        /* Find all records with same key */
        int best = i;
        int j = i + 1;
        
        while (j < list->count && compare_records(&list->records[i], &list->records[j]) == 0) {
            /* Same position - keep the one with highest score */
            if (list->records[j].score > list->records[best].score) {
                list->records[best].is_duplicate = 1;
                best = j;
            } else {
                list->records[j].is_duplicate = 1;
            }
            j++;
        }
        
        i = j;
    }
}

/* Write SAM record with modified FLAG */
static int write_sam_record(const char *in_buf, char *out_buf, 
                     unsigned long *out_pos, unsigned long out_capacity, 
                     sam_record_t *rec) {
    unsigned long needed;
    char flag_str[32];
    int flag_str_len;
    
    /* Update FLAG if duplicate */
    int new_flag = rec->flag;
    if (rec->is_duplicate) {
        new_flag |= BAM_FDUP;
    }
    
    flag_str_len = snprintf(flag_str, sizeof(flag_str), "%d", new_flag);
    
    /* Calculate space needed */
    needed = rec->line_len - rec->flag_len + flag_str_len + 1;  /* +1 for newline */
    
    if (*out_pos + needed > out_capacity)
        return -1;

    /* Copy part before FLAG */
    unsigned long before_flag = rec->flag_offset - rec->line_offset;
    memcpy(out_buf + *out_pos, in_buf + rec->line_offset, before_flag);
    *out_pos += before_flag;
    
    /* Write new FLAG */
    memcpy(out_buf + *out_pos, flag_str, flag_str_len);
    *out_pos += flag_str_len;
    
    /* Copy part after FLAG */
    unsigned long after_flag_offset = rec->flag_offset + rec->flag_len;
    unsigned long after_flag_len = rec->line_offset + rec->line_len - after_flag_offset;
    memcpy(out_buf + *out_pos, in_buf + after_flag_offset, after_flag_len);
    *out_pos += after_flag_len;
    
    /* Add newline */
    out_buf[(*out_pos)++] = '\n';

    return 0;
}

/* Core markdup function */
static int markdup_core(const char *in_buf, char *out_buf, 
                 unsigned long size, unsigned long out_buf_capacity, 
                 unsigned long *out_size) {
    record_list_t *list = NULL;
    int ret = -1;
    
    // 初始化 out_size 为 0，防止使用未初始化的值
    if (out_size) {
        *out_size = 0;
    }
    
    if (!in_buf || !out_buf || size == 0) {
        return -1;
    }

    /* Initialize record list */
    list = record_list_init();
    if (!list) {
        return -1;
    }

    /* First pass: read all records */
    unsigned long pos = 0;
    while (pos < size) {
        // Skip headers and empty lines
        while (pos < size && (in_buf[pos] == '\n' || in_buf[pos] == '\r' || in_buf[pos] == '@')) {
            if (in_buf[pos] == '@') {
                // Skip entire header line
                while (pos < size && in_buf[pos] != '\n') pos++;
                if (pos < size) pos++;
            } else {
                pos++;
            }
        }
        
        if (pos >= size) break;
        
        unsigned long line_start = pos;
        while (pos < size && in_buf[pos] != '\n' && in_buf[pos] != '\r') pos++;
        int line_len = pos - line_start;
        if (pos < size && in_buf[pos] == '\n') pos++;
        
        if (line_len == 0) continue;
        
        if (list->count >= list->capacity) {
            int new_cap = list->capacity * 2;
            sam_record_t *new_recs = realloc(list->records, new_cap * sizeof(sam_record_t));
            if (!new_recs) {
                goto cleanup;
            }
            list->records = new_recs;
            list->capacity = new_cap;
        }
        
        sam_record_t *rec = &list->records[list->count];
        if (parse_sam_line_markdup(in_buf + line_start, line_start, line_len, rec, &list->ref_map) < 0) {
            continue;
        }
        
        list->count++;
    }

    /* Mark duplicates using sort-based algorithm */
    mark_duplicates_sorted(list);
    
    /* Second pass: write output */
    unsigned long out_pos = 0;
    
    /* Copy header first */
    unsigned long hdr_pos = 0;
    while (hdr_pos < size) {
        if (in_buf[hdr_pos] != '@')
            break;
        
        unsigned long line_start = hdr_pos;
        while (hdr_pos < size && in_buf[hdr_pos] != '\n')
            hdr_pos++;
        
        if (hdr_pos < size) hdr_pos++;
        
        unsigned long len = hdr_pos - line_start;
        if (out_pos + len > out_buf_capacity) {
            goto cleanup;
        }
        
        memcpy(out_buf + out_pos, in_buf + line_start, len);
        out_pos += len;
    }

    /* Write records */
    for (int i = 0; i < list->count; i++) {
        if (write_sam_record(in_buf, out_buf, &out_pos, out_buf_capacity, &list->records[i]) < 0) {
            goto cleanup;
        }
    }
    
    if (out_size)
        *out_size = out_pos;

    ret = 0;
    
cleanup:
    record_list_free(list);
    return ret;
}

// ==================== 从核入口函数 ====================

// 从核入口：每个 CPE 负责 paras[_PEN] 这一份
// 注意：函数名不带 slave_ 前缀，编译器会自动添加
void sam_process_cpe(SamProcessPara paras[64])
{
    int id = _PEN;
    SamProcessPara *para = &paras[id];

    char         *in_buf  = para->in_buf;
    char         *out_buf = para->out_buf;
    unsigned long size    = para->size;
    unsigned long out_buf_capacity = para->out_buf_capacity;
    int           mode    = para->mode;

    if (!in_buf || !out_buf || size == 0) {
        return;
    }
    
    // 如果没有设置 out_buf_capacity，默认等于 size
    if (out_buf_capacity == 0) {
        out_buf_capacity = size;
    }

    // 根据模式选择处理流程
    if (mode == MODE_SORT_ONLY) {
        // 仅排序
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
            if (out_pos + lines[i].len > size) break;
            memcpy(out_buf + out_pos,
                   in_buf  + lines[i].start,
                   (unsigned long)lines[i].len);
            out_pos += lines[i].len;
        }

        *(para->out_size) = out_pos;
        if (lines) free(lines);
        
    } else if (mode == MODE_MARKDUP_ONLY) {
        // 仅去重
        int ret = markdup_core(in_buf, out_buf, size, out_buf_capacity, para->out_size);
        if (ret != 0) {
            // 如果失败，确保 out_size 为 0
            *(para->out_size) = 0;
        }
        
    } else if (mode == MODE_ALL) {
        // 先排序再去重
        // 第一步：排序到 out_buf
        LineInfo *lines = 0;
        int n_lines = parse_sam_lines(in_buf, size, &lines);

        if (n_lines > 1 && lines) {
            quicksort_lineinfo(lines, 0, n_lines - 1);
        }

        unsigned long out_pos = 0;
        int i;
        for (i = 0; i < n_lines; ++i) {
            if (lines[i].len == 0) continue;
            if (out_pos + lines[i].len > size) break;
            memcpy(out_buf + out_pos,
                   in_buf  + lines[i].start,
                   (unsigned long)lines[i].len);
            out_pos += lines[i].len;
        }
        if (lines) free(lines);
        
        // 第二步：从 out_buf 去重回 in_buf，再复制回 out_buf
        // 注意：这里需要临时交换 in/out buffer
        // in_buf 和 out_buf 的容量相同，都是 out_buf_capacity
        unsigned long sorted_size = out_pos;
        unsigned long markdup_size = 0;
        int ret = markdup_core(out_buf, in_buf, sorted_size, out_buf_capacity, &markdup_size);
        
        // 复制回 out_buf
        if (ret == 0 && markdup_size > 0 && markdup_size <= out_buf_capacity) {
            memcpy(out_buf, in_buf, markdup_size);
            *(para->out_size) = markdup_size;
        } else {
            // 失败时，至少保留排序结果
            *(para->out_size) = sorted_size;
        }
    }
}
