/*  Simplified markdup for Sunway CPE
    
    Extracts core duplicate marking algorithm without htslib dependencies.
    Memory-optimized version using sorting instead of hash table.
    
    Based on samtools bam_markdup.c
    Copyright (C) 2017-2023 Genome Research Ltd.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <slave.h>

/* Keep consistent with definition in main.c on MPE side */
typedef struct {
    char *in_buf;
    char *out_buf;
    unsigned long size;
    unsigned long *out_size;
} SamMarkPara;

/* Buffer reader for CPE-side SAM parsing */
typedef struct {
    const char *buf;
    unsigned long size;
    unsigned long pos;
} sam_buf_reader_t;

/* Compact SAM record - only store what we need
 * Memory: ~32 bytes per record (vs 68 bytes before)
 */
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

/* Reference name to ID mapping */
#define MAX_REFS 256

typedef struct {
    char names[MAX_REFS][256];
    int count;
} ref_map_t;

typedef struct {
    sam_record_t *records;
    int count;
    int capacity;
    ref_map_t ref_map;
} record_list_t;

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

/* Get or create reference ID from name */
int get_ref_id(ref_map_t *map, const char *rname, int len) {
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
record_list_t *record_list_init(void) {
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
void record_list_free(record_list_t *list) {
    if (!list) return;
    free(list->records);
    free(list);
}

/* Calculate quality score from quality string */
int64_t calc_score(const char *qual, int len) {
    int64_t score = 0;
    for (int i = 0; i < len; i++) {
        int q = qual[i] - 33;  /* Phred+33 */
        if (q > 15) q = 15;    /* Cap at 15 as per markdup spec */
        if (q > 0) score += q;
    }
    return score;
}

/* Parse one line of SAM - extract only fields needed for duplicate detection */
int parse_sam_line(const char *line, unsigned long line_offset, int len, 
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

/* Read one SAM line from buffer, skip headers */
int read_sam_line(sam_buf_reader_t *r, const char **line_start, 
                  int *line_len, unsigned long *line_offset) {
    if (!r || r->pos >= r->size) return -1;
    
    /* Skip empty lines and headers */
    while (r->pos < r->size) {
        while (r->pos < r->size && (r->buf[r->pos] == '\n' || r->buf[r->pos] == '\r'))
            r->pos++;
        
        if (r->pos >= r->size) return -1;
        
        /* Skip header lines */
        if (r->buf[r->pos] == '@') {
            while (r->pos < r->size && r->buf[r->pos] != '\n')
                r->pos++;
            continue;
        }
        
        break;
    }
    
    if (r->pos >= r->size) return -1;
    
    *line_offset = r->pos;
    *line_start = r->buf + r->pos;
    unsigned long start = r->pos;
    
    while (r->pos < r->size && r->buf[r->pos] != '\n' && r->buf[r->pos] != '\r')
        r->pos++;
    
    *line_len = r->pos - start;
    
    if (r->pos < r->size && r->buf[r->pos] == '\n')
        r->pos++;
    
    return (*line_len > 0) ? 0 : -1;
}

/* Comparison function for sorting records by position */
int compare_records(const void *a, const void *b) {
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
void mark_duplicates_sorted(record_list_t *list) {
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
int write_sam_record(const char *in_buf, char *out_buf, 
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

/* Core markdup function for CPE */
int markdup_core(const char *in_buf, char *out_buf, 
                 unsigned long size, unsigned long *out_size) {
    sam_buf_reader_t reader;
    record_list_t *list = NULL;
    const char *line;
    int line_len;
    unsigned long line_offset;
    int ret = -1;
    
    if (!in_buf || !out_buf || size == 0) {
        return -1;
    }

    reader.buf = in_buf;
    reader.size = size;
    reader.pos = 0;
    
    /* Initialize record list */
    list = record_list_init();
    if (!list) {
        return -1;
    }

    /* First pass: read all records */
    while (read_sam_line(&reader, &line, &line_len, &line_offset) == 0) {
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
        if (parse_sam_line(line, line_offset, line_len, rec, &list->ref_map) < 0) {
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
        if (out_pos + len > size) {
            goto cleanup;
        }
        
        memcpy(out_buf + out_pos, in_buf + line_start, len);
        out_pos += len;
    }

    /* Write records */
    for (int i = 0; i < list->count; i++) {
        if (write_sam_record(in_buf, out_buf, &out_pos, size, &list->records[i]) < 0) {
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

/* CPE entry function */
void slave_sam_markdup_cpe(void *arg)
{
    SamMarkPara *paras = (SamMarkPara *)arg;
    int id = _PEN;

    if (!paras)
        return;

    if (id < 0 || id >= 64)
        return;

    /* Each CPE processes its corresponding para */
    SamMarkPara *p = &paras[id];

    if (!p->in_buf || !p->out_buf || p->size == 0) {
        if (p->out_size)
            *(p->out_size) = 0;
        return;
    }

    /* Call core markdup function */
    markdup_core(p->in_buf, p->out_buf, p->size, p->out_size);
}
