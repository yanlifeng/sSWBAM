// sam_sort_para.h
// 公共参数结构，主核和从核都会使用

#ifndef SAM_SORT_PARA_H
#define SAM_SORT_PARA_H

typedef struct {
    char   *in_buf;       // 输入 SAM buffer
    char   *out_buf;      // 输出 SAM buffer（排序后的结果）
    unsigned long size;   // 输入有效字节数
    unsigned long *out_size; // CPE 处理完后写回输出长度
} SamSortPara;

#endif // SAM_SORT_PARA_H

