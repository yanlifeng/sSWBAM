// sam_process_para.h
// 公共参数结构，主核和从核都会使用

#ifndef SAM_PROCESS_PARA_H
#define SAM_PROCESS_PARA_H

// 处理模式
#define MODE_SORT_ONLY      1
#define MODE_MARKDUP_ONLY   2
#define MODE_ALL            3  // sort + markdup

typedef struct {
    char   *in_buf;       // 输入 SAM buffer
    char   *out_buf;      // 输出 SAM buffer（处理后的结果）
    unsigned long size;   // 输入有效字节数
    unsigned long out_buf_capacity; // 输出 buffer 容量（可能大于 size）
    unsigned long *out_size; // CPE 处理完后写回输出长度
    int     mode;         // 处理模式：MODE_SORT_ONLY, MODE_MARKDUP_ONLY, MODE_ALL
} SamProcessPara;

#endif // SAM_PROCESS_PARA_H

