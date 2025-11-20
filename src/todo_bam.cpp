// swbam_test.cpp
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <mutex>

// =============== Sunway Para & CPE 声明（你的工程里应已存在） ===============
typedef struct {
    char   *in_buffer;
    char   *out_buffer;
    size_t  in_size;
    size_t *out_size;
    int     level;     // 压缩等级，解压时可以忽略
} Para;

// 由 athread/SLAVE_FUN 生成的从核函数
extern "C" void slave_compressfunc(Para paras[64]);
extern "C" void slave_decompressfunc(Para paras[64]);

// 你的全局互斥锁（防止多个 CPU 线程同时调用 athread）
extern std::mutex globalMutex;

// 块大小：必须和 from核解压使用的 BLOCK_SIZE 一致
static const size_t SWBAM_BLOCK_SIZE = 2 * 1024 * 1024; // 2MB 例子

// =============== 简化版 header / record 结构体 ===============

struct swbam_hdr_t {
    int32_t n_ref;
    char  **ref_name;   // [n_ref]
    int32_t *ref_len;   // [n_ref]
    char   *text;
    int32_t l_text;
};

struct swbam1_t {
    int32_t tid;
    int32_t pos;
    int32_t mate_tid;
    int32_t mate_pos;
    int32_t tlen;
    uint16_t flag;
    uint8_t  mapq;

    uint32_t l_qname;
    uint32_t l_cigar;
    uint32_t l_seq;
    uint32_t l_qual;
    uint32_t l_aux;

    char    *qname;
    char    *cigar;
    char    *seq;
    char    *qual;
    uint8_t *aux;
};

// =============== swbam_hdr / swbam1 基本操作 ===============

swbam_hdr_t* swbam_hdr_init(int32_t n_ref) {
    swbam_hdr_t* h = (swbam_hdr_t*)std::calloc(1, sizeof(swbam_hdr_t));
    h->n_ref = n_ref;
    if (n_ref > 0) {
        h->ref_name = (char**)std::calloc(n_ref, sizeof(char*));
        h->ref_len  = (int32_t*)std::calloc(n_ref, sizeof(int32_t));
    }
    h->text = nullptr;
    h->l_text = 0;
    return h;
}

void swbam_hdr_destroy(swbam_hdr_t* h) {
    if (!h) return;
    if (h->ref_name) {
        for (int i = 0; i < h->n_ref; ++i) {
            std::free(h->ref_name[i]);
        }
        std::free(h->ref_name);
    }
    if (h->ref_len) std::free(h->ref_len);
    if (h->text) std::free(h->text);
    std::free(h);
}

swbam1_t* swbam1_init() {
    swbam1_t* b = (swbam1_t*)std::calloc(1, sizeof(swbam1_t));
    b->tid = -1;
    b->mate_tid = -1;
    return b;
}

void swbam1_destroy(swbam1_t* b) {
    if (!b) return;
    std::free(b->qname);
    std::free(b->cigar);
    std::free(b->seq);
    std::free(b->qual);
    std::free(b->aux);
    std::free(b);
}

swbam1_t* swbam1_dup(const swbam1_t* src) {
    if (!src) return nullptr;
    swbam1_t* b = swbam1_init();
    *b = *src; // 结构体浅拷贝
    // 把指针字段改成深拷贝
    b->qname = (char*)std::malloc(src->l_qname);
    std::memcpy(b->qname, src->qname, src->l_qname);
    b->cigar = (char*)std::malloc(src->l_cigar);
    std::memcpy(b->cigar, src->cigar, src->l_cigar);
    b->seq = (char*)std::malloc(src->l_seq);
    std::memcpy(b->seq, src->seq, src->l_seq);
    b->qual = (char*)std::malloc(src->l_qual);
    std::memcpy(b->qual, src->qual, src->l_qual);
    if (src->l_aux > 0) {
        b->aux = (uint8_t*)std::malloc(src->l_aux);
        std::memcpy(b->aux, src->aux, src->l_aux);
    } else {
        b->aux = nullptr;
    }
    return b;
}

// =============== 压缩/解压包装（调用 Sunway CPE） ===============

bool swbam_compress_block(const uint8_t* in, size_t in_size,
                          std::vector<uint8_t>& out,
                          int level)
{
    Para paras[64];
    size_t out_sizes[64] = {0};
    std::memset(paras, 0, sizeof(paras));

    // 输出缓冲区预估：略放大一点
    size_t bound = in_size * 2 + 1024;
    out.resize(bound);

    paras[0].in_buffer  = (char*)in;
    paras[0].out_buffer = (char*)out.data();
    paras[0].in_size    = in_size;
    paras[0].out_size   = &out_sizes[0];
    paras[0].level      = level;

    for (int i = 1; i < 64; ++i) {
        paras[i].in_buffer = nullptr;
        paras[i].in_size   = 0;
        paras[i].out_buffer= nullptr;
        paras[i].out_size  = &out_sizes[i];
        paras[i].level     = level;
    }

    {
        std::lock_guard<std::mutex> guard(globalMutex);
        __real_athread_spawn((void*)slave_compressfunc, paras, 1);
        athread_join();
    }

    if (out_sizes[0] == 0) return false;
    out.resize(out_sizes[0]);
    return true;
}

bool swbam_decompress_block(std::ifstream& fp,
                            uint64_t offset,
                            size_t comp_size,
                            uint8_t* out_buf,
                            size_t& out_size)
{
    std::vector<uint8_t> comp_buf(comp_size);

    fp.seekg(offset, std::ios::beg);
    fp.read((char*)comp_buf.data(), comp_size);
    if (!fp) return false;

    Para paras[64];
    size_t out_sizes[64] = {0};
    std::memset(paras, 0, sizeof(paras));

    paras[0].in_buffer  = (char*)comp_buf.data();
    paras[0].out_buffer = (char*)out_buf;
    paras[0].in_size    = comp_size;
    paras[0].out_size   = &out_sizes[0];

    for (int i = 1; i < 64; ++i) {
        paras[i].in_buffer = nullptr;
        paras[i].in_size   = 0;
        paras[i].out_buffer= nullptr;
        paras[i].out_size  = &out_sizes[i];
    }

    {
        std::lock_guard<std::mutex> guard(globalMutex);
        __real_athread_spawn((void*)slave_decompressfunc, paras, 1);
        athread_join();
    }

    out_size = out_sizes[0];
    return true;
}

// =============== header 序列化/反序列化 ===============

bool swbam_write_header(std::ofstream& fp, const swbam_hdr_t* h) {
    const char magic[8] = { 'S','W','B','A','M', 1, 0, 0 };
    fp.write(magic, 8);
    int32_t n_ref = h->n_ref;
    int32_t l_text = h->l_text;
    fp.write((char*)&n_ref, sizeof(int32_t));
    fp.write((char*)&l_text, sizeof(int32_t));
    if (l_text > 0 && h->text) {
        fp.write(h->text, l_text);
    }
    for (int i = 0; i < n_ref; ++i) {
        int32_t name_len = (int32_t)std::strlen(h->ref_name[i]);
        fp.write((char*)&name_len, sizeof(int32_t));
        fp.write(h->ref_name[i], name_len);
        fp.write((char*)&h->ref_len[i], sizeof(int32_t));
    }
    return (bool)fp;
}

bool swbam_read_header(std::ifstream& fp, swbam_hdr_t*& hdr_out) {
    char magic[8];
    fp.read(magic, 8);
    if (!fp) return false;
    if (std::memcmp(magic, "SWBAM", 5) != 0) {
        fprintf(stderr, "Not a SWBAM file.\n");
        return false;
    }
    int32_t n_ref = 0, l_text = 0;
    fp.read((char*)&n_ref, sizeof(int32_t));
    fp.read((char*)&l_text, sizeof(int32_t));
    if (!fp) return false;

    swbam_hdr_t* h = swbam_hdr_init(n_ref);
    h->l_text = l_text;
    if (l_text > 0) {
        h->text = (char*)std::malloc(l_text);
        fp.read(h->text, l_text);
    }

    for (int i = 0; i < n_ref; ++i) {
        int32_t name_len = 0;
        fp.read((char*)&name_len, sizeof(int32_t));
        h->ref_name[i] = (char*)std::malloc(name_len + 1);
        fp.read(h->ref_name[i], name_len);
        h->ref_name[i][name_len] = '\0';
        fp.read((char*)&h->ref_len[i], sizeof(int32_t));
    }

    if (!fp) {
        swbam_hdr_destroy(h);
        return false;
    }
    hdr_out = h;
    return true;
}

// =============== record 序列化/反序列化（块内） ===============

size_t swbam1_serialize(const swbam1_t* b, uint8_t* buf, size_t buf_len) {
    size_t need = sizeof(int32_t)*5 + sizeof(uint16_t) + 2 + sizeof(uint32_t)*5
                  + b->l_qname + b->l_cigar + b->l_seq + b->l_qual + b->l_aux;
    if (need > buf_len) return 0;

    uint8_t* p = buf;
    auto W32 = [&](int32_t v) { std::memcpy(p, &v, 4); p += 4; };
    auto W16 = [&](uint16_t v){ std::memcpy(p, &v, 2); p += 2; };
    auto W8  = [&](uint8_t  v){ *p++ = v; };

    W32(b->tid);
    W32(b->pos);
    W32(b->mate_tid);
    W32(b->mate_pos);
    W32(b->tlen);
    W16(b->flag);
    W8(b->mapq);
    W8(0); // reserved

    auto W32u = [&](uint32_t v){ std::memcpy(p, &v, 4); p += 4; };
    W32u(b->l_qname);
    W32u(b->l_cigar);
    W32u(b->l_seq);
    W32u(b->l_qual);
    W32u(b->l_aux);

    std::memcpy(p, b->qname, b->l_qname); p += b->l_qname;
    std::memcpy(p, b->cigar, b->l_cigar); p += b->l_cigar;
    std::memcpy(p, b->seq,   b->l_seq);   p += b->l_seq;
    std::memcpy(p, b->qual,  b->l_qual);  p += b->l_qual;
    if (b->l_aux > 0) {
        std::memcpy(p, b->aux, b->l_aux); p += b->l_aux;
    }
    return (size_t)(p - buf);
}

size_t swbam1_deserialize(uint8_t* buf, size_t buf_len, swbam1_t*& b_out) {
    if (buf_len < sizeof(int32_t)*5 + sizeof(uint16_t) + 2 + sizeof(uint32_t)*5)
        return 0;
    uint8_t* p = buf;
    auto R32 = [&](int32_t& v){ std::memcpy(&v, p, 4); p += 4; };
    auto R16 = [&](uint16_t& v){ std::memcpy(&v, p, 2); p += 2; };
    auto R8  = [&](uint8_t&  v){ v = *p++; };

    swbam1_t* b = swbam1_init();
    R32(b->tid);
    R32(b->pos);
    R32(b->mate_tid);
    R32(b->mate_pos);
    R32(b->tlen);
    R16(b->flag);
    R8(b->mapq);
    uint8_t reserved; R8(reserved);

    auto R32u = [&](uint32_t& v){ std::memcpy(&v, p, 4); p += 4; };
    R32u(b->l_qname);
    R32u(b->l_cigar);
    R32u(b->l_seq);
    R32u(b->l_qual);
    R32u(b->l_aux);

    size_t need = sizeof(int32_t)*5 + sizeof(uint16_t) + 2 + sizeof(uint32_t)*5
                  + b->l_qname + b->l_cigar + b->l_seq + b->l_qual + b->l_aux;
    if (need > buf_len) {
        swbam1_destroy(b);
        return 0;
    }

    b->qname = (char*)std::malloc(b->l_qname);
    std::memcpy(b->qname, p, b->l_qname); p += b->l_qname;

    b->cigar = (char*)std::malloc(b->l_cigar);
    std::memcpy(b->cigar, p, b->l_cigar); p += b->l_cigar;

    b->seq = (char*)std::malloc(b->l_seq);
    std::memcpy(b->seq, p, b->l_seq); p += b->l_seq;

    b->qual = (char*)std::malloc(b->l_qual);
    std::memcpy(b->qual, p, b->l_qual); p += b->l_qual;

    if (b->l_aux > 0) {
        b->aux = (uint8_t*)std::malloc(b->l_aux);
        std::memcpy(b->aux, p, b->l_aux); p += b->l_aux;
    } else {
        b->aux = nullptr;
    }

    b_out = b;
    return need;
}

// =============== swbam reader/writer 状态结构 ===============

struct SwbamIndex {
    size_t blocknum = 0;
    std::vector<size_t> comp_sizes;
    std::vector<uint64_t> offsets; // 相对于 data_start_offset 的偏移
};

struct swbam_reader_t {
    std::ifstream fp;
    swbam_hdr_t*  hdr = nullptr;

    SwbamIndex index;
    uint64_t   data_start_offset = 0;

    uint8_t*   ublock = nullptr;
    size_t     ublock_size = 0;

    std::vector<swbam1_t*> cur_recs;
    size_t     rec_idx = 0;
    size_t     cur_block = 0;

    swbam_reader_t(const char* path) : fp(path, std::ios::binary) {}
    ~swbam_reader_t() {
        if (hdr) swbam_hdr_destroy(hdr);
        if (ublock) std::free(ublock);
        for (auto* r : cur_recs) swbam1_destroy(r);
    }
};

struct swbam_writer_t {
    std::ofstream fp;
    swbam_hdr_t*  hdr = nullptr;
    int           level = 1;

    std::vector<size_t> comp_sizes;

    uint8_t* ublock = nullptr;
    size_t   ublock_size = 0;

    std::vector<swbam1_t*> recs_buf;

    swbam_writer_t(const char* path) : fp(path, std::ios::binary) {}
    ~swbam_writer_t() {
        if (hdr) swbam_hdr_destroy(hdr);
        if (ublock) std::free(ublock);
        for (auto* r : recs_buf) swbam1_destroy(r);
    }
};

// =============== 读取末尾 index（压缩块大小 + 偏移） ===============

bool swbam_read_index(const char* path, SwbamIndex& idx,
                      uint64_t& data_start_offset)
{
    std::ifstream iff(path, std::ios::binary | std::ios::ate);
    if (!iff.is_open()) return false;
    size_t file_size = (size_t)iff.tellg();
    if (file_size < sizeof(size_t)) return false;

    // blocknum 在文件末尾
    iff.seekg(-(long)sizeof(size_t), std::ios::end);
    size_t blocknum = 0;
    iff.read((char*)&blocknum, sizeof(size_t));
    idx.blocknum = blocknum;
    if (blocknum == 0) {
        data_start_offset = file_size; // 只有头，没有数据块
        return true;
    }

    size_t index_bytes = blocknum * sizeof(size_t);
    iff.seekg(-(long)(index_bytes + sizeof(size_t)), std::ios::end);

    idx.comp_sizes.resize(blocknum);
    size_t total_comp = 0;
    for (size_t i = 0; i < blocknum; ++i) {
        size_t bs = 0;
        iff.read((char*)&bs, sizeof(size_t));
        idx.comp_sizes[i] = bs;
        total_comp += bs;
    }

    uint64_t tail_bytes = index_bytes + sizeof(size_t);
    data_start_offset = file_size - total_comp - tail_bytes;

    idx.offsets.resize(blocknum);
    uint64_t off = 0;
    for (size_t i = 0; i < blocknum; ++i) {
        idx.offsets[i] = off;
        off += idx.comp_sizes[i];
    }

    return true;
}

// =============== swbam 打开/关闭/读写 ===============

swbam_writer_t* swbam_open_writer(const char* path,
                                  const swbam_hdr_t* hdr,
                                  int level)
{
    swbam_writer_t* w = new swbam_writer_t(path);
    if (!w->fp.is_open()) {
        delete w;
        return nullptr;
    }
    w->hdr = swbam_hdr_init(hdr->n_ref);
    // header 深拷贝（简单写一下）
    w->hdr->l_text = hdr->l_text;
    if (hdr->l_text > 0) {
        w->hdr->text = (char*)std::malloc(hdr->l_text);
        std::memcpy(w->hdr->text, hdr->text, hdr->l_text);
    }
    for (int i = 0; i < hdr->n_ref; ++i) {
        size_t len = std::strlen(hdr->ref_name[i]);
        w->hdr->ref_name[i] = (char*)std::malloc(len + 1);
        std::memcpy(w->hdr->ref_name[i], hdr->ref_name[i], len + 1);
        w->hdr->ref_len[i] = hdr->ref_len[i];
    }

    if (!swbam_write_header(w->fp, w->hdr)) {
        delete w;
        return nullptr;
    }

    w->level = level;
    w->ublock = (uint8_t*)std::malloc(SWBAM_BLOCK_SIZE);
    w->ublock_size = 0;
    return w;
}

bool swbam_writer_flush_block(swbam_writer_t* w) {
    if (w->ublock_size == 0) return true;
    std::vector<uint8_t> comp;
    if (!swbam_compress_block(w->ublock, w->ublock_size, comp, w->level))
        return false;
    w->fp.write((char*)comp.data(), comp.size());
    w->comp_sizes.push_back(comp.size());
    w->ublock_size = 0;
    for (auto* r : w->recs_buf) swbam1_destroy(r);
    w->recs_buf.clear();
    return (bool)w->fp;
}

int swbam_write1(swbam_writer_t* w, const swbam1_t* b) {
    swbam1_t* copy = swbam1_dup(b);
    if (!copy) return -1;
    w->recs_buf.push_back(copy);

    size_t written = 0;
    uint8_t* buf = w->ublock;
    written = 0;
    for (size_t i = 0; i < w->recs_buf.size(); ++i) {
        size_t ret = swbam1_serialize(w->recs_buf[i],
                                      buf + written,
                                      SWBAM_BLOCK_SIZE - written);
        if (ret == 0) {
            // 最后一条装不下了：去掉最后一条，flush 旧块，再单独写
            w->recs_buf.pop_back();
            swbam1_destroy(copy);

            written = 0;
            for (size_t j = 0; j < w->recs_buf.size(); ++j) {
                size_t r2 = swbam1_serialize(w->recs_buf[j],
                                             buf + written,
                                             SWBAM_BLOCK_SIZE - written);
                if (r2 == 0) return -1;
                written += r2;
            }
            w->ublock_size = written;
            if (!swbam_writer_flush_block(w)) return -1;

            // 新块写这条
            swbam1_t* copy2 = swbam1_dup(b);
            if (!copy2) return -1;
            w->recs_buf.push_back(copy2);
            size_t r3 = swbam1_serialize(copy2, w->ublock, SWBAM_BLOCK_SIZE);
            if (r3 == 0) return -1;
            w->ublock_size = r3;
            return 0;
        } else {
            written += ret;
        }
    }
    w->ublock_size = written;
    return 0;
}

void swbam_close_writer(swbam_writer_t* w) {
    if (!w) return;
    if (w->ublock_size > 0) {
        swbam_writer_flush_block(w);
    }
    // 写 index
    for (size_t sz : w->comp_sizes) {
        w->fp.write((char*)&sz, sizeof(size_t));
    }
    size_t blocknum = w->comp_sizes.size();
    w->fp.write((char*)&blocknum, sizeof(size_t));
    delete w;
}

swbam_reader_t* swbam_open_reader(const char* path) {
    swbam_reader_t* r = new swbam_reader_t(path);
    if (!r->fp.is_open()) {
        delete r;
        return nullptr;
    }

    // 先读 header
    if (!swbam_read_header(r->fp, r->hdr)) {
        delete r;
        return nullptr;
    }

    // 读 index
    if (!swbam_read_index(path, r->index, r->data_start_offset)) {
        delete r;
        return nullptr;
    }

    r->ublock = (uint8_t*)std::malloc(SWBAM_BLOCK_SIZE);
    r->ublock_size = 0;
    r->cur_block = 0;
    r->rec_idx = 0;
    return r;
}

bool swbam_load_next_block(swbam_reader_t* r) {
    for (auto* x : r->cur_recs) swbam1_destroy(x);
    r->cur_recs.clear();
    r->rec_idx = 0;

    if (r->cur_block >= r->index.blocknum) {
        return false;
    }

    size_t comp_size = r->index.comp_sizes[r->cur_block];
    uint64_t off = r->data_start_offset + r->index.offsets[r->cur_block];

    size_t out_size = 0;
    if (!swbam_decompress_block(r->fp, off, comp_size, r->ublock, out_size))
        return false;
    r->ublock_size = out_size;

    size_t used = 0;
    while (used < r->ublock_size) {
        swbam1_t* rec = nullptr;
        size_t consumed = swbam1_deserialize(r->ublock + used,
                                             r->ublock_size - used,
                                             rec);
        if (consumed == 0) break;
        used += consumed;
        r->cur_recs.push_back(rec);
    }

    r->cur_block++;
    r->rec_idx = 0;
    return !r->cur_recs.empty();
}

int swbam_read1(swbam_reader_t* r, swbam1_t* b) {
    while (true) {
        if (r->rec_idx < r->cur_recs.size()) {
            swbam1_t* src = r->cur_recs[r->rec_idx++];
            swbam1_t* tmp = swbam1_dup(src);
            *b = *tmp;
            // 把 tmp 的指针直接“塞”给 b，再释放临时结构体壳
            std::free(tmp);
            return 1;
        }
        if (!swbam_load_next_block(r)) {
            return 0; // EOF
        }
    }
}

swbam_hdr_t* swbam_get_header(swbam_reader_t* r) {
    return r ? r->hdr : nullptr;
}

void swbam_close_reader(swbam_reader_t* r) {
    if (!r) return;
    delete r;
}

// =============== main：从参数读文件名，简单测试 ===============
//
// 用法：
//   1) 生成一个测试文件：
//      ./swbam_test write test.sw.bam
//   2) 读这个文件并打印：
//      ./swbam_test read test.sw.bam
//
int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "Usage:\n"
                     "  %s write <output.sw.bam>\n"
                     "  %s read  <input.sw.bam>\n",
                     argv[0], argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    const char* path = argv[2];

    if (mode == "write") {
        // 构造一个假的 header：1 条参考序列 chr1，长度 1e6
        swbam_hdr_t* hdr = swbam_hdr_init(1);
        const char* name = "chr1";
        hdr->ref_name[0] = (char*)std::malloc(std::strlen(name) + 1);
        std::strcpy(hdr->ref_name[0], name);
        hdr->ref_len[0] = 1000000;
        hdr->text = nullptr;
        hdr->l_text = 0;

        swbam_writer_t* w = swbam_open_writer(path, hdr, 3);
        if (!w) {
            std::fprintf(stderr, "Failed to open writer for %s\n", path);
            swbam_hdr_destroy(hdr);
            return 1;
        }

        // 构造几条假记录
        for (int i = 0; i < 10; ++i) {
            swbam1_t* b = swbam1_init();
            b->tid = 0;
            b->pos = i * 100;
            b->mate_tid = -1;
            b->mate_pos = -1;
            b->tlen = 0;
            b->flag = 0;
            b->mapq = 60;

            std::string qn = "read_" + std::to_string(i);
            std::string cg = "100M";
            std::string sq = "ACGTACGTAC"; // just toy
            std::string qu = "IIIIIIIIII";

            b->l_qname = qn.size();
            b->l_cigar = cg.size();
            b->l_seq   = sq.size();
            b->l_qual  = qu.size();
            b->l_aux   = 0;

            b->qname = (char*)std::malloc(b->l_qname);
            std::memcpy(b->qname, qn.data(), b->l_qname);
            b->cigar = (char*)std::malloc(b->l_cigar);
            std::memcpy(b->cigar, cg.data(), b->l_cigar);
            b->seq = (char*)std::malloc(b->l_seq);
            std::memcpy(b->seq, sq.data(), b->l_seq);
            b->qual = (char*)std::malloc(b->l_qual);
            std::memcpy(b->qual, qu.data(), b->l_qual);
            b->aux = nullptr;

            if (swbam_write1(w, b) != 0) {
                std::fprintf(stderr, "swbam_write1 failed at %d\n", i);
            }
            swbam1_destroy(b);
        }

        swbam_close_writer(w);
        swbam_hdr_destroy(hdr);
        std::fprintf(stderr, "Wrote test sw.bam to %s\n", path);
    } else if (mode == "read") {
        swbam_reader_t* r = swbam_open_reader(path);
        if (!r) {
            std::fprintf(stderr, "Failed to open reader for %s\n", path);
            return 1;
        }
        swbam_hdr_t* hdr = swbam_get_header(r);
        std::fprintf(stderr, "Header: n_ref=%d\n", hdr->n_ref);
        for (int i = 0; i < hdr->n_ref; ++i) {
            std::fprintf(stderr, "  ref[%d]: %s len=%d\n",
                         i, hdr->ref_name[i], hdr->ref_len[i]);
        }

        swbam1_t* b = swbam1_init();
        int cnt = 0;
        while (true) {
            int ret = swbam_read1(r, b);
            if (ret <= 0) break;
            std::string qn(b->qname, b->l_qname);
            std::string cg(b->cigar, b->l_cigar);
            std::string sq(b->seq,   b->l_seq);
            std::string qu(b->qual,  b->l_qual);
            std::printf("REC %d: tid=%d pos=%d mapq=%d qname=%s cigar=%s seq=%s qual=%s\n",
                        cnt, b->tid, b->pos, b->mapq,
                        qn.c_str(), cg.c_str(), sq.c_str(), qu.c_str());
            ++cnt;
            // 清理 b 的指针，下次再重用 b 需要先 free
            std::free(b->qname); b->qname = nullptr; b->l_qname = 0;
            std::free(b->cigar); b->cigar = nullptr; b->l_cigar = 0;
            std::free(b->seq);   b->seq   = nullptr; b->l_seq   = 0;
            std::free(b->qual);  b->qual  = nullptr; b->l_qual  = 0;
            std::free(b->aux);   b->aux   = nullptr; b->l_aux   = 0;
        }
        swbam1_destroy(b);
        swbam_close_reader(r);
        std::fprintf(stderr, "Read %d records from %s\n", cnt, path);
    } else {
        std::fprintf(stderr, "Unknown mode: %s (use 'write' or 'read')\n", mode.c_str());
        return 1;
    }

    return 0;
}
