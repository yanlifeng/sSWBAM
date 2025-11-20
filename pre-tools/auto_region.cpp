// auto_region.cpp
// 用法：
//   g++ -O3 -std=gnu++11 -fopenmp auto_region.cpp -o auto_region
//   ./auto_region ref.fa in.sam out_dir
//
// 功能：
//   1. 从 ref.fa 中解析出 chr 名和长度（只保留 chr1-22, chrX, chrY）。
//   2. 整个 in.sam 读入内存 sam_buf，一次遍历：
//        - 收集 header 行；
//        - 解析每个 alignment，记录 chr_id, pos, offset, len；
//        - 按 chr/bin 累积“字节权重” bin_weight，用于估计每个区间的 SAM 字节数。
//   3. 对每个 chr，用 bin_weight + 目标大小（默认 2 MB）切成若干 region[start,end]，不跨 chr。
//   4. 构建 per-chr 的 record 索引列表 chr_rec_indices[chr]。
//   5. 使用 OpenMP 按 chr 并行：每个线程只处理一个 chr，把属于它的记录写到相应 region 的 SAM 文件里：
//        out_dir/chrX_start_end.sam
//      每个 region 文件会写完整 SAM header。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// ---------------- 可调参数 ----------------
static const double TARGET_REGION_MB = 64.0;         // 目标每个 region 文件大小（MB）
static const int    BIN_SIZE        = 1000;       // 深度离散 bin 大小（bp）

// ---------------- 计时辅助 ----------------
static double now_ms()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

// ---------------- 结构体定义 ----------------
struct Region {
    int64_t start;   // 1-based, inclusive
    int64_t end;     // inclusive
};

struct ChrInfo {
    std::string name;
    int64_t     length;

    int64_t     num_bins;
    std::vector<double> bin_weight;   // 每个 bin 累积的字节数

    std::vector<Region> regions;      // 最终切出来的 region 列表
};

// SAM 中的每一条 alignment 记录（不含 header）
struct SamRecord {
    int     chr_id;     // chrs 下标
    int64_t pos;        // 1-based POS
    size_t  offset;     // 在 sam_buf 中的行起始 offset
    size_t  len;        // 该行长度（包括 '\n'）
};

// ---------------- 工具：判断是否为需要的 chr 名 ----------------
// 只保留：chr1-22, chrX, chrY
static bool is_wanted_chr_name(const std::string& name)
{
    if (name == "chrX" || name == "chrY") return true;
    if (name.size() <= 3) return false;
    if (name.compare(0, 3, "chr") != 0) return false;

    int v = 0;
    int i = 3;
    if (i >= (int)name.size()) return false;
    for (; i < (int)name.size(); ++i) {
        char c = name[i];
        if (c < '0' || c > '9') return false;
        v = v * 10 + (c - '0');
        if (v > 22) return false;  // 超过 22 不要
    }
    if (v >= 1 && v <= 22) return true;
    return false;
}

// ---------------- FASTA 解析：只保留 chr1-22, chrX, chrY ----------------
static bool parse_fasta(const std::string& fasta_path,
                        std::vector<ChrInfo>& chrs,
                        std::unordered_map<std::string,int>& chr_index)
{
    FILE* fp = std::fopen(fasta_path.c_str(), "rb");
    if (!fp) {
        std::fprintf(stderr, "Failed to open fasta: %s (%s)\n",
                     fasta_path.c_str(), std::strerror(errno));
        return false;
    }

    char*  line = NULL;
    size_t cap  = 0;

    ChrInfo cur;
    bool    in_seq       = false;
    bool    keep_current = false;
    size_t  skipped_chrs = 0;

    while (true) {
        ssize_t n = getline(&line, &cap, fp);
        if (n < 0) break;

        // 去掉行尾 \n/\r
        if (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) {
            while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) {
                line[--n] = '\0';
            }
        }
        if (n == 0) continue;

        if (line[0] == '>') {
            if (in_seq && keep_current) {
                chrs.push_back(cur);
            }
            cur = ChrInfo();
            in_seq = true;

            int i = 1;
            while (line[i] == ' ' || line[i] == '\t') ++i;
            int start = i;
            while (line[i] != '\0' && line[i] != ' ' && line[i] != '\t') ++i;
            int end = i;

            cur.name.assign(line + start, end - start);
            cur.length = 0;

            keep_current = is_wanted_chr_name(cur.name);
            if (!keep_current) skipped_chrs++;
        } else {
            if (!in_seq) continue;
            if (!keep_current) continue;
            for (int i = 0; i < n; ++i) {
                char c = line[i];
                if (c != '\n' && c != '\r' && c != ' ' && c != '\t')
                    cur.length++;
            }
        }
    }

    if (in_seq && keep_current) {
        chrs.push_back(cur);
    }

    std::fclose(fp);
    if (line) std::free(line);

    if (chrs.empty()) {
        std::fprintf(stderr,
                     "No wanted chromosomes (chr1-22, chrX, chrY) found in fasta: %s\n",
                     fasta_path.c_str());
        return false;
    }

    for (size_t i = 0; i < chrs.size(); ++i) {
        ChrInfo& c = chrs[i];
        if (c.length <= 0) {
            c.num_bins = 0;
        } else {
            c.num_bins = (c.length + BIN_SIZE - 1) / BIN_SIZE;
        }
        c.bin_weight.assign(c.num_bins, 0.0);
        chr_index[c.name] = (int)i;
    }

    std::fprintf(stderr,
                 "Parsed FASTA: kept %zu chromosomes (chr1-22, chrX, chrY), skipped %zu others\n",
                 chrs.size(), skipped_chrs);
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::fprintf(stderr, "  chr[%zu]: %s len=%ld bins=%ld\n",
                     i, chrs[i].name.c_str(),
                     (long)chrs[i].length,
                     (long)chrs[i].num_bins);
    }
    return true;
}

// ---------------- SAM 行解析：RNAME + POS ----------------
static bool parse_sam_rname_pos(const char* line,
                                size_t len,       // 不含 '\n'
                                const char*& rname_start,
                                size_t& rname_len,
                                int64_t& pos_out)
{
    if (len == 0) return false;
    if (line[0] == '@') return false;

    int  field = 0;
    int  i     = 0;
    int  start = 0;
    int  end   = 0;

    const char* rname_s = NULL;
    int         rname_l = 0;
    const char* pos_s   = NULL;
    int         pos_l   = 0;

    for (i = 0; i <= (int)len; ++i) {
        char c = (i == (int)len) ? '\t' : line[i];

        if (c == '\t') {
            end = i;
            if (field == 2) {
                rname_s = line + start;
                rname_l = end - start;
            } else if (field == 3) {
                pos_s = line + start;
                pos_l = end - start;
                break;
            }
            field++;
            start = i + 1;
        }
    }

    if (!rname_s || !pos_s || pos_l <= 0) return false;

    int64_t value = 0;
    int     neg   = 0;
    int     idx   = 0;
    if (pos_s[0] == '-') {
        neg = 1;
        idx++;
        if (idx >= pos_l) return false;
    }
    for (; idx < pos_l; ++idx) {
        char c = pos_s[idx];
        if (c < '0' || c > '9') return false;
        value = value * 10 + (c - '0');
    }
    if (neg) value = -value;

    rname_start = rname_s;
    rname_len   = (size_t)rname_l;
    pos_out     = value;
    return true;
}

// ---------------- 整个 SAM 读入内存并解析 ----------------
//
// 输出：
//   sam_buf/sam_size        : 整个 SAM 文件内容
//   header_lines            : SAM header 行（含 '\n'）
//   records                 : 所有 chr1-22/X/Y 的 alignment 记录
//   同时更新 chrs[].bin_weight[bin] += line_len_bytes
//
static bool load_and_parse_sam(const std::string& sam_path,
                               char*& sam_buf,
                               size_t& sam_size,
                               const std::unordered_map<std::string,int>& chr_index,
                               std::vector<ChrInfo>& chrs,
                               std::vector<std::string>& header_lines,
                               std::vector<SamRecord>& records)
{
    struct stat st;
    if (stat(sam_path.c_str(), &st) != 0) {
        std::fprintf(stderr, "stat failed for SAM: %s (%s)\n",
                     sam_path.c_str(), std::strerror(errno));
        return false;
    }

    if (st.st_size <= 0) {
        std::fprintf(stderr, "Empty SAM file: %s\n", sam_path.c_str());
        return false;
    }

    sam_size = (size_t)st.st_size;
    sam_buf  = (char*)std::malloc(sam_size);
    if (!sam_buf) {
        std::fprintf(stderr, "malloc sam_buf failed, size=%zu\n", sam_size);
        return false;
    }

    FILE* fp = std::fopen(sam_path.c_str(), "rb");
    if (!fp) {
        std::fprintf(stderr, "Failed to open SAM: %s (%s)\n",
                     sam_path.c_str(), std::strerror(errno));
        std::free(sam_buf);
        sam_buf = nullptr;
        return false;
    }

    size_t nread = std::fread(sam_buf, 1, sam_size, fp);
    std::fclose(fp);
    if (nread != sam_size) {
        std::fprintf(stderr, "fread SAM incomplete: expect=%zu got=%zu\n",
                     sam_size, nread);
        std::free(sam_buf);
        sam_buf = nullptr;
        return false;
    }

    // 解析行
    size_t i = 0;
    int64_t total_reads = 0;
    int64_t used_reads  = 0;

    while (i < sam_size) {
        size_t line_start = i;
        while (i < sam_size && sam_buf[i] != '\n') ++i;
        size_t line_end = i;          // 不含 '\n'
        if (i < sam_size && sam_buf[i] == '\n') ++i;
        size_t line_len = i - line_start;  // 行长度（含 '\n'，最后一行没换行也行）

        if (line_len == 0) continue;

        char*  line_ptr  = sam_buf + line_start;
        size_t text_len  = line_end - line_start;

        if (line_ptr[0] == '@') {
            // header 行，直接拷一份
            header_lines.emplace_back(line_ptr, line_len);
            continue;
        }

        total_reads++;

        const char* rname_s = NULL;
        size_t      rname_l = 0;
        int64_t     pos     = 0;
        if (!parse_sam_rname_pos(line_ptr, text_len, rname_s, rname_l, pos)) {
            continue;
        }

        std::string rname(rname_s, rname_l);
        auto it = chr_index.find(rname);
        if (it == chr_index.end()) {
            // 不在 chr1-22/X/Y 内
            continue;
        }
        int chr_id = it->second;
        ChrInfo& c = chrs[chr_id];

        if (pos <= 0 || pos > c.length) {
            continue;
        }

        // 更新该 chr 的 bin_weight
        if (c.num_bins > 0) {
            int bin_idx = (int)((pos - 1) / BIN_SIZE);
            if (bin_idx < 0) bin_idx = 0;
            if (bin_idx >= (int)c.num_bins) bin_idx = (int)c.num_bins - 1;
            c.bin_weight[bin_idx] += (double)line_len;
        }

        // 记录 SamRecord
        SamRecord rec;
        rec.chr_id = chr_id;
        rec.pos    = pos;
        rec.offset = line_start;
        rec.len    = line_len;
        records.push_back(rec);
        used_reads++;
    }

    std::fprintf(stderr,
                 "SAM loaded into memory: size=%.3f MB, total_reads=%ld, used_reads=%ld\n",
                 sam_size / 1024.0 / 1024.0,
                 (long)total_reads, (long)used_reads);
    return true;
}

static void build_regions_for_chr(ChrInfo& c, double target_bytes)
{
    if (c.length <= 0) return;

    if (c.num_bins == 0 || c.bin_weight.empty()) {
        // 没有 bin 信息，整个 chr 一个 region
        Region r; r.start = 1; r.end = c.length;
        c.regions.push_back(r);
        return;
    }

    double   accum_bytes   = 0.0;     // 当前 region 已累积的“字节”
    int64_t  region_start  = 1;       // 当前 region 的起始坐标（1-based）
    int64_t  chr_len       = c.length;

    for (int64_t b = 0; b < c.num_bins; ++b) {
        double w = c.bin_weight[b];

        int64_t bin_start_pos = b * (int64_t)BIN_SIZE + 1;
        int64_t bin_end_pos   = (b + 1) * (int64_t)BIN_SIZE;
        if (bin_end_pos > chr_len) bin_end_pos = chr_len;

        if (region_start > chr_len)
            break;

        // 如果把这个 bin 加进来就 >= target_bytes，
        // 那就直接在这个 bin 的末尾结束一个 region，允许略微超过 target_bytes。
        if (accum_bytes + w >= target_bytes) {
            Region r;
            r.start = region_start;
            r.end   = bin_end_pos;  // 这个 bin 的尾巴作为 region 结束
            c.regions.push_back(r);

            region_start = bin_end_pos + 1;
            accum_bytes  = 0.0;     // 下一个 region 从下一个位置重新累计
        } else {
            // 还没到 target，继续累计
            accum_bytes += w;
        }
    }

    // chr 尾巴部分，如果还有没覆盖的区间，就作为最后一个 region
    if (region_start <= chr_len) {
        Region r;
        r.start = region_start;
        r.end   = chr_len;
        c.regions.push_back(r);
    }

    std::fprintf(stderr, "  chr %s: regions=%zu\n",
                 c.name.c_str(), c.regions.size());
}


// ---------------- 对每个 chr 划分 region ----------------
static void build_regions_for_chr_init(ChrInfo& c, double target_bytes)
{
    if (c.length <= 0) return;

    if (c.num_bins == 0 || c.bin_weight.empty()) {
        // 没有 bin 信息，简单一个 region 覆盖全部
        Region r; r.start = 1; r.end = c.length;
        c.regions.push_back(r);
        return;
    }

    double   accum_bytes   = 0.0;
    int64_t  region_start  = 1;
    int64_t  chr_len       = c.length;

    for (int64_t b = 0; b < c.num_bins; ++b) {
        double w = c.bin_weight[b];
        accum_bytes += w;

        int64_t bin_start_pos = b * (int64_t)BIN_SIZE + 1;
        int64_t bin_end_pos   = (b + 1) * (int64_t)BIN_SIZE;
        if (bin_end_pos > chr_len) bin_end_pos = chr_len;

        while (accum_bytes >= target_bytes) {
            Region r;
            r.start = region_start;
            r.end   = bin_end_pos;  // 简单地切在 bin 边界
            c.regions.push_back(r);

            region_start = bin_end_pos + 1;
            if (region_start > chr_len) {
                accum_bytes = 0.0;
                break;
            }

            accum_bytes -= target_bytes;
        }

        if (region_start > chr_len) break;
    }

    if (region_start <= chr_len) {
        Region r;
        r.start = region_start;
        r.end   = chr_len;
        c.regions.push_back(r);
    }

    std::fprintf(stderr, "  chr %s: regions=%zu\n",
                 c.name.c_str(), c.regions.size());
}

// ---------------- 按 chr 构建 per-chr record 列表 ----------------
static void build_chr_record_indices(const std::vector<SamRecord>& records,
                                     size_t n_chr,
                                     std::vector< std::vector<int> >& chr_rec_indices)
{
    chr_rec_indices.clear();
    chr_rec_indices.resize(n_chr);
    for (size_t i = 0; i < records.size(); ++i) {
        const SamRecord& rec = records[i];
        if (rec.chr_id < 0 || rec.chr_id >= (int)n_chr) continue;
        chr_rec_indices[rec.chr_id].push_back((int)i);
    }

    for (size_t c = 0; c < n_chr; ++c) {
        std::fprintf(stderr,
                     "  chr_rec_indices[%zu] = %zu records\n",
                     c, chr_rec_indices[c].size());
    }
}

// ---------------- 并行按 chr split：从 sam_buf 直接写 region 文件（每个region一块buffer） ----------------
static bool split_sam_by_regions_from_memory(const char* sam_buf,
                                             const std::string& out_dir,
                                             const std::vector<ChrInfo>& chrs,
                                             const std::vector<std::string>& header_lines,
                                             const std::vector<SamRecord>& records,
                                             const std::vector< std::vector<int> >& chr_rec_indices)
{
    int n_chr = (int)chrs.size();
    if (n_chr == 0) {
        std::fprintf(stderr, "No chromosomes to split.\n");
        return true;
    }

    int global_ok = 1;

    // 按 chr 并行
    #pragma omp parallel for schedule(dynamic)
    for (int chr_id = 0; chr_id < n_chr; ++chr_id) {
        const ChrInfo& c = chrs[chr_id];
        const std::vector<int>& idx_list = chr_rec_indices[chr_id];

        if (c.regions.empty() || idx_list.empty()) {
            continue;
        }

        size_t n_region = c.regions.size();
        // 每个 region 一块 buffer：保存属于该 region 的 SamRecord 下标
        std::vector< std::vector<int> > reg_rec_ids(n_region);

        // 先按 region 把 record 分桶
        for (size_t k = 0; k < idx_list.size(); ++k) {
            int rec_id = idx_list[k];
            const SamRecord& rec = records[rec_id];

            int64_t pos = rec.pos;
            if (pos <= 0 || pos > c.length) continue;

            // 二分查找该 pos 属于哪个 region
            int left = 0;
            int right = (int)n_region - 1;
            int found = -1;

            while (left <= right) {
                int mid = (left + right) / 2;
                const Region& rg = c.regions[mid];
                if (pos < rg.start) {
                    right = mid - 1;
                } else if (pos > rg.end) {
                    left = mid + 1;
                } else {
                    found = mid;
                    break;
                }
            }

            if (found >= 0) {
                reg_rec_ids[found].push_back(rec_id);
            }
            // 找不到的就直接忽略（理论上不会太多）
        }

        // 再按 region 输出文件
        int64_t written_reads = 0;

        for (size_t r_idx = 0; r_idx < n_region; ++r_idx) {
            const std::vector<int>& rec_ids = reg_rec_ids[r_idx];
            if (rec_ids.empty()) continue;  // 这个 region 没有 read，跳过

            const Region& rg = c.regions[r_idx];

            char path[4096];
            std::snprintf(path, sizeof(path),
                          "%s/%s_%ld_%ld.sam",
                          out_dir.c_str(), c.name.c_str(),
                          (long)rg.start, (long)rg.end);

            FILE* out = std::fopen(path, "wb");
            if (!out) {
                std::fprintf(stderr,
                             "[OMP] Failed to open region SAM: %s (%s)\n",
                             path, std::strerror(errno));
                // 简单处理一下错误标记
                #pragma omp critical
                {
                    global_ok = 0;
                }
                continue;
            }

            // 写 header
            for (const auto& hline : header_lines) {
                std::fwrite(hline.data(), 1, hline.size(), out);
            }

            // 写属于这个 region 的所有 record（按它们在原 SAM 里的顺序）
            for (size_t t = 0; t < rec_ids.size(); ++t) {
                int rec_id = rec_ids[t];
                const SamRecord& rec = records[rec_id];
                std::fwrite(sam_buf + rec.offset, 1, rec.len, out);
                written_reads++;
            }

            std::fclose(out);
        }

        std::fprintf(stderr,
                     "[OMP] chr %s done, regions=%zu, written_reads=%ld\n",
                     c.name.c_str(), c.regions.size(), (long)written_reads);
    }

    return global_ok != 0;
}


// ---------------- 并行按 chr split：从 sam_buf 直接写 region 文件 ----------------
static bool split_sam_by_regions_from_memory_init(const char* sam_buf,
                                             const std::string& out_dir,
                                             const std::vector<ChrInfo>& chrs,
                                             const std::vector<std::string>& header_lines,
                                             const std::vector<SamRecord>& records,
                                             const std::vector< std::vector<int> >& chr_rec_indices)
{
    int n_chr = (int)chrs.size();
    if (n_chr == 0) {
        std::fprintf(stderr, "No chromosomes to split.\n");
        return true;
    }

    int global_ok = 1;

    struct RegionFile {
        FILE* fp;
        bool  header_written;
        RegionFile() : fp(nullptr), header_written(false) {}
    };

    // 按 chr 并行
    #pragma omp parallel for schedule(dynamic)
    for (int chr_id = 0; chr_id < n_chr; ++chr_id) {
        const ChrInfo& c = chrs[chr_id];
        const std::vector<int>& idx_list = chr_rec_indices[chr_id];

        if (c.regions.empty() || idx_list.empty()) {
            continue;
        }

        std::vector<RegionFile> region_files(c.regions.size());

        auto close_region_files = [&]() {
            for (size_t i = 0; i < region_files.size(); ++i) {
                if (region_files[i].fp) {
                    std::fclose(region_files[i].fp);
                    region_files[i].fp = nullptr;
                }
            }
        };

        size_t region_idx = 0;
        int64_t written_reads = 0;

        // 假设同一个 chr 内记录大致按 pos 升序
        for (size_t k = 0; k < idx_list.size(); ++k) {
            int rec_id = idx_list[k];
            const SamRecord& rec = records[rec_id];

            int64_t pos = rec.pos;
            if (pos <= 0 || pos > c.length) continue;

            while (region_idx < c.regions.size() &&
                   pos > c.regions[region_idx].end) {
                region_idx++;
            }
            if (region_idx >= c.regions.size()) {
                // 后面的记录都在最后一个 region 之后，直接结束
                break;
            }
            const Region& rg = c.regions[region_idx];
            if (pos < rg.start || pos > rg.end) {
                // 理论上不会出现（pos 应在当前 region 内），出现就略过
                continue;
            }

            RegionFile& rf = region_files[region_idx];
            if (!rf.fp) {
                char path[4096];
                std::snprintf(path, sizeof(path),
                              "%s/%s_%ld_%ld.sam",
                              out_dir.c_str(), c.name.c_str(),
                              (long)rg.start, (long)rg.end);
                rf.fp = std::fopen(path, "wb");
                if (!rf.fp) {
                    std::fprintf(stderr,
                                 "[OMP] Failed to open region SAM: %s (%s)\n",
                                 path, std::strerror(errno));
                    global_ok = 0;
                    continue;
                }
                // 写 header
                for (const auto& hline : header_lines) {
                    std::fwrite(hline.data(), 1, hline.size(), rf.fp);
                }
                rf.header_written = true;
            }

            std::fwrite(sam_buf + rec.offset, 1, rec.len, rf.fp);
            written_reads++;
        }

        close_region_files();

        std::fprintf(stderr,
                     "[OMP] chr %s done, regions=%zu, written_reads=%ld\n",
                     c.name.c_str(), c.regions.size(), (long)written_reads);
    }

    return global_ok != 0;
}

// ---------------- main ----------------
int main(int argc, char** argv)
{
    if (argc < 4) {
        std::fprintf(stderr,
                     "Usage: %s <ref.fa> <in.sam> <out_dir>\n"
                     "Example:\n"
                     "  %s ref.fa input.sam out_regions\n",
                     argv[0], argv[0]);
        return 1;
    }

    std::string fasta_path = argv[1];
    std::string sam_path   = argv[2];
    std::string out_dir    = argv[3];

    // 创建输出目录（如果不存在）
    struct stat st;
    if (stat(out_dir.c_str(), &st) != 0) {
        if (mkdir(out_dir.c_str(), 0755) != 0) {
            std::fprintf(stderr, "Failed to create out_dir: %s (%s)\n",
                         out_dir.c_str(), std::strerror(errno));
            return 1;
        }
    } else {
        if (!S_ISDIR(st.st_mode)) {
            std::fprintf(stderr, "out_dir exists and is not a directory: %s\n",
                         out_dir.c_str());
            return 1;
        }
    }

    // 1. 解析 FASTA（只保留 chr1-22, chrX, chrY）
    std::vector<ChrInfo> chrs;
    std::unordered_map<std::string,int> chr_index;
    if (!parse_fasta(fasta_path, chrs, chr_index)) {
        return 1;
    }

    // 2. 整个 SAM 读入内存并解析
    char* sam_buf = nullptr;
    size_t sam_size = 0;
    std::vector<std::string> header_lines;
    std::vector<SamRecord>   records;

    double t_load0 = now_ms();
    if (!load_and_parse_sam(sam_path, sam_buf, sam_size,
                            chr_index, chrs,
                            header_lines, records)) {
        return 1;
    }
    double t_load1 = now_ms();
    std::fprintf(stderr, "Load & parse SAM time: %.3f ms\n", t_load1 - t_load0);

    // 3. 对每个 chr 划分 region
    double target_bytes = TARGET_REGION_MB * 1024.0 * 1024.0;
    std::fprintf(stderr, "Target region size: %.1f MB (%.0f bytes)\n",
                 TARGET_REGION_MB, target_bytes);

    double t_reg0 = now_ms();
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::fprintf(stderr, "Build regions for chr %s ...\n",
                     chrs[i].name.c_str());
        build_regions_for_chr(chrs[i], target_bytes);
    }
    double t_reg1 = now_ms();
    std::fprintf(stderr, "Region building time: %.3f ms\n", t_reg1 - t_reg0);

    // 统计总 region 数
    size_t total_regions = 0;
    for (size_t i = 0; i < chrs.size(); ++i)
        total_regions += chrs[i].regions.size();
    std::fprintf(stderr, "Total regions: %zu\n", total_regions);

    // 4. 构建 per-chr record 列表
    std::vector< std::vector<int> > chr_rec_indices;
    build_chr_record_indices(records, chrs.size(), chr_rec_indices);

    // 5. 并行按 chr 切分：直接从 sam_buf 写 region
    double t_split0 = now_ms();
    std::fprintf(stderr, "Using OpenMP to split SAM by chromosome from memory.\n");
    bool ok = split_sam_by_regions_from_memory(sam_buf, out_dir,
                                               chrs, header_lines,
                                               records, chr_rec_indices);
    double t_split1 = now_ms();
    std::fprintf(stderr, "Split time: %.3f ms\n", t_split1 - t_split0);

    if (sam_buf) std::free(sam_buf);

    return ok ? 0 : 1;
}

