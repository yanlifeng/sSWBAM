// split_from_region.cpp
//
// 用法：
//   g++ -O3 -std=gnu++11 split_from_region.cpp -o split_from_region
//   ./split_from_region region.txt all.sam out_regions_sam
//
// 功能：
//   1. 从 region.txt 读取若干 region：每行格式为
//        chr  start  end
//      例如：
//        chr1  1  1000000
//        chr1  1000001  2000000
//      支持空行和以 '#' 开头的注释行。
//   2. 要求 region 总数 < 3000，否则退出。
//   3. 为每个 region 分配一个 2MB 的内存 buffer。
//   4. 扫描 all.sam：
//        - 收集 header 行（以 '@' 开头）；
//        - 解析每条对齐记录的 RNAME 和 POS；
//        - 根据 (chr, pos) 找到所属 region，将该行放入对应 buffer；
//        - buffer 满了则 flush 到文件，然后继续装；
//      输出文件命名为： out_dir/chr_start_end.sam
//      每个 region 文件在第一次写入前会先写入完整 SAM header。

#define _GNU_SOURCE
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <errno.h>
#include <algorithm>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cstring>
#include <algorithm>
#include <sys/time.h>
#include <iostream>


static double now_ms()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

static const size_t REGION_BUF_SIZE = 512u * 1024u;
static const size_t MAX_REGION_NUM  = 3000;

// ------------- 简单结构体 -------------

struct Region {
    std::string chr;
    long long   start;   // 1-based
    long long   end;     // inclusive

    std::string out_path;

    std::vector<char> buffer;
    size_t            used;
    bool              header_written;

    Region() : start(0), end(0), used(0), header_written(false) {}
};

void set_nofile_limit(rlim_t target_nofile) {
    struct rlimit rl;
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) {
        std::cerr << "[WARN] getrlimit failed: " << strerror(errno) << "\n";
        return;
    }

    std::cout << "[INFO] Current RLIMIT_NOFILE soft=" << rl.rlim_cur
              << " hard=" << rl.rlim_max << "\n";

    // 如果当前 soft 已经够大，就不用改了，直接返回
    if (target_nofile <= rl.rlim_cur) {
        std::cout << "[INFO] target_nofile (" << target_nofile
                  << ") <= current soft limit, no change.\n";
        return;
    }

    // 想调得比 hard 还大，就 clamp 到 hard
    if (target_nofile > rl.rlim_max) {
        std::cerr << "[WARN] target_nofile(" << target_nofile
                  << ") > hard limit(" << rl.rlim_max
                  << "), will set to hard limit.\n";
        target_nofile = rl.rlim_max;
    }

    rl.rlim_cur = target_nofile;
    if (setrlimit(RLIMIT_NOFILE, &rl) != 0) {
        std::cerr << "[WARN] setrlimit failed: " << strerror(errno) << "\n";
    } else {
        std::cout << "[INFO] RLIMIT_NOFILE soft set to " << target_nofile << "\n";
    }
}



// ------------- 工具函数：字符串 trim -------------

static void trim(std::string &s) {
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n'))
        ++i;
    size_t j = s.size();
    while (j > i && (s[j-1] == ' ' || s[j-1] == '\t' || s[j-1] == '\r' || s[j-1] == '\n'))
        --j;
    if (i == 0 && j == s.size()) return;
    s.assign(s.begin() + i, s.begin() + j);
}

// ------------- 解析 region.txt -------------

static bool load_regions(const std::string &region_file,
                         const std::string &out_dir,
                         std::vector<Region> &regions)
{
    FILE *fp = std::fopen(region_file.c_str(), "rb");
    if (!fp) {
        std::fprintf(stderr, "Failed to open region file: %s (%s)\n",
                     region_file.c_str(), std::strerror(errno));
        return false;
    }

    char *line = nullptr;
    size_t cap = 0;
    long long line_no = 0;

    while (true) {
        ssize_t n = getline(&line, &cap, fp);
        if (n < 0) break;
        ++line_no;

        std::string s(line, (size_t)n);
        trim(s);
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        // 拆分为 tokens
        std::vector<std::string> tokens;
        size_t i = 0;
        while (i < s.size()) {
            while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
            if (i >= s.size()) break;
            size_t j = i;
            while (j < s.size() && s[j] != ' ' && s[j] != '\t') ++j;
            tokens.emplace_back(s.substr(i, j - i));
            i = j;
        }
        if (tokens.size() < 3) {
            std::fprintf(stderr,
                         "Bad line in region file (need at least 3 columns): line %lld: %s\n",
                         line_no, s.c_str());
            std::fclose(fp);
            if (line) std::free(line);
            return false;
        }

        Region r;
        r.chr = tokens[0];
        char *endp = nullptr;
        r.start = std::strtoll(tokens[1].c_str(), &endp, 10);
        if (*endp != '\0') {
            std::fprintf(stderr, "Bad start in line %lld: %s\n",
                         line_no, tokens[1].c_str());
            std::fclose(fp);
            if (line) std::free(line);
            return false;
        }
        r.end = std::strtoll(tokens[2].c_str(), &endp, 10);
        if (*endp != '\0') {
            std::fprintf(stderr, "Bad end in line %lld: %s\n",
                         line_no, tokens[2].c_str());
            std::fclose(fp);
            if (line) std::free(line);
            return false;
        }
        if (r.start <= 0 || r.start > r.end) {
            std::fprintf(stderr, "Invalid region [%lld,%lld] at line %lld\n",
                         r.start, r.end, line_no);
            std::fclose(fp);
            if (line) std::free(line);
            return false;
        }

        // 构造输出文件名： out_dir/chr_start_end.sam
        char path[4096];
        std::snprintf(path, sizeof(path), "%s/%s_%lld_%lld.sam",
                      out_dir.c_str(), r.chr.c_str(), r.start, r.end);
        r.out_path = path;

        r.buffer.resize(REGION_BUF_SIZE);
        r.used = 0;
        r.header_written = false;

        regions.push_back(r);
        if (regions.size() >= MAX_REGION_NUM) {
            std::fprintf(stderr,
                         "Too many regions (>= %zu). Adjust MAX_REGION_NUM or region.txt\n",
                         (size_t)MAX_REGION_NUM);
            std::fclose(fp);
            if (line) std::free(line);
            return false;
        }
    }

    std::fclose(fp);
    if (line) std::free(line);

    std::fprintf(stderr, "Loaded %zu regions from %s\n",
                 regions.size(), region_file.c_str());
    return true;
}

// ------------- SAM RNAME+POS 解析 -------------

static bool parse_sam_rname_pos(const char* line,
                                size_t len,       // 不含换行
                                const char*& rname_start,
                                size_t& rname_len,
                                long long& pos_out)
{
    if (len == 0) return false;
    if (line[0] == '@') return false; // header

    int    field = 0;
    size_t i     = 0;
    size_t start = 0;
    size_t end   = 0;

    const char* rname_s = nullptr;
    size_t      rname_l = 0;
    const char* pos_s   = nullptr;
    size_t      pos_l   = 0;

    for (i = 0; i <= len; ++i) {
        char c = (i == len) ? '\t' : line[i];

        if (c == '\t') {
            end = i;
            if (field == 2) {
                rname_s = line + start;
                rname_l = end - start;
            } else if (field == 3) {
                pos_s = line + start;
                pos_l = end - start;
                break;   // 后面字段不关心
            }
            field++;
            start = i + 1;
        }
    }

    if (!rname_s || !pos_s || pos_l == 0) return false;

    // 解析 POS
    long long value = 0;
    int       neg   = 0;
    size_t    idx   = 0;
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
    rname_len   = rname_l;
    pos_out     = value;
    return true;
}

// ------------- 用于 region 查找的 per-chr 索引 -------------

// 构造： chr -> 对应 region 下标列表（按 start 升序）
static void build_chr_region_index(const std::vector<Region> &regions,
                                   std::unordered_map<std::string, std::vector<int> > &chr2regs)
{
    for (size_t i = 0; i < regions.size(); ++i) {
        const Region &r = regions[i];
        chr2regs[r.chr].push_back((int)i);
    }

    // 每个 chr 对应 region 已在 region.txt 中按顺序给出的话，可以不排序；
    // 为稳妥起见这里还是按 start 排一下。
    for (auto &kv : chr2regs) {
        std::vector<int> &idxs = kv.second;
        std::sort(idxs.begin(), idxs.end(), [&](int a, int b) {
            const Region &ra = regions[a];
            const Region &rb = regions[b];
            if (ra.start != rb.start) return ra.start < rb.start;
            return ra.end < rb.end;
        });
    }
}

// 根据 chr + pos 找到对应 region 的下标（在 regions 向量中的 index）。
// 返回 -1 表示没找到。
static int find_region_for_pos(const std::string &chr,
                               long long pos,
                               const std::vector<Region> &regions,
                               const std::unordered_map<std::string, std::vector<int> > &chr2regs)
{
    auto it = chr2regs.find(chr);
    if (it == chr2regs.end()) return -1;

    const std::vector<int> &idxs = it->second;
    int left = 0;
    int right = (int)idxs.size() - 1;
    int found = -1;

    while (left <= right) {
        int mid = (left + right) / 2;
        const Region &r = regions[idxs[mid]];
        if (pos < r.start) {
            right = mid - 1;
        } else if (pos > r.end) {
            left = mid + 1;
        } else {
            found = idxs[mid];
            break;
        }
    }
    return found;
}

// ------------- 把一个 region 的 buffer flush 到文件 -------------

static bool flush_region_buffer(Region &r,
                                const std::vector<std::string> &header_lines)
{
    if (r.used == 0) return true; // 没数据就不用写

    FILE *fp = std::fopen(r.out_path.c_str(), r.header_written ? "ab" : "wb");
    if (!fp) {
        std::fprintf(stderr,
                     "Failed to open region file for writing: %s (%s)\n",
                     r.out_path.c_str(), std::strerror(errno));
        return false;
    }

    // 第一次写入时需要 header
    if (!r.header_written) {
        for (const auto &h : header_lines) {
            std::fwrite(h.data(), 1, h.size(), fp);
        }
        r.header_written = true;
    }

    std::fwrite(r.buffer.data(), 1, r.used, fp);
    std::fclose(fp);

    r.used = 0;
    return true;
}

// ------------- 主逻辑：按 region.txt split SAM -------------

static bool split_by_regions(const std::string &sam_path,
                             std::vector<Region> &regions,
                             const std::string &out_dir)
{
    // 先构造 chr -> region 下标列表
    std::unordered_map<std::string, std::vector<int> > chr2regs;
    build_chr_region_index(regions, chr2regs);

    FILE *fp = std::fopen(sam_path.c_str(), "rb");
    if (!fp) {
        std::fprintf(stderr, "Failed to open SAM: %s (%s)\n",
                     sam_path.c_str(), std::strerror(errno));
        return false;
    }

    char  *line = nullptr;
    size_t cap  = 0;

    std::vector<std::string> header_lines;

    long long total_records   = 0;
    long long assigned_records = 0;

    while (true) {
        ssize_t n = getline(&line, &cap, fp);
        if (n < 0) break;
        if (n == 0) continue;

        size_t line_len = (size_t)n;
        size_t text_len = line_len;
        while (text_len > 0 &&
              (line[text_len-1] == '\n' || line[text_len-1] == '\r')) {
            --text_len;
        }

        if (text_len == 0) continue;

        if (line[0] == '@') {
            // header 直接缓存起来，后面每个 region 文件写一次
            header_lines.emplace_back(line, line_len);
            continue;
        }

        total_records++;

        const char *rname_s = nullptr;
        size_t      rname_l = 0;
        long long   pos     = 0;
        if (!parse_sam_rname_pos(line, text_len, rname_s, rname_l, pos)) {
            // 解析失败就丢掉
            continue;
        }

        if (pos <= 0) {
            // unmapped 或非法 POS，暂时不写
            continue;
        }

        std::string chr(rname_s, rname_l);
        int ridx = find_region_for_pos(chr, pos, regions, chr2regs);
        if (ridx < 0) {
            // 说明这个 chr/pos 不在任何 region 里，略过
            continue;
        }

        Region &r = regions[ridx];

        // 如果 buffer 放不下这一行，先 flush
        if (r.used + line_len > r.buffer.size()) {
            if (!flush_region_buffer(r, header_lines)) {
                std::fprintf(stderr,
                             "Flush failed for region file: %s\n",
                             r.out_path.c_str());
                std::fclose(fp);
                if (line) std::free(line);
                return false;
            }
        }

        // 再把这一行放进去
        if (line_len > r.buffer.size()) {
            // 特殊情况：单行比 buffer 还大，直接写一次到文件中（附加 header）
            if (!r.header_written) {
                // 写 header + 这行
                FILE *fp2 = std::fopen(r.out_path.c_str(), "wb");
                if (!fp2) {
                    std::fprintf(stderr,
                                 "Failed to open region file (large line): %s (%s)\n",
                                 r.out_path.c_str(), std::strerror(errno));
                    std::fclose(fp);
                    if (line) std::free(line);
                    return false;
                }
                for (const auto &h : header_lines) {
                    std::fwrite(h.data(), 1, h.size(), fp2);
                }
                std::fwrite(line, 1, line_len, fp2);
                std::fclose(fp2);
                r.header_written = true;
            } else {
                FILE *fp2 = std::fopen(r.out_path.c_str(), "ab");
                if (!fp2) {
                    std::fprintf(stderr,
                                 "Failed to open region file (large line, append): %s (%s)\n",
                                 r.out_path.c_str(), std::strerror(errno));
                    std::fclose(fp);
                    if (line) std::free(line);
                    return false;
                }
                std::fwrite(line, 1, line_len, fp2);
                std::fclose(fp2);
            }
        } else {
            // 正常情况：复制到 buffer
            std::memcpy(r.buffer.data() + r.used, line, line_len);
            r.used += line_len;
        }

        assigned_records++;
    }

    std::fclose(fp);
    if (line) std::free(line);

    // 把每个 region 剩余的数据 flush 一次
    for (size_t i = 0; i < regions.size(); ++i) {
        Region &r = regions[i];
        if (r.used > 0) {
            if (!flush_region_buffer(r, header_lines)) {
                std::fprintf(stderr,
                             "Final flush failed for region file: %s\n",
                             r.out_path.c_str());
                return false;
            }
        }
    }

    std::fprintf(stderr,
                 "Split done. total_records=%lld, assigned_records=%lld\n",
                 total_records, assigned_records);
    return true;
}

// ------------- main -------------

int main(int argc, char **argv)
{
    //volatile int aa = 0;
    //fprintf(stderr, "aa %d\n", aa);
    //while(aa == 0) {
    //    
    //}
    //fprintf(stderr, "aa %d\n", aa);

    
    if (argc < 4) {
        std::fprintf(stderr,
                     "Usage: %s <region.txt> <all.sam> <out_dir>\n"
                     "Example:\n"
                     "  %s region.txt all.sam out_regions_sam\n",
                     argv[0], argv[0]);
        return 1;
    }

    std::string region_file = argv[1];
    std::string sam_file    = argv[2];
    std::string out_dir     = argv[3];

    // 创建输出目录（若不存在）
    struct stat st;
    if (stat(out_dir.c_str(), &st) != 0) {
        if (mkdir(out_dir.c_str(), 0755) != 0) {
            std::fprintf(stderr,
                         "Failed to create out_dir %s (%s)\n",
                         out_dir.c_str(), std::strerror(errno));
            return 1;
        }
    } else {
        if (!S_ISDIR(st.st_mode)) {
            std::fprintf(stderr,
                         "out_dir exists and is not a directory: %s\n",
                         out_dir.c_str());
            return 1;
        }
    }
    double t0 = now_ms();

    std::vector<Region> regions;
    if (!load_regions(region_file, out_dir, regions)) {
        return 1;
    }

    set_nofile_limit(regions.size() + 128);
    fprintf(stderr, "load_regions %.2f ms\n", now_ms() - t0);

    t0 = now_ms();
    if (regions.empty()) {
        std::fprintf(stderr, "No regions loaded from %s\n", region_file.c_str());
        return 1;
    }

    if (!split_by_regions(sam_file, regions, out_dir)) {
        return 1;
    }
    fprintf(stderr, "split_by_regions %.2f ms\n", now_ms() - t0);

    return 0;
}

