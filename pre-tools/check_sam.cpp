// check_sam.cpp
// 用法：
//   g++ -O3 -std=gnu++11 check_sam.cpp -o check_sam
//   ./check_sam sam_dir
//
// 功能：
//   1) 对 sam_dir 中的每个 SAM 文件：
//        - 从文件名解析 chrName, start, end：
//            例：chr13_80350001_87900000.sam
//                chr10_42163648_43168944_1708.sam.sorted.sw.sam
//            取第一个 ".sam" 前的部分，按 '_' 拆出前三段：
//                tokens[0]=chrName, tokens[1]=start, tokens[2]=end
//        - 扫描 SAM 内容，解析每条记录的 RNAME, POS；
//
//        - 只打印**有问题**的文件：
//             - RNAME != chrName，或
//             - POS 不在 [start,end]
//          正常通过的文件不打印。
//
//   2) 扫描结束后，在 sam_dir 里生成 region_auto.txt：
//        每行：chr  start  end
//        只包含能从文件名解析出的 region（自动去重）。
//
// ----------------------------------------------------------

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_set>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

// 简单辅助：判断字符串是否以 suffix 结尾（目前没用到，可留着）
static bool ends_with(const std::string& s, const std::string& suf) {
    if (s.size() < suf.size()) return false;
    return std::memcmp(s.data() + s.size() - suf.size(), suf.data(), suf.size()) == 0;
}

// 从文件名解析 chrName, start, end
// 规则：
//   1. 取第一个 ".sam" 之前的部分 base：例如
//        "chr13_80350001_87900000.sam"                   → base="chr13_80350001_87900000"
//        "chr10_42163648_43168944_1708.sam.sorted.sw.sam"→ base="chr10_42163648_43168944_1708"
//   2. 以 '_' 分割 base，要求至少 3 段，
//      tokens[0] = chrName, tokens[1] = start, tokens[2] = end
static bool parse_filename_region(const std::string& fname,
                                  std::string& chrName,
                                  long long& start,
                                  long long& end)
{
    // 取第一个 ".sam" 的位置
    size_t pos = fname.find(".sam");
    if (pos == std::string::npos) {
        return false;
    }
    std::string base = fname.substr(0, pos);

    // 按 '_' 分割
    std::vector<std::string> tokens;
    size_t i = 0;
    while (i < base.size()) {
        size_t j = i;
        while (j < base.size() && base[j] != '_') ++j;
        tokens.push_back(base.substr(i, j - i));
        if (j == base.size()) break;
        i = j + 1;
    }
    if (tokens.size() < 3) {
        return false;
    }

    chrName = tokens[0];
    char* endp = nullptr;

    start = std::strtoll(tokens[1].c_str(), &endp, 10);
    if (*endp != '\0') return false;

    end = std::strtoll(tokens[2].c_str(), &endp, 10);
    if (*endp != '\0') return false;

    if (start > end || start <= 0) return false;
    return true;
}

// 解析 SAM 一行的 RNAME + POS
// 输入: line(不一定以 '\0' 结束)，len(不含 '\n')。
// 输出: rname_start/rname_len, pos_out
// 返回: true=解析成功且不是 header；false=失败或 header。
static bool parse_sam_rname_pos(const char* line,
                                size_t len,
                                const char*& rname_start,
                                size_t& rname_len,
                                long long& pos_out)
{
    if (len == 0) return false;
    if (line[0] == '@') return false; // header 行

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
                break;  // 后面字段不关心
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

// 检查单个 SAM 文件
// 返回：true = 该文件通过（没有 bad_chr/bad_range）；false = 有问题，需要在主调那边计数。
static bool check_one_sam(const std::string& path, const std::string& fname)
{
    std::string chrName;
    long long region_start = 0;
    long long region_end   = 0;

    if (!parse_filename_region(fname, chrName, region_start, region_end)) {
        // 无法从文件名解析 region，就跳过检查（不判为“错误文件”）
        // 如果你想在这里报警，可以解开注释：
        // std::fprintf(stderr,
        //              "[SKIP] %s: cannot parse chr/start/end from filename\n",
        //              fname.c_str());
        return true;
    }

    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) {
        std::fprintf(stderr,
                     "[ERROR] Failed to open %s (%s)\n",
                     path.c_str(), std::strerror(errno));
        return false;
    }

    char*  line = nullptr;
    size_t cap  = 0;

    long long total_records    = 0;
    long long checked_records  = 0;
    long long bad_chr          = 0;
    long long bad_range        = 0;
    long long unmapped_or_zero = 0;

    const int MAX_PRINT_ERR = 10;
    int printed_err = 0;

    while (true) {
        ssize_t n = getline(&line, &cap, fp);
        if (n < 0) break;
        if (n == 0) continue;

        size_t len      = (size_t)n;
        size_t text_len = len;
        while (text_len > 0 &&
              (line[text_len-1] == '\n' || line[text_len-1] == '\r')) {
            --text_len;
        }
        if (text_len == 0) continue;

        if (line[0] == '@') {
            continue; // header
        }

        total_records++;

        const char* rname_s = nullptr;
        size_t      rname_l = 0;
        long long   pos     = 0;
        if (!parse_sam_rname_pos(line, text_len, rname_s, rname_l, pos)) {
            // 解析失败，记到 unmapped_or_zero
            unmapped_or_zero++;
            if (printed_err < MAX_PRINT_ERR) {
                std::fprintf(stderr,
                             "  [WARN] %s: failed to parse RNAME/POS, line: %.80s\n",
                             fname.c_str(), line);
                printed_err++;
            }
            continue;
        }

        std::string rname(rname_s, rname_l);

        if (rname == "*" || pos <= 0) {
            unmapped_or_zero++;
            continue;
        }

        checked_records++;

        if (rname != chrName) {
            bad_chr++;
            if (printed_err < MAX_PRINT_ERR) {
                std::fprintf(stderr,
                             "  [ERR-CHR] %s: RNAME=%s POS=%lld (expect chr=%s [%lld,%lld])\n",
                             fname.c_str(), rname.c_str(), pos,
                             chrName.c_str(), region_start, region_end);
                printed_err++;
            }
            continue;
        }

        if (pos < region_start || pos > region_end) {
            bad_range++;
            if (printed_err < MAX_PRINT_ERR) {
                std::fprintf(stderr,
                             "  [ERR-RANGE] %s: RNAME=%s POS=%lld not in [%lld,%lld]\n",
                             fname.c_str(), rname.c_str(), pos,
                             region_start, region_end);
                printed_err++;
            }
        }
    }

    if (line) std::free(line);
    std::fclose(fp);

    // 如果没有 chr/范围错误，就认为该文件通过，啥也不打印
    if (bad_chr == 0 && bad_range == 0) {
        return true;
    }

    // 有问题的文件，打印一份总览
    std::fprintf(stderr,
                 "[FAIL] %s (chr=%s [%lld,%lld]):\n"
                 "  total_records     = %lld\n"
                 "  checked_records   = %lld (mapped, parsed OK)\n"
                 "  bad_chr           = %lld\n"
                 "  bad_range         = %lld\n"
                 "  unmapped_or_zero  = %lld\n",
                 fname.c_str(),
                 chrName.c_str(), region_start, region_end,
                 total_records,
                 checked_records,
                 bad_chr,
                 bad_range,
                 unmapped_or_zero);

    return false;
}

// 用于 region_auto.txt 的结构体
struct RegionInfo {
    std::string chr;
    long long   start;
    long long   end;
};

// ---------------- main ----------------

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr,
                     "Usage: %s <sam_dir>\n"
                     "Example:\n"
                     "  %s out_regions\n",
                     argv[0], argv[0]);
        return 1;
    }

    std::string dir_path = argv[1];

    DIR* dp = opendir(dir_path.c_str());
    if (!dp) {
        std::fprintf(stderr,
                     "Failed to open directory %s (%s)\n",
                     dir_path.c_str(), std::strerror(errno));
        return 1;
    }

    struct dirent* ent;
    int file_count   = 0;
    int fail_count   = 0;

    // 为生成 region_auto.txt 准备：去重用 set，顺序用 vector
    std::vector<RegionInfo> region_list;
    std::unordered_set<std::string> region_set; // key = "chr\tstart\tend"

    while ((ent = readdir(dp)) != nullptr) {
        const char* name = ent->d_name;
        if (std::strcmp(name, ".") == 0 || std::strcmp(name, "..") == 0)
            continue;

        std::string fname(name);
        // 只处理包含 ".sam" 的文件
        if (fname.find(".sam") == std::string::npos) continue;

        std::string full_path = dir_path;
        if (!full_path.empty() && full_path.back() != '/')
            full_path.push_back('/');
        full_path += fname;

        struct stat st;
        if (stat(full_path.c_str(), &st) != 0) {
            std::fprintf(stderr,
                         "  [SKIP] %s: stat failed (%s)\n",
                         full_path.c_str(), std::strerror(errno));
            continue;
        }
        if (!S_ISREG(st.st_mode)) {
            continue;
        }

        // 先尝试从文件名解析 region，用于 region_auto.txt
        std::string chrName;
        long long rstart = 0, rend = 0;
        if (parse_filename_region(fname, chrName, rstart, rend)) {
            char keybuf[256];
            std::snprintf(keybuf, sizeof(keybuf),
                          "%s\t%lld\t%lld",
                          chrName.c_str(), rstart, rend);
            std::string key(keybuf);
            if (region_set.find(key) == region_set.end()) {
                region_set.insert(key);
                RegionInfo info;
                info.chr   = chrName;
                info.start = rstart;
                info.end   = rend;
                region_list.push_back(info);
            }
        }

        file_count++;
        bool ok = check_one_sam(full_path, fname);
        if (!ok) {
            fail_count++;
        }
    }

    closedir(dp);

    // 生成 region_auto.txt
    std::string region_auto_path = "region_auto.txt";

    FILE* fp_out = std::fopen(region_auto_path.c_str(), "wb");
    if (!fp_out) {
        std::fprintf(stderr,
                     "Failed to write region_auto.txt at %s (%s)\n",
                     region_auto_path.c_str(), std::strerror(errno));
        // 不影响检查结果，只是提示一下
    } else {
        for (size_t i = 0; i < region_list.size(); ++i) {
            const RegionInfo& r = region_list[i];
            std::fprintf(fp_out, "%s\t%lld\t%lld\n",
                         r.chr.c_str(), r.start, r.end);
        }
        std::fclose(fp_out);
        std::fprintf(stderr,
                     "region_auto.txt written: %s (regions=%zu)\n",
                     region_auto_path.c_str(), region_list.size());
    }

    std::fprintf(stderr,
                 "Checked %d SAM files in directory: %s\n"
                 "  Failed files: %d\n",
                 file_count, dir_path.c_str(), fail_count);
    return 0;
}

