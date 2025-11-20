#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cstdlib>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cstring>

static const int NUM_REGIONS = 384 * 16;             // 3072
static const size_t BUF_SIZE = 4ULL * 1024 * 1024;  // 4MB

struct ChromInfo {
    std::string name;
    uint64_t length = 0;
    uint64_t offset = 0; // 全局起始坐标（0-based）
};

struct RegionMeta {
    std::string chr;
    uint64_t start_pos = 0; // 1-based
    uint64_t end_pos   = 0; // 1-based, inclusive
};

// 尝试把 NOFILE 软限制设到 target_nofile（如果超出硬限制，会被 clamp）
void set_nofile_limit(rlim_t target_nofile) {
    struct rlimit rl;
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) {
        std::cerr << "[WARN] getrlimit failed: " << strerror(errno) << "\n";
        return;
    }

    std::cout << "[INFO] Current RLIMIT_NOFILE soft=" << rl.rlim_cur
              << " hard=" << rl.rlim_max << "\n";

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

// 粗糙判断是不是我们关心的 chr1-22,X,Y
bool is_target_chrom(const std::string &name) {
    if (name.rfind("chr", 0) == 0) {
        std::string tail = name.substr(3);
        if (tail == "X" || tail == "Y") return true;
        // 数字 1-22
        char *endptr = nullptr;
        long v = std::strtol(tail.c_str(), &endptr, 10);
        if (*endptr == '\0' && v >= 1 && v <= 22) return true;
    }
    return false;
}

// 读取 fasta，统计 chr1-22,X,Y 的长度，并为它们分配 offset
uint64_t load_reference_lengths(
    const std::string &fa_path,
    std::unordered_map<std::string, ChromInfo> &chroms,
    std::vector<ChromInfo> &chrom_order)
{
    std::ifstream fin(fa_path);
    if (!fin) {
        std::cerr << "[ERROR] Failed to open fasta: " << fa_path << "\n";
        std::exit(1);
    }

    std::string line;
    std::string current_chr;
    uint64_t total_len = 0;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // 解析 header：>chr1 AC:... -> name="chr1"
            std::string header = line.substr(1);

            // 去掉行尾空白
            while (!header.empty() && (header.back() == '\r' || header.back() == ' ' || header.back() == '\t'))
                header.pop_back();
            // 去掉行首空白
            size_t p = 0;
            while (p < header.size() && (header[p] == ' ' || header[p] == '\t'))
                ++p;
            header = header.substr(p);

            // 取第一个空格/Tab 之前的部分
            size_t sp = header.find_first_of(" \t");
            std::string name = (sp == std::string::npos) ? header : header.substr(0, sp);

            current_chr.clear();
            if (is_target_chrom(name)) {
                current_chr = name;
                if (chroms.find(name) == chroms.end()) {
                    ChromInfo ci;
                    ci.name = name;
                    ci.length = 0;
                    ci.offset = 0; // 后面统一赋值
                    chroms[name] = ci;
                    chrom_order.push_back(ci); // 先 push，长度后面再同步更新
                }
            }
        } else {
            if (!current_chr.empty()) {
                // 累加当前 chr 的长度（忽略换行）
                uint64_t add = 0;
                for (char c : line) {
                    if (c == '\n' || c == '\r') continue;
                    ++add;
                }
                chroms[current_chr].length += add;
            }
        }
    }
    fin.close();

    // 重新按 chrom_order 同步长度，并计算 offset
    uint64_t offset = 0;
    for (auto &ci : chrom_order) {
        auto it = chroms.find(ci.name);
        if (it == chroms.end()) continue;
        it->second.offset = offset;
        ci.length = it->second.length;
        ci.offset = offset;
        offset += it->second.length;
    }
    total_len = offset;

    std::cout << "[INFO] Reference total length (chr1-22,X,Y) = " << total_len << "\n";
    for (auto &ci : chrom_order) {
        std::cout << "  " << ci.name << ": len=" << ci.length
                  << " offset=" << ci.offset << "\n";
    }

    return total_len;
}

// 把全局 0-based 坐标转成 (chr_index, 1-based pos)
bool global_to_chr_pos(const std::vector<ChromInfo> &chrom_order,
                       uint64_t global_pos,
                       int &chr_index,
                       uint64_t &chr_pos_1based)
{
    // chrom 数量很少，直接线性扫描就行
    for (size_t i = 0; i < chrom_order.size(); ++i) {
        const auto &ci = chrom_order[i];
        if (global_pos >= ci.offset && global_pos < ci.offset + ci.length) {
            chr_index = static_cast<int>(i);
            chr_pos_1based = (global_pos - ci.offset) + 1;
            return true;
        }
    }
    return false;
}

// 根据 (chrom, pos) 计算全局坐标，再映射到 region id
// pos 是 1-based SAM POS
int coord_to_region(const std::unordered_map<std::string, ChromInfo> &chroms,
                    const std::string &rname, int pos,
                    uint64_t region_size, int num_regions)
{
    auto it = chroms.find(rname);
    if (it == chroms.end()) return -1;
    if (pos <= 0) return -1;

    uint64_t global_pos = it->second.offset + static_cast<uint64_t>(pos - 1);
    uint64_t rid = global_pos / region_size;
    if (rid >= static_cast<uint64_t>(num_regions)) {
        rid = num_regions - 1;
    }
    return static_cast<int>(rid);
}

// 简单创建目录（如果不存在）
void make_dir(const std::string &dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) == 0) {
        if (S_ISDIR(st.st_mode)) return;
        std::cerr << "[ERROR] " << dir << " exists and is not a directory.\n";
        std::exit(1);
    }
    if (mkdir(dir.c_str(), 0755) != 0) {
        std::cerr << "[ERROR] mkdir(" << dir << ") failed: " << strerror(errno) << "\n";
        std::exit(1);
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <ref.fa> <aln.sam> <out_dir>\n";
        return 1;
    }

    std::string fa_path = argv[1];
    std::string sam_path = argv[2];
    std::string out_dir  = argv[3];

    // 1) 尝试将 NOFILE 软限制设到 ~NUM_REGIONS（留点余量）
    set_nofile_limit(NUM_REGIONS + 128);

    // 2) 读取 reference，计算 chr1-22,X,Y 长度和全局 offset
    std::unordered_map<std::string, ChromInfo> chroms;
    std::vector<ChromInfo> chrom_order;
    uint64_t total_len = load_reference_lengths(fa_path, chroms, chrom_order);
    if (total_len == 0) {
        std::cerr << "[ERROR] total_len == 0, check fasta / chrom names.\n";
        return 1;
    }

    // 3) 计算 region 大小（向上取整）
    uint64_t region_size = (total_len + NUM_REGIONS - 1) / NUM_REGIONS;
    std::cout << "[INFO] NUM_REGIONS=" << NUM_REGIONS
              << " region_size=" << region_size << "\n";

    // 4) 为每个 region 预先算一个“代表区间”用于文件命名
    std::vector<RegionMeta> regions(NUM_REGIONS);
    for (int rid = 0; rid < NUM_REGIONS; ++rid) {
        uint64_t global_start = static_cast<uint64_t>(rid) * region_size;
        uint64_t global_end_excl = global_start + region_size;
        if (global_end_excl > total_len) global_end_excl = total_len;
        if (global_start >= total_len) { // 保险起见
            global_start = total_len - 1;
            global_end_excl = total_len;
        }
        uint64_t global_end = global_end_excl - 1; // 0-based

        int chr_idx_s = -1, chr_idx_e = -1;
        uint64_t pos_s = 0, pos_e = 0;
        if (!global_to_chr_pos(chrom_order, global_start, chr_idx_s, pos_s)) {
            // 不太可能发生，直接标记为未知
            regions[rid].chr = "unknown";
            regions[rid].start_pos = 0;
            regions[rid].end_pos   = 0;
            continue;
        }
        if (!global_to_chr_pos(chrom_order, global_end, chr_idx_e, pos_e)) {
            // 同上
            regions[rid].chr = "unknown";
            regions[rid].start_pos = 0;
            regions[rid].end_pos   = 0;
            continue;
        }

        if (chr_idx_s == chr_idx_e) {
            // region 完全落在同一条染色体上
            regions[rid].chr = chrom_order[chr_idx_s].name;
            regions[rid].start_pos = pos_s;
            regions[rid].end_pos   = pos_e;
        } else {
            // region 跨染色体：用起始染色体，end 写到这一条染色体的末尾
            const auto &ci = chrom_order[chr_idx_s];
            regions[rid].chr = ci.name;
            regions[rid].start_pos = pos_s;
            regions[rid].end_pos   = ci.length; // 到该 chr 末尾
        }
    }
	
// === 新增：把所有 region 信息写到 region_info.txt ===
    {
        std::string info_path = "./region_info_6k.txt";
        std::ofstream info_out(info_path);
        if (!info_out) {
            std::cerr << "[ERROR] Failed to open " << info_path << " for write.\n";
            return 1;
        }
        for (int i = 0; i < NUM_REGIONS; ++i) {
            const auto &rm = regions[i];
            // 一行三个 word：chr start end
            info_out << rm.chr << " " << rm.start_pos << " " << rm.end_pos << "\n";
        }
        info_out.close();
        std::cout << "[INFO] Wrote region info to " << info_path << "\n";
    }

    // 5) 准备输出目录及 per-region 文件和 buffer
    make_dir(out_dir);

    std::vector<std::ofstream> region_files(NUM_REGIONS);
    std::vector<std::string>  buffers(NUM_REGIONS);
    for (int i = 0; i < NUM_REGIONS; ++i) {
        const auto &rm = regions[i];
        // 文件名格式：chr5_1_10000_88.sam
        // 若 chr == "unknown"，也照样写出来，方便调试
        std::string fname = out_dir + "/" + rm.chr + "_" +
                            std::to_string(rm.start_pos) + "_" +
                            std::to_string(rm.end_pos) + "_" +
                            std::to_string(i) + ".sam";

        region_files[i].open(fname, std::ios::out | std::ios::binary);
        if (!region_files[i]) {
            std::cerr << "[ERROR] Failed to open " << fname << " for write.\n";
            return 1;
        }
        buffers[i].reserve(BUF_SIZE);
    }
    std::cout << "[INFO] Opened " << NUM_REGIONS << " region files.\n";

    // 6) 读取 SAM，按 region 分发
    std::ifstream sam_in(sam_path);
    if (!sam_in) {
        std::cerr << "[ERROR] Failed to open SAM: " << sam_path << "\n";
        return 1;
    }

    std::string line;
    uint64_t total_reads = 0;
    uint64_t mapped_reads = 0;
    uint64_t unmapped_reads = 0;

    while (std::getline(sam_in, line)) {
        if (line.empty()) continue;

        // SAM header：这里直接跳过
        if (line[0] == '@') {
            continue;
        }

        ++total_reads;

        // 解析 RNAME 和 POS
        // col1: QNAME
        // col2: FLAG
        // col3: RNAME
        // col4: POS
        size_t t1 = line.find('\t');
        if (t1 == std::string::npos) { ++unmapped_reads; continue; }
        size_t t2 = line.find('\t', t1 + 1);
        if (t2 == std::string::npos) { ++unmapped_reads; continue; }
        size_t t3 = line.find('\t', t2 + 1);
        if (t3 == std::string::npos) { ++unmapped_reads; continue; }
        size_t t4 = line.find('\t', t3 + 1);
        if (t4 == std::string::npos) t4 = line.size();

        std::string rname = line.substr(t2 + 1, t3 - (t2 + 1));
        std::string pos_str = line.substr(t3 + 1, t4 - (t3 + 1));
        int pos = std::atoi(pos_str.c_str());
        if (rname == "*" || pos <= 0) {
            ++unmapped_reads;
            continue;
        }

        int rid = coord_to_region(chroms, rname, pos, region_size, NUM_REGIONS);
        if (rid < 0) {
            ++unmapped_reads;
            continue;
        }
        ++mapped_reads;

        // 往对应 region buffer 写一行（加 '\n'）
        std::string &buf = buffers[rid];
        size_t need = line.size() + 1; // +1 for '\n'
        if (buf.size() + need > BUF_SIZE) {
            // flush
            region_files[rid].write(buf.data(), buf.size());
            buf.clear();
        }
        buf.append(line);
        buf.push_back('\n');
    }

    sam_in.close();

    // flush 所有剩余 buffer
    for (int i = 0; i < NUM_REGIONS; ++i) {
        if (!buffers[i].empty()) {
            region_files[i].write(buffers[i].data(), buffers[i].size());
        }
        region_files[i].close();
    }

    std::cout << "[INFO] Done.\n";
    std::cout << "  total_reads    = " << total_reads << "\n";
    std::cout << "  mapped_reads   = " << mapped_reads << "\n";
    std::cout << "  unmapped_reads = " << unmapped_reads << "\n";

    return 0;
}

