# sSWBAM

A simple parallel SAM sorting tool for the Sunway platform.

## 工作流程

### 第一步：在 x86 平台上预处理（需要大内存和 OpenMP 并行）

1. **自动划分区域并生成 SAM 文件**
   ```bash
   # 在 x86 机器上运行 auto_region
   ./pre-tools/auto_region <input.sam> <output_dir>
   ```
   - 这一步会自动分析 SAM 文件，按染色体和位置划分区域
   - 生成划分好的 SAM 文件到 `output_dir`
   - 需要较大内存，使用 OpenMP 并行加速

2. **检查划分结果并生成区域配置文件**
   ```bash
   # 检查生成的 SAM 文件并生成 region_auto.txt
   ./pre-tools/check_sam <output_dir>
   ```
   - 验证划分的正确性
   - 生成 `region_auto.txt` 配置文件，记录每个区域的信息

### 第二步：在 Sunway 平台上排序

1. **根据区域配置文件划分 SAM**
   ```bash
   # 在 Sunway 上运行 split_from_region
   ./pre-tools/split_from_region <region_auto.txt> <input.sam> <out_regions_sam>
   ```
   - 根据 `region_auto.txt` 将输入 SAM 文件划分到 `out_regions_sam` 目录

2. **编译 Sunway SAM 排序工具**
   ```bash
   make
   ```
   - 生成可执行文件 `sw_sam_sort`

3. **运行 Sunway 并行排序**
   ```bash
   ./sw_sam_sort <out_regions_sam> <out_regions_sorted_sam>
   ```
   - 使用 Sunway 从核（CPE）并行排序
   - 每个 SAM 文件按 RNAME（染色体）和 POS（位置）排序
   - 输出排序后的文件到 `out_regions_sorted_sam` 目录

## 目录结构

```
sSWBAM/
├── pre-tools/          # x86 预处理工具
│   ├── auto_region.cpp      # 自动区域划分
│   ├── check_sam.cpp        # 检查并生成配置
│   └── split_from_region.cpp # 根据配置划分
├── src/                # Sunway 主核代码
│   └── main.c               # 主核入口和文件 I/O
├── slave/              # Sunway 从核代码
│   ├── sam_sort_para.h      # 参数结构定义
│   └── slave.c              # 从核排序逻辑
└── Makefile            # Sunway 编译配置
```

## 编译说明

### x86 预处理工具
```bash
cd pre-tools
g++ -O3 -fopenmp auto_region.cpp -o auto_region
g++ -O3 check_sam.cpp -o check_sam
g++ -O3 split_from_region.cpp -o split_from_region
```

### Sunway 排序工具
```bash
make
```

## 注意事项

- x86 预处理步骤需要足够内存来加载整个 SAM 文件
- Sunway 排序工具使用 64 个 CPE 并行处理，每个 CPE 处理一个 SAM 文件
- 单个 SAM 文件大小不应超过 100MB（可在 `src/main.c` 中调整 `MAX_BUF_SIZE`）
- 输出文件名格式：`<input_filename>.sorted.sw.sam`
