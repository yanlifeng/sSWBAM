# sSWBAM

A simple parallel SAM sorting and duplicate marking tool for the Sunway platform.

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

2. **编译 Sunway SAM 处理工具**
   ```bash
   make
   ```
   - 生成可执行文件 `sw_sam_process`
   - 支持排序、去重、以及排序+去重的完整流程

3. **运行 Sunway 并行处理（排序 + 去重）**
   ```bash
   # 完整流程：排序 + 标记重复（推荐）
   ./sw_sam_process --all <out_regions_sam> <out_regions_processed>
   
   # 或者分步执行：
   # 步骤 1：仅排序
   ./sw_sam_process --sort <out_regions_sam> <out_regions_sorted>
   
   # 步骤 2：仅标记重复（需要输入已排序的文件）
   ./sw_sam_process --markdup <out_regions_sorted> <out_regions_marked>
   ```
   - 使用 Sunway 从核（CPE）并行处理，每个 CPE 处理一个文件
   - `--all`: 先按 RNAME（染色体）+ POS（位置）排序，再标记重复序列
   - `--sort`: 仅排序
   - `--markdup`: 仅标记重复（输入必须已排序）
   - 输出处理后的文件到指定目录
   - **输出目录**：如果不存在会自动创建，如果存在会清空后使用
   
   **输出文件命名规则**：
   - `--sort` 模式：`input.sam` → `input.sorted.sam`
   - `--markdup` 模式：`input.sam` → `input.markdup.sam`
   - `--all` 模式：`input.sam` → `input.sorted.markdup.sam`

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

### Sunway 处理工具
```bash
make
```

## 处理模式说明

### `--all` 模式（推荐）
- 完整流程：先排序再标记重复
- 适用于需要完整处理的场景
- 在从核内部完成排序和去重，效率最高

### `--sort` 模式
- 仅对 SAM 文件按 RNAME（染色体）+ POS（位置）排序
- 适用于只需要排序的场景

### `--markdup` 模式
- 仅标记重复序列（PCR duplicates）
- **要求输入文件必须已经排序**
- 标记重复后，FLAG 字段会添加 1024（0x400）标志

## 注意事项

- x86 预处理步骤需要足够内存来加载整个 SAM 文件
- Sunway 处理工具使用 64 个 CPE 并行处理，每个 CPE 处理一个 SAM 文件
- 单个 SAM 文件大小不应超过 100MB（可在 `src/main.c` 中调整 `MAX_BUF_SIZE`）
- `--all` 模式会在从核内部完成排序和去重，无需中间文件
- 去重算法基于位置和质量分数，保留质量最高的序列
