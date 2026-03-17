# Easy Immunofluorescence Analysis / 简易免疫荧光分析工具

Batch immunofluorescence (IF) image analysis tool with automatic channel detection, colocalization analysis, bleed-through correction, and statistical comparison.

批量免疫荧光图像分析工具，支持自动通道检测、共定位分析、串色校正和组间统计比较。

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/github/v/release/zhumiao-cloud/easy-immunofluorance-analysis)](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/releases)

---

## Features / 功能特性

| Feature / 功能 | Description / 说明 |
|----------------|-------------------|
| **Automatic Channel Detection** / 自动通道检测 | Auto-detect DAPI/488/594/647 channels from metadata, filenames, or heuristics<br>从 OME-TIFF 元数据、文件名或启发式算法自动识别通道 |
| **Flexible Input Modes** / 灵活输入模式 | Multi-channel single files or sample folders with single-channel files<br>支持多通道单文件，或单通道文件组成的样本文件夹 |
| **Nuclei Segmentation** / 细胞核分割 | Watershed-based segmentation with configurable parameters<br>基于 Watershed 的分割，参数可配置 |
| **Intensity Analysis** / 强度定量 | Background-subtracted signal intensity and positive area quantification<br>背景扣除后的信号强度和阳性面积分数 |
| **Colocalization Analysis** / 共定位分析 | Pearson r, Manders' coefficients, overlap coefficient<br>Pearson 相关系数、Manders' M1/M2、重叠系数 |
| **Bleed-through Correction** / 串色校正 | Automatic or manual crosstalk coefficient estimation<br>自动或手动估计串色系数，消除通道间串扰 |
| **Debris Filtering** / 碎点过滤 | Area-based filtering of small objects after thresholding<br>阈值后按面积过滤小物体，减少背景噪声 |
| **Statistical Comparison** / 组间统计 | Two-group testing (Welch t-test / Mann-Whitney U) with FDR correction<br>两组比较，支持 FDR 校正 |
| **Visualization** / 可视化输出 | Automated QC overlays, channel diagnostics, summary plots<br>自动生成 QC 叠加图、通道诊断图、汇总统计图 |

---

## Installation / 安装

### Requirements / 环境要求

- Python >= 3.8
- Recommended: Use conda environment / 推荐使用 conda 环境

### Dependencies / 安装依赖

```bash
pip install numpy pandas scipy scikit-image opencv-python-headless matplotlib tifffile pillow
```

Or use requirements.txt / 或使用 requirements.txt：
```bash
pip install -r requirements.txt
```

---

## Quick Start / 快速开始

### 1. Prepare Your Data / 准备数据

**Mode 1: Multi-channel single file** (Recommended, e.g., Zeiss/Thermo TIFF)
**模式一：多通道单文件**（推荐，如 Zeiss/Thermo 导出的 TIFF）

```
data/
├── group1/              # 实验组1
│   ├── image1.tif       # 4-channel TIFF: DAPI/488/594/647
│   └── image2.tif
└── group2/              # 实验组2
    ├── image1.tif
    └── ...
```

**Mode 2: Sample folders with single-channel files**
**模式二：单通道文件组成的样本文件夹**

```
data/
├── group1/
│   ├── sample1/
│   │   ├── DAPI.tif
│   │   ├── 488.tif
│   │   └── 594.tif
│   └── sample2/
│       └── ...
└── group2/
    └── ...
```

Supported formats / 支持格式: `.tif`, `.tiff`, `.png`, `.jpg`, `.jpeg`, `.bmp`

### 2. Configure the Script / 配置脚本

Edit the **User Settings** section at the top of `if_analysis.py`:
编辑脚本顶部的 **User Settings** 区域：

```python
# Input/output paths / 输入输出路径
GROUP1_DIR = r"/path/to/group1"
GROUP2_DIR = r"/path/to/group2"
GROUP1_NAME = "Control"         # 图表中组1的显示名称
GROUP2_NAME = "Treatment"       # 图表中组2的显示名称
OUTPUT_DIR = r"/path/to/output"

# Channel settings / 通道设置
DAPI_CHANNEL = "DAPI"
ANALYSIS_CHANNELS = ["488", "594", "647"]

# Manual channel order (when auto-detection fails)
# 手动指定通道顺序（当自动检测失败时使用）
FILE_CHANNEL_ORDER = None  # or ["DAPI", "488", "594", "647"]
```

### 3. Run Analysis / 运行分析

```bash
python if_analysis.py
```

Or use command line arguments to override config:
或使用命令行参数覆盖配置：

```bash
python if_analysis.py \
    --group1-dir /path/to/group1 \
    --group2-dir /path/to/group2 \
    --group1-name Control \
    --group2-name Treatment \
    --intensity-channels 488 647 \
    --coloc-pairs 488:647
```

---

## Configuration Guide / 配置指南

### Channel Configuration / 通道配置

```python
# Input discovery mode / 输入文件发现模式
INPUT_DISCOVERY_MODE = "auto"   # auto / sample_folders / image_files

# Auto channel detection settings / 自动通道检测设置
AUTO_DETECT_CHANNELS = True       # Enable auto-detection / 启用自动检测
STRICT_CHANNEL_DETECTION = True   # Error on failure (False uses fallback)
                                  # 检测失败时报错（False 则使用兜底顺序）
FILE_CHANNEL_ORDER = None         # Manual channel order override
                                  # 手动指定通道顺序

# Filename mapping (when filenames don't contain channel names)
# 文件名映射（当文件名不直接包含通道名时使用）
FILENAME_CHANNEL_MAP = {
    "CH1": "DAPI",
    "CH2": "488",
    "CH3": "594",
    "CH4": "647",
}

# Filename role mapping (ignore overlay/merge files)
# 文件角色映射（忽略 overlay/merge 文件）
FILENAME_ROLE_MAP = {
    "OVERLAY": "overlay",
    "MERGE": "overlay",
}
```

### Threshold Settings / 阳性阈值设置

```python
POSITIVE_THRESHOLD_RULES = {
    "DAPI": {"method": "otsu", "scale": 1.0, "min_value": 0.0},
    "594": {"method": "otsu", "scale": 1.7, "min_value": 30.0},
    "488": {"method": "otsu", "scale": 1.1, "min_value": 30.0},
    "647": {"method": "otsu", "scale": 1.3, "min_value": 30.0},
}
```

**Threshold methods / 阈值方法说明：**

| Method / 方法 | Description / 说明 | Use case / 适用场景 |
|---------------|-------------------|---------------------|
| `otsu` | Otsu auto-threshold × scale factor<br>Otsu 自动阈值 × scale 系数 | Stable signal/background contrast<br>信号/背景对比度稳定 |
| `manual` | Fixed threshold value<br>固定阈值 | Known optimal threshold<br>已知最佳阈值 |
| `percentile` | Percentile-based threshold<br>分位数阈值 | Variable signal intensity<br>信号强度变化大 |

### Background Subtraction / 背景扣除设置

```python
# Global default background percentile
# 全局默认背景分位数
BACKGROUND_PERCENTILE = 3.0

# Per-channel overrides / 按通道覆盖
BACKGROUND_PERCENTILE_RULES = {
    "DAPI": 1.0,   # DAPI uses lower background estimation
    "488": 3.0,
    "594": 3.0,
    "647": 3.0,
}
```

### Debris Filtering (Positive Object Area) / 碎点过滤

```python
# Global default (0 = disabled)
# 全局默认（0 表示关闭）
MIN_POSITIVE_OBJECT_AREA = 0

# Per-channel overrides / 按通道覆盖
POSITIVE_OBJECT_AREA_RULES = {
    "DAPI": 0,     # No filtering for DAPI
    "488": 100,    # Filter objects < 100 pixels for 488
    "594": 200,    # Filter objects < 200 pixels for 594
    "647": 100,
}
```

### Nuclei Segmentation / 细胞核分割设置

```python
# Nucleus segmentation parameters / 核分割参数
MIN_NUCLEUS_AREA = 0.5          # Minimum nucleus area in pixels
GAUSSIAN_BLUR_SIZE = 7          # Gaussian blur kernel size (odd number)
MIN_PEAK_DISTANCE = 1           # Minimum distance between watershed seeds
MASK_DILATION_RADIUS = 0        # ROI dilation radius beyond nuclei
```

**Recommended for 60X objective / 60X 物镜推荐参数：**
```python
MIN_NUCLEUS_AREA = 2200         # Remove debris
GAUSSIAN_BLUR_SIZE = 7          # Smooth noise
MIN_PEAK_DISTANCE = 25          # Prevent oversegmentation
MASK_DILATION_RADIUS = 18       # Include cytoplasm region
```

### Bleed-through Correction / 串色校正设置

```python
# Enable bleed-through correction / 启用串色校正
ENABLE_BLEEDTHROUGH_CORRECTION = True

# Define bleed-through direction: target <- source
# 定义串色方向：目标通道 <- 来源通道
BLEEDTHROUGH_SOURCE_MAP = {
    "488": "DAPI",  # Subtract DAPI bleed-through from 488
}

# Automatic estimation parameters / 自动估计参数
BLEEDTHROUGH_DEFAULT_RULE = {
    "mode": "auto",                    # auto / manual
    "estimate_mask": "nuclei",         # Estimation mask: nuclei/roi/all
    "ratio_percentile": 20.0,          # Low percentile for coefficient
    "source_threshold_percentile": 75.0,  # Source strong signal percentile
    "min_pixels": 100,                 # Minimum pixels for estimation
    "max_coefficient": 3.0,            # Maximum allowed coefficient
}

# Manual coefficients (override auto) / 手动指定系数
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    # "488": 0.15,  # Bleed-through coefficient for 488 <- DAPI
}
```

### Visualization / 可视化设置

```python
# Channel pseudocolors (RGB) / 通道伪彩色
CHANNEL_COLORS = {
    "DAPI": (0, 0, 255),      # Blue / 蓝色
    "488": (0, 255, 0),       # Green / 绿色
    "594": (255, 0, 0),       # Red / 红色
    "647": (255, 0, 178),     # Magenta / 品红
}

# Display modes / 显示模式
RAW_CHANNEL_DISPLAY_MODE = "color"       # raw: color / gray
POSITIVE_PREVIEW_DISPLAY_MODE = "color"  # preview: colormap / color / gray
POSITIVE_PREVIEW_COLORMAP = "inferno"    # Heatmap colormap

# Output control / 输出控制
SAVE_CHANNEL_DIAGNOSTICS = True
SAVE_QC_OVERLAYS = True
SAVE_BLEEDTHROUGH_DIAGNOSTICS = True
SAVE_FIGURE_PDF = True
SAVE_FIGURE_PNG = True
FIGURE_DPI = 500
```

---

## Output Structure / 输出结构

```
output/
└── Group1_vs_Group2_20260318_143022/
    ├── per_image_results.csv       # Per-image metrics / 单图像指标
    ├── group_statistics.csv        # Inter-group statistics / 组间统计
    ├── failed_images.csv           # Failed image log / 失败日志
    ├── analysis_report.txt         # Text report / 文本报告
    ├── run_config.json             # Config snapshot / 配置快照
    ├── summary_plots.pdf           # Summary plots (PDF)
    ├── summary_plots.png           # Summary plots (PNG)
    ├── qc_overlays/                # QC overlays / QC 叠加图
    │   └── ...
    ├── channel_diagnostics/        # Channel diagnostics / 通道诊断图
    │   └── ...
    ├── bleedthrough/               # Bleed-through diagnostics / 串色诊断图
    │   └── ...
    └── merged_denoised/            # Denoised merged images / 去噪合并图
        └── ...
```

---

## Results Interpretation / 结果解读

### Per-Image Metrics / 单图像指标

| Metric / 指标 | Description / 说明 |
|---------------|-------------------|
| `cell_count` | Number of detected nuclei / 检测到的细胞核数量 |
| `cell_density_per_mp` | Cell density per million pixels / 细胞密度（每百万像素） |
| `mean_nucleus_area_px` | Mean nucleus area in pixels / 平均核面积（像素） |
| `{channel}_roi_mean_bgsub` | ROI mean intensity (bg-subtracted) / ROI 平均强度（背景扣除后） |
| `{channel}_roi_integrated_bgsub` | ROI integrated intensity / ROI 累积强度 |
| `{channel}_integrated_per_cell` | Integrated intensity per cell / 每细胞累积强度 |
| `{channel}_norm_to_DAPI` | Intensity normalized to DAPI / 相对于 DAPI 的归一化强度 |
| `{channel}_positive_area_fraction` | Positive signal area fraction / 阳性信号面积比例 |

### Colocalization Metrics / 共定位指标

| Metric / 指标 | Range / 范围 | Description / 说明 |
|---------------|--------------|-------------------|
| `{A}_VS_{B}_pearson_r` | -1 ~ 1 | Pearson correlation coefficient |
| `{A}_VS_{B}_manders_{A}` | 0 ~ 1 | Manders' M1 (fraction of A overlapping B) |
| `{A}_VS_{B}_manders_{B}` | 0 ~ 1 | Manders' M2 (fraction of B overlapping A) |
| `{A}_VS_{B}_overlap_coeff` | 0 ~ 1 | Overlap coefficient |
| `{A}_VS_{B}_coloc_fraction` | 0 ~ 1 | Colocalized pixel fraction |

### Bleed-through Records / 串色校正记录

| Metric / 指标 | Description / 说明 |
|---------------|-------------------|
| `{channel}_bleedthrough_applied` | Whether correction was applied / 是否应用校正 |
| `{channel}_bleedthrough_source` | Source channel of crosstalk / 串色来源通道 |
| `{channel}_bleedthrough_coefficient` | Estimated coefficient / 估计的串色系数 |
| `{channel}_bleedthrough_mode` | auto / manual / 自动或手动模式 |

---

## Command Line Interface / 命令行参数

```bash
python if_analysis.py [options / 选项]

# Input/Output / 输入输出
--group1-dir PATH               Group 1 input directory / 组1输入目录
--group2-dir PATH               Group 2 input directory / 组2输入目录
--group1-name NAME              Group 1 display name / 组1名称
--group2-name NAME              Group 2 display name / 组2名称
--output-dir PATH               Output directory / 输出目录
--input-discovery-mode MODE     auto / sample_folders / image_files

# Channel Settings / 通道设置
--file-channel-order CH1 CH2... Manual channel order / 手动指定通道顺序
--fallback-channel-order ...    Fallback order / 兜底通道顺序
--no-auto-detect                Disable auto-detection / 关闭自动检测
--not-strict-channel-detection  Allow fallback on failure / 允许兜底
--dapi-channel NAME             DAPI channel name / DAPI 通道名
--intensity-channels CH1...     Channels for intensity analysis / 强度分析通道
--coloc-pairs A:B C:D...        Colocalization pairs / 共定位通道对
--all-coloc-pairs               Auto all pairwise colocalization / 自动两两配对
--no-auto-coloc-pairs           Disable auto pairing / 关闭自动配对

# Filename Mapping / 文件名映射
--filename-channel-map "CH1:DAPI CH2:488"
--filename-role-map "OVERLAY:ignore"
--recursive-scan                Recursive subdirectory scan / 递归扫描

# Image Processing / 图像处理参数
--background-percentile FLOAT   Background subtraction percentile
--min-nucleus-area INT          Minimum nucleus area / 最小核面积
--min-positive-object-area INT  Min positive object area / 最小阳性面积
--positive-object-area-rules "488:100 594:200"
--gaussian-blur-size INT        Gaussian blur kernel size / 高斯核大小
--min-peak-distance INT         Watershed seed distance / 种子最小距离
--mask-dilation-radius INT      ROI dilation radius / ROI 扩张半径

# Display Settings / 显示设置
--raw-channel-display-mode MODE     gray / color
--positive-preview-display-mode MODE  colormap / color / gray
--positive-preview-colormap NAME    inferno / magma / gray

# Feature Toggles / 功能开关
--no-channel-diagnostics        Disable channel diagnostics / 不保存通道诊断
--no-bleedthrough-correction    Disable bleed-through correction / 关闭串色校正
--no-bleedthrough-diagnostics   Disable bleed-through diagnostics
--bleedthrough-488-manual FLOAT Manual 488 coefficient / 488 手动系数
--no-qc                         Disable QC overlays / 不保存 QC 图
```

---

## Troubleshooting / 常见问题

### Channel Detection Errors / 通道检测失败

**Symptom / 症状：** `ValueError: Channel detection is ambiguous...`

**Solution / 解决：** Manually specify channel order
```python
FILE_CHANNEL_ORDER = ["DAPI", "488", "594", "647"]
# Or via CLI / 或使用命令行
python if_analysis.py --file-channel-order DAPI 488 594 647
```

### Oversegmentation / 细胞核过度分割

**Symptom / 症状：** One nucleus split into multiple fragments

**Solution / 解决：** Increase `MIN_PEAK_DISTANCE`
```python
MIN_PEAK_DISTANCE = 25  # Or higher
```

### High Background Noise / 背景噪声过多

**Symptom / 症状：** Too many debris in positive regions

**Solution / 解决：** Increase debris filtering area
```python
POSITIVE_OBJECT_AREA_RULES = {
    "488": 200,   # Increase filter threshold
    "594": 300,
}
```

### Excessive Bleed-through Correction / 串色校正过度

**Symptom / 症状：** Signal over-subtracted

**Solution / 解决：** Use manual coefficient or adjust parameters
```python
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    "488": 0.10,  # Manually specify coefficient
}
# Or adjust estimation / 或调整估计参数
BLEEDTHROUGH_DEFAULT_RULE["ratio_percentile"] = 10.0
```

### Missing tifffile / 缺少 tifffile

```bash
pip install tifffile
```

---

## Citation / 引用

If you use this tool in your research, please cite:
如果您在研究中使用了本工具，请引用：

```bibtex
@software{easy_immunofluorescence_analysis,
  title = {Easy Immunofluorescence Analysis: Batch IF Image Analysis Tool},
  author = {zhumiao-cloud},
  year = {2026},
  url = {https://github.com/zhumiao-cloud/easy-immunofluorance-analysis}
}
```

---

## License / 许可证

MIT License - see [LICENSE](LICENSE) file for details.

---

## Changelog / 更新日志

See [CHANGELOG.md](CHANGELOG.md) for version history.
详见 [CHANGELOG.md](CHANGELOG.md)。

---

## Contributing / 贡献与反馈

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue.
欢迎提交 Pull Request 或 Issue。
