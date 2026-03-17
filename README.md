# Easy Immunofluorescence Analysis

Easy batch immunofluorescence (IF) image analysis tool with automatic channel detection, colocalization analysis, bleed-through correction, and statistical comparison.

一个简单易用的免疫荧光图像批量分析工具，支持自动通道识别、共定位分析、串色校正和统计学比较。

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## ✨ Features | 功能特性

- 🔬 **Automatic Channel Detection** - 自动识别 DAPI/488/594/647 通道
- 🧬 **Nuclei Segmentation** - 基于 DAPI 的细胞核分割（分水岭算法）
- 📊 **Intensity Analysis** - 信号强度和阳性面积分数分析
- 🔗 **Colocalization Analysis** - 共定位分析（Pearson 相关系数、Manders 系数、Overlap 系数）
- 🎨 **Bleed-through Correction** - 多通道串色校正（自动/手动系数估计）
- 🖼️ **Flexible Visualization** - 灵活的可视化（伪彩色/colormap/灰度）
- 📈 **Statistical Comparison** - 两组统计学比较（Welch t-test / Mann-Whitney U）
- ✅ **QC Visualization** - 自动生成 QC 诊断图和通道诊断图

---

## 📦 Installation | 安装

### Requirements | 环境要求

```bash
Python >= 3.8
```

### Install Dependencies | 安装依赖

```bash
pip install numpy pandas scipy scikit-image opencv-python-headless matplotlib tifffile pillow
```

Or install from requirements.txt:
```bash
pip install -r requirements.txt
```

---

## 🚀 Quick Start | 快速开始

### 1. Prepare Your Images | 准备图像

将你的免疫荧光图像放在两个文件夹中（例如实验组和对照组）：

**多通道图像模式**（单个文件包含多个通道）：
```
data/
├── group1/          # 对照组
│   ├── image1.tif   # 4通道 TIFF: DAPI/488/594/647
│   ├── image2.tif
│   └── ...
└── group2/          # 实验组
    ├── image1.tif
    └── ...
```

**样本文件夹模式**（每个通道单独文件）：
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

支持的格式：`.tif`, `.tiff`, `.png`, `.jpg`, `.jpeg`, `.bmp`

### 2. Configure the Script | 配置脚本

编辑脚本顶部的 **User Settings**（用户配置区）：

```python
# ============================================================================
# User Settings: edit this section
# ============================================================================

# 输入输出路径
GROUP1_DIR = r"/path/to/your/group1"
GROUP2_DIR = r"/path/to/your/group2"
GROUP1_NAME = "Control"
GROUP2_NAME = "Treatment"
OUTPUT_DIR = r"/path/to/output"

# 通道配置
DAPI_CHANNEL = "DAPI"
ANALYSIS_CHANNELS = ["488", "594", "647"]

# 手动指定通道顺序（可选）
FILE_CHANNEL_ORDER = None  # 或 ["647", "488", "DAPI"]

# 共定位分析
AUTO_ALL_PAIRWISE_STATS = True  # 自动所有两两配对
```

### 3. Run the Analysis | 运行分析

```bash
python if_analysis.py
```

Or with command line arguments:
```bash
python if_analysis.py     --group1-dir /path/to/group1     --group2-dir /path/to/group2     --group1-name Control     --group2-name Treatment     --intensity-channels 488 647     --coloc-pairs 488:647
```

---

## ⚙️ Configuration Guide | 配置详解

### Channel Configuration | 通道配置

| Parameter | Description | Example |
|-----------|-------------|---------|
| `FILE_CHANNEL_ORDER` | 手动指定通道顺序 | `["647", "488", "DAPI"]` |
| `DAPI_CHANNEL` | DAPI 通道名称 | `"DAPI"` |
| `ANALYSIS_CHANNELS` | 需要分析的通道 | `["488", "594", "647"]` |
| `AUTO_ALL_PAIRWISE_STATS` | 自动两两共定位分析 | `True` / `False` |
| `MANUAL_COLOCALIZATION_PAIRS` | 手动指定共定位对 | `[("488", "594")]` |

### Threshold & Filtering Settings | 阈值与过滤设置

```python
# 阳性阈值规则（支持 Otsu、Manual、Percentile 三种方法）
POSITIVE_THRESHOLD_RULES = {
    "DAPI": {"method": "otsu", "scale": 1.0, "min_value": 0.0},
    "594": {"method": "otsu", "scale": 1.7, "min_value": 30.0},
    "488": {"method": "otsu", "scale": 1.1, "min_value": 30.0},
    "647": {"method": "otsu", "scale": 1.3, "min_value": 30.0},
}

# 阳性区域最小面积（去除碎点）
MIN_POSITIVE_OBJECT_AREA = 0  # 全局默认值
POSITIVE_OBJECT_AREA_RULES = {  # 按通道覆盖
    "DAPI": 0,
    "488": 100,
    "594": 200,
    "647": 100,
}
```

### Background Subtraction | 背景扣除

```python
# 背景估计分位数（全局默认值）
BACKGROUND_PERCENTILE = 3.0

# 按通道覆盖（DAPI 可设更低）
BACKGROUND_PERCENTILE_RULES = {
    "DAPI": 1.0,
    "488": 3.0,
    "594": 3.0,
    "647": 3.0,
}
```

### Bleed-through Correction | 串色校正

```python
ENABLE_BLEEDTHROUGH_CORRECTION = True

# 定义 target <- source 的串色关系
BLEEDTHROUGH_SOURCE_MAP = {
    "488": DAPI_CHANNEL,  # 488 从 DAPI 扣除
    # "594": DAPI_CHANNEL,
    # "647": "488",
}

# 自动估计参数
BLEEDTHROUGH_DEFAULT_RULE = {
    "mode": "auto",                    # "auto" 或 "manual"
    "estimate_mask": "nuclei",         # 估计区域："nuclei"/"roi"/"all"
    "ratio_percentile": 20.0,          # 取较低分位避免过度扣除
    "source_threshold_percentile": 75.0,  # 只在 source 较强像素上估计
    "min_pixels": 100,
    "max_coefficient": 3.0,
}

# 手动指定系数（覆盖自动估计）
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    # "488": 0.15,
}
```

### Visualization Settings | 可视化设置

```python
# 通道伪彩色（RGB）
CHANNEL_COLORS = {
    "DAPI": (0, 0, 255),      # 蓝色
    "488": (0, 255, 0),       # 绿色
    "594": (255, 0, 0),       # 红色
    "647": (255, 0, 178),     # 品红
}

# 原始单通道显示方式："color"（伪彩色）或 "gray"（灰度）
RAW_CHANNEL_DISPLAY_MODE = "color"

# 阳性预览显示方式："colormap" / "color" / "gray"
POSITIVE_PREVIEW_DISPLAY_MODE = "color"

# colormap 名称（仅在 colormap 模式下使用）
POSITIVE_PREVIEW_COLORMAP = "inferno"
```

### Nuclei Segmentation Parameters | 细胞核分割参数

```python
MIN_NUCLEUS_AREA = 0.5           # 最小细胞核面积（像素）
GAUSSIAN_BLUR_SIZE = 7           # 高斯模糊核大小
MIN_PEAK_DISTANCE = 1            # 分水岭种子最小距离
MASK_DILATION_RADIUS = 0         # ROI 相对 nuclei 向外扩张的像素半径
```

> 💡 **物镜参数建议**：
> 
> | 参数 | 20X | 60X |
> |------|-----|-----|
> | `MIN_NUCLEUS_AREA` | 50 | 2200 |
> | `MIN_PEAK_DISTANCE` | 6 | 25 |
> | `MASK_DILATION_RADIUS` | 6 | 18 |
> | `BACKGROUND_PERCENTILE` | 3 | 5 |

---

## 📁 Output Structure | 输出结构

```
output/
└── Group1_vs_Group2_20260316_143022/
    ├── per_image_results.csv       # 每张图像的详细结果
    ├── group_statistics.csv        # 组间统计分析
    ├── failed_images.csv           # 分析失败的图像
    ├── analysis_report.txt         # 文本报告（含配置记录）
    ├── run_config.json             # 运行配置（JSON格式）
    ├── summary_plots.pdf           # 汇总图表
    ├── summary_plots.png
    ├── qc_overlays/                # QC 质控图
    │   ├── image1_qc.png
    │   └── ...
    ├── channel_diagnostics/        # 通道诊断图
    │   ├── image1_channels.png
    │   └── ...
    ├── merged_denoised/            # 降噪后的合并图
    │   └── ...
    └── channel_diagnostics/bleedthrough/  # 串色校正诊断图
        └── ...
```

---

## 📊 Results Interpretation | 结果解读

### Per-Image Metrics | 单图像指标

| Metric | Description |
|--------|-------------|
| `cell_count` | 细胞核数量 |
| `cell_density_per_mp` | 细胞密度（每百万像素） |
| `mean_nucleus_area_px` | 平均细胞核面积（像素） |
| `roi_area_fraction` | ROI 占图像面积比例 |
| `{channel}_roi_mean_bgsub` | ROI 内背景扣除后的平均强度 |
| `{channel}_roi_integrated_bgsub` | ROI 内背景扣除后的总强度 |
| `{channel}_positive_area_fraction` | 阳性区域比例（经阈值和面积过滤） |
| `{channel}_norm_to_DAPI` | 相对于 DAPI 的归一化强度 |
| `{channel}_integrated_per_cell` | 每个细胞的平均强度 |

### Colocalization Metrics | 共定位指标

| Metric | Range | Description |
|--------|-------|-------------|
| `pearson_r` | -1 to 1 | Pearson 相关系数（线性相关性） |
| `manders_M1` | 0 to 1 | Manders M1 系数（Channel A 重叠于 B 的比例） |
| `manders_M2` | 0 to 1 | Manders M2 系数（Channel B 重叠于 A 的比例） |
| `overlap_coeff` | 0 to 1 | 重叠系数（Overlap coefficient） |
| `coloc_fraction` | 0 to 1 | 共定位像素占联合阳性区域的比例 |

### Statistics | 统计学结果

- **Normality Test**: Shapiro-Wilk (n >= 3)
- **Group Comparison**: 
  - 正态分布：Welch t-test
  - 非正态：Mann-Whitney U test
- **Multiple Testing Correction**: Benjamini-Hochberg FDR
- **Effect Size**: Hedges' g (t-test) / Rank-biserial r (Mann-Whitney)

---

## 🛠️ Advanced Usage | 高级用法

### Command Line Arguments | 命令行参数

```bash
python if_analysis.py     --group1-dir PATH     --group2-dir PATH     --group1-name NAME     --group2-name NAME     --output-dir PATH     --file-channel-order 647 488 DAPI     --intensity-channels 488 647     --coloc-pairs 488:647     --background-percentile 5     --min-nucleus-area 2200     --min-positive-object-area 200     --positive-object-area-rules 488:100 594:200     --gaussian-blur-size 7     --min-peak-distance 25     --mask-dilation-radius 18     --raw-channel-display-mode color     --positive-preview-display-mode colormap     --positive-preview-colormap inferno     --no-qc     --no-channel-diagnostics     --no-bleedthrough-correction     --recursive-scan
```

### Input Discovery Modes | 输入发现模式

```python
INPUT_DISCOVERY_MODE = "auto"  # 自动检测输入结构
```

- `"auto"` - 自动识别：检测到单通道图则使用样本文件夹模式，否则使用直接图像模式
- `"sample_folders"` - 只分析子文件夹（每个子文件夹包含多个单通道文件）
- `"image_files"` - 只分析直接图像文件（多通道 TIFF 等）

### Custom Channel Names | 自定义通道名

如果你的文件名使用不同的通道命名（如 CH1、CH2）：

```python
FILENAME_CHANNEL_MAP = {
    "CH1": "DAPI",
    "CH2": "488",
    "CH3": "594",
    "CH4": "647",
}
```

忽略 overlay/merge 文件：
```python
FILENAME_ROLE_MAP = {
    "OVERLAY": "overlay",
    "MERGE": "overlay",
    "RGB": "overlay",
}
```

支持的通道别名：
- **DAPI**: DAPI, Hoechst, Hoechst33342, 405, Blue, Nuclei
- **488**: 488, FITC, GFP, EGFP, Green, Alexa488, AF488
- **594**: 594, TRITC, Cy3, TexasRed, TXRed, RFP, 561, 568, Alexa594, AF594
- **647**: 647, Cy5, APC, Alexa647, AF647, FarRed, 640

---

## 🐛 Troubleshooting | 常见问题

### Q: 通道识别错误
**A**: 手动设置 `FILE_CHANNEL_ORDER`:
```python
FILE_CHANNEL_ORDER = ["594", "488", "DAPI"]
```
或检查 `channel_diagnostics/` 中的诊断图查看实际通道顺序。

### Q: 细胞核分割不准确
**A**: 调整分割参数：
```python
MIN_NUCLEUS_AREA = 2200        # 增大以去除小噪声
GAUSSIAN_BLUR_SIZE = 9         # 增大以平滑噪声
MIN_PEAK_DISTANCE = 30         # 增大以避免过度分割
```

### Q: 594 假阳性过高
**A**: 提高阈值并增加最小面积过滤：
```python
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 2.0, "min_value": 50.0}
}
POSITIVE_OBJECT_AREA_RULES = {
    "594": 200  # 去除小碎点
}
```

### Q: 串色校正过度/不足
**A**: 切换到手动模式或调整参数：
```python
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    "488": 0.12  # 手动指定系数
}
# 或调整自动估计参数
BLEEDTHROUGH_DEFAULT_RULE["ratio_percentile"] = 10.0  # 更保守的估计
```

### Q: 缺少 tifffile
**A**: 
```bash
pip install tifffile
```

---

## 📝 Version History | 版本历史

- **v1.0.6** - 新增阳性区域面积过滤、按通道背景分位数、可视化模式选择、改进串色校正
- **v1.0.5** - 修复通道识别、增强 QC 输出、添加串色校正

---

## 📄 Citation | 引用

If you use this tool in your research, please cite:

```bibtex
@software{easy_immunofluorescence_analysis,
  title = {Easy Immunofluorescence Analysis: Batch Immunofluorescence Image Analysis Tool},
  author = {zhumiao-cloud},
  year = {2026},
  url = {https://github.com/zhumiao-cloud/easy-immunofluorance-analysis}
}
```

---

## 📜 License | 许可证

MIT License - see [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing | 贡献

Contributions are welcome! Please feel free to submit a Pull Request.

---

## 📧 Contact | 联系方式

For questions or suggestions, please open an [Issue](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues).
