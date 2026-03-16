# Easy IF Analysis

Easy batch immunofluorescence (IF) image analysis tool with automatic channel detection, colocalization analysis, and statistical comparison.

一个简单易用的免疫荧光图像批量分析工具，支持自动通道识别、共定位分析和统计学比较。

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## ✨ Features | 功能特性

- 🔬 **Automatic Channel Detection** - 自动识别 DAPI/488/594/647 通道
- 🧬 **Nuclei Segmentation** - 基于 DAPI 的细胞核分割（分水岭算法）
- 📊 **Intensity Analysis** - 信号强度和阳性面积分数分析
- 🔗 **Colocalization Analysis** - 共定位分析（Pearson 相关系数、Manders 系数）
- 🎨 **Bleed-through Correction** - DAPI 串色校正
- 📈 **Statistical Comparison** - 两组统计学比较（t-test / Mann-Whitney U）
- 🖼️ **QC Visualization** - 自动生成 QC 诊断图

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

```
data/
├── group1/          # 对照组
│   ├── image1.tif
│   ├── image2.tif
│   └── ...
└── group2/          # 实验组
    ├── image1.tif
    ├── image2.tif
    └── ...
```

支持的格式：`.tif`, `.tiff`, `.png`, `.jpg`, `.bmp`

### 2. Configure the Script | 配置脚本

编辑脚本顶部的 **Configuration Section**（用户配置区）：

```python
# ============================================================================
# 用户配置区：绝大多数情况下，只需要修改这一段
# ============================================================================

GROUP1_DIR = r"/path/to/your/group1"
GROUP2_DIR = r"/path/to/your/group2"
GROUP1_NAME = "Control"
GROUP2_NAME = "Treatment"
OUTPUT_DIR = r"/path/to/output"

# 手动指定文件中的通道顺序（可选）
# 如果不确定，设为 None，让脚本自动检测
FILE_CHANNEL_ORDER = None  # 或 ["647", "488", "DAPI"]

# 需要分析的通道
INTENSITY_CHANNELS = ["488", "594"]  # 或 ["488", "647"]

# 共定位分析对
COLOCALIZATION_PAIRS = [("488", "594")]  # 或 [("488", "647")]
```

### 3. Run the Analysis | 运行分析

```bash
python if_analysis.py
```

Or with command line arguments:
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

## ⚙️ Configuration Guide | 配置详解

### Channel Configuration | 通道配置

| Parameter | Description | Example |
|-----------|-------------|---------|
| `FILE_CHANNEL_ORDER` | 手动指定通道顺序 | `["647", "488", "DAPI"]` |
| `DAPI_CHANNEL` | DAPI 通道名称 | `"DAPI"` |
| `INTENSITY_CHANNELS` | 需要分析强度的通道 | `["488", "594"]` |
| `COLOCALIZATION_PAIRS` | 共定位分析对 | `[("488", "594")]` |

### Threshold Settings | 阈值设置

```python
# 自定义阳性阈值（可选）
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 1.6, "min_value": 30.0}
}
```

Methods available:
- `otsu`: Otsu 自动阈值 × scale
- `manual`: 手动指定阈值
- `percentile`: 基于分位数的阈值

### Bleed-through Correction | 串色校正

```python
ENABLE_BLEEDTHROUGH_CORRECTION = True

BLEEDTHROUGH_RULES = {
    "488": {
        "source": "DAPI",           # 从 DAPI 扣除
        "mode": "auto",             # 自动估计系数
        "estimate_mask": "nuclei",  # 在细胞核区域估计
    }
}
```

### Analysis Parameters | 分析参数

```python
BACKGROUND_PERCENTILE = 3      # 背景估计分位数
MIN_NUCLEUS_AREA = 50          # 最小细胞核面积（像素）
GAUSSIAN_BLUR_SIZE = 7         # 高斯模糊核大小
MIN_PEAK_DISTANCE = 10         # 分水岭种子最小距离
MASK_DILATION_RADIUS = 6       # ROI 扩张半径
```

> 💡 **60X 物镜建议参数**：
> - `MIN_NUCLEUS_AREA = 2200`
> - `MIN_PEAK_DISTANCE = 25`
> - `MASK_DILATION_RADIUS = 18`
> - `BACKGROUND_PERCENTILE = 5`

---

## 📁 Output Structure | 输出结构

```
output/
└── Group1_vs_Group2_20260316_143022/
    ├── per_image_results.csv       # 每张图像的详细结果
    ├── group_statistics.csv        # 组间统计分析
    ├── failed_images.csv           # 分析失败的图像
    ├── analysis_report.txt         # 文本报告
    ├── run_config.json             # 运行配置记录
    ├── summary_plots.pdf           # 汇总图表
    ├── summary_plots.png
    ├── qc_overlays/                # QC 质控图
    │   ├── image1_qc.png
    │   └── ...
    └── channel_diagnostics/        # 通道诊断图
        ├── image1_channels.png
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
| `{channel}_roi_mean_bgsub` | ROI 内背景扣除后的平均强度 |
| `{channel}_positive_area_fraction` | 阳性区域比例 |
| `{channel}_norm_to_DAPI` | 相对于 DAPI 的归一化强度 |

### Colocalization Metrics | 共定位指标

| Metric | Range | Description |
|--------|-------|-------------|
| `pearson_r` | -1 to 1 | Pearson 相关系数 |
| `manders_M1` | 0 to 1 | Manders M1 系数（Channel A 与 B 重叠比例） |
| `manders_M2` | 0 to 1 | Manders M2 系数（Channel B 与 A 重叠比例） |
| `overlap_coeff` | 0 to 1 | 重叠系数 |
| `coloc_fraction` | 0 to 1 | 共定位像素比例 |

---

## 🛠️ Advanced Usage | 高级用法

### Command Line Arguments | 命令行参数

```bash
python if_analysis.py \
    --group1-dir PATH \
    --group2-dir PATH \
    --group1-name NAME \
    --group2-name NAME \
    --output-dir PATH \
    --file-channel-order 647 488 DAPI \
    --intensity-channels 488 647 \
    --coloc-pairs 488:647 \
    --background-percentile 5 \
    --min-nucleus-area 2200 \
    --gaussian-blur-size 7 \
    --min-peak-distance 25 \
    --mask-dilation-radius 18 \
    --no-qc \
    --recursive-scan
```

### Custom Channel Names | 自定义通道名

如果你的通道名称不同（如使用 Cy3、Cy5），可以在配置中指定：

```python
canonicalize_channel_label("Cy3")   # 返回 "594"
canonicalize_channel_label("Cy5")   # 返回 "647"
```

脚本支持的别名：
- **DAPI**: DAPI, Hoechst, 405, Blue, Nuclei
- **488**: 488, FITC, GFP, Green, Alexa488
- **594**: 594, TRITC, Cy3, TexasRed, RFP, 561, 568
- **647**: 647, Cy5, APC, Alexa647, FarRed, 640

---

## 🐛 Troubleshooting | 常见问题

### Q: 通道识别错误
**A**: 手动设置 `FILE_CHANNEL_ORDER`:
```python
FILE_CHANNEL_ORDER = ["594", "488", "DAPI"]
```

### Q: 细胞核分割不准确
**A**: 调整分割参数：
```python
MIN_NUCLEUS_AREA = 2200        # 增大以去除小噪声
GAUSSIAN_BLUR_SIZE = 9         # 增大以平滑噪声
MIN_PEAK_DISTANCE = 30         # 增大以避免过度分割
```

### Q: 594 假阳性过高
**A**: 提高阈值：
```python
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 2.0, "min_value": 50.0}
}
```

### Q: 缺少 tifffile
**A**: 
```bash
pip install tifffile
```

---

## 📄 Citation | 引用

If you use this tool in your research, please cite:

```bibtex
@software{easy_if_analysis,
  title = {Easy IF Analysis: Batch Immunofluorescence Image Analysis Tool},
  author = {Your Name},
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
