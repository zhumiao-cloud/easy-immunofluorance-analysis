# Easy Immunofluorescence Analysis

批量免疫荧光图像分析工具，支持自动通道检测、共定位分析、串色校正和组间统计比较。

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 功能特性

| 功能 | 说明 |
|------|------|
| **自动通道检测** | 从 OME-TIFF 元数据、文件名或启发式算法自动识别 DAPI/488/594/647 通道 |
| **灵活输入模式** | 支持多通道单文件，或单通道文件组成的样本文件夹 |
| **细胞核分割** | 基于 Watershed 的细胞核分割，参数可配置 |
| **强度定量** | 背景扣除后的信号强度和阳性面积分数 |
| **共定位分析** | Pearson 相关系数、Manders' M1/M2、重叠系数 |
| **串色校正** | 自动或手动估计串色系数，消除通道间串扰 |
| **碎点过滤** | 阈值后按面积过滤小物体，减少背景噪声 |
| **组间统计** | 两组比较（Welch t-test / Mann-Whitney U），支持 FDR 校正 |
| **可视化输出** | 自动生成 QC 叠加图、通道诊断图、汇总统计图 |

---

## 安装

### 环境要求

- Python >= 3.8
- 推荐：使用 conda 环境

### 安装依赖

```bash
pip install numpy pandas scipy scikit-image opencv-python-headless matplotlib tifffile pillow
```

或使用 requirements.txt：
```bash
pip install -r requirements.txt
```

---

## 快速开始

### 1. 准备数据

**模式一：多通道单文件**（推荐，如 Zeiss/Thermo 导出的 TIFF）
```
data/
├── group1/              # 实验组1
│   ├── image1.tif       # 4通道 TIFF (DAPI/488/594/647)
│   └── image2.tif
└── group2/              # 实验组2
    ├── image1.tif
    └── ...
```

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

### 2. 配置脚本

编辑脚本顶部的 **User Settings** 区域：

```python
# 输入输出路径
GROUP1_DIR = r"/path/to/group1"
GROUP2_DIR = r"/path/to/group2"
GROUP1_NAME = "Control"
GROUP2_NAME = "Treatment"
OUTPUT_DIR = r"/path/to/output"

# 通道设置
DAPI_CHANNEL = "DAPI"
ANALYSIS_CHANNELS = ["488", "594", "647"]

# 手动指定通道顺序（当自动检测失败时使用）
FILE_CHANNEL_ORDER = None  # 或 ["DAPI", "488", "594", "647"]
```

### 3. 运行分析

```bash
python if_analysis.py
```

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

## 详细配置指南

### 通道配置

```python
# 输入文件发现模式
INPUT_DISCOVERY_MODE = "auto"  # auto / sample_folders / image_files

# 自动通道检测设置
AUTO_DETECT_CHANNELS = True       # 启用自动检测
STRICT_CHANNEL_DETECTION = True   # 检测失败时报错（False 则使用兜底顺序）
FILE_CHANNEL_ORDER = None         # 手动指定通道顺序

# 文件名映射（当文件名不直接包含通道名时使用）
FILENAME_CHANNEL_MAP = {
    "CH1": "DAPI",
    "CH2": "488",
    "CH3": "594",
    "CH4": "647",
}

# 文件角色映射（忽略 overlay/merge 文件）
FILENAME_ROLE_MAP = {
    "OVERLAY": "overlay",
    "MERGE": "overlay",
}
```

### 阳性阈值设置

```python
POSITIVE_THRESHOLD_RULES = {
    "DAPI": {"method": "otsu", "scale": 1.0, "min_value": 0.0},
    "594": {"method": "otsu", "scale": 1.7, "min_value": 30.0},
    "488": {"method": "otsu", "scale": 1.1, "min_value": 30.0},
    "647": {"method": "otsu", "scale": 1.3, "min_value": 30.0},
}
```

**阈值方法说明：**
| 方法 | 说明 | 适用场景 |
|------|------|----------|
| `otsu` | Otsu 自动阈值 × scale 系数 | 信号/背景对比度稳定 |
| `manual` | 固定阈值 | 已知最佳阈值 |
| `percentile` | 分位数阈值 | 信号强度变化大 |

### 背景扣除设置

```python
# 全局默认背景分位数
BACKGROUND_PERCENTILE = 3.0

# 按通道覆盖
BACKGROUND_PERCENTILE_RULES = {
    "DAPI": 1.0,   # DAPI 使用更低的背景估计
    "488": 3.0,
    "594": 3.0,
    "647": 3.0,
}
```

### 碎点过滤（Positive Object Area）

```python
# 全局默认（0 表示关闭）
MIN_POSITIVE_OBJECT_AREA = 0

# 按通道覆盖
POSITIVE_OBJECT_AREA_RULES = {
    "DAPI": 0,     # DAPI 不过滤
    "488": 100,    # 488 通道过滤 < 100 像素的碎点
    "594": 200,    # 594 通道过滤 < 200 像素的碎点
    "647": 100,
}
```

### 细胞核分割设置

```python
# 核分割参数
MIN_NUCLEUS_AREA = 0.5          # 最小核面积（像素）
GAUSSIAN_BLUR_SIZE = 7          # 高斯平滑核大小（奇数）
MIN_PEAK_DISTANCE = 1           # Watershed 种子点最小距离
MASK_DILATION_RADIUS = 0        # ROI 相对核向外扩张半径
```

**60X 物镜推荐参数：**
```python
MIN_NUCLEUS_AREA = 2200         # 去除碎片
GAUSSIAN_BLUR_SIZE = 7          # 平滑噪声
MIN_PEAK_DISTANCE = 25          # 防止过度分割
MASK_DILATION_RADIUS = 18       # 包含细胞质区域
```

### 串色校正设置

```python
# 启用串色校正
ENABLE_BLEEDTHROUGH_CORRECTION = True

# 定义串色方向：target <- source
BLEEDTHROUGH_SOURCE_MAP = {
    "488": "DAPI",  # 488 通道扣除来自 DAPI 的串色
    # "594": "DAPI",
    # "647": "488",
}

# 自动估计参数
BLEEDTHROUGH_DEFAULT_RULE = {
    "mode": "auto",                    # auto / manual
    "estimate_mask": "nuclei",         # 估计 mask: nuclei / roi / all
    "ratio_percentile": 20.0,          # 低分位数估计系数
    "source_threshold_percentile": 75.0,  # source 强信号分位数
    "min_pixels": 100,                 # 最小估计像素数
    "max_coefficient": 3.0,            # 最大允许系数
}

# 手动指定系数（覆盖 auto）
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    # "488": 0.15,  # 488 <- DAPI 的串色系数
}
```

### 可视化设置

```python
# 通道伪彩色（RGB）
CHANNEL_COLORS = {
    "DAPI": (0, 0, 255),      # 蓝色
    "488": (0, 255, 0),       # 绿色
    "594": (255, 0, 0),       # 红色
    "647": (255, 0, 178),     # 品红
}

# 显示模式
RAW_CHANNEL_DISPLAY_MODE = "color"       # 原始通道: color / gray
POSITIVE_PREVIEW_DISPLAY_MODE = "color"  # 阳性预览: colormap / color / gray
POSITIVE_PREVIEW_COLORMAP = "inferno"    # 热力图 colormap

# 输出控制
SAVE_CHANNEL_DIAGNOSTICS = True          # 保存通道诊断图
SAVE_QC_OVERLAYS = True                  # 保存 QC 叠加图
SAVE_BLEEDTHROUGH_DIAGNOSTICS = True     # 保存串色校正诊断图
SAVE_FIGURE_PDF = True                   # 保存汇总 PDF
SAVE_FIGURE_PNG = True                   # 保存汇总 PNG
FIGURE_DPI = 500                         # 输出分辨率
```

---

## 输出结构

```
output/
└── Group1_vs_Group2_20260318_143022/
    ├── per_image_results.csv       # 每张图像的详细指标
    ├── group_statistics.csv        # 组间统计分析
    ├── failed_images.csv           # 失败图像日志
    ├── analysis_report.txt         # 文本报告
    ├── run_config.json             # 配置快照
    ├── summary_plots.pdf           # 汇总图（PDF）
    ├── summary_plots.png           # 汇总图（PNG）
    ├── qc_overlays/                # QC 叠加图
    │   ├── Group1__image1_qc.png
    │   └── ...
    ├── channel_diagnostics/        # 通道诊断图
    │   ├── Group1__image1_channels.png
    │   └── ...
    ├── bleedthrough/               # 串色校正诊断图
    │   └── ...
    └── merged_denoised/            # 去噪后的合并图
        └── ...
```

---

## 结果解读

### 单图像指标

| 指标 | 说明 |
|------|------|
| `cell_count` | 检测到的细胞核数量 |
| `cell_density_per_mp` | 细胞密度（每百万像素） |
| `mean_nucleus_area_px` | 平均核面积（像素） |
| `{channel}_roi_mean_bgsub` | ROI 平均强度（背景扣除后） |
| `{channel}_roi_integrated_bgsub` | ROI 累积强度（背景扣除后） |
| `{channel}_integrated_per_cell` | 每细胞累积强度 |
| `{channel}_norm_to_DAPI` | 相对于 DAPI 的归一化强度 |
| `{channel}_positive_area_fraction` | 阳性信号面积比例 |

### 共定位指标

| 指标 | 范围 | 说明 |
|------|------|------|
| `{A}_VS_{B}_pearson_r` | -1 ~ 1 | Pearson 相关系数 |
| `{A}_VS_{B}_manders_{A}` | 0 ~ 1 | Manders' M1（A 中与 B 重叠的比例） |
| `{A}_VS_{B}_manders_{B}` | 0 ~ 1 | Manders' M2（B 中与 A 重叠的比例） |
| `{A}_VS_{B}_overlap_coeff` | 0 ~ 1 | 重叠系数 |
| `{A}_VS_{B}_coloc_fraction` | 0 ~ 1 | 共定位像素比例 |

### 串色校正记录

| 指标 | 说明 |
|------|------|
| `{channel}_bleedthrough_applied` | 是否应用校正 |
| `{channel}_bleedthrough_source` | 串色来源通道 |
| `{channel}_bleedthrough_coefficient` | 估计的串色系数 |
| `{channel}_bleedthrough_mode` | auto / manual |

---

## 命令行参数

```bash
python if_analysis.py [选项]

# 输入输出
--group1-dir PATH               组1输入目录
--group2-dir PATH               组2输入目录
--group1-name NAME              组1名称
--group2-name NAME              组2名称
--output-dir PATH               输出目录
--input-discovery-mode MODE     输入模式 (auto/sample_folders/image_files)

# 通道设置
--file-channel-order CH1 CH2... 手动指定通道顺序
--fallback-channel-order ...    兜底通道顺序
--no-auto-detect                关闭自动通道检测
--not-strict-channel-detection  允许检测失败后使用兜底顺序
--dapi-channel NAME             DAPI 通道名
--intensity-channels CH1 CH2... 强度分析通道
--coloc-pairs A:B C:D ...       共定位通道对
--all-coloc-pairs               自动所有两两配对
--no-auto-coloc-pairs           关闭自动配对

# 文件名映射
--filename-channel-map "CH1:DAPI CH2:488"
--filename-role-map "OVERLAY:ignore MERGE:ignore"
--recursive-scan                递归扫描子目录

# 图像处理参数
--background-percentile FLOAT   背景扣除分位数
--min-nucleus-area INT          最小核面积
--min-positive-object-area INT  阈值后最小阳性面积
--positive-object-area-rules "488:100 594:200"
--gaussian-blur-size INT        高斯平滑核大小
--min-peak-distance INT         Watershed 种子最小距离
--mask-dilation-radius INT      ROI 扩张半径

# 显示设置
--raw-channel-display-mode MODE     gray / color
--positive-preview-display-mode MODE  colormap / color / gray
--positive-preview-colormap NAME    inferno / magma / gray

# 功能开关
--no-channel-diagnostics        不保存通道诊断图
--no-bleedthrough-correction    关闭串色校正
--no-bleedthrough-diagnostics   不保存串色诊断图
--bleedthrough-488-manual FLOAT 488 通道手动串色系数
--no-qc                         不保存 QC 叠加图
```

---

## 常见问题

### 通道检测失败

**症状：** `ValueError: 通道仍然存在歧义...`

**解决：** 手动指定通道顺序
```python
FILE_CHANNEL_ORDER = ["DAPI", "488", "594", "647"]
# 或使用命令行
python if_analysis.py --file-channel-order DAPI 488 594 647
```

### 细胞核过度分割

**症状：** 一个核被分成多个

**解决：** 增加 `MIN_PEAK_DISTANCE`
```python
MIN_PEAK_DISTANCE = 25  # 或更大
```

### 背景噪声过多

**症状：** 阳性区域有太多碎点

**解决：** 增加碎点过滤面积
```python
POSITIVE_OBJECT_AREA_RULES = {
    "488": 200,   # 增加过滤阈值
    "594": 300,
}
```

### 串色校正过度

**症状：** 信号被过度扣除

**解决：** 使用手动系数或调整参数
```python
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    "488": 0.10,  # 手动指定系数
}
# 或调整估计参数
BLEEDTHROUGH_DEFAULT_RULE["ratio_percentile"] = 10.0  # 更保守的估计
```

### 缺少 tifffile

```bash
pip install tifffile
```

---

## 引用

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

## 许可证

MIT License - 详见 [LICENSE](LICENSE) 文件。

---

## 更新日志

详见 [CHANGELOG.md](CHANGELOG.md)。

---

## 贡献与反馈

欢迎提交 Pull Request 或 Issue。
