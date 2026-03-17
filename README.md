# Easy Immunofluorescence Analysis | 简易免疫荧光分析工具

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/github/v/release/zhumiao-cloud/easy-immunofluorance-analysis)](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/releases)

> **English**: A streamlined Python tool for batch immunofluorescence (IF) image analysis with automated workflows, from raw images to publication-ready statistics.
> 
> **中文**：一个简化的 Python 工具，用于批量免疫荧光 (IF) 图像分析，实现从原始图像到可发表统计结果的自动化流程。

**Features | 功能特性**：Automatic Channel Detection (自动通道识别) • Nuclei Segmentation (细胞核分割) • Colocalization Analysis (共定位分析) • Bleed-through Correction (串色校正) • Statistical Comparison (统计比较)

---

## 🎯 Why This Tool? | 为什么选择这个工具？

| Feature | Easy IF | CellProfiler | ImageJ/FIJI |
|---------|---------|--------------|-------------|
| Setup time | **Minutes / 分钟** | Hours / 小时 | Hours / 小时 |
| Batch processing | ✅ Native / 原生支持 | ⚠️ Pipelines | ⚠️ Macros / 宏 |
| Auto channel detection | ✅ Built-in / 内置 | ❌ Manual / 手动 | ❌ Manual / 手动 |
| Colocalization | ✅ Automated / 自动 | ⚠️ Plugins / 插件 | ⚠️ Plugins / 插件 |
| Bleed-through correction | ✅ Auto + Manual / 自动+手动 | ❌ | ❌ |
| Statistical output | ✅ Ready-to-publish / 直接可用 | ⚠️ Export needed | ⚠️ Export needed |
| Command line interface | ✅ Yes / 支持 | ⚠️ Limited / 有限 | ❌ |

**Ideal for | 适用对象**：Researchers who need reproducible, quantitative IF analysis without spending days setting up complex pipelines.

需要可重复、定量 IF 分析的研究人员，无需花费数天设置复杂的分析流程。

---

## 🚀 Quick Start | 快速开始

### 1. Install | 安装

```bash
# Clone the repository / 克隆仓库
git clone https://github.com/zhumiao-cloud/easy-immunofluorance-analysis.git
cd easy-immunofluorance-analysis

# Install dependencies / 安装依赖
pip install -r requirements.txt
```

**Requirements | 环境要求**：Python 3.8+, 4GB+ RAM recommended for large images / 推荐 4GB+ 内存处理大图像

### 2. Prepare Data | 准备数据

Organize your images in two folders (control vs. treatment) / 将图像整理到两个文件夹中（对照组 vs 实验组）：

```
experiment/
├── control/              # Group 1 / 组1
│   ├── image1.tif        # 4-channel TIFF: DAPI/488/594/647
│   └── image2.tif
└── treatment/            # Group 2 / 组2
    ├── image1.tif
    └── image2.tif
```

**Supports | 支持格式**：`.tif`, `.tiff`, `.png`, `.jpg` (single or multi-channel / 单通道或多通道)

### 3. Configure & Run | 配置并运行

Edit the **User Settings** section in `if_analysis.py` / 编辑 `if_analysis.py` 中的 **用户设置** 部分：

```python
GROUP1_DIR = r"/path/to/experiment/control"
GROUP2_DIR = r"/path/to/experiment/treatment"
GROUP1_NAME = "Control"      # 对照组
GROUP2_NAME = "Treatment"    # 实验组
OUTPUT_DIR = r"/path/to/output"
```

**Run | 运行**：
```bash
python if_analysis.py
```

**Or use command line** (no editing needed / 无需编辑)：
```bash
python if_analysis.py \
    --group1-dir /path/to/control \
    --group2-dir /path/to/treatment \
    --group1-name Control \
    --group2-name Treatment \
    --intensity-channels 488 594 647
```

---

## 📊 What You Get | 输出结果

### Automated Outputs | 自动输出

```
output/
└── Control_vs_Treatment_20260317_143022/
    ├── per_image_results.csv       # All metrics per image / 每张图像的所有指标
    ├── group_statistics.csv        # Statistical comparison / 统计比较
    ├── summary_plots.pdf           # Publication-ready figures / 可直接发表的图表
    ├── analysis_report.txt         # Complete analysis log / 完整分析日志
    ├── qc_overlays/                # Quality control images / 质控图像
    ├── channel_diagnostics/        # Channel detection verification / 通道检测验证
    └── merged_denoised/            # Denoised composite images / 降噪合成图像
```

---

## ⚙️ Configuration | 配置

### Key Parameters | 关键参数

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FILE_CHANNEL_ORDER` | Auto-detect / 自动检测 | Force channel order: `["DAPI", "488", "594", "647"]` / 强制通道顺序 |
| `POSITIVE_THRESHOLD_RULES` | Otsu-based / Otsu算法 | Threshold method per channel / 每通道阈值方法 |
| `MIN_NUCLEUS_AREA` | 0.5 px | Minimum nucleus size / 最小细胞核大小 |
| `MASK_DILATION_RADIUS` | 0 px | Expand ROI beyond nuclei / ROI 相对细胞核扩展半径 |
| `ENABLE_BLEEDTHROUGH_CORRECTION` | True / 开启 | Correct spectral bleed-through / 校正光谱串色 |

### Microscope-Specific Presets | 显微镜专用预设

**20X Objective | 20X 物镜** (standard tissue sections / 标准组织切片)：
```python
MIN_NUCLEUS_AREA = 50
MIN_PEAK_DISTANCE = 6
MASK_DILATION_RADIUS = 6
```

**60X/63X Objective | 60X/63X 物镜** (high-resolution imaging / 高分辨率成像)：
```python
MIN_NUCLEUS_AREA = 2200
MIN_PEAK_DISTANCE = 25
MASK_DILATION_RADIUS = 18
BACKGROUND_PERCENTILE = 5
```

See [Configuration Guide](docs/CONFIGURATION.md) for detailed documentation.

详细文档请参见 [配置指南](docs/CONFIGURATION.md)。

---

## 🔬 Analysis Capabilities | 分析能力

### 1. Intensity Quantification | 强度定量
- **Background subtraction / 背景扣除**: Percentile-based estimation / 基于百分位数估计
- **Signal metrics / 信号指标**: Mean, integrated intensity per ROI/cell / ROI/细胞的平均、总强度
- **Positive area / 阳性区域**: Fraction of cells/ROI above threshold / 高于阈值的细胞/ROI比例
- **Normalization / 归一化**: To DAPI or total cell count / 相对于 DAPI 或总细胞数

### 2. Colocalization Analysis | 共定位分析
Automatically computed for all channel pairs / 自动计算所有通道对：
- **Pearson r**: Linear correlation (-1 to 1) / 线性相关系数
- **Manders M1/M2**: Fraction of overlapping signal / 重叠信号比例
- **Overlap coefficient**: Spatial overlap measure / 空间重叠度量
- **Colocalization fraction**: Percentage of colocalized pixels / 共定位像素百分比

### 3. Bleed-through Correction | 串色校正
- **Auto-estimation / 自动估计**: Calculates correction coefficients from nuclei regions / 从细胞核区域计算校正系数
- **Manual override / 手动覆盖**: Specify exact coefficients when needed / 需要时指定精确系数
- **Validation / 验证**: Diagnostic plots show before/after correction / 诊断图显示校正前后

### 4. Statistical Comparison | 统计比较
Between-group comparison with automatic test selection / 组间比较，自动选择检验方法：
- **Normality / 正态性**: Shapiro-Wilk test / Shapiro-Wilk 检验
- **Comparison / 比较**: Welch t-test (normal / 正态) or Mann-Whitney U (non-normal / 非正态)
- **Multiple testing / 多重检验**: Benjamini-Hochberg FDR correction / Benjamini-Hochberg FDR 校正
- **Effect size / 效应量**: Hedges' g or Rank-biserial r

---

## 🐛 Common Issues & Solutions | 常见问题与解决方案

<details>
<summary><strong>Q: Incorrect channel detection / 通道识别错误</strong></summary>

Manually specify channel order / 手动指定通道顺序：
```python
FILE_CHANNEL_ORDER = ["DAPI", "488", "594", "647"]
```
Check `channel_diagnostics/image_channels.png` to verify / 检查验证。
</details>

<details>
<summary><strong>Q: Over/under-segmented nuclei / 细胞核分割过度/不足</strong></summary>

Adjust based on your objective / 根据物镜调整：
```python
# Too many small fragments / 碎片过多 → increase area / 增大面积
MIN_NUCLEUS_AREA = 2200    # For 60X

# Merged nuclei / 细胞核合并 → increase peak distance / 增大峰值距离
MIN_PEAK_DISTANCE = 30
```
</details>

<details>
<summary><strong>Q: High false positives in marker channel / 标记通道假阳性过高</strong></summary>

Increase threshold and add area filtering / 提高阈值并增加面积过滤：
```python
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 2.0, "min_value": 50.0}
}
POSITIVE_OBJECT_AREA_RULES = {"594": 200}  # Remove debris / 去除碎点
```
</details>

<details>
<summary><strong>Q: Excessive bleed-through correction / 串色校正过度</strong></summary>

Switch to manual or adjust estimation / 切换到手动或调整估计：
```python
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {"488": 0.12}
# OR / 或者 more conservative auto-estimation / 更保守的自动估计：
BLEEDTHROUGH_DEFAULT_RULE["ratio_percentile"] = 10.0
```
</details>

See [Troubleshooting Guide](docs/TROUBLESHOOTING.md) for more solutions.

更多解决方案请参见 [故障排除指南](docs/TROUBLESHOOTING.md)。

---

## 📖 Documentation | 文档

- [Installation Guide | 安装指南](docs/INSTALLATION.md)
- [Configuration Reference | 配置参考](docs/CONFIGURATION.md)
- [Output Format Specification | 输出格式说明](docs/OUTPUT_FORMAT.md)
- [Troubleshooting | 故障排除](docs/TROUBLESHOOTING.md)
- [API Documentation | API 文档](docs/API.md) (for programmatic use / 程序化使用)

---

## 🧪 Citation | 引用

If you use this tool in your research, please cite:

如果在研究中使用此工具，请引用：

```bibtex
@software{easy_immunofluorescence_analysis,
  title = {Easy Immunofluorescence Analysis: Batch Immunofluorescence Image Analysis Tool},
  author = {zhumiao-cloud},
  year = {2026},
  url = {https://github.com/zhumiao-cloud/easy-immunofluorance-analysis}
}
```

---

## 🗺️ Roadmap | 路线图

- [ ] GUI interface for non-programmers / 非程序员图形界面
- [ ] Support for 3D/Z-stack analysis / 3D/Z-stack 分析支持
- [ ] Plugin system for custom metrics / 自定义指标插件系统
- [ ] Integration with OMERO/image databases / OMERO/图像数据库集成
- [ ] Batch processing on HPC clusters / HPC 集群批量处理

---

## 🤝 Contributing | 贡献

Contributions welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

欢迎贡献！请参见 [CONTRIBUTING.md](CONTRIBUTING.md) 了解指南。

**Ways to contribute | 贡献方式**：
- Report bugs or request features via [Issues](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues)
- Submit improvements via [Pull Requests](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/pulls)
- Share example data or use cases / 分享示例数据或用例
- Improve documentation / 改进文档

---

## 📜 License | 许可证

MIT License - see [LICENSE](LICENSE) for details.

MIT 许可证 - 详情参见 [LICENSE](LICENSE)。

---

## 💬 Support | 支持

- **Bug reports / 错误报告**: [Open an Issue](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues/new)
- **Questions / 问题**: [Discussions](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/discussions)
- **Email / 邮件**: Open an issue for private inquiries / 如需私人咨询请开 Issue

---

<p align="center">Made with 🔬 for the research community | 为科研社区打造</p>
