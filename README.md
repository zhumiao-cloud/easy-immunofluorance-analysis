# Easy Immunofluorescence Analysis

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/github/v/release/zhumiao-cloud/easy-immunofluorance-analysis)](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/releases)

> A streamlined Python tool for batch immunofluorescence (IF) image analysis with automated workflows, from raw images to publication-ready statistics.

**Features**: Automatic Channel Detection • Nuclei Segmentation • Colocalization Analysis • Bleed-through Correction • Statistical Comparison

---

## 🎯 Why This Tool?

| Feature | Easy IF | CellProfiler | ImageJ/FIJI |
|---------|---------|--------------|-------------|
| Setup time | **Minutes** | Hours | Hours |
| Batch processing | ✅ Native | ⚠️ Pipelines | ⚠️ Macros |
| Auto channel detection | ✅ Built-in | ❌ Manual | ❌ Manual |
| Colocalization | ✅ Automated | ⚠️ Plugins | ⚠️ Plugins |
| Bleed-through correction | ✅ Auto + Manual | ❌ | ❌ |
| Statistical output | ✅ Ready-to-publish | ⚠️ Export needed | ⚠️ Export needed |
| Command line interface | ✅ Yes | ⚠️ Limited | ❌ |

**Ideal for**: Researchers who need reproducible, quantitative IF analysis without spending days setting up complex pipelines.

---

## 🚀 Quick Start (5 minutes)

### 1. Install

```bash
# Clone the repository
git clone https://github.com/zhumiao-cloud/easy-immunofluorance-analysis.git
cd easy-immunofluorance-analysis

# Install dependencies
pip install -r requirements.txt
```

**Requirements**: Python 3.8+, 4GB+ RAM recommended for large images

### 2. Prepare Data

Organize your images in two folders (control vs. treatment):

```
experiment/
├── control/              # Group 1
│   ├── image1.tif        # 4-channel TIFF (DAPI/488/594/647)
│   └── image2.tif
└── treatment/            # Group 2
    ├── image1.tif
    └── image2.tif
```

Supports: `.tif`, `.tiff`, `.png`, `.jpg` (single or multi-channel)

### 3. Configure & Run

Edit the **User Settings** section in `if_analysis.py`:

```python
GROUP1_DIR = r"/path/to/experiment/control"
GROUP2_DIR = r"/path/to/experiment/treatment"
GROUP1_NAME = "Control"
GROUP2_NAME = "Treatment"
OUTPUT_DIR = r"/path/to/output"
```

Run:
```bash
python if_analysis.py
```

**Or use command line** (no editing needed):
```bash
python if_analysis.py \
    --group1-dir /path/to/control \
    --group2-dir /path/to/treatment \
    --group1-name Control \
    --group2-name Treatment \
    --intensity-channels 488 594 647
```

---

## 📊 What You Get

### Automated Outputs

```
output/
└── Control_vs_Treatment_20260317_143022/
    ├── per_image_results.csv       # All metrics per image
    ├── group_statistics.csv        # Statistical comparison
    ├── summary_plots.pdf           # Publication-ready figures
    ├── analysis_report.txt         # Complete analysis log
    ├── qc_overlays/                # Quality control images
    ├── channel_diagnostics/        # Channel detection verification
    └── merged_denoised/            # Denoised composite images
```

### Example Output Preview

<!-- Add example images here when available -->
<!-- ![Summary Plots](docs/images/summary_example.png) -->
<!-- ![QC Overlay](docs/images/qc_example.png) -->

---

## ⚙️ Configuration

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FILE_CHANNEL_ORDER` | Auto-detect | Force channel order: `["DAPI", "488", "594", "647"]` |
| `POSITIVE_THRESHOLD_RULES` | Otsu-based | Threshold method per channel |
| `MIN_NUCLEUS_AREA` | 0.5 px | Minimum nucleus size |
| `MASK_DILATION_RADIUS` | 0 px | Expand ROI beyond nuclei |
| `ENABLE_BLEEDTHROUGH_CORRECTION` | True | Correct spectral bleed-through |

### Microscope-Specific Presets

**20X Objective** (standard tissue sections):
```python
MIN_NUCLEUS_AREA = 50
MIN_PEAK_DISTANCE = 6
MASK_DILATION_RADIUS = 6
```

**60X/63X Objective** (high-resolution imaging):
```python
MIN_NUCLEUS_AREA = 2200
MIN_PEAK_DISTANCE = 25
MASK_DILATION_RADIUS = 18
BACKGROUND_PERCENTILE = 5
```

See [Configuration Guide](docs/CONFIGURATION.md) for detailed parameter documentation.

---

## 🔬 Analysis Capabilities

### 1. Intensity Quantification
- **Background subtraction**: Percentile-based estimation
- **Signal metrics**: Mean, integrated intensity per ROI/cell
- **Positive area**: Fraction of cells/ROI above threshold
- **Normalization**: To DAPI or total cell count

### 2. Colocalization Analysis
Automatically computed for all channel pairs:
- **Pearson r**: Linear correlation (-1 to 1)
- **Manders M1/M2**: Fraction of overlapping signal
- **Overlap coefficient**: Spatial overlap measure
- **Colocalization fraction**: Percentage of colocalized pixels

### 3. Bleed-through Correction
- **Auto-estimation**: Calculates correction coefficients from nuclei regions
- **Manual override**: Specify exact coefficients when needed
- **Validation**: Diagnostic plots show before/after correction

### 4. Statistical Comparison
Between-group comparison with automatic test selection:
- **Normality**: Shapiro-Wilk test
- **Comparison**: Welch t-test (normal) or Mann-Whitney U (non-normal)
- **Multiple testing**: Benjamini-Hochberg FDR correction
- **Effect size**: Hedges' g or Rank-biserial r

---

## 🐛 Common Issues & Solutions

<details>
<summary><strong>Incorrect channel detection</strong></summary>

Manually specify channel order:
```python
FILE_CHANNEL_ORDER = ["DAPI", "488", "594", "647"]
```
Check `channel_diagnostics/image_channels.png` to verify.
</details>

<details>
<summary><strong>Over/under-segmented nuclei</strong></summary>

Adjust based on your objective:
```python
# Too many small fragments → increase area
MIN_NUCLEUS_AREA = 2200    # For 60X

# Merged nuclei → increase peak distance
MIN_PEAK_DISTANCE = 30
```
</details>

<details>
<summary><strong>High false positives in marker channel</strong></summary>

Increase threshold and add area filtering:
```python
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 2.0, "min_value": 50.0}
}
POSITIVE_OBJECT_AREA_RULES = {"594": 200}  # Remove debris
```
</details>

<details>
<summary><strong>Excessive bleed-through correction</strong></summary>

Switch to manual or adjust estimation:
```python
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {"488": 0.12}
# OR more conservative auto-estimation:
BLEEDTHROUGH_DEFAULT_RULE["ratio_percentile"] = 10.0
```
</details>

See [Troubleshooting Guide](docs/TROUBLESHOOTING.md) for more solutions.

---

## 📖 Documentation

- [Installation Guide](docs/INSTALLATION.md)
- [Configuration Reference](docs/CONFIGURATION.md)
- [Output Format Specification](docs/OUTPUT_FORMAT.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [API Documentation](docs/API.md) (for programmatic use)

---

## 🧪 Citation

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

## 🗺️ Roadmap

- [ ] GUI interface for non-programmers
- [ ] Support for 3D/Z-stack analysis
- [ ] Plugin system for custom metrics
- [ ] Integration with OMERO/image databases
- [ ] Batch processing on HPC clusters

---

## 🤝 Contributing

Contributions welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

**Ways to contribute**:
- Report bugs or request features via [Issues](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues)
- Submit improvements via [Pull Requests](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/pulls)
- Share example data or use cases
- Improve documentation

---

## 📜 License

MIT License - see [LICENSE](LICENSE) for details.

---

## 💬 Support

- **Bug reports**: [Open an Issue](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues/new)
- **Questions**: [Discussions](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/discussions)
- **Email**: Open an issue for private inquiries

---

<p align="center">Made with 🔬 for the research community</p>
