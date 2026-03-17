# Easy Immunofluorescence Analysis

Easy batch immunofluorescence (IF) image analysis tool with automatic channel detection, colocalization analysis, bleed-through correction, and statistical comparison.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## ✨ Features

- 🔬 **Automatic Channel Detection** - Automatic identification of DAPI/488/594/647 channels
- 🧬 **Nuclei Segmentation** - DAPI-based nuclei segmentation (watershed algorithm)
- 📊 **Intensity Analysis** - Signal intensity and positive area fraction analysis
- 🔗 **Colocalization Analysis** - Pearson correlation coefficient, Manders coefficients, Overlap coefficient
- 🎨 **Bleed-through Correction** - Multi-channel bleed-through correction (auto/manual coefficient estimation)
- 🖼️ **Flexible Visualization** - Pseudocolor/colormap/grayscale visualization modes
- 📈 **Statistical Comparison** - Two-group statistical comparison (Welch t-test / Mann-Whitney U)
- ✅ **QC Visualization** - Automatic generation of QC diagnostic and channel diagnostic images

---

## 📦 Installation

### Requirements

```bash
Python >= 3.8
```

### Install Dependencies

```bash
pip install numpy pandas scipy scikit-image opencv-python-headless matplotlib tifffile pillow
```

Or install from requirements.txt:
```bash
pip install -r requirements.txt
```

---

## 🚀 Quick Start

### 1. Prepare Your Images

Place your immunofluorescence images in two folders (e.g., experimental group and control group):

**Multi-channel Image Mode** (single file contains multiple channels):
```
data/
├── group1/          # Control group
│   ├── image1.tif   # 4-channel TIFF: DAPI/488/594/647
│   ├── image2.tif
│   └── ...
└── group2/          # Experimental group
    ├── image1.tif
    └── ...
```

**Sample Folder Mode** (each channel as separate file):
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

Supported formats: `.tif`, `.tiff`, `.png`, `.jpg`, `.jpeg`, `.bmp`

### 2. Configure the Script

Edit the **User Settings** section at the top of the script:

```python
# ============================================================================
# User Settings: edit this section
# ============================================================================

# Input/Output paths
GROUP1_DIR = r"/path/to/your/group1"
GROUP2_DIR = r"/path/to/your/group2"
GROUP1_NAME = "Control"
GROUP2_NAME = "Treatment"
OUTPUT_DIR = r"/path/to/output"

# Channel configuration
DAPI_CHANNEL = "DAPI"
ANALYSIS_CHANNELS = ["488", "594", "647"]

# Manually specify channel order (optional)
FILE_CHANNEL_ORDER = None  # Or ["647", "488", "DAPI"]

# Colocalization analysis
AUTO_ALL_PAIRWISE_STATS = True  # Automatically generate all pairwise comparisons
```

### 3. Run the Analysis

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

## ⚙️ Configuration Guide

### Channel Configuration

| Parameter | Description | Example |
|-----------|-------------|---------|
| `FILE_CHANNEL_ORDER` | Manually specify channel order | `["647", "488", "DAPI"]` |
| `DAPI_CHANNEL` | DAPI channel name | `"DAPI"` |
| `ANALYSIS_CHANNELS` | Channels to analyze | `["488", "594", "647"]` |
| `AUTO_ALL_PAIRWISE_STATS` | Automatic pairwise colocalization analysis | `True` / `False` |
| `MANUAL_COLOCALIZATION_PAIRS` | Manually specify colocalization pairs | `[("488", "594")]` |

### Threshold & Filtering Settings

```python
# Positive threshold rules (supports Otsu, Manual, Percentile methods)
POSITIVE_THRESHOLD_RULES = {
    "DAPI": {"method": "otsu", "scale": 1.0, "min_value": 0.0},
    "594": {"method": "otsu", "scale": 1.7, "min_value": 30.0},
    "488": {"method": "otsu", "scale": 1.1, "min_value": 30.0},
    "647": {"method": "otsu", "scale": 1.3, "min_value": 30.0},
}

# Minimum positive object area (to remove debris)
MIN_POSITIVE_OBJECT_AREA = 0  # Global default
POSITIVE_OBJECT_AREA_RULES = {  # Per-channel overrides
    "DAPI": 0,
    "488": 100,
    "594": 200,
    "647": 100,
}
```

### Background Subtraction

```python
# Background estimation percentile (global default)
BACKGROUND_PERCENTILE = 3.0

# Per-channel overrides (DAPI can be set lower)
BACKGROUND_PERCENTILE_RULES = {
    "DAPI": 1.0,
    "488": 3.0,
    "594": 3.0,
    "647": 3.0,
}
```

### Bleed-through Correction

```python
ENABLE_BLEEDTHROUGH_CORRECTION = True

# Define target <- source bleed-through relationships
BLEEDTHROUGH_SOURCE_MAP = {
    "488": DAPI_CHANNEL,  # Subtract DAPI from 488
    # "594": DAPI_CHANNEL,
    # "647": "488",
}

# Auto-estimation parameters
BLEEDTHROUGH_DEFAULT_RULE = {
    "mode": "auto",                    # "auto" or "manual"
    "estimate_mask": "nuclei",         # Estimation region: "nuclei"/"roi"/"all"
    "ratio_percentile": 20.0,          # Use lower percentile to avoid over-correction
    "source_threshold_percentile": 75.0,  # Only estimate on strong source pixels
    "min_pixels": 100,
    "max_coefficient": 3.0,
}

# Manual coefficient overrides
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    # "488": 0.15,
}
```

### Visualization Settings

```python
# Channel pseudocolors (RGB)
CHANNEL_COLORS = {
    "DAPI": (0, 0, 255),      # Blue
    "488": (0, 255, 0),       # Green
    "594": (255, 0, 0),       # Red
    "647": (255, 0, 178),     # Magenta
}

# Raw single-channel display: "color" (pseudocolor) or "gray" (grayscale)
RAW_CHANNEL_DISPLAY_MODE = "color"

# Positive preview display: "colormap" / "color" / "gray"
POSITIVE_PREVIEW_DISPLAY_MODE = "color"

# Colormap name (only used in colormap mode)
POSITIVE_PREVIEW_COLORMAP = "inferno"
```

### Nuclei Segmentation Parameters

```python
MIN_NUCLEUS_AREA = 0.5           # Minimum nucleus area (pixels)
GAUSSIAN_BLUR_SIZE = 7           # Gaussian blur kernel size
MIN_PEAK_DISTANCE = 1            # Minimum distance between watershed seeds
MASK_DILATION_RADIUS = 0         # ROI dilation radius relative to nuclei
```

> 💡 **Objective Lens Parameter Recommendations**:
> 
> | Parameter | 20X | 60X |
> |-----------|-----|-----|
> | `MIN_NUCLEUS_AREA` | 50 | 2200 |
> | `MIN_PEAK_DISTANCE` | 6 | 25 |
> | `MASK_DILATION_RADIUS` | 6 | 18 |
> | `BACKGROUND_PERCENTILE` | 3 | 5 |

---

## 📁 Output Structure

```
output/
└── Group1_vs_Group2_20260316_143022/
    ├── per_image_results.csv       # Detailed results per image
    ├── group_statistics.csv        # Inter-group statistical analysis
    ├── failed_images.csv           # Images that failed analysis
    ├── analysis_report.txt         # Text report (with configuration log)
    ├── run_config.json             # Runtime configuration (JSON format)
    ├── summary_plots.pdf           # Summary plots
    ├── summary_plots.png
    ├── qc_overlays/                # QC quality control images
    │   ├── image1_qc.png
    │   └── ...
    ├── channel_diagnostics/        # Channel diagnostic images
    │   ├── image1_channels.png
    │   └── ...
    ├── merged_denoised/            # Denoised merged images
    │   └── ...
    └── channel_diagnostics/bleedthrough/  # Bleed-through correction diagnostics
        └── ...
```

---

## 📊 Results Interpretation

### Per-Image Metrics

| Metric | Description |
|--------|-------------|
| `cell_count` | Nuclei count |
| `cell_density_per_mp` | Cell density (per million pixels) |
| `mean_nucleus_area_px` | Mean nucleus area (pixels) |
| `roi_area_fraction` | ROI area fraction relative to image |
| `{channel}_roi_mean_bgsub` | Mean intensity in ROI after background subtraction |
| `{channel}_roi_integrated_bgsub` | Total intensity in ROI after background subtraction |
| `{channel}_positive_area_fraction` | Positive area fraction (after threshold and area filtering) |
| `{channel}_norm_to_DAPI` | Normalized intensity relative to DAPI |
| `{channel}_integrated_per_cell` | Average intensity per cell |

### Colocalization Metrics

| Metric | Range | Description |
|--------|-------|-------------|
| `pearson_r` | -1 to 1 | Pearson correlation coefficient (linear correlation) |
| `manders_M1` | 0 to 1 | Manders M1 coefficient (fraction of Channel A overlapping with B) |
| `manders_M2` | 0 to 1 | Manders M2 coefficient (fraction of Channel B overlapping with A) |
| `overlap_coeff` | 0 to 1 | Overlap coefficient |
| `coloc_fraction` | 0 to 1 | Colocalized pixels as fraction of union of positive regions |

### Statistics

- **Normality Test**: Shapiro-Wilk (n >= 3)
- **Group Comparison**: 
  - Normal distribution: Welch t-test
  - Non-normal: Mann-Whitney U test
- **Multiple Testing Correction**: Benjamini-Hochberg FDR
- **Effect Size**: Hedges' g (t-test) / Rank-biserial r (Mann-Whitney)

---

## 🛠️ Advanced Usage

### Command Line Arguments

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
    --min-positive-object-area 200 \
    --positive-object-area-rules 488:100 594:200 \
    --gaussian-blur-size 7 \
    --min-peak-distance 25 \
    --mask-dilation-radius 18 \
    --raw-channel-display-mode color \
    --positive-preview-display-mode colormap \
    --positive-preview-colormap inferno \
    --no-qc \
    --no-channel-diagnostics \
    --no-bleedthrough-correction \
    --recursive-scan
```

### Input Discovery Modes

```python
INPUT_DISCOVERY_MODE = "auto"  # Auto-detect input structure
```

- `"auto"` - Automatic detection: uses sample folder mode if single-channel images detected, otherwise direct image mode
- `"sample_folders"` - Only analyze subfolders (each containing multiple single-channel files)
- `"image_files"` - Only analyze direct image files (multi-channel TIFF, etc.)

### Custom Channel Names

If your files use different channel naming (e.g., CH1, CH2):

```python
FILENAME_CHANNEL_MAP = {
    "CH1": "DAPI",
    "CH2": "488",
    "CH3": "594",
    "CH4": "647",
}
```

Ignore overlay/merge files:
```python
FILENAME_ROLE_MAP = {
    "OVERLAY": "overlay",
    "MERGE": "overlay",
    "RGB": "overlay",
}
```

Supported channel aliases:
- **DAPI**: DAPI, Hoechst, Hoechst33342, 405, Blue, Nuclei
- **488**: 488, FITC, GFP, EGFP, Green, Alexa488, AF488
- **594**: 594, TRITC, Cy3, TexasRed, TXRed, RFP, 561, 568, Alexa594, AF594
- **647**: 647, Cy5, APC, Alexa647, AF647, FarRed, 640

---

## 🐛 Troubleshooting

### Q: Incorrect channel detection
**A**: Manually set `FILE_CHANNEL_ORDER`:
```python
FILE_CHANNEL_ORDER = ["594", "488", "DAPI"]
```
Or check the diagnostic images in `channel_diagnostics/` to see actual channel order.

### Q: Inaccurate nuclei segmentation
**A**: Adjust segmentation parameters:
```python
MIN_NUCLEUS_AREA = 2200        # Increase to remove small noise
GAUSSIAN_BLUR_SIZE = 9         # Increase to smooth noise
MIN_PEAK_DISTANCE = 30         # Increase to prevent over-segmentation
```

### Q: High 594 false positives
**A**: Increase threshold and add minimum area filtering:
```python
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 2.0, "min_value": 50.0}
}
POSITIVE_OBJECT_AREA_RULES = {
    "594": 200  # Remove small debris
}
```

### Q: Over/under-correction of bleed-through
**A**: Switch to manual mode or adjust parameters:
```python
BLEEDTHROUGH_MANUAL_COEFFICIENTS = {
    "488": 0.12  # Manually specify coefficient
}
# Or adjust auto-estimation parameters
BLEEDTHROUGH_DEFAULT_RULE["ratio_percentile"] = 10.0  # More conservative estimation
```

### Q: Missing tifffile
**A**: 
```bash
pip install tifffile
```

---

## 📝 Version History

- **v1.0.6** - Added positive area filtering, per-channel background percentiles, visualization mode selection, improved bleed-through correction
- **v1.0.5** - Fixed channel detection, enhanced QC output, added bleed-through correction

---

## 📄 Citation

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

## 📜 License

MIT License - see [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## 📧 Contact

For questions or suggestions, please open an [Issue](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues).
