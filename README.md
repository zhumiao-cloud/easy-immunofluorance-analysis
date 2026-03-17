# Easy Immunofluorescence Analysis

Easy batch immunofluorescence (IF) image analysis tool with automatic channel detection, colocalization analysis, bleed-through correction, and statistical comparison.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/github/v/release/zhumiao-cloud/easy-immunofluorance-analysis)](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/releases)

---

## Features

- Automatic Channel Detection - Auto-detect DAPI/488/594/647 channels from metadata or heuristics
- Nuclei Segmentation - Watershed-based nuclei segmentation with configurable parameters
- Intensity Analysis - Signal intensity and positive area fraction quantification
- Colocalization Analysis - Pearson correlation, Manders' coefficients, overlap coefficient
- Bleed-through Correction - Automatic or manual coefficient estimation for crosstalk removal
- Flexible Input Modes - Multi-channel images or sample folders with single-channel files
- Statistical Comparison - Two-group statistical testing (Welch t-test / Mann-Whitney U)
- QC Visualization - Automated QC overlays and channel diagnostics

---

## Installation

### Requirements

```bash
Python >= 3.8
```

### Dependencies

```bash
pip install numpy pandas scipy scikit-image opencv-python-headless matplotlib tifffile pillow
```

Or install from requirements.txt:
```bash
pip install -r requirements.txt
```

---

## Quick Start

### 1. Prepare Your Images

Organize images in two folders (e.g., control vs. treatment):

**Multi-channel mode** (single file with all channels):
```
data/
├── group1/
│   ├── image1.tif   # 4-channel TIFF: DAPI/488/594/647
│   └── image2.tif
└── group2/
    ├── image1.tif
    └── ...
```

**Sample folder mode** (one file per channel):
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

Edit the **User Settings** section at the top of `if_analysis.py`:

```python
GROUP1_DIR = r"/path/to/your/group1"
GROUP2_DIR = r"/path/to/your/group2"
GROUP1_NAME = "Control"
GROUP2_NAME = "Treatment"
OUTPUT_DIR = r"/path/to/output"

DAPI_CHANNEL = "DAPI"
ANALYSIS_CHANNELS = ["488", "594", "647"]

# Manual channel order (optional, for overriding auto-detection)
FILE_CHANNEL_ORDER = None  # or ["647", "488", "DAPI"]
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

## Configuration Guide

### Channel Configuration

| Parameter | Description | Example |
|-----------|-------------|---------|
| `FILE_CHANNEL_ORDER` | Manual channel order override | `["647", "488", "DAPI"]` |
| `DAPI_CHANNEL` | Nuclear channel name | `"DAPI"` |
| `ANALYSIS_CHANNELS` | Channels to analyze | `["488", "594", "647"]` |
| `AUTO_ALL_PAIRWISE_STATS` | Auto-generate all pairwise coloc | `True` / `False` |

### Threshold Settings

```python
# Positive threshold rules (Otsu, Manual, or Percentile methods)
CHANNEL_SETTINGS = {
    "DAPI": {"positive_threshold": {"method": "otsu", "scale": 1.0, "min_value": 0.0}},
    "594": {"positive_threshold": {"method": "otsu", "scale": 1.7, "min_value": 30.0}},
    "488": {"positive_threshold": {"method": "otsu", "scale": 1.1, "min_value": 30.0}},
    "647": {"positive_threshold": {"method": "otsu", "scale": 1.3, "min_value": 30.0}},
}
```

Available methods:
- `otsu`: Otsu auto-threshold × scale factor
- `manual`: Fixed threshold value
- `percentile`: Percentile-based threshold

### Bleed-through Correction

```python
ENABLE_BLEEDTHROUGH_CORRECTION = True

BLEEDTHROUGH_SOURCE_MAP = {
    "488": "DAPI",  # Correct 488 bleed-through from DAPI
}

BLEEDTHROUGH_DEFAULT_RULE = {
    "mode": "auto",              # "auto" or "manual"
    "estimate_mask": "nuclei",   # "nuclei", "roi", or "all"
    "ratio_percentile": 20.0,
    "source_threshold_percentile": 75.0,
}
```

### Nucleus Segmentation Parameters

```python
NUCLEUS_SEGMENTATION_SETTINGS = {
    "background_percentile": 0.0,    # Background subtraction percentile
    "gaussian_blur_size": 0,         # Gaussian blur kernel (0=off)
    "threshold": {"method": "manual", "value": 0.0},
    "min_mask_object_area": 0,       # Remove small objects
    "min_mask_hole_area": 0,         # Fill small holes
    "opening_radius": 0,             # Morphological opening
    "closing_radius": 0,             # Morphological closing
    "min_peak_distance": 0,          # Watershed seed distance
    "min_nucleus_area": 0,           # Final area filter
    "mask_dilation_radius": 0,       # ROI dilation radius
}
```

> **Recommended for 60X objective**:
> - `min_nucleus_area = 2200`
> - `min_peak_distance = 25`
> - `mask_dilation_radius = 18`

---

## Output Structure

```
output/
└── Group1_vs_Group2_20260318_143022/
    ├── per_image_results.csv       # Per-image detailed metrics
    ├── group_statistics.csv        # Inter-group statistical analysis
    ├── failed_images.csv           # Failed image log
    ├── analysis_report.txt         # Text summary report
    ├── run_config.json             # Configuration snapshot
    ├── summary_plots.pdf           # Summary visualizations
    ├── summary_plots.png
    ├── qc_overlays/                # QC overlay images
    │   ├── image1_qc.png
    │   └── ...
    └── channel_diagnostics/        # Channel diagnostic images
        ├── image1_channels.png
        └── ...
```

---

## Results Interpretation

### Per-Image Metrics

| Metric | Description |
|--------|-------------|
| `cell_count` | Number of detected nuclei |
| `cell_density_per_mp` | Cell density (per million pixels) |
| `mean_nucleus_area_px` | Mean nucleus area (pixels) |
| `{channel}_roi_mean_bgsub` | ROI mean intensity (background-subtracted) |
| `{channel}_roi_integrated_bgsub` | ROI integrated intensity (background-subtracted) |
| `{channel}_positive_area_fraction` | Positive signal area fraction |
| `{channel}_norm_to_DAPI` | Intensity normalized to DAPI |

### Colocalization Metrics

| Metric | Range | Description |
|--------|-------|-------------|
| `pearson_r` | -1 to 1 | Pearson correlation coefficient |
| `manders_{A}` | 0 to 1 | Manders' M1 (fraction of A overlapping B) |
| `manders_{B}` | 0 to 1 | Manders' M2 (fraction of B overlapping A) |
| `overlap_coeff` | 0 to 1 | Overlap coefficient |
| `coloc_fraction` | 0 to 1 | Colocalized pixel fraction |

---

## Command Line Interface

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
    --recursive-scan \
    --filename-channel-map "CH1:DAPI,CH2:488,CH3:594"
```

---

## Troubleshooting

### Channel Detection Errors
Set manual channel order:
```python
FILE_CHANNEL_ORDER = ["594", "488", "DAPI"]
```

### Inaccurate Nuclei Segmentation
Adjust segmentation parameters:
```python
NUCLEUS_SEGMENTATION_SETTINGS = {
    "min_nucleus_area": 2200,      # Increase to remove debris
    "gaussian_blur_size": 9,       # Increase to smooth noise
    "min_peak_distance": 30,       # Increase to prevent oversegmentation
}
```

### High False Positives
Increase the threshold:
```python
CHANNEL_SETTINGS = {
    "594": {"positive_threshold": {"method": "otsu", "scale": 2.0, "min_value": 50.0}}
}
```

### Missing tifffile
```bash
pip install tifffile
```

---

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{easy_immunofluorescence_analysis,
  title = {Easy Immunofluorescence Analysis: Batch IF Image Analysis Tool},
  author = {zhumiao-cloud},
  year = {2026},
  url = {https://github.com/zhumiao-cloud/easy-immunofluorance-analysis}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## Contact

For questions or suggestions, please open an [Issue](https://github.com/zhumiao-cloud/easy-immunofluorance-analysis/issues).
