#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch immunofluorescence analysis for direct multichannel images and
sample-folder inputs made of single-channel files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import warnings
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from datetime import datetime
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage as ndi, stats
from skimage import feature, filters, measure, morphology, segmentation

try:
    import tifffile
except Exception:
    tifffile = None

try:
    from PIL import Image, ImageDraw, ImageFont
except Exception:
    Image = None
    ImageDraw = None
    ImageFont = None


# ============================================================================
# User Settings: edit this section
# ============================================================================

# Paths
GROUP1_DIR = r"/Users/yemingzhu/Downloads/课题文件/原始数据胸腺/免疫荧光/all_in_one/all_in_one_analysis/cd11c_cd31_cxcr4/3m"  # 组1输入目录
GROUP2_DIR = r"/Users/yemingzhu/Downloads/课题文件/原始数据胸腺/免疫荧光/all_in_one/all_in_one_analysis/cd11c_cd31_cxcr4/22m"  # 组2输入目录
GROUP1_NAME = "3m"  # 图表和结果表里的组1名称
GROUP2_NAME = "22m"  # 图表和结果表里的组2名称
OUTPUT_DIR = r"/Users/yemingzhu/Downloads/课题文件/原始数据胸腺/免疫荧光/if_analysis_results"  # 输出总目录

# Input structure
INPUT_DISCOVERY_MODE = "auto"  # auto / sample_folders / image_files
IMAGE_PATTERNS = ["*.tif", "*.tiff", "*.png", "*.jpg", "*.jpeg", "*.bmp"]  # 会被扫描的图像后缀
RECURSIVE_SCAN = False  # True 时递归扫描子目录

AUTO_DETECT_CHANNELS = True  # 是否尝试自动识别通道顺序
STRICT_CHANNEL_DETECTION = True  # True 时检测不明确就报错，不静默兜底
FILE_CHANNEL_ORDER: Optional[List[str]] = None  # 单个多通道文件的手动通道顺序；None 表示自动判断

FILENAME_ROLE_MAP = {  # 用于在样本文件夹模式下忽略 overlay / merge 一类文件
    "OVERLAY": "overlay",
    "MERGE": "overlay",
    "MERGED": "overlay",
    "COMPOSITE": "overlay",
    "RGB": "overlay",
}

FILENAME_CHANNEL_MAP = {  # 当文件名不直接写 DAPI/488/594/647 时使用
    "CH1": "DAPI",
    "CH2": "488",
    "CH3": "594",
    "CH4": "647",
}

# Channels and statistics
DAPI_CHANNEL = "DAPI"  # 核通道名称
ANALYSIS_CHANNELS = ["488", "594", "647"]  # 进入强度和共定位分析的通道
AUTO_ALL_PAIRWISE_STATS = True  # True 时自动做 ANALYSIS_CHANNELS 的两两配对
MANUAL_COLOCALIZATION_PAIRS: List[Tuple[str, str]] = []  # 仅在 AUTO_ALL_PAIRWISE_STATS=False 时使用

# Quantification
BACKGROUND_PERCENTILE = 3.0  # 背景估计分位数的全局默认值
BACKGROUND_PERCENTILE_RULES: Dict[str, float] = {  # 按通道覆盖背景分位数；DAPI 可单独设低一点
    "DAPI": 1.0,
    "488": 3.0,
    "594": 3.0,
    "647": 3.0,
}
MIN_NUCLEUS_AREA = 0.5  # 仅用于 DAPI 核分割；运行时会转成 >=1 的整数
MIN_POSITIVE_OBJECT_AREA = 0  # 阳性阈值后碎点过滤的全局默认面积；0 表示关闭
POSITIVE_OBJECT_AREA_RULES: Dict[str, int] = {  # 按通道覆盖碎点过滤面积；DAPI 可与 marker 分开设置
    "DAPI": 0,
    "488": 100,
    "594": 200,
    "647": 100,
}
GAUSSIAN_BLUR_SIZE = 7  # DAPI 分割前的高斯平滑核大小
MIN_PEAK_DISTANCE = 1  # watershed 种子点最小距离
MASK_DILATION_RADIUS = 0  # ROI 相对 nuclei 向外扩张的像素半径

POSITIVE_THRESHOLD_RULES = {  # 各通道阳性阈值规则；若包含 DAPI，也会影响 DAPI 的 split/merged positive preview
    "DAPI": {"method": "otsu", "scale": 1.0, "min_value": 0.0},
    "594": {"method": "otsu", "scale": 1.7, "min_value": 30.0},
    "488": {"method": "otsu", "scale": 1.1, "min_value": 30.0},
    "647": {"method": "otsu", "scale": 1.3, "min_value": 30.0},
}

# Bleed-through correction
ENABLE_BLEEDTHROUGH_CORRECTION = True  # 是否启用串色扣除

BLEEDTHROUGH_SOURCE_MAP = {  # target <- source
    "488": DAPI_CHANNEL,
    # "594": DAPI_CHANNEL,
    # "647": "488",
}

BLEEDTHROUGH_DEFAULT_RULE: Dict[str, Any] = {  # 自动估计串色系数时的默认参数
    "mode": "auto",
    "estimate_mask": "nuclei",
    "ratio_percentile": 20.0,
    "source_threshold_percentile": 75.0,
    "min_pixels": 100,
    "max_coefficient": 3.0,
}

BLEEDTHROUGH_MANUAL_COEFFICIENTS: Dict[str, float] = {}  # 需要手动系数时在这里填 target:coefficient
BLEEDTHROUGH_RULE_OVERRIDES: Dict[str, Dict[str, Any]] = {}  # 只在少数通道要覆盖默认规则时填写

# Visualization and outputs
CHANNEL_COLORS: Dict[str, Tuple[int, int, int]] = {  # split preview 与 merged 共用这套颜色
    "DAPI": (0, 0, 255),
    "488": (0, 255, 0),
    "594": (255, 0, 0),
    "647": (255, 0, 178),
}
RAW_CHANNEL_DISPLAY_MODE = "color"           # 原始单通道面板显示方式: color / gray
POSITIVE_PREVIEW_DISPLAY_MODE = "color"      # 阳性预览显示方式: colormap / color / gray
POSITIVE_PREVIEW_COLORMAP = "inferno"        # 仅在 POSITIVE_PREVIEW_DISPLAY_MODE="colormap" 时使用
SAVE_CHANNEL_DIAGNOSTICS = True              # 是否保存 channel_diagnostics
SAVE_QC_OVERLAYS = True                      # 是否保存 qc_overlays
SAVE_BLEEDTHROUGH_DIAGNOSTICS = True         # 是否保存串色校正诊断图
SAVE_FIGURE_PDF = True                       # 是否保存汇总 PDF
SAVE_FIGURE_PNG = True                       # 是否保存汇总 PNG
FIGURE_DPI = 500                             # 输出图分辨率

# Derived from the settings above. Usually do not edit.
INTENSITY_CHANNELS = list(ANALYSIS_CHANNELS)
AUTO_COLOCALIZATION_PAIRS = bool(AUTO_ALL_PAIRWISE_STATS)
COLOCALIZATION_PAIRS = list(MANUAL_COLOCALIZATION_PAIRS)
FALLBACK_CHANNEL_ORDER = [DAPI_CHANNEL] + [x for x in INTENSITY_CHANNELS if x != DAPI_CHANNEL]


# ============================================================================
# Internal constants: usually do not edit
# ============================================================================

SIGNIFICANCE_LEVELS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]
DEFAULT_GROUP_COLORS = ("#4472C4", "#ED7D31","#A747BE")
MIN_PIXELS_FOR_COLOC = 10
ALLOWED_CHANNEL_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_\-+.]+$")
CANONICAL_CHANNELS = ("DAPI", "488", "594", "647")
CHANNEL_PSEUDOCOLOR_RGB = {
    "DAPI": (0, 0, 255),
    "488": (0, 255, 0),
    "594": (255, 0, 0),
    "647": (255, 0, 178),
}


# ============================================================================
# 数据结构
# ============================================================================

@dataclass
class AnalysisConfig:
    group1_dir: Path
    group2_dir: Path
    group1_name: str
    group2_name: str
    output_dir: Path

    input_discovery_mode: str
    file_channel_order: Optional[List[str]]
    auto_detect_channels: bool
    strict_channel_detection: bool
    save_channel_diagnostics: bool
    fallback_channel_order: List[str]
    positive_threshold_rules: Dict[str, Dict[str, Any]]

    enable_bleedthrough_correction: bool
    save_bleedthrough_diagnostics: bool
    bleedthrough_rules: Dict[str, Dict[str, Any]]

    dapi_channel: str
    intensity_channels: List[str]
    auto_colocalization_pairs: bool
    colocalization_pairs: List[Tuple[str, str]]
    filename_role_map: Dict[str, str]
    filename_channel_map: Dict[str, str]

    image_patterns: List[str] = field(default_factory=lambda: ["*.tif", "*.tiff", "*.png", "*.jpg", "*.jpeg", "*.bmp"])
    recursive_scan: bool = False
    background_percentile: float = 10.0
    background_percentile_rules: Dict[str, float] = field(default_factory=dict)
    min_nucleus_area: int = 50
    min_positive_object_area: int = 0
    positive_object_area_rules: Dict[str, int] = field(default_factory=dict)
    gaussian_blur_size: int = 5
    min_peak_distance: int = 6
    mask_dilation_radius: int = 6
    channel_colors: Dict[str, Tuple[int, int, int]] = field(default_factory=dict)
    raw_channel_display_mode: str = "gray"
    positive_preview_display_mode: str = "colormap"
    positive_preview_colormap: str = "inferno"
    save_qc_overlays: bool = True
    save_figure_pdf: bool = True
    save_figure_png: bool = True
    figure_dpi: int = 220


@dataclass
class ImageLoadResult:
    image: np.ndarray
    source_kind: str
    axes: str
    photometric: str
    raw_channel_names: List[str]
    notes: List[str] = field(default_factory=list)


@dataclass
class ChannelResolution:
    channel_order: List[str]
    channel_to_index: Dict[str, int]
    method: str
    notes: List[str]


@dataclass
class BleedthroughCorrectionInfo:
    applied: bool
    target_channel: str
    source_channel: str
    mode: str
    coefficient: float
    note: str = ""
    source_background: float = float("nan")


@dataclass
class PositiveThresholdInfo:
    channel_name: str
    threshold: float
    base_threshold: float
    method: str
    note: str = ""


@dataclass
class ImageAnalysisResult:
    values: Dict[str, Any]
    roi_mask: np.ndarray
    nuclei_labels: np.ndarray


@dataclass
class SampleFolderContents:
    sample_dir: Path
    channel_files: Dict[str, Path]
    duplicate_channel_files: Dict[str, List[Path]] = field(default_factory=dict)
    overlay_files: List[Path] = field(default_factory=list)
    unclassified_files: List[Path] = field(default_factory=list)


@dataclass
class AnalysisInput:
    path: Path
    source_mode: str
    sample_contents: Optional[SampleFolderContents] = None


# ============================================================================
# 工具函数
# ============================================================================


def normalize_channel_name(name: str) -> str:
    value = str(name).strip()
    if not value:
        raise ValueError("通道名称不能为空。")
    return value.upper()


def simplify_token(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", str(text).upper())


def canonicalize_channel_label(name: str) -> Optional[str]:
    token = simplify_token(name)
    if not token:
        return None

    direct_map = {
        "DAPI": "DAPI",
        "HOECHST": "DAPI",
        "HOECHST33342": "DAPI",
        "H33342": "DAPI",
        "NUC": "DAPI",
        "NUCLEI": "DAPI",
        "NUCLEUS": "DAPI",
        "BLUE": "DAPI",
        "405": "DAPI",
        "488": "488",
        "AF488": "488",
        "ALEXA488": "488",
        "ALEXAFLUOR488": "488",
        "FITC": "488",
        "GFP": "488",
        "EGFP": "488",
        "GREEN": "488",
        "594": "594",
        "AF594": "594",
        "ALEXA594": "594",
        "ALEXAFLUOR594": "594",
        "TRITC": "594",
        "CY3": "594",
        "TEXASRED": "594",
        "TXRED": "594",
        "RFP": "594",
        "561": "594",
        "568": "594",
        "647": "647",
        "AF647": "647",
        "ALEXA647": "647",
        "ALEXAFLUOR647": "647",
        "CY5": "647",
        "FARRED": "647",
        "APC": "647",
        "640": "647",
    }
    if token in direct_map:
        return direct_map[token]

    if "DAPI" in token or "HOECHST" in token or "NUC" in token or "405" in token:
        return "DAPI"
    if "488" in token or "FITC" in token or "GFP" in token or "GREEN" in token:
        return "488"
    if "594" in token or "561" in token or "568" in token or "CY3" in token or "TRITC" in token or "TEXASRED" in token or "RFP" in token:
        return "594"
    if "647" in token or "640" in token or "CY5" in token or "FARRED" in token or "APC" in token:
        return "647"
    return None


def ensure_odd(value: int) -> int:
    value = int(max(1, value))
    return value if value % 2 == 1 else value + 1


def natural_sort_key(value: Any) -> List[Any]:
    text = str(value)
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def unique_in_order(items: Iterable[Any]) -> List[Any]:
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def sanitize_filename(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    text = text.strip("._")
    return text or "output"


def build_output_basename(group_name: str, image_path: Path) -> str:
    group_token = sanitize_filename(group_name)
    image_token = sanitize_filename(image_path.stem)
    digest = hashlib.md5(str(image_path.resolve()).encode("utf-8")).hexdigest()[:8]
    return f"{group_token}__{image_token}__{digest}"


def normalize_filename_role_map(raw_map: Dict[str, str]) -> Dict[str, str]:
    normalized: Dict[str, str] = {}
    role_aliases = {
        "OVERLAY": "overlay",
        "MERGE": "overlay",
        "MERGED": "overlay",
        "COMPOSITE": "overlay",
        "RGB": "overlay",
        "IGNORE": "ignore",
        "IGNORED": "ignore",
        "SKIP": "ignore",
        "EXCLUDE": "ignore",
    }
    for raw_token, raw_role in (raw_map or {}).items():
        token = simplify_token(raw_token)
        if not token:
            raise ValueError("FILENAME_ROLE_MAP 中存在空关键字。")
        role_key = simplify_token(raw_role)
        if role_key not in role_aliases:
            raise ValueError(
                f"FILENAME_ROLE_MAP 中存在不支持的角色: {raw_role}。"
                " 可选值: overlay / ignore"
            )
        normalized[token] = role_aliases[role_key]
    return normalized


def match_filename_role(text: str, filename_role_map: Dict[str, str]) -> Optional[str]:
    raw_text = str(text).upper()
    simplified_text = simplify_token(text)
    for token, role_name in sorted(filename_role_map.items(), key=lambda item: (-len(item[0]), item[0])):
        if re.search(rf"(?<![A-Z0-9]){re.escape(token)}(?![A-Z0-9])", raw_text):
            return role_name
        if simplified_text == token or simplified_text.startswith(token) or simplified_text.endswith(token):
            return role_name
    return None


def infer_file_role_from_path(path: Path, filename_role_map: Optional[Dict[str, str]] = None) -> Optional[str]:
    candidates = [
        path.stem,
        path.name,
        path.stem.replace("_", " "),
        path.stem.replace("-", " "),
    ]
    for candidate in candidates:
        if filename_role_map:
            role = match_filename_role(candidate, filename_role_map)
            if role is not None:
                return role
    return None


def is_overlay_like_path(path: Path, filename_role_map: Optional[Dict[str, str]] = None) -> bool:
    role = infer_file_role_from_path(path, filename_role_map=filename_role_map)
    if role in {"overlay", "ignore"}:
        return True
    token = simplify_token(path.stem)
    overlay_tokens = ("OVERLAY", "MERGE", "MERGED", "COMPOSITE", "RGB")
    return any(item in token for item in overlay_tokens)


def normalize_filename_channel_map(raw_map: Dict[str, str]) -> Dict[str, str]:
    normalized: Dict[str, str] = {}
    for raw_token, raw_channel in raw_map.items():
        token = simplify_token(raw_token)
        if not token:
            raise ValueError("FILENAME_CHANNEL_MAP 中存在空关键字。")
        channel_name = normalize_channel_name(raw_channel)
        if channel_name not in CANONICAL_CHANNELS:
            raise ValueError(f"FILENAME_CHANNEL_MAP 中存在不支持的通道名: {raw_channel}")
        normalized[token] = channel_name
    return normalized


def match_filename_channel(text: str, filename_channel_map: Dict[str, str]) -> Optional[str]:
    raw_text = str(text).upper()
    simplified_text = simplify_token(text)
    for token, channel_name in sorted(filename_channel_map.items(), key=lambda item: (-len(item[0]), item[0])):
        if re.search(rf"(?<![A-Z0-9]){re.escape(token)}(?![A-Z0-9])", raw_text):
            return channel_name
        if simplified_text == token or simplified_text.startswith(token) or simplified_text.endswith(token):
            return channel_name
    return None


def infer_channel_from_path(path: Path, filename_channel_map: Optional[Dict[str, str]] = None) -> Optional[str]:
    candidates = [
        path.stem,
        path.name,
        path.stem.replace("_", " "),
        path.stem.replace("-", " "),
    ]
    for candidate in candidates:
        if filename_channel_map:
            mapped = match_filename_channel(candidate, filename_channel_map)
            if mapped is not None:
                return mapped
        canonical = canonicalize_channel_label(candidate)
        if canonical is not None:
            return canonical
    return None


def build_all_channel_pairs(channel_names: Sequence[str]) -> List[Tuple[str, str]]:
    normalized = [normalize_channel_name(name) for name in channel_names]
    return [(a, b) for a, b in combinations(normalized, 2)]


def normalize_colocalization_pairs(raw_pairs: Sequence[Tuple[str, str]]) -> List[Tuple[str, str]]:
    normalized: List[Tuple[str, str]] = []
    seen = set()

    for channel_a, channel_b in raw_pairs:
        a_name = normalize_channel_name(channel_a)
        b_name = normalize_channel_name(channel_b)
        if a_name == b_name:
            raise ValueError(f"共定位通道对不能是同一个通道: {a_name}")
        key = tuple(sorted((a_name, b_name)))
        if key in seen:
            continue
        seen.add(key)
        normalized.append((a_name, b_name))

    return normalized


def clean_numeric(values: Sequence[Any]) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64).ravel()
    if arr.size == 0:
        return np.asarray([], dtype=np.float64)
    return arr[np.isfinite(arr)]


def sample_std(values: Sequence[Any]) -> float:
    arr = clean_numeric(values)
    if arr.size < 2:
        return float("nan")
    return float(np.std(arr, ddof=1))


def sem(values: Sequence[Any]) -> float:
    arr = clean_numeric(values)
    if arr.size < 2:
        return float("nan")
    return float(sample_std(arr) / math.sqrt(arr.size))


def iqr(values: Sequence[Any]) -> float:
    arr = clean_numeric(values)
    if arr.size == 0:
        return float("nan")
    q1, q3 = np.percentile(arr, [25, 75])
    return float(q3 - q1)


def significance_label(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "NA"
    for threshold, label in SIGNIFICANCE_LEVELS:
        if p_value < threshold:
            return label
    return "ns"


def format_p_value(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "NA"
    if p_value < 1e-4:
        return f"{p_value:.2e}"
    return f"{p_value:.4f}"


def make_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def array_to_uint8(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float64)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.zeros(arr.shape, dtype=np.uint8)
    arr = arr.copy()
    arr[~finite] = 0.0
    min_value = float(np.min(arr))
    max_value = float(np.max(arr))
    if max_value <= min_value:
        return np.zeros(arr.shape, dtype=np.uint8)
    scaled = (arr - min_value) / (max_value - min_value)
    return np.clip(np.round(scaled * 255.0), 0, 255).astype(np.uint8)


def robust_normalize_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    if arr.size == 0:
        return np.zeros(arr.shape, dtype=np.uint8)
    p1 = float(np.percentile(arr, 1))
    p99 = float(np.percentile(arr, 99))
    if p99 <= p1:
        return array_to_uint8(arr)
    scaled = np.clip((arr - p1) / (p99 - p1), 0, 1)
    return np.clip(np.round(scaled * 255.0), 0, 255).astype(np.uint8)


def robust_otsu_threshold(values: Sequence[Any]) -> float:
    arr = clean_numeric(values)
    if arr.size == 0:
        return 0.0
    if float(np.max(arr)) <= float(np.min(arr)):
        return float(np.max(arr))
    try:
        return float(filters.threshold_otsu(arr))
    except Exception:
        return float(np.median(arr))


def make_thresholded_preview(signal: np.ndarray, threshold: float) -> np.ndarray:
    return make_thresholded_preview_with_cleanup(signal, threshold, min_area=0)


def normalize_raw_channel_display_mode(mode: str) -> str:
    token = simplify_token(mode)
    aliases = {
        "GRAY": "gray",
        "GREY": "gray",
        "GRAYSCALE": "gray",
        "GREYSCALE": "gray",
        "COLOR": "color",
        "COLOUR": "color",
        "RGB": "color",
        "PSEUDOCOLOR": "color",
        "PSEUDO": "color",
    }
    if token not in aliases:
        raise ValueError(f"不支持的 RAW_CHANNEL_DISPLAY_MODE: {mode}。可选值: gray / color")
    return aliases[token]


def normalize_positive_preview_display_mode(mode: str) -> str:
    token = simplify_token(mode)
    aliases = {
        "GRAY": "gray",
        "GREY": "gray",
        "GRAYSCALE": "gray",
        "GREYSCALE": "gray",
        "COLOR": "color",
        "COLOUR": "color",
        "RGB": "color",
        "PSEUDOCOLOR": "color",
        "PSEUDO": "color",
        "COLORMAP": "colormap",
        "CMAP": "colormap",
        "HEATMAP": "colormap",
    }
    if token not in aliases:
        raise ValueError(
            f"不支持的 POSITIVE_PREVIEW_DISPLAY_MODE: {mode}。"
            " 可选值: colormap / color / gray"
        )
    return aliases[token]


def filter_positive_mask(mask: np.ndarray, min_area: int) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    min_area = int(max(0, min_area))
    if min_area <= 1 or not np.any(mask):
        return mask
    return np.asarray(remove_small_objects_compat(mask, min_area), dtype=bool)


def make_positive_mask(signal: np.ndarray, threshold: float, min_area: int = 0) -> np.ndarray:
    arr = np.asarray(signal, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    threshold = float(max(0.0, threshold))
    mask = arr > threshold
    return filter_positive_mask(mask, min_area=min_area)


def make_thresholded_preview_with_cleanup(signal: np.ndarray, threshold: float, min_area: int = 0) -> np.ndarray:
    arr = np.asarray(signal, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    mask = make_positive_mask(arr, threshold, min_area=min_area)
    out = arr.copy()
    out[~mask] = 0.0
    return out


def positive_fraction_in_roi(signal: np.ndarray, roi_mask: np.ndarray, threshold: float, min_area: int = 0) -> float:
    roi_mask = np.asarray(roi_mask, dtype=bool)
    positive_mask = make_positive_mask(signal, threshold, min_area=min_area)
    if np.any(roi_mask):
        return float(np.mean(positive_mask[roi_mask]))
    return float(np.mean(positive_mask))


def apply_colormap_to_u8(gray_u8: np.ndarray, cmap_name: str = "gray") -> np.ndarray:
    gray_u8 = np.asarray(gray_u8, dtype=np.uint8)
    if cmap_name.lower() in {"gray", "grey"}:
        return np.repeat(gray_u8[..., None], 3, axis=2)
    cmap = matplotlib.colormaps.get_cmap(cmap_name)
    rgba = cmap(gray_u8.astype(np.float32) / 255.0)
    return np.clip(np.round(rgba[..., :3] * 255.0), 0, 255).astype(np.uint8)


def ensure_rgb_uint8(image: np.ndarray, cmap_name: str = "gray") -> np.ndarray:
    arr = np.asarray(image)
    if arr.ndim == 2:
        gray_u8 = arr if arr.dtype == np.uint8 else robust_normalize_for_display(arr)
        return apply_colormap_to_u8(gray_u8, cmap_name=cmap_name)
    if arr.ndim == 3 and arr.shape[-1] == 3:
        if arr.dtype == np.uint8:
            return arr
        arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
        if float(np.max(arr)) <= 1.0:
            return np.clip(np.round(arr * 255.0), 0, 255).astype(np.uint8)
        return np.clip(np.round(arr), 0, 255).astype(np.uint8)
    raise ValueError(f"不支持的图像形状: {arr.shape}")


def resolve_channel_color(
    channel_name: str,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> Tuple[int, int, int]:
    channel = canonicalize_channel_label(channel_name) or normalize_channel_name(channel_name)
    colors = channel_colors or CHANNEL_PSEUDOCOLOR_RGB
    if channel in colors:
        color = colors[channel]
        return tuple(int(np.clip(float(x), 0, 255)) for x in color)
    return tuple(int(x) for x in CHANNEL_PSEUDOCOLOR_RGB.get(channel, (255, 255, 255)))


def colorize_channel_u8(
    gray_u8: np.ndarray,
    channel_name: str,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> np.ndarray:
    gray_u8 = np.asarray(gray_u8, dtype=np.uint8)
    rgb = np.asarray(resolve_channel_color(channel_name, channel_colors=channel_colors), dtype=np.float32) / 255.0
    scaled = (gray_u8.astype(np.float32) / 255.0)[..., None] * rgb[None, None, :]
    return np.clip(np.round(scaled * 255.0), 0, 255).astype(np.uint8)


def render_channel_panel(
    image: np.ndarray,
    channel_name: str,
    display_mode: str,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> np.ndarray:
    display_mode = normalize_raw_channel_display_mode(display_mode)
    gray_u8 = robust_normalize_for_display(image)
    if display_mode == "gray":
        return apply_colormap_to_u8(gray_u8, cmap_name="gray")
    return colorize_channel_u8(gray_u8, channel_name=channel_name, channel_colors=channel_colors)


def render_positive_preview_panel(
    image: np.ndarray,
    channel_name: str,
    display_mode: str,
    cmap_name: str,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> np.ndarray:
    display_mode = normalize_positive_preview_display_mode(display_mode)
    gray_u8 = robust_normalize_for_display(image)
    if display_mode == "gray":
        return apply_colormap_to_u8(gray_u8, cmap_name="gray")
    if display_mode == "color":
        return colorize_channel_u8(gray_u8, channel_name=channel_name, channel_colors=channel_colors)
    return apply_colormap_to_u8(gray_u8, cmap_name=cmap_name)


@lru_cache(maxsize=16)
def get_overlay_font(font_size: int) -> Any:
    if ImageFont is None:
        return None
    candidates = [
        "/System/Library/Fonts/PingFang.ttc",
        "/System/Library/Fonts/STHeiti Light.ttc",
        "/Library/Fonts/Arial Unicode.ttf",
    ]
    for font_path in candidates:
        try:
            if Path(font_path).exists():
                return ImageFont.truetype(font_path, int(font_size))
        except Exception:
            continue
    try:
        return ImageFont.load_default()
    except Exception:
        return None


def draw_multiline_text(
    canvas: np.ndarray,
    text: str,
    x: int,
    y: int,
    font_scale: float = 0.42,
    color: Tuple[int, int, int] = (245, 245, 245),
    thickness: int = 1,
    line_gap: int = 6,
) -> None:
    if Image is not None and ImageDraw is not None:
        pil_image = Image.fromarray(canvas)
        draw = ImageDraw.Draw(pil_image)
        font = get_overlay_font(max(12, int(22 * font_scale)))
        cursor_y = int(y)
        for raw_line in str(text).splitlines():
            line = raw_line.strip()
            if line:
                draw.text((int(x), cursor_y), line, fill=color, font=font)
                try:
                    bbox = draw.textbbox((int(x), cursor_y), line, font=font)
                    line_height = max(14, int(bbox[3] - bbox[1]))
                except Exception:
                    line_height = max(14, int(22 * font_scale))
            else:
                line_height = max(14, int(22 * font_scale))
            cursor_y += line_height + int(line_gap)
        canvas[:] = np.asarray(pil_image)
        return

    cursor_y = int(y)
    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if line:
            cv2.putText(
                canvas,
                line,
                (int(x), cursor_y),
                cv2.FONT_HERSHEY_SIMPLEX,
                float(font_scale),
                color,
                int(thickness),
                cv2.LINE_AA,
            )
        cursor_y += int(18 * font_scale) + int(line_gap)


def save_lossless_panel_grid(
    out_path: Path,
    title: str,
    rows: Sequence[Sequence[Tuple[str, np.ndarray]]],
    panel_gap: int = 24,
    title_band_height: int = 58,
    header_height: int = 72,
    bg_color: Tuple[int, int, int] = (18, 18, 18),
) -> Path:
    if not rows or not rows[0]:
        raise ValueError("rows 不能为空。")

    nrows = len(rows)
    ncols = max(len(row) for row in rows)
    sample_image = ensure_rgb_uint8(rows[0][0][1])
    panel_h, panel_w = sample_image.shape[:2]

    canvas_h = header_height + panel_gap + nrows * (title_band_height + panel_h) + (nrows * panel_gap)
    canvas_w = panel_gap + ncols * panel_w + (ncols * panel_gap)
    canvas = np.full((canvas_h, canvas_w, 3), bg_color, dtype=np.uint8)

    draw_multiline_text(canvas, title, x=panel_gap, y=30, font_scale=0.52, thickness=1, line_gap=8)

    for row_idx, row in enumerate(rows):
        y0 = header_height + panel_gap + row_idx * (title_band_height + panel_h + panel_gap)
        for col_idx, (caption, image) in enumerate(row):
            x0 = panel_gap + col_idx * (panel_w + panel_gap)
            title_box = canvas[y0:y0 + title_band_height, x0:x0 + panel_w]
            title_box[:] = (28, 28, 28)
            draw_multiline_text(title_box, caption, x=10, y=24, font_scale=0.40, thickness=1, line_gap=5)
            canvas[y0 + title_band_height:y0 + title_band_height + panel_h, x0:x0 + panel_w] = ensure_rgb_uint8(image)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    bgr = cv2.cvtColor(canvas, cv2.COLOR_RGB2BGR)
    if not cv2.imwrite(str(out_path), bgr, [cv2.IMWRITE_PNG_COMPRESSION, 0]):
        raise IOError(f"无法保存图片: {out_path}")
    return out_path


def save_rgb_image(out_path: Path, image: np.ndarray) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rgb = ensure_rgb_uint8(image)
    bgr = cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)
    if not cv2.imwrite(str(out_path), bgr, [cv2.IMWRITE_PNG_COMPRESSION, 0]):
        raise IOError(f"无法保存图片: {out_path}")
    return out_path


def build_merged_preview(
    channel_images: Dict[str, np.ndarray],
    positive_thresholds: Optional[Dict[str, PositiveThresholdInfo]] = None,
    thresholded: bool = False,
    min_positive_object_area: int = 0,
    positive_object_area_rules: Optional[Dict[str, int]] = None,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> np.ndarray:
    if not channel_images:
        return np.zeros((1, 1, 3), dtype=np.uint8)

    sample = next(iter(channel_images.values()))
    sample = np.asarray(sample)
    merged = np.zeros((sample.shape[0], sample.shape[1], 3), dtype=np.uint8)

    ordered_channels = [name for name in CANONICAL_CHANNELS if name in channel_images]
    ordered_channels.extend(
        name for name in channel_images.keys()
        if normalize_channel_name(name) not in {normalize_channel_name(x) for x in ordered_channels}
    )

    for channel_name in ordered_channels:
        preview = np.asarray(channel_images[channel_name], dtype=np.float64)
        threshold_info = (positive_thresholds or {}).get(normalize_channel_name(channel_name))
        if thresholded and threshold_info is not None:
            channel_min_area = resolve_positive_object_area(
                channel_name,
                default_min_area=min_positive_object_area,
                area_rules=positive_object_area_rules,
            )
            preview = make_thresholded_preview_with_cleanup(
                preview,
                threshold_info.threshold,
                min_area=channel_min_area,
            )
        colored = colorize_channel_u8(
            robust_normalize_for_display(preview),
            channel_name=channel_name,
            channel_colors=channel_colors,
        )
        merged = np.maximum(merged, colored)

    return merged


def build_denoised_merged_image(
    image_path: Path,
    group_name: str,
    channel_images: Dict[str, np.ndarray],
    merged_dir: Path,
    min_positive_object_area: int = 0,
    positive_object_area_rules: Optional[Dict[str, int]] = None,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> Optional[Path]:
    if not channel_images:
        return None
    make_output_dir(merged_dir)
    merged = build_merged_preview(
        channel_images,
        thresholded=False,
        min_positive_object_area=min_positive_object_area,
        positive_object_area_rules=positive_object_area_rules,
        channel_colors=channel_colors,
    )
    out_path = merged_dir / f"{build_output_basename(group_name, image_path)}_merged_denoised.png"
    return save_rgb_image(out_path, merged)


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=np.float64)
    out = np.full(p.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(p)
    if not np.any(valid):
        return out
    pv = p[valid]
    n = pv.size
    order = np.argsort(pv)
    ranks = np.arange(1, n + 1, dtype=np.float64)
    adjusted_sorted = pv[order] * n / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0, 1)
    restored = np.empty_like(pv)
    restored[order] = adjusted_sorted
    out[valid] = restored
    return out


def remove_small_objects_compat(mask: np.ndarray, min_size: int) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return morphology.remove_small_objects(mask, min_size=int(min_size))
        except TypeError:
            return morphology.remove_small_objects(mask, max_size=int(min_size))


def remove_small_holes_compat(mask: np.ndarray, area_threshold: int) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return morphology.remove_small_holes(mask, area_threshold=int(area_threshold))
        except TypeError:
            return morphology.remove_small_holes(mask, max_size=int(area_threshold))


def requested_channels_from_config(config: AnalysisConfig) -> List[str]:
    return unique_in_order(
        [config.dapi_channel]
        + list(config.intensity_channels)
        + [x for pair in config.colocalization_pairs for x in pair]
    )


# ============================================================================
# 图像读取与通道检测
# ============================================================================


def discover_images(folder: Path, patterns: Sequence[str], recursive: bool = False) -> List[Path]:
    if not folder.exists():
        raise FileNotFoundError(f"文件夹不存在: {folder}")
    if not folder.is_dir():
        raise NotADirectoryError(f"不是文件夹: {folder}")

    found: List[Path] = []
    for pattern in patterns:
        if recursive:
            found.extend(folder.rglob(pattern))
        else:
            found.extend(folder.glob(pattern))

    unique_paths = sorted({path.resolve() for path in found if path.is_file()}, key=natural_sort_key)
    return [Path(path) for path in unique_paths]


def inspect_sample_folder(
    sample_dir: Path,
    patterns: Sequence[str],
    filename_channel_map: Optional[Dict[str, str]] = None,
    filename_role_map: Optional[Dict[str, str]] = None,
) -> SampleFolderContents:
    files = discover_images(sample_dir, patterns, recursive=False)
    channel_files: Dict[str, Path] = {}
    duplicate_channel_files: Dict[str, List[Path]] = {}
    overlay_files: List[Path] = []
    unclassified_files: List[Path] = []

    for file_path in files:
        if is_overlay_like_path(file_path, filename_role_map=filename_role_map):
            overlay_files.append(file_path)
            continue

        channel_name = infer_channel_from_path(file_path, filename_channel_map=filename_channel_map)
        if channel_name is None:
            unclassified_files.append(file_path)
            continue

        if channel_name in channel_files:
            duplicate_channel_files.setdefault(channel_name, [channel_files[channel_name]]).append(file_path)
            continue

        channel_files[channel_name] = file_path

    return SampleFolderContents(
        sample_dir=sample_dir,
        channel_files=channel_files,
        duplicate_channel_files=duplicate_channel_files,
        overlay_files=overlay_files,
        unclassified_files=unclassified_files,
    )


def looks_like_sample_folder(
    sample_dir: Path,
    patterns: Sequence[str],
    dapi_channel: str,
    filename_channel_map: Optional[Dict[str, str]] = None,
    filename_role_map: Optional[Dict[str, str]] = None,
) -> Optional[SampleFolderContents]:
    contents = inspect_sample_folder(
        sample_dir,
        patterns,
        filename_channel_map=filename_channel_map,
        filename_role_map=filename_role_map,
    )
    non_dapi_channels = [name for name in contents.channel_files if name != normalize_channel_name(dapi_channel)]
    if normalize_channel_name(dapi_channel) in contents.channel_files and len(non_dapi_channels) >= 1:
        return contents
    return None


def discover_analysis_inputs(folder: Path, config: AnalysisConfig) -> List[AnalysisInput]:
    if not folder.exists():
        raise FileNotFoundError(f"文件夹不存在: {folder}")
    if not folder.is_dir():
        raise NotADirectoryError(f"不是文件夹: {folder}")

    mode = str(config.input_discovery_mode or "auto").strip().lower()
    candidate_dirs = folder.rglob("*") if config.recursive_scan else folder.iterdir()
    sample_inputs: List[AnalysisInput] = []
    sample_dirs: List[Path] = []

    if mode in {"auto", "sample_folders"}:
        for candidate in sorted((path for path in candidate_dirs if path.is_dir()), key=natural_sort_key):
            contents = looks_like_sample_folder(
                candidate,
                config.image_patterns,
                config.dapi_channel,
                filename_channel_map=config.filename_channel_map,
                filename_role_map=config.filename_role_map,
            )
            if contents is None:
                continue
            sample_inputs.append(AnalysisInput(path=candidate, source_mode="sample_folder", sample_contents=contents))
            sample_dirs.append(candidate.resolve())

    if mode == "sample_folders":
        return sorted(sample_inputs, key=lambda item: natural_sort_key(item.path))

    direct_files = discover_images(folder, config.image_patterns, recursive=config.recursive_scan)
    direct_inputs: List[AnalysisInput] = []
    for file_path in direct_files:
        if mode == "auto" and any(sample_dir in file_path.resolve().parents for sample_dir in sample_dirs):
            continue
        direct_inputs.append(AnalysisInput(path=file_path, source_mode="image_file"))

    if mode == "image_files":
        all_inputs = direct_inputs
    else:
        all_inputs = sample_inputs + direct_inputs
    all_inputs = sorted(all_inputs, key=lambda item: natural_sort_key(item.path))
    return all_inputs


def extract_ome_channel_names(ome_xml: str) -> List[str]:
    if not ome_xml:
        return []
    try:
        root = ET.fromstring(ome_xml)
    except Exception:
        return []

    namespace = {}
    if root.tag.startswith("{") and "}" in root.tag:
        namespace["ome"] = root.tag.split("}", 1)[0].strip("{")
        image_nodes = root.findall(".//ome:Image", namespace)
    else:
        image_nodes = root.findall(".//Image")

    if not image_nodes:
        image_nodes = [root]

    for image_node in image_nodes:
        if namespace:
            pixels_node = image_node.find(".//ome:Pixels", namespace)
            channel_nodes = pixels_node.findall(".//ome:Channel", namespace) if pixels_node is not None else []
        else:
            pixels_node = image_node.find(".//Pixels")
            channel_nodes = pixels_node.findall(".//Channel") if pixels_node is not None else []

        if channel_nodes:
            out = []
            for channel in channel_nodes:
                name = (
                    channel.attrib.get("Name")
                    or channel.attrib.get("Fluor")
                    or channel.attrib.get("ID")
                    or ""
                ).strip()
                out.append(name or "")
            return out
    return []


def _drop_singleton_axes(arr: np.ndarray, axes: str) -> Tuple[np.ndarray, str]:
    arr = np.asarray(arr)
    axes = str(axes or "")
    if not axes or len(axes) != arr.ndim:
        return np.squeeze(arr), axes

    keep = [i for i, size in enumerate(arr.shape) if size != 1]
    if len(keep) == arr.ndim:
        return arr, axes
    if not keep:
        return np.squeeze(arr), ""
    arr2 = np.squeeze(arr)
    axes2 = "".join(axes[i] for i in keep)
    return arr2, axes2


def _standardize_to_hwc(arr: np.ndarray, axes: str) -> Tuple[np.ndarray, str]:
    arr = np.asarray(arr)
    axes = str(axes or "")

    if arr.ndim == 2:
        return arr[..., np.newaxis], axes or "YX"

    if arr.ndim != 3:
        raise ValueError(f"暂不支持 {arr.ndim} 维图像，形状为 {arr.shape}。请先导出单张 2D 图像。")

    if axes and len(axes) == 3:
        axes_upper = axes.upper()
        if "Y" in axes_upper and "X" in axes_upper and ("C" in axes_upper or "S" in axes_upper):
            channel_axis = axes_upper.index("C") if "C" in axes_upper else axes_upper.index("S")
            arr = np.moveaxis(arr, channel_axis, -1)
            return arr, axes_upper

    if arr.shape[-1] <= 8 and arr.shape[0] > 8 and arr.shape[1] > 8:
        return arr, axes
    if arr.shape[0] <= 8 and arr.shape[1] > 8 and arr.shape[2] > 8:
        return np.moveaxis(arr, 0, -1), axes

    raise ValueError(
        f"无法自动识别通道轴，图像形状为 {arr.shape}。仅支持 HxWxC 或 CxHxW。"
    )


def load_multichannel_image(path: Path) -> ImageLoadResult:
    suffix = path.suffix.lower()

    if suffix in {".tif", ".tiff"}:
        if tifffile is None:
            raise ModuleNotFoundError(
                "当前环境缺少 tifffile，无法读取 TIFF。请先运行：python -m pip install tifffile"
            )

        notes: List[str] = []
        with tifffile.TiffFile(str(path)) as tif:
            series = tif.series[0]
            axes = str(getattr(series, "axes", "") or "")
            raw = series.asarray()
            raw, axes = _drop_singleton_axes(raw, axes)

            photometric = ""
            try:
                page0 = tif.pages[0]
                photometric_obj = getattr(page0, "photometric", None)
                photometric = getattr(photometric_obj, "name", str(photometric_obj or ""))
            except Exception:
                photometric = ""

            raw_channel_names = extract_ome_channel_names(getattr(tif, "ome_metadata", "") or "")
            if raw_channel_names:
                notes.append(f"检测到 OME/TIFF metadata 通道名: {raw_channel_names}")

        raw = np.asarray(raw)
        raw = np.squeeze(raw)
        image, axes = _standardize_to_hwc(raw, axes)

        return ImageLoadResult(
            image=image.astype(np.float64, copy=False),
            source_kind="tiff",
            axes=axes,
            photometric=str(photometric or ""),
            raw_channel_names=list(raw_channel_names),
            notes=notes,
        )

    image = cv2.imread(str(path), cv2.IMREAD_UNCHANGED)
    if image is None:
        raise ValueError(f"无法读取图像: {path}")

    image = np.asarray(image)
    image = np.squeeze(image)
    image, axes = _standardize_to_hwc(image, "")
    return ImageLoadResult(
        image=image.astype(np.float64, copy=False),
        source_kind="cv2",
        axes=axes,
        photometric="",
        raw_channel_names=[],
        notes=["非 TIFF 图像使用 cv2 读取；3 通道时按 BGR 颜色语义解释。"],
    )


def load_single_channel_image(path: Path) -> Tuple[np.ndarray, List[str]]:
    load_result = load_multichannel_image(path)
    image = np.asarray(load_result.image, dtype=np.float64)
    if image.ndim != 3:
        raise ValueError(f"{path.name} 不是可识别的单通道图像，形状为 {image.shape}。")

    n_channels = int(image.shape[-1])
    notes = list(load_result.notes)
    if n_channels == 1:
        return image[..., 0], notes

    if n_channels in (3, 4):
        notes.append(f"{path.name} 为 {n_channels} 通道单标记图，按 max(R,G,B[,A]) 压缩为灰度后参与组装。")
        return np.max(image[..., : min(n_channels, 3)], axis=-1), notes

    raise ValueError(
        f"{path.name} 读出 {n_channels} 个通道，无法作为单独 marker 图像自动组装。"
    )


def load_sample_folder_multichannel(sample_dir: Path, config: AnalysisConfig) -> Tuple[ImageLoadResult, SampleFolderContents]:
    contents = inspect_sample_folder(
        sample_dir,
        config.image_patterns,
        filename_channel_map=config.filename_channel_map,
        filename_role_map=config.filename_role_map,
    )
    if contents.duplicate_channel_files:
        duplicate_text = "; ".join(
            f"{channel} -> {[path.name for path in paths]}"
            for channel, paths in sorted(contents.duplicate_channel_files.items())
        )
        raise ValueError(f"{sample_dir.name} 中检测到重复通道文件: {duplicate_text}")

    requested_channels = requested_channels_from_config(config)
    missing_channels = [name for name in requested_channels if name not in contents.channel_files]
    if missing_channels:
        available = sorted(contents.channel_files)
        raise ValueError(
            f"{sample_dir.name} 缺少请求通道 {missing_channels}。"
            f" 当前识别到的单通道文件只有: {available}。"
        )

    ordered_channels = [
        name
        for name in unique_in_order(requested_channels + list(CANONICAL_CHANNELS))
        if name in contents.channel_files
    ]

    arrays: List[np.ndarray] = []
    notes: List[str] = [f"由样本文件夹组装多通道图像: {sample_dir}"]
    reference_shape: Optional[Tuple[int, int]] = None

    for channel_name in ordered_channels:
        channel_path = contents.channel_files[channel_name]
        gray_image, image_notes = load_single_channel_image(channel_path)
        gray_image = np.asarray(gray_image, dtype=np.float64)
        if gray_image.ndim != 2:
            raise ValueError(f"{channel_path.name} 压缩后不是 2D 图像，形状为 {gray_image.shape}。")

        if reference_shape is None:
            reference_shape = gray_image.shape
        elif gray_image.shape != reference_shape:
            raise ValueError(
                f"{sample_dir.name} 中各通道图像尺寸不一致："
                f"期望 {reference_shape}，但 {channel_path.name} 为 {gray_image.shape}。"
            )

        arrays.append(gray_image)
        notes.append(f"{channel_name} <- {channel_path.name}")
        notes.extend(image_notes)

    if contents.overlay_files:
        notes.append(
            "忽略 overlay/composite/映射忽略文件: "
            + ", ".join(path.name for path in sorted(contents.overlay_files, key=natural_sort_key))
        )

    image = np.stack(arrays, axis=-1)
    return ImageLoadResult(
        image=image.astype(np.float64, copy=False),
        source_kind="sample_folder",
        axes="YXC",
        photometric="",
        raw_channel_names=list(ordered_channels),
        notes=notes,
    ), contents


def score_channel_for_dapi(
    channel: np.ndarray,
    min_area: int,
    gaussian_blur_size: int,
    min_peak_distance: int,
) -> Dict[str, float]:
    channel = np.asarray(channel, dtype=np.float64)
    channel = np.nan_to_num(channel, nan=0.0, posinf=0.0, neginf=0.0)

    bgsub = np.clip(channel - np.percentile(channel, 5), 0, None)
    u8 = array_to_uint8(bgsub)
    if np.max(u8) <= 0:
        return {
            "score": 0.0,
            "cell_count": 0.0,
            "mask_fraction": 1.0,
            "mean_area": 0.0,
        }

    if gaussian_blur_size > 1:
        smooth = cv2.GaussianBlur(u8, (gaussian_blur_size, gaussian_blur_size), 0)
    else:
        smooth = u8

    threshold = robust_otsu_threshold(smooth.ravel())
    mask = smooth > threshold
    cleanup_threshold = max(8, min_area // 6)
    mask = remove_small_objects_compat(mask, cleanup_threshold)
    mask = remove_small_holes_compat(mask, cleanup_threshold)
    mask = morphology.opening(mask, morphology.disk(1))
    mask = morphology.closing(mask, morphology.disk(1))

    if not np.any(mask):
        return {
            "score": 0.0,
            "cell_count": 0.0,
            "mask_fraction": 1.0,
            "mean_area": 0.0,
        }

    distance = ndi.distance_transform_edt(mask)
    coordinates = feature.peak_local_max(
        distance,
        labels=mask,
        min_distance=max(1, min_peak_distance // 2),
        exclude_border=False,
    )
    markers = np.zeros(mask.shape, dtype=np.int32)
    if coordinates.size > 0:
        for idx, (row, col) in enumerate(coordinates, start=1):
            markers[row, col] = idx
        markers = ndi.label(markers > 0)[0]
    else:
        markers = ndi.label(mask)[0]

    labels = segmentation.watershed(-distance, markers, mask=mask)

    relaxed_min_area = max(20, min_area // 4)
    areas = [
        float(region.area)
        for region in measure.regionprops(labels)
        if region.area >= relaxed_min_area
    ]
    cell_count = len(areas)
    mean_area = float(np.mean(areas)) if areas else 0.0
    mask_fraction = float(np.mean(mask))

    if cell_count <= 0:
        score = 0.0
    else:
        target = max(float(relaxed_min_area), 1.0)
        area_ratio = mean_area / target if mean_area > 0 else 0.0
        area_penalty = math.exp(-abs(math.log(max(area_ratio, 1e-6))))
        score = float(cell_count * area_penalty / (1.0 + 8.0 * mask_fraction))

    return {
        "score": score,
        "cell_count": float(cell_count),
        "mask_fraction": mask_fraction,
        "mean_area": mean_area,
    }


def detect_dapi_index(image: np.ndarray, config: AnalysisConfig) -> Tuple[int, List[Dict[str, float]]]:
    details: List[Dict[str, float]] = []
    for idx in range(image.shape[-1]):
        metrics = score_channel_for_dapi(
            image[..., idx],
            min_area=config.min_nucleus_area,
            gaussian_blur_size=config.gaussian_blur_size,
            min_peak_distance=config.min_peak_distance,
        )
        metrics["index"] = float(idx)
        details.append(metrics)

    if not details:
        raise ValueError("图像没有可用通道。")

    ranked = sorted(details, key=lambda x: x["score"], reverse=True)
    best = ranked[0]
    best_idx = int(best["index"])
    return best_idx, ranked


def build_full_channel_order(
    n_channels: int,
    assignments: Dict[int, str],
) -> List[str]:
    order: List[str] = []
    used_names = set()

    for idx in range(n_channels):
        if idx in assignments:
            name = normalize_channel_name(assignments[idx])
        else:
            name = f"CH{idx}"
        if name in used_names:
            name = f"{name}_{idx}"
        used_names.add(name)
        order.append(name)

    return order


def resolve_channel_mapping(load_result: ImageLoadResult, config: AnalysisConfig, image_path: Path) -> ChannelResolution:
    n_channels = int(load_result.image.shape[-1])
    requested = requested_channels_from_config(config)
    notes: List[str] = list(load_result.notes)

    if n_channels <= 0:
        raise ValueError(f"{image_path.name} 没有检测到通道。")

    if config.file_channel_order:
        manual = [normalize_channel_name(x) for x in config.file_channel_order]
        if len(manual) > n_channels:
            photometric_upper = str(load_result.photometric or "").upper()
            axes_upper = str(load_result.axes or "").upper()
            looks_like_rgb_overlay = (
                n_channels == 3
                and (
                    "OVERLAY" in image_path.name.upper()
                    or "RGB" in photometric_upper
                    or "S" in axes_upper
                )
            )
            overlay_hint = ""
            if looks_like_rgb_overlay:
                overlay_hint = (
                    " 该文件当前读出来是 3 通道 RGB/overlay 合成图，不是 4 个独立原始通道。"
                    " 如果样本实际做了 DAPI/488/594/647 四通道共染，请导出原始多通道 TIFF/OME-TIFF"
                    "（或每个通道单独灰度图），不要使用 Overlay 图。"
                    " Overlay 中 594 和 647 进入显示用颜色后，无法再可靠拆开。"
                )
            raise ValueError(
                f"{image_path.name} 共有 {n_channels} 个通道，但 FILE_CHANNEL_ORDER 配置了 {len(manual)} 个。"
                f"{overlay_hint}"
            )
        assignments = {idx: name for idx, name in enumerate(manual)}
        order = build_full_channel_order(n_channels, assignments)
        channel_to_index = {name: idx for idx, name in enumerate(order)}
        missing = [name for name in requested if name not in channel_to_index]
        if missing:
            raise ValueError(
                f"{image_path.name} 的手动 FILE_CHANNEL_ORDER 未覆盖以下请求通道: {missing}。"
            )
        notes.append(f"使用手动 FILE_CHANNEL_ORDER: {manual}")
        return ChannelResolution(
            channel_order=order,
            channel_to_index=channel_to_index,
            method="manual",
            notes=notes,
        )

    if not config.auto_detect_channels:
        raise ValueError(
            f"{image_path.name} 未提供 FILE_CHANNEL_ORDER，且 AUTO_DETECT_CHANNELS=False。"
        )

    assignments: Dict[int, str] = {}
    methods: List[str] = []

    # 1) 先尝试 metadata
    if load_result.raw_channel_names:
        for idx, raw_name in enumerate(load_result.raw_channel_names[:n_channels]):
            canonical = canonicalize_channel_label(raw_name)
            if canonical and canonical not in assignments.values():
                assignments[idx] = canonical
                notes.append(f"metadata: index {idx} -> {canonical} (raw='{raw_name}')")
        if assignments:
            methods.append("metadata")

    # 2) RGB/BGR 合成图按颜色语义自动映射
    photometric_upper = str(load_result.photometric or "").upper()
    is_tiff_rgb = (load_result.source_kind == "tiff") and ("S" in str(load_result.axes).upper() or "RGB" in photometric_upper)
    is_cv2_color = (load_result.source_kind == "cv2") and (n_channels in (3, 4))

    if is_tiff_rgb and n_channels >= 3:
        # tifffile 读 RGB TIFF 时，index 0/1/2 通常是 R/G/B
        rgb_map = ["594", "488", "DAPI", "ALPHA"]
        for idx in range(min(n_channels, len(rgb_map))):
            if idx not in assignments and rgb_map[idx] != "ALPHA":
                assignments[idx] = rgb_map[idx]
        methods.append("rgb_tiff_semantics")
        notes.append("检测到 RGB TIFF：按 RGB -> [594, 488, DAPI] 映射。")

    elif is_cv2_color and n_channels >= 3:
        # cv2 读普通彩图时是 BGR
        bgr_map = ["DAPI", "488", "594", "ALPHA"]
        for idx in range(min(n_channels, len(bgr_map))):
            if idx not in assignments and bgr_map[idx] != "ALPHA":
                assignments[idx] = bgr_map[idx]
        methods.append("bgr_cv2_semantics")
        notes.append("检测到 cv2 彩图：按 BGR -> [DAPI, 488, 594] 映射。")

    # 3) 如果 DAPI 还没解决，启发式找 DAPI
    if config.dapi_channel not in assignments.values():
        dapi_idx, ranked = detect_dapi_index(load_result.image, config)
        assignments[dapi_idx] = config.dapi_channel
        methods.append("dapi_heuristic")
        rank_text = ", ".join(
            [
                f"idx {int(item['index'])}: score={item['score']:.3f}, count={item['cell_count']:.0f}, mask={item['mask_fraction']:.3f}"
                for item in ranked
            ]
        )
        notes.append(f"DAPI 启发式评分 -> {rank_text}")
        notes.append(f"启发式选择 DAPI = index {dapi_idx}")

    assigned_names = set(assignments.values())
    missing_requested = [name for name in requested if name not in assigned_names]

    if missing_requested:
        remaining_indices = [idx for idx in range(n_channels) if idx not in assignments]
        if len(missing_requested) == 1 and len(remaining_indices) >= 1:
            assignments[remaining_indices[0]] = missing_requested[0]
            notes.append(
                f"只剩 1 个未分配请求通道，自动补齐: index {remaining_indices[0]} -> {missing_requested[0]}"
            )
            methods.append("single_remaining_fill")
        elif config.strict_channel_detection:
            raise ValueError(
                f"{image_path.name} 的通道仍然存在歧义：已确定 {assignments}，"
                f"但还缺少 {missing_requested}。"
                " 这通常说明文件是没有 metadata 的原始灰度通道堆栈，无法仅凭颜色自动区分 488/594。"
                " 请手动设置 FILE_CHANNEL_ORDER，例如 ['594', '488', 'DAPI']。"
            )
        else:
            fallback = [normalize_channel_name(x) for x in config.fallback_channel_order]
            fallback = [x for x in fallback if x not in assigned_names]
            remaining_indices = [idx for idx in range(n_channels) if idx not in assignments]
            for idx, name in zip(remaining_indices, fallback):
                assignments[idx] = name
            methods.append("fallback_order")
            notes.append(f"使用 FALLBACK_CHANNEL_ORDER 兜底: {config.fallback_channel_order}")

    order = build_full_channel_order(n_channels, assignments)
    channel_to_index = {name: idx for idx, name in enumerate(order)}
    missing = [name for name in requested if name not in channel_to_index]
    if missing:
        raise ValueError(
            f"{image_path.name} 自动检测完成后仍缺少请求通道: {missing}。"
        )

    return ChannelResolution(
        channel_order=order,
        channel_to_index=channel_to_index,
        method=" + ".join(unique_in_order(methods)) if methods else "auto",
        notes=notes,
    )


# ============================================================================
# 配置校验
# ============================================================================


def validate_config(config: AnalysisConfig) -> AnalysisConfig:
    mode = str(config.input_discovery_mode or "auto").strip().lower()
    mode_aliases = {
        "auto": "auto",
        "sample_folders": "sample_folders",
        "sample_folder": "sample_folders",
        "folders": "sample_folders",
        "folder": "sample_folders",
        "image_files": "image_files",
        "image_file": "image_files",
        "images": "image_files",
        "files": "image_files",
    }
    if mode not in mode_aliases:
        raise ValueError(
            f"INPUT_DISCOVERY_MODE 不支持: {config.input_discovery_mode}。"
            " 可选值: auto / sample_folders / image_files"
        )
    config.input_discovery_mode = mode_aliases[mode]

    if config.file_channel_order:
        config.file_channel_order = [normalize_channel_name(x) for x in config.file_channel_order]
        if len(config.file_channel_order) != len(set(config.file_channel_order)):
            raise ValueError(f"FILE_CHANNEL_ORDER 中存在重复通道名: {config.file_channel_order}")
        for name in config.file_channel_order:
            if not ALLOWED_CHANNEL_NAME_PATTERN.match(name):
                raise ValueError(f"不支持的通道名称格式: {name}")

    config.filename_role_map = normalize_filename_role_map(config.filename_role_map)
    config.filename_channel_map = normalize_filename_channel_map(config.filename_channel_map)
    config.fallback_channel_order = [normalize_channel_name(x) for x in config.fallback_channel_order]
    config.positive_threshold_rules = _normalize_positive_threshold_rules(config.positive_threshold_rules)
    config.background_percentile_rules = _normalize_background_percentile_rules(config.background_percentile_rules)
    config.positive_object_area_rules = _normalize_positive_object_area_rules(config.positive_object_area_rules)
    config.bleedthrough_rules = _normalize_bleedthrough_rules(config.bleedthrough_rules)
    config.dapi_channel = normalize_channel_name(config.dapi_channel)
    config.intensity_channels = [normalize_channel_name(x) for x in config.intensity_channels]
    config.colocalization_pairs = normalize_colocalization_pairs(config.colocalization_pairs)
    config.channel_colors = _normalize_channel_colors(config.channel_colors)
    config.raw_channel_display_mode = normalize_raw_channel_display_mode(config.raw_channel_display_mode)
    config.positive_preview_display_mode = normalize_positive_preview_display_mode(config.positive_preview_display_mode)

    if config.auto_colocalization_pairs:
        config.colocalization_pairs = build_all_channel_pairs(config.intensity_channels)

    for channel_name, rule in config.positive_threshold_rules.items():
        method = str(rule.get("method", "otsu")).strip().lower()
        if method not in {"otsu", "manual", "percentile"}:
            raise ValueError(f"不支持的 positive threshold method: {channel_name} -> {method}")
        rule["method"] = method
        if "scale" in rule:
            rule["scale"] = float(max(0.0, rule["scale"]))
        if "offset" in rule:
            rule["offset"] = float(rule["offset"])
        if "min_value" in rule:
            rule["min_value"] = float(max(0.0, rule["min_value"]))
        if method == "manual":
            rule["value"] = float(max(0.0, rule.get("value", 0.0)))
        if method == "percentile":
            rule["percentile"] = float(np.clip(rule.get("percentile", 99.0), 0.0, 100.0))

    for target_name, rule in config.bleedthrough_rules.items():
        source_name = normalize_channel_name(rule.get("source", "")) if rule.get("source", "") else ""
        if source_name and target_name == source_name:
            raise ValueError(f"串色校正中 target 与 source 不能相同: {target_name}")
        mode = str(rule.get("mode", "auto")).strip().lower()
        if mode not in {"auto", "manual"}:
            raise ValueError(f"不支持的 bleed-through mode: {mode}")
        estimate_mask = str(rule.get("estimate_mask", "nuclei")).strip().lower()
        if estimate_mask not in {"nuclei", "roi", "all"}:
            raise ValueError(f"不支持的 estimate_mask: {estimate_mask}")
        if "ratio_percentile" in rule:
            rule["ratio_percentile"] = float(np.clip(rule["ratio_percentile"], 0.0, 100.0))
        if "source_threshold_percentile" in rule:
            rule["source_threshold_percentile"] = float(np.clip(rule["source_threshold_percentile"], 0.0, 100.0))
        if "max_coefficient" in rule:
            rule["max_coefficient"] = float(max(0.0, rule["max_coefficient"]))
        if "min_pixels" in rule:
            rule["min_pixels"] = int(max(1, rule["min_pixels"]))
        if mode == "manual":
            rule["coefficient"] = float(max(0.0, rule.get("coefficient", 0.0)))

    requested = requested_channels_from_config(config)
    if len(requested) != len(set(requested)):
        raise ValueError(f"请求通道中存在重复项: {requested}")

    if config.dapi_channel not in requested:
        raise ValueError("DAPI_CHANNEL 不在请求通道列表中。")

    config.gaussian_blur_size = ensure_odd(config.gaussian_blur_size)
    config.min_nucleus_area = int(max(1, config.min_nucleus_area))
    config.min_positive_object_area = int(max(0, config.min_positive_object_area))
    config.min_peak_distance = int(max(1, config.min_peak_distance))
    config.mask_dilation_radius = int(max(0, config.mask_dilation_radius))
    config.background_percentile = float(np.clip(config.background_percentile, 0.0, 100.0))
    try:
        matplotlib.colormaps.get_cmap(config.positive_preview_colormap)
    except Exception as exc:
        raise ValueError(f"不支持的 POSITIVE_PREVIEW_COLORMAP: {config.positive_preview_colormap}") from exc

    config.group1_dir = Path(config.group1_dir)
    config.group2_dir = Path(config.group2_dir)
    config.output_dir = Path(config.output_dir)
    return config


# ============================================================================
# 图像分析
# ============================================================================


def segment_nuclei(
    dapi: np.ndarray,
    min_area: int,
    gaussian_blur_size: int,
    min_peak_distance: int,
) -> Tuple[np.ndarray, np.ndarray, int, float]:
    dapi = np.asarray(dapi, dtype=np.float64)
    dapi = np.nan_to_num(dapi, nan=0.0, posinf=0.0, neginf=0.0)
    if dapi.ndim != 2:
        raise ValueError(f"DAPI 通道必须是 2D，当前形状: {dapi.shape}")

    dapi_bg = np.clip(dapi - np.percentile(dapi, 5), 0, None)
    dapi_u8 = array_to_uint8(dapi_bg)

    if gaussian_blur_size > 1:
        dapi_smooth = cv2.GaussianBlur(dapi_u8, (gaussian_blur_size, gaussian_blur_size), 0)
    else:
        dapi_smooth = dapi_u8

    if int(np.max(dapi_smooth)) <= 0:
        empty_labels = np.zeros(dapi.shape, dtype=np.int32)
        empty_mask = np.zeros(dapi.shape, dtype=bool)
        return empty_labels, empty_mask, 0, 0.0

    threshold = robust_otsu_threshold(dapi_smooth.ravel())
    mask = dapi_smooth > threshold
    cleanup_threshold = max(3, (min_area // 2) - 1)
    mask = remove_small_objects_compat(mask, cleanup_threshold)
    mask = remove_small_holes_compat(mask, cleanup_threshold)
    mask = morphology.opening(mask, morphology.disk(1))
    mask = morphology.closing(mask, morphology.disk(1))

    if not np.any(mask):
        empty_labels = np.zeros(dapi.shape, dtype=np.int32)
        empty_mask = np.zeros(dapi.shape, dtype=bool)
        return empty_labels, empty_mask, 0, 0.0

    distance = ndi.distance_transform_edt(mask)
    coordinates = feature.peak_local_max(
        distance,
        labels=mask,
        min_distance=min_peak_distance,
        exclude_border=False,
    )

    markers = np.zeros(mask.shape, dtype=np.int32)
    if coordinates.size > 0:
        for idx, (row, col) in enumerate(coordinates, start=1):
            markers[row, col] = idx
        markers = ndi.label(markers > 0)[0]
    else:
        markers = ndi.label(mask)[0]

    labels = segmentation.watershed(-distance, markers, mask=mask)

    filtered = np.zeros_like(labels, dtype=np.int32)
    areas: List[float] = []
    next_label = 1
    for region in measure.regionprops(labels):
        if region.area >= min_area:
            filtered[labels == region.label] = next_label
            next_label += 1
            areas.append(float(region.area))

    cell_count = len(areas)
    mean_area = float(np.mean(areas)) if areas else 0.0
    return filtered, filtered > 0, cell_count, mean_area


def dilate_mask(mask: np.ndarray, radius: int) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    if radius <= 0 or not np.any(mask):
        return mask
    return morphology.dilation(mask, morphology.disk(radius))


def estimate_background(channel: np.ndarray, roi_mask: np.ndarray, percentile: float) -> float:
    channel = np.asarray(channel, dtype=np.float64)
    channel = np.nan_to_num(channel, nan=0.0, posinf=0.0, neginf=0.0)
    roi_mask = np.asarray(roi_mask, dtype=bool)

    outside = channel[~roi_mask]
    if outside.size >= 100:
        return float(np.percentile(outside, percentile))
    return float(np.percentile(channel, percentile))


def prepare_signal(channel: np.ndarray, roi_mask: np.ndarray, background_percentile: float) -> Tuple[np.ndarray, float]:
    bg = estimate_background(channel, roi_mask, background_percentile)
    signal = np.asarray(channel, dtype=np.float64) - bg
    signal = np.clip(signal, 0, None)
    signal = np.nan_to_num(signal, nan=0.0, posinf=0.0, neginf=0.0)
    return signal, bg


def _normalize_positive_threshold_rules(rules: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    normalized: Dict[str, Dict[str, Any]] = {}
    for channel_name, rule in (rules or {}).items():
        if rule is None:
            continue
        channel = normalize_channel_name(channel_name)
        if isinstance(rule, (int, float)):
            normalized[channel] = {"method": "manual", "value": float(rule)}
            continue
        if not isinstance(rule, dict):
            raise ValueError(f"POSITIVE_THRESHOLD_RULES[{channel_name}] 必须是 dict 或数字。")
        normalized[channel] = dict(rule)
    return normalized


def _normalize_positive_object_area_rules(rules: Dict[str, Any]) -> Dict[str, int]:
    normalized: Dict[str, int] = {}
    for channel_name, value in (rules or {}).items():
        channel = normalize_channel_name(channel_name)
        normalized[channel] = int(max(0, value))
    return normalized


def resolve_positive_object_area(
    channel_name: str,
    default_min_area: int,
    area_rules: Optional[Dict[str, int]] = None,
) -> int:
    channel = normalize_channel_name(channel_name)
    rules = area_rules or {}
    return int(max(0, rules.get(channel, default_min_area)))


def _normalize_background_percentile_rules(rules: Dict[str, Any]) -> Dict[str, float]:
    normalized: Dict[str, float] = {}
    for channel_name, value in (rules or {}).items():
        channel = normalize_channel_name(channel_name)
        normalized[channel] = float(np.clip(value, 0.0, 100.0))
    return normalized


def resolve_background_percentile(
    channel_name: str,
    default_percentile: float,
    percentile_rules: Optional[Dict[str, float]] = None,
) -> float:
    channel = normalize_channel_name(channel_name)
    rules = percentile_rules or {}
    return float(np.clip(rules.get(channel, default_percentile), 0.0, 100.0))


def _normalize_channel_colors(raw_map: Dict[str, Any]) -> Dict[str, Tuple[int, int, int]]:
    normalized: Dict[str, Tuple[int, int, int]] = {}
    for channel_name, raw_color in (raw_map or {}).items():
        channel = normalize_channel_name(channel_name)
        if isinstance(raw_color, str):
            text = str(raw_color).strip()
            if text.startswith("#"):
                text = text[1:]
            if len(text) != 6 or not re.fullmatch(r"[0-9A-Fa-f]{6}", text):
                raise ValueError(f"CHANNEL_COLORS[{channel_name}] 颜色格式错误: {raw_color}")
            color = tuple(int(text[i:i + 2], 16) for i in (0, 2, 4))
        elif isinstance(raw_color, (list, tuple)) and len(raw_color) == 3:
            color = tuple(int(np.clip(float(x), 0, 255)) for x in raw_color)
        else:
            raise ValueError(
                f"CHANNEL_COLORS[{channel_name}] 必须是 RGB 三元组或 #RRGGBB 字符串。"
            )
        normalized[channel] = color
    return normalized


def _normalize_bleedthrough_rules(rules: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    normalized: Dict[str, Dict[str, Any]] = {}
    for target_name, rule in (rules or {}).items():
        if rule is None:
            continue
        target = normalize_channel_name(target_name)
        rule_copy = dict(rule)
        if "source" in rule_copy and rule_copy["source"] is not None:
            rule_copy["source"] = normalize_channel_name(rule_copy["source"])
        rule_copy["mode"] = str(rule_copy.get("mode", "auto")).strip().lower()
        rule_copy["estimate_mask"] = str(rule_copy.get("estimate_mask", "nuclei")).strip().lower()
        normalized[target] = rule_copy
    return normalized


def build_bleedthrough_rules_from_simple_config(
    source_map: Dict[str, str],
    default_rule: Optional[Dict[str, Any]] = None,
    manual_coefficients: Optional[Dict[str, float]] = None,
    rule_overrides: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Dict[str, Any]]:
    rules: Dict[str, Dict[str, Any]] = {}
    default_rule = dict(default_rule or {})
    manual_coefficients = manual_coefficients or {}

    for target_name, source_name in (source_map or {}).items():
        target = normalize_channel_name(target_name)
        source = normalize_channel_name(source_name)
        rule = dict(default_rule)
        rule["source"] = source

        manual_key = next(
            (key for key in manual_coefficients if normalize_channel_name(key) == target),
            None,
        )
        if manual_key is not None:
            rule["mode"] = "manual"
            rule["coefficient"] = float(max(0.0, manual_coefficients[manual_key]))

        rules[target] = rule

    for target_name, raw_override in (rule_overrides or {}).items():
        if raw_override is None:
            continue
        target = normalize_channel_name(target_name)
        override = dict(raw_override)
        rule = dict(rules.get(target, {}))
        rule.update(override)
        if "source" not in rule or rule.get("source") in {None, ""}:
            raise ValueError(
                f"BLEEDTHROUGH_RULE_OVERRIDES[{target_name}] 缺少 source。"
                " 如果该通道不在 BLEEDTHROUGH_SOURCE_MAP 中，请在 override 里显式写 source。"
            )
        rules[target] = rule

    return rules


def resolve_positive_threshold(
    channel_name: str,
    signal: np.ndarray,
    roi_mask: np.ndarray,
    config: AnalysisConfig,
) -> PositiveThresholdInfo:
    channel = normalize_channel_name(channel_name)
    roi_mask = np.asarray(roi_mask, dtype=bool)
    signal = np.asarray(signal, dtype=np.float64)

    roi_values = signal[roi_mask] if np.any(roi_mask) else signal.ravel()
    values = clean_numeric(roi_values)
    if values.size == 0:
        return PositiveThresholdInfo(
            channel_name=channel,
            threshold=0.0,
            base_threshold=0.0,
            method="empty",
            note="ROI 内没有可用像素，阈值记为 0。",
        )

    rule = dict(config.positive_threshold_rules.get(channel, {}))
    method = str(rule.get("method", "otsu")).strip().lower()
    scale = float(rule.get("scale", 1.0))
    offset = float(rule.get("offset", 0.0))
    min_value = float(rule.get("min_value", 0.0))

    if method == "manual":
        base_threshold = float(rule.get("value", 0.0))
        threshold = base_threshold
        note = f"manual={base_threshold:.4g}"
    elif method == "percentile":
        percentile = float(np.clip(rule.get("percentile", 99.0), 0.0, 100.0))
        base_threshold = float(np.percentile(values, percentile))
        threshold = (base_threshold * scale) + offset
        note = f"percentile={percentile:.1f} | base={base_threshold:.4g} | scale={scale:.3g} | offset={offset:.3g}"
    else:
        base_threshold = robust_otsu_threshold(values)
        threshold = (base_threshold * scale) + offset
        note = f"otsu={base_threshold:.4g} | scale={scale:.3g} | offset={offset:.3g}"
        method = "otsu"

    threshold = float(max(0.0, threshold, min_value))
    if min_value > 0:
        note = f"{note} | min={min_value:.3g}"

    return PositiveThresholdInfo(
        channel_name=channel,
        threshold=threshold,
        base_threshold=float(base_threshold),
        method=method,
        note=note,
    )


def choose_bleedthrough_mask(
    mask_name: str,
    nuclei_mask: np.ndarray,
    roi_mask: np.ndarray,
    shape: Tuple[int, int],
) -> np.ndarray:
    mode = str(mask_name or "nuclei").strip().lower()
    if mode == "roi":
        mask = np.asarray(roi_mask, dtype=bool)
    elif mode == "all":
        mask = np.ones(shape, dtype=bool)
    else:
        mask = np.asarray(nuclei_mask, dtype=bool)
    return mask


def estimate_bleedthrough_coefficient(
    target_signal: np.ndarray,
    source_signal: np.ndarray,
    estimation_mask: np.ndarray,
    ratio_percentile: float = 20.0,
    source_threshold_percentile: float = 75.0,
    min_pixels: int = 100,
    max_coefficient: float = 3.0,
) -> Tuple[float, str]:
    target_signal = np.asarray(target_signal, dtype=np.float64)
    source_signal = np.asarray(source_signal, dtype=np.float64)
    estimation_mask = np.asarray(estimation_mask, dtype=bool)

    valid = (
        estimation_mask
        & np.isfinite(target_signal)
        & np.isfinite(source_signal)
        & (source_signal > 0)
    )
    valid_count = int(np.count_nonzero(valid))
    if valid_count < max(20, min_pixels // 2):
        return 0.0, f"有效估计像素过少: {valid_count}"

    source_values = source_signal[valid]
    source_threshold_percentile = float(np.clip(source_threshold_percentile, 0.0, 100.0))
    ratio_percentile = float(np.clip(ratio_percentile, 0.0, 100.0))
    src_threshold = float(np.percentile(source_values, source_threshold_percentile))
    strong = valid & (source_signal >= src_threshold)
    strong_count = int(np.count_nonzero(strong))

    if strong_count < min_pixels:
        fallback_percentile = 50.0 if source_threshold_percentile > 50.0 else source_threshold_percentile
        src_threshold = float(np.percentile(source_values, fallback_percentile))
        strong = valid & (source_signal >= src_threshold)
        strong_count = int(np.count_nonzero(strong))

    if strong_count < max(20, min_pixels // 2):
        return 0.0, f"source 强像素不足: {strong_count}"

    ratio = target_signal[strong] / np.maximum(source_signal[strong], 1e-9)
    ratio = ratio[np.isfinite(ratio)]
    ratio = ratio[ratio >= 0]

    if ratio.size < max(20, min_pixels // 2):
        return 0.0, f"可用 ratio 像素不足: {ratio.size}"

    low_cut = float(np.percentile(ratio, ratio_percentile))
    selected = ratio[ratio <= low_cut]
    if selected.size == 0:
        coefficient = low_cut
    else:
        coefficient = float(np.median(selected))

    coefficient = float(np.clip(coefficient, 0.0, max_coefficient))
    note = (
        f"auto estimate on {strong_count} strong-source pixels | "
        f"ratio_percentile={ratio_percentile:.1f} | "
        f"source_threshold={source_threshold_percentile:.1f} | "
        f"coef={coefficient:.6g}"
    )
    return coefficient, note


def apply_bleedthrough_correction(
    target_channel_name: str,
    raw_channels: Dict[str, np.ndarray],
    roi_mask: np.ndarray,
    nuclei_mask: np.ndarray,
    background_percentile: float,
    background_percentile_rules: Optional[Dict[str, float]],
    bleedthrough_rules: Dict[str, Dict[str, Any]],
) -> Tuple[np.ndarray, float, BleedthroughCorrectionInfo, np.ndarray, np.ndarray]:
    target_name = normalize_channel_name(target_channel_name)
    target_raw = np.asarray(raw_channels[target_name], dtype=np.float64)
    target_bg_percentile = resolve_background_percentile(
        target_name,
        default_percentile=background_percentile,
        percentile_rules=background_percentile_rules,
    )
    target_signal, target_bg = prepare_signal(target_raw, roi_mask=roi_mask, background_percentile=target_bg_percentile)

    default_info = BleedthroughCorrectionInfo(
        applied=False,
        target_channel=target_name,
        source_channel="",
        mode="none",
        coefficient=0.0,
        note="未启用该通道的串色校正。",
        source_background=float("nan"),
    )
    zero_source = np.zeros_like(target_signal, dtype=np.float64)

    rule = (bleedthrough_rules or {}).get(target_name)
    if not rule:
        return target_signal, target_bg, default_info, target_signal, zero_source

    source_name = normalize_channel_name(rule.get("source", ""))
    if not source_name:
        info = BleedthroughCorrectionInfo(
            applied=False,
            target_channel=target_name,
            source_channel="",
            mode=str(rule.get("mode", "auto")).strip().lower(),
            coefficient=0.0,
            note="串色校正规则缺少 source。",
            source_background=float("nan"),
        )
        return target_signal, target_bg, info, target_signal, zero_source

    if source_name == target_name:
        info = BleedthroughCorrectionInfo(
            applied=False,
            target_channel=target_name,
            source_channel=source_name,
            mode=str(rule.get("mode", "auto")).strip().lower(),
            coefficient=0.0,
            note="source 与 target 相同，已跳过串色校正。",
            source_background=float("nan"),
        )
        return target_signal, target_bg, info, target_signal, zero_source

    if source_name not in raw_channels:
        info = BleedthroughCorrectionInfo(
            applied=False,
            target_channel=target_name,
            source_channel=source_name,
            mode=str(rule.get("mode", "auto")).strip().lower(),
            coefficient=0.0,
            note=f"source 通道不存在: {source_name}",
            source_background=float("nan"),
        )
        return target_signal, target_bg, info, target_signal, zero_source

    source_raw = np.asarray(raw_channels[source_name], dtype=np.float64)
    source_bg_percentile = resolve_background_percentile(
        source_name,
        default_percentile=background_percentile,
        percentile_rules=background_percentile_rules,
    )
    source_signal, source_bg = prepare_signal(source_raw, roi_mask=roi_mask, background_percentile=source_bg_percentile)

    mode = str(rule.get("mode", "auto")).strip().lower()
    if mode == "manual":
        coefficient = float(max(0.0, rule.get("coefficient", 0.0)))
        note = f"manual coefficient={coefficient:.6g}"
    else:
        estimation_mask = choose_bleedthrough_mask(
            mask_name=str(rule.get("estimate_mask", "nuclei")),
            nuclei_mask=nuclei_mask,
            roi_mask=roi_mask,
            shape=target_signal.shape,
        )
        coefficient, note = estimate_bleedthrough_coefficient(
            target_signal=target_signal,
            source_signal=source_signal,
            estimation_mask=estimation_mask,
            ratio_percentile=float(rule.get("ratio_percentile", 20.0)),
            source_threshold_percentile=float(rule.get("source_threshold_percentile", 75.0)),
            min_pixels=int(rule.get("min_pixels", 100)),
            max_coefficient=float(rule.get("max_coefficient", 3.0)),
        )

    source_contribution = coefficient * source_signal
    corrected_signal = np.clip(target_signal - source_contribution, 0, None)
    info = BleedthroughCorrectionInfo(
        applied=bool(coefficient > 0),
        target_channel=target_name,
        source_channel=source_name,
        mode=mode,
        coefficient=float(coefficient),
        note=note,
        source_background=float(source_bg),
    )
    return corrected_signal, target_bg, info, target_signal, source_contribution


def compute_intensity_metrics(
    channel_name: str,
    raw_channel: np.ndarray,
    signal: np.ndarray,
    roi_mask: np.ndarray,
    dapi_signal: np.ndarray,
    cell_count: int,
    dapi_name: str,
    background_value: float,
    positive_threshold: PositiveThresholdInfo,
    min_positive_object_area: int,
) -> Dict[str, Any]:
    roi_mask = np.asarray(roi_mask, dtype=bool)
    channel_name = normalize_channel_name(channel_name)
    dapi_name = normalize_channel_name(dapi_name)

    full_mean_raw = float(np.mean(raw_channel))
    full_mean_bgsub = float(np.mean(signal))
    roi_pixel_count = int(np.count_nonzero(roi_mask))

    if roi_pixel_count > 0:
        roi_values = signal[roi_mask]
        roi_mean_bgsub = float(np.mean(roi_values))
        roi_integrated_bgsub = float(np.sum(roi_values))
        positive_area_fraction = positive_fraction_in_roi(
            signal,
            roi_mask,
            positive_threshold.threshold,
            min_area=min_positive_object_area,
        )
    else:
        roi_mean_bgsub = float("nan")
        roi_integrated_bgsub = 0.0
        positive_area_fraction = float("nan")

    dapi_roi_integrated_bgsub = float(np.sum(dapi_signal[roi_mask])) if roi_pixel_count > 0 else 0.0
    norm_to_dapi = (
        roi_integrated_bgsub / dapi_roi_integrated_bgsub if dapi_roi_integrated_bgsub > 0 else float("nan")
    )
    integrated_per_cell = roi_integrated_bgsub / cell_count if cell_count > 0 else float("nan")

    return {
        f"{channel_name}_background": background_value,
        f"{channel_name}_full_mean_raw": full_mean_raw,
        f"{channel_name}_full_mean_bgsub": full_mean_bgsub,
        f"{channel_name}_roi_mean_bgsub": roi_mean_bgsub,
        f"{channel_name}_roi_integrated_bgsub": roi_integrated_bgsub,
        f"{channel_name}_norm_to_{dapi_name}": norm_to_dapi,
        f"{channel_name}_integrated_per_cell": integrated_per_cell,
        f"{channel_name}_positive_threshold": float(positive_threshold.threshold),
        f"{channel_name}_positive_threshold_base": float(positive_threshold.base_threshold),
        f"{channel_name}_positive_threshold_method": positive_threshold.method,
        f"{channel_name}_positive_threshold_note": positive_threshold.note,
        f"{channel_name}_positive_area_fraction": positive_area_fraction,
    }


def compute_colocalization(
    channel_a_name: str,
    signal_a: np.ndarray,
    channel_b_name: str,
    signal_b: np.ndarray,
    roi_mask: np.ndarray,
    threshold_a: PositiveThresholdInfo,
    threshold_b: PositiveThresholdInfo,
    min_positive_object_area_a: int,
    min_positive_object_area_b: int,
) -> Dict[str, float]:
    a_name = normalize_channel_name(channel_a_name)
    b_name = normalize_channel_name(channel_b_name)
    signal_a = np.asarray(signal_a, dtype=np.float64)
    signal_b = np.asarray(signal_b, dtype=np.float64)
    roi_mask = np.asarray(roi_mask, dtype=bool)

    if not np.any(roi_mask):
        roi_mask = np.ones(signal_a.shape, dtype=bool)

    mask_a = roi_mask & make_positive_mask(signal_a, threshold_a.threshold, min_area=min_positive_object_area_a)
    mask_b = roi_mask & make_positive_mask(signal_b, threshold_b.threshold, min_area=min_positive_object_area_b)
    union_mask = mask_a | mask_b
    intersection_mask = mask_a & mask_b

    if int(np.count_nonzero(union_mask)) >= MIN_PIXELS_FOR_COLOC:
        vec_a = signal_a[union_mask]
        vec_b = signal_b[union_mask]
        if np.std(vec_a) > 0 and np.std(vec_b) > 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pearson_r, _ = stats.pearsonr(vec_a, vec_b)
            pearson_r = float(pearson_r)
        else:
            pearson_r = float("nan")

        numerator = float(np.sum(vec_a * vec_b))
        denominator = float(np.sqrt(np.sum(vec_a ** 2) * np.sum(vec_b ** 2)))
        overlap_coeff = numerator / denominator if denominator > 0 else float("nan")
    else:
        pearson_r = float("nan")
        overlap_coeff = float("nan")

    denom_a = float(np.sum(signal_a[mask_a]))
    denom_b = float(np.sum(signal_b[mask_b]))
    inter_a = float(np.sum(signal_a[intersection_mask]))
    inter_b = float(np.sum(signal_b[intersection_mask]))
    manders_m1 = inter_a / denom_a if denom_a > 0 else float("nan")
    manders_m2 = inter_b / denom_b if denom_b > 0 else float("nan")
    coloc_fraction = (
        float(np.count_nonzero(intersection_mask)) / float(np.count_nonzero(union_mask))
        if np.any(union_mask)
        else float("nan")
    )

    prefix = f"{a_name}_VS_{b_name}"
    return {
        f"{prefix}_pearson_r": pearson_r,
        f"{prefix}_manders_{a_name}": manders_m1,
        f"{prefix}_manders_{b_name}": manders_m2,
        f"{prefix}_overlap_coeff": overlap_coeff,
        f"{prefix}_coloc_fraction": coloc_fraction,
        f"{prefix}_threshold_{a_name}": float(threshold_a.threshold),
        f"{prefix}_threshold_{b_name}": float(threshold_b.threshold),
    }


def build_channel_diagnostic_figure(
    image_path: Path,
    group_name: str,
    image: np.ndarray,
    resolution: ChannelResolution,
    signals: Dict[str, np.ndarray],
    positive_thresholds: Dict[str, PositiveThresholdInfo],
    roi_mask: np.ndarray,
    config: AnalysisConfig,
    diagnostics_dir: Path,
) -> Optional[Path]:
    if image.ndim != 3:
        return None

    make_output_dir(diagnostics_dir)
    n_channels = image.shape[-1]
    raw_channels: Dict[str, np.ndarray] = {
        name: image[..., idx]
        for idx, name in enumerate(resolution.channel_order)
        if idx < image.shape[-1]
    }
    raw_merged = build_merged_preview(raw_channels, thresholded=False, channel_colors=config.channel_colors)
    denoised_merged = build_merged_preview(
        signals,
        positive_thresholds=positive_thresholds,
        thresholded=True,
        min_positive_object_area=config.min_positive_object_area,
        positive_object_area_rules=config.positive_object_area_rules,
        channel_colors=config.channel_colors,
    )

    top_row: List[Tuple[str, np.ndarray]] = [("Raw merged", raw_merged)]
    bottom_row: List[Tuple[str, np.ndarray]] = [("Thresholded merged", denoised_merged)]

    for idx in range(n_channels):
        name = resolution.channel_order[idx]
        top_row.append(
            (
                f"idx {idx}\n{name}",
                render_channel_panel(
                    image[..., idx],
                    channel_name=name,
                    display_mode=config.raw_channel_display_mode,
                    channel_colors=config.channel_colors,
                ),
            )
        )

        threshold_info = positive_thresholds.get(name)
        if name in signals and threshold_info is not None:
            channel_min_area = resolve_positive_object_area(
                name,
                default_min_area=config.min_positive_object_area,
                area_rules=config.positive_object_area_rules,
            )
            preview = make_thresholded_preview_with_cleanup(
                signals[name],
                threshold_info.threshold,
                min_area=channel_min_area,
            )
            positive_fraction = positive_fraction_in_roi(
                signals[name],
                roi_mask,
                threshold_info.threshold,
                min_area=channel_min_area,
            )
            bottom_row.append(
                (
                    f"{name} positive preview\nthr={threshold_info.threshold:.4g} | roi+={positive_fraction:.2%}",
                    render_positive_preview_panel(
                        preview,
                        channel_name=name,
                        display_mode=config.positive_preview_display_mode,
                        cmap_name=config.positive_preview_colormap,
                        channel_colors=config.channel_colors,
                    ),
                )
            )
        else:
            bottom_row.append((f"{name}\nnot analyzed", np.zeros((*image.shape[:2], 3), dtype=np.uint8)))

    method_text = resolution.method or "auto"
    out_path = diagnostics_dir / f"{build_output_basename(group_name, image_path)}_channels.png"
    return save_lossless_panel_grid(
        out_path=out_path,
        title=f"[{group_name}] {image_path.name}\nChannel diagnostics | row1=raw | row2=positive preview | {method_text}",
        rows=[top_row, bottom_row],
    )


def build_qc_figure(
    image_path: Path,
    group_name: str,
    dapi: np.ndarray,
    roi_mask: np.ndarray,
    nuclei_labels: np.ndarray,
    raw_channels: Dict[str, np.ndarray],
    signals: Dict[str, np.ndarray],
    corrections: Dict[str, BleedthroughCorrectionInfo],
    positive_thresholds: Dict[str, PositiveThresholdInfo],
    config: AnalysisConfig,
    resolution: ChannelResolution,
    qc_dir: Path,
) -> Optional[Path]:
    if not config.save_qc_overlays:
        return None

    make_output_dir(qc_dir)
    dapi_u8 = robust_normalize_for_display(dapi)
    boundary = segmentation.find_boundaries(nuclei_labels, mode="outer")

    non_dapi_channels = [
        c for c in unique_in_order(
            config.intensity_channels + [x for p in config.colocalization_pairs for x in p]
        )
        if c != config.dapi_channel
    ]
    if not non_dapi_channels:
        non_dapi_channels = [config.dapi_channel]

    dapi_idx = resolution.channel_to_index.get(config.dapi_channel, -1)
    qc_panels: List[Tuple[str, np.ndarray]] = [
        (
            f"{config.dapi_channel} raw (idx {dapi_idx})",
            render_channel_panel(
                dapi,
                channel_name=config.dapi_channel,
                display_mode=config.raw_channel_display_mode,
                channel_colors=config.channel_colors,
            ),
        )
    ]

    overlay = np.dstack([dapi_u8, dapi_u8, dapi_u8]).astype(np.uint8)
    overlay[boundary] = np.array([255, 0, 0], dtype=np.uint8)
    roi_outline = segmentation.find_boundaries(roi_mask, mode="outer")
    overlay[roi_outline] = np.array([0, 255, 0], dtype=np.uint8)
    qc_panels.append(("Nuclei / ROI QC", overlay))

    merged_channels: Dict[str, np.ndarray] = {config.dapi_channel: dapi}

    for channel_name in non_dapi_channels:
        raw = raw_channels.get(channel_name, dapi)
        signal = signals.get(channel_name, dapi)
        threshold_info = positive_thresholds.get(channel_name)
        channel_min_area = resolve_positive_object_area(
            channel_name,
            default_min_area=config.min_positive_object_area,
            area_rules=config.positive_object_area_rules,
        )
        preview = (
            make_thresholded_preview_with_cleanup(
                signal,
                threshold_info.threshold,
                min_area=channel_min_area,
            )
            if threshold_info is not None
            else signal
        )
        positive_fraction = (
            positive_fraction_in_roi(
                signal,
                roi_mask,
                threshold_info.threshold,
                min_area=channel_min_area,
            )
            if threshold_info is not None
            else float("nan")
        )

        channel_idx = resolution.channel_to_index.get(channel_name, -1)
        qc_panels.append(
            (
                f"{channel_name} raw (idx {channel_idx})",
                render_channel_panel(
                    raw,
                    channel_name=channel_name,
                    display_mode=config.raw_channel_display_mode,
                    channel_colors=config.channel_colors,
                ),
            )
        )
        correction_info = corrections.get(channel_name)
        if correction_info and correction_info.applied:
            threshold_text = f"{threshold_info.threshold:.4g}" if threshold_info is not None else "NA"
            preview_title = (
                f"{channel_name} positive preview\nthr={threshold_text} | roi+={positive_fraction:.2%}"
            )
        elif threshold_info is not None:
            preview_title = (
                f"{channel_name} positive preview\nthr={threshold_info.threshold:.4g} | roi+={positive_fraction:.2%}"
            )
        else:
            preview_title = f"{channel_name} bg-subtracted"
        qc_panels.append(
            (
                preview_title,
                render_positive_preview_panel(
                    preview,
                    channel_name=channel_name,
                    display_mode=config.positive_preview_display_mode,
                    cmap_name=config.positive_preview_colormap,
                    channel_colors=config.channel_colors,
                ),
            )
        )
        merged_channels[channel_name] = preview

    merged_thresholded = build_merged_preview(
        signals,
        positive_thresholds=positive_thresholds,
        thresholded=True,
        min_positive_object_area=config.min_positive_object_area,
        positive_object_area_rules=config.positive_object_area_rules,
        channel_colors=config.channel_colors,
    )
    qc_panels.append(("Merged QC (thresholded)", merged_thresholded))

    method_text = resolution.method or "auto"
    out_path = qc_dir / f"{build_output_basename(group_name, image_path)}_qc.png"
    return save_lossless_panel_grid(
        out_path=out_path,
        title=f"[{group_name}] {image_path.name}\n{method_text} | order={resolution.channel_order}",
        rows=[qc_panels],
    )


def build_bleedthrough_diagnostic_figure(
    image_path: Path,
    group_name: str,
    channel_name: str,
    raw_channel: np.ndarray,
    uncorrected_signal: np.ndarray,
    source_channel_name: str,
    source_contribution: np.ndarray,
    corrected_signal: np.ndarray,
    info: BleedthroughCorrectionInfo,
    raw_channel_display_mode: str,
    channel_colors: Optional[Dict[str, Tuple[int, int, int]]],
    diagnostics_dir: Path,
) -> Optional[Path]:
    if not info.applied:
        return None

    out_dir = diagnostics_dir / "bleedthrough"
    make_output_dir(out_dir)

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    axes[0].imshow(
        render_channel_panel(
            raw_channel,
            channel_name=channel_name,
            display_mode=raw_channel_display_mode,
            channel_colors=channel_colors,
        )
    )
    axes[0].set_title(f"{channel_name} raw")
    axes[0].axis("off")

    axes[1].imshow(
        render_channel_panel(
            uncorrected_signal,
            channel_name=channel_name,
            display_mode=raw_channel_display_mode,
            channel_colors=channel_colors,
        )
    )
    axes[1].set_title(f"{channel_name} bg-subtracted")
    axes[1].axis("off")

    axes[2].imshow(
        render_channel_panel(
            source_contribution,
            channel_name=source_channel_name or info.source_channel or channel_name,
            display_mode=raw_channel_display_mode,
            channel_colors=channel_colors,
        )
    )
    axes[2].set_title(f"{info.source_channel} contribution\n(k={info.coefficient:.3g})")
    axes[2].axis("off")

    axes[3].imshow(
        render_channel_panel(
            corrected_signal,
            channel_name=channel_name,
            display_mode=raw_channel_display_mode,
            channel_colors=channel_colors,
        )
    )
    axes[3].set_title(f"{channel_name} corrected")
    axes[3].axis("off")

    fig.suptitle(
        f"[{group_name}] {image_path.name}\nBleed-through correction: {channel_name} <- {info.source_channel} | {info.mode}",
        fontsize=11,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.90])
    out_path = out_dir / f"{build_output_basename(group_name, image_path)}_{sanitize_filename(channel_name)}_bleedthrough.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out_path


def analyze_loaded_image(
    image_path: Path,
    group_name: str,
    config: AnalysisConfig,
    qc_dir: Path,
    diagnostics_dir: Path,
    merged_dir: Path,
    load_result: ImageLoadResult,
    extra_result_fields: Optional[Dict[str, Any]] = None,
) -> ImageAnalysisResult:
    image = load_result.image
    resolution = resolve_channel_mapping(load_result, config, image_path=image_path)

    channel_arrays: Dict[str, np.ndarray] = {
        name: image[..., idx].astype(np.float64, copy=False)
        for name, idx in resolution.channel_to_index.items()
        if idx < image.shape[-1]
    }

    dapi = channel_arrays[config.dapi_channel]
    nuclei_labels, nuclei_mask, cell_count, mean_nucleus_area = segment_nuclei(
        dapi=dapi,
        min_area=config.min_nucleus_area,
        gaussian_blur_size=config.gaussian_blur_size,
        min_peak_distance=config.min_peak_distance,
    )
    roi_mask = dilate_mask(nuclei_mask, config.mask_dilation_radius)

    height, width = dapi.shape
    total_pixels = int(height * width)
    cell_density_per_mp = float(cell_count / total_pixels * 1_000_000) if total_pixels > 0 else float("nan")

    signals: Dict[str, np.ndarray] = {}
    uncorrected_signals: Dict[str, np.ndarray] = {}
    source_contributions: Dict[str, np.ndarray] = {}
    backgrounds: Dict[str, float] = {}
    corrections: Dict[str, BleedthroughCorrectionInfo] = {}
    positive_thresholds: Dict[str, PositiveThresholdInfo] = {}

    all_needed_channels = requested_channels_from_config(config)
    bleedthrough_rules = config.bleedthrough_rules if config.enable_bleedthrough_correction else {}

    for channel_name in all_needed_channels:
        if channel_name in bleedthrough_rules:
            signal, bg, correction_info, uncorrected_signal, source_contribution = apply_bleedthrough_correction(
                target_channel_name=channel_name,
                raw_channels=channel_arrays,
                roi_mask=roi_mask,
                nuclei_mask=nuclei_mask,
                background_percentile=config.background_percentile,
                background_percentile_rules=config.background_percentile_rules,
                bleedthrough_rules=bleedthrough_rules,
            )
        else:
            channel_bg_percentile = resolve_background_percentile(
                channel_name,
                default_percentile=config.background_percentile,
                percentile_rules=config.background_percentile_rules,
            )
            signal, bg = prepare_signal(
                channel_arrays[channel_name],
                roi_mask=roi_mask,
                background_percentile=channel_bg_percentile,
            )
            correction_info = BleedthroughCorrectionInfo(
                applied=False,
                target_channel=channel_name,
                source_channel="",
                mode="none",
                coefficient=0.0,
                note="未对该通道应用串色校正。",
                source_background=float("nan"),
            )
            uncorrected_signal = signal
            source_contribution = np.zeros_like(signal, dtype=np.float64)

        signals[channel_name] = signal
        backgrounds[channel_name] = bg
        corrections[channel_name] = correction_info
        uncorrected_signals[channel_name] = uncorrected_signal
        source_contributions[channel_name] = source_contribution

    for channel_name in all_needed_channels:
        positive_thresholds[channel_name] = resolve_positive_threshold(
            channel_name=channel_name,
            signal=signals[channel_name],
            roi_mask=roi_mask,
            config=config,
        )

    dapi_signal = signals[config.dapi_channel]
    dapi_roi_pixel_count = int(np.count_nonzero(roi_mask))
    dapi_roi_integrated_bgsub = float(np.sum(dapi_signal[roi_mask])) if dapi_roi_pixel_count > 0 else 0.0
    dapi_roi_mean_bgsub = float(np.mean(dapi_signal[roi_mask])) if dapi_roi_pixel_count > 0 else float("nan")

    result: Dict[str, Any] = {
        "image_name": image_path.stem,
        "image_file": image_path.name,
        "image_path": str(image_path),
        "image_output_id": build_output_basename(group_name, image_path),
        "group": group_name,
        "height_px": height,
        "width_px": width,
        "image_area_px": total_pixels,
        "n_channels": int(image.shape[-1]),
        "cell_count": int(cell_count),
        "cell_density_per_mp": cell_density_per_mp,
        "mean_nucleus_area_px": float(mean_nucleus_area),
        f"{config.dapi_channel}_background": backgrounds[config.dapi_channel],
        f"{config.dapi_channel}_roi_mean_bgsub": dapi_roi_mean_bgsub,
        f"{config.dapi_channel}_roi_integrated_bgsub": dapi_roi_integrated_bgsub,
        "roi_area_fraction": float(np.mean(roi_mask)),
        "channel_detection_method": resolution.method,
        "resolved_channel_order": " | ".join(f"{idx}:{name}" for idx, name in enumerate(resolution.channel_order)),
        "channel_detection_notes": " || ".join(resolution.notes),
        "tiff_axes": load_result.axes,
        "tiff_photometric": load_result.photometric,
    }
    if extra_result_fields:
        result.update(extra_result_fields)

    for canonical_name in CANONICAL_CHANNELS:
        result[f"channel_index_{canonical_name}"] = resolution.channel_to_index.get(canonical_name, np.nan)

    for channel_name, correction_info in corrections.items():
        result[f"{channel_name}_bleedthrough_applied"] = bool(correction_info.applied)
        result[f"{channel_name}_bleedthrough_source"] = correction_info.source_channel
        result[f"{channel_name}_bleedthrough_mode"] = correction_info.mode
        result[f"{channel_name}_bleedthrough_coefficient"] = float(correction_info.coefficient)
        result[f"{channel_name}_bleedthrough_note"] = correction_info.note
        if channel_name in positive_thresholds:
            threshold_info = positive_thresholds[channel_name]
            result[f"{channel_name}_positive_threshold"] = float(threshold_info.threshold)
            result[f"{channel_name}_positive_threshold_base"] = float(threshold_info.base_threshold)
            result[f"{channel_name}_positive_threshold_method"] = threshold_info.method
            result[f"{channel_name}_positive_threshold_note"] = threshold_info.note

    for channel_name in config.intensity_channels:
        channel_min_area = resolve_positive_object_area(
            channel_name,
            default_min_area=config.min_positive_object_area,
            area_rules=config.positive_object_area_rules,
        )
        result.update(
            compute_intensity_metrics(
                channel_name=channel_name,
                raw_channel=channel_arrays[channel_name],
                signal=signals[channel_name],
                roi_mask=roi_mask,
                dapi_signal=dapi_signal,
                cell_count=cell_count,
                dapi_name=config.dapi_channel,
                background_value=backgrounds[channel_name],
                positive_threshold=positive_thresholds[channel_name],
                min_positive_object_area=channel_min_area,
            )
        )

    for channel_a, channel_b in config.colocalization_pairs:
        min_area_a = resolve_positive_object_area(
            channel_a,
            default_min_area=config.min_positive_object_area,
            area_rules=config.positive_object_area_rules,
        )
        min_area_b = resolve_positive_object_area(
            channel_b,
            default_min_area=config.min_positive_object_area,
            area_rules=config.positive_object_area_rules,
        )
        result.update(
            compute_colocalization(
                channel_a_name=channel_a,
                signal_a=signals[channel_a],
                channel_b_name=channel_b,
                signal_b=signals[channel_b],
                roi_mask=roi_mask,
                threshold_a=positive_thresholds[channel_a],
                threshold_b=positive_thresholds[channel_b],
                min_positive_object_area_a=min_area_a,
                min_positive_object_area_b=min_area_b,
            )
        )

    merged_path = build_denoised_merged_image(
        image_path=image_path,
        group_name=group_name,
        channel_images=signals,
        merged_dir=merged_dir,
        min_positive_object_area=config.min_positive_object_area,
        positive_object_area_rules=config.positive_object_area_rules,
        channel_colors=config.channel_colors,
    )
    result["merged_denoised_path"] = str(merged_path) if merged_path is not None else ""

    qc_path = build_qc_figure(
        image_path=image_path,
        group_name=group_name,
        dapi=dapi,
        roi_mask=roi_mask,
        nuclei_labels=nuclei_labels,
        raw_channels=channel_arrays,
        signals=signals,
        corrections=corrections,
        positive_thresholds=positive_thresholds,
        config=config,
        resolution=resolution,
        qc_dir=qc_dir,
    )
    result["qc_overlay_path"] = str(qc_path) if qc_path is not None else ""

    diag_path = None
    if config.save_channel_diagnostics:
        diag_path = build_channel_diagnostic_figure(
            image_path=image_path,
            group_name=group_name,
            image=image,
            resolution=resolution,
            signals=signals,
            positive_thresholds=positive_thresholds,
            roi_mask=roi_mask,
            config=config,
            diagnostics_dir=diagnostics_dir,
        )
    result["channel_diagnostic_path"] = str(diag_path) if diag_path is not None else ""

    bleedthrough_paths: List[str] = []
    if config.enable_bleedthrough_correction and config.save_bleedthrough_diagnostics:
        for channel_name, correction_info in corrections.items():
            bt_path = build_bleedthrough_diagnostic_figure(
                image_path=image_path,
                group_name=group_name,
                channel_name=channel_name,
                raw_channel=channel_arrays[channel_name],
                uncorrected_signal=uncorrected_signals[channel_name],
                source_channel_name=correction_info.source_channel,
                source_contribution=source_contributions[channel_name],
                corrected_signal=signals[channel_name],
                info=correction_info,
                raw_channel_display_mode=config.raw_channel_display_mode,
                channel_colors=config.channel_colors,
                diagnostics_dir=diagnostics_dir,
            )
            if bt_path is not None:
                bleedthrough_paths.append(str(bt_path))
    result["bleedthrough_diagnostic_paths"] = " | ".join(bleedthrough_paths)

    return ImageAnalysisResult(values=result, roi_mask=roi_mask, nuclei_labels=nuclei_labels)


def analyze_single_image(
    image_path: Path,
    group_name: str,
    config: AnalysisConfig,
    qc_dir: Path,
    diagnostics_dir: Path,
    merged_dir: Path,
) -> ImageAnalysisResult:
    load_result = load_multichannel_image(image_path)
    return analyze_loaded_image(
        image_path=image_path,
        group_name=group_name,
        config=config,
        qc_dir=qc_dir,
        diagnostics_dir=diagnostics_dir,
        merged_dir=merged_dir,
        load_result=load_result,
        extra_result_fields={
            "input_mode": "image_file",
            "component_channel_files": str(image_path),
            "ignored_overlay_files": "",
        },
    )


def analyze_sample_folder(
    sample_dir: Path,
    group_name: str,
    config: AnalysisConfig,
    qc_dir: Path,
    diagnostics_dir: Path,
    merged_dir: Path,
    sample_contents: Optional[SampleFolderContents] = None,
) -> ImageAnalysisResult:
    load_result, inspected_contents = load_sample_folder_multichannel(sample_dir, config)
    contents = sample_contents or inspected_contents
    component_text = " | ".join(
        f"{channel}:{path.name}"
        for channel, path in sorted(contents.channel_files.items(), key=lambda item: natural_sort_key(item[0]))
    )
    ignored_overlay_text = " | ".join(path.name for path in sorted(contents.overlay_files, key=natural_sort_key))
    return analyze_loaded_image(
        image_path=sample_dir,
        group_name=group_name,
        config=config,
        qc_dir=qc_dir,
        diagnostics_dir=diagnostics_dir,
        merged_dir=merged_dir,
        load_result=load_result,
        extra_result_fields={
            "input_mode": "sample_folder",
            "component_channel_files": component_text,
            "ignored_overlay_files": ignored_overlay_text,
            "sample_dir": str(sample_dir),
        },
    )


# ============================================================================
# 统计分析
# ============================================================================


def safe_shapiro(values: Sequence[Any]) -> float:
    arr = clean_numeric(values)
    if arr.size < 3:
        return float("nan")
    if np.allclose(arr, arr[0]):
        return float("nan")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return float(stats.shapiro(arr).pvalue)
    except Exception:
        return float("nan")


def hedges_g(group1: Sequence[Any], group2: Sequence[Any]) -> float:
    x1 = clean_numeric(group1)
    x2 = clean_numeric(group2)
    n1, n2 = x1.size, x2.size
    if n1 < 2 or n2 < 2:
        return float("nan")
    sd1 = sample_std(x1)
    sd2 = sample_std(x2)
    if not np.isfinite(sd1) or not np.isfinite(sd2):
        return float("nan")
    pooled_num = (n1 - 1) * (sd1 ** 2) + (n2 - 1) * (sd2 ** 2)
    pooled_den = n1 + n2 - 2
    if pooled_den <= 0:
        return float("nan")
    pooled_sd = math.sqrt(pooled_num / pooled_den) if pooled_num > 0 else 0.0
    if pooled_sd == 0:
        return float("nan")
    d = (float(np.mean(x1)) - float(np.mean(x2))) / pooled_sd
    correction = 1.0
    if (n1 + n2) > 3:
        correction = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    return float(d * correction)


def rank_biserial_from_u(u_stat: float, n1: int, n2: int) -> float:
    if n1 <= 0 or n2 <= 0:
        return float("nan")
    return float(2.0 * u_stat / (n1 * n2) - 1.0)


def compare_two_groups(metric_name: str, group1: Sequence[Any], group2: Sequence[Any]) -> Dict[str, Any]:
    x1 = clean_numeric(group1)
    x2 = clean_numeric(group2)
    n1, n2 = x1.size, x2.size

    result: Dict[str, Any] = {
        "metric": metric_name,
        "n_group1": int(n1),
        "n_group2": int(n2),
        "mean_group1": float(np.mean(x1)) if n1 > 0 else float("nan"),
        "mean_group2": float(np.mean(x2)) if n2 > 0 else float("nan"),
        "median_group1": float(np.median(x1)) if n1 > 0 else float("nan"),
        "median_group2": float(np.median(x2)) if n2 > 0 else float("nan"),
        "sd_group1": sample_std(x1),
        "sd_group2": sample_std(x2),
        "sem_group1": sem(x1),
        "sem_group2": sem(x2),
        "iqr_group1": iqr(x1),
        "iqr_group2": iqr(x2),
        "shapiro_p_group1": safe_shapiro(x1),
        "shapiro_p_group2": safe_shapiro(x2),
        "test": "NA",
        "statistic": float("nan"),
        "p_value_raw": float("nan"),
        "effect_size": float("nan"),
        "effect_size_name": "NA",
    }

    if n1 == 0 or n2 == 0:
        result["test"] = "no_valid_data"
        result["significance_raw"] = "NA"
        return result

    if n1 < 2 or n2 < 2:
        result["test"] = "insufficient_n"
        result["significance_raw"] = "NA"
        return result

    normal_1 = bool(np.isfinite(result["shapiro_p_group1"]) and result["shapiro_p_group1"] > 0.05)
    normal_2 = bool(np.isfinite(result["shapiro_p_group2"]) and result["shapiro_p_group2"] > 0.05)

    if normal_1 and normal_2:
        test_res = stats.ttest_ind(x1, x2, equal_var=False, nan_policy="omit")
        result["test"] = "Welch_t_test"
        result["statistic"] = float(test_res.statistic)
        result["p_value_raw"] = float(test_res.pvalue)
        result["effect_size"] = hedges_g(x1, x2)
        result["effect_size_name"] = "Hedges_g"
    else:
        try:
            test_res = stats.mannwhitneyu(x1, x2, alternative="two-sided", method="auto")
        except TypeError:
            test_res = stats.mannwhitneyu(x1, x2, alternative="two-sided")
        result["test"] = "Mann_Whitney_U"
        result["statistic"] = float(test_res.statistic)
        result["p_value_raw"] = float(test_res.pvalue)
        result["effect_size"] = rank_biserial_from_u(float(test_res.statistic), n1, n2)
        result["effect_size_name"] = "Rank_biserial_r"

    result["significance_raw"] = significance_label(result["p_value_raw"])
    return result


def apply_fdr_correction(stats_df: pd.DataFrame) -> pd.DataFrame:
    stats_df = stats_df.copy()
    if "p_value_raw" not in stats_df.columns:
        stats_df["p_value_fdr_bh"] = np.nan
        stats_df["significance_fdr"] = "NA"
        return stats_df

    raw_p = pd.to_numeric(stats_df["p_value_raw"], errors="coerce").to_numpy(dtype=np.float64)
    adjusted = bh_fdr(raw_p)
    stats_df["p_value_fdr_bh"] = adjusted
    stats_df["significance_fdr"] = [significance_label(x) for x in adjusted]
    return stats_df


# ============================================================================
# 可视化与报告
# ============================================================================


def metric_specs(config: AnalysisConfig) -> List[Tuple[str, str, str]]:
    specs: List[Tuple[str, str, str]] = [
        ("cell_count", "Cell count", "count"),
        ("cell_density_per_mp", "Cell density", "cells / MP"),
        ("mean_nucleus_area_px", "Mean nucleus area", "px^2"),
        ("roi_area_fraction", "ROI area fraction", "fraction"),
    ]

    for channel in config.intensity_channels:
        specs.extend([
            (f"{channel}_roi_mean_bgsub", f"{channel} ROI mean", "a.u."),
            (f"{channel}_integrated_per_cell", f"{channel} per cell", "a.u. / cell"),
            (f"{channel}_norm_to_{config.dapi_channel}", f"{channel}/{config.dapi_channel}", "ratio"),
            (f"{channel}_positive_area_fraction", f"{channel} positive area", "fraction"),
        ])

    for channel_a, channel_b in config.colocalization_pairs:
        prefix = f"{channel_a}_VS_{channel_b}"
        specs.extend([
            (f"{prefix}_pearson_r", f"{channel_a} vs {channel_b} Pearson r", "r"),
            (f"{prefix}_manders_{channel_a}", f"Manders {channel_a}", "coefficient"),
            (f"{prefix}_manders_{channel_b}", f"Manders {channel_b}", "coefficient"),
            (f"{prefix}_overlap_coeff", f"{channel_a} vs {channel_b} overlap", "coefficient"),
            (f"{prefix}_coloc_fraction", f"{channel_a} vs {channel_b} coloc fraction", "fraction"),
        ])

    return specs


def safe_ylim(values: Sequence[float]) -> Tuple[float, float]:
    arr = clean_numeric(values)
    if arr.size == 0:
        return 0.0, 1.0
    min_v = float(np.min(arr))
    max_v = float(np.max(arr))
    if max_v == min_v:
        delta = max(1.0, abs(max_v) * 0.2)
        return min_v - delta, max_v + delta
    margin = (max_v - min_v) * 0.25
    return min_v - margin * 0.4, max_v + margin


def _boxplot_compat(ax: Any, box_data: Sequence[np.ndarray], labels: Sequence[str]) -> Any:
    try:
        return ax.boxplot(box_data, tick_labels=list(labels), patch_artist=True, widths=0.55)
    except TypeError:
        return ax.boxplot(box_data, labels=list(labels), patch_artist=True, widths=0.55)


def create_summary_figure(
    per_image_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    config: AnalysisConfig,
    output_base: Path,
) -> Dict[str, Path]:
    stats_lookup = stats_df.set_index("metric").to_dict(orient="index") if not stats_df.empty else {}
    specs = [(col, title, ylabel) for col, title, ylabel in metric_specs(config) if col in per_image_df.columns]
    if not specs:
        return {}

    n_metrics = len(specs)
    ncols = 3
    nrows = int(math.ceil(n_metrics / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.2 * nrows))
    axes = np.atleast_1d(axes).ravel()

    rng = np.random.default_rng(42)
    group1_df = per_image_df[per_image_df["group"] == config.group1_name]
    group2_df = per_image_df[per_image_df["group"] == config.group2_name]

    for ax, (metric, title, ylabel) in zip(axes, specs):
        data1 = clean_numeric(group1_df[metric].to_numpy()) if metric in group1_df.columns else np.array([])
        data2 = clean_numeric(group2_df[metric].to_numpy()) if metric in group2_df.columns else np.array([])

        if data1.size == 0 and data2.size == 0:
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            ax.axis("off")
            continue

        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylabel)
        ax.set_xticks([1, 2], labels=[config.group1_name, config.group2_name])
        ax.grid(axis="y", alpha=0.25)

        if data1.size > 0 and data2.size > 0:
            box_data = [data1, data2]
            box = _boxplot_compat(ax, box_data, [config.group1_name, config.group2_name])
            for patch, color in zip(box["boxes"], DEFAULT_GROUP_COLORS):
                patch.set_facecolor(color)
                patch.set_alpha(0.50)
            for median in box["medians"]:
                median.set_color("black")
                median.set_linewidth(1.5)

        if data1.size > 0:
            x1 = rng.normal(1.0, 0.045, size=data1.size)
            ax.scatter(x1, data1, s=36, alpha=0.8, zorder=3, color=DEFAULT_GROUP_COLORS[0], edgecolors="white", linewidths=0.5)
        if data2.size > 0:
            x2 = rng.normal(2.0, 0.045, size=data2.size)
            ax.scatter(x2, data2, s=36, alpha=0.8, zorder=3, color=DEFAULT_GROUP_COLORS[1], edgecolors="white", linewidths=0.5)

        combined = np.concatenate([data1, data2]) if data1.size > 0 and data2.size > 0 else (data1 if data1.size > 0 else data2)
        y_min, y_max = safe_ylim(combined)
        ax.set_ylim(y_min, y_max)

        if data1.size == 0 or data2.size == 0:
            ax.text(
                1.5,
                y_max - (y_max - y_min) * 0.05,
                "Insufficient valid data for group comparison",
                ha="center",
                va="top",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.7, edgecolor="lightgray"),
            )
            continue

        if metric in stats_lookup:
            stat_row = stats_lookup[metric]
            sig = stat_row.get("significance_fdr", "NA")
            p_raw = stat_row.get("p_value_raw", np.nan)
            p_fdr = stat_row.get("p_value_fdr_bh", np.nan)
            test_name = str(stat_row.get("test", "NA")).replace("_", " ")
            text_y = y_max - (y_max - y_min) * 0.05
            annotation = f"{sig} | raw p={format_p_value(p_raw)}\nFDR p={format_p_value(p_fdr)} | {test_name}"
            ax.text(1.5, text_y, annotation, ha="center", va="top", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.7, edgecolor="lightgray"))

    for ax in axes[n_metrics:]:
        ax.axis("off")

    fig.suptitle(
        f"IF analysis summary: {config.group1_name} vs {config.group2_name}",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    outputs: Dict[str, Path] = {}
    if config.save_figure_png:
        png_path = output_base.with_suffix(".png")
        fig.savefig(png_path, dpi=config.figure_dpi, bbox_inches="tight")
        outputs["png"] = png_path
    if config.save_figure_pdf:
        pdf_path = output_base.with_suffix(".pdf")
        fig.savefig(pdf_path, dpi=config.figure_dpi, bbox_inches="tight")
        outputs["pdf"] = pdf_path
    plt.close(fig)
    return outputs


def write_text_report(
    report_path: Path,
    config: AnalysisConfig,
    group1_inputs: List[AnalysisInput],
    group2_inputs: List[AnalysisInput],
    per_image_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    errors_df: pd.DataFrame,
) -> None:
    lines: List[str] = []
    lines.append("Immunofluorescence Analysis Report")
    lines.append("=" * 72)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Group 1: {config.group1_name} | folder={config.group1_dir}")
    lines.append(f"Group 2: {config.group2_name} | folder={config.group2_dir}")
    lines.append(f"Inputs found: {config.group1_name}={len(group1_inputs)}, {config.group2_name}={len(group2_inputs)}")
    lines.append(f"Successful analyses: {config.group1_name}={(per_image_df['group'] == config.group1_name).sum()}, {config.group2_name}={(per_image_df['group'] == config.group2_name).sum()}")
    lines.append(f"Input discovery mode: {config.input_discovery_mode}")
    lines.append(f"Manual FILE_CHANNEL_ORDER: {config.file_channel_order if config.file_channel_order else 'AUTO'}")
    lines.append(f"AUTO_DETECT_CHANNELS: {config.auto_detect_channels}")
    lines.append(f"STRICT_CHANNEL_DETECTION: {config.strict_channel_detection}")
    lines.append(f"Fallback channel order: {config.fallback_channel_order}")
    lines.append(f"Positive threshold rules: {json.dumps(config.positive_threshold_rules, ensure_ascii=False)}")
    lines.append(f"Background percentile default: {config.background_percentile}")
    lines.append(f"Background percentile rules: {json.dumps(config.background_percentile_rules, ensure_ascii=False)}")
    lines.append(f"Min positive object area (default): {config.min_positive_object_area}")
    lines.append(f"Positive object area rules: {json.dumps(config.positive_object_area_rules, ensure_ascii=False)}")
    lines.append(f"ENABLE_BLEEDTHROUGH_CORRECTION: {config.enable_bleedthrough_correction}")
    lines.append(f"SAVE_BLEEDTHROUGH_DIAGNOSTICS: {config.save_bleedthrough_diagnostics}")
    lines.append(f"Bleedthrough rules: {json.dumps(config.bleedthrough_rules, ensure_ascii=False)}")
    lines.append(f"DAPI channel: {config.dapi_channel}")
    lines.append(f"Intensity channels: {config.intensity_channels}")
    lines.append(f"Auto colocalization pairs: {config.auto_colocalization_pairs}")
    lines.append(f"Colocalization pairs: {config.colocalization_pairs}")
    lines.append(f"Channel colors: {json.dumps(config.channel_colors, ensure_ascii=False)}")
    lines.append(f"Raw channel display mode: {config.raw_channel_display_mode}")
    lines.append(f"Positive preview display mode: {config.positive_preview_display_mode}")
    lines.append(f"Positive preview colormap: {config.positive_preview_colormap}")
    lines.append(f"Filename role map: {json.dumps(config.filename_role_map, ensure_ascii=False)}")
    lines.append(f"Filename channel map: {json.dumps(config.filename_channel_map, ensure_ascii=False)}")
    lines.append("")
    lines.append("Statistics policy")
    lines.append("- Normality test: Shapiro-Wilk when n >= 3 and data are non-constant")
    lines.append("- Group test: Welch t-test if both groups pass normality; otherwise Mann-Whitney U")
    lines.append("- Multiple testing correction: Benjamini-Hochberg FDR")
    lines.append("- SD: sample SD (ddof=1)")
    lines.append("")

    if not per_image_df.empty and "resolved_channel_order" in per_image_df.columns:
        lines.append("Per-image channel mapping")
        lines.append("-" * 72)
        for _, row in per_image_df[["group", "image_file", "image_output_id", "input_mode", "channel_detection_method", "resolved_channel_order"]].iterrows():
            lines.append(
                f"[{row['group']}] {row['image_file']} ({row['image_output_id']}): "
                f"mode={row['input_mode']} | method={row['channel_detection_method']} | {row['resolved_channel_order']}"
            )
        lines.append("")

    if not stats_df.empty:
        lines.append("Metric summary")
        lines.append("-" * 72)
        for _, row in stats_df.iterrows():
            lines.append(
                f"{row['metric']}: "
                f"{config.group1_name}={row['mean_group1']:.4g}±{row['sd_group1']:.4g} (n={int(row['n_group1'])}), "
                f"{config.group2_name}={row['mean_group2']:.4g}±{row['sd_group2']:.4g} (n={int(row['n_group2'])}), "
                f"test={row['test']}, raw p={format_p_value(row['p_value_raw'])}, "
                f"FDR p={format_p_value(row['p_value_fdr_bh'])}, sig={row['significance_fdr']}, "
                f"{row['effect_size_name']}={row['effect_size']:.4g}"
            )
        lines.append("")

    if not errors_df.empty:
        lines.append("Skipped / failed images")
        lines.append("-" * 72)
        for _, row in errors_df.iterrows():
            lines.append(f"[{row['group']}] {row['image_path']}: {row['error']}")
        lines.append("")

    lines.append("Interpretation note")
    lines.append("- 当前脚本默认以“图像”为统计单位。若一只动物/一个样本包含多个视野，严格统计应先按生物学重复汇总，再做组间比较。")

    report_path.write_text("\n".join(lines), encoding="utf-8")


# ============================================================================
# 主流程
# ============================================================================


def run_analysis(config: AnalysisConfig) -> Dict[str, Any]:
    config = validate_config(config)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_tag = f"{sanitize_filename(config.group1_name)}_vs_{sanitize_filename(config.group2_name)}_{timestamp}"
    run_dir = config.output_dir / run_tag
    make_output_dir(run_dir)
    qc_dir = run_dir / "qc_overlays"
    diagnostics_dir = run_dir / "channel_diagnostics"
    merged_dir = run_dir / "merged_denoised"
    if config.save_qc_overlays:
        make_output_dir(qc_dir)
    if config.save_channel_diagnostics:
        make_output_dir(diagnostics_dir)
    make_output_dir(merged_dir)

    group1_inputs = discover_analysis_inputs(config.group1_dir, config)
    group2_inputs = discover_analysis_inputs(config.group2_dir, config)

    if len(group1_inputs) == 0:
        raise RuntimeError(f"组 {config.group1_name} 文件夹中未找到可分析输入: {config.group1_dir}")
    if len(group2_inputs) == 0:
        raise RuntimeError(f"组 {config.group2_name} 文件夹中未找到可分析输入: {config.group2_dir}")

    print("=" * 88)
    print("Easy IF batch analysis (fixed)")
    print("=" * 88)
    print(f"Group 1: {config.group1_name} | inputs={len(group1_inputs)} | folder={config.group1_dir}")
    print(f"Group 2: {config.group2_name} | inputs={len(group2_inputs)} | folder={config.group2_dir}")
    print(f"Input discovery mode   : {config.input_discovery_mode}")
    print(f"Manual FILE_CHANNEL_ORDER: {config.file_channel_order if config.file_channel_order else 'AUTO'}")
    print(f"Auto detect channels   : {config.auto_detect_channels}")
    print(f"Strict detection       : {config.strict_channel_detection}")
    print(f"Positive thresholds    : {json.dumps(config.positive_threshold_rules, ensure_ascii=False)}")
    print(f"Background pct default : {config.background_percentile}")
    print(f"Background pct rules   : {json.dumps(config.background_percentile_rules, ensure_ascii=False)}")
    print(f"Min positive obj area  : {config.min_positive_object_area}")
    print(f"Positive obj area rules: {json.dumps(config.positive_object_area_rules, ensure_ascii=False)}")
    print(f"Intensity channels     : {config.intensity_channels}")
    print(f"Auto coloc pairs       : {config.auto_colocalization_pairs}")
    print(f"Colocalization pairs   : {config.colocalization_pairs}")
    print(f"Channel colors         : {json.dumps(config.channel_colors, ensure_ascii=False)}")
    print(f"Raw channel display    : {config.raw_channel_display_mode}")
    print(f"Positive preview mode  : {config.positive_preview_display_mode}")
    print(f"Positive preview cmap  : {config.positive_preview_colormap}")
    print(f"Filename role map      : {json.dumps(config.filename_role_map, ensure_ascii=False)}")
    print(f"Filename channel map   : {json.dumps(config.filename_channel_map, ensure_ascii=False)}")
    print(f"Output                 : {run_dir}")

    all_rows: List[Dict[str, Any]] = []
    errors: List[Dict[str, str]] = []

    for group_name, input_list in ((config.group1_name, group1_inputs), (config.group2_name, group2_inputs)):
        print("-" * 88)
        print(f"Analyzing group: {group_name}")
        for index, analysis_input in enumerate(input_list, start=1):
            input_path = analysis_input.path
            try:
                if analysis_input.source_mode == "sample_folder":
                    result = analyze_sample_folder(
                        sample_dir=input_path,
                        group_name=group_name,
                        config=config,
                        qc_dir=qc_dir,
                        diagnostics_dir=diagnostics_dir,
                        merged_dir=merged_dir,
                        sample_contents=analysis_input.sample_contents,
                    )
                else:
                    result = analyze_single_image(
                        image_path=input_path,
                        group_name=group_name,
                        config=config,
                        qc_dir=qc_dir,
                        diagnostics_dir=diagnostics_dir,
                        merged_dir=merged_dir,
                    )
                all_rows.append(result.values)
                print(
                    f"[{group_name}] {index:>3d}/{len(input_list):<3d} {input_path.name} | "
                    f"mode={result.values.get('input_mode', 'NA')} | "
                    f"cells={result.values['cell_count']} | "
                    f"roi_area={result.values['roi_area_fraction']:.3f} | "
                    f"mapping={result.values['resolved_channel_order']}"
                )
            except Exception as exc:
                errors.append({
                    "group": group_name,
                    "image_path": str(input_path),
                    "error": str(exc),
                })
                print(f"[ERROR] [{group_name}] {input_path.name}: {exc}")

    per_image_df = pd.DataFrame(all_rows)
    errors_df = pd.DataFrame(errors)

    if per_image_df.empty:
        raise RuntimeError("没有任何图像成功完成分析，请检查输入图像格式和通道配置。")

    if (per_image_df["group"] == config.group1_name).sum() == 0 or (per_image_df["group"] == config.group2_name).sum() == 0:
        raise RuntimeError(
            "至少有一组没有成功分析到任何图像，无法进行组间统计。"
            " 请检查 failed_images.csv 中的错误信息。"
        )

    metrics_to_test = [col for col, _, _ in metric_specs(config) if col in per_image_df.columns]
    stats_rows = []
    for metric in metrics_to_test:
        g1 = per_image_df.loc[per_image_df["group"] == config.group1_name, metric].to_numpy()
        g2 = per_image_df.loc[per_image_df["group"] == config.group2_name, metric].to_numpy()
        stats_rows.append(compare_two_groups(metric, g1, g2))

    stats_df = pd.DataFrame(stats_rows)
    if not stats_df.empty:
        stats_df = apply_fdr_correction(stats_df)
        stats_df = stats_df.sort_values(by=["p_value_fdr_bh", "p_value_raw", "metric"], na_position="last").reset_index(drop=True)

    per_image_csv = run_dir / "per_image_results.csv"
    stats_csv = run_dir / "group_statistics.csv"
    error_csv = run_dir / "failed_images.csv"
    report_txt = run_dir / "analysis_report.txt"
    config_json = run_dir / "run_config.json"
    summary_base = run_dir / "summary_plots"

    per_image_df.to_csv(per_image_csv, index=False)
    stats_df.to_csv(stats_csv, index=False)
    errors_df.to_csv(error_csv, index=False)
    write_text_report(report_txt, config, group1_inputs, group2_inputs, per_image_df, stats_df, errors_df)
    config_json.write_text(json.dumps(asdict(config), ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    figure_paths = create_summary_figure(per_image_df, stats_df, config, summary_base)

    print("=" * 88)
    print("Analysis finished")
    print(f"Per-image results: {per_image_csv}")
    print(f"Group statistics : {stats_csv}")
    print(f"Error log        : {error_csv}")
    print(f"Report           : {report_txt}")
    if figure_paths:
        for key, value in figure_paths.items():
            print(f"Summary figure ({key}): {value}")
    print(f"Merged denoised   : {merged_dir}")
    if config.save_qc_overlays:
        print(f"QC overlays      : {qc_dir}")
    if config.save_channel_diagnostics:
        print(f"Channel diagnostics: {diagnostics_dir}")

    return {
        "run_dir": run_dir,
        "per_image_results": per_image_csv,
        "group_statistics": stats_csv,
        "failed_images": error_csv,
        "report": report_txt,
        "config": config_json,
        "figures": figure_paths,
        "per_image_df": per_image_df,
        "stats_df": stats_df,
        "errors_df": errors_df,
    }


# ============================================================================
# 命令行接口
# ============================================================================


def parse_coloc_pairs(raw_pairs: Sequence[str]) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for item in raw_pairs:
        parts = re.split(r"[:/,|]", item)
        parts = [x.strip() for x in parts if x.strip()]
        if len(parts) != 2:
            raise ValueError(f"共定位参数格式错误: {item}，请使用 488:594 这样的格式。")
        out.append((parts[0], parts[1]))
    return out


def parse_optional_channel_order(raw: Optional[Sequence[str]]) -> Optional[List[str]]:
    if raw is None:
        return None
    values = [str(x).strip() for x in raw if str(x).strip()]
    if not values:
        return None
    if len(values) == 1 and values[0].upper() in {"AUTO", "NONE", "NULL"}:
        return None
    return values


def parse_filename_channel_map(raw_items: Sequence[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in raw_items:
        parts = re.split(r"[:=]", str(item), maxsplit=1)
        parts = [x.strip() for x in parts if x.strip()]
        if len(parts) != 2:
            raise ValueError(f"文件名通道映射格式错误: {item}，请使用 CH1:DAPI 这样的格式。")
        out[parts[0]] = parts[1]
    return out


def parse_filename_role_map(raw_items: Sequence[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in raw_items:
        parts = re.split(r"[:=]", str(item), maxsplit=1)
        parts = [x.strip() for x in parts if x.strip()]
        if len(parts) != 2:
            raise ValueError(f"文件名角色映射格式错误: {item}，请使用 CH0:overlay 这样的格式。")
        out[parts[0]] = parts[1]
    return out


def parse_channel_int_map(raw_items: Sequence[str], item_name: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for item in raw_items:
        parts = re.split(r"[:=]", str(item), maxsplit=1)
        parts = [x.strip() for x in parts if x.strip()]
        if len(parts) != 2:
            raise ValueError(f"{item_name} 格式错误: {item}，请使用 488:100 这样的格式。")
        out[parts[0]] = int(float(parts[1]))
    return out


def build_default_config() -> AnalysisConfig:
    file_channel_order = globals().get("FILE_CHANNEL_ORDER", None)
    if file_channel_order is not None:
        file_channel_order = list(file_channel_order)
    bleedthrough_rules = build_bleedthrough_rules_from_simple_config(
        source_map=json.loads(json.dumps(BLEEDTHROUGH_SOURCE_MAP)),
        default_rule=json.loads(json.dumps(BLEEDTHROUGH_DEFAULT_RULE)),
        manual_coefficients=json.loads(json.dumps(BLEEDTHROUGH_MANUAL_COEFFICIENTS)),
        rule_overrides=json.loads(json.dumps(BLEEDTHROUGH_RULE_OVERRIDES)),
    )
    return AnalysisConfig(
        group1_dir=Path(GROUP1_DIR),
        group2_dir=Path(GROUP2_DIR),
        group1_name=GROUP1_NAME,
        group2_name=GROUP2_NAME,
        output_dir=Path(OUTPUT_DIR),

        input_discovery_mode=str(INPUT_DISCOVERY_MODE),
        file_channel_order=file_channel_order,
        auto_detect_channels=bool(AUTO_DETECT_CHANNELS),
        strict_channel_detection=bool(STRICT_CHANNEL_DETECTION),
        save_channel_diagnostics=bool(SAVE_CHANNEL_DIAGNOSTICS),
        fallback_channel_order=list(FALLBACK_CHANNEL_ORDER),
        positive_threshold_rules=json.loads(json.dumps(POSITIVE_THRESHOLD_RULES)),

        enable_bleedthrough_correction=bool(ENABLE_BLEEDTHROUGH_CORRECTION),
        save_bleedthrough_diagnostics=bool(SAVE_BLEEDTHROUGH_DIAGNOSTICS),
        bleedthrough_rules=bleedthrough_rules,

        dapi_channel=DAPI_CHANNEL,
        intensity_channels=list(INTENSITY_CHANNELS),
        auto_colocalization_pairs=bool(AUTO_COLOCALIZATION_PAIRS),
        colocalization_pairs=list(COLOCALIZATION_PAIRS),
        filename_role_map=json.loads(json.dumps(FILENAME_ROLE_MAP)),
        filename_channel_map=json.loads(json.dumps(FILENAME_CHANNEL_MAP)),

        image_patterns=list(IMAGE_PATTERNS),
        recursive_scan=RECURSIVE_SCAN,
        background_percentile=BACKGROUND_PERCENTILE,
        background_percentile_rules=json.loads(json.dumps(BACKGROUND_PERCENTILE_RULES)),
        min_nucleus_area=MIN_NUCLEUS_AREA,
        min_positive_object_area=MIN_POSITIVE_OBJECT_AREA,
        positive_object_area_rules=json.loads(json.dumps(POSITIVE_OBJECT_AREA_RULES)),
        gaussian_blur_size=GAUSSIAN_BLUR_SIZE,
        min_peak_distance=MIN_PEAK_DISTANCE,
        mask_dilation_radius=MASK_DILATION_RADIUS,
        channel_colors=json.loads(json.dumps(CHANNEL_COLORS)),
        raw_channel_display_mode=str(RAW_CHANNEL_DISPLAY_MODE),
        positive_preview_display_mode=str(POSITIVE_PREVIEW_DISPLAY_MODE),
        positive_preview_colormap=str(POSITIVE_PREVIEW_COLORMAP),
        save_qc_overlays=SAVE_QC_OVERLAYS,
        save_figure_pdf=SAVE_FIGURE_PDF,
        save_figure_png=SAVE_FIGURE_PNG,
        figure_dpi=FIGURE_DPI,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Easy batch immunofluorescence analysis (fixed + configurable bleed-through correction). 也可直接修改脚本顶部配置区后运行。"
    )
    parser.add_argument("--group1-dir", type=str, help="组1图片文件夹")
    parser.add_argument("--group2-dir", type=str, help="组2图片文件夹")
    parser.add_argument("--group1-name", type=str, help="组1名称")
    parser.add_argument("--group2-name", type=str, help="组2名称")
    parser.add_argument("--output-dir", type=str, help="输出文件夹")
    parser.add_argument("--input-discovery-mode", type=str, help="输入模式：auto / sample_folders / image_files")

    parser.add_argument("--file-channel-order", nargs="+", help="手动指定文件通道顺序；填 AUTO 表示自动检测")
    parser.add_argument("--fallback-channel-order", nargs="+", help="关闭 strict 时的兜底顺序")
    parser.add_argument("--no-auto-detect", action="store_true", help="关闭自动通道检测")
    parser.add_argument("--not-strict-channel-detection", action="store_true", help="允许自动检测失败后按 fallback 兜底")
    parser.add_argument("--no-channel-diagnostics", action="store_true", help="不保存原始通道诊断图")
    parser.add_argument("--no-bleedthrough-correction", action="store_true", help="关闭串色校正")
    parser.add_argument("--no-bleedthrough-diagnostics", action="store_true", help="不保存串色诊断图")
    parser.add_argument("--bleedthrough-488-manual", type=float, help="手动指定 488 <- DAPI 的扣除系数，会把 488 规则改为 manual 模式")

    parser.add_argument("--dapi-channel", type=str, help="DAPI 通道名")
    parser.add_argument("--intensity-channels", nargs="+", help="做强度分析的通道")
    parser.add_argument("--coloc-pairs", nargs="+", help="共定位通道对，如 488:594 594:647")
    parser.add_argument("--filename-role-map", nargs="+", help="文件名关键字与文件角色映射，如 CH0:overlay Merge:ignore")
    parser.add_argument("--filename-channel-map", nargs="+", help="文件名关键字与通道名映射，如 CH1:DAPI CH2:488")
    parser.add_argument("--all-coloc-pairs", action="store_true", help="自动对所有 intensity channels 生成两两共定位/相关性统计")
    parser.add_argument("--no-auto-coloc-pairs", action="store_true", help="关闭自动两两配对，仅使用 --coloc-pairs 或顶部配置中的手动列表")
    parser.add_argument("--recursive-scan", action="store_true", help="递归扫描子文件夹")
    parser.add_argument("--background-percentile", type=float, help="背景分位数")
    parser.add_argument("--min-nucleus-area", type=int, help="最小核面积")
    parser.add_argument("--min-positive-object-area", type=int, help="阈值后最小阳性区域面积；用于去除 marker 碎点")
    parser.add_argument("--positive-object-area-rules", nargs="+", help="按通道覆盖阳性碎点面积，如 488:50 594:120")
    parser.add_argument("--gaussian-blur-size", type=int, help="DAPI 平滑核大小")
    parser.add_argument("--min-peak-distance", type=int, help="watershed 种子最小距离")
    parser.add_argument("--mask-dilation-radius", type=int, help="DAPI 掩膜扩张像素")
    parser.add_argument("--raw-channel-display-mode", type=str, help="原始单通道显示方式：gray / color")
    parser.add_argument("--positive-preview-display-mode", type=str, help="阳性预览显示方式：colormap / color / gray")
    parser.add_argument("--positive-preview-colormap", type=str, help="阳性预览 colormap，如 inferno / magma / gray")
    parser.add_argument("--no-qc", action="store_true", help="不保存 QC 叠加图")
    return parser


def apply_cli_overrides(config: AnalysisConfig, args: argparse.Namespace) -> AnalysisConfig:
    if args.group1_dir:
        config.group1_dir = Path(args.group1_dir)
    if args.group2_dir:
        config.group2_dir = Path(args.group2_dir)
    if args.group1_name:
        config.group1_name = args.group1_name
    if args.group2_name:
        config.group2_name = args.group2_name
    if args.output_dir:
        config.output_dir = Path(args.output_dir)
    if args.input_discovery_mode:
        config.input_discovery_mode = str(args.input_discovery_mode)

    cli_order = parse_optional_channel_order(args.file_channel_order)
    if args.file_channel_order is not None:
        config.file_channel_order = cli_order

    if args.fallback_channel_order:
        config.fallback_channel_order = list(args.fallback_channel_order)
    if args.no_auto_detect:
        config.auto_detect_channels = False
    if args.not_strict_channel_detection:
        config.strict_channel_detection = False
    if args.no_channel_diagnostics:
        config.save_channel_diagnostics = False
    if args.dapi_channel:
        config.dapi_channel = args.dapi_channel
    if args.intensity_channels:
        config.intensity_channels = list(args.intensity_channels)
    if args.no_bleedthrough_correction:
        config.enable_bleedthrough_correction = False
    if args.no_bleedthrough_diagnostics:
        config.save_bleedthrough_diagnostics = False
    if args.bleedthrough_488_manual is not None:
        config.bleedthrough_rules.setdefault("488", {"source": config.dapi_channel})
        config.bleedthrough_rules["488"]["source"] = config.dapi_channel
        config.bleedthrough_rules["488"]["mode"] = "manual"
        config.bleedthrough_rules["488"]["coefficient"] = float(args.bleedthrough_488_manual)
    if args.filename_role_map:
        config.filename_role_map = parse_filename_role_map(args.filename_role_map)
    if args.filename_channel_map:
        config.filename_channel_map = parse_filename_channel_map(args.filename_channel_map)
    if args.all_coloc_pairs:
        config.auto_colocalization_pairs = True
        config.colocalization_pairs = []
    if args.no_auto_coloc_pairs:
        config.auto_colocalization_pairs = False
    if args.coloc_pairs:
        config.auto_colocalization_pairs = False
        config.colocalization_pairs = parse_coloc_pairs(args.coloc_pairs)
    if args.recursive_scan:
        config.recursive_scan = True
    if args.background_percentile is not None:
        config.background_percentile = float(args.background_percentile)
    if args.min_nucleus_area is not None:
        config.min_nucleus_area = int(args.min_nucleus_area)
    if args.min_positive_object_area is not None:
        config.min_positive_object_area = int(args.min_positive_object_area)
    if args.positive_object_area_rules:
        config.positive_object_area_rules = parse_channel_int_map(
            args.positive_object_area_rules,
            item_name="positive object area rules",
        )
    if args.gaussian_blur_size is not None:
        config.gaussian_blur_size = int(args.gaussian_blur_size)
    if args.min_peak_distance is not None:
        config.min_peak_distance = int(args.min_peak_distance)
    if args.mask_dilation_radius is not None:
        config.mask_dilation_radius = int(args.mask_dilation_radius)
    if args.raw_channel_display_mode:
        config.raw_channel_display_mode = str(args.raw_channel_display_mode)
    if args.positive_preview_display_mode:
        config.positive_preview_display_mode = str(args.positive_preview_display_mode)
    if args.positive_preview_colormap:
        config.positive_preview_colormap = str(args.positive_preview_colormap)
    if args.no_qc:
        config.save_qc_overlays = False
    return config


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    config = build_default_config()
    config = apply_cli_overrides(config, args)
    run_analysis(config)


if __name__ == "__main__":
    main()
