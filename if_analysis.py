#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Easy batch immunofluorescence (IF) analysis script - fixed version.

这版重点修复了原脚本中的两个核心问题：

1) 通道识别
   - 原脚本对 TIFF 直接按 FILE_CHANNEL_ORDER 绑定，容易把 RGB TIFF / 多通道 TIFF 的通道顺序搞错。
   - 现在优先使用：
       a. 手动 FILE_CHANNEL_ORDER（如果你明确知道顺序）
       b. TIFF/OME 元数据中的通道名
       c. RGB/BGR 合成图的颜色语义自动映射
          - cv2 读 PNG/JPG/BMP: BGR -> [DAPI, 488, 594]
          - tifffile 读 RGB TIFF  : RGB -> [594, 488, DAPI]
       d. 仅对 DAPI 做启发式自动识别
   - 如果是“没有元数据的原始灰度平面堆栈”，脚本无法可靠地区分 488 和 594。
     这种情况下默认严格报错，避免静默分析错通道。

2) QC 输出
   - 现在会保存：
       a. 标准 QC 图（含 DAPI raw、ROI QC、首个分析通道 raw / bg-subtracted、Merged QC）
       b. channel_diagnostics：每个原始 index 的单独 raw 图 + 自动映射结果
   - 这样能直接审查“原始第 0/1/2 通道分别是什么”。

为什么不能单纯“按颜色自动筛选通道”？
-----------------------------------
因为显微图像里常见的是两种完全不同的数据：

A. 原始多通道 TIFF / OME-TIFF
   - 底层通常是多个灰度平面（plane 0, plane 1, plane 2 ...）
   - “蓝/绿/红”只是软件显示时人为指定的伪彩
   - 文件本身很多时候并没有把“这是 DAPI / 这是 488 / 这是 594”写死到像素里
   - 如果没有 metadata，就无法只靠“颜色”知道哪个平面是哪个 marker

B. 已经合成好的 RGB/BGR 彩图
   - 这类图像才真的有 R/G/B 三个颜色通道
   - 这时可以按蓝= DAPI、绿=488、红=594 这样映射
   - 但如果 594 和 647 已经一起混进红色，就再也无法分开

所以：
- “颜色自动检测”只适用于 RGB/BGR 合成图
- “原始多通道灰度堆栈”必须依赖 metadata，或者手动给顺序
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
# 用户配置区：绝大多数情况下，只需要修改这一段
# ============================================================================

GROUP1_DIR = r"/Users/yemingzhu/Downloads/课题文件/原始数据胸腺/免疫荧光/all_in_one/all_in_one_analysis/p21/3m"
GROUP2_DIR = r"/Users/yemingzhu/Downloads/课题文件/原始数据胸腺/免疫荧光/all_in_one/all_in_one_analysis/p21/22m"
GROUP1_NAME = "3m"
GROUP2_NAME = "22m"
OUTPUT_DIR = r"/Users/yemingzhu/Downloads/课题文件/原始数据胸腺/免疫荧光/if_analysis_results"

# 手动指定文件中的通道顺序（按读入后的 index 顺序写）
# - 如果你明确知道真实顺序，直接填写，例如：
#     FILE_CHANNEL_ORDER = ["594", "488", "DAPI"]
# - 如果不确定，设为 None，让脚本自动检测 metadata / RGB 颜色语义
FILE_CHANNEL_ORDER: Optional[List[str]] = None #  如果染色为 647和 488，这里要改一下["647", "488", "DAPI"]

# 自动通道检测
AUTO_DETECT_CHANNELS = True
STRICT_CHANNEL_DETECTION = True  # 模糊时直接报错，避免静默错配
SAVE_CHANNEL_DIAGNOSTICS = True  # 保存每个 raw index 的诊断图
FALLBACK_CHANNEL_ORDER = ["DAPI", "488", "594", "647"]  # 仅在关闭 strict 时作为兜底顺序

# 阳性阈值
# 作用：
# - 控制 positive_area_fraction / 共定位等“什么算阳性”
# - 同时影响 channel_diagnostics 和 QC 图里的阈值化预览
# 默认：
# - 每个通道都使用 ROI 内 bg-subtracted 信号的 Otsu 自动阈值
# 如果 594 假阳性偏多，可提高阈值，例如：
# POSITIVE_THRESHOLD_RULES = {
#     "594": {"method": "otsu", "scale": 1.8, "min_value": 30.0}
# }
# 可选 method:
# - "otsu"      : threshold = otsu * scale + offset，再与 min_value 取较大
# - "manual"    : threshold = value
# - "percentile": threshold = ROI 内 percentile 分位数，再叠加 scale / offset / min_value
#POSITIVE_THRESHOLD_RULES: Dict[str, Dict[str, Any]] = {}  #空字典

# 修改594阈值，建议采用 Otsu 自动阈值乘以 1.8 倍缩放（scale）并设定 30 为下限，使 594 信号避免将弱噪声误判为阳性，
# 若噪声仍偏多可上调 scale 至 2.0–2.2，若真信号被过度压制则下调至 1.4–1.6，亦可直接通过 {"method": "manual", "value": 80} 完全手动指定阈值。
POSITIVE_THRESHOLD_RULES = {
    "594": {"method": "otsu", "scale": 1.6, "min_value": 30.0}  
}


# 串色 / bleed-through 校正
# 适用场景：
# - 现在通道顺序已经正确
# - 但 488 原始 plane 中仍混入了 DAPI 结构（看起来像 "DAPI + 488"）
# 核心思想：
#   corrected_488 = max((488 - bg488) - k * (DAPI - bgDAPI), 0)
# 这不是“按颜色筛选”，而是按线性串色模型扣除 DAPI 对 488 的贡献。
ENABLE_BLEEDTHROUGH_CORRECTION = True
SAVE_BLEEDTHROUGH_DIAGNOSTICS = True

BLEEDTHROUGH_RULES: Dict[str, Dict[str, Any]] = {
    "488": {
        "source": "DAPI",                # 从哪个通道扣除串色
        "mode": "auto",                  # "auto" 或 "manual"
        "coefficient": 0.18,             # mode="manual" 时使用；auto 模式下会忽略
        "estimate_mask": "nuclei",       # "nuclei" / "roi" / "all"
        "ratio_percentile": 20.0,        # auto 模式：取较低分位，尽量避免过度扣除
        "source_threshold_percentile": 75.0,  # auto 模式：只在 source 较强像素上估计
        "min_pixels": 100,               # auto 模式：最少用于估计的像素数
        "max_coefficient": 3.0,          # auto 模式：k 的上限，避免异常值
    }
}

# DAPI 及需要分析的通道
DAPI_CHANNEL = "DAPI"
INTENSITY_CHANNELS = ["488", "594"]  # 如果要同时分析 488/594/647，可改成 ["488", "594", "647"]
AUTO_COLOCALIZATION_PAIRS = True     # True: 自动对 INTENSITY_CHANNELS 生成全部两两组合
COLOCALIZATION_PAIRS: List[Tuple[str, str]] = []  # 仅在 AUTO_COLOCALIZATION_PAIRS=False 时生效

# 文件扫描
IMAGE_PATTERNS = ["*.tif", "*.tiff", "*.png", "*.jpg", "*.jpeg", "*.bmp"]
RECURSIVE_SCAN = False

# 分析参数
BACKGROUND_PERCENTILE = 3   # 背景估计分位数；值越大，通常扣掉的背景越多,60X 的话可填 5
MIN_NUCLEUS_AREA = 1     # 细胞核最小面积；小于这个值的区域会被当作噪声去掉，60X的话可填 2200
GAUSSIAN_BLUR_SIZE = 7      # 分割前的高斯模糊核大小；越大越平滑，但也更容易丢细节
MIN_PEAK_DISTANCE = 0.2      # 分水岭分割时种子点的最小距离；越大越不容易把挨得近的核分开 ，60X 的话可填 25
MASK_DILATION_RADIUS = 0.4   # 以细胞核为基础向外扩张的像素半径，用来近似细胞 ROI ，60X 的话可填 18
SAVE_QC_OVERLAYS = True     # 是否保存 QC 质控图；True 保存，False 不保存

# 输出
SAVE_FIGURE_PDF = True
SAVE_FIGURE_PNG = True
FIGURE_DPI = 500


# ============================================================================
# 一般不需要改动的参数
# ============================================================================

SIGNIFICANCE_LEVELS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]
DEFAULT_GROUP_COLORS = ("#4472C4", "#ED7D31")
MIN_PIXELS_FOR_COLOC = 10
ALLOWED_CHANNEL_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_\-+.]+$")
CANONICAL_CHANNELS = ("DAPI", "488", "594", "647")


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

    image_patterns: List[str] = field(default_factory=lambda: ["*.tif", "*.tiff", "*.png", "*.jpg", "*.jpeg", "*.bmp"])
    recursive_scan: bool = False
    background_percentile: float = 10.0
    min_nucleus_area: int = 50
    gaussian_blur_size: int = 5
    min_peak_distance: int = 6
    mask_dilation_radius: int = 6
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
    arr = np.asarray(signal, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    threshold = float(max(0.0, threshold))
    out = arr.copy()
    out[out <= threshold] = 0.0
    return out


def make_positive_mask(signal: np.ndarray, threshold: float) -> np.ndarray:
    arr = np.asarray(signal, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    threshold = float(max(0.0, threshold))
    return arr > threshold


def positive_fraction_in_roi(signal: np.ndarray, roi_mask: np.ndarray, threshold: float) -> float:
    roi_mask = np.asarray(roi_mask, dtype=bool)
    positive_mask = make_positive_mask(signal, threshold)
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


def build_merged_preview(
    channel_images: Dict[str, np.ndarray],
    positive_thresholds: Optional[Dict[str, PositiveThresholdInfo]] = None,
    thresholded: bool = False,
) -> np.ndarray:
    if not channel_images:
        return np.zeros((1, 1, 3), dtype=np.uint8)

    sample = next(iter(channel_images.values()))
    sample = np.asarray(sample)
    merged = np.zeros((sample.shape[0], sample.shape[1], 3), dtype=np.uint8)

    dapi = channel_images.get("DAPI")
    if dapi is not None:
        merged[..., 2] = np.maximum(merged[..., 2], robust_normalize_for_display(dapi))

    if "488" in channel_images:
        green = channel_images["488"]
        if thresholded and positive_thresholds and "488" in positive_thresholds:
            green = make_thresholded_preview(green, positive_thresholds["488"].threshold)
        merged[..., 1] = np.maximum(merged[..., 1], robust_normalize_for_display(green))

    if "594" in channel_images:
        ch594 = channel_images["594"]
        if thresholded and positive_thresholds and "594" in positive_thresholds:
            ch594 = make_thresholded_preview(ch594, positive_thresholds["594"].threshold)
        merged[..., 0] = np.maximum(merged[..., 0], robust_normalize_for_display(ch594))

    if "647" in channel_images:
        ch647 = channel_images["647"]
        if thresholded and positive_thresholds and "647" in positive_thresholds:
            ch647 = make_thresholded_preview(ch647, positive_thresholds["647"].threshold)
        ch647_u8 = robust_normalize_for_display(ch647)
        merged[..., 0] = np.maximum(merged[..., 0], ch647_u8)
        merged[..., 2] = np.maximum(merged[..., 2], np.clip(np.round(ch647_u8 * 0.7), 0, 255).astype(np.uint8))

    return merged


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
            raise ValueError(
                f"{image_path.name} 共有 {n_channels} 个通道，但 FILE_CHANNEL_ORDER 配置了 {len(manual)} 个。"
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
    if config.file_channel_order:
        config.file_channel_order = [normalize_channel_name(x) for x in config.file_channel_order]
        if len(config.file_channel_order) != len(set(config.file_channel_order)):
            raise ValueError(f"FILE_CHANNEL_ORDER 中存在重复通道名: {config.file_channel_order}")
        for name in config.file_channel_order:
            if not ALLOWED_CHANNEL_NAME_PATTERN.match(name):
                raise ValueError(f"不支持的通道名称格式: {name}")

    config.fallback_channel_order = [normalize_channel_name(x) for x in config.fallback_channel_order]
    config.positive_threshold_rules = _normalize_positive_threshold_rules(config.positive_threshold_rules)
    config.bleedthrough_rules = _normalize_bleedthrough_rules(config.bleedthrough_rules)
    config.dapi_channel = normalize_channel_name(config.dapi_channel)
    config.intensity_channels = [normalize_channel_name(x) for x in config.intensity_channels]
    config.colocalization_pairs = normalize_colocalization_pairs(config.colocalization_pairs)

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
    config.min_peak_distance = int(max(1, config.min_peak_distance))
    config.mask_dilation_radius = int(max(0, config.mask_dilation_radius))
    config.background_percentile = float(np.clip(config.background_percentile, 0.0, 100.0))

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
    bleedthrough_rules: Dict[str, Dict[str, Any]],
) -> Tuple[np.ndarray, float, BleedthroughCorrectionInfo, np.ndarray, np.ndarray]:
    target_name = normalize_channel_name(target_channel_name)
    target_raw = np.asarray(raw_channels[target_name], dtype=np.float64)
    target_signal, target_bg = prepare_signal(target_raw, roi_mask=roi_mask, background_percentile=background_percentile)

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
    source_signal, source_bg = prepare_signal(source_raw, roi_mask=roi_mask, background_percentile=background_percentile)

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
        positive_area_fraction = float(np.mean(roi_values > positive_threshold.threshold))
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
) -> Dict[str, float]:
    a_name = normalize_channel_name(channel_a_name)
    b_name = normalize_channel_name(channel_b_name)
    signal_a = np.asarray(signal_a, dtype=np.float64)
    signal_b = np.asarray(signal_b, dtype=np.float64)
    roi_mask = np.asarray(roi_mask, dtype=bool)

    if not np.any(roi_mask):
        roi_mask = np.ones(signal_a.shape, dtype=bool)

    mask_a = roi_mask & (signal_a > threshold_a.threshold)
    mask_b = roi_mask & (signal_b > threshold_b.threshold)
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
    raw_merged = build_merged_preview(raw_channels, thresholded=False)
    denoised_merged = build_merged_preview(signals, positive_thresholds=positive_thresholds, thresholded=True)

    top_row: List[Tuple[str, np.ndarray]] = [("Raw merged", raw_merged)]
    bottom_row: List[Tuple[str, np.ndarray]] = [("Thresholded merged", denoised_merged)]

    for idx in range(n_channels):
        name = resolution.channel_order[idx]
        top_row.append((f"idx {idx}\n{name}", robust_normalize_for_display(image[..., idx])))

        threshold_info = positive_thresholds.get(name)
        if name in signals and threshold_info is not None:
            preview = make_thresholded_preview(signals[name], threshold_info.threshold)
            positive_fraction = positive_fraction_in_roi(signals[name], roi_mask, threshold_info.threshold)
            bottom_row.append(
                (
                    f"{name} positive preview\nthr={threshold_info.threshold:.4g} | roi+={positive_fraction:.2%}",
                    apply_colormap_to_u8(robust_normalize_for_display(preview), cmap_name="inferno"),
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
        (f"{config.dapi_channel} raw (idx {dapi_idx})", dapi_u8)
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
        preview = (
            make_thresholded_preview(signal, threshold_info.threshold)
            if threshold_info is not None
            else signal
        )
        positive_fraction = (
            positive_fraction_in_roi(signal, roi_mask, threshold_info.threshold)
            if threshold_info is not None
            else float("nan")
        )

        channel_idx = resolution.channel_to_index.get(channel_name, -1)
        qc_panels.append((f"{channel_name} raw (idx {channel_idx})", robust_normalize_for_display(raw)))
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
        qc_panels.append((preview_title, apply_colormap_to_u8(robust_normalize_for_display(preview), cmap_name="inferno")))
        merged_channels[channel_name] = preview

    qc_panels.append(("Merged QC (thresholded)", build_merged_preview(merged_channels, thresholded=False)))

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
    diagnostics_dir: Path,
) -> Optional[Path]:
    if not info.applied:
        return None

    out_dir = diagnostics_dir / "bleedthrough"
    make_output_dir(out_dir)

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    axes[0].imshow(robust_normalize_for_display(raw_channel), cmap="gray")
    axes[0].set_title(f"{channel_name} raw")
    axes[0].axis("off")

    axes[1].imshow(robust_normalize_for_display(uncorrected_signal), cmap="gray")
    axes[1].set_title(f"{channel_name} bg-subtracted")
    axes[1].axis("off")

    axes[2].imshow(robust_normalize_for_display(source_contribution), cmap="gray")
    axes[2].set_title(f"{info.source_channel} contribution\n(k={info.coefficient:.3g})")
    axes[2].axis("off")

    axes[3].imshow(robust_normalize_for_display(corrected_signal), cmap="gray")
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


def analyze_single_image(
    image_path: Path,
    group_name: str,
    config: AnalysisConfig,
    qc_dir: Path,
    diagnostics_dir: Path,
) -> ImageAnalysisResult:
    load_result = load_multichannel_image(image_path)
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
                bleedthrough_rules=bleedthrough_rules,
            )
        else:
            signal, bg = prepare_signal(
                channel_arrays[channel_name],
                roi_mask=roi_mask,
                background_percentile=config.background_percentile,
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
            )
        )

    for channel_a, channel_b in config.colocalization_pairs:
        result.update(
            compute_colocalization(
                channel_a_name=channel_a,
                signal_a=signals[channel_a],
                channel_b_name=channel_b,
                signal_b=signals[channel_b],
                roi_mask=roi_mask,
                threshold_a=positive_thresholds[channel_a],
                threshold_b=positive_thresholds[channel_b],
            )
        )

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
                diagnostics_dir=diagnostics_dir,
            )
            if bt_path is not None:
                bleedthrough_paths.append(str(bt_path))
    result["bleedthrough_diagnostic_paths"] = " | ".join(bleedthrough_paths)

    return ImageAnalysisResult(values=result, roi_mask=roi_mask, nuclei_labels=nuclei_labels)


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
    group1_images: List[Path],
    group2_images: List[Path],
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
    lines.append(f"Images found: {config.group1_name}={len(group1_images)}, {config.group2_name}={len(group2_images)}")
    lines.append(f"Successful analyses: {config.group1_name}={(per_image_df['group'] == config.group1_name).sum()}, {config.group2_name}={(per_image_df['group'] == config.group2_name).sum()}")
    lines.append(f"Manual FILE_CHANNEL_ORDER: {config.file_channel_order if config.file_channel_order else 'AUTO'}")
    lines.append(f"AUTO_DETECT_CHANNELS: {config.auto_detect_channels}")
    lines.append(f"STRICT_CHANNEL_DETECTION: {config.strict_channel_detection}")
    lines.append(f"Fallback channel order: {config.fallback_channel_order}")
    lines.append(f"Positive threshold rules: {json.dumps(config.positive_threshold_rules, ensure_ascii=False)}")
    lines.append(f"ENABLE_BLEEDTHROUGH_CORRECTION: {config.enable_bleedthrough_correction}")
    lines.append(f"SAVE_BLEEDTHROUGH_DIAGNOSTICS: {config.save_bleedthrough_diagnostics}")
    lines.append(f"Bleedthrough rules: {json.dumps(config.bleedthrough_rules, ensure_ascii=False)}")
    lines.append(f"DAPI channel: {config.dapi_channel}")
    lines.append(f"Intensity channels: {config.intensity_channels}")
    lines.append(f"Auto colocalization pairs: {config.auto_colocalization_pairs}")
    lines.append(f"Colocalization pairs: {config.colocalization_pairs}")
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
        for _, row in per_image_df[["group", "image_file", "image_output_id", "channel_detection_method", "resolved_channel_order"]].iterrows():
            lines.append(
                f"[{row['group']}] {row['image_file']} ({row['image_output_id']}): "
                f"method={row['channel_detection_method']} | {row['resolved_channel_order']}"
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
    if config.save_qc_overlays:
        make_output_dir(qc_dir)
    if config.save_channel_diagnostics:
        make_output_dir(diagnostics_dir)

    group1_images = discover_images(config.group1_dir, config.image_patterns, recursive=config.recursive_scan)
    group2_images = discover_images(config.group2_dir, config.image_patterns, recursive=config.recursive_scan)

    if len(group1_images) == 0:
        raise RuntimeError(f"组 {config.group1_name} 文件夹中未找到图像: {config.group1_dir}")
    if len(group2_images) == 0:
        raise RuntimeError(f"组 {config.group2_name} 文件夹中未找到图像: {config.group2_dir}")

    print("=" * 88)
    print("Easy IF batch analysis (fixed)")
    print("=" * 88)
    print(f"Group 1: {config.group1_name} | images={len(group1_images)} | folder={config.group1_dir}")
    print(f"Group 2: {config.group2_name} | images={len(group2_images)} | folder={config.group2_dir}")
    print(f"Manual FILE_CHANNEL_ORDER: {config.file_channel_order if config.file_channel_order else 'AUTO'}")
    print(f"Auto detect channels   : {config.auto_detect_channels}")
    print(f"Strict detection       : {config.strict_channel_detection}")
    print(f"Positive thresholds    : {json.dumps(config.positive_threshold_rules, ensure_ascii=False)}")
    print(f"Intensity channels     : {config.intensity_channels}")
    print(f"Auto coloc pairs       : {config.auto_colocalization_pairs}")
    print(f"Colocalization pairs   : {config.colocalization_pairs}")
    print(f"Output                 : {run_dir}")

    all_rows: List[Dict[str, Any]] = []
    errors: List[Dict[str, str]] = []

    for group_name, image_list in ((config.group1_name, group1_images), (config.group2_name, group2_images)):
        print("-" * 88)
        print(f"Analyzing group: {group_name}")
        for index, image_path in enumerate(image_list, start=1):
            try:
                result = analyze_single_image(
                    image_path=image_path,
                    group_name=group_name,
                    config=config,
                    qc_dir=qc_dir,
                    diagnostics_dir=diagnostics_dir,
                )
                all_rows.append(result.values)
                print(
                    f"[{group_name}] {index:>3d}/{len(image_list):<3d} {image_path.name} | "
                    f"cells={result.values['cell_count']} | "
                    f"roi_area={result.values['roi_area_fraction']:.3f} | "
                    f"mapping={result.values['resolved_channel_order']}"
                )
            except Exception as exc:
                errors.append({
                    "group": group_name,
                    "image_path": str(image_path),
                    "error": str(exc),
                })
                print(f"[ERROR] [{group_name}] {image_path.name}: {exc}")

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
    write_text_report(report_txt, config, group1_images, group2_images, per_image_df, stats_df, errors_df)
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


def build_default_config() -> AnalysisConfig:
    return AnalysisConfig(
        group1_dir=Path(GROUP1_DIR),
        group2_dir=Path(GROUP2_DIR),
        group1_name=GROUP1_NAME,
        group2_name=GROUP2_NAME,
        output_dir=Path(OUTPUT_DIR),

        file_channel_order=list(FILE_CHANNEL_ORDER) if FILE_CHANNEL_ORDER else None,
        auto_detect_channels=bool(AUTO_DETECT_CHANNELS),
        strict_channel_detection=bool(STRICT_CHANNEL_DETECTION),
        save_channel_diagnostics=bool(SAVE_CHANNEL_DIAGNOSTICS),
        fallback_channel_order=list(FALLBACK_CHANNEL_ORDER),
        positive_threshold_rules=json.loads(json.dumps(POSITIVE_THRESHOLD_RULES)),

        enable_bleedthrough_correction=bool(ENABLE_BLEEDTHROUGH_CORRECTION),
        save_bleedthrough_diagnostics=bool(SAVE_BLEEDTHROUGH_DIAGNOSTICS),
        bleedthrough_rules=json.loads(json.dumps(BLEEDTHROUGH_RULES)),

        dapi_channel=DAPI_CHANNEL,
        intensity_channels=list(INTENSITY_CHANNELS),
        auto_colocalization_pairs=bool(AUTO_COLOCALIZATION_PAIRS),
        colocalization_pairs=list(COLOCALIZATION_PAIRS),

        image_patterns=list(IMAGE_PATTERNS),
        recursive_scan=RECURSIVE_SCAN,
        background_percentile=BACKGROUND_PERCENTILE,
        min_nucleus_area=MIN_NUCLEUS_AREA,
        gaussian_blur_size=GAUSSIAN_BLUR_SIZE,
        min_peak_distance=MIN_PEAK_DISTANCE,
        mask_dilation_radius=MASK_DILATION_RADIUS,
        save_qc_overlays=SAVE_QC_OVERLAYS,
        save_figure_pdf=SAVE_FIGURE_PDF,
        save_figure_png=SAVE_FIGURE_PNG,
        figure_dpi=FIGURE_DPI,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Easy batch immunofluorescence analysis (fixed + DAPI bleed-through correction). 也可直接修改脚本顶部配置区后运行。"
    )
    parser.add_argument("--group1-dir", type=str, help="组1图片文件夹")
    parser.add_argument("--group2-dir", type=str, help="组2图片文件夹")
    parser.add_argument("--group1-name", type=str, help="组1名称")
    parser.add_argument("--group2-name", type=str, help="组2名称")
    parser.add_argument("--output-dir", type=str, help="输出文件夹")

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
    parser.add_argument("--all-coloc-pairs", action="store_true", help="自动对所有 intensity channels 生成两两共定位/相关性统计")
    parser.add_argument("--no-auto-coloc-pairs", action="store_true", help="关闭自动两两配对，仅使用 --coloc-pairs 或顶部配置中的手动列表")
    parser.add_argument("--recursive-scan", action="store_true", help="递归扫描子文件夹")
    parser.add_argument("--background-percentile", type=float, help="背景分位数")
    parser.add_argument("--min-nucleus-area", type=int, help="最小核面积")
    parser.add_argument("--gaussian-blur-size", type=int, help="DAPI 平滑核大小")
    parser.add_argument("--min-peak-distance", type=int, help="watershed 种子最小距离")
    parser.add_argument("--mask-dilation-radius", type=int, help="DAPI 掩膜扩张像素")
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
    if args.no_bleedthrough_correction:
        config.enable_bleedthrough_correction = False
    if args.no_bleedthrough_diagnostics:
        config.save_bleedthrough_diagnostics = False
    if args.bleedthrough_488_manual is not None:
        config.bleedthrough_rules.setdefault("488", {"source": config.dapi_channel})
        config.bleedthrough_rules["488"]["source"] = config.dapi_channel
        config.bleedthrough_rules["488"]["mode"] = "manual"
        config.bleedthrough_rules["488"]["coefficient"] = float(args.bleedthrough_488_manual)

    if args.dapi_channel:
        config.dapi_channel = args.dapi_channel
    if args.intensity_channels:
        config.intensity_channels = list(args.intensity_channels)
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
    if args.gaussian_blur_size is not None:
        config.gaussian_blur_size = int(args.gaussian_blur_size)
    if args.min_peak_distance is not None:
        config.min_peak_distance = int(args.min_peak_distance)
    if args.mask_dilation_radius is not None:
        config.mask_dilation_radius = int(args.mask_dilation_radius)
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
