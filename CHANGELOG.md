# Changelog

## [1.0.7] - 2025-03-18

### 新增功能

- **样本文件夹模式增强**
  - 支持从样本文件夹（sample folders）自动组装多通道图像
  - 每个通道单独文件，自动识别并合并
  - 支持文件名映射（FILENAME_CHANNEL_MAP）处理非标准命名
  - 自动忽略 overlay/merge/composite/RGB 合成图

- **文件名角色映射**
  - 新增 FILENAME_ROLE_MAP 配置，可标记 overlay/merge 文件为忽略
  - 避免合成图干扰单通道文件识别

- **递归扫描选项**
  - 新增 RECURSIVE_SCAN 参数，支持递归扫描子目录

- **显示模式自定义**
  - RAW_CHANNEL_DISPLAY_MODE: 原始通道显示（color/gray）
  - POSITIVE_PREVIEW_DISPLAY_MODE: 阳性预览显示（colormap/color/gray）
  - POSITIVE_PREVIEW_COLORMAP: 自定义阳性预览色图

### 改进

- **通道识别进一步优化**
  - 增强 DAPI 启发式评分算法
  - 改进 RGB/BGR 合成图的通道语义映射

- **核分割配置细化**
  - NUCLEUS_SEGMENTATION_SETTINGS 新增更多可调参数
  - 支持 background_percentile、gaussian_blur_size 等预处理
  - 支持 min_mask_object_area、min_mask_hole_area 形态学过滤
  - 支持 opening_radius、closing_radius 形态学操作

- **阳性阈值配置通道化**
  - CHANNEL_SETTINGS 支持按通道独立配置
  - 每个通道可独立设置 background_percentile、positive_min_area、positive_threshold、color

- **命令行接口扩展**
  - 新增 --recursive-scan 递归扫描参数
  - 新增 --filename-role-map、--filename-channel-map 文件名映射参数
  - 新增 --raw-display-mode、--positive-display-mode 显示模式参数

## [1.0.4] - 2025-03-17

### 改进

- **通道识别增强**
  - 优先使用 TIFF/OME 元数据中的通道名
  - 支持 RGB TIFF (R/G/B -> 594/488/DAPI) 和 BGR 彩图的颜色语义自动映射
  - 新增 DAPI 启发式自动识别（基于形态学特征评分）
  - 添加 STRICT_CHANNEL_DETECTION 模式，模糊时直接报错避免静默错配

- **QC 输出增强**
  - 新增 channel_diagnostics 诊断图，显示每个原始 index 的 raw 图和自动映射结果
  - 可直接审查原始通道分配是否正确

- **阈值配置更灵活**
  - 支持三种阈值方法：Otsu（自动）、manual（手动）、percentile（分位数）
  - 可配置 scale、offset、min_value 等参数

- **串色校正改进**
  - 支持自动/手动两种模式
  - 自动模式可配置估计 mask、ratio_percentile 等参数
  - 新增串色诊断图保存

- **命令行接口**
  - 新增 --file-channel-order、--bleedthrough-488-manual 等参数
  - 支持 --all-coloc-pairs 自动生成所有共定位对
