# Changelog

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
