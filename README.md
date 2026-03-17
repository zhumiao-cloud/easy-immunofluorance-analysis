# Easy Immunofluorescence Analysis

批量免疫荧光(IF)图像分析工具，支持自动通道检测、共定位分析和统计比较。

## 功能特性

- **自动通道识别**: 支持 TIFF/OME 元数据、RGB/BGR 颜色语义、启发式 DAPI 检测
- **细胞核分割**: 基于分水岭算法的自动细胞核识别  
- **信号强度分析**: 背景扣除、阳性面积分数计算
- **共定位分析**: Pearson 相关系数、Manders 系数、Overlap 系数
- **串色校正**: 自动/手动模式去除 DAPI 对 488 通道的串色
- **QC 诊断图**: 通道映射验证、阈值预览、分割结果可视化
- **两组统计比较**: Welch t-test 或 Mann-Whitney U，自动正态性检验，FDR 校正

## 使用方法

编辑脚本顶部的 用户配置区，设置输入输出路径和分析参数，然后运行：

```bash
python if_analysis.py
```

## 输出文件

- per_image_results.csv - 每张图像的详细分析结果
- group_statistics.csv - 两组间统计比较  
- summary_plots.pdf/png - 汇总图表
- analysis_report.txt - 文本报告
- qc_overlays/ - QC 叠加图
- channel_diagnostics/ - 通道诊断图

## 依赖

```bash
pip install numpy pandas scipy scikit-image opencv-python matplotlib tifffile pillow
```

## 更新日志

查看 [CHANGELOG.md](CHANGELOG.md) 了解版本更新详情。

## 许可证

MIT License
