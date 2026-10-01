# =============================================================================
# 文件名：build_AF_3GEO_SeuratObj.R
#
# 目的：
#   仅整合 GSE255612、GSE198204、GSE238242 中的人类 RNA/snRNA 数据。
#   每个数据集先独立完成 metadata 整理、基础 QC、DoubletFinder、DecontX，
#   再把所有保留下来的细胞合并到一个 Seurat v5 对象并用 Harmony 去批次整合。
#
# 重要原则：
#   1. 不保留 GSE198204 的小鼠 Sham/TAC 数据。
#   2. 不把 GSE238242 的 ATAC 模态放入本脚本；这里只处理 RNA/GEX。
#   3. 每个样本独立做 QC、DoubletFinder、DecontX，避免把不同 donor/library 混在一起估算。
#   4. DecontX 以 contamination <= 0.30 作为保留阈值，即污染比例超过 30% 的细胞删除。
#   5. DoubletFinder 的 expected doublet rate 根据 10x recovered-cell/multiplet-rate 表插值估算，
#      并进一步根据样本内初步聚类做 homotypic doublet 调整。
#   6. 最终 Harmony 以每个样本/library 为 batch 单位；AF/Control 不作为 batch 变量。
#   7. 不自动安装包，不自动改输入文件，不吞掉报错。任何关键结构不符合预期直接停止。
#
# GSE255612 特别说明：
#   GEO 当前 supplementary 提供的是 processed expression matrix；这不是 DecontX 所要求的
#   原始整数 counts，也不应当伪装成 DoubletFinder 的原始输入。
#   因此本脚本要求先从该数据集 SRA 原始 reads 用统一 Cell Ranger 流程重建 34 个 donor
#   的 filtered_feature_bc_matrix，并放在：
#
#   GSE255612/raw_counts_by_sample/<sample_id>/
#       matrix.mtx.gz
#       barcodes.tsv.gz
#       features.tsv.gz
#
#   本脚本不会下载 FASTQ，也不会运行 Cell Ranger；这里只消费已经重建好的原始 counts。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 加载依赖包并固定 Seurat v5 数据结构
# -----------------------------------------------------------------------------
# Seurat/SeuratObject：构建、预处理、合并和整合对象。
# Matrix：稀疏矩阵。
# DoubletFinder：双细胞识别。
# SingleCellExperiment + decontX：环境 RNA 污染估计和校正。
# harmony：Seurat v5 HarmonyIntegration 的底层依赖。
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(DoubletFinder)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(decontX)
  library(harmony)
})

# 本脚本按 Seurat v5 / Assay5 编写。
stopifnot(packageVersion("Seurat") >= "5.0.0")
stopifnot(packageVersion("SeuratObject") >= "5.0.0")
options(Seurat.object.assay.version = "v5")


# -----------------------------------------------------------------------------
# 1. 设置输入目录、输出目录和可编辑参数
# -----------------------------------------------------------------------------
# 默认从 AF原始数据 目录启动 R；如脚本放在别处，只修改 base_dir。
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

gse255_dir <- file.path(base_dir, "GSE255612")
gse198_dir <- file.path(base_dir, "GSE198204")
gse238_dir <- file.path(base_dir, "GSE238242")

gse255_raw_root <- file.path(gse255_dir, "raw_counts_by_sample")

# 输出目录只保存 R 对象和最终 metadata，不修改三个原始 GEO 文件夹。
output_dir <- file.path(base_dir, "SeuratObj_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 临时解压目录使用 R 会话临时目录，脚本末尾统一删除。
work_dir <- file.path(tempdir(), "AF_3GEO_human_pipeline")
unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# 基础 QC 参数集中放在这里，便于你按项目需要直接修改。
# GSE255612、GSE238242 为 snRNA / Multiome GEX，因此线粒体阈值设得更严格。
# GSE198204 为 scRNA，因此 percent.mt 默认放宽到 20%。
qc_parameters <- data.frame(
  dataset = c("GSE255612", "GSE198204", "GSE238242"),
  min_features = c(200, 200, 200),
  max_features = c(8000, 8000, 8000),
  min_counts = c(500, 500, 500),
  max_percent_mt = c(10, 20, 10),
  stringsAsFactors = FALSE
)

# DoubletFinder 固定使用 30 个 PC；pK 由每个样本自己的 sweep 最大 BCmetric 自动确定。
df_pcs <- 1:30

df_pN <- 0.25

# DecontX 30% 阈值：污染比例大于 0.30 的细胞删除。
decontx_max_contamination <- 0.30

# 每个样本初步聚类仅用于 DoubletFinder homotypic adjustment 和 DecontX 的 z 标签。
precluster_resolution <- 0.5

# 最终跨样本整合使用 3000 个高变基因和 30 个 PC。
integration_nfeatures <- 3000
integration_pcs <- 1:30


# -----------------------------------------------------------------------------
# 2. 定义 10x multiplet rate 估算函数
# -----------------------------------------------------------------------------
# 10x Chromium 3' v3.1 官方表：
# recovered cells 500/1000/.../10000 对应约 0.4%/0.8%/.../8.0% multiplet rate。
# 对实际保留下来的细胞数做线性插值；超过表格范围时固定使用边界值。
tenx_doublet_rate <- function(n_cells) {
  recovered_cells <- c(
    500, 1000, 2000, 3000, 4000,
    5000, 6000, 7000, 8000, 9000, 10000
  )
  multiplet_rate <- c(
    0.004, 0.008, 0.016, 0.024, 0.032,
    0.040, 0.048, 0.056, 0.064, 0.072, 0.080
  )
  clipped_cells <- pmin(pmax(n_cells, 500), 10000)
  as.numeric(
    approx(
      x = recovered_cells,
      y = multiplet_rate,
      xout = clipped_cells,
      rule = 2
    )$y
  )
}


# -----------------------------------------------------------------------------
# 3. 定义统一 metadata 字段写入函数
# -----------------------------------------------------------------------------
# 每个 Seurat 对象都统一具有以下标准字段：
# dataset / GSE / GSM / BioSample / sample_id / donor_id / sample_uid / batch_id
# species / tissue / condition_original / AF_status / sex / age_years / cell_type
# platform / library_strategy / technology / genome_build / source_matrix_state
#
# 对于数据库本身提供的更多 cell-level metadata，后面分别以 gse255_ / gse238_ 前缀完整保留。
add_standard_metadata <- function(
  obj,
  dataset,
  GSM,
  BioSample,
  sample_id,
  donor_id,
  tissue,
  condition_original,
  AF_status,
  sex,
  age_years,
  cell_type,
  platform,
  library_strategy,
  technology,
  genome_build,
  source_matrix_state
) {
  obj$dataset <- dataset
  obj$GSE <- dataset
  obj$GSM <- GSM
  obj$BioSample <- BioSample
  obj$sample_id <- sample_id
  obj$donor_id <- donor_id
  obj$sample_uid <- paste(dataset, sample_id, sep = "__")
  obj$batch_id <- paste(dataset, sample_id, sep = "__")
  obj$species <- "Homo sapiens"
  obj$tissue <- tissue
  obj$condition_original <- condition_original
  obj$AF_status <- AF_status
  obj$sex <- sex
  obj$age_years <- age_years
  obj$cell_type <- cell_type
  obj$platform <- platform
  obj$library_strategy <- library_strategy
  obj$technology <- technology
  obj$genome_build <- genome_build
  obj$source_matrix_state <- source_matrix_state
  obj
}


# -----------------------------------------------------------------------------
# 4. 定义基础 QC + DoubletFinder + DecontX 的统一单样本函数
# -----------------------------------------------------------------------------
# 输入必须是单个 donor/library 的原始整数 counts Seurat 对象。
# 流程顺序固定：
#   原始 QC 指标 -> 基础 QC -> 初步降维聚类 -> DoubletFinder -> 删除 doublet
#   -> DecontX -> 删除 contamination >30% -> 用 DecontX 校正后的整数 counts 重建对象。
run_qc_doubletfinder_decontx <- function(obj, dataset_name) {
  # 读取该数据库的固定 QC 参数。
  qc_row <- qc_parameters[
    match(dataset_name, qc_parameters$dataset),
    ,
    drop = FALSE
  ]
  stopifnot(nrow(qc_row) == 1)

  # 验证当前输入确实存在 counts layer，并且非零值都是整数。
  raw_counts <- SeuratObject::LayerData(
    object = obj,
    assay = "RNA",
    layer = "counts"
  )
  stopifnot(ncol(raw_counts) == ncol(obj))
  stopifnot(all(raw_counts@x >= 0))
  stopifnot(all(raw_counts@x == floor(raw_counts@x)))

  # 计算常规 QC 指标；这里只记录，不做任何隐式修正。
  obj$percent.mt <- Seurat::PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )
  obj$percent.ribo <- Seurat::PercentageFeatureSet(
    obj,
    pattern = "^RP[SL]"
  )

  # 保存过滤前的基础 QC 指标，便于以后回溯。
  obj$raw_nCount_RNA <- obj$nCount_RNA
  obj$raw_nFeature_RNA <- obj$nFeature_RNA
  obj$raw_percent_mt <- obj$percent.mt
  obj$raw_percent_ribo <- obj$percent.ribo

  # 按顶部 qc_parameters 的显式阈值生成 qc_pass。
  obj$qc_pass <-
    obj$nFeature_RNA >= qc_row$min_features &
    obj$nFeature_RNA <= qc_row$max_features &
    obj$nCount_RNA >= qc_row$min_counts &
    obj$percent.mt <= qc_row$max_percent_mt

  # 只保留基础 QC 通过的细胞。
  obj <- subset(
    x = obj,
    subset = qc_pass
  )
  stopifnot(ncol(obj) > 100)

  # DoubletFinder 需要先完成常规标准化、高变基因、PCA 和粗聚类。
  obj <- Seurat::NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  obj <- Seurat::FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  obj <- Seurat::ScaleData(
    obj,
    features = Seurat::VariableFeatures(obj),
    verbose = FALSE
  )
  obj <- Seurat::RunPCA(
    obj,
    features = Seurat::VariableFeatures(obj),
    npcs = max(df_pcs),
    verbose = FALSE
  )
  obj <- Seurat::FindNeighbors(
    obj,
    dims = df_pcs,
    verbose = FALSE
  )
  obj <- Seurat::FindClusters(
    obj,
    resolution = precluster_resolution,
    verbose = FALSE
  )

  # 对该样本独立做 pN-pK 参数扫描。
  sweep_list <- DoubletFinder::paramSweep(
    seu = obj,
    PCs = df_pcs,
    sct = FALSE,
    num.cores = 1
  )
  sweep_stats <- DoubletFinder::summarizeSweep(
    sweep_list,
    GT = FALSE
  )
  bcmvn <- DoubletFinder::find.pK(sweep_stats)

  # 取 BCmetric 最大值对应的 pK，完全由该样本自身数据确定。
  best_pK <- as.numeric(
    as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]
    )
  )
  stopifnot(length(best_pK) == 1)
  stopifnot(!is.na(best_pK))

  # 根据 10x multiplet 表估计 expected doublets。
  expected_rate <- tenx_doublet_rate(ncol(obj))
  expected_doublets <- round(expected_rate * ncol(obj))

  # 使用样本内粗聚类估计 homotypic doublet 比例。
  homotypic_prop <- DoubletFinder::modelHomotypic(
    annotations = obj$seurat_clusters
  )
  expected_doublets_adjusted <- max(
    1,
    round(
      expected_doublets * (1 - homotypic_prop)
    )
  )

  # 正式运行 DoubletFinder。
  obj <- DoubletFinder::doubletFinder(
    seu = obj,
    PCs = df_pcs,
    pN = df_pN,
    pK = best_pK,
    nExp = expected_doublets_adjusted,
    reuse.pANN = FALSE,
    sct = FALSE,
    annotations = obj$seurat_clusters
  )

  # DoubletFinder 会生成一个 DF.classifications_* 列；严格要求只有一个。
  df_class_column <- grep(
    pattern = "^DF.classifications_",
    x = colnames(obj@meta.data),
    value = TRUE
  )
  stopifnot(length(df_class_column) == 1)

  # 保存标准化 DoubletFinder 输出字段。
  obj$doublet_status <- obj@meta.data[[df_class_column]]
  obj$doublet_expected_rate <- expected_rate
  obj$doublet_expected_n <- expected_doublets
  obj$doublet_expected_n_homotypic_adjusted <- expected_doublets_adjusted
  obj$doublet_pK <- best_pK

  # 只保留 DoubletFinder 判定的 Singlet。
  obj <- subset(
    x = obj,
    subset = doublet_status == "Singlet"
  )
  stopifnot(ncol(obj) > 100)

  # 提取当前 singlet 的原始整数 counts，转换为 SingleCellExperiment。
  singlet_counts <- SeuratObject::LayerData(
    object = obj,
    assay = "RNA",
    layer = "counts"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = singlet_counts)
  )

  # 直接使用 DoubletFinder 前的样本内粗聚类作为 DecontX 的 z 标签。
  # DecontX 文档允许使用用户提供的 cluster labels，避免再次自动聚类。
  sce <- decontX::decontX(
    x = sce,
    assayName = "counts",
    z = as.character(obj$seurat_clusters),
    seed = 12345,
    verbose = TRUE
  )

  # 读取每个细胞的污染比例并写回 Seurat metadata。
  obj$decontX_contamination <- as.numeric(
    SummarizedExperiment::colData(sce)$decontX_contamination
  )
  obj$decontX_pass_30pct <-
    obj$decontX_contamination <= decontx_max_contamination

  # 只保留污染比例不超过 30% 的细胞。
  obj <- subset(
    x = obj,
    subset = decontX_pass_30pct
  )
  stopifnot(ncol(obj) > 100)

  # DecontX 输出可以是非整数；按照 decontX 官方说明四舍五入成整数 counts，
  # 作为后续 Seurat counts 层使用。
  decont_counts <- round(
    decontX::decontXcounts(sce)
  )
  decont_counts <- decont_counts[
    ,
    Seurat::Cells(obj),
    drop = FALSE
  ]
  decont_counts <- as(decont_counts, "dgCMatrix")
  stopifnot(all(decont_counts@x >= 0))
  stopifnot(all(decont_counts@x == floor(decont_counts@x)))

  # 保存原始 QC 指标到 metadata，然后去掉旧的 nCount/nFeature 字段，
  # 使新对象的 nCount_RNA / nFeature_RNA 由 DecontX 校正后的 counts 自动重新计算。
  final_meta <- obj@meta.data
  final_meta$preDecontX_nCount_RNA <- final_meta$nCount_RNA
  final_meta$preDecontX_nFeature_RNA <- final_meta$nFeature_RNA
  final_meta$nCount_RNA <- NULL
  final_meta$nFeature_RNA <- NULL

  # 用校正后的 counts 重新构建干净对象，并完整继承 metadata。
  final_obj <- SeuratObject::CreateSeuratObject(
    counts = decont_counts,
    assay = "RNA",
    project = unique(final_meta$sample_uid),
    meta.data = final_meta,
    min.cells = 0,
    min.features = 0
  )

  # 重新计算 DecontX 后的线粒体和核糖体比例，便于后续核查。
  final_obj$postDecontX_percent_mt <- Seurat::PercentageFeatureSet(
    final_obj,
    pattern = "^MT-"
  )
  final_obj$postDecontX_percent_ribo <- Seurat::PercentageFeatureSet(
    final_obj,
    pattern = "^RP[SL]"
  )

  # 给每个最终单样本对象生成 LogNormalize data layer，便于单独检查。
  final_obj <- Seurat::NormalizeData(
    final_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  final_obj
}


# -----------------------------------------------------------------------------
# 5. GSE255612：读取官方 cell-level metadata，并准备 34 个人类 donor 表
# -----------------------------------------------------------------------------
# 当前 GEO supplementary 的 metadata 文件继续用于：
# donor / AF 分组 / sex / 作者 cell type 等注释。
# RNA counts 本身必须来自 SRA 重建后的 raw_counts_by_sample。
gse255_meta_file <- file.path(
  gse255_dir,
  "GSE255612_AF_snRNA_MetaData.txt.gz"
)
stopifnot(file.exists(gse255_meta_file))
stopifnot(dir.exists(gse255_raw_root))

gse255_native_meta <- read.delim(
  gzfile(gse255_meta_file),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# GEO metadata 第一条数据行为字段 TYPE 描述，不是真实细胞，固定删除。
gse255_native_meta <- gse255_native_meta[-1, , drop = FALSE]

# 严格要求存在后续匹配所需的核心列。
stopifnot(
  all(
    c(
      "NAME",
      "biosample_id",
      "donor_id",
      "sex",
      "af",
      "cell_type_assigned"
    ) %in% colnames(gse255_native_meta)
  )
)

# 34 个 GEO Sample 与 donor 的固定一一对应关系。
# gem_group 用于把每个 donor 原始 10x barcode 的 -1 后缀映射到聚合 metadata 的 GEM group。
gse255_samples <- data.frame(
  sample_id = c(
    "1279_3n", "1296_1n", "1334_3n", "1357_2n", "1360_1n", "1365_1n",
    "1369_3n", "1377_3n", "1440_1n", "1465_1n", "1488_2n", "1490_3n",
    "1498_1n", "1513_1n", "1515_1n", "1540_3n", "1543_1n", "1558_1n",
    "1561_1n", "1570_1n", "1582_1n", "1588_1n", "1591_1n", "1593_1n",
    "1603_1n", "1610_1n", "1623_1n", "1626_1n", "1696_1n", "1698_3n",
    "1741_1n", "1758_1n", "1762_1n", "1789_1n"
  ),
  GSM = paste0("GSM", 8076040:8076073),
  gem_group = seq_len(34),
  stringsAsFactors = FALSE
)
gse255_samples$donor_id <- paste0("P", gse255_samples$sample_id)

# case/control 到统一 AF_status 的固定映射。
gse255_af_map <- c(case = "AF", control = "Control")

# 逐 donor 读取原始整数 counts，并插入该 donor 的官方 cell-level metadata。
gse255_obj_list <- lapply(
  seq_len(nrow(gse255_samples)),
  function(index) {
    sample_row <- gse255_samples[index, , drop = FALSE]
    sample_dir <- file.path(
      gse255_raw_root,
      sample_row$sample_id
    )
    stopifnot(dir.exists(sample_dir))

    counts <- Seurat::Read10X(
      data.dir = sample_dir,
      gene.column = 2,
      unique.features = TRUE,
      strip.suffix = FALSE
    )
    counts <- as(counts, "dgCMatrix")

    # 当前 donor 的官方 metadata。
    native_sample_meta <- gse255_native_meta[
      gse255_native_meta$biosample_id == sample_row$donor_id,
      ,
      drop = FALSE
    ]
    stopifnot(nrow(native_sample_meta) > 0)

    # Cell Ranger aggr/GEM group 规则：单样本原始 barcode 末尾 -1 替换成该 donor 的 GEM group。
    aggregated_cell_id <- sub(
      pattern = "-1$",
      replacement = paste0("-", sample_row$gem_group),
      x = colnames(counts)
    )
    native_index <- match(
      aggregated_cell_id,
      native_sample_meta$NAME
    )
    stopifnot(!anyNA(native_index))
    native_sample_meta <- native_sample_meta[
      native_index,
      ,
      drop = FALSE
    ]
    stopifnot(
      identical(
        native_sample_meta$NAME,
        aggregated_cell_id
      )
    )

    # 给细胞名加入 dataset + sample 前缀，保证全局唯一。
    new_cell_names <- paste0(
      "GSE255612__",
      sample_row$sample_id,
      "__",
      colnames(counts)
    )
    colnames(counts) <- new_cell_names

    obj <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_row$sample_id,
      min.cells = 0,
      min.features = 0
    )

    # 把 GSE255612 原始 cell-level metadata 全部保留，并统一加 gse255_ 前缀防止字段冲突。
    native_meta_for_seurat <- native_sample_meta
    colnames(native_meta_for_seurat) <- paste0(
      "gse255_",
      make.names(colnames(native_meta_for_seurat))
    )
    rownames(native_meta_for_seurat) <- new_cell_names
    obj <- Seurat::AddMetaData(
      object = obj,
      metadata = native_meta_for_seurat
    )

    # 从数据库原始 metadata 中提取标准字段。
    condition_original <- unname(
      gse255_af_map[native_sample_meta$af]
    )
    stopifnot(!anyNA(condition_original))

    obj <- add_standard_metadata(
      obj = obj,
      dataset = "GSE255612",
      GSM = sample_row$GSM,
      BioSample = NA_character_,
      sample_id = sample_row$sample_id,
      donor_id = sample_row$donor_id,
      tissue = "Left atrium",
      condition_original = condition_original,
      AF_status = condition_original,
      sex = as.character(native_sample_meta$sex),
      age_years = NA_real_,
      cell_type = as.character(native_sample_meta$cell_type_assigned),
      platform = "Illumina NovaSeq 6000",
      library_strategy = "RNA-Seq",
      technology = "10x Genomics single-nucleus RNA-seq",
      genome_build = "GRCh38",
      source_matrix_state = "raw_integer_counts_rebuilt_from_SRA"
    )

    # 逐 donor 独立做基础 QC、DoubletFinder、DecontX。
    run_qc_doubletfinder_decontx(
      obj = obj,
      dataset_name = "GSE255612"
    )
  }
)
names(gse255_obj_list) <- gse255_samples$sample_id


# -----------------------------------------------------------------------------
# 6. GSE198204：只读取 6 个人类左心耳样本
# -----------------------------------------------------------------------------
# GEO 当前 human 样本：3 个 SR + 3 个 AF。
# mouse Sham/TAC 不进入本脚本。
gse198_tar <- file.path(
  gse198_dir,
  "GSE198204_RAW.tar"
)
stopifnot(file.exists(gse198_tar))

gse198_extract_dir <- file.path(
  work_dir,
  "GSE198204_extracted"
)
dir.create(gse198_extract_dir, recursive = TRUE, showWarnings = FALSE)
untar(gse198_tar, exdir = gse198_extract_dir)

# 下表只写公开 GEO 记录已经确认的信息；没有公开确认的信息明确写 NA，不猜测。
gse198_samples <- data.frame(
  file_prefix = c(
    "GSM5940688_LA_SR_01",
    "GSM5940689_LA_SR_02",
    "GSM6552875_LA_SR_03",
    "GSM5940691_LA_AF_01",
    "GSM6552876_LA_AF_02",
    "GSM6552877_LA_AF_03"
  ),
  GSM = c(
    "GSM5940688",
    "GSM5940689",
    "GSM6552875",
    "GSM5940691",
    "GSM6552876",
    "GSM6552877"
  ),
  BioSample = c(
    "SAMN26540411",
    NA,
    "SAMN30711027",
    "SAMN26540408",
    "SAMN30711026",
    "SAMN30711025"
  ),
  sample_id = c(
    "LA_SR_01",
    "LA_SR_02",
    "LA_SR_03",
    "LA_AF_01",
    "LA_AF_02",
    "LA_AF_03"
  ),
  condition_original = c(
    "SR", "SR", "SR",
    "AF", "AF", "AF"
  ),
  AF_status = c(
    "Control", "Control", "Control",
    "AF", "AF", "AF"
  ),
  age_years = c(
    55,
    NA,
    42,
    55,
    40,
    62
  ),
  stringsAsFactors = FALSE
)

# 逐样本读取 GEO 提供的标准 10x 原始 gene-count matrix。
gse198_obj_list <- lapply(
  seq_len(nrow(gse198_samples)),
  function(index) {
    sample_row <- gse198_samples[index, , drop = FALSE]

    counts <- Seurat::ReadMtx(
      mtx = file.path(
        gse198_extract_dir,
        paste0(sample_row$file_prefix, "_matrix.mtx.gz")
      ),
      cells = file.path(
        gse198_extract_dir,
        paste0(sample_row$file_prefix, "_barcodes.tsv.gz")
      ),
      features = file.path(
        gse198_extract_dir,
        paste0(sample_row$file_prefix, "_features.tsv.gz")
      ),
      cell.column = 1,
      feature.column = 2,
      unique.features = TRUE,
      strip.suffix = FALSE
    )
    counts <- as(counts, "dgCMatrix")

    new_cell_names <- paste0(
      "GSE198204__",
      sample_row$sample_id,
      "__",
      colnames(counts)
    )
    colnames(counts) <- new_cell_names

    obj <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_row$sample_id,
      min.cells = 0,
      min.features = 0
    )

    # GSE198204 没有额外的 cell-level metadata 文件，
    # 因此这里把 GEO sample record 披露的样本级信息标准化后复制到该样本全部细胞。
    obj <- add_standard_metadata(
      obj = obj,
      dataset = "GSE198204",
      GSM = sample_row$GSM,
      BioSample = sample_row$BioSample,
      sample_id = sample_row$sample_id,
      donor_id = sample_row$sample_id,
      tissue = "Left atrial appendage",
      condition_original = sample_row$condition_original,
      AF_status = sample_row$AF_status,
      sex = NA_character_,
      age_years = sample_row$age_years,
      cell_type = NA_character_,
      platform = "Illumina NovaSeq 6000",
      library_strategy = "RNA-Seq",
      technology = "10x Genomics Single Cell 3' v3 scRNA-seq",
      genome_build = "GRCh38",
      source_matrix_state = "GEO_raw_gene_counts_matrix"
    )

    run_qc_doubletfinder_decontx(
      obj = obj,
      dataset_name = "GSE198204"
    )
  }
)
names(gse198_obj_list) <- gse198_samples$sample_id


# -----------------------------------------------------------------------------
# 7. GSE238242：只读取 7 个人类 RNA/GEX 样本，并完整插入作者 metadata
# -----------------------------------------------------------------------------
# GSE238242 是 human left atrial appendage 10x Multiome。
# 本脚本只处理 RNA：4 个 SR + 3 个 AF；ATAC 文件完全不动。
gse238_tar <- file.path(
  gse238_dir,
  "GSE238242_RAW.tar"
)
gse238_meta_file <- file.path(
  gse238_dir,
  "GSE238242_snAF.metadata.tsv.gz"
)
stopifnot(file.exists(gse238_tar))
stopifnot(file.exists(gse238_meta_file))

gse238_extract_dir <- file.path(
  work_dir,
  "GSE238242_extracted"
)
dir.create(gse238_extract_dir, recursive = TRUE, showWarnings = FALSE)
untar(gse238_tar, exdir = gse238_extract_dir)

# 读取作者公开的完整 cell-level metadata。
gse238_native_meta <- read.delim(
  gzfile(gse238_meta_file),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  all(
    c("cell_type", "sample", "sex", "Rhythm") %in%
      colnames(gse238_native_meta)
  )
)

# 7 个 RNA 样本；age/sex 来自研究公开 sample metadata。
gse238_samples <- data.frame(
  sample_id = c(
    "CF69", "CF77", "CF89", "CF91",
    "CF93", "CF97", "CF102"
  ),
  GSM = c(
    "GSM7660990", "GSM7660991", "GSM7660992", "GSM7660993",
    "GSM7660994", "GSM7660995", "GSM7660996"
  ),
  condition_original = c(
    "SR", "SR", "SR", "SR",
    "AF", "AF", "AF"
  ),
  AF_status = c(
    "Control", "Control", "Control", "Control",
    "AF", "AF", "AF"
  ),
  age_years = c(
    51, 78, 69, 57,
    80, 58, 75
  ),
  sex_expected = c(
    "F", "M", "F", "M",
    "M", "M", "M"
  ),
  counts_file = c(
    "GSM7660990_CF69RNA.counts.tsv.gz",
    "GSM7660991_CF77RNA.counts.tsv.gz",
    "GSM7660992_CF89RNA.counts.tsv.gz",
    "GSM7660993_CF91RNA.counts.tsv.gz",
    "GSM7660994_CF93RNA.counts.tsv.gz",
    "GSM7660995_CF97RNA.counts.tsv.gz",
    "GSM7660996_CF102RNA.counts.tsv.gz"
  ),
  stringsAsFactors = FALSE
)

# 逐 donor 读取 gene × cell TSV counts，并精确匹配作者 metadata。
gse238_obj_list <- lapply(
  seq_len(nrow(gse238_samples)),
  function(index) {
    sample_row <- gse238_samples[index, , drop = FALSE]

    dense_counts <- read.delim(
      gzfile(
        file.path(
          gse238_extract_dir,
          sample_row$counts_file
        )
      ),
      header = TRUE,
      row.names = 1,
      sep = "\t",
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    counts <- Matrix::Matrix(
      as.matrix(dense_counts),
      sparse = TRUE
    )
    counts <- as(counts, "dgCMatrix")
    rownames(counts) <- make.unique(rownames(counts))
    rm(dense_counts)
    invisible(gc())

    # 只取该 donor 的作者 metadata。
    native_sample_meta <- gse238_native_meta[
      gse238_native_meta$sample == sample_row$sample_id,
      ,
      drop = FALSE
    ]
    stopifnot(nrow(native_sample_meta) > 0)

    # counts 和作者 metadata 的末尾批次编号不同，先去掉最后一个 -数字 后缀再一一匹配。
    native_sample_meta$barcode_core <- sub(
      "-[0-9]+$",
      "",
      rownames(native_sample_meta)
    )
    count_barcode_core <- sub(
      "-[0-9]+$",
      "",
      colnames(counts)
    )
    native_index <- match(
      count_barcode_core,
      native_sample_meta$barcode_core
    )
    stopifnot(!anyNA(native_index))
    native_sample_meta <- native_sample_meta[
      native_index,
      ,
      drop = FALSE
    ]
    stopifnot(
      all(
        native_sample_meta$Rhythm == sample_row$condition_original
      )
    )

    # 全局唯一细胞名。
    new_cell_names <- paste0(
      "GSE238242__",
      sample_row$sample_id,
      "__",
      colnames(counts)
    )
    colnames(counts) <- new_cell_names

    obj <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_row$sample_id,
      min.cells = 0,
      min.features = 0
    )

    # 完整保留 GSE238242 原始 cell-level metadata，统一加 gse238_ 前缀。
    native_meta_for_seurat <- native_sample_meta
    native_meta_for_seurat$barcode_core <- NULL
    colnames(native_meta_for_seurat) <- paste0(
      "gse238_",
      make.names(colnames(native_meta_for_seurat))
    )
    rownames(native_meta_for_seurat) <- new_cell_names
    obj <- Seurat::AddMetaData(
      object = obj,
      metadata = native_meta_for_seurat
    )

    # 写入跨数据库统一标准字段。
    obj <- add_standard_metadata(
      obj = obj,
      dataset = "GSE238242",
      GSM = sample_row$GSM,
      BioSample = NA_character_,
      sample_id = sample_row$sample_id,
      donor_id = sample_row$sample_id,
      tissue = "Left atrial appendage",
      condition_original = as.character(native_sample_meta$Rhythm),
      AF_status = sample_row$AF_status,
      sex = as.character(native_sample_meta$sex),
      age_years = sample_row$age_years,
      cell_type = as.character(native_sample_meta$cell_type),
      platform = "Illumina NovaSeq 6000",
      library_strategy = "RNA-Seq",
      technology = "10x Genomics Single Cell Multiome Gene Expression",
      genome_build = "GRCh38",
      source_matrix_state = "GEO_raw_RNA_counts"
    )

    # 作者 metadata 中 sex 必须与公开 sample-level sex 一致。
    standardized_sex <- toupper(substr(obj$sex, 1, 1))
    stopifnot(
      all(
        standardized_sex == sample_row$sex_expected
      )
    )

    run_qc_doubletfinder_decontx(
      obj = obj,
      dataset_name = "GSE238242"
    )
  }
)
names(gse238_obj_list) <- gse238_samples$sample_id


# -----------------------------------------------------------------------------
# 8. 分别生成三个数据库清洗后的 human RNA Seurat 对象
# -----------------------------------------------------------------------------
# merge.data = FALSE：只合并 counts，后续统一重新 NormalizeData。
# Seurat v5 会保留每个样本独立 layer，后续 HarmonyIntegration 正是利用这些 layer 做整合。
gse255_obj <- merge(
  x = gse255_obj_list[[1]],
  y = gse255_obj_list[-1],
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "GSE255612_human_clean"
)

gse198_obj <- merge(
  x = gse198_obj_list[[1]],
  y = gse198_obj_list[-1],
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "GSE198204_human_clean"
)

gse238_obj <- merge(
  x = gse238_obj_list[[1]],
  y = gse238_obj_list[-1],
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "GSE238242_human_clean"
)

# 保存三个数据库各自完成 QC + DoubletFinder + DecontX 后的对象。
saveRDS(
  gse255_obj,
  file.path(output_dir, "GSE255612_human_QC_DoubletFinder_DecontX.rds"),
  compress = FALSE
)
saveRDS(
  gse198_obj,
  file.path(output_dir, "GSE198204_human_QC_DoubletFinder_DecontX.rds"),
  compress = FALSE
)
saveRDS(
  gse238_obj,
  file.path(output_dir, "GSE238242_human_QC_DoubletFinder_DecontX.rds"),
  compress = FALSE
)


# -----------------------------------------------------------------------------
# 9. 合并三个数据库，并用 Seurat v5 HarmonyIntegration 去批次
# -----------------------------------------------------------------------------
# 三个对象都已经是 DecontX 校正后的整数 counts。
# 最终 merge 继续保留样本级 layer，不把样本来源提前压扁。
seurobj <- merge(
  x = gse255_obj,
  y = list(
    gse198_obj,
    gse238_obj
  ),
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "AF_3GEO_human"
)

# 最终对象必须只含人类数据，并且三个 GSE 都存在。
stopifnot(all(seurobj$species == "Homo sapiens"))
stopifnot(
  setequal(
    unique(seurobj$dataset),
    c("GSE255612", "GSE198204", "GSE238242")
  )
)

# 各样本 layer 分别做 LogNormalize。
seurobj <- Seurat::NormalizeData(
  seurobj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# 在多 layer 对象上寻找跨样本共同高变基因。
seurobj <- Seurat::FindVariableFeatures(
  seurobj,
  selection.method = "vst",
  nfeatures = integration_nfeatures,
  verbose = FALSE
)

# 使用共同高变基因 ScaleData 和 PCA。
seurobj <- Seurat::ScaleData(
  seurobj,
  features = Seurat::VariableFeatures(seurobj),
  verbose = FALSE
)
seurobj <- Seurat::RunPCA(
  seurobj,
  features = Seurat::VariableFeatures(seurobj),
  npcs = max(integration_pcs),
  verbose = FALSE
)

# Seurat v5 HarmonyIntegration：每个样本/library 的 layer 作为 batch 单位。
# AF_status 不作为 batch 字段，因此不会直接把疾病分组当作技术因素消除。
seurobj <- Seurat::IntegrateLayers(
  object = seurobj,
  method = Seurat::HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

# 在 Harmony 空间构建邻居、聚类和 UMAP。
seurobj <- Seurat::FindNeighbors(
  seurobj,
  reduction = "harmony",
  dims = integration_pcs,
  verbose = FALSE
)
seurobj <- Seurat::FindClusters(
  seurobj,
  resolution = 0.5,
  cluster.name = "harmony_clusters",
  verbose = FALSE
)
seurobj <- Seurat::RunUMAP(
  seurobj,
  reduction = "harmony",
  dims = integration_pcs,
  reduction.name = "umap.harmony",
  verbose = FALSE
)

# 完成整合后重新 JoinLayers，方便后续 FindMarkers/差异表达等标准 Seurat 操作。
seurobj[["RNA"]] <- SeuratObject::JoinLayers(
  seurobj[["RNA"]]
)

# 明确记录最终对象处理状态。
seurobj$integration_method <- "Seurat_v5_HarmonyIntegration"
seurobj$decontX_threshold <- decontx_max_contamination
seurobj$pipeline_state <- "QC_DoubletFinder_DecontX30_Harmony_integrated"


# -----------------------------------------------------------------------------
# 10. 导出最终 Seurat 对象和标准化 metadata
# -----------------------------------------------------------------------------
# 最终主对象。
saveRDS(
  seurobj,
  file.path(
    output_dir,
    "AF_3GEO_human_QC_DoubletFinder_DecontX30_Harmony_seurobj.rds"
  ),
  compress = FALSE
)

# 导出完整 metadata：
# 一方面包含跨数据库统一标准字段，另一方面保留 GSE255612 / GSE238242 的原始作者字段。
metadata_output_file <- file.path(
  output_dir,
  "AF_3GEO_human_standardized_metadata.tsv.gz"
)
metadata_connection <- gzfile(
  metadata_output_file,
  open = "wt"
)
write.table(
  seurobj@meta.data,
  file = metadata_connection,
  sep = "\t",
  quote = FALSE,
  col.names = NA,
  row.names = TRUE
)
close(metadata_connection)


# -----------------------------------------------------------------------------
# 11. 运行结束后的结构核查输出
# -----------------------------------------------------------------------------
# 这些语句只打印统计，不再修改数据。
print(seurobj)
print(table(seurobj$dataset))
print(table(seurobj$AF_status))
print(table(seurobj$dataset, seurobj$AF_status))
print(table(seurobj$dataset, seurobj$sample_id))
print(table(seurobj$doublet_status))
print(summary(seurobj$decontX_contamination))
print(SeuratObject::Layers(seurobj[["RNA"]]))
print(Seurat::Reductions(seurobj))

cat(
  "\n最终 Seurat 对象：\n",
  file.path(
    output_dir,
    "AF_3GEO_human_QC_DoubletFinder_DecontX30_Harmony_seurobj.rds"
  ),
  "\n\n最终 metadata：\n",
  metadata_output_file,
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# 12. 清理本脚本产生的临时解压目录
# -----------------------------------------------------------------------------
# 只删除 R tempdir 下的临时文件；三个 GEO 原始目录和正式输出完全不动。
unlink(
  work_dir,
  recursive = TRUE,
  force = TRUE
)
