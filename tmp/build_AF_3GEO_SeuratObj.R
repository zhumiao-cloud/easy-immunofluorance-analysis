# =============================================================================
# 文件名：build_AF_3GEO_SeuratObj.R
#
# 目的：仅处理 GSE255612、GSE198204、GSE238242 的人类 RNA/snRNA 数据。
# 流程：分别插入 metadata -> 基础 QC -> DoubletFinder -> DecontX(30%)
#       -> 三库合并 -> Seurat v5 Harmony 去批次整合 -> 输出 seurobj + metadata。
#
# 重要说明：
# 1. 不保留 GSE198204 小鼠数据；不处理 GSE238242 ATAC。
# 2. 每个 donor/library 独立做 QC、DoubletFinder、DecontX。
# 3. DecontX 的“30%”解释为：保留 contamination <= 0.30 的细胞。
# 4. GSE255612 GEO supplementary 是 processed expression，不可伪装成 raw counts。
#    本脚本要求先从 SRA 用统一 Cell Ranger 流程重建 34 个 donor 的整数 counts，目录：
#      GSE255612/raw_counts_by_sample/<sample_id>/
#        matrix.mtx.gz
#        barcodes.tsv.gz
#        features.tsv.gz
# 5. 本脚本不下载 FASTQ、不跑 Cell Ranger、不自动安装包、不吞报错。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 加载依赖
# -----------------------------------------------------------------------------
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

stopifnot(packageVersion("Seurat") >= "5.0.0")
stopifnot(packageVersion("SeuratObject") >= "5.0.0")
options(Seurat.object.assay.version = "v5")


# -----------------------------------------------------------------------------
# 1. 路径与统一参数
# -----------------------------------------------------------------------------
# 默认从 AF原始数据 目录启动 R；如脚本放在其他目录，仅修改 base_dir。
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
gse255_dir <- file.path(base_dir, "GSE255612")
gse198_dir <- file.path(base_dir, "GSE198204")
gse238_dir <- file.path(base_dir, "GSE238242")
gse255_raw_root <- file.path(gse255_dir, "raw_counts_by_sample")

output_dir <- file.path(base_dir, "SeuratObj_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# TAR 只解压到 R 临时目录，脚本结束后删除；原始 GEO 目录完全不动。
work_dir <- file.path(tempdir(), "AF_3GEO_human_pipeline")
unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# 基础 QC 阈值集中管理，便于按你的标准直接修改。
qc_parameters <- data.frame(
  dataset = c("GSE255612", "GSE198204", "GSE238242"),
  min_features = c(200, 200, 200),
  max_features = c(8000, 8000, 8000),
  min_counts = c(500, 500, 500),
  max_percent_mt = c(10, 20, 10),
  stringsAsFactors = FALSE
)

df_pcs <- 1:30
df_pN <- 0.25
precluster_resolution <- 0.5
decontx_max_contamination <- 0.30
integration_nfeatures <- 3000
integration_pcs <- 1:30


# -----------------------------------------------------------------------------
# 2. 10x expected doublet rate
# -----------------------------------------------------------------------------
# 依据 10x Chromium 3' v3.1 recovered cells / multiplet rate 表线性插值。
tenx_doublet_rate <- function(n_cells) {
  recovered_cells <- c(
    500, 1000, 2000, 3000, 4000, 5000,
    6000, 7000, 8000, 9000, 10000
  )
  multiplet_rate <- c(
    0.004, 0.008, 0.016, 0.024, 0.032, 0.040,
    0.048, 0.056, 0.064, 0.072, 0.080
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
# 3. 统一标准 metadata 字段
# -----------------------------------------------------------------------------
# 三个数据库最终统一具有：
# dataset/GSE/GSM/BioSample/sample_id/donor_id/sample_uid/batch_id/species/tissue
# condition_original/AF_status/sex/age_years/cell_type/platform/library_strategy
# technology/genome_build/source_matrix_state。
# 数据库原始 cell-level metadata 另外完整保留，并增加 gse255_ / gse238_ 前缀。
add_standard_metadata <- function(
  obj, dataset, GSM, BioSample, sample_id, donor_id, tissue,
  condition_original, AF_status, sex, age_years, cell_type,
  platform, library_strategy, technology, genome_build, source_matrix_state
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
# 4. 单样本：基础 QC -> DoubletFinder -> DecontX
# -----------------------------------------------------------------------------
run_qc_doubletfinder_decontx <- function(obj, dataset_name) {
  # 读取该数据库预先写死的 QC 参数。
  qc_row <- qc_parameters[
    match(dataset_name, qc_parameters$dataset),
    ,
    drop = FALSE
  ]
  stopifnot(nrow(qc_row) == 1)

  # DecontX 和 DoubletFinder 必须基于原始整数 counts；这里做硬性检查。
  raw_counts <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
  stopifnot(ncol(raw_counts) == ncol(obj))
  stopifnot(all(raw_counts@x >= 0))
  stopifnot(all(raw_counts@x == floor(raw_counts@x)))

  # 基础 QC 指标。
  obj$percent.mt <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent.ribo <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj$raw_nCount_RNA <- obj$nCount_RNA
  obj$raw_nFeature_RNA <- obj$nFeature_RNA
  obj$raw_percent_mt <- obj$percent.mt
  obj$raw_percent_ribo <- obj$percent.ribo

  # 显式 QC 判定。
  obj$qc_pass <-
    obj$nFeature_RNA >= qc_row$min_features &
    obj$nFeature_RNA <= qc_row$max_features &
    obj$nCount_RNA >= qc_row$min_counts &
    obj$percent.mt <= qc_row$max_percent_mt

  obj <- subset(obj, subset = qc_pass)
  stopifnot(ncol(obj) > 100)

  # DoubletFinder 前的标准 Seurat 预处理和粗聚类。
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
  obj <- Seurat::FindNeighbors(obj, dims = df_pcs, verbose = FALSE)
  obj <- Seurat::FindClusters(
    obj,
    resolution = precluster_resolution,
    verbose = FALSE
  )

  # 每个样本独立做 pK sweep，取 BCmetric 最大值对应的 pK。
  sweep_list <- DoubletFinder::paramSweep(
    seu = obj,
    PCs = df_pcs,
    sct = FALSE,
    num.cores = 1
  )
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep_stats)
  best_pK <- as.numeric(
    as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]
    )
  )
  stopifnot(length(best_pK) == 1)
  stopifnot(!is.na(best_pK))

  # 10x expected doublets + homotypic adjustment。
  expected_rate <- tenx_doublet_rate(ncol(obj))
  expected_doublets <- round(expected_rate * ncol(obj))
  homotypic_prop <- DoubletFinder::modelHomotypic(obj$seurat_clusters)
  expected_doublets_adjusted <- max(
    1,
    round(expected_doublets * (1 - homotypic_prop))
  )

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

  df_class_column <- grep(
    "^DF.classifications_",
    colnames(obj@meta.data),
    value = TRUE
  )
  stopifnot(length(df_class_column) == 1)

  obj$doublet_status <- obj@meta.data[[df_class_column]]
  obj$doublet_expected_rate <- expected_rate
  obj$doublet_expected_n <- expected_doublets
  obj$doublet_expected_n_homotypic_adjusted <- expected_doublets_adjusted
  obj$doublet_pK <- best_pK

  # 删除 DoubletFinder doublet。
  obj <- subset(obj, subset = doublet_status == "Singlet")
  stopifnot(ncol(obj) > 100)

  # DecontX 使用当前 singlet 原始 counts；粗聚类标签作为 z。
  singlet_counts <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = singlet_counts)
  )
  sce <- decontX::decontX(
    x = sce,
    assayName = "counts",
    z = as.character(obj$seurat_clusters),
    seed = 12345,
    verbose = TRUE
  )

  # 写入 DecontX contamination，并按 <=30% 保留。
  obj$decontX_contamination <- as.numeric(
    SummarizedExperiment::colData(sce)$decontX_contamination
  )
  obj$decontX_pass_30pct <-
    obj$decontX_contamination <= decontx_max_contamination
  obj <- subset(obj, subset = decontX_pass_30pct)
  stopifnot(ncol(obj) > 100)

  # DecontX 结果允许非整数；按照官方说明四舍五入为整数 counts 后再进入 Seurat。
  decont_counts <- round(decontX::decontXcounts(sce))
  decont_counts <- decont_counts[, Seurat::Cells(obj), drop = FALSE]
  decont_counts <- as(decont_counts, "dgCMatrix")
  stopifnot(all(decont_counts@x >= 0))
  stopifnot(all(decont_counts@x == floor(decont_counts@x)))

  # 旧 nCount/nFeature 先另存为 preDecontX_*，新对象自动按校正 counts 重算。
  final_meta <- obj@meta.data
  final_meta$preDecontX_nCount_RNA <- final_meta$nCount_RNA
  final_meta$preDecontX_nFeature_RNA <- final_meta$nFeature_RNA
  final_meta$nCount_RNA <- NULL
  final_meta$nFeature_RNA <- NULL

  final_obj <- SeuratObject::CreateSeuratObject(
    counts = decont_counts,
    assay = "RNA",
    project = unique(final_meta$sample_uid),
    meta.data = final_meta,
    min.cells = 0,
    min.features = 0
  )

  # DecontX 后重新计算 QC 指标。
  final_obj$postDecontX_percent_mt <- Seurat::PercentageFeatureSet(
    final_obj,
    pattern = "^MT-"
  )
  final_obj$postDecontX_percent_ribo <- Seurat::PercentageFeatureSet(
    final_obj,
    pattern = "^RP[SL]"
  )

  # 给单样本对象生成 data layer，便于单独检查。
  final_obj <- Seurat::NormalizeData(
    final_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  final_obj
}


# -----------------------------------------------------------------------------
# 5. GSE255612：34 个 human donor + 完整官方 cell-level metadata
# -----------------------------------------------------------------------------
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

# 第一条数据为 TYPE/group/numeric 字段说明，不是真实细胞。
gse255_native_meta <- gse255_native_meta[-1, , drop = FALSE]

# 这是官方文件实际核心列名。
stopifnot(
  all(
    c(
      "NAME", "biosample_id", "donor_id", "cell_type",
      "sex", "af", "n_umi", "n_genes", "cellranger_percent_mito"
    ) %in% colnames(gse255_native_meta)
  )
)

# 34 个 GEO sample 与 sample_id 的固定对应关系。
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
  stringsAsFactors = FALSE
)
gse255_samples$donor_id <- paste0("P", gse255_samples$sample_id)
gse255_af_map <- c(case = "AF", control = "Control")

gse255_obj_list <- lapply(
  seq_len(nrow(gse255_samples)),
  function(index) {
    sample_row <- gse255_samples[index, , drop = FALSE]
    sample_dir <- file.path(gse255_raw_root, sample_row$sample_id)
    stopifnot(dir.exists(sample_dir))

    # SRA/Cell Ranger 重建后的单 donor 原始整数 10x counts。
    counts <- Seurat::Read10X(
      data.dir = sample_dir,
      gene.column = 2,
      unique.features = TRUE,
      strip.suffix = FALSE
    )
    counts <- as(counts, "dgCMatrix")

    # 官方 metadata 的 biosample_id 实际为 1279_3n 这类 sample_id。
    native_sample_meta <- gse255_native_meta[
      gse255_native_meta$biosample_id == sample_row$sample_id,
      ,
      drop = FALSE
    ]
    stopifnot(nrow(native_sample_meta) > 0)
    stopifnot(all(native_sample_meta$donor_id == sample_row$donor_id))

    # 官方 NAME 例如 AAAC...-1-0；去掉最后一个 -数字 后恢复单样本 10x barcode AAAC...-1。
    native_sample_meta$raw_barcode <- sub(
      "-[0-9]+$",
      "",
      native_sample_meta$NAME
    )
    native_index <- match(
      colnames(counts),
      native_sample_meta$raw_barcode
    )
    stopifnot(!anyNA(native_index))
    native_sample_meta <- native_sample_meta[native_index, , drop = FALSE]
    stopifnot(
      identical(
        native_sample_meta$raw_barcode,
        colnames(counts)
      )
    )

    new_cell_names <- paste0(
      "GSE255612__", sample_row$sample_id, "__", colnames(counts)
    )
    colnames(counts) <- new_cell_names

    obj <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_row$sample_id,
      min.cells = 0,
      min.features = 0
    )

    # 完整保留作者 metadata；raw_barcode 也保留用于追溯。
    native_meta_for_seurat <- native_sample_meta
    colnames(native_meta_for_seurat) <- paste0(
      "gse255_",
      make.names(colnames(native_meta_for_seurat))
    )
    rownames(native_meta_for_seurat) <- new_cell_names
    obj <- Seurat::AddMetaData(obj, native_meta_for_seurat)

    # 标准字段直接从官方 metadata 提取。
    condition_standard <- unname(gse255_af_map[native_sample_meta$af])
    stopifnot(!anyNA(condition_standard))

    obj <- add_standard_metadata(
      obj = obj,
      dataset = "GSE255612",
      GSM = sample_row$GSM,
      BioSample = NA_character_,
      sample_id = sample_row$sample_id,
      donor_id = sample_row$donor_id,
      tissue = "Left atrium",
      condition_original = as.character(native_sample_meta$af),
      AF_status = condition_standard,
      sex = as.character(native_sample_meta$sex),
      age_years = NA_real_,
      cell_type = as.character(native_sample_meta$cell_type),
      platform = "Illumina NovaSeq 6000",
      library_strategy = "RNA-Seq",
      technology = "10x Genomics single-nucleus RNA-seq",
      genome_build = "GRCh38",
      source_matrix_state = "raw_integer_counts_rebuilt_from_SRA"
    )

    run_qc_doubletfinder_decontx(obj, "GSE255612")
  }
)
names(gse255_obj_list) <- gse255_samples$sample_id


# -----------------------------------------------------------------------------
# 6. GSE198204：只取 6 个人类 Left atrial appendage 样本
# -----------------------------------------------------------------------------
gse198_tar <- file.path(gse198_dir, "GSE198204_RAW.tar")
stopifnot(file.exists(gse198_tar))

gse198_extract_dir <- file.path(work_dir, "GSE198204_extracted")
dir.create(gse198_extract_dir, recursive = TRUE, showWarnings = FALSE)
untar(gse198_tar, exdir = gse198_extract_dir)

# 只写当前已从 GEO 核实的信息；未公开确认的信息明确使用 NA，不推测。
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
    "GSM5940688", "GSM5940689", "GSM6552875",
    "GSM5940691", "GSM6552876", "GSM6552877"
  ),
  BioSample = c(
    "SAMN26540411", NA, "SAMN30711027",
    "SAMN26540408", "SAMN30711026", "SAMN30711025"
  ),
  sample_id = c(
    "LA_SR_01", "LA_SR_02", "LA_SR_03",
    "LA_AF_01", "LA_AF_02", "LA_AF_03"
  ),
  condition_original = c("SR", "SR", "SR", "AF", "AF", "AF"),
  AF_status = c("Control", "Control", "Control", "AF", "AF", "AF"),
  age_years = c(55, NA, 42, 55, 40, 62),
  stringsAsFactors = FALSE
)

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
    colnames(counts) <- paste0(
      "GSE198204__", sample_row$sample_id, "__", colnames(counts)
    )

    obj <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_row$sample_id,
      min.cells = 0,
      min.features = 0
    )

    # 本数据集没有额外 cell-level metadata 文件，GEO sample-level 信息复制到该样本全部细胞。
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

    run_qc_doubletfinder_decontx(obj, "GSE198204")
  }
)
names(gse198_obj_list) <- gse198_samples$sample_id


# -----------------------------------------------------------------------------
# 7. GSE238242：7 个人类 Multiome RNA/GEX donor + 完整作者 metadata
# -----------------------------------------------------------------------------
gse238_tar <- file.path(gse238_dir, "GSE238242_RAW.tar")
gse238_meta_file <- file.path(
  gse238_dir,
  "GSE238242_snAF.metadata.tsv.gz"
)
stopifnot(file.exists(gse238_tar))
stopifnot(file.exists(gse238_meta_file))

gse238_extract_dir <- file.path(work_dir, "GSE238242_extracted")
dir.create(gse238_extract_dir, recursive = TRUE, showWarnings = FALSE)
untar(gse238_tar, exdir = gse238_extract_dir)

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
  all(c("cell_type", "sample", "sex", "Rhythm") %in%
        colnames(gse238_native_meta))
)

gse238_samples <- data.frame(
  sample_id = c("CF69", "CF77", "CF89", "CF91", "CF93", "CF97", "CF102"),
  GSM = c(
    "GSM7660990", "GSM7660991", "GSM7660992", "GSM7660993",
    "GSM7660994", "GSM7660995", "GSM7660996"
  ),
  condition_original = c("SR", "SR", "SR", "SR", "AF", "AF", "AF"),
  AF_status = c("Control", "Control", "Control", "Control", "AF", "AF", "AF"),
  age_years = c(51, 78, 69, 57, 80, 58, 75),
  sex_expected = c("F", "M", "F", "M", "M", "M", "M"),
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

gse238_obj_list <- lapply(
  seq_len(nrow(gse238_samples)),
  function(index) {
    sample_row <- gse238_samples[index, , drop = FALSE]

    # GEO TSV 第一列是 gene，其余列为细胞 barcode；转换为稀疏 counts。
    dense_counts <- read.delim(
      gzfile(file.path(gse238_extract_dir, sample_row$counts_file)),
      header = TRUE,
      row.names = 1,
      sep = "\t",
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    counts <- Matrix::Matrix(as.matrix(dense_counts), sparse = TRUE)
    counts <- as(counts, "dgCMatrix")
    rownames(counts) <- make.unique(rownames(counts))
    rm(dense_counts)
    invisible(gc())

    # 只取当前 donor 的作者原始 metadata。
    native_sample_meta <- gse238_native_meta[
      gse238_native_meta$sample == sample_row$sample_id,
      ,
      drop = FALSE
    ]
    stopifnot(nrow(native_sample_meta) > 0)

    # counts 与 metadata 的末尾批次编号不同，去掉最后一个 -数字 后按 barcode 主体匹配。
    native_sample_meta$barcode_core <- sub(
      "-[0-9]+$",
      "",
      rownames(native_sample_meta)
    )
    count_barcode_core <- sub("-[0-9]+$", "", colnames(counts))
    native_index <- match(count_barcode_core, native_sample_meta$barcode_core)
    stopifnot(!anyNA(native_index))
    native_sample_meta <- native_sample_meta[native_index, , drop = FALSE]
    stopifnot(all(native_sample_meta$Rhythm == sample_row$condition_original))

    new_cell_names <- paste0(
      "GSE238242__", sample_row$sample_id, "__", colnames(counts)
    )
    colnames(counts) <- new_cell_names

    obj <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_row$sample_id,
      min.cells = 0,
      min.features = 0
    )

    # 完整保留作者原始 cell-level metadata，统一前缀避免字段冲突。
    native_meta_for_seurat <- native_sample_meta
    native_meta_for_seurat$barcode_core <- NULL
    colnames(native_meta_for_seurat) <- paste0(
      "gse238_",
      make.names(colnames(native_meta_for_seurat))
    )
    rownames(native_meta_for_seurat) <- new_cell_names
    obj <- Seurat::AddMetaData(obj, native_meta_for_seurat)

    # 统一标准字段。
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

    # 公开 sample-level sex 与 cell-level metadata 必须一致。
    standardized_sex <- toupper(substr(obj$sex, 1, 1))
    stopifnot(all(standardized_sex == sample_row$sex_expected))

    run_qc_doubletfinder_decontx(obj, "GSE238242")
  }
)
names(gse238_obj_list) <- gse238_samples$sample_id


# -----------------------------------------------------------------------------
# 8. 分别生成三个数据库的 human cleaned 对象
# -----------------------------------------------------------------------------
# merge.data=FALSE：只合并校正后的 counts；后面统一重新 NormalizeData。
gse255_obj <- merge(
  gse255_obj_list[[1]],
  y = gse255_obj_list[-1],
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "GSE255612_human_clean"
)
gse198_obj <- merge(
  gse198_obj_list[[1]],
  y = gse198_obj_list[-1],
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "GSE198204_human_clean"
)
gse238_obj <- merge(
  gse238_obj_list[[1]],
  y = gse238_obj_list[-1],
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "GSE238242_human_clean"
)

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
# 9. 三库合并 + Seurat v5 HarmonyIntegration
# -----------------------------------------------------------------------------
# Seurat v5 merge 保留每个 sample/library 的独立 layer；HarmonyIntegration 据此去批次。
seurobj <- merge(
  gse255_obj,
  y = list(gse198_obj, gse238_obj),
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "AF_3GEO_human"
)

stopifnot(all(seurobj$species == "Homo sapiens"))
stopifnot(
  setequal(
    unique(seurobj$dataset),
    c("GSE255612", "GSE198204", "GSE238242")
  )
)

# 各 layer 分别标准化，再共同选高变基因、Scale、PCA。
seurobj <- Seurat::NormalizeData(
  seurobj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
seurobj <- Seurat::FindVariableFeatures(
  seurobj,
  selection.method = "vst",
  nfeatures = integration_nfeatures,
  verbose = FALSE
)
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

# AF/Control 不进入 batch correction；batch 来自每个独立 sample layer。
seurobj <- Seurat::IntegrateLayers(
  object = seurobj,
  method = Seurat::HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

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

# 整合完成后合并 RNA layers，方便后续标准差异表达流程。
seurobj[["RNA"]] <- SeuratObject::JoinLayers(seurobj[["RNA"]])
seurobj$integration_method <- "Seurat_v5_HarmonyIntegration"
seurobj$decontX_threshold <- decontx_max_contamination
seurobj$pipeline_state <- "QC_DoubletFinder_DecontX30_Harmony_integrated"


# -----------------------------------------------------------------------------
# 10. 输出最终对象和完整标准化 metadata
# -----------------------------------------------------------------------------
final_rds <- file.path(
  output_dir,
  "AF_3GEO_human_QC_DoubletFinder_DecontX30_Harmony_seurobj.rds"
)
saveRDS(seurobj, final_rds, compress = FALSE)

metadata_output_file <- file.path(
  output_dir,
  "AF_3GEO_human_standardized_metadata.tsv.gz"
)
metadata_connection <- gzfile(metadata_output_file, open = "wt")
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
# 11. 结束时打印核查信息
# -----------------------------------------------------------------------------
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
  "\n最终 Seurat 对象：\n", final_rds,
  "\n\n最终 metadata：\n", metadata_output_file,
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# 12. 清理临时解压文件
# -----------------------------------------------------------------------------
unlink(work_dir, recursive = TRUE, force = TRUE)
