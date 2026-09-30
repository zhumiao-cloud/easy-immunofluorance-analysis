# =============================================================================
# 文件名：build_AF_3GEO_SeuratObj.R
#
# 目的：
#   将 GSE255612、GSE198204、GSE238242 整理为 Seurat v5 对象，
#   并生成三套人类 RNA/snRNA 数据的汇总对象 seurobj。
#
# 使用方式：
#   1. 把本脚本放在“AF原始数据”目录下。
#   2. 在 R / RStudio 中把工作目录切换到该目录，然后 source 本脚本；
#      或在 shell 中进入该目录后执行：
#      Rscript build_AF_3GEO_SeuratObj.R
#   3. 如脚本不放在数据根目录，请手动修改第 1 节的 base_dir。
#
# 数据处理原则：
#   - 不进行 QC，不删除任何细胞，不按 nFeature/nCount/线粒体比例过滤。
#   - 不执行 Harmony、CCA、RPCA、IntegrateLayers 等批次校正。
#   - GSE255612 是作者公开的 Scanpy processed real-valued expression matrix，
#     只作为 Seurat v5 data layer 使用，不伪装成原始 UMI counts。
#   - GSE198204 与 GSE238242 的 RNA 是 counts；保留 counts，同时生成 LogNormalize data。
#   - GSE198204 的小鼠 Sham/TAC 单独保存，不与人数据混合。
#   - GSE238242 的 ATAC 不加入三数据库 RNA 汇总对象；原始 ATAC 文件保持不动。
#   - 本脚本没有条件兜底或异常吞噬逻辑；任何结构不符合预期都会直接停止。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 加载依赖包并固定 Seurat v5 数据结构
# -----------------------------------------------------------------------------
# Seurat：创建、标准化、合并 Seurat 对象。
# SeuratObject：直接操作 Assay5 和 layer。
# Matrix：读取 Matrix Market 稀疏矩阵并保持稀疏存储。
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# GSE255612 需要 data-only Assay5，因此要求 Seurat/SeuratObject >= 5.0。
stopifnot(packageVersion("Seurat") >= "5.0.0")
stopifnot(packageVersion("SeuratObject") >= "5.0.0")

# 明确让新建 RNA assay 使用 Seurat v5 Assay5。
options(Seurat.object.assay.version = "v5")


# -----------------------------------------------------------------------------
# 1. 设置数据根目录、临时目录和输出目录
# -----------------------------------------------------------------------------
# 默认认为本脚本就在 AF原始数据 目录，并且从该目录启动 R。
# 如需从其他位置运行，只修改下面这一行即可。
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# 三个 GEO 数据集所在目录。
gse255_dir <- file.path(base_dir, "GSE255612")
gse198_dir <- file.path(base_dir, "GSE198204")
gse238_dir <- file.path(base_dir, "GSE238242")

# 临时解压目录放在 R 会话临时目录中，不写入或删除三个原始 GEO 目录。
work_dir <- file.path(tempdir(), "AF_3GEO_Seurat_build")
unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# 所有 RDS 结果统一写入 AF原始数据/SeuratObj_outputs。
output_dir <- file.path(base_dir, "SeuratObj_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. 函数：读取 GSE198204 的一个标准 10x RNA 样本
# -----------------------------------------------------------------------------
# 输入：
#   extracted_dir       RAW.tar 解压目录。
#   file_prefix         GEO 文件前缀。
#   sample_id           简化样本名。
#   condition_original  SR / AF / Sham / TAC。
#   species_group       Human / Mouse。
#
# 输出：
#   一个不做 QC、不删除细胞的 Seurat v5 RNA counts 对象。
read_gse198_sample <- function(
  extracted_dir,
  file_prefix,
  sample_id,
  condition_original,
  species_group
) {
  # 按 GEO 实际文件结构读取 matrix + barcodes + features。
  counts <- Seurat::ReadMtx(
    mtx = file.path(
      extracted_dir,
      paste0(file_prefix, "_matrix.mtx.gz")
    ),
    cells = file.path(
      extracted_dir,
      paste0(file_prefix, "_barcodes.tsv.gz")
    ),
    features = file.path(
      extracted_dir,
      paste0(file_prefix, "_features.tsv.gz")
    ),
    cell.column = 1,
    feature.column = 2,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  # 添加“数据库 + 样本”前缀，避免不同数据库出现相同 10x barcode。
  new_cell_names <- paste0(
    "GSE198204__",
    sample_id,
    "__",
    colnames(counts)
  )
  colnames(counts) <- new_cell_names

  # 建立固定的样本级 metadata；不依赖表达阈值。
  meta <- data.frame(
    dataset = rep("GSE198204", ncol(counts)),
    sample_id = rep(sample_id, ncol(counts)),
    condition_original = rep(condition_original, ncol(counts)),
    species_group = rep(species_group, ncol(counts)),
    source_matrix_state = rep("raw_10x_counts", ncol(counts)),
    stringsAsFactors = FALSE,
    row.names = colnames(counts)
  )

  # 统一人类 AF/SR 标签；鼠 Sham/TAC 在 AF_status 中保持 NA。
  status_map <- c(
    AF = "AF",
    SR = "Control",
    Sham = NA_character_,
    TAC = NA_character_
  )
  meta$AF_status <- unname(
    status_map[meta$condition_original]
  )

  # counts 列名与 metadata 行名必须严格逐一一致。
  stopifnot(identical(colnames(counts), rownames(meta)))

  # 创建 Seurat 对象；min.cells/min.features 均为 0，不删除任何 feature/cell。
  obj <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = sample_id,
    meta.data = meta,
    min.cells = 0,
    min.features = 0
  )

  return(obj)
}


# -----------------------------------------------------------------------------
# 3. 函数：读取 GSE238242 的一个 RNA counts.tsv.gz 样本
# -----------------------------------------------------------------------------
# 输入：
#   counts_file      一个 donor 的 gene × cell RNA counts TSV。
#   global_meta      GEO 提供的 GSE238242_snAF.metadata.tsv.gz。
#   sample_id        CF69 / CF77 / CF89 / CF91 / CF93 / CF97 / CF102。
#   expected_rhythm  SR / AF。
#
# 输出：
#   一个带作者 cell_type、sex、Rhythm 等 metadata 的 Seurat RNA counts 对象。
read_gse238_rna_sample <- function(
  counts_file,
  global_meta,
  sample_id,
  expected_rhythm
) {
  # 文件第一列是基因名，其余列是细胞 barcode；直接把第一列设为 row.names。
  dense_counts <- read.delim(
    gzfile(counts_file),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # 转成稀疏 dgCMatrix，避免后续对象持续占用密集矩阵内存。
  counts <- Matrix::Matrix(
    as.matrix(dense_counts),
    sparse = TRUE
  )
  counts <- as(counts, "dgCMatrix")

  # 重复 gene symbol 只做唯一化命名，不删除任何 feature。
  rownames(counts) <- make.unique(rownames(counts))

  # 密集临时矩阵完成转换后立即释放。
  rm(dense_counts)
  invisible(gc())

  # 保存 counts 文件中的原始 barcode。
  original_barcodes <- colnames(counts)

  # 固定取当前 donor 的全局 metadata。
  sample_meta <- global_meta[
    global_meta$sample == sample_id,
    ,
    drop = FALSE
  ]

  # GEO 全局 metadata barcode 含聚合数字后缀；
  # 在当前 donor 内去掉末尾“-数字”后再与 counts barcode 一一对应。
  sample_meta$barcode_core <- sub(
    "-[0-9]+$",
    "",
    rownames(sample_meta)
  )
  count_barcode_core <- sub(
    "-[0-9]+$",
    "",
    original_barcodes
  )

  # 按 barcode 主体构建固定索引。
  meta_index <- match(
    count_barcode_core,
    sample_meta$barcode_core
  )

  # 每个 counts barcode 都必须找到 metadata，且 Rhythm 必须与样本表一致。
  stopifnot(!anyNA(meta_index))
  sample_meta <- sample_meta[
    meta_index,
    ,
    drop = FALSE
  ]
  stopifnot(all(sample_meta$Rhythm == expected_rhythm))

  # 添加“数据库 + donor”前缀，保证跨数据库 cell name 唯一。
  new_cell_names <- paste0(
    "GSE238242__",
    sample_id,
    "__",
    original_barcodes
  )
  colnames(counts) <- new_cell_names
  rownames(sample_meta) <- new_cell_names

  # 添加三数据库统一 metadata 字段，同时保留作者原始 metadata。
  sample_meta$dataset <- "GSE238242"
  sample_meta$sample_id <- sample_id
  sample_meta$condition_original <- as.character(sample_meta$Rhythm)
  sample_meta$species_group <- "Human"
  sample_meta$source_matrix_state <- "raw_RNA_counts"

  status_map <- c(
    AF = "AF",
    SR = "Control"
  )
  sample_meta$AF_status <- unname(
    status_map[sample_meta$condition_original]
  )

  # 删除仅用于内部匹配的临时字段。
  sample_meta$barcode_core <- NULL

  # counts 与 metadata 必须严格按同一细胞顺序排列。
  stopifnot(identical(colnames(counts), rownames(sample_meta)))

  # 创建 Seurat 对象，不设置任何 QC 过滤阈值。
  obj <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = sample_id,
    meta.data = sample_meta,
    min.cells = 0,
    min.features = 0
  )

  return(obj)
}


# -----------------------------------------------------------------------------
# 4. GSE255612：读取作者公开的 processed snRNA expression
# -----------------------------------------------------------------------------
# 已核对的数据结构：
#   Matrix：36600 genes × 179697 cells
#   barcodes：179697 行
#   genes：36600 行；第 1 列 Ensembl ID，第 2 列 gene symbol
#   metadata：第 1 个数据行是 TYPE 描述行，之后才是 179697 个细胞
#
# Matrix Market 文件含非整数实数，因此不作为 UMI counts 使用。
gse255_matrix_file <- file.path(
  gse255_dir,
  "GSE255612_AF_snRNA_Matrix_V1.mtx.gz"
)
gse255_barcode_file <- file.path(
  gse255_dir,
  "GSE255612_AF_snRNA_Processed_Expression_Matrix_barcodes_V1.tsv.gz"
)
gse255_gene_file <- file.path(
  gse255_dir,
  "GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz"
)
gse255_meta_file <- file.path(
  gse255_dir,
  "GSE255612_AF_snRNA_MetaData.txt.gz"
)

# 四个必须文件缺失时直接停止。
stopifnot(file.exists(gse255_matrix_file))
stopifnot(file.exists(gse255_barcode_file))
stopifnot(file.exists(gse255_gene_file))
stopifnot(file.exists(gse255_meta_file))

# 读取 Matrix Market 稀疏矩阵。
gse255_matrix <- Matrix::readMM(
  gzfile(gse255_matrix_file)
)
gse255_matrix <- as(
  gse255_matrix,
  "dgCMatrix"
)

# 读取无表头 barcode 文件。
gse255_barcodes <- read.delim(
  gzfile(gse255_barcode_file),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

# 读取无表头 gene 文件。
gse255_genes <- read.delim(
  gzfile(gse255_gene_file),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

# 按已核对的官方实际规模做严格检查。
stopifnot(nrow(gse255_matrix) == 36600)
stopifnot(ncol(gse255_matrix) == 179697)
stopifnot(nrow(gse255_genes) == 36600)
stopifnot(nrow(gse255_barcodes) == 179697)
stopifnot(ncol(gse255_genes) >= 2)

# 第 2 列 gene symbol 作为 feature name；重复名称只唯一化，不删 feature。
rownames(gse255_matrix) <- make.unique(
  as.character(gse255_genes[[2]])
)

# barcode 文件第 1 列作为细胞名。
colnames(gse255_matrix) <- as.character(
  gse255_barcodes[[1]]
)

# 读取作者提供的 cell-level metadata。
gse255_meta <- read.delim(
  gzfile(gse255_meta_file),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 第一个数据行是 TYPE/group/numeric 的字段说明，不是细胞，固定删除。
gse255_meta <- gse255_meta[
  -1,
  ,
  drop = FALSE
]

# NAME 为原始细胞 barcode。
rownames(gse255_meta) <- gse255_meta$NAME

# 严格按表达矩阵 cell 顺序重排 metadata。
gse255_meta <- gse255_meta[
  colnames(gse255_matrix),
  ,
  drop = FALSE
]

# matrix 与 metadata 必须一一对应。
stopifnot(nrow(gse255_meta) == 179697)
stopifnot(
  identical(
    rownames(gse255_meta),
    colnames(gse255_matrix)
  )
)

# 添加数据库前缀，避免合并后三套数据 cell name 冲突。
gse255_new_cell_names <- paste0(
  "GSE255612__",
  colnames(gse255_matrix)
)
colnames(gse255_matrix) <- gse255_new_cell_names
rownames(gse255_meta) <- gse255_new_cell_names

# 添加统一 metadata 字段。
gse255_meta$dataset <- "GSE255612"
gse255_meta$sample_id <- as.character(
  gse255_meta$biosample_id
)

gse255_condition_map <- c(
  case = "AF",
  control = "NF"
)
gse255_status_map <- c(
  case = "AF",
  control = "Control"
)

gse255_meta$condition_original <- unname(
  gse255_condition_map[gse255_meta$af]
)
gse255_meta$AF_status <- unname(
  gse255_status_map[gse255_meta$af]
)
gse255_meta$species_group <- "Human"
gse255_meta$source_matrix_state <- "GEO_processed_Scanpy_expression"

# case/control 必须全部成功映射。
stopifnot(!anyNA(gse255_meta$condition_original))
stopifnot(!anyNA(gse255_meta$AF_status))

# processed expression 只进入 data layer，不人为生成 counts layer。
gse255_assay <- SeuratObject::CreateAssay5Object(
  data = gse255_matrix
)

# 创建 GSE255612 Seurat 对象。
gse255_obj <- SeuratObject::CreateSeuratObject(
  gse255_assay,
  assay = "RNA",
  project = "GSE255612",
  meta.data = gse255_meta
)

# 单独保存该数据库对象。
saveRDS(
  gse255_obj,
  file = file.path(
    output_dir,
    "GSE255612_processed_RNA_Seurat.rds"
  ),
  compress = FALSE
)


# -----------------------------------------------------------------------------
# 5. GSE198204：解压 8 个标准 10x 样本并创建 Seurat 对象
# -----------------------------------------------------------------------------
# Human：
#   SR：LA_SR_01、LA_SR_02、LA_SR_03
#   AF：LA_AF_01、LA_AF_02、LA_AF_03
#
# Mouse：
#   Sham_RFP_cells、TAC_RFP_cells
#
# 人与鼠分别保存，鼠数据不进入最终三数据库 human seurobj。
gse198_tar <- file.path(
  gse198_dir,
  "GSE198204_RAW.tar"
)
stopifnot(file.exists(gse198_tar))

# 创建临时解压目录。
gse198_extract_dir <- file.path(
  work_dir,
  "GSE198204_extracted"
)
dir.create(
  gse198_extract_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 解压 supplementary TAR；原始 tar 不修改。
untar(
  gse198_tar,
  exdir = gse198_extract_dir
)

# 固定列出实际 8 个 GEO 样本及分组。
gse198_samples <- data.frame(
  file_prefix = c(
    "GSM5940688_LA_SR_01",
    "GSM5940689_LA_SR_02",
    "GSM6552875_LA_SR_03",
    "GSM5940691_LA_AF_01",
    "GSM6552876_LA_AF_02",
    "GSM6552877_LA_AF_03",
    "GSM5940692_Sham_RFP_cells",
    "GSM5940693_TAC_RFP_cells"
  ),
  sample_id = c(
    "LA_SR_01",
    "LA_SR_02",
    "LA_SR_03",
    "LA_AF_01",
    "LA_AF_02",
    "LA_AF_03",
    "Sham_RFP_cells",
    "TAC_RFP_cells"
  ),
  condition_original = c(
    "SR",
    "SR",
    "SR",
    "AF",
    "AF",
    "AF",
    "Sham",
    "TAC"
  ),
  species_group = c(
    "Human",
    "Human",
    "Human",
    "Human",
    "Human",
    "Human",
    "Mouse",
    "Mouse"
  ),
  stringsAsFactors = FALSE
)

# 按固定样本表逐一建立样本级 Seurat 对象。
gse198_obj_list <- lapply(
  seq_len(nrow(gse198_samples)),
  function(index) {
    read_gse198_sample(
      extracted_dir = gse198_extract_dir,
      file_prefix = gse198_samples$file_prefix[index],
      sample_id = gse198_samples$sample_id[index],
      condition_original = gse198_samples$condition_original[index],
      species_group = gse198_samples$species_group[index]
    )
  }
)
names(gse198_obj_list) <- gse198_samples$sample_id

# 前 6 个固定为人类 AF/SR 样本，合并为一个 human RNA 对象。
gse198_human_obj <- Reduce(
  function(x, y) {
    merge(
      x = x,
      y = y,
      merge.data = FALSE,
      merge.dr = FALSE,
      project = "GSE198204_human"
    )
  },
  gse198_obj_list[1:6]
)

# 对 human counts 进行 LogNormalize；counts 原样保留。
gse198_human_obj <- Seurat::NormalizeData(
  object = gse198_human_obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# 后 2 个固定为小鼠 Sham/TAC；单独合并，不进入 human 汇总对象。
gse198_mouse_obj <- Reduce(
  function(x, y) {
    merge(
      x = x,
      y = y,
      merge.data = FALSE,
      merge.dr = FALSE,
      project = "GSE198204_mouse"
    )
  },
  gse198_obj_list[7:8]
)

# 分别保存人和鼠对象。
saveRDS(
  gse198_human_obj,
  file = file.path(
    output_dir,
    "GSE198204_human_RNA_Seurat.rds"
  ),
  compress = FALSE
)

saveRDS(
  gse198_mouse_obj,
  file = file.path(
    output_dir,
    "GSE198204_mouse_RNA_Seurat.rds"
  ),
  compress = FALSE
)


# -----------------------------------------------------------------------------
# 6. GSE238242：读取 7 个 human snRNA counts 样本
# -----------------------------------------------------------------------------
# 该数据集为 10x Multiome。
# 最终三数据库共同模态只取 RNA：
#   SR：CF69、CF77、CF89、CF91
#   AF：CF93、CF97、CF102
#
# 对应 ATAC 文件不删除、不修改，也不放入 RNA assay。
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

# 创建临时解压目录。
gse238_extract_dir <- file.path(
  work_dir,
  "GSE238242_extracted"
)
dir.create(
  gse238_extract_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 解压 TAR；原始 RNA/ATAC supplementary 文件不修改。
untar(
  gse238_tar,
  exdir = gse238_extract_dir
)

# 读取作者提供的全局 metadata；第一列作为 cell row.names。
gse238_meta <- read.delim(
  gzfile(gse238_meta_file),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 核对后续匹配所需的核心字段。
stopifnot(
  all(
    c(
      "cell_type",
      "sample",
      "sex",
      "Rhythm"
    ) %in% colnames(gse238_meta)
  )
)

# 固定列出 7 个 RNA 文件与分组。
gse238_samples <- data.frame(
  sample_id = c(
    "CF69",
    "CF77",
    "CF89",
    "CF91",
    "CF93",
    "CF97",
    "CF102"
  ),
  expected_rhythm = c(
    "SR",
    "SR",
    "SR",
    "SR",
    "AF",
    "AF",
    "AF"
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

# 逐 donor 创建 RNA Seurat 对象。
gse238_obj_list <- lapply(
  seq_len(nrow(gse238_samples)),
  function(index) {
    read_gse238_rna_sample(
      counts_file = file.path(
        gse238_extract_dir,
        gse238_samples$counts_file[index]
      ),
      global_meta = gse238_meta,
      sample_id = gse238_samples$sample_id[index],
      expected_rhythm = gse238_samples$expected_rhythm[index]
    )
  }
)
names(gse238_obj_list) <- gse238_samples$sample_id

# 合并 7 个 donor，保留各自来源的 counts layer。
gse238_rna_obj <- Reduce(
  function(x, y) {
    merge(
      x = x,
      y = y,
      merge.data = FALSE,
      merge.dr = FALSE,
      project = "GSE238242_RNA"
    )
  },
  gse238_obj_list
)

# 对原始 RNA counts 执行 LogNormalize；不做 QC/聚类/细胞删除。
gse238_rna_obj <- Seurat::NormalizeData(
  object = gse238_rna_obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# 保存 GSE238242 RNA 对象。
saveRDS(
  gse238_rna_obj,
  file = file.path(
    output_dir,
    "GSE238242_human_RNA_Seurat.rds"
  ),
  compress = FALSE
)


# -----------------------------------------------------------------------------
# 7. 合并三套 human RNA/snRNA 数据，生成最终 seurobj
# -----------------------------------------------------------------------------
# 最终对象包含：
#   GSE255612：18 AF + 16 Control donor，processed snRNA expression
#   GSE198204：3 AF + 3 SR donor，raw counts + LogNormalize data
#   GSE238242：3 AF + 4 SR donor，raw counts + LogNormalize data
#
# Seurat v5 merge 会把不同来源表达信息保留在分层 layer 中。
# 这里不执行 JoinLayers，也不执行 batch correction。
# 原因是 GSE255612 与另外两套数据库公开的数据层级不同。
seurobj <- merge(
  x = gse255_obj,
  y = list(
    gse198_human_obj,
    gse238_rna_obj
  ),
  merge.data = TRUE,
  merge.dr = FALSE,
  project = "AF_3GEO_human_RNA"
)

# 明确记录最终对象目前只是汇总，并未做批次校正。
seurobj$three_GEO_merge_state <- "merged_not_batch_corrected"

# 最终对象只能包含人类细胞。
stopifnot(
  all(
    seurobj$species_group == "Human"
  )
)

# 最终 metadata 必须同时存在三个数据库来源。
stopifnot(
  setequal(
    unique(seurobj$dataset),
    c(
      "GSE255612",
      "GSE198204",
      "GSE238242"
    )
  )
)

# 保存最终三数据库 human RNA/snRNA 汇总对象。
saveRDS(
  seurobj,
  file = file.path(
    output_dir,
    "AF_3GEO_human_RNA_seurobj.rds"
  ),
  compress = FALSE
)


# -----------------------------------------------------------------------------
# 8. 打印最终对象摘要，并清理脚本产生的临时解压目录
# -----------------------------------------------------------------------------
# 仅打印结构，不执行 QC、降维、聚类或批次校正。
print(seurobj)
print(table(seurobj$dataset))
print(
  table(
    seurobj$dataset,
    seurobj$AF_status,
    useNA = "always"
  )
)
print(
  SeuratObject::Layers(
    seurobj[["RNA"]]
  )
)

# 打印最终 RDS 路径。
cat(
  "\n最终对象已写入：\n",
  file.path(
    output_dir,
    "AF_3GEO_human_RNA_seurobj.rds"
  ),
  "\n",
  sep = ""
)

# 删除本脚本在 R 临时目录中创建的解压文件；
# 三个 GEO 原始目录和 SeuratObj_outputs 完全不动。
unlink(
  work_dir,
  recursive = TRUE,
  force = TRUE
)
