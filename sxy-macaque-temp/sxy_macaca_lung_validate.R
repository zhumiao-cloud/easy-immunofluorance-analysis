#!/usr/bin/env Rscript
# 食蟹猴肺单细胞/单核数据的 10x 原始计数校验。
# 输入：10x三件套目录、样本名、验证结果目录。
# 不重排矩阵、不进行归一化、不按基因符号合并重复项。
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
tenx <- normalizePath(args[[1]], mustWork = TRUE)
sample <- args[[2]]
out <- args[[3]]
featurefile <- file.path(tenx, "features.tsv.gz")
barcodefile <- file.path(tenx, "barcodes.tsv.gz")
matrixfile <- file.path(tenx, "matrix.mtx.gz")
stopifnot(all(file.exists(c(featurefile, barcodefile, matrixfile))))
# 第1列是稳定的基因ID；第2列是显示名；第3列是10x feature type。
features <- read.delim(gzfile(featurefile), header = FALSE, sep = "\t",
                       quote = "", comment.char = "", fill = FALSE,
                       stringsAsFactors = FALSE, check.names = FALSE)
barcodes <- readLines(gzfile(barcodefile), warn = FALSE)
stopifnot(ncol(features) == 3L, nrow(features) > 10000L, length(barcodes) > 0L)
stopifnot(!anyNA(features), !anyDuplicated(features[[1]]),
          all(nzchar(features[[1]])), all(nzchar(features[[2]])),
          all(features[[3]] == "Gene Expression"))
stopifnot(!anyNA(barcodes), !anyDuplicated(barcodes), all(nzchar(barcodes)))
# Matrix Market 的行号对应 features 行序，列号对应 barcodes 行序。
# 从原始MEX文件独立读取，避免仅凭维度相同错误判断对齐。
disk <- as(Matrix::readMM(gzfile(matrixfile)), "dgCMatrix")
stopifnot(identical(nrow(disk), nrow(features)),
          identical(ncol(disk), length(barcodes)))
stopifnot(all(is.finite(disk@x)), all(disk@x >= 0),
          all(disk@x == floor(disk@x)))
dimnames(disk) <- list(features[[1]], barcodes)
# 使用 gene.column=1 保留基因ID，不使用可能重复的 gene symbol。
m <- Seurat::Read10X(data.dir = tenx, gene.column = 1L,
                     unique.features = FALSE)
if (is.list(m)) {
  stopifnot(identical(names(m), "Gene Expression"))
  m <- m[["Gene Expression"]]
}
stopifnot(inherits(m, "dgCMatrix"))
stopifnot(identical(rownames(m), features[[1]]),
          identical(colnames(m), barcodes),
          identical(dim(m), dim(disk)))
# 精确比较非零坐标和值，任何行列错位或者表达量变化均终止。
stopifnot(identical(disk@p, m@p), identical(disk@i, m@i),
          identical(disk@x, m@x))
# 允许零表达的基因和低UMI细胞保留；下游质量控制另做。
obj <- Seurat::CreateSeuratObject(counts = m, project = sample,
                                  min.cells = 0, min.features = 0)
stopifnot(identical(rownames(obj), rownames(m)),
          identical(colnames(obj), colnames(m)))
obj$source_accession <- "CRA002577"
obj$organism <- "Macaca fascicularis"
obj$tissue <- "lung"
obj$sample_id <- sample
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rds <- file.path(out, "seurat_unfiltered_after_cellranger.rds")
saveRDS(obj, file = rds, compress = TRUE)
check <- readRDS(rds)
stopifnot(identical(rownames(check), rownames(obj)),
          identical(colnames(check), colnames(obj)))
summary <- data.frame(
  sample = sample,
  genes = nrow(m),
  barcodes = ncol(m),
  nonzero_entries = length(m@x),
  total_umis = sum(m@x),
  gene_id_order = "PASS",
  barcode_order = "PASS",
  mex_sparse_coordinates_and_counts = "PASS",
  seurat_rds_roundtrip = "PASS")
write.table(summary, file.path(out, "validation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("READ10X_SEURAT_PASS", sample, "GENES", nrow(m), "CELLS", ncol(m),
    "NNZ", length(m@x), "TOTAL_UMIS", sum(m@x), "\n")
