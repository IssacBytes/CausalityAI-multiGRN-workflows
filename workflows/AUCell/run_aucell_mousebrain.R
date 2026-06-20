#!/usr/bin/env Rscript
## =============================================================================
## AUCell 运行工作流 / AUCell run workflow
##
## 官方 AUCell vignette 的完整流程：GSE60361 小鼠脑单细胞 → 排名 → AUC → 阈值分配
## → 用论文给的真实细胞标签做混淆矩阵验证。
##
##   数据   : GSE60361（Zeisel et al. 2015，3005 细胞）—— 从 GEO 下载，自动缓存
##   基因集 : AUCell 包自带 inst/examples/geneSignatures.gmt（脑细胞类型签名）
##   标签   : AUCell 包自带 inst/examples/mouseBrain_cellLabels.tsv（真实细胞类型）
##   依赖   : AUCell, GSEABase, Matrix；下载需 GEOquery + data.table；阈值需 mixtools
##
## 用法（项目根目录下）：
##   Rscript workflows/AUCell/run_aucell_mousebrain.R
##
## 产出（写入 results/，已被 .gitignore 忽略）：
##   results/aucell_mousebrain_AUC.csv          AUC 分数（基因集 × 细胞）
##   results/aucell_scores_long.csv             统一长表（CLAUDE.md §7；含 binary/threshold）
##   results/aucell_confusion_matrix.csv        效度验证：分配 vs 真实细胞类型
##
## 注：首次运行从 GEO 下载 GSE60361（约几十 MB）并缓存到 data/；之后自动读缓存、不再联网。
## =============================================================================

## ---- 可移植路径（CLAUDE.md §5）----
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- "--file="
script_arg <- args[startsWith(args, file_arg)]
script_dir <- if (length(script_arg) > 0)
  dirname(normalizePath(sub(file_arg, "", script_arg[1]), winslash = "/")) else getwd()
project_root <- normalizePath(file.path(script_dir, "../.."), winslash = "/")
setwd(project_root)

local_lib <- Sys.getenv("MULTIGRN_RLIB", unset = file.path(project_root, "r-lib"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(AUCell)
  library(GSEABase)
  library(Matrix)
})

set.seed(123)

data_dir    <- "data"
results_dir <- "results"
dir.create(data_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 数据：自动缓存（首次下载，之后读缓存）----
cacheFile <- file.path(data_dir, "exprMatrix_MouseBrain.RData")
if (file.exists(cacheFile)) {
  message("读取本地缓存: ", cacheFile)
  load(cacheFile)                                  # 载入 mouseBrainExprMatrix
} else {
  message("首次运行：从 GEO 下载 GSE60361 ...")
  suppressPackageStartupMessages({ library(GEOquery); library(data.table) })
  geoFile <- tempfile(fileext = ".txt.gz")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60361/suppl/GSE60361_C1-3005-Expression.txt.gz",
    destfile = geoFile)
  exprMatrix <- fread(geoFile, sep = "\t")
  geneNames  <- unname(unlist(exprMatrix[, 1, with = FALSE]))
  exprMatrix <- as.matrix(exprMatrix[, -1, with = FALSE])
  rownames(exprMatrix) <- geneNames
  exprMatrix <- exprMatrix[unique(rownames(exprMatrix)), ]
  file.remove(geoFile)
  mouseBrainExprMatrix <- as(exprMatrix, "dgCMatrix")
  save(mouseBrainExprMatrix, file = cacheFile)     # 缓存，下次自动复用
}
exprMatrix <- mouseBrainExprMatrix
message(sprintf("表达矩阵: %d 基因 x %d 细胞", nrow(exprMatrix), ncol(exprMatrix)))

## ---- 基因集：包自带脑签名 ----
gmtFile  <- file.path(system.file("examples", package = "AUCell"), "geneSignatures.gmt")
geneSets <- getGmt(gmtFile)
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix))            # 过滤到数据中存在的基因
geneSets <- setGeneSetNames(geneSets,
              paste(names(geneSets), " (", nGenes(geneSets), "g)", sep = ""))

## ---- 三步打分 ----
message("AUCell: buildRankings → calcAUC → exploreThresholds ...")
cells_rankings   <- AUCell_buildRankings(exprMatrix, plotStats = FALSE, verbose = FALSE)
cells_AUC        <- AUCell_calcAUC(geneSets, cells_rankings, verbose = FALSE)
set.seed(333)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = FALSE,
                                             assign = TRUE, verbose = FALSE)

auc_mat   <- getAUC(cells_AUC)
selThr    <- getThresholdSelected(cells_assignment)                  # 各签名选定阈值

## ---- 效度验证：混淆矩阵（分配细胞 vs 真实细胞类型）----
labFile    <- file.path(system.file("examples", package = "AUCell"), "mouseBrain_cellLabels.tsv")
cellLabels <- read.table(labFile, row.names = 1, header = TRUE, sep = "\t")
cellTypes  <- unique(cellLabels[, "level1class"])
confMat <- t(sapply(cells_assignment, function(x)
  table(cellLabels[x$assignment, "level1class"])[cellTypes]))
colnames(confMat) <- cellTypes
confMat[is.na(confMat)] <- 0

## ---- 导出 ----
write.csv(auc_mat, file.path(results_dir, "aucell_mousebrain_AUC.csv"))
write.csv(confMat, file.path(results_dir, "aucell_confusion_matrix.csv"))

# 统一长表（§7）：AUCell 能填 binary（是否过阈值）与 threshold
assignedSets <- intersect(rownames(auc_mat), names(selThr))
long <- do.call(rbind, lapply(assignedSets, function(gs) {
  sc  <- auc_mat[gs, ]
  thr <- selThr[[gs]]
  data.frame(method = "AUCell", cell = colnames(auc_mat), signature = gs,
             score = as.numeric(sc),
             binary = as.integer(sc > thr),                          # 过阈值=1
             threshold = thr, row.names = NULL)
}))
write.csv(long, file.path(results_dir, "aucell_scores_long.csv"), row.names = FALSE)

## ---- 摘要 ----
cat("\n=== AUCell 完成 / done ===\n")
cat(sprintf("AUC 矩阵: %d 签名 x %d 细胞\n", nrow(auc_mat), ncol(auc_mat)))
cat("各签名分配的活跃细胞数:\n")
print(sapply(cells_assignment, function(x) length(x$assignment)))
cat("\n效度验证 — 混淆矩阵（行=AUCell 按签名分配的细胞, 列=真实细胞类型）:\n")
print(confMat)
cat("\n结果已写入:\n")
cat("  ", file.path(results_dir, "aucell_mousebrain_AUC.csv"), "\n")
cat("  ", file.path(results_dir, "aucell_scores_long.csv"), "\n")
cat("  ", file.path(results_dir, "aucell_confusion_matrix.csv"), "\n")
