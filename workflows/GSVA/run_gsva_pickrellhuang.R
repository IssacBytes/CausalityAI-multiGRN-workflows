#!/usr/bin/env Rscript
## =============================================================================
## GSVA 运行工作流 / GSVA run workflow
##
## 用 GSVA::gsvaParam() 在官方 GSVAdata 真实数据上跑 GSVA method，重点演示 `kcdf`
## 随数据尺度的选择：同一批样本的微阵列（连续值→Gaussian）与 RNA-seq（计数→Poisson）。
##
##   数据   : GSVAdata::commonPickrellHuang（微阵列 RMA eset + RNA-seq 计数 eset）
##   基因集 : MSigDB C2 canonical 通路（KEGG / REACTOME / BIOCARTA）
##   依赖   : GSVA, GSVAdata, GSEABase, Biobase, BiocParallel
##
## 用法（项目根目录下）：
##   Rscript workflows/GSVA/run_gsva_pickrellhuang.R
##
## 产出（写入 results/，已被 .gitignore 忽略）：
##   results/gsva_huang_microarray_scores.csv   微阵列 GSVA 分数（通路 × 样本）
##   results/gsva_pickrell_rnaseq_scores.csv    RNA-seq GSVA 分数（通路 × 样本）
##   results/gsva_scores_long.csv               统一长表（CLAUDE.md §7 schema，微阵列）
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
  library(GSVA)
  library(GSVAdata)
  library(GSEABase)
  library(Biobase)
  library(BiocParallel)
})

set.seed(123)
BPPARAM <- if (.Platform$OS.type == "windows")
  SnowParam(workers = 1) else MulticoreParam(workers = 2)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 官方数据（随 GSVAdata 安装，无需联网）----
data(commonPickrellHuang)   # huangArrayRMAnoBatchCommon_eset / pickrellCountsArgonneCQNcommon_eset
data(c2BroadSets)

canonicalC2BroadSets <- c2BroadSets[c(
  grep("^KEGG",     names(c2BroadSets)),
  grep("^REACTOME", names(c2BroadSets)),
  grep("^BIOCARTA", names(c2BroadSets)))]

message(sprintf("基因集: %d 条 canonical 通路 | 样本: %d",
                length(canonicalC2BroadSets), ncol(huangArrayRMAnoBatchCommon_eset)))

## ---- GSVA 打分：★ kcdf 随数据尺度选 ----
## annotation = NullIdentifier()：eset 行名已是 Entrez（与 c2BroadSets 直配），
## 禁用不必要的 ID 映射，避免 AnnotationIdentifier 报错（见 GSVA 学习稿 Ch4.3）。

# 微阵列 RMA（连续值）→ kcdf 默认 "Gaussian"
message("GSVA on 微阵列 (kcdf=Gaussian) ...")
esmicro <- gsva(gsvaParam(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets,
                          minSize = 5, maxSize = 500,         # kcdf 取默认 Gaussian
                          annotation = NullIdentifier()),
                BPPARAM = BPPARAM, verbose = FALSE)

# RNA-seq 整数计数 → kcdf = "Poisson"
message("GSVA on RNA-seq (kcdf=Poisson) ...")
esrnaseq <- gsva(gsvaParam(pickrellCountsArgonneCQNcommon_eset, canonicalC2BroadSets,
                           minSize = 5, maxSize = 500, kcdf = "Poisson",
                           annotation = NullIdentifier()),
                 BPPARAM = BPPARAM, verbose = FALSE)

m_micro <- exprs(esmicro)
m_rna   <- exprs(esrnaseq)

## ---- 跨平台一致性（官方 vignette 的验证：通路分数应跨平台高相关）----
common_paths <- intersect(rownames(m_micro), rownames(m_rna))
pwy_corr <- sapply(common_paths, function(p)
  cor(m_micro[p, ], m_rna[p, ], method = "spearman"))

## ---- 导出 ----
write.csv(m_micro, file.path(results_dir, "gsva_huang_microarray_scores.csv"))
write.csv(m_rna,   file.path(results_dir, "gsva_pickrell_rnaseq_scores.csv"))

long <- data.frame(
  method    = "GSVA",
  cell      = rep(colnames(m_micro), each = nrow(m_micro)),   # 样本 ID（微阵列）
  signature = rep(rownames(m_micro), times = ncol(m_micro)),
  score     = as.numeric(m_micro),
  binary    = NA_integer_,
  threshold = NA_real_,
  row.names = NULL)
write.csv(long, file.path(results_dir, "gsva_scores_long.csv"), row.names = FALSE)

## ---- 摘要 ----
cat("\n=== GSVA 完成 / done ===\n")
cat(sprintf("微阵列(Gaussian): %d 通路 x %d 样本 | 分数范围 [%.3f, %.3f]\n",
            nrow(m_micro), ncol(m_micro), min(m_micro), max(m_micro)))
cat(sprintf("RNA-seq(Poisson): %d 通路 x %d 样本 | 分数范围 [%.3f, %.3f]\n",
            nrow(m_rna), ncol(m_rna), min(m_rna), max(m_rna)))
cat(sprintf("跨平台通路分数 Spearman 相关（%d 条共有通路）: 中位数 %.3f\n",
            length(pwy_corr), median(pwy_corr)))
cat("→ 中位为正 = 两平台通路分数总体同向（GSVA 把基因表达聚合到通路层面；\n")
cat("  官方 vignette 的论点是通路层面的跨平台一致性与基因层面相当，并非完美相关）。\n")
cat("\n结果已写入:\n")
cat("  ", file.path(results_dir, "gsva_huang_microarray_scores.csv"), "\n")
cat("  ", file.path(results_dir, "gsva_pickrell_rnaseq_scores.csv"), "\n")
cat("  ", file.path(results_dir, "gsva_scores_long.csv"), "\n")
