#!/usr/bin/env Rscript
## =============================================================================
## ssGSEA 运行工作流 / ssGSEA run workflow
##
## ssGSEA = GSVA 包的 `ssgsea` method（见根目录 CLAUDE.md §1）——本脚本用
## GSVA::ssgseaParam() 在官方 GSVAdata 真实数据上跑通 ssGSEA，并导出分数。
##
##   数据   : GSVAdata::commonPickrellHuang（淋巴母细胞系；这里用微阵列 RMA eset）
##   基因集 : MSigDB C2 canonical 通路（KEGG / REACTOME / BIOCARTA）
##   依赖   : GSVA, GSVAdata, GSEABase, Biobase, BiocParallel（均 Bioconductor）
##
## 用法（项目根目录下）：
##   Rscript workflows/GSVA/run_ssgsea_pickrellhuang.R
##
## 产出（写入 results/，已被 .gitignore 忽略）：
##   results/ssgsea_pickrellhuang_scores.csv   宽表：通路 × 样本
##   results/ssgsea_scores_long.csv            统一长表（CLAUDE.md §7 schema）
## =============================================================================

## ---- 可移植路径（CLAUDE.md §5）：从脚本位置推导 project_root，不硬编码 ----
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- "--file="
script_arg <- args[startsWith(args, file_arg)]
script_dir <- if (length(script_arg) > 0)
  dirname(normalizePath(sub(file_arg, "", script_arg[1]), winslash = "/")) else getwd()
project_root <- normalizePath(file.path(script_dir, "../.."), winslash = "/")
setwd(project_root)

## R 库路径：优先环境变量，回退项目内 r-lib/
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

## OS 自适应并行后端（CLAUDE.md §5）
BPPARAM <- if (.Platform$OS.type == "windows")
  SnowParam(workers = 1) else MulticoreParam(workers = 2)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 官方数据（随 GSVAdata 安装，无需联网下载）----
data(commonPickrellHuang)   # huangArrayRMAnoBatchCommon_eset（微阵列 RMA，连续值）
data(c2BroadSets)           # MSigDB C2 curated gene sets（Entrez ID）

exprSet <- huangArrayRMAnoBatchCommon_eset

## MSigDB C2 canonical 通路：KEGG / REACTOME / BIOCARTA
canonicalC2BroadSets <- c2BroadSets[c(
  grep("^KEGG",     names(c2BroadSets)),
  grep("^REACTOME", names(c2BroadSets)),
  grep("^BIOCARTA", names(c2BroadSets)))]

message(sprintf("数据: %d 基因 x %d 样本 | 基因集: %d 条 canonical 通路",
                nrow(exprSet), ncol(exprSet), length(canonicalC2BroadSets)))

## ---- ssGSEA 打分 ----
## annotation = NullIdentifier()：eset 行名已是 Entrez ID（与 c2BroadSets 一致），
## 禁用不必要的 ID 映射，避免 "no method for coercing this S4 class to a vector"
## 报错（详见 workflows/GSVA/GSVA_学习手册.Rmd Ch4.3）。
ssgseaPar <- ssgseaParam(exprSet, canonicalC2BroadSets,
                         minSize    = 5,        # 过滤过小的通路
                         maxSize    = 500,      # 过滤过大的通路
                         alpha      = 0.25,     # 随机游走尾部权重（Barbie 2009 默认）
                         normalize  = TRUE,     # 按极差归一化
                         annotation = NullIdentifier())

es     <- gsva(ssgseaPar, BPPARAM = BPPARAM, verbose = FALSE)
es_mat <- exprs(es)

message(sprintf("ssGSEA 分数矩阵: %d 通路 x %d 样本", nrow(es_mat), ncol(es_mat)))

## ---- 导出：宽表 + 跨方法统一长表（CLAUDE.md §7）----
write.csv(es_mat, file.path(results_dir, "ssgsea_pickrellhuang_scores.csv"))

long <- data.frame(
  method    = "ssGSEA",
  cell      = rep(colnames(es_mat), each = nrow(es_mat)),   # 样本 ID
  signature = rep(rownames(es_mat), times = ncol(es_mat)),  # 通路名
  score     = as.numeric(es_mat),
  binary    = NA_integer_,                                  # 连续分，无阈值二值化
  threshold = NA_real_,
  row.names = NULL)
write.csv(long, file.path(results_dir, "ssgsea_scores_long.csv"), row.names = FALSE)

## ---- 摘要 ----
cat("\n=== ssGSEA 完成 / done ===\n")
cat(sprintf("通路 x 样本: %d x %d | 分数范围: [%.3f, %.3f]\n",
            nrow(es_mat), ncol(es_mat), min(es_mat), max(es_mat)))
cat("分数预览（前 3 通路 x 前 3 样本）:\n")
print(round(es_mat[1:3, 1:3], 3))
cat("\n结果已写入:\n")
cat("  ", file.path(results_dir, "ssgsea_pickrellhuang_scores.csv"), "\n")
cat("  ", file.path(results_dir, "ssgsea_scores_long.csv"), "\n")
