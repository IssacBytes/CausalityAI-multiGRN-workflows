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

## ---- 可视化：性别特异基因集的 ssGSEA 分数（官方 GSVA vignette 招牌验证图）----
## Pickrell/Huang 是 HapMap 淋巴母细胞系，无明显表型分组（热图意义不大）；
## 但有真实"性别"表型。官方用性别特异基因集验证：MSY（男性 Y 染色体基因）应在
## 男性高、XiE（X 失活逃逸基因）应在女性高。这是干净、可解读的展示级验证图。
fig_paths <- character(0)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  figs_dir <- file.path(results_dir, "figures")
  dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
  cjk_font <- if (.Platform$OS.type == "windows") "Microsoft YaHei" else "PingFang SC"  # §5
  qz <- capabilities("aqua")

  data(genderGenesEntrez)                                   # msYgenesEntrez / XiEgenesEntrez
  sexSets <- GeneSetCollection(
    GeneSet(as.character(msYgenesEntrez), geneIdType = EntrezIdentifier(), setName = "MSY"),
    GeneSet(as.character(XiEgenesEntrez), geneIdType = EntrezIdentifier(), setName = "XiE"))
  es_sex <- gsva(ssgseaParam(exprSet, sexSets, minSize = 2, maxSize = 500,
                             alpha = 0.25, normalize = TRUE,
                             annotation = NullIdentifier()), verbose = FALSE)
  sm <- exprs(es_sex)

  dfs <- data.frame(
    geneset = factor(rep(rownames(sm), times = ncol(sm)),
                     levels = c("MSY", "XiE"),
                     labels = c("MSY\n(男性 Y 染色体)", "XiE\n(X 失活逃逸)")),
    score   = as.numeric(sm),
    gender  = factor(rep(exprSet$Gender, each = nrow(sm)),
                     levels = c("Male", "Female"), labels = c("男性", "女性")))

  pal <- c(`男性` = "#2C7FB8", `女性` = "#D95F62")
  p <- ggplot(dfs, aes(geneset, score, fill = gender)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9,
                 position = position_dodge(0.66)) +
    scale_fill_manual(values = pal, name = "性别") +
    labs(title = "ssGSEA：性别特异基因集分数（按性别）",
         x = NULL, y = "ssGSEA 富集分数") +
    theme_classic(base_size = 13, base_family = cjk_font) +
    theme(plot.title  = element_text(face = "bold", size = 15),
          axis.text.x = element_text(size = 11, lineheight = 0.9),
          legend.position = "top",
          legend.title = element_text(face = "bold"))

  fp <- file.path(figs_dir, "ssgsea_sex_signature.png")
  ggsave(fp, p, width = 7, height = 5.2, dpi = 300,
         device = if (qz) "png" else NULL, type = if (qz) "quartz" else NULL)
  fig_paths <- c(fig_paths, fp)
} else {
  message("未安装 ggplot2，跳过可视化（不影响打分与导出）")
}

## ---- 摘要 ----
cat("\n=== ssGSEA 完成 / done ===\n")
cat(sprintf("通路 x 样本: %d x %d | 分数范围: [%.3f, %.3f]\n",
            nrow(es_mat), ncol(es_mat), min(es_mat), max(es_mat)))
cat("分数预览（前 3 通路 x 前 3 样本）:\n")
print(round(es_mat[1:3, 1:3], 3))
cat("\n结果已写入:\n")
cat("  ", file.path(results_dir, "ssgsea_pickrellhuang_scores.csv"), "\n")
cat("  ", file.path(results_dir, "ssgsea_scores_long.csv"), "\n")
for (fp in fig_paths) cat("  ", fp, "  ← 图\n")
