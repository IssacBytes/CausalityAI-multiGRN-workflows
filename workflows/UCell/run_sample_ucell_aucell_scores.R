## Minimal example for calling local UCell and AUCell packages.
## Run from the project root with:
##   Rscript workflows/UCell/run_sample_ucell_aucell_scores.R

# 可移植 R 库路径（CLAUDE.md §5）：优先环境变量，回退项目内 r-lib/
local_lib <- Sys.getenv("MULTIGRN_RLIB", unset = file.path(getwd(), "r-lib"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

library(UCell)
library(AUCell)
library(BiocParallel)

data(sample.matrix, package = "UCell")

gene_sets <- list(
  Tcell = c("CD2", "CD3E", "CD3D"),
  Myeloid = c("SPI1", "FCER1G", "CSF1R")
)

## UCell: directly scores a genes x cells expression matrix.
ucell_scores <- ScoreSignatures_UCell(
  sample.matrix,
  features = gene_sets,
  BPPARAM = SnowParam(workers = 1)
)

cat("\nUCell scores:\n")
print(head(ucell_scores))

## AUCell: first builds per-cell rankings, then calculates AUC scores.
set.seed(1)
rankings <- AUCell_buildRankings(
  as.matrix(sample.matrix),
  plotStats = FALSE,
  verbose = FALSE
)

auc <- AUCell_calcAUC(
  gene_sets,
  rankings,
  aucMaxRank = ceiling(0.05 * nrow(sample.matrix)),
  verbose = FALSE
)

auc_scores <- getAUC(auc)

cat("\nAUCell scores:\n")
print(auc_scores[, 1:3])

## ---- 可视化：UCell 分数按签名的小提琴图（官方 UCell vignette 用法）----
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("reshape2", quietly = TRUE)) {
  library(ggplot2)
  figs_dir <- file.path("results", "figures")
  dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
  cjk_font <- if (.Platform$OS.type == "windows") "Microsoft YaHei" else "PingFang SC"  # §5
  qz <- capabilities("aqua")

  dfu <- reshape2::melt(as.matrix(ucell_scores),
                        varnames = c("cell", "signature"), value.name = "UCell")
  dfu$signature <- sub("_UCell$", "", dfu$signature)        # 去掉列名后缀

  p <- ggplot(dfu, aes(signature, UCell, fill = signature)) +
    geom_violin(scale = "width", alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.9) +
    geom_jitter(width = 0.08, size = 0.5, alpha = 0.5) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(title = "UCell：各签名分数分布",
         x = "签名", y = "UCell 分数 (0–1)") +
    theme_bw(base_size = 11, base_family = cjk_font) +
    theme(plot.title = element_text(face = "bold"))

  fp <- file.path(figs_dir, "ucell_score_violin.png")
  ggsave(fp, p, width = 6, height = 4.5, dpi = 300,
         device = if (qz) "png" else NULL, type = if (qz) "quartz" else NULL)
  cat("\n图已写入:", fp, "  ← UCell 分数小提琴图\n")

  ## ---- 可视化2：降维特征图，按 UCell 分数着色（官方 FeaturePlot 等价图）----
  ## sample.matrix 无现成降维坐标，且仅 30 细胞（UMAP/tSNE 不适用）→ 用 PCA 算 2D 嵌入。
  logm <- log1p(as.matrix(sample.matrix))                     # log 归一化
  vg   <- order(apply(logm, 1, var), decreasing = TRUE)[seq_len(min(1000, nrow(logm)))]
  pca  <- prcomp(t(logm[vg, ]), scale. = FALSE)
  emb  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     cell = colnames(sample.matrix))
  dfp  <- merge(dfu, emb, by = "cell")

  p2 <- ggplot(dfp[order(dfp$UCell), ], aes(PC1, PC2, color = UCell)) +
    geom_point(size = 2.6, alpha = 0.9) +
    facet_wrap(~ signature) +
    scale_color_gradient(low = "grey85", high = "#D7301F", name = "UCell\n分数") +
    labs(title = "UCell：PCA 嵌入上的签名活性",
         x = "PC1", y = "PC2") +
    theme_minimal(base_size = 11, base_family = cjk_font) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  fp2 <- file.path(figs_dir, "ucell_pca_featureplot.png")
  ggsave(fp2, p2, width = 7.5, height = 4, dpi = 300,
         device = if (qz) "png" else NULL, type = if (qz) "quartz" else NULL)
  cat("图已写入:", fp2, "  ← UCell PCA 特征图\n")
}
