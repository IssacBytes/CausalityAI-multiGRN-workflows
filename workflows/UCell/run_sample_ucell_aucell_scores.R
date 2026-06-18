## Minimal example for calling local UCell and AUCell packages.
## Run from the project root with:
## F:/R-4.6.0/bin/Rscript.exe workflows/UCell/run_sample_ucell_aucell_scores.R

.libPaths(c("F:/WorkSpace/Causality+AI-multiGRN/r-lib", .libPaths()))

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
