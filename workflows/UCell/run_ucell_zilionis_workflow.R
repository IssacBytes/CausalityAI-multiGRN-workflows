#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args[startsWith(args, file_arg)]
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub(file_arg, "", script_arg[1]), winslash = "/"))
} else {
  getwd()
}
project_root <- normalizePath(file.path(script_dir, "../.."), winslash = "/")
setwd(project_root)

# 可移植 R 库路径（CLAUDE.md §5）：优先环境变量，回退项目内 r-lib/
local_lib <- Sys.getenv("MULTIGRN_RLIB", unset = file.path(project_root, "r-lib"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(UCell)
  library(BiocParallel)
  library(scRNAseq)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

set.seed(123)

BPPARAM <- SnowParam(workers = 1)
data_dir <- "data"
results_dir <- "results"
cache_file <- file.path(data_dir, "ucell_zilionis_immune_5000.rds")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

make_long_scores <- function(score_matrix) {
  data.frame(
    cell = rep(rownames(score_matrix), times = ncol(score_matrix)),
    signature = rep(colnames(score_matrix), each = nrow(score_matrix)),
    score = as.numeric(as.matrix(score_matrix)),
    row.names = NULL
  )
}

make_gene_overlap <- function(signatures, expr_mat) {
  rows <- lapply(names(signatures), function(signature_name) {
    genes <- signatures[[signature_name]]
    clean_genes <- sub("[+-]$", "", genes)
    data.frame(
      signature = signature_name,
      gene = genes,
      clean_gene = clean_genes,
      direction = ifelse(grepl("-$", genes), "negative", "positive"),
      present_in_matrix = clean_genes %in% rownames(expr_mat),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

load_or_create_zilionis_cache <- function(cache_file) {
  if (file.exists(cache_file)) {
    message("Reading cached data: ", cache_file)
    return(readRDS(cache_file))
  }

  message("Downloading ZilionisLungData() through scRNAseq...")
  lung <- scRNAseq::ZilionisLungData()
  immune <- lung$Used & lung$used_in_NSCLC_immune
  lung <- lung[, immune]
  lung <- lung[, seq_len(min(5000, ncol(lung)))]

  assay_name <- if ("counts" %in% SummarizedExperiment::assayNames(lung)) {
    "counts"
  } else {
    SummarizedExperiment::assayNames(lung)[1]
  }

  exp_mat <- Matrix::Matrix(
    SummarizedExperiment::assay(lung, assay_name),
    sparse = TRUE
  )
  colnames(exp_mat) <- paste0(colnames(exp_mat), seq_len(ncol(exp_mat)))

  cached <- list(
    exp.mat = exp_mat,
    colData = as.data.frame(SummarizedExperiment::colData(lung)),
    assay_name = assay_name,
    source = "scRNAseq::ZilionisLungData; immune cells; first 5000 cells"
  )
  saveRDS(cached, cache_file)
  message("Saved cache: ", cache_file)
  cached
}

cached <- load_or_create_zilionis_cache(cache_file)
exp.mat <- cached$exp.mat

stopifnot(inherits(exp.mat, "dgCMatrix"))
stopifnot(ncol(exp.mat) <= 5000)

signatures <- list(
  Tcell = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC"),
  Myeloid = c("CD14", "LYZ", "CSF1R", "FCER1G", "SPI1", "LCK-"),
  NK = c("KLRD1", "NCR1", "NKG7", "CD3D-", "CD3E-"),
  Plasma_cell = c("MZB1", "DERL3", "CD19-")
)

gene_overlap <- make_gene_overlap(signatures, exp.mat)
write.csv(
  gene_overlap,
  file.path(results_dir, "ucell_zilionis_gene_overlap.csv"),
  row.names = FALSE
)

message("Running ScoreSignatures_UCell()...")
ucell_scores <- ScoreSignatures_UCell(
  exp.mat,
  features = signatures,
  BPPARAM = BPPARAM
)

write.csv(ucell_scores, file.path(results_dir, "ucell_zilionis_scores.csv"))
saveRDS(ucell_scores, file.path(results_dir, "ucell_zilionis_scores.rds"))
write.csv(
  make_long_scores(ucell_scores),
  file.path(results_dir, "ucell_zilionis_scores_long.csv"),
  row.names = FALSE
)

message("Comparing maxRank values...")
scores_800 <- ScoreSignatures_UCell(
  exp.mat,
  features = signatures,
  maxRank = 800,
  BPPARAM = BPPARAM
)
scores_1500 <- ScoreSignatures_UCell(
  exp.mat,
  features = signatures,
  maxRank = 1500,
  BPPARAM = BPPARAM
)
scores_3000 <- ScoreSignatures_UCell(
  exp.mat,
  features = signatures,
  maxRank = 3000,
  BPPARAM = BPPARAM
)

compare_maxrank <- data.frame(
  cell = rownames(scores_1500),
  Tcell_800 = as.numeric(scores_800[, "Tcell_UCell"]),
  Tcell_1500 = as.numeric(scores_1500[, "Tcell_UCell"]),
  Tcell_3000 = as.numeric(scores_3000[, "Tcell_UCell"]),
  row.names = NULL
)
write.csv(
  compare_maxrank,
  file.path(results_dir, "ucell_zilionis_compare_maxRank.csv"),
  row.names = FALSE
)

message("Precomputing ranks and rescoring...")
ranks <- StoreRankings_UCell(
  exp.mat,
  maxRank = 1500,
  BPPARAM = BPPARAM
)
saveRDS(ranks, file.path(results_dir, "ucell_zilionis_ranks.rds"))

scores_from_ranks <- ScoreSignatures_UCell(
  features = signatures,
  precalc.ranks = ranks,
  BPPARAM = BPPARAM
)
write.csv(
  scores_from_ranks,
  file.path(results_dir, "ucell_zilionis_scores_from_ranks.csv")
)

message("Running SingleCellExperiment workflow...")
sce <- SingleCellExperiment(list(counts = exp.mat))
sce <- ScoreSignatures_UCell(
  sce,
  features = signatures,
  assay = "counts",
  BPPARAM = BPPARAM
)
saveRDS(sce, file.path(results_dir, "ucell_zilionis_sce_with_scores.rds"))

score_delta <- max(abs(as.matrix(ucell_scores) - as.matrix(scores_from_ranks)))
summary_table <- data.frame(
  metric = c(
    "n_genes",
    "n_cells",
    "n_signatures",
    "max_abs_delta_direct_vs_precalc_ranks",
    "sce_has_ucell_altExp"
  ),
  value = c(
    nrow(exp.mat),
    ncol(exp.mat),
    length(signatures),
    score_delta,
    "UCell" %in% SingleCellExperiment::altExpNames(sce)
  )
)
write.csv(
  summary_table,
  file.path(results_dir, "ucell_zilionis_workflow_summary.csv"),
  row.names = FALSE
)

print(summary_table)
message("Done. Results saved in: ", normalizePath(results_dir))
