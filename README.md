# CausalityAI multiGRN Workflows

This repository collects reproducible workflow notes and scripts for learning and comparing three gene signature scoring methods used in this research project:

- UCell
- AUCell
- ssGSEA

The official source repositories are managed as Git submodules, so their upstream provenance and commit pins stay explicit.

## Clone

```bash
git clone --recurse-submodules <repo-url>
```

If the repository was cloned without submodules:

```bash
git submodule update --init --recursive
```

## Repository Layout

- `UCell/`: official UCell source repository as a submodule.
- `AUCell/`: official AUCell source repository as a submodule.
- `ssGSEA-gpmodule/`: official ssGSEA GenePattern module as a submodule.
- `examples/`: small runnable examples for package calls.
- `scripts/`: batch scripts for reproducible workflows.
- `styles/`: shared document styling.
- `*.Rmd`: bilingual workflow notebooks.

Generated outputs are intentionally ignored:

- `data/`
- `results/`
- rendered `*.html`, `*.pdf`, `*.tex`, and `*_files/`
- project-local R library `r-lib/`

## R Setup

This project uses a project-local R library:

```r
.libPaths(c("F:/WorkSpace/Causality+AI-multiGRN/r-lib", .libPaths()))
```

Install required packages into that library before running the workflows:

```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c(
  "UCell",
  "AUCell",
  "scRNAseq",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "BiocParallel",
  "Matrix",
  "ggplot2"
))
```

In this workspace, `UCell` and `AUCell` were installed from the local submodule sources.

## Run UCell Workflows

Small official UCell example:

```powershell
& 'F:\R-4.6.0\bin\Rscript.exe' examples\run_signature_scores.R
```

Official Zilionis dataset workflow:

```powershell
& 'F:\R-4.6.0\bin\Rscript.exe' scripts\run_ucell_zilionis_workflow.R
```

The first full-data run downloads `scRNAseq::ZilionisLungData()`, filters immune cells, keeps the first 5000 cells, and caches the processed object at:

```text
data/ucell_zilionis_immune_5000.rds
```

Subsequent runs reuse that cache.

## Notes

- Use `ScoreSignatures_UCell()` for matrix and SingleCellExperiment workflows.
- Use `StoreRankings_UCell()` when repeatedly testing many signatures on the same expression matrix.
- Seurat and SmoothKNN workflows are intentionally left for a later step to avoid heavy dependencies in this first repository version.
