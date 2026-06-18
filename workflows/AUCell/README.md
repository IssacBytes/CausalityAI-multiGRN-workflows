# AUCell workflows

Learning and reproducibility workflows for AUCell (v1.25.2).

Based on the official AUCell source at the root-level `AUCell` submodule
(`aertslab/AUCell`). Reference: Aibar et al., 2017, *Nature Methods* (SCENIC).

## Files

- `AUCell_完整学习手册.Rmd` — bilingual (Chinese/code) R Markdown notebook.
  Covers the full AUCell workflow from source-code level:
  algorithm internals, all function parameters, threshold methods, and
  end-to-end reproducible example. Renders to HTML (floating TOC, code
  folding) or PDF via `xelatex`.
- `AUCell_完整学习手册.R` — plain R script version of the same material
  (1 892 lines, 10 chapters, dense inline comments). Useful for stepping
  through individual lines without rendering.

## Chapters

| # | Topic |
|---|-------|
| 0 | Package dependencies & installation |
| 1 | AUC algorithm internals (trapezoid integral, source-code walkthrough) |
| 2 | Supported input formats (matrix / dgCMatrix / SCE / ExpressionSet) |
| 3 | Gene-set operations (`nGenes`, `subsetGeneSets`, `setGeneSetNames`) |
| 4 | `AUCell_buildRankings()` — all 7 parameters |
| 5 | `AUCell_calcAUC()` — `aucMaxRank` deep-dive |
| 6 | `AUCell_run()` — block-processing one-liner |
| 7 | Threshold exploration (4 methods) + cell assignment + binarisation |
| 8 | `aucellResults` object operations (subset, cbind, orderAUC) |
| 9 | Visualisation (`plotHist`, `plotTSNE`, `plotEmb_rgb`) |
| 10 | Full end-to-end example + FAQ + quick-reference table |

## Render

```r
# HTML (recommended for interactive exploration)
rmarkdown::render("workflows/AUCell/AUCell_完整学习手册.Rmd",
                  output_format = "html_document",
                  output_dir    = "results/AUCell")

# PDF
rmarkdown::render("workflows/AUCell/AUCell_完整学习手册.Rmd",
                  output_format = "pdf_document",
                  output_dir    = "results/AUCell")
```

Rendered outputs (`*.html`, `*.pdf`) are excluded by `.gitignore`.
