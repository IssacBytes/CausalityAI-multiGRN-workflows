# UCell workflows

This folder contains the learning and reproducibility workflows for UCell.

## Files

- `UCell_zilionis_workflow_zh_en.Rmd`: bilingual R Markdown notebook for the official Zilionis data workflow.
- `run_ucell_zilionis_workflow.R`: batch script version of the same workflow.
- `run_sample_ucell_aucell_scores.R`: small smoke-test example using UCell's bundled `sample.matrix`.

## Run

From the project root:

```powershell
& 'F:\R-4.6.0\bin\Rscript.exe' workflows\UCell\run_ucell_zilionis_workflow.R
```

The first run downloads and caches the processed Zilionis data under `data/`.
Generated outputs are written under `results/`.
