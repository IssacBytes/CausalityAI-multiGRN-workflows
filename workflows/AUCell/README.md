# AUCell workflows

Learning workflow for AUCell (v1.25.2), grounded entirely in the **official
AUCell vignette + bundled data** — no fabricated/simulated matrices.

Reference: Aibar et al., 2017, *Nature Methods* (SCENIC); official vignette
`vignette("AUCell")`. Source tracked as the root-level `AUCell` submodule
(`aertslab/AUCell`).

## Files

- `AUCell_完整学习手册.Rmd` — 类型 A 学习稿，**唯一学习载体**。中英对照标题 +
  中文正文。源码级讲解全部函数参数与调用调度流程，**全程使用官方真实数据**：
  - 表达矩阵：GSE60361 小鼠脑（Zeisel et al., 2015，官方 vignette 下载）
  - 基因集：包自带 `inst/examples/geneSignatures.gmt`（脑细胞类型签名）
  - t-SNE / 标签：包自带 `cellsTsne.RData` / `mouseBrain_cellLabels.tsv`
  - 导出 HTML / PDF（xelatex）供学习；逐块可测（默认 `eval=FALSE`）。
- `run_aucell_mousebrain.R` — **AUCell 运行脚本**（测试/复现用）。官方 GSE60361 小鼠脑
  完整流程：下载（自动缓存到 `data/`）→ buildRankings → calcAUC → 阈值分配 →
  **混淆矩阵效度验证**（对照包自带真实细胞标签）。结果写 `results/`（AUC 宽表 + §7 长表
  + 混淆矩阵）。可移植路径（§5）。运行：`Rscript workflows/AUCell/run_aucell_mousebrain.R`
  （首次联网下载约几十 MB，之后读缓存）。

## Chapters

| # | Topic |
|---|-------|
| 0 | 依赖与安装 |
| 1 | AUCell 是什么 + 三步工作流（直觉层）|
| 2 | 准备官方数据（GSE60361 下载 + 缓存 + 子采样）|
| 3 | 基因集准备（官方 gmt + subsetGeneSets/nGenes/setGeneSetNames + 对照集）|
| 4 | `AUCell_buildRankings()` 逐参数（真实 plotStats 图）|
| 5 | `AUCell_calcAUC()` 逐参数（`aucMaxRank` + 4 种基因集输入格式）|
| 6 | 一步法 `AUCell_run()` |
| 7 | 阈值探索（真实 oligodendrocyte 双峰 vs 随机单峰；4 种方法；手动调整）|
| 8 | 结果对象 + 可视化（包自带真实 t-SNE 着色 + `AUCell_plotTSNE`）|
| 9 | 效度验证（细胞标签混淆矩阵）+ 为什么用 AUCell（vs 均值）|
| 10 | 速查表 + FAQ |
| 附录 A | AUC 积分推导（选读，超出学习必需深度）|

> 注：按项目规范（CLAUDE.md §3），学习稿只用 `.Rmd`；不再保留平行 `.R` 学习副本。

## Render

```r
rmarkdown::render("workflows/AUCell/AUCell_完整学习手册.Rmd",
                  output_format = "html_document", output_dir = "results/AUCell")
rmarkdown::render("workflows/AUCell/AUCell_完整学习手册.Rmd",
                  output_format = "pdf_document",  output_dir = "results/AUCell")
```

> 真正运行需联网下载 GSE60361（约几十 MB，首次后缓存）。渲染产物
> (`*.html`/`*.pdf`)、`data/`、`results/` 按 `.gitignore` 不入库。
