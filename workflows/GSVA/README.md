# GSVA workflows

Learning and reproducibility workflows for the **GSVA** package (v2.x).

Based on the GSVA source under the local `GSVA` directory (Bioconductor
`GSVA`). References: Hänzelmann et al., 2013, *BMC Bioinformatics* (GSVA);
Barbie et al., 2009, *Nature* (ssGSEA).

> **范围说明 / Scope**：GSVA 包是四方法统一框架（`gsva` / `ssgsea` / `plage` /
> `zscore`）。本项目**只聚焦 `gsva` 与 `ssgsea` 两个 method**（通路/基因集分数）；
> `plage` 与 `zscore` 不在范围内。`ssGSEA` 即 GSVA 包内的一个 method，故无需单独
> 立项，远端 `ssGSEA-gpmodule` submodule 仅作历史 provenance。详见根目录 `CLAUDE.md` §1。

**全程使用官方真实数据**（`GSVAdata`）——非模拟数据，遵守 `CLAUDE.md` §3 数据铁律。

## Files

- `GSVA_学习手册.Rmd` — 类型 A 学习稿，**唯一学习载体**（中英对照标题 + 中文正文 + 代码注释）。
  源码级讲解 `gsva` 与 `ssgsea` 的全部参数、param-object API、调用调度流程。
  真实数据：`commonPickrellHuang`（微阵列 + RNA-seq，演示 `kcdf`）、`gbm_VerhaakEtAl`
  （分子亚型复现）、`c2BroadSets`（MSigDB C2 通路）。**导出 HTML / PDF（xelatex）供学习**，逐块可测。
- `run_gsva_pickrellhuang.R` — **GSVA 运行脚本**（测试/复现用）。用 `gsvaParam` 在官方
  `commonPickrellHuang` 上跑，演示 `kcdf`（微阵列→Gaussian、RNA-seq 计数→Poisson）+
  跨平台一致性检查；结果写 `results/`（两平台宽表 + §7 长表）。离线可跑。
  运行：`Rscript workflows/GSVA/run_gsva_pickrellhuang.R`
- `run_ssgsea_pickrellhuang.R` — **ssGSEA 运行脚本**（测试/复现用）。ssGSEA 即 GSVA 包的
  `ssgsea` method（CLAUDE.md §1），用官方 `commonPickrellHuang` 微阵列 eset + C2 canonical
  通路打分，结果写 `results/`（宽表 + §7 统一长表）。可移植路径（§5），离线可跑。
  运行：`Rscript workflows/GSVA/run_ssgsea_pickrellhuang.R`

## Chapters

| # | Topic |
|---|-------|
| 0 | 依赖与安装（含 `GSVAdata` / `edgeR`）|
| 1 | GSVA 与 ssGSEA 算法直觉（两种富集分数）+ 与 AUCell/UCell 区别 |
| 2 | 准备官方数据（commonPickrellHuang / gbm / c2BroadSets）|
| 3 | param-object API 范式（`gsva(param)`）+ 基因集准备 |
| 4 | GSVA 方法 `gsvaParam` 逐参数（`kcdf` 在真实微阵列 vs RNA-seq 上演示 + tau/maxDiff/...）|
| 5 | ssGSEA 方法 `ssgseaParam` 逐参数（`alpha`/`normalize`）|
| 6 | 调度流程与返回值结构 |
| 7 | 官方真实示例（跨平台 kcdf 一致性 + GBM 分子亚型复现 + 跨方法长表导出）|
| 8 | 速查表 + gsva vs ssgsea 选择 + FAQ + 附录（KS 游走直觉）|

## Render

```r
rmarkdown::render("workflows/GSVA/GSVA_学习手册.Rmd",
                  output_format = "html_document",
                  output_dir    = "results/GSVA")
```

Rendered outputs (`*.html`, `*.pdf`) are excluded by `.gitignore`.
