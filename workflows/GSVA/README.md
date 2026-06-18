# GSVA workflows

Learning and reproducibility workflows for the **GSVA** package (v2.x).

Based on the GSVA source under the local `GSVA` directory (Bioconductor
`GSVA`). References: Hänzelmann et al., 2013, *BMC Bioinformatics* (GSVA);
Barbie et al., 2009, *Nature* (ssGSEA).

> **范围说明 / Scope**：GSVA 包是四方法统一框架（`gsva` / `ssgsea` / `plage` /
> `zscore`）。本项目**只聚焦 `gsva` 与 `ssgsea` 两个 method**（通路/基因集分数）；
> `plage` 与 `zscore` 不在范围内。`ssGSEA` 即 GSVA 包内的一个 method，故无需单独
> 立项，远端 `ssGSEA-gpmodule` submodule 仅作历史 provenance。详见根目录 `CLAUDE.md` §1。

## Files

- `GSVA_学习手册.Rmd` — 类型 A 学习稿（中英对照标题 + 中文正文 + 代码注释）。
  源码级讲解 `gsva` 与 `ssgsea` 的全部参数、param-object API、调用调度流程，
  含端到端示例与速查表。可渲染 HTML / PDF（xelatex），可逐块测试。
- `GSVA_学习手册.R` — 同内容的纯脚本版，便于逐行 step。

## Chapters

| # | Topic |
|---|-------|
| 0 | 依赖与安装 |
| 1 | GSVA 与 ssGSEA 算法直觉（两种富集分数）+ 与 AUCell/UCell 的区别 |
| 2 | param-object API 范式（`gsva(param)`，新旧 API 区别）|
| 3 | 数据输入 + `kcdf` 与数据类型的匹配 |
| 4 | GSVA 方法：`gsvaParam` 逐参数（`kcdf`/`tau`/`maxDiff`/`absRanking`/`sparse`/`filterRows`）|
| 5 | ssGSEA 方法：`ssgseaParam` 逐参数（`alpha`/`normalize` + NA 策略）|
| 6 | 调度流程与返回值结构 |
| 7 | 端到端示例（gsva + ssgsea 同数据对比）+ 跨方法长表导出 + FAQ |
| 8 | 速查表 + gsva vs ssgsea 选择建议 |

## Render

```r
rmarkdown::render("workflows/GSVA/GSVA_学习手册.Rmd",
                  output_format = "html_document",
                  output_dir    = "results/GSVA")
```

Rendered outputs (`*.html`, `*.pdf`) are excluded by `.gitignore`.
