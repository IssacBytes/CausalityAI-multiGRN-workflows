# UCell workflows

UCell 的学习与复现工作流。表达数据与签名**全部为官方真实数据**（`sample.matrix` /
`scRNAseq::ZilionisLungData()` + 官方 vignette 签名），符合 `CLAUDE.md` §3 数据铁律。

参考：Andreatta & Carmona, 2021, *CSBJ* (UCell)；官方 `vignette("UCell")`。

## Files

- `UCell_学习手册.Rmd` — **类型 A 学习稿**（唯一学习载体）。源码级讲解 UCell 全部参数、
  正负基因集、调用调度流程，含端到端示例与速查表。真实数据用包自带 `data(sample.matrix)`。
  导出 HTML / PDF（xelatex）供学习，逐块可测。
- `UCell_zilionis_workflow_zh_en.Rmd` — **类型 B 复现 workflow**（雏形）。用官方
  `scRNAseq::ZilionisLungData()`（约 5000 免疫细胞）走完整流程，中英对照。
- `run_ucell_zilionis_workflow.R` — 上述 workflow 的批处理脚本版。
- `run_sample_ucell_aucell_scores.R` — 用包自带 `sample.matrix` 的冒烟测试（UCell + AUCell）。

## 学习稿章节（Type A）

| # | Topic |
|---|-------|
| 0 | 依赖与安装 |
| 1 | UCell 是什么 + 算法直觉（Mann-Whitney U）+ 与 AUCell/GSVA 区别 |
| 2 | 准备官方数据（`sample.matrix` + 官方签名）|
| 3 | 基因签名格式 + 正负基因（`+/-`、`w_neg`）|
| 4 | `ScoreSignatures_UCell()` 逐参数（`maxRank`★ / `w_neg`★ / chunk.size / ...）|
| 5 | 预计算排名 `StoreRankings_UCell()` + `precalc.ranks` |
| 6 | Seurat / SCE 集成（`AddModuleScore_UCell()`）|
| 7 | kNN 平滑 `SmoothKNN()` |
| 8 | 端到端示例 + 跨方法长表 + 速查 + FAQ + 附录（U 统计量）|

## Run

```r
# 学习稿渲染（无需联网，sample.matrix 随包自带）
rmarkdown::render("workflows/UCell/UCell_学习手册.Rmd",
                  output_format = "html_document", output_dir = "results/UCell")

# 复现 workflow 批处理（首次会下载并缓存 Zilionis 数据到 data/）
# Rscript workflows/UCell/run_ucell_zilionis_workflow.R
```

> 路径已按 `CLAUDE.md` §5 改为可移植（动态 `project_root` + 环境变量 `MULTIGRN_RLIB`），
> 不再硬编码机器路径。渲染产物 / `data/` / `results/` 按 `.gitignore` 不入库。
