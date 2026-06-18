# CLAUDE.md — 项目总指导（Project Charter）

> 本文件是本仓库所有 Claude 协作者（多台机器、多分支）的**唯一权威规范**。
> 任何新建文档、脚本、提交前，先读本文件。冲突时以本文件为准。
> 维护者：以「AUCell 总指导」分支起草；修改需在 PR 描述里说明理由。

---

## 1. 项目目标 / What the user needs

研究者正在搭建**多 GRN（基因调控网络）基因集打分方法**的学习与对比框架，最终用于自己的因果 + AI 分析。

三种核心打分方法（均为 R / Bioconductor / S4）：

| 方法 | 算法核心 | 典型场景 |
|------|---------|---------|
| **UCell** | Mann-Whitney U 统计量（基于排名）| scRNA-seq，大数据友好，支持 Seurat/SCE |
| **AUCell** | AUC 曲线下面积（基于排名）| scRNA-seq，逐细胞活性打分 |
| **GSVA**（含 ssGSEA）| KS 随机游走 / 富集分数 | bulk / 单细胞 / 微阵列；本项目只用 `gsva` + `ssgsea` 两个 method |

> ✅ **第三方法已定（2026-06）= GSVA 包**。经查本地源码 `GSVA/R/`，GSVA 包是**四方法统一框架**（`gsva` / `ssgsea` / `plage` / `zscore`），原文 `AllClasses.R:103`："GSVA implements four single-sample gene set analysis methods: PLAGE, combined z-scores, ssGSEA, and GSVA"。统一入口 `gsva(param)`，参数对象有 `gsvaParam` / `ssgseaParam` / `plageParam` / `zscoreParam`。
> 因此：
> - **ssGSEA 是 GSVA 包的一个 method**，无需单独立项；远端 `ssGSEA-gpmodule` submodule（GSEA-MSigDB 的 GenePattern 独立实现）降为**历史 provenance**，不再单独写 workflow。
> - **本项目只聚焦 `gsva` 与 `ssgsea` 两个 method**（用户只需通路/基因集分数）；`plage` 与 `zscore` **不在范围内**，学习稿中至多一句带过。

用户的真实诉求分两层：
1. **先学懂**：每个方法的源码调用方式、每个参数的含义与取值流程 —— 这样才能自己做分析。
2. **再对比/整合**：三种方法在同一数据、同一基因集上的打分结果可横向比较。

因此每个方法的交付物必须同时服务这两层（见 §3）。

---

## 2. 角色分工与协作 / How we work

- **一台机器 = 一个方法 = 一个分支**。分支命名：`claude-<method>`（小写），例：`claude-aucell`、`claude-ucell`、`claude-gsva`。
- 各分支独立推进自己负责的方法，**不要跨分支改别人的方法目录**，除非是为对齐本规范做的最小改动（并在 PR 说明）。
- `main` 是集成分支。每个方法成熟后通过 **PR 合入 main**，不直接 push main。
- 跨机器共享的唯一真相是远端仓库；本地差异（路径、R 库位置）一律通过 §5 的可移植写法吸收，**不得**把机器特定路径提交进版本库。

---

## 3. 产出内容规范 / Deliverables

每个方法有**两套面向不同目的的产出**，对应用户的两类真实需求。两者并存、互补。

### 类型 A —— 学习理解稿（Comprehension / Study）
**目的**：给用户自己 + 未来科研组成员上手用，为分析真实实验数据做准备。**学懂参数 + 学会函数调度流程**。

- **深度天花板（最关键的约束）**：
  - ✅ 讲透：每个参数的语义/类型/默认值/取值流程；函数之间如何串联调度（buildRankings→calcAUC→threshold 这类 **workflow orchestration**）；S4 不同输入类型走哪个分支；返回对象怎么取用。
  - ✅ 算法只讲到「**直觉 + 输入输出契约**」层面，让人知道某参数"在哪一步、起什么作用、调大调小有何后果"。
  - 🚫 **不往下走到内核**：不逐行复刻算法、不展开 C/Fortran/数值推导、不做"重写这个包"级别的源码考古。源码引用仅用于定位"参数在哪生效"，点到为止。
  - 判断准则：读者读完应能**正确调用与编排**，而非**重新实现**。
- **🔴 数据来源（铁律 / 全局项目点）**：学习稿的**实质示例必须使用官方源的真实数据**——即该包**官方 vignette 所用的真实数据集 + 包自带数据**，**严禁**用自造的模拟矩阵（`matrix(rnorm/rpois(...))` 之类）冒充。
  - 仅在演示**纯 API 机制**（如对象构造、参数签名）时可用模拟数据，且必须**照抄官方 toy 片段**并注明"官方示例数据，结果无意义"，不得自己另编一套。
  - 各方法的官方数据落点：
    - **AUCell** → GSE60361 小鼠脑（vignette 下载）+ 包自带 `geneSignatures.gmt` / `cellsTsne.RData` / `mouseBrain_cellLabels.tsv`
    - **GSVA** → `GSVAdata`：`commonPickrellHuang`（Huang 微阵列 + Pickrell RNA-seq 计数）+ `c2BroadSets` / `c7.immunesigdb`
    - **UCell** → `scRNAseq::ZilionisLungData()`（官方 vignette 数据）
  - 价值：唯有真实数据才能展示真实分布（如 AUCell 的双峰阈值）与**效度验证**（混淆矩阵对照已知细胞类型），模拟数据"output is meaningless"。
- **载体：仅 `.Rmd`**。学习稿就是这份 `.Rmd`——逐块可测，并**导出 PDF / HTML 供学习**（这是它的主要用途）。`echo=TRUE`（代码即教学内容）；默认 `eval=FALSE`（展示真实官方代码，不强制联网跑）。**重清楚，不重排版**。
- **不需要平行的 `.R` 学习副本**。`.R` 文件在本项目里一律是**测试/运行脚本**（如 Type B 的 `run_*.R`、冒烟测试），不是学习文档。学习只认 `.Rmd` → PDF/HTML。
- 结构：函数签名 → 逐参数 → 调用/分发流程（图或表）→ 返回值结构 → 最小可跑示例 → 速查表。
- **模板参照：`workflows/AUCell/AUCell_完整学习手册.Rmd`**（已按上述铁律重写：全程 GSE60361 官方数据、积分推导降为附录、无 `.R` 副本）。后续方法照此标杆。

### 类型 B —— 报告稿（Report / Publication）
**目的**：对外展示级产出——组会、论文附录、科研报告。**美观、科研风格、LaTeX 排版、可视化丰富**。

- **呈现标准（最关键的约束）**：
  - 真实公开数据上的完整分析 + 高质量图（ggplot2 统一主题）+ 结果表 + 文字解读结论。
  - LaTeX 科研排版：标题页、目录、图/表编号与题注、统一配色字体、参考文献。
  - **`echo=FALSE` 默认隐藏代码**，正文只呈现图、表、结论叙述（与学习稿相反）。计算细节如需展示，收进附录代码块。
  - md → **PDF（xelatex + 中文字体）为主**，HTML 备用。
- 形态：`.Rmd` 报告 + 配套 `run_*.R` 批处理（保证可复现，下载缓存→`data/`，结果→`results/`）。
- 模板参照：`workflows/UCell/UCell_zilionis_workflow_zh_en.Rmd`（现为雏形，需按 §4 报告排版标准升级到"美观"档）。

### 一句话区分
> **学习稿**＝把代码摊开教你怎么调（`echo=TRUE`、模拟数据、重清楚）；
> **报告稿**＝把代码藏起来给你看结果（`echo=FALSE`、真实数据、重美观）。

### 每方法目录结构（标准）
```
workflows/<Method>/
  README.md                              # 索引：文件清单 + 渲染/运行命令
  <Method>_学习手册.Rmd                  # 类型 A 学习稿（唯一学习载体，导出 PDF/HTML）
  <Method>_<dataset>_report_zh_en.Rmd    # 类型 B 报告稿（既有 *_workflow_* 命名予以保留）
  run_<method>_<dataset>.R               # 测试/运行脚本（Type B 批处理、冒烟测试等）
```

---

## 4. 风格规范 / Style

### 语言
- **章节标题**：中英对照，格式 `## 2. 计算 AUC / Calculate AUC`。
- **正文叙述**：中文为主，关键术语首次出现附英文（如「曲线下面积（AUC）」）。
- **代码注释**：中文为主；函数签名、参数名、报错信息保留英文原文。
- 面向「类型 B」复现文档时，关键结论句给出中英双句（参照 UCell 现有写法）。

### 文档结构
- 编号体系：类型 A 用「章（Chapter/Ch）」，类型 B 用「节（`## 1.` `## 2.`）」。
- 每个函数讲解固定顺序：**签名 → 逐参数 → 内部算法（源码行号）→ 返回值结构 → 示例**。
- 引用源码必须标注**文件名 + 行号**，例：`源码 02_calcAUC.R Line 275-282`。

### YAML（HTML + PDF 双输出，必须可同时渲染）
```yaml
output:
  html_document: { toc: true, toc_depth: 3, toc_float: true, code_folding: show, df_print: paged }
  pdf_document:  { toc: true, toc_depth: 3, latex_engine: xelatex }   # 中文必须 xelatex
```
- 共享样式放 `styles/`（已有 `styles/ucell_pdf.css`）；新方法复用，不各搞一套。
- 代码块全部具名（`{r calcAUC-demo}`），便于 RStudio chunk 导航与逐块测试。

### 代码约定
- 缩进 2 空格；管道按项目现有风格；变量 `snake_case`，函数沿用上游包命名。
- 每个参数单独一行注释说明（类型/默认/含义），重点参数标 `★`。
- 随机性必须 `set.seed()` 固定并注明原因。

### 报告稿专属（类型 B 才需要）

**代码占比 —— 原则：展示旋钮，隐藏接线（show the knobs, hide the wiring）**
- 全局 `echo=FALSE, message=FALSE, warning=FALSE`。
- **正文裸代码 ≤ 全文篇幅约 10%，且只能是"决策性代码"**：即定义方法本身的核心函数调用 + 确切参数值（1–5 行，放进「方法/参数」标注框）。例：`AUCell_calcAUC(geneSets, rankings, aucMaxRank = 1500, normAUC = TRUE)`。
- 其余全部隐藏：数据读取、reshape、绘图代码、阈值循环=「接线」，藏掉或收进末尾**附录/复现章节**；完整可跑代码交给配套 `run_*.R`。
- 版面由 **图 > 表 > 文字解读** 主导，代码是配角。这是上限不是目标——能压到只剩参数行更好。
- 理由：读者要的是"做了什么、发现了什么"，可复现性靠①正文写清参数 + ②附可跑脚本来保证，而非把代码塞进叙述。

**可视化 —— 按"回答什么问题"选图（哪种合适用哪种）**

| 图类型 | 回答的问题 | 何时用 | 实现 |
|--------|-----------|--------|------|
| 密度/直方/山脊图 | 单签名打分分布、是否双峰、阈值依据 | 阈值类方法必备；多签名并列用 ridgeline | `ggridges` |
| 降维特征图(UMAP/t-SNE，连续渐变) | 活性空间分布 | 单细胞头号必备 | ggplot+viridis |
| 二值/分配图(embedding) | 哪些细胞过阈值 | 做了细胞分配时 | ggplot 两色 |
| RGB 叠加(2–3 签名) | 多签名共定位 | 可选 | `plotEmb_rgb` |
| 小提琴/箱线(按细胞类型分组) | 效度验证：签名是否在预期型里高 | 必备验证图 | ggplot violin |
| 聚类热图(细胞×签名，标准化+注释条) | 全局结构、签名共现 | 签名多时 | `ComplexHeatmap` |
| 跨方法散点+Spearman/相关热图 | 三方法是否一致 | 对比报告 payoff 图 | ggplot+`ggpubr` |

- **两张工作马**：降维特征图 + 按细胞类型的小提琴图，几乎每份报告都该有。
- **审美红线**：连续量用 **viridis**（感知均匀、色盲友好），分类用固定定性色板；**禁用** 3D、双 Y 轴、饼图、jet/rainbow 色阶、无 alpha 的过度重叠散点（重叠改 `alpha` 或 `geom_hex`）。
- 统一 `ggplot2` + 共享主题 `theme_multigrn()`（放 `styles/`，定义字体/配色/留白）；每张图 `fig.cap=` 必填、自动编号。

**表格**：`knitr::kable()` + `kableExtra`（PDF 走 `booktabs`），不贴原始 `print()`。

**LaTeX 排版**：共享 preamble 放 `styles/`（如 `multigrn_preamble.tex`），经 YAML `includes: in_header:` 引入；统一中文字体（`xelatex` + `Noto Serif CJK` / `PingFang` / `SimSun` 按可用性回退）；启用 `toc`、图表编号、`booktabs`。

**图片质量**：`dpi=300`，矢量优先（PDF 图用 PDF/SVG 设备），`fig.align="center"`。

---

## 5. 路径与环境 / Portability（强制）

**禁止提交任何机器特定的绝对路径。** 已知风险：UCell 脚本里的 `F:/WorkSpace/...r-lib`、本机的 `/Users/jzyserver/...`。

统一写法：
```r
# project_root 从脚本/文档位置动态推导，绝不硬编码
# .Rmd 中：
knitr::opts_knit$set(root.dir = rprojroot::find_root(rprojroot::has_file(".gitmodules")))
# .R 脚本中（参照现有 run_ucell 写法）：从 --file= 推导 ../..

# R 库路径：优先环境变量，回退到项目内 r-lib/（已被 .gitignore 忽略）
local_lib <- Sys.getenv("MULTIGRN_RLIB", unset = file.path(getwd(), "r-lib"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# 并行后端 OS 自适应（Windows 用 Snow，Mac/Linux 用 Multicore）
BPPARAM <- if (.Platform$OS.type == "windows")
  BiocParallel::SnowParam(workers = 1) else BiocParallel::MulticoreParam(workers = 2)
```
- 所有产出文件路径相对 `project_root`；下载数据 → `data/`，结果 → `results/`（均已被 `.gitignore` 忽略）。

---

## 6. Git 与 submodule 纪律

- **三个上游包是 submodule**（`AUCell` / `UCell` / `ssGSEA-gpmodule`）：只读、只 pin 版本，**绝不在 submodule 内部改动或提交**。我们的代码一律放 `workflows/`。
- 提交粒度：一个逻辑改动一个 commit；信息用中文或英文均可，**结尾附** `Co-Authored-By: Claude <模型名> <noreply@anthropic.com>`。
- 渲染产物（`*.html`/`*.pdf`/`*.tex`/`*_files/`）、`data/`、`results/`、`r-lib/` 已在 `.gitignore`，**不提交**。
- 根目录草稿不入库；正式版只放 `workflows/<Method>/`（`.gitignore` 已忽略根级 `/AUCell_完整学习手册.*`）。
- 工作流：`claude-<method>` 分支 → push → 开 PR 合入 `main`。

---

## 7. 跨方法可比性 / Comparison schema（整合目标）

三种方法的打分结果必须能并到一张长表里横向比较。**统一输出 schema**（已在 `run_sample_ucell_aucell_scores.R` 雏形基础上扩展）：

```r
# 标准长表：每行 = 一个 (方法, 细胞, 基因集) 的打分
data.frame(
  method    = character(),  # "UCell" | "AUCell" | "GSVA" | ...
  cell      = character(),  # 细胞 ID
  signature = character(),  # 基因集名
  score     = numeric(),    # 连续打分
  binary    = integer(),    # 可选：阈值二值化 0/1（无则 NA）
  threshold = numeric()     # 可选：所用阈值（无则 NA）
)
```
- 每个方法的「类型 B」workflow 末尾导出该长表到 `results/<method>_scores_long.csv`。
- 对比用同一**公开数据集**与同一**基因集 .gmt**（待 §8 选定基准数据/基因集）。
- 评分尺度不同（UCell U 统计量 vs AUCell AUC vs GSVA 富集分数），对比时按方法内分布标准化后再比，且在文档中明确说明各分数的语义不可直接等同。

---

## 8. 当前状态与待办 / Status board

| 方法 | 类型A 学习稿 | 类型B 报告稿 | 分支 | 备注 |
|------|:---:|:---:|------|------|
| AUCell | ✅ 官方数据版 | ⬜ 待做 | `claude-aucell-rewrite`(PR#3) | 已按数据铁律重写为 GSE60361 真实数据；旧模拟版被替换 |
| UCell | ⬜ 待做 | 🟡 雏形 | `claude-ucell` | Zilionis 文档已在 main，需按 §4 升级到"美观报告"档 |
| GSVA | ⚠️ 需重写 | ⬜ | `claude-gsva`(已合并) | 现 main 上是**模拟数据版**，违反数据铁律；待用 `GSVAdata`(Pickrell/Huang + c2BroadSets) 重写 |

**优先级**：① GSVA 学习稿按数据铁律用 `GSVAdata` 重写（与 AUCell 对齐）；② 补 AUCell↔UCell 对称缺口（AUCell 加报告稿、UCell 加学习稿）；③ 各方法报告稿。

**待用户决策**（阻塞项）：
1. ~~第三个方法是 GSVA 还是 ssGSEA~~ → **已定 = GSVA 包，只做 `gsva`+`ssgsea` 两 method**（见 §1）。
2. 跨方法对比选用的**基准公开数据集**与**基准基因集 .gmt**（建议沿用 UCell 已用的 Zilionis，以便直接对齐）。
3. 仓库名不一致：远端 `CausalityAI-multiGRN-workflows` vs 本地路径 `Causality+AI+mulitiGRN` vs UCell 脚本里的 `Causality+AI-multiGRN`——仅影响硬编码路径，按 §5 改造后即无关，但请知悉。
