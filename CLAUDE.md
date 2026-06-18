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
| **GSVA / ssGSEA** | KS 随机游走 | bulk / 单细胞 / 微阵列 |

> ⚠️ **第三方法待澄清**：远端 submodule 与 README 把第三个方法定为 **ssGSEA**（`GSEA-MSigDB/ssGSEA-gpmodule`），但本地源码目录是 **GSVA**（`GSVA` 包），且根目录同时存在 GSVA 与 ssGSEA 两篇文献 PDF。两者算法相关但不同。**在用户明确前，GSVA 与 ssGSEA 视为各自独立的待办方法，不要混为一谈。** 见 §8。

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

每个方法在 `workflows/<Method>/` 下，**最终应包含以下两类互补文档**：

### 类型 A —— 学习手册（Learning Manual）
**目的**：源码级、逐参数的参考书。读者读完能独立调用每个函数。
- 文件名：`<Method>_学习手册.Rmd`（+ 同名 `.R` 纯脚本版，便于逐行 step）
- 内容必含：算法原理（含源码级公式/伪代码 + 源文件行号引用）、所有导出函数的**全部参数**（类型/默认值/含义/取值建议）、S4 方法分发差异、对象操作、可视化、端到端示例、函数速查表。
- 数据：可用模拟数据；默认 `eval=FALSE`（参考型，避免渲染时强依赖安装）。
- 模板参照：`workflows/AUCell/AUCell_完整学习手册.Rmd`（10 章结构）。

### 类型 B —— 可复现 Workflow（Reproducible Workflow）
**目的**：在**真实公开数据**上跑通完整流程，结果可复现、可对比。
- 文件名：`<Method>_<dataset>_workflow_zh_en.Rmd` + 配套 `run_<method>_<dataset>_workflow.R`（批处理版）
- 数据：公开数据集（如 UCell 用 `scRNAseq::ZilionisLungData()`）；下载缓存进 `data/`，结果写 `results/`。
- 默认 `eval=TRUE`，能 `Rscript` 一键跑完。
- 必含**参数实跑探索**（在真实数据上演示关键参数如何影响结果）。
- 模板参照：`workflows/UCell/UCell_zilionis_workflow_zh_en.Rmd`。

### 每方法目录结构（标准）
```
workflows/<Method>/
  README.md                              # 索引：文件清单 + 渲染/运行命令
  <Method>_学习手册.Rmd                  # 类型 A
  <Method>_学习手册.R                    # 类型 A 脚本版
  <Method>_<dataset>_workflow_zh_en.Rmd  # 类型 B
  run_<method>_<dataset>_workflow.R      # 类型 B 批处理
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

| 方法 | 类型A 学习手册 | 类型B 复现workflow | 分支 | 备注 |
|------|:---:|:---:|------|------|
| AUCell | ✅ 已完成 | ⬜ 待做 | `claude-aucell` | 手册 10 章已推送 |
| UCell | ⬜ 待做 | ✅ 已完成(Zilionis) | `claude-ucell` | 已在 main |
| GSVA / ssGSEA | ⬜ | ⬜ | `claude-gsva`? | **方法身份待用户澄清，见 §1** |

**优先级**：先补齐 AUCell↔UCell 的对称缺口（AUCell 加类型 B、UCell 加类型 A），再启动第三方法。

**待用户决策**（阻塞项）：
1. 第三个方法到底是 **GSVA** 还是 **ssGSEA**？还是两个都要？这决定 submodule 与目录结构。
2. 跨方法对比选用的**基准公开数据集**与**基准基因集 .gmt**（建议沿用 UCell 已用的 Zilionis，以便直接对齐）。
3. 仓库名不一致：远端 `CausalityAI-multiGRN-workflows` vs 本地路径 `Causality+AI+mulitiGRN` vs UCell 脚本里的 `Causality+AI-multiGRN`——仅影响硬编码路径，按 §5 改造后即无关，但请知悉。
