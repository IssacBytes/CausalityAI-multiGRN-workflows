# AGENTS.md

本仓库的**唯一权威协作规范是 [`CLAUDE.md`](CLAUDE.md)**。
所有 AI 协作者（Claude、Codex 等）开始任何工作前，请先完整阅读 `CLAUDE.md` 并严格遵守。

> The single source of truth for collaboration in this repo is **[`CLAUDE.md`](CLAUDE.md)**.
> Every AI collaborator (Claude, Codex, …) must read `CLAUDE.md` in full and follow it before doing any work.

## 快速要点 / Quick pointers（细则一律以 CLAUDE.md 为准）

- **分支 / Branch**：一台机器 = 一个方法 = 一个 `claude-<method>` 分支；通过 PR 合入 `main`，不直接 push `main`。（Codex 协作者请用 `codex-<method>` 命名以便区分。）
- **两类产出 / Two deliverables**（见 CLAUDE.md §3）：
  - 学习稿（Comprehension）—— 讲透参数 + 函数调度流程，**不深入内核**；`echo=TRUE`、模拟数据、重清楚。
  - 报告稿（Report）—— 美观科研 LaTeX、可视化丰富、真实数据；`echo=FALSE`、md→PDF；正文裸代码 ≤ ~10%，只留决策性参数调用。
- **路径可移植 / Portability**（§5）：禁止提交任何机器特定绝对路径（如 `F:/...`、`/Users/...`）。
- **submodule 纪律**（§6）：`AUCell`/`UCell`/`ssGSEA-gpmodule` 只读，绝不在其内部改动；我们的代码只放 `workflows/`。
- **不提交产物**：`*.html`/`*.pdf`/`data/`/`results/`/`r-lib/` 已被 `.gitignore` 忽略。
- **待决事项**（§8）：第三个方法是 GSVA 还是 ssGSEA、基准数据集与基因集，未澄清前不要擅自决定。
