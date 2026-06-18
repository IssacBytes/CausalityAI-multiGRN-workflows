################################################################################
##                                                                            ##
##            AUCell 完整学习手册 —— 源码级参数解析与调用流程                    ##
##                                                                            ##
##  本文件以可运行 R 脚本形式，系统讲解 AUCell 所有函数的参数含义、              ##
##  调用流程、内部算法原理，以及不同输入类型的行为差异。                         ##
##                                                                            ##
##  参考源码路径：Causality+AI+mulitiGRN/AUCell/R/                            ##
##  参考文献：Aibar et al., 2017, Nature Methods (SCENIC)                     ##
##                                                                            ##
##  章节目录：                                                                 ##
##    Ch0  包依赖与安装                                                        ##
##    Ch1  核心算法原理（AUC 数学推导）                                         ##
##    Ch2  数据输入格式                                                        ##
##    Ch3  基因集操作函数                                                      ##
##    Ch4  Step1: AUCell_buildRankings()                                      ##
##    Ch5  Step2: AUCell_calcAUC()                                            ##
##    Ch6  一步法: AUCell_run()                                               ##
##    Ch7  Step3: 阈值探索与细胞分配                                            ##
##    Ch8  aucellResults 对象操作                                              ##
##    Ch9  可视化函数                                                          ##
##    Ch10 完整端到端示例                                                      ##
##                                                                            ##
################################################################################


#### ===== Chapter 0: 包依赖与安装 ===== ####
# 源文件：AUCell/DESCRIPTION, AUCell/NAMESPACE
#
# AUCell 的依赖分三类：
#   Imports   — 必须安装，内部调用
#   Suggests  — 可选，增强功能
#   Enhances  — 可选，并行加速

# ---------- 0.1 核心依赖（Imports）---------- #
# 以下包在 AUCell 安装时会自动安装
#
#   data.table       高效内存数据操作（排除 shift 函数以避免命名冲突）
#   DelayedArray     延迟计算框架，支持超大矩阵不全量读入内存
#   DelayedMatrixStats 对 DelayedArray 的矩阵统计（核心用于 colRanks）
#   GSEABase         GeneSet / GeneSetCollection S4 对象系统
#   Matrix           稀疏矩阵（dgCMatrix），单细胞数据的标准格式
#   methods          R 的 S4 类/泛型系统
#   mixtools         混合高斯分布拟合（用于阈值计算）
#   R.utils          工具函数（capitalize 等）
#   stats            cor/quantile/density/hclust 等统计函数
#   SummarizedExperiment  aucellResults 的父类容器
#   BiocGenerics     泛型函数（cbind/rbind 重载）
#   graphics / grDevices / utils  基础绘图与工具

# ---------- 0.2 可选增强（Suggests / Enhances）---------- #
# 并行计算（显著加速 buildRankings 和 calcAUC）
#   doMC / doRNG / doParallel / foreach / doSNOW
#
# 可视化增强
#   Rtsne      t-SNE 降维
#   dynamicTreeCut  层次聚类动态剪枝（orderAUC 需要）
#   R2HTML     生成 HTML 报告（AUCell_plotTSNE asPNG=TRUE 需要）
#   NMF        颜色调色板
#
# 数据格式支持
#   Biobase    ExpressionSet 对象支持
#   SingleCellExperiment  单细胞标准容器

# ---------- 0.3 安装代码 ---------- #
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 安装 AUCell 主包（含所有 Imports）
BiocManager::install("AUCell")

# 安装可选增强包
BiocManager::install(c(
  "doMC", "doRNG", "doParallel", "foreach",   # 并行
  "dynamicTreeCut",                            # orderAUC
  "Rtsne",                                     # t-SNE
  "SingleCellExperiment",                      # SCE 支持
  "Biobase"                                    # ExpressionSet 支持
))

# 加载
library(AUCell)
library(GSEABase)   # GeneSet/GeneSetCollection 操作
library(Matrix)     # 稀疏矩阵 dgCMatrix


#### ===== Chapter 1: 核心算法原理 ===== ####
# 源文件：AUCell/R/02_calcAUC.R（函数 .auc, .AUC.geneSet_norm, .AUC.geneSet_old）
#
# ── 1.1 方法论：为什么基于排序？ ──────────────────────────────────────────────
#
# AUCell 不直接用表达量，而是先把每个细胞内的基因按表达量从高到低排序。
# 这样做的好处：
#   1. 独立于归一化方式（UMI/TPM/FPKM 均可）
#   2. 独立于测序深度（每个细胞内部排序自成体系）
#   3. 对噪声鲁棒（单个基因的绝对值不重要，相对位置才重要）
#
# ── 1.2 AUC 的直觉理解 ────────────────────────────────────────────────────────
#
# 问题：在一个细胞中，目标基因集（如神经元标记基因）是否"集体排在前面"？
#
# 方法：只看排名前 aucMaxRank 位的基因，统计有多少属于目标基因集。
#   → 如果该基因集的基因大量出现在"前 5%"，说明它在这个细胞中高度活跃
#
# ── 1.3 AUC 的数学计算（梯形积分）─────────────────────────────────────────────
#
# 源码（.auc 函数，02_calcAUC.R Line 275-282）：
#
#   .auc <- function(oneRanking, aucThreshold, maxAUC) {
#     x <- unlist(oneRanking)           # 该细胞中，目标基因集各基因的排名
#     x <- sort(x[x < aucThreshold])   # 只保留排名 < aucMaxRank 的基因
#     y <- seq_along(x)                 # 累计计数 1, 2, 3, ...
#     sum(diff(c(x, aucThreshold)) * y) / maxAUC  # 梯形积分
#   }
#
# 几何意义：
#   - x 轴：基因排名（1 到 aucMaxRank）
#   - y 轴：到达该排名时，累计发现了多少目标基因
#   - 计算折线下方的面积（梯形法）
#   - 面积越大 → 目标基因越集中在排名前端 → AUC 越高
#
# 例子（5 个目标基因，aucMaxRank=20）：
#   假设目标基因排名为 [2, 5, 8, 12, 19]
#   x = c(2, 5, 8, 12, 19)
#   y = c(1, 2, 3,  4,  5)
#   面积 = (5-2)*1 + (8-5)*2 + (12-8)*3 + (19-12)*4 + (20-19)*5 = 3+6+12+28+5 = 54
#   maxAUC（归一化分母）= 见下方
#
# ── 1.4 归一化（normAUC=TRUE 的作用）─────────────────────────────────────────
#
# 源码（.AUC.geneSet_norm，Line 249-271）：
#
#   # 理论最大面积（所有目标基因都排在最前面时）：
#   x_th <- sort(1:nrow(gSetRanks))          # 假设目标基因排名 1,2,...,n
#   x_th <- x_th[x_th < aucThreshold]
#   y_th <- seq_along(x_th)
#   maxAUC <- sum(diff(c(x_th, aucThreshold)) * y_th)
#
# 归一化后，AUC 值 ∈ [0, 1]：
#   AUC = 1.0  → 所有目标基因都在排名前 aucMaxRank 位（最理想）
#   AUC = 0.0  → 没有任何目标基因出现在排名前 aucMaxRank 位
#   AUC = 0.5  → 中等活跃
#
# normAUC=FALSE 时（.AUC.geneSet_old）：
#   maxAUC = aucThreshold * nrow(gSetRanks)  # 简单相乘，不考虑基因集大小
#   → 不同大小的基因集 AUC 不可比，通常不用
#
# ── 1.5 与 GSEA/ssGSEA 的区别 ────────────────────────────────────────────────
#
# GSEA：
#   - 基于 Kolmogorov-Smirnov 统计量
#   - 考虑全排序，使用"游走"的最大偏差
#   - 设计用于样本间比较，需要排列检验
#
# ssGSEA（Single-sample GSEA）：
#   - 对 GSEA 的扩展，可对单个样本评分
#   - 使用差分面积统计量
#
# AUCell：
#   - 纯粹基于排名的梯形积分（ROC 曲线下面积的变体）
#   - 只看"前 X%"的基因（aucMaxRank），计算速度快
#   - 专为单细胞设计：每个细胞独立排序，不受批次效应影响
#   - 无需 p 值，直接输出 0-1 的活跃度分数


#### ===== Chapter 2: 数据输入格式 ===== ####
# 源文件：AUCell/R/01_buildRankings.R（S4 方法分发部分）
#
# AUCell 的核心函数（buildRankings、calcAUC、run）通过 S4 泛型实现多态，
# 自动识别不同的输入类型。支持的格式如下：
#
# ── 2.1 支持的表达矩阵类型 ────────────────────────────────────────────────────

## 格式一：普通 matrix（密集矩阵）
# 行 = 基因（需有行名），列 = 细胞（需有列名）
# 小数据集使用，内存占用与矩阵大小成正比
exprMatrix_dense <- matrix(
  data   = sample(0:10, 200 * 50, replace = TRUE),  # 模拟表达值（整数）
  nrow   = 200,                  # 200 个基因
  ncol   = 50,                   # 50 个细胞
  dimnames = list(
    paste0("Gene", 1:200),       # 行名 = 基因名（必须！）
    paste0("Cell", 1:50)         # 列名 = 细胞ID（必须！）
  )
)

## 格式二：稀疏矩阵 dgCMatrix（强烈推荐）
# 单细胞数据中大部分基因表达值为 0，稀疏格式节省 90%+ 内存
# buildRankings 检测到 dgCMatrix 后自动启用 splitByBlocks=TRUE
library(Matrix)
exprMatrix_sparse <- as(exprMatrix_dense, "dgCMatrix")
# 验证：
class(exprMatrix_sparse)   # "dgCMatrix"
object.size(exprMatrix_dense) / object.size(exprMatrix_sparse)  # 通常大几倍

## 格式三：ExpressionSet（Biobase 的老式容器）
# 内部自动调用 Biobase::exprs(exprMat) 提取矩阵
library(Biobase)
eset <- ExpressionSet(assayData = exprMatrix_dense)
# AUCell_buildRankings(eset) 等同于 AUCell_buildRankings(exprs(eset))

## 格式四：SummarizedExperiment / SingleCellExperiment（现代 Bioc 标准）
library(SummarizedExperiment)
sce <- SummarizedExperiment(
  assays = list(counts = exprMatrix_dense)
)
# 默认使用第一个 assay（通常是 counts）
# 若有多个 assay，用 assayName 参数指定：
#   AUCell_buildRankings(sce, assayName = "counts")
#   AUCell_run(sce, geneSets, assayName = "logcounts")

# ── 2.2 输入要求汇总 ─────────────────────────────────────────────────────────
#
#  要求                  说明
#  ──────────────────────────────────────────────────────────────────────
#  行 = 基因             rownames 必须是基因名（用于与基因集匹配）
#  列 = 细胞             colnames 必须是细胞 ID（用于结果矩阵的列名）
#  值 = 非负数           原始计数、UMI、TPM 均可，不要求预先归一化
#  行名唯一              重复行名会导致排序结果混乱
#  列名唯一              cbind 合并时会检查列名是否重复
#
# ── 2.3 实际场景：从 Seurat 对象提取 ────────────────────────────────────────
#
# # 方法 A：直接用 RNA assay 的稀疏矩阵
# exprMat <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
# # exprMat 是 dgCMatrix，可直接传给 AUCell_buildRankings
#
# # 方法 B：转为普通矩阵（小数据集时）
# exprMat <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))


#### ===== Chapter 3: 基因集操作函数 ===== ####
# 源文件：AUCell/R/aux_methods.R, aux_getSetNames.R
# 相关格式：GSEABase::GeneSet, GSEABase::GeneSetCollection
#
# AUCell 使用 GSEABase 的对象系统管理基因集。
# 理解两种对象的区别非常重要：
#   GeneSet         —— 单个基因集（一组基因 + 名称 + 元信息）
#   GeneSetCollection —— 多个 GeneSet 的容器

# ── 3.1 从 .gmt 文件读取基因集（最常用）─────────────────────────────────────
#
# .gmt 格式（Tab 分隔）：
#   基因集名称  描述/来源  Gene1  Gene2  Gene3  ...
#   （每行一个基因集，基因数不固定）

library(GSEABase)

# 包内自带示例 gmt 文件
gmtFile <- system.file("examples", "geneSignatures.gmt",
                       package = "AUCell")

# getGmt() 返回 GeneSetCollection
geneSets <- getGmt(gmtFile)
geneSets  # 打印摘要，显示有多少个基因集
#   GeneSetCollection
#     names: Oligodendrocyte_Cahoy, Neuron_Cahoy, ...
#     unique identifiers: ABCA8, ABHD2, ...
#     types in collection: NULL

# 查看基因集名称
names(geneSets)  # 字符向量：所有基因集的名字

# 查看某一个具体的 GeneSet
geneSets[["Oligodendrocyte_Cahoy"]]   # 单个 GeneSet 对象
GSEABase::geneIds(geneSets[["Oligodendrocyte_Cahoy"]])  # 提取基因名向量

# ── 3.2 手动创建基因集 ────────────────────────────────────────────────────────

# 方式 A：字符向量（单基因集，最简单）
geneSet_simple <- c("Gene1", "Gene5", "Gene12", "Gene30")
# 直接传给 calcAUC，会自动包装为 list(geneSet = ...)

# 方式 B：命名列表（多基因集）
geneSets_list <- list(
  NeuronMarkers   = c("SYP", "NEFL", "MAP2", "TUBB3"),
  AstroMarkers    = c("GFAP", "AQP4", "S100B", "VIM"),
  OligoMarkers    = c("MBP", "MOG", "PLP1", "CNP")
)

# 方式 C：GSEABase GeneSet 对象
geneSet_obj <- GeneSet(
  c("GFAP", "AQP4", "S100B"),
  setName    = "AstrocyteSignature",
  setIdentifier = "AS_001"
)

# 方式 D：GeneSetCollection 对象
gsc <- GeneSetCollection(list(
  GeneSet(c("SYP", "MAP2"), setName = "Neuron"),
  GeneSet(c("GFAP", "AQP4"), setName = "Astro")
))

# ── 3.3 nGenes()：查看基因集大小 ─────────────────────────────────────────────
# 源码：aux_methods.R，两个 S4 方法
#
# S4 方法分发：
#   nGenes,GeneSet-method          → length(geneIds(geneSet))
#   nGenes,GeneSetCollection-method → 命名数值向量

# 单个 GeneSet
nGenes(geneSet_obj)  # 整数标量：3

# GeneSetCollection
nGenes(geneSets)     # 命名向量：每个基因集的基因数
#   Oligodendrocyte_Cahoy  Neuron_Cahoy  ...
#            1238              1022      ...

# ── 3.4 subsetGeneSets()：过滤为表达矩阵中存在的基因 ─────────────────────────
# 源码：aux_methods.R，S4 方法 subsetGeneSets,GeneSetCollection-method
#
# 必须步骤！基因集中的基因 ID 必须与表达矩阵的行名一致。
# 如果不过滤，calcAUC 会警告"XX% 的基因不在数据中"。
#
# 内部实现：对每个基因集，取 geneIds(gs) ∩ geneNames（用 GSEABase & 运算符）

# 先构造一个含基因名的向量（通常是 rownames(exprMatrix)）
available_genes <- paste0("Gene", 1:200)  # 示例

# 过滤后基因集只保留 available_genes 中存在的基因
geneSets_filtered <- subsetGeneSets(geneSets, available_genes)

# 对比过滤前后的基因集大小
cbind(
  before = nGenes(geneSets),
  after  = nGenes(geneSets_filtered)
)

# ── 3.5 setGeneSetNames()：重命名基因集 ──────────────────────────────────────
# 源码：aux_methods.R，S4 方法 setGeneSetNames,GeneSetCollection-method
#
# 用途：在名称中附加基因数等元信息，方便后续阅图
# 约束：newNames 长度必须等于 geneSets 的数量

# 在名称后附加基因数
newNames <- paste0(names(geneSets_filtered),
                   " (", nGenes(geneSets_filtered), "g)")
geneSets_named <- setGeneSetNames(geneSets_filtered, newNames)
names(geneSets_named)
#   "Oligodendrocyte_Cahoy (1238g)"  "Neuron_Cahoy (1022g)"  ...

# ── 3.6 getSetNames()：用正则模式提取基因集名称 ──────────────────────────────
# 源码：aux_getSetNames.R
# 函数签名：getSetNames(aucMat, patterns, startChr="^", endChr=" |_")
#
# 用途：当有很多基因集时，快速筛选出感兴趣的子集
# 内部实现：grep(paste0(startChr, pattern, endChr), rownames(aucMat))

# 假设已经有了 cells_AUC（第 5 章产出）
# 找出名称中包含 "Oligodendrocyte" 的基因集
# gs_oligo <- getSetNames(cells_AUC, "Oligodendrocyte")

# 支持多个 pattern（向量化）
# gs_neuro_oligo <- getSetNames(cells_AUC, c("Neuron", "Oligodendrocyte"))

# ── 3.7 实用小贴士 ────────────────────────────────────────────────────────────
#
# 1. 基因 ID 匹配：基因集和表达矩阵必须用同一套 ID 系统
#    （如 Gene Symbol vs Ensembl ID vs Entrez ID）
#    不一致时需要先做 ID 转换（用 AnnotationDbi / biomaRt）
#
# 2. 基因集大小建议：
#    < 10 基因   → AUC 值不稳定，容易出现 AUC=0
#    10~100 基因 → 可用，但方差较大
#    100~2000 基因 → 理想范围，稳定可靠
#    > 2000 基因  → 趋向于覆盖管家基因，特异性下降
#
# 3. calcAUC 的自动过滤规则（源码 02_calcAUC.R Line 191-203）：
#    如果某个基因集中 ≥ 80% 的基因不在排序中
#    → 该基因集被自动排除，并输出 warning
#    如果所有基因集都 ≥ 80% 缺失 → 直接报错，提示检查 ID 匹配


#### ===== Chapter 4: Step 1 — AUCell_buildRankings() ===== ####
# 源文件：AUCell/R/01_buildRankings.R
#
# 目的：为每个细胞建立基因表达排序（从高到低）。
# 这是 AUCell 流程的第一步，产出物是后续所有计算的基础。
#
# ── 4.1 函数签名（完整） ──────────────────────────────────────────────────────
#
#   AUCell_buildRankings(
#     exprMat,                   # 表达矩阵（见 Ch2）
#     featureType   = "genes",   # 行的类型名称（仅用于命名，不影响计算）
#     plotStats     = TRUE,      # 是否绘制每细胞基因检测数分布图
#     splitByBlocks = FALSE,     # 是否分块并行处理
#     BPPARAM       = NULL,      # BiocParallel 并行参数
#     keepZeroesAsNA = FALSE,    # 零值处理方式
#     verbose       = TRUE,      # 显示进度消息
#     nCores        = NULL,      # 已废弃，改用 BPPARAM
#     mctype        = NULL,      # 已废弃
#     ...                        # 透传给内部函数
#   )
#
# ── 4.2 S4 方法分发（输入类型自动识别）──────────────────────────────────────
#
# 源码在 01_buildRankings.R 中定义了 4 个 S4 方法：
#
#   exprMat 类型              内部行为
#   ──────────────────────────────────────────────────────────────
#   matrix                   默认方法，splitByBlocks 遵从参数设置
#   dgCMatrix（稀疏矩阵）     强制 splitByBlocks=TRUE，启用块处理
#   ExpressionSet            自动调用 Biobase::exprs(exprMat) 提取矩阵
#   SummarizedExperiment /   自动调用 assay(exprMat) 或 assays(exprMat)[[assayName]]
#   SingleCellExperiment

# ── 4.3 参数详解 ──────────────────────────────────────────────────────────────

## 参数 1: exprMat（必需）
# 类型：matrix / dgCMatrix / ExpressionSet / SummarizedExperiment
# 格式：行 = 基因（有行名），列 = 细胞（有列名），值 = 表达量（非负）
# 推荐：优先使用 dgCMatrix（单细胞数据标准稀疏格式）

## 参数 2: featureType = "genes"
# 类型：character（长度 1）
# 作用：给排序对象命名，如"genes"/"peaks"/"regions"
# 重要性：★☆☆ 仅用于提示信息和结果对象的描述，不影响任何计算
# 场景：做 ATAC-seq 的 peak 排序时可改为 featureType="peaks"
featureType_example <- "genes"  # 保持默认即可

## 参数 3: plotStats = TRUE
# 类型：logical
# 作用：调用内部 plotGeneCount() 函数，绘制"每个细胞中检测到的基因数"分布图
# 重要性：★★★ 必看！用于确定合适的 aucMaxRank
#
# 输出内容：
#   - 箱线图（上）：基因检测数的整体分布
#   - 直方图（下）：各细胞的基因检测数频率分布
#   - 控制台打印 6 个分位数：min / 1% / 5% / 10% / 50% / 100%
#
# 如何利用：
#   看这个图上的"中位数"或"25%分位数"
#   aucMaxRank 必须小于这个值，否则大量细胞的排名数据不足
#   例：若中位数是 3000 个基因，则 aucMaxRank = 3000 * 5% = 150 是安全的

## 参数 4: splitByBlocks = FALSE（dgCMatrix 自动=TRUE）
# 类型：logical
# 作用：是否将矩阵分块后再处理（可与 BPPARAM 配合实现并行）
# 内部实现：调用 DelayedArray::blockApply()，将矩阵按列分块
#
# 关键规则（源码 Line ~85）：
#   - 输入为 dgCMatrix 时：自动设置为 TRUE，无论用户传入什么值
#   - 输入为 matrix 时：遵从用户设置（默认 FALSE）
#   - splitByBlocks=TRUE + BPPARAM 多核 = 最大并行加速
#
# 何时手动设置 TRUE：
#   - 输入是普通 matrix 但数据量很大，想手动并行
#   - 内存受限时，分块避免一次性读入

## 参数 5: BPPARAM = NULL
# 类型：BiocParallelParam 对象 或 NULL
# 作用：控制并行计算方式（仅在 splitByBlocks=TRUE 时生效）
# NULL → 单核串行（默认）
#
# 可用的 BPPARAM 类型（需要 BiocParallel 包）：
library(BiocParallel)
#
# 多核（推荐 Linux/Mac）
bpparam_mc  <- MulticoreParam(workers = 4)   # 4 个核
#
# 多进程（Windows 兼容）
bpparam_sp  <- SnowParam(workers = 4)        # 4 个进程
#
# 串行（调试用）
bpparam_seq <- SerialParam()
#
# 使用示例（大数据集）：
# cells_rankings <- AUCell_buildRankings(
#   exprMatrix_sparse,
#   plotStats     = TRUE,
#   BPPARAM       = MulticoreParam(4)
# )

## 参数 6: keepZeroesAsNA = FALSE
# 类型：logical
# 作用：控制零表达量基因在排序中的处理方式
#
# FALSE（默认）：
#   零值基因被随机打乱，分布在排序的末尾区域
#   → 体现"这些基因确实没有表达，但排名顺序未知"
#   → 对 AUC 计算结果有轻微随机性（每次运行可能微小差异）
#
# TRUE（严格模式）：
#   零值基因的排名设为 NA，不参与 AUC 计算
#   → 更保守：只统计真正有表达的基因
#   → 适合：数据质量高、零值可信的情况
#
# 建议：大多数情况保持 FALSE（默认），对结果影响很小

## 参数 7: verbose = TRUE
# 类型：logical
# 作用：控制是否在控制台显示进度消息（正在处理第几块等）
# FALSE：静默运行，适合批处理脚本

# ── 4.4 内部算法（源码级） ────────────────────────────────────────────────────
#
# 核心计算在 .AUCell_buildRankings() 内部函数中：
#
#   Step 1: 将表达矩阵取负号（-exprMat）
#           → 使 colRanks() 的排序从"高到低"变成"低到高"（rank 1 = 最高表达）
#
#   Step 2: 调用 DelayedMatrixStats::colRanks(-exprMat, ties.method="random")
#           → 对每一列（细胞）独立计算各行（基因）的排名
#           → ties.method="random"：相同值的基因排名随机打散（避免系统性偏差）
#
#   Step 3: 将排名矩阵封装为 aucellResults 对象
#           → 内部存储在 SummarizedExperiment 的 assays[["ranking"]] 中

# ── 4.5 返回值结构 ────────────────────────────────────────────────────────────
#
# 返回：aucellResults 对象（S4 类，继承自 SummarizedExperiment）
#
# 内容：
#   assayNames(result) → "ranking"（只有一个 assay）
#   dim(result)        → c(nGenes, nCells)（基因数 × 细胞数）
#   rownames(result)   → 基因名（与输入 exprMat 相同）
#   colnames(result)   → 细胞 ID（与输入 exprMat 相同）
#
# 访问排名矩阵：
#   getRanking(result)  ← 推荐，有类型检查
#   assay(result, "ranking")  ← 直接访问（SummarizedExperiment 方式）
#
# 额外槽：
#   result@nGenesDetected → 每个细胞检测到的基因数（分位数向量）

# ── 4.6 完整使用示例 ──────────────────────────────────────────────────────────

# 构建示例数据
set.seed(42)
nGenes <- 5000
nCells <- 300
exprMat <- matrix(
  rpois(nGenes * nCells, lambda = 0.3),   # 模拟单细胞稀疏数据（泊松分布）
  nrow = nGenes, ncol = nCells,
  dimnames = list(paste0("Gene", 1:nGenes),
                  paste0("Cell", 1:nCells))
)
exprMat_sparse <- as(exprMat, "dgCMatrix")

## 基础用法（会自动绘制 plotStats 图）
cells_rankings <- AUCell_buildRankings(
  exprMat   = exprMat_sparse,  # 推荐稀疏矩阵
  plotStats = TRUE,            # 必看！观察基因检测数分布
  verbose   = TRUE
)

## 查看结果
cells_rankings           # 打印 aucellResults 摘要
class(cells_rankings)    # "aucellResults"
dim(cells_rankings)      # c(5000, 300)
assayNames(cells_rankings)  # "ranking"

## 查看排名矩阵的一小块（前3基因，前4细胞）
getRanking(cells_rankings)[1:3, 1:4]
#   Cell1  Cell2  Cell3  Cell4
# Gene1  2341   1823   3011   456    （排名值，越小越靠前）
# Gene2   789   4231    123  2901
# Gene3  3412    567   2134   789

## 生产场景（并行 + 静默）
# cells_rankings <- AUCell_buildRankings(
#   exprMat      = exprMat_sparse,
#   plotStats    = FALSE,          # 批处理时关闭绘图
#   BPPARAM      = MulticoreParam(8),
#   verbose      = FALSE
# )


#### ===== Chapter 5: Step 2 — AUCell_calcAUC() ===== ####
# 源文件：AUCell/R/02_calcAUC.R
#
# 目的：对每个细胞、每个基因集，计算 AUC 活跃度分数（0-1）。
# 输入：Step1 的排序（cells_rankings）+ 基因集列表
# 输出：AUC 矩阵（基因集 × 细胞）
#
# ── 5.1 函数签名（完整） ──────────────────────────────────────────────────────
#
#   AUCell_calcAUC(
#     geneSets,                              # 基因集（多种格式，见下）
#     rankings,                              # buildRankings 的输出
#     nCores     = 1,                        # 计算核数
#     normAUC    = TRUE,                     # 是否归一化 AUC 最大值为 1
#     aucMaxRank = ceiling(0.05*nrow(rankings)),  # ★ 最关键参数
#     verbose    = TRUE
#   )
#
# ── 5.2 参数详解 ──────────────────────────────────────────────────────────────

## 参数 1: geneSets（必需）
# 支持 4 种类型，S4 泛型自动分发：
#
#  输入类型              内部转换                     适用场景
#  ──────────────────────────────────────────────────────────────────
#  character 向量        list(geneSet = 向量)         单个基因集，快速测试
#  named list            直接使用                     最常用，多基因集
#  GeneSet               setNames(list(geneIds), setName)  GSEABase 单对象
#  GeneSetCollection     GSEABase::geneIds(geneSets)  GSEABase 集合对象

# 示例：4 种格式等价写法
gs_char <- c("Gene1", "Gene5", "Gene20")   # 字符向量（单基因集）

gs_list <- list(                             # 命名列表（多基因集）
  Set_A = c("Gene1", "Gene5", "Gene20"),
  Set_B = c("Gene3", "Gene8", "Gene15", "Gene22")
)

gs_GeneSet <- GeneSet(
  c("Gene1", "Gene5", "Gene20"),
  setName = "Set_A"
)

gs_GeneSetCollection <- GeneSetCollection(list(
  GeneSet(c("Gene1", "Gene5", "Gene20"), setName = "Set_A"),
  GeneSet(c("Gene3", "Gene8", "Gene15"), setName = "Set_B")
))

## 参数 2: rankings（必需）
# 类型：aucellResults 对象，内含 "ranking" assay
# 必须是 AUCell_buildRankings() 的输出
# 内部验证（源码 Line 113-123）：
#   - 检查 class 是否为 "aucellResults"
#   - 检查 assayNames 是否为 "ranking"
#   - 若是旧版 matrixWrapper → 报错提示用 updateAucellResults() 升级
#   - 验证通过后，调用 getRanking() 提取裸矩阵进行计算

## ★ 参数 3: aucMaxRank = ceiling(0.05 * nrow(rankings))  ← 最重要参数！
# 类型：integer（正整数）
# 默认值：基因总数的 5%（如 5000 基因 → aucMaxRank = 250）
#
# 含义：只看排名前 aucMaxRank 位的基因，在这个"窗口"里统计目标基因的密度
#
# 几何意义：
#   AUC 积分的上限（x 轴截止点）
#   aucMaxRank 越大 → 看更多基因 → 对低表达基因更敏感，但噪声也更多
#   aucMaxRank 越小 → 只看高表达基因 → 更严格，假阳性更少，但可能漏掉部分信号
#
# 选择建议（三步法）：
#   Step A: 先运行 buildRankings(plotStats=TRUE)，看基因检测数分布图
#   Step B: 找到"大多数细胞"的基因检测数（例如 50% 分位数 = 2000 基因）
#   Step C: aucMaxRank 应该 << 2000（如 5% × 2000 = 100）
#
# 典型值范围：
#   1% - 5%  → 保守，适合高质量、高深度数据
#   5% - 10% → 推荐默认值，平衡灵敏度与特异性
#   10%-20%  → 宽松，适合基因检测数高的 bulk-like 数据
#
# 警告机制（源码 Line 108-111）：
#   aucMaxRank < 300 → 输出 warning（提醒基因数可能太少）
#   aucMaxRank <= 0  → 直接 stop（错误）
#
# 实际设置示例：
nGenes_total <- nrow(cells_rankings)   # 总基因数 5000
aucMaxRank_5pct  <- ceiling(0.05 * nGenes_total)   # 250（推荐默认）
aucMaxRank_10pct <- ceiling(0.10 * nGenes_total)   # 500（宽松）
aucMaxRank_1pct  <- ceiling(0.01 * nGenes_total)   # 50（严格）

## 参数 4: normAUC = TRUE
# 类型：logical
# 作用：是否将 AUC 最大值归一化为 1
#
# TRUE（默认，.AUC.geneSet_norm）：
#   maxAUC 计算方式：假设所有目标基因都排在最前面时的理论最大面积
#   → AUC 值 ∈ [0, 1]，不同大小的基因集可直接比较
#   → 推荐！
#
# FALSE（.AUC.geneSet_old）：
#   maxAUC = aucThreshold × n目标基因（简单乘积）
#   → AUC 值不归一化，较大的基因集倾向于有更高的绝对值
#   → 不推荐用于基因集间比较

## 参数 5: nCores = 1
# 类型：integer
# 注意：源码中多核支持（doMC/doSNOW）已保留但不完全稳定
# 更推荐用 AUCell_run() 的 BPPARAM 参数实现并行

## 参数 6: verbose = TRUE
# 同 buildRankings，控制消息输出

# ── 5.3 内部计算流程（源码级） ────────────────────────────────────────────────
#
# .AUCell_calcAUC() 内部执行顺序（02_calcAUC.R）：
#
#   1. 验证输入
#      - geneSets 必须是 named list
#      - 移除空基因集（length <= 0）
#      - 验证 rankings 是合法的 aucellResults
#      - 提取 getRanking(rankings) 为裸矩阵
#
#   2. 选择 AUC 计算函数
#      - normAUC=TRUE  → .AUC.geneSet_norm（归一化版本）
#      - normAUC=FALSE → .AUC.geneSet_old（非归一化版本）
#
#   3. 对每个基因集调用 .AUC.geneSet_norm(geneSet, rankings, aucMaxRank, name)
#
#      内部逻辑（.AUC.geneSet_norm，Line 249-271）：
#        a. geneSet <- unique(geneSet)            # 去重
#        b. geneSet <- geneSet[geneSet %in% rownames(rankings)]  # 过滤不在数据中的基因
#        c. gSetRanks <- rankings[geneSet, ]      # 提取目标基因的排名子矩阵
#        d. 计算 maxAUC（理论上界）
#        e. apply(gSetRanks, 2, .auc, aucMaxRank, maxAUC)  # 对每个细胞计算 AUC
#        f. 返回 AUC 向量 + missing 计数 + nGenes
#
#      .auc() 函数（Line 275-282）—— AUC 的核心：
#        x <- sort(oneRanking[oneRanking < aucMaxRank])  # 目标基因中排名 < 阈值的
#        y <- seq_along(x)                               # 累计计数
#        AUC = sum(diff(c(x, aucMaxRank)) * y) / maxAUC  # 梯形积分
#
#   4. 过滤缺失基因超过 80% 的基因集（warning）
#   5. 移除 missing/nGenes 辅助列
#   6. 封装为 aucellResults（assay = "AUC"）

# ── 5.4 返回值结构 ────────────────────────────────────────────────────────────
#
# 返回：aucellResults 对象
#   assayNames(result) → "AUC"
#   dim(result)        → c(nGeneSets, nCells)
#   rownames(result)   → 基因集名称
#   colnames(result)   → 细胞 ID
#   值域：[0, 1]（normAUC=TRUE 时）
#
# 访问 AUC 矩阵：
#   getAUC(cells_AUC)  ← 推荐

# ── 5.5 完整使用示例 ──────────────────────────────────────────────────────────

# 准备示例基因集（用表达矩阵中存在的基因名）
set.seed(1)
geneSets_demo <- list(
  SignatureA = sample(paste0("Gene", 1:5000), 100),   # 随机 100 基因
  SignatureB = sample(paste0("Gene", 1:5000), 50),    # 随机 50 基因
  SignatureC = paste0("Gene", 1:30)                    # 前 30 个基因
)

## 基础调用（默认参数）
cells_AUC <- AUCell_calcAUC(
  geneSets   = geneSets_demo,
  rankings   = cells_rankings,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),  # 明确指定
  normAUC    = TRUE,
  verbose    = TRUE
)

## 查看结果
cells_AUC
dim(cells_AUC)       # c(3, 300)  → 3个基因集 × 300个细胞
assayNames(cells_AUC) # "AUC"

## 提取 AUC 矩阵
auc_matrix <- getAUC(cells_AUC)
class(auc_matrix)    # "matrix"
auc_matrix[, 1:5]   # 前 5 个细胞的 AUC 值

## 查看 AUC 值的分布
summary(auc_matrix["SignatureA", ])
# Min. 1st Qu. Median  Mean 3rd Qu.  Max.
# 0.00   0.03   0.05  0.05   0.07   0.15

## 不同 aucMaxRank 的对比（重要！）
cells_AUC_5pct  <- AUCell_calcAUC(geneSets_demo, cells_rankings,
                                   aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
cells_AUC_10pct <- AUCell_calcAUC(geneSets_demo, cells_rankings,
                                   aucMaxRank = ceiling(0.10 * nrow(cells_rankings)))
# 对比两种设置下 AUC 值的分布差异
# 较大的 aucMaxRank → AUC 均值通常更高（窗口更宽，更容易捕捉到目标基因）


#### ===== Chapter 6: 一步法 — AUCell_run() ===== ####
# 源文件：AUCell/R/aux_AUCell_run_byBlocks.R
#
# 等价于：AUCell_buildRankings() + AUCell_calcAUC() 的组合
# 区别：run() 使用 DelayedArray::blockApply() 分块串联，
#        不保存中间排序结果，直接输出 AUC
#
# ── 6.1 函数签名（完整） ──────────────────────────────────────────────────────
#
#   AUCell_run(
#     exprMat,                              # 表达矩阵（同 buildRankings）
#     geneSets,                             # 基因集（同 calcAUC）
#     featureType    = 'genes',             # 行类型名称（同 buildRankings）
#     keepZeroesAsNA = FALSE,               # 零值处理（同 buildRankings）
#     normAUC        = TRUE,                # AUC 归一化（同 calcAUC）
#     aucMaxRank     = ceiling(0.05*nrow(exprMat)),  # 关键参数（同 calcAUC）
#     BPPARAM        = NULL,                # 并行参数
#     assayName      = NULL,                # 仅 SE 对象需要指定 assay 名
#     ...
#   )
#
# ── 6.2 内部实现（源码级）────────────────────────────────────────────────────
#
# .AUCell_run()（块调度器，Line 120-138）：
#   1. 将 exprMat 包装为 DelayedArray
#   2. 使用 colAutoGrid(exprMat) 自动划分列方向的分块网格
#   3. blockApply() 对每个块调用 .AUCell_run_internal()
#   4. 所有块的 AUC 结果通过 cbind() 拼回
#
# .AUCell_run_internal()（单块处理，Line 141-148）：
#   AUCell_buildRankings(block, splitByBlocks=FALSE, plotStats=FALSE, verbose=FALSE)
#   ↓
#   AUCell_calcAUC(geneSets, ranked, normAUC, aucMaxRank, verbose=FALSE)
#
# 关键区别：
#   buildRankings 分步 → 可以保留 rankings 对象，复用于多个基因集
#   run() 一步         → 不保存 rankings，基因集固定，内存更节省
#
# ── 6.3 S4 方法分发 ───────────────────────────────────────────────────────────
#
#   exprMat 类型                内部处理
#   ──────────────────────────────────────────────────────────────
#   matrix                     直接 blockApply 分块
#   dgCMatrix                  转为 DelayedArray 后 blockApply
#   SummarizedExperiment       先提取 assay（assayName 指定）再 blockApply

# ── 6.4 参数特有说明 ──────────────────────────────────────────────────────────

## 参数: assayName = NULL（仅 SE 对象）
# 当 exprMat 是 SummarizedExperiment 或 SingleCellExperiment 时：
#   NULL → 使用第一个 assay（有多个时会 warning）
#   "counts" / "logcounts" 等 → 明确指定 assay 名称
#
# 内部验证（源码 Line 100-112）：
#   is.null(assayName) + 多 assay → warning，使用第一个
#   assayName 不存在 → stop

# ── 6.5 AUCell_run() vs 分步调用的选择 ──────────────────────────────────────
#
#   场景                          推荐方式
#   ────────────────────────────────────────────────────────────
#   需要探索不同基因集             buildRankings → calcAUC × N
#   （rankings 重复使用）          （保存 rankings 一次，多次 calcAUC）
#
#   基因集固定，追求内存效率        AUCell_run()
#   （不需要中间 rankings）
#
#   超大矩阵（>10万细胞）           AUCell_run() + BPPARAM 并行
#   内存有限
#
#   需要检查每细胞基因数分布        buildRankings(plotStats=TRUE) 分步
#   （确认 aucMaxRank 合理性）

# ── 6.6 完整使用示例 ──────────────────────────────────────────────────────────

## 等价用法 A：分步（推荐用于探索阶段）
cells_rankings_A <- AUCell_buildRankings(
  exprMat_sparse, plotStats = TRUE, verbose = TRUE
)
cells_AUC_A <- AUCell_calcAUC(
  geneSets_demo, cells_rankings_A,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings_A))
)

## 等价用法 B：一步法（推荐用于生产流程）
cells_AUC_B <- AUCell_run(
  exprMat    = exprMat_sparse,
  geneSets   = geneSets_demo,
  aucMaxRank = ceiling(0.05 * nrow(exprMat_sparse)),  # 注意：基于 exprMat 行数
  normAUC    = TRUE,
  verbose    = FALSE
)

## 验证两种方式结果一致（精度差异来自随机打平）
all.equal(
  getAUC(cells_AUC_A),
  getAUC(cells_AUC_B),
  tolerance = 0.01   # 因 ties.method="random"，允许小误差
)

## 并行一步法（大数据集）
# cells_AUC_parallel <- AUCell_run(
#   exprMat    = exprMat_sparse,
#   geneSets   = geneSets_demo,
#   aucMaxRank = 500,
#   BPPARAM    = MulticoreParam(8)
# )

## SummarizedExperiment 输入
# sce <- SingleCellExperiment(assays = list(counts = exprMat_sparse))
# cells_AUC_sce <- AUCell_run(
#   exprMat   = sce,
#   geneSets  = geneSets_demo,
#   assayName = "counts"          # 指定用哪个 assay
# )


#### ===== Chapter 7: Step 3 — 阈值探索与细胞分配 ===== ####
# 源文件：
#   AUCell/R/03_exploreThresholds.R        —— 主函数
#   AUCell/R/priv_auc.assignmnetThreshold_v6.R  —— 4种阈值算法
#   AUCell/R/aux_AUCell_assignCells.R      —— 手动阈值分配
#   AUCell/R/aux_binarizeAUC.R             —— 二值化
#   AUCell/R/aux_thresholds.R              —— 结果提取辅助函数
#
# 目的：确定每个基因集的"活跃/非活跃"分界点（AUC 阈值），
#       并将细胞分类为"该基因集活跃"或"非活跃"。

# ── 7.1 AUCell_exploreThresholds()：自动探索阈值 ────────────────────────────
#
# 函数签名（完整）：
#
#   AUCell_exploreThresholds(
#     cellsAUC,              # AUC 矩阵（calcAUC 的输出）
#     thrP            = 0.01,  # 异常值概率阈值
#     smallestPopPercent = 0.25, # 最小预期群体比例
#     plotHist        = TRUE,  # 是否绘制诊断直方图
#     densAdjust      = 2,     # 密度曲线平滑度
#     assignCells     = FALSE, # 是否在返回值中包含细胞分配
#     nBreaks         = 100,   # 直方图柱数
#     nCores          = 1,
#     verbose         = TRUE
#   )

## 参数详解

## thrP = 0.01
# 类型：numeric（0 到 1 之间）
# 含义：用于异常值检测的概率阈值（显著性水平）
# 具体作用：
#   对于全局分布方法（Global_k1）：
#     阈值 = qnorm(1 - thrP/nCells, mean, sd)
#     → 在 nCells 个细胞中，期望有 thrP × nCells 个异常值
#     → 例：thrP=0.01，nCells=3000 → 期望约 30 个假阳性
#   对于混合模型方法（L_k2, R_k3）：
#     thrP 用于在混合分布中找各子分布的边界
#
# 调整建议：
#   thrP 越小 → 阈值越高 → 更严格 → 假阳性少，但可能漏掉真实信号
#   thrP 越大 → 阈值越低 → 更宽松 → 灵敏度高，但假阳性多
#   默认 0.01 适合大多数情况

## smallestPopPercent = 0.25
# 类型：numeric（0 到 1 之间）
# 含义：预期最小活跃细胞群体的比例
# 作用：用于过滤不合理的双峰检测
#   如果检测到的第二个峰只包含 < 25% 的细胞，认为不是真实的活跃群体
# 调整建议：
#   稀有细胞类型（如 < 5% 的细胞）→ 需要降低此值，如 0.05
#   常见细胞类型（如 Neuron 占 50%）→ 默认值合适

## plotHist = TRUE
# 类型：logical
# 作用：为每个基因集绘制带阈值线的 AUC 分布直方图
# 直方图颜色编码：
#   蓝色（dodgerblue4）  → 通过选定阈值的细胞
#   灰色（slategray2）   → 未通过阈值的细胞
# 阈值线颜色：
#   蓝色粗实线  → selected（最终选定的阈值）
#   蓝色细虚线  → minimumDens（密度曲线拐点）
#   灰色虚线    → Global_k1（全局正态分布）
#   红色虚线    → L_k2（二分布左侧）
#   粉色虚线    → R_k3（三分布右侧）

## densAdjust = 2
# 类型：numeric（> 0）
# 含义：density() 函数的 adjust 参数，控制核密度估计的平滑程度
# 调整建议：
#   densAdjust < 1 → 更锯齿，捕捉细节峰（噪声也多）
#   densAdjust = 1 → 标准平滑
#   densAdjust = 2 → 更平滑（默认），减少假双峰
#   densAdjust > 3 → 过度平滑，可能错过真实双峰

## assignCells = FALSE
# 类型：logical
# 作用：是否在返回值中包含每个基因集的细胞分配结果
# FALSE → 只返回阈值信息（更快，不存储细胞列表）
# TRUE  → 返回值包含 $assignment（通过阈值的细胞 ID 向量）
# 建议：探索阶段可用 FALSE 快速查看阈值，确定后再用 TRUE 或 assignCells()

# ── 7.2 四种阈值方法详解（源码：priv_auc.assignmnetThreshold_v6.R）────────────
#
# 对每个基因集，代码会尝试计算以下 4 种阈值，
# 然后选择其中"最高的非全局"阈值作为 selected

# 方法 1: minimumDens（蓝色）—— 密度曲线拐点法
#
#   原理：计算 AUC 分布的核密度曲线（density(auc, adjust=densAdjust)），
#         找到曲线的局部最小值（inflection point）
#   实现：diff(sign(diff(density$y)))  → 找到二阶导符号变化的点
#   适用：双峰分布（"低活跃"和"高活跃"两群细胞）
#   局限：单峰分布时找不到有效拐点
#   跳过条件（skipSmallDens）：
#     如果零值或接近零值的细胞太多（前 10% AUC ≈ 0），跳过
#     如果检测到的拐点位于分布右侧更高的峰，跳过（避免假阳性）

# 方法 2: L_k2（红色）—— 双分布混合模型左侧
#
#   原理：拟合 2 个高斯分布的混合模型（mixtools::normalmixEM(k=2)），
#         取左侧（低均值）那个分布，计算其右尾阈值
#   阈值 = qnorm(1 - thrP/nCells, mean_left, sd_left)
#   含义：在"非活跃"细胞群中，只有 thrP/nCells 的概率超过此值
#   适用：明确的双峰分布
#   跳过条件（skipRed）：
#     分布太偏斜（skewness 不对称）
#     两个分布高度重叠（无法区分）

# 方法 3: R_k3（粉色）—— 三分布混合模型右侧
#
#   原理：拟合 3 个高斯分布的混合模型（mixtools::normalmixEM(k=3)），
#         取右侧（高均值）那个分布，计算其左尾阈值
#   阈值 = qnorm(thrP, mean_right, sd_right)
#   含义：在"活跃"细胞群中，只有 thrP 的概率低于此值
#   适用：三峰或复杂分布（如中间有过渡群体）

# 方法 4: Global_k1（灰色）—— 全局正态分布
#
#   原理：将所有细胞的 AUC 值视为单一正态分布，
#         基于均值 + 标准差计算异常值阈值
#   阈值 = qnorm(1 - thrP/nCells, mean_all, sd_all)
#   含义：在整体分布中，位于右尾极端的 thrP 比例细胞
#   通常作为"兜底"：双峰不明确时使用
#   注意：在大多数情况下此方法倾向于被其他方法覆盖

# ── 7.3 selected 的自动选择逻辑 ─────────────────────────────────────────────
#
# 源码逻辑（priv_auc.assignmnetThreshold_v6.R）：
#
#   候选阈值顺序（从高到低排序）：
#     minimumDens > L_k2 > R_k3 > Global_k1
#
#   选择规则：
#     1. 优先选择最高的"非全局"阈值（minimumDens/L_k2/R_k3 之一）
#     2. 如果所有非全局方法都不可用，使用 Global_k1
#     3. 每种方法都有"跳过条件"，不满足时自动跳过
#
# ★ 重要提醒：
#   自动选择 ≠ 正确答案！
#   必须人工查看直方图，理解数据分布后再决定是否接受自动阈值。

# ── 7.4 返回值结构（详细）────────────────────────────────────────────────────
#
# 返回：命名列表（每个基因集一个元素）
#
# 结构：
#   thresholds                         ← 顶层列表
#   └─ $SignatureA                     ← 基因集名称
#      ├─ $aucThr                      ← 阈值信息
#      │  ├─ $selected                 ← 数值标量：自动选定的阈值
#      │  ├─ $thresholds               ← 数据框：所有方法的阈值和对应细胞数
#      │  │     method   threshold  nCells
#      │  │     minimumDens  0.15    120
#      │  │     L_k2         0.12    145
#      │  │     R_k3         0.18     89
#      │  │     Global_k1   0.20     45
#      │  └─ $comment                  ← 字符串：说明为什么选择该方法
#      └─ $assignment                  ← 仅 assignCells=TRUE 时
#            字符向量：通过阈值的细胞名称

# ── 7.5 完整使用示例 ──────────────────────────────────────────────────────────

set.seed(333)   # 固定随机数种子（混合模型初始化有随机性）

cells_assignment <- AUCell_exploreThresholds(
  cellsAUC          = cells_AUC,
  thrP              = 0.01,
  smallestPopPercent = 0.25,
  plotHist          = TRUE,    # 打开，查看直方图！
  densAdjust        = 2,
  assignCells       = TRUE,    # 同时获取细胞分配
  verbose           = TRUE
)

## 查看结果结构
str(cells_assignment, max.level = 3)
names(cells_assignment)          # 所有基因集名称

## 提取某个基因集的信息
gs_result <- cells_assignment[["SignatureA"]]
gs_result$aucThr$selected        # 自动选定的阈值（数值）
gs_result$aucThr$thresholds      # 所有方法的阈值对比
gs_result$assignment             # 通过阈值的细胞名称（字符向量）
length(gs_result$assignment)     # 有多少细胞被认为"活跃"

## 批量提取所有阈值（辅助函数）
thresholds_vec <- getThresholdSelected(cells_assignment)
thresholds_vec  # 命名数值向量
#   SignatureA  SignatureB  SignatureC
#        0.15        0.12        0.20

## 批量提取所有细胞分配
assignments_list <- getAssignments(cells_assignment)
assignments_list$SignatureA  # 字符向量：该基因集活跃的细胞 ID

# ── 7.6 直方图分布形态解读 ───────────────────────────────────────────────────
#
# 形态 1：标准双峰（最理想）
#   外观：左侧低峰 + 右侧高峰，中间有明显谷底
#   解读：明确的两群细胞（活跃 vs 非活跃）
#   阈值策略：minimumDens（谷底）通常是最佳阈值
#   代表：如 Oligodendrocyte 标记在包含神经元的数据中
#
# 形态 2：单峰带右尾
#   外观：左侧高峰，右侧有较长的尾部，少数高 AUC 值
#   解读：大多数细胞不活跃，少数异常高活跃
#   阈值策略：L_k2 或 Global_k1（异常值检测）
#   代表：如 Microglia 在多细胞类型数据中
#
# 形态 3：全部偏右（高 AUC）
#   外观：几乎所有细胞 AUC 都很高，没有明显低值群体
#   解读：该基因集在几乎所有细胞中都活跃（可能是管家基因集）
#   阈值策略：慎用！可能无法有效区分细胞
#   处理：考虑换一个更特异性的基因集
#
# 形态 4：接近随机（中间平坦）
#   外观：AUC 值均匀分布在 0-0.1 之间，类似正态分布
#   解读：该基因集对这个数据集没有区分能力（如随机基因集）
#   处理：该基因集可能不适合此数据
#
# 参考：包自带的 geneSignatures.gmt 中包含随机和管家基因集作为对照

# ── 7.7 手动调整阈值：AUCell_assignCells() ──────────────────────────────────
# 源文件：aux_AUCell_assignCells.R
#
# 函数签名：
#   AUCell_assignCells(
#     cellsAUC,       # AUC 矩阵（aucellResults）
#     thresholds,     # 命名数值向量（名称 = 基因集，值 = 阈值）
#     nCores = 1
#   )
#
# 适用场景：
#   - 查看 exploreThresholds 的直方图后，认为自动阈值不合理
#   - 需要将自动阈值中某些值改为自定义值
#   - 有先验知识（如论文中报告的阈值）

## 获取自动阈值，修改其中一个
new_thresholds <- getThresholdSelected(cells_assignment)
new_thresholds   # 先查看当前自动阈值

# 手动修改某个基因集的阈值
new_thresholds["SignatureA"] <- 0.20   # 提高阈值，更严格

# 也可以手动创建阈值向量
manual_thresholds <- c(
  SignatureA = 0.20,
  SignatureB = 0.15,
  SignatureC = 0.25
)

# 应用新阈值
new_assignments <- AUCell_assignCells(
  cellsAUC   = cells_AUC,
  thresholds = manual_thresholds,
  nCores     = 1
)

# 返回值结构与 exploreThresholds 完全相同
new_assignments$SignatureA$aucThr$selected   # = 0.20（你设定的值）
new_assignments$SignatureA$assignment        # 通过新阈值的细胞

# ── 7.8 二值化 AUC 矩阵：binarizeAUC() ──────────────────────────────────────
# 源文件：aux_binarizeAUC.R
#
# 函数签名：
#   binarizeAUC(auc, thresholds)
#
# 功能：将连续的 AUC 矩阵转化为二值矩阵
#   1 = 细胞通过该基因集的阈值（活跃）
#   0 = 细胞未通过阈值（非活跃）
#
# 内部实现：
#   1. getThresholdSelected(thresholds)    提取每个基因集的阈值
#   2. 对每个基因集：AUC > threshold → 1，否则 0
#   3. reshape2::melt() 长转宽 → table() 转为矩阵

binary_matrix <- binarizeAUC(
  auc        = cells_AUC,
  thresholds = cells_assignment   # 传入 exploreThresholds 的结果
)

dim(binary_matrix)        # nCells × nGeneSets
class(binary_matrix[1,1]) # numeric（0 或 1）

# 查看每个基因集中有多少活跃细胞
colSums(binary_matrix)
#   SignatureA  SignatureB  SignatureC
#          89          67         134


#### ===== Chapter 8: aucellResults 对象操作 ===== ####
# 源文件：AUCell/R/class_aucellResults.R
#         AUCell/R/aux_orderAUC.R
#         AUCell/R/aux_updateAucellResults.R
#
# aucellResults 是 AUCell 的核心 S4 类：
#   继承：SummarizedExperiment
#   扩展槽：nGenesDetected（每细胞基因数的分位数）
#   包含两种 assay 之一：
#     "ranking"  ← buildRankings() 的输出
#     "AUC"      ← calcAUC() 的输出

# ── 8.1 提取核心矩阵 ─────────────────────────────────────────────────────────

## getAUC()：提取 AUC 矩阵
auc_mat <- getAUC(cells_AUC)
class(auc_mat)     # "matrix"
dim(auc_mat)       # c(nGeneSets, nCells)
# 注意：如果对 ranking 对象调用 getAUC() 会报错：
# "This object does not contain an AUC matrix."

## getRanking()：提取排名矩阵
rank_mat <- getRanking(cells_rankings)
class(rank_mat)    # "DelayedMatrix" 或 "matrix"（取决于输入格式）
dim(rank_mat)      # c(nGenes, nCells)
# 注意：如果对 AUC 对象调用 getRanking() 会报错

# ── 8.2 标准维度操作 ─────────────────────────────────────────────────────────

dim(cells_AUC)          # c(3, 300)
nrow(cells_AUC)         # 3（基因集数）
ncol(cells_AUC)         # 300（细胞数）
rownames(cells_AUC)     # 基因集名称
colnames(cells_AUC)     # 细胞 ID
assayNames(cells_AUC)   # "AUC"

# ── 8.3 子集操作 ─────────────────────────────────────────────────────────────
# aucellResults 继承了 SummarizedExperiment 的 [ 方法

## 按基因集行子集（第一维）
cells_AUC[1, ]            # 第一个基因集，所有细胞 → aucellResults（1行）
cells_AUC["SignatureA", ] # 按名称子集
cells_AUC[c("SignatureA","SignatureB"), ]  # 多个基因集

## 按细胞列子集（第二维）
cells_AUC[, 1:100]        # 前 100 个细胞
cells_AUC[, c("Cell1","Cell5","Cell10")]  # 按名称子集

## 同时子集
cells_AUC[1:2, 1:100]    # 前 2 个基因集 × 前 100 个细胞

## 对 ranking 对象的子集
cells_rankings[1:10, 1:5]   # 前 10 个基因 × 前 5 个细胞

# ── 8.4 cbind()：合并多个 aucellResults（按细胞）────────────────────────────
# 源码：class_aucellResults.R，S4 方法 cbind,aucellResults-method
#
# 用途：将多个批次/样本的结果合并
# 约束（源码严格验证）：
#   1. 所有对象必须是同类型（同为 "AUC" 或同为 "ranking"）
#   2. 行数必须相同（同样的基因集或基因）
#   3. 行名必须相同（顺序可不同，但名称必须一一对应）
#   4. 列名不能有重复（细胞 ID 必须唯一）
#
# 内部实现：
#   - 对 assay 矩阵做 cbind（列合并）
#   - 对 nGenesDetected 做 c() 合并

## 示例：合并两个批次
set.seed(1)
exprMat_batch1 <- matrix(
  rpois(5000 * 150, 0.3), 5000, 150,
  dimnames = list(paste0("Gene", 1:5000), paste0("B1_Cell", 1:150))
)
exprMat_batch2 <- matrix(
  rpois(5000 * 150, 0.3), 5000, 150,
  dimnames = list(paste0("Gene", 1:5000), paste0("B2_Cell", 1:150))
)

rankings_b1 <- AUCell_buildRankings(as(exprMat_batch1, "dgCMatrix"),
                                    plotStats = FALSE, verbose = FALSE)
rankings_b2 <- AUCell_buildRankings(as(exprMat_batch2, "dgCMatrix"),
                                    plotStats = FALSE, verbose = FALSE)

auc_b1 <- AUCell_calcAUC(geneSets_demo, rankings_b1, verbose = FALSE)
auc_b2 <- AUCell_calcAUC(geneSets_demo, rankings_b2, verbose = FALSE)

## 合并
auc_combined <- cbind(auc_b1, auc_b2)
dim(auc_combined)      # c(3, 300)（150 + 150 个细胞）
colnames(auc_combined) # "B1_Cell1" ... "B2_Cell150"（全部细胞 ID）

## 常见错误示例：
# cbind(auc_b1, auc_b2[1:2, ])  → ERROR：行数不同
# cbind(rankings_b1, auc_b1)    → ERROR：类型不同（ranking vs AUC）
# 若两批次有重复细胞名 → ERROR：细胞 ID 重复

# rankings 也可以 cbind（用于分批次处理大数据集）
rankings_combined <- cbind(rankings_b1, rankings_b2)
dim(rankings_combined)  # c(5000, 300)

# ── 8.5 orderAUC()：按相似性排序基因集 ─────────────────────────────────────
# 源文件：aux_orderAUC.R
#
# 函数签名：orderAUC(auc)
#
# 目的：对多个基因集按照其 AUC 模式的相似性重新排序，
#       使得功能相似的基因集相邻，便于热图可视化
#
# 内部算法：
#   1. 过滤 AUC 值恒定（SD = 0）的基因集
#   2. 计算基因集间的 Spearman 相关系数矩阵：cor(t(AUC), method="spearman")
#   3. 转换为距离矩阵：1 - cor_matrix
#   4. 层次聚类：hclust(dist, method="ward.D2")
#   5. 动态剪枝：dynamicTreeCut::cutreeDynamic()
#   6. 按聚类标签（降序）+ 树状图顺序排列

sorted_gs_names <- orderAUC(cells_AUC)
sorted_gs_names  # 返回排序后的基因集名称向量

## 按排序结果重排 AUC 对象
cells_AUC_sorted <- cells_AUC[sorted_gs_names, ]

# ── 8.6 updateAucellResults()：升级旧版对象 ──────────────────────────────────
# 源文件：aux_updateAucellResults.R
#
# 函数签名：
#   updateAucellResults(oldAucObject, objectType = "AUC")
#
# 用途：将 AUCell 旧版（< 1.17）的 matrixWrapper 对象转换为现代 aucellResults
# 参数 objectType：
#   "AUC"     → 创建含 "AUC" assay 的 aucellResults
#   "ranking" → 创建含 "ranking" assay 的 aucellResults

# oldAUC <- readRDS("old_aucell_v1.rds")          # 读取旧版结果
# new_auc <- updateAucellResults(oldAUC, "AUC")   # 升级

# ── 8.7 show() 方法：打印摘要 ────────────────────────────────────────────────
# 当直接输入对象名时，自动调用 show()

cells_AUC
# 输出：
# class: aucellResults
# dim: 3 300
# assay names: AUC
# AUC (Gene sets x Cells):
#               Cell1  Cell2  Cell3  ...
#   SignatureA  0.041  0.063  0.029  ...
#   SignatureB  0.055  0.031  0.078  ...
#   SignatureC  0.112  0.088  0.095  ...


#### ===== Chapter 9: 可视化函数 ===== ####
# 源文件：
#   AUCell/R/aux_plotGeneCount.R    —— 基因检测数分布
#   AUCell/R/aux_AUCell_plotHist.R  —— AUC 直方图
#   AUCell/R/aux_AUCell_plotTSNE.R  —— 综合 t-SNE 可视化
#   AUCell/R/priv_plots.R           —— 内部绘图函数
#   AUCell/R/aux_plotEmb_rgb.R      —— RGB 嵌入图

# ── 9.1 plotGeneCount()：基因检测数分布 ──────────────────────────────────────
# 源文件：aux_plotGeneCount.R
#
# 函数签名：
#   plotGeneCount(exprMat, plotStats = TRUE, verbose = TRUE)
#
# 目的：在 buildRankings 之前快速了解数据质量，为 aucMaxRank 选择提供依据
# 内部实现：
#   colSums(exprMat > 0)    → 每个细胞检测到的基因数（非零基因数）
#   对稀疏矩阵：Matrix::colSums()
#
# 参数：
#   exprMat    → 表达矩阵
#   plotStats  → TRUE：绘制箱线图（上）+ 直方图（下）
#   verbose    → TRUE：打印 6 个分位数
#
# 返回值（invisible）：命名数值向量，6 个分位数
#   names: "min" / "1%" / "5%" / "10%" / "50%" / "100%"

## 示例
gene_stats <- plotGeneCount(exprMat_sparse, plotStats = TRUE)
gene_stats
# min  1%   5%  10%  50% 100%
#   0  12   45   89  312  891
# → 中位数 312 个基因 → aucMaxRank 建议 = 312 × 5% ≈ 16（此例数据极稀疏）

# ── 9.2 AUCell_plotHist()：AUC 分布直方图 ────────────────────────────────────
# 源文件：aux_AUCell_plotHist.R
#
# 函数签名：
#   AUCell_plotHist(
#     cellsAUC,                        # aucellResults 或 matrix
#     aucThr    = max(cellsAUC),        # 阈值线位置（默认最大值）
#     nBreaks   = 100,                  # 直方图柱数
#     onColor   = "dodgerblue4",        # 通过阈值的柱颜色
#     offColor  = "slategray2",         # 未通过阈值的柱颜色
#     ...                               # 传给 hist()
#   )
#
# 逻辑（.auc_plotHist 内部）：
#   1. hist() 计算分组，但不立即绘图
#   2. 对每个柱：若右边界 < aucThr → offColor，否则 → onColor
#   3. 重新 plot(hist_obj, col = bar_colors, ...)
#
# 返回值（invisible）：list，每个基因集一个 hist 对象

## 为单个基因集绘制直方图（子集后传入）
par(mfrow = c(1, 3))  # 并排显示 3 个基因集
AUCell_plotHist(
  cellsAUC  = cells_AUC["SignatureA", ],   # 子集单行
  aucThr    = 0.15,                          # 自定义阈值线
  nBreaks   = 50,
  onColor   = "dodgerblue4",
  offColor  = "slategray2",
  main      = "SignatureA"                   # 透传给 hist()
)
abline(v = 0.15, col = "red", lwd = 2)      # 手动画阈值线（更清晰）

## 为所有基因集批量绘制
par(mfrow = c(1, nrow(cells_AUC)))
hist_list <- AUCell_plotHist(
  cellsAUC = cells_AUC,                          # 所有基因集
  aucThr   = getThresholdSelected(cells_assignment)  # 各自的阈值
)

# ── 9.3 AUCell_plotTSNE()：综合 t-SNE 可视化 ─────────────────────────────────
# 源文件：aux_AUCell_plotTSNE.R
#
# 函数签名（完整）：
#   AUCell_plotTSNE(
#     tSNE,                                        # 2D 嵌入坐标（必需）
#     exprMat        = NULL,                       # 表达矩阵（用于 expression 图）
#     cellsAUC       = NULL,                       # AUC 矩阵
#     thresholds     = NULL,                       # 阈值（NULL 则自动计算）
#     reorderGeneSets = FALSE,                     # 按相似性排序基因集
#     cex            = 1,                          # 点大小
#     alphaOn        = 1,                          # 活跃细胞透明度
#     alphaOff       = 0.2,                        # 非活跃细胞透明度
#     borderColor    = adjustcolor("lightgray", alpha.f = 0.1),
#     offColor       = "lightgray",                # 非活跃细胞颜色
#     plots          = c("histogram","binaryAUC","AUC","expression"),
#     exprCols       = c("goldenrod1","darkorange","brown"),  # 表达量渐变色
#     asPNG          = FALSE,                      # 是否保存为 PNG 文件
#     ...
#   )

## 参数 tSNE（必需）
# 类型：matrix，2 列，行名 = 细胞 ID
# 可以是 t-SNE / UMAP / PCA 的任意 2D 坐标
# 行名必须与 cells_AUC 的列名一致

## 参数 plots：控制生成哪些类型的图（可组合）
#
# "histogram"  → AUC 分布直方图（调用 .auc_plotHist）
# "binaryAUC"  → t-SNE 上二值着色（蓝色=活跃，灰色=非活跃）
#                内部调用 .auc_plotBinaryTsne()
# "AUC"        → t-SNE 上连续 AUC 值渐变着色
#                内部调用 .auc_plotGradientTsne()
# "expression" → t-SNE 上该基因集关键 TF 的表达量着色
#                需要 exprMat 参数
#                注意：从基因集名称提取 TF 名（有命名约定限制）

## 参数 thresholds：3 种传入方式
# NULL          → 自动调用 AUCell_exploreThresholds() 计算
# named list    → exploreThresholds() 的完整输出（直接传入）
# named numeric → 手动阈值向量（名称=基因集，值=阈值）

## 参数 alphaOn / alphaOff：透明度控制
# alphaOn  = 1.0  → 活跃细胞完全不透明（鲜艳）
# alphaOff = 0.2  → 非活跃细胞半透明（灰淡背景）
# 调整建议：细胞数很多时，降低 alphaOff 使背景更淡

## 参数 reorderGeneSets = FALSE
# TRUE → 先调用 orderAUC() 按相似性排序基因集，再逐个绘图
# FALSE → 按原顺序绘图

## 参数 asPNG = FALSE
# TRUE → 在当前目录生成 PNG 文件 + HTML 报告（需要 R2HTML 包）
# FALSE → 直接在 R 图形设备绘制

## 示例（使用模拟 t-SNE 坐标）
set.seed(42)
tSNE_demo <- matrix(
  rnorm(300 * 2), nrow = 300, ncol = 2,
  dimnames = list(paste0("Cell", 1:300), c("tSNE_1", "tSNE_2"))
)

# 只绘制直方图和二值 t-SNE（最常用组合）
par(mfrow = c(nrow(cells_AUC), 2))
AUCell_plotTSNE(
  tSNE       = tSNE_demo,
  cellsAUC   = cells_AUC,
  thresholds = cells_assignment,    # 传入已计算的阈值
  plots      = c("histogram", "binaryAUC"),
  alphaOn    = 1,
  alphaOff   = 0.1,
  offColor   = "lightgray"
)

# 完整 4 图（需要表达矩阵用于 expression 图）
# par(mfrow = c(nrow(cells_AUC), 4))
# AUCell_plotTSNE(
#   tSNE       = tSNE_demo,
#   exprMat    = exprMat_sparse,     # 提供表达矩阵
#   cellsAUC   = cells_AUC,
#   thresholds = cells_assignment,
#   plots      = c("histogram", "binaryAUC", "AUC", "expression"),
# )

# ── 9.4 plotEmb_rgb()：3 基因集 RGB 编码可视化 ──────────────────────────────
# 源文件：aux_plotEmb_rgb.R
#
# 函数签名：
#   plotEmb_rgb(
#     aucMat,                    # AUC 矩阵（getAUC 的输出 或 aucellResults）
#     embedding,                 # 2D 嵌入坐标（矩阵，行=细胞）
#     geneSetsByCol,             # 基因集到 RGB 通道的映射
#     aucType        = "AUC",    # "AUC"（连续）或 "binary"（二值）
#     aucMaxContrast = 0.8,      # AUC 对比度调整
#     offColor       = "#c0c0c030",  # 全"关闭"细胞的颜色（带透明度）
#     showPlot       = TRUE,
#     showLegend     = TRUE,
#     ...
#   )
#
# 目的：同时在一张图上展示 3 个基因集的活跃情况
#   红色通道（R）→ 基因集 1
#   绿色通道（G）→ 基因集 2
#   蓝色通道（B）→ 基因集 3
#   混合颜色  →  多个基因集同时活跃（如 R+G = 黄色）

## 参数 geneSetsByCol：两种格式
# 格式 A：字符向量（≤ 3 个元素）→ 按 R/G/B 顺序自动分配
geneSetsByCol_vec <- c("SignatureA", "SignatureB", "SignatureC")

# 格式 B：命名列表 → 明确指定哪个基因集映射到哪个颜色通道
geneSetsByCol_list <- list(
  red   = "SignatureA",    # 神经元 → 红色
  green = "SignatureB",    # 胶质 → 绿色
  blue  = "SignatureC"     # 免疫 → 蓝色
)

# 可以将多个基因集分组到一个通道（取平均 AUC）
geneSetsByCol_grouped <- list(
  red   = c("SignatureA", "SignatureB"),  # 多个基因集 → 红色
  green = "SignatureC"                    # 单个 → 绿色
  # 省略 blue → 不显示蓝色通道
)

## 参数 aucType
# "AUC" / "continuous" / "cont" → 用连续 AUC 值映射颜色强度
# "binary" / "bin"              → 先二值化（>0 = 活跃），再映射
#
# 连续模式内部：
#   channel_value = mean(AUC[genesets_in_channel, cell]) / aucMaxContrast
#   → 除以 aucMaxContrast 后 clip 到 [0, 1]
#
# 二值模式内部：
#   先用 binarizeAUC() 得到 0/1，再取该通道基因集的平均值
#   → 值为 0（全非活跃）到 1（全活跃）之间

## 参数 aucMaxContrast = 0.8
# 调整颜色对比度（连续模式才生效）
# AUC 除以 aucMaxContrast 后 clip 到 [0, 1]
# 0.8 → AUC ≥ 0.8 的细胞颜色饱和（最亮）
# 降低（如 0.3）→ 整体颜色更明亮（适合 AUC 值普遍较低的数据）
# 提高（如 1.0）→ 只有最高 AUC 细胞才颜色鲜艳

## 参数 offColor = "#c0c0c030"
# 所有通道都接近 0 的细胞显示此颜色
# 默认：半透明浅灰色（#c0c0c0 = silver，30 = alpha 约 0.19）
# 建议：保持半透明，使活跃细胞的颜色更突出

## 示例
auc_mat <- getAUC(cells_AUC)
par(mfrow = c(1, 2))

# 连续模式
plotEmb_rgb(
  aucMat        = auc_mat,
  embedding     = tSNE_demo,
  geneSetsByCol = geneSetsByCol_list,
  aucType       = "AUC",
  aucMaxContrast = 0.8,
  offColor      = "#c0c0c030",
  showPlot      = TRUE,
  showLegend    = TRUE,
  main          = "RGB AUC (连续)"
)

# 二值模式
plotEmb_rgb(
  aucMat        = auc_mat,
  embedding     = tSNE_demo,
  geneSetsByCol = geneSetsByCol_list,
  aucType       = "binary",
  offColor      = "#c0c0c030",
  main          = "RGB AUC (二值)"
)

## 返回值（invisible）：字符向量，每个细胞一个 hex 颜色码
# 带 attributes：
#   attr(result, "red")   → 红色通道的基因集名称
#   attr(result, "green") → 绿色通道的基因集名称
#   attr(result, "blue")  → 蓝色通道的基因集名称

# ── 9.5 plotTsne_cellProps()：细胞元数据着色 ─────────────────────────────────
# 源文件：priv_plots.R
#
# 函数签名：
#   plotTsne_cellProps(
#     tSNE,
#     cellInfo,                        # 细胞元数据（data.frame，行=细胞）
#     colVars    = NULL,                # 分类变量的颜色映射列表
#     cex        = 1,
#     sub        = "",
#     gradientCols = c("yellow","orange","red"),
#     showLegend = TRUE
#   )
#
# 用途：在 t-SNE 上展示已知的细胞注释信息（如细胞类型、cluster）
# 自动判断变量类型：
#   character / factor → 分类着色（rainbow 调色板 或 colVars 中的自定义颜色）
#   numeric           → 渐变着色（gradientCols 调色板）

## 示例
cellInfo_demo <- data.frame(
  CellType  = sample(c("Neuron","Astro","Oligo"), 300, replace = TRUE),
  nGenes    = round(rnorm(300, mean = 2000, sd = 500)),
  row.names = paste0("Cell", 1:300)
)

# 自定义颜色（分类变量）
colVars <- list(
  CellType = c(
    Neuron = "orange",
    Astro  = "cornflowerblue",
    Oligo  = "forestgreen"
  )
)

par(mfrow = c(1, 2))
plotTsne_cellProps(
  tSNE     = tSNE_demo,
  cellInfo = cellInfo_demo,
  colVars  = colVars,
  showLegend = TRUE
)


#### ===== Chapter 10: 完整端到端示例 ===== ####
# 使用包内自带数据，完整演示从原始数据到结果可视化的全流程

# ── 10.1 准备数据 ─────────────────────────────────────────────────────────────

library(AUCell)
library(GSEABase)
library(Matrix)

## 构建模拟单细胞表达矩阵
# 参数设定：5种细胞类型，每种60个细胞，共1500个基因
set.seed(123)
nCells_per_type <- 60
cell_types <- c("Neuron", "Astro", "Oligo", "Microglia", "OPC")
nTypes <- length(cell_types)
nGenes_total <- 1500

# 每种细胞类型的标记基因（各 50 个）
marker_genes <- lapply(seq_len(nTypes), function(i) {
  paste0("Gene", ((i-1)*50 + 1):(i*50))
})
names(marker_genes) <- cell_types

# 生成表达矩阵：标记基因在对应细胞中高表达
exprMat_full <- matrix(
  rpois(nGenes_total * (nCells_per_type * nTypes), lambda = 0.2),
  nrow = nGenes_total,
  ncol = nCells_per_type * nTypes,
  dimnames = list(
    paste0("Gene", 1:nGenes_total),
    paste0(rep(cell_types, each = nCells_per_type), "_",
           rep(1:nCells_per_type, nTypes))
  )
)

# 增加标记基因表达量（模拟细胞类型特异性）
for (i in seq_len(nTypes)) {
  cell_cols <- which(grepl(paste0("^", cell_types[i], "_"), colnames(exprMat_full)))
  gene_rows <- which(rownames(exprMat_full) %in% marker_genes[[i]])
  exprMat_full[gene_rows, cell_cols] <- rpois(
    length(gene_rows) * length(cell_cols), lambda = 5
  )
}

# 转为稀疏矩阵
exprMat_sp <- as(exprMat_full, "dgCMatrix")
cat("表达矩阵维度：", nrow(exprMat_sp), "基因 ×", ncol(exprMat_sp), "细胞\n")

# 读取基因集（用包内自带 gmt 文件 + 自定义标记基因集）
gmtFile <- system.file("examples", "geneSignatures.gmt", package = "AUCell")
geneSets_pkg <- getGmt(gmtFile)
cat("包内基因集：", length(geneSets_pkg), "个\n")

# 补充自定义基因集（使用模拟数据的标记基因）
geneSets_custom <- GeneSetCollection(lapply(seq_len(nTypes), function(i) {
  GeneSet(marker_genes[[i]], setName = paste0(cell_types[i], "_markers"))
}))

# 合并基因集（以自定义为主）
all_geneSets <- geneSets_custom

# 过滤：只保留存在于表达矩阵中的基因
all_geneSets <- subsetGeneSets(all_geneSets, rownames(exprMat_sp))
cat("过滤后基因集大小：\n")
print(cbind(nGenes = nGenes(all_geneSets)))

# ── 10.2 Step 1: 构建排序 ─────────────────────────────────────────────────────

cat("\n=== Step 1: 构建基因排序 ===\n")
cells_rankings_full <- AUCell_buildRankings(
  exprMat   = exprMat_sp,
  plotStats = TRUE,         # ★ 必看图！确认 aucMaxRank 选择
  verbose   = TRUE
)

# 从图中观察：50% 分位数 ≈ X 个基因
# aucMaxRank 的默认 5% 是否合理？

cat("排序对象维度：", dim(cells_rankings_full), "\n")

# ── 10.3 Step 2: 计算 AUC ─────────────────────────────────────────────────────

cat("\n=== Step 2: 计算 AUC 分数 ===\n")

# 根据 Step 1 的图确定 aucMaxRank
# 假设 50% 分位数 ≈ 250 个基因，则 5% ≈ 13，适当放宽到 10%
aucMaxRank_chosen <- ceiling(0.10 * nrow(cells_rankings_full))
cat("使用 aucMaxRank =", aucMaxRank_chosen, "\n")

cells_AUC_full <- AUCell_calcAUC(
  geneSets   = all_geneSets,
  rankings   = cells_rankings_full,
  aucMaxRank = aucMaxRank_chosen,
  normAUC    = TRUE,
  verbose    = TRUE
)

cat("AUC 矩阵维度：", dim(cells_AUC_full), "\n")

# 快速查看 AUC 分布
auc_full <- getAUC(cells_AUC_full)
cat("\nAUC 值统计（各基因集）：\n")
print(round(t(apply(auc_full, 1, summary)), 3))

# ── 10.4 Step 3: 探索阈值 ─────────────────────────────────────────────────────

cat("\n=== Step 3: 探索阈值 ===\n")

set.seed(333)
cells_assignment_full <- AUCell_exploreThresholds(
  cellsAUC          = cells_AUC_full,
  thrP              = 0.01,
  smallestPopPercent = 0.25,
  plotHist          = TRUE,    # ★ 仔细看每个直方图！
  densAdjust        = 2,
  assignCells       = TRUE
)

# 查看自动选定的阈值
auto_thresholds <- getThresholdSelected(cells_assignment_full)
cat("\n自动选定阈值：\n")
print(round(auto_thresholds, 4))

# 查看每个基因集的细胞分配数
n_assigned <- sapply(getAssignments(cells_assignment_full), length)
cat("\n各基因集活跃细胞数：\n")
print(n_assigned)

# ── 10.5 手动调整阈值（如有需要）─────────────────────────────────────────────

cat("\n=== 手动调整阈值（示例）===\n")

# 假设看图后认为 Neuron_markers 的阈值太低，手动提高
adjusted_thresholds <- auto_thresholds
# adjusted_thresholds["Neuron_markers"] <- 0.30  # 按需调整

final_assignments <- AUCell_assignCells(
  cellsAUC   = cells_AUC_full,
  thresholds = adjusted_thresholds
)

# ── 10.6 结果可视化 ──────────────────────────────────────────────────────────

cat("\n=== 可视化结果 ===\n")

# 生成模拟 UMAP/t-SNE 坐标（实际分析中用真实的降维结果）
set.seed(42)
tsne_coords <- matrix(
  rnorm(ncol(exprMat_sp) * 2),
  nrow = ncol(exprMat_sp), ncol = 2,
  dimnames = list(colnames(exprMat_sp), c("UMAP_1", "UMAP_2"))
)

# 图 A：AUC 直方图（每个基因集一张）
par(mfrow = c(2, 3))
AUCell_plotHist(
  cellsAUC = cells_AUC_full,
  aucThr   = auto_thresholds
)

# 图 B：t-SNE 上的二值着色
par(mfrow = c(2, 3))
AUCell_plotTSNE(
  tSNE       = tsne_coords,
  cellsAUC   = cells_AUC_full,
  thresholds = final_assignments,
  plots      = c("binaryAUC"),
  cex        = 0.8,
  alphaOn    = 1,
  alphaOff   = 0.05
)

# 图 C：RGB 多基因集（只展示前 3 个）
auc_full_mat <- getAUC(cells_AUC_full)
plotEmb_rgb(
  aucMat        = auc_full_mat,
  embedding     = tsne_coords,
  geneSetsByCol = list(
    red   = rownames(auc_full_mat)[1],
    green = rownames(auc_full_mat)[2],
    blue  = rownames(auc_full_mat)[3]
  ),
  aucType        = "AUC",
  aucMaxContrast = 0.8,
  main           = "三基因集 RGB 叠加"
)

# ── 10.7 结果解读与常见问题 ──────────────────────────────────────────────────

# Q1: AUC 全是 0 或接近 0？
#   → 检查基因集和表达矩阵是否用同一套基因 ID
#   → 检查 aucMaxRank 是否太小（从 plotGeneCount 的图来判断）
#   → 代码：
#   sum(marker_genes[["Neuron"]] %in% rownames(exprMat_sp))
#   # 应该 > 0，且比例越高越好

# Q2: 阈值直方图是单峰，没有明显分界？
#   → 该基因集可能对此数据无区分能力
#   → 尝试：subsetGeneSets() 重新过滤 → 减少背景基因
#   → 或：考虑换一个更特异的基因集

# Q3: 多样本怎么处理？
#   → 分别运行 calcAUC，再用 cbind() 合并（见 Ch8）
#   → 注意：不要混合不同数据集的 rankings，因为排名是细胞内相对的

# Q4: 如何保存结果？
saveRDS(cells_AUC_full,        "cells_AUC.rds")
saveRDS(cells_assignment_full, "cells_assignment.rds")
write.csv(getAUC(cells_AUC_full), "AUC_matrix.csv")

# Q5: 如何与 Seurat 结合？
# library(Seurat)
# # 从 Seurat 提取计数矩阵
# exprMat_seurat <- GetAssayData(seurat_obj, slot = "counts")
# # 运行 AUCell
# cells_AUC_seurat <- AUCell_run(exprMat_seurat, all_geneSets)
# # 将 AUC 值加回 Seurat（每个基因集作为一个 meta feature）
# for (gs in rownames(cells_AUC_seurat)) {
#   seurat_obj[[gs]] <- getAUC(cells_AUC_seurat)[gs, colnames(seurat_obj)]
# }


################################################################################
##
##  附录：函数快速参考卡
##
##  函数                      关键参数                     返回值
##  ─────────────────────────────────────────────────────────────────────────
##  AUCell_buildRankings()    exprMat, aucMaxRank*,        aucellResults(ranking)
##                            plotStats, BPPARAM
##
##  AUCell_calcAUC()          geneSets, rankings,          aucellResults(AUC)
##                            aucMaxRank★, normAUC
##
##  AUCell_run()              exprMat, geneSets,           aucellResults(AUC)
##                            aucMaxRank★, BPPARAM
##
##  AUCell_exploreThresholds() cellsAUC, thrP,             named list
##                            plotHist, assignCells        (aucThr + assignment)
##
##  AUCell_assignCells()      cellsAUC, thresholds         named list
##
##  binarizeAUC()             auc, thresholds              binary matrix
##
##  getAUC()                  aucellResults                matrix
##  getRanking()              aucellResults                matrix
##  getThresholdSelected()    thresholds_list              named numeric
##  getAssignments()          thresholds_list              named list of cells
##
##  AUCell_plotHist()         cellsAUC, aucThr             invisible(hist list)
##  AUCell_plotTSNE()         tSNE, cellsAUC, thresholds   list + plots
##  plotEmb_rgb()             aucMat, embedding,           invisible(color vec)
##                            geneSetsByCol, aucType
##  plotGeneCount()           exprMat                      named numeric(6分位数)
##  orderAUC()                auc                          character vector
##
##  subsetGeneSets()          geneSets, geneNames          GeneSetCollection
##  nGenes()                  geneSets                     integer/vector
##  setGeneSetNames()         geneSets, newNames           GeneSetCollection
##  getSetNames()             aucMat, patterns             named character
##  cbind(aucellResults,...)  多个 aucellResults            aucellResults（合并）
##  updateAucellResults()     oldObject, objectType        aucellResults
##
##  ★ aucMaxRank 是全流程最重要的参数，必须结合 plotGeneCount 图来选择
##
################################################################################
