#!/usr/bin/env Rscript

# Stage 2 deliverable pipeline: analysis and preliminary figures/tables
# Keeps original script internals untouched.

source("deliverable_pipeline/_shared.R", encoding = "UTF-8")

scripts <- c(
  "02_analysis/cox回归分析.R",
  "02_analysis/亚组与交互作用分析.R",
  "02_analysis/限制性立方样条分析.R",
  "02_analysis/人群归因分析.R",
  "02_analysis/ICD-10 Chapter Enrichment Analysis for CKM Disease Communities.R",
  "02_analysis/Node Role Classification (Guimerà-Amaral framework).R",
  "02_analysis/Normalized Connectivity Index (NCI) Matrix across 16 Disease Communities.R",
  "02_analysis/Bootstrap C-statistic and ΔC 95% Confidence Intervals.R",
  "02_analysis/Top1-share distribution with τ = 0.40 threshold justification.R",
  "02_analysis/Leiden Resolution Sweep.R",
  "02_analysis/社区ICD章节富集分析.R",
  "02_analysis/社区间标准化连接强度热图.R",
  "02_analysis/桥接疾病的节点角色分类.R"
)

log_step("Stage 2: Analysis started")
invisible(lapply(scripts, run_script))
log_step("Stage 2: Analysis completed")
