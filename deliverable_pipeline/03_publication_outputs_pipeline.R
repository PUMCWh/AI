#!/usr/bin/env Rscript

# Stage 3 deliverable pipeline: publication-grade tables and figures
# Keeps all original plotting parameters unchanged by executing source scripts as-is.

source("deliverable_pipeline/_shared.R", encoding = "UTF-8")

scripts <- c(
  "03_publication_outputs/生成正文表1.R",
  "03_publication_outputs/生成正文表2.R",
  "03_publication_outputs/生成正文图3.R",
  "03_publication_outputs/生成出版级森林图（基于cox回归）.R",
  "03_publication_outputs/STROBE Flow Diagram for Figure.R",
  "03_publication_outputs/Fig S3 Top 100 diseases by Multimorbidity Centrality (MMC), circular bar chart.R",
  "03_publication_outputs/出版级分辨率扫描结果图.R",
  "03_publication_outputs/Fig 3 Kaplan-Meier 生存曲线 (Publication grade BMC Medicine 投稿版.R",
  "03_publication_outputs/Generate Supplementary Tables S2 - S9 in publication format.R",
  "03_publication_outputs/# Step 15-16-17 三连发 — 补表 S8、补表 S9、补图 S7.R"
)

log_step("Stage 3: Publication outputs started")
invisible(lapply(scripts, run_script))
log_step("Stage 3: Publication outputs completed")
