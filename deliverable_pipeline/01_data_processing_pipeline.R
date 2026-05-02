#!/usr/bin/env Rscript

# Stage 1 deliverable pipeline: data preparation and processing
# NOTE: plotting/statistical parameters remain unchanged because this pipeline
# simply executes original stage scripts in reproducible order.

source("deliverable_pipeline/_shared.R", encoding = "UTF-8")

scripts <- c(
  "01_data_processing/数据准备.R",
  "01_data_processing/偏相关网络构建.R",
  "01_data_processing/社区发现.R",
  "01_data_processing/社区映射回个体.R",
  "01_data_processing/主导社区分配.R",
  "01_data_processing/QC_Step4_快速检查.R"
)

log_step("Stage 1: Data processing started")
invisible(lapply(scripts, run_script))
log_step("Stage 1: Data processing completed")
