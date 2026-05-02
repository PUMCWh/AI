#!/usr/bin/env Rscript
# ==============================================================================
# Generate Table 1: Baseline characteristics (Lancet-style)
#
# Columns: All / Female / Male
# Rows:    Age (mean SD) / Age groups / CKM domains / Baseline dx /
#          Follow-up / Deaths / Mortality rate
#
# Output: Table1_baseline_characteristics.xlsx
# ==============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)       # for parquet
  library(openxlsx)
})

# =========================
# Paths
# =========================
COX_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
PARQUET <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\analysis_CKM_strict.parquet"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\manuscript_tables"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) 加载 cox 数据集 (患者级)
# =========================
cat("==> Step 1) 加载 cox 数据集...\n")
cox <- fread(COX_CSV,
             select = c("pid", "age_at_index", "sex_binary",
                        "time_years", "event", "baseline_unique_codes",
                        "dom_community_unw"))
cat(sprintf("   原始 cohort: %d 行\n", nrow(cox)))

# 排除 NONE
cox <- cox[dom_community_unw != "NONE"]
cat(sprintf("   排除 NONE 后: %d 行\n", nrow(cox)))

# =========================
# 2) 加载 parquet, 提取每个患者的 CKM_C/K/M
# =========================
cat("\n==> Step 2) 加载 parquet 数据 (CKM 域标识)...\n")

# parquet 是诊断级, 一个患者多条; 但 CKM_C/K/M 是患者级常量
# 用 arrow 高效读取去重
pq <- arrow::open_dataset(PARQUET)
ckm <- as.data.table(pq |>
                       dplyr::select(pid, CKM_C, CKM_K, CKM_M, n_ckm_domains) |>
                       dplyr::distinct() |>
                       dplyr::collect())

# 防御: 同一 pid 应该只有一行 (常量), 万一有重复则去重
ckm <- unique(ckm, by = "pid")
cat(sprintf("   唯一患者数 (parquet): %d\n", nrow(ckm)))

# Merge into cox
cox <- merge(cox, ckm, by = "pid", all.x = TRUE)
cox[is.na(CKM_C),         CKM_C := 0L]
cox[is.na(CKM_K),         CKM_K := 0L]
cox[is.na(CKM_M),         CKM_M := 0L]
cox[is.na(n_ckm_domains), n_ckm_domains := 0L]

cat(sprintf("   合并后: %d 行 (CKM_C=%d, CKM_K=%d, CKM_M=%d, >=2 域=%d)\n",
            nrow(cox),
            sum(cox$CKM_C == 1), sum(cox$CKM_K == 1),
            sum(cox$CKM_M == 1), sum(cox$n_ckm_domains >= 2)))

# =========================
# 3) 派生 age groups
# =========================
cox[, age_group := fcase(
  age_at_index < 45,                       "18-44 years",
  age_at_index >= 45 & age_at_index < 65,  "45-64 years",
  age_at_index >= 65,                       ">=65 years"
)]
cox[, age_group := factor(age_group,
                          levels = c("18-44 years", "45-64 years",
                                     ">=65 years"))]

# Sex labels (sex_binary: 0/1; 通常 0=Female, 1=Male, 但需要确认)
# 从 数据准备.R: line 128 — sex 列从原始 char 来, sex_binary 在 Step 4.5 派生
# 我们假设 sex_binary == 1 是 Male (常用约定)
cox[, sex_label := fifelse(sex_binary == 1L, "Male", "Female")]

# =========================
# 4) 计算 Table 1 各 cell
# =========================
cat("\n==> Step 3) 计算 Table 1...\n")

# 通用辅助
fmt_n_pct <- function(n, total) {
  if (total == 0) return("0 (0.0)")
  sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * n / total)
}
fmt_mean_sd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}
fmt_median_iqr <- function(x, dig = 0) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", dig, "f (%.", dig, "f-%.", dig, "f)"),
          q[1], q[2], q[3])
}

# 一个函数算所有 metrics 给定 subset
make_col <- function(dt) {
  N <- nrow(dt)
  c(
    n_header           = sprintf("%s", format(N, big.mark = ",")),
    age_mean_sd        = fmt_mean_sd(dt$age_at_index),
    age_18_44          = fmt_n_pct(sum(dt$age_group == "18-44 years"), N),
    age_45_64          = fmt_n_pct(sum(dt$age_group == "45-64 years"), N),
    age_65_plus        = fmt_n_pct(sum(dt$age_group == ">=65 years"), N),
    ckm_C              = fmt_n_pct(sum(dt$CKM_C == 1L), N),
    ckm_K              = fmt_n_pct(sum(dt$CKM_K == 1L), N),
    ckm_M              = fmt_n_pct(sum(dt$CKM_M == 1L), N),
    ckm_2plus          = fmt_n_pct(sum(dt$n_ckm_domains >= 2L), N),
    baseline_dx        = fmt_median_iqr(dt$baseline_unique_codes, 0),
    followup           = fmt_median_iqr(dt$time_years, 2),
    deaths             = fmt_n_pct(sum(dt$event == 1L), N),
    mort_rate          = sprintf("%.2f",
                                 1000 * sum(dt$event == 1L) / sum(dt$time_years))
  )
}

# 计算 3 列
col_all    <- make_col(cox)
col_female <- make_col(cox[sex_label == "Female"])
col_male   <- make_col(cox[sex_label == "Male"])

# 验证 sex 假设 (女性占比通常 50-60%)
n_f <- sum(cox$sex_label == "Female")
n_m <- sum(cox$sex_label == "Male")
cat(sprintf("   Sex 分布: Female=%s (%.1f%%), Male=%s (%.1f%%)\n",
            format(n_f, big.mark = ","), 100 * n_f / nrow(cox),
            format(n_m, big.mark = ","), 100 * n_m / nrow(cox)))
if (n_f / nrow(cox) < 0.30 || n_f / nrow(cox) > 0.70) {
  cat("   提示: 女性占比异常, 检查 sex_binary 编码方向 (0/1 vs Female/Male)\n")
}

# =========================
# 5) 组装 Table 1
# =========================
header_row <- c("Characteristic", "All", "Female", "Male")
data_rows <- list(
  c("Total participants, n",                col_all["n_header"],   col_female["n_header"], col_male["n_header"]),
  c("Age, mean (SD), years",                col_all["age_mean_sd"], col_female["age_mean_sd"], col_male["age_mean_sd"]),
  c("Age group, n (%)",                     "",  "",  ""),
  c("  18-44 years",                        col_all["age_18_44"], col_female["age_18_44"], col_male["age_18_44"]),
  c("  45-64 years",                        col_all["age_45_64"], col_female["age_45_64"], col_male["age_45_64"]),
  c("  >=65 years",                         col_all["age_65_plus"], col_female["age_65_plus"], col_male["age_65_plus"]),
  c("CKM domains, n (%)",                   "",  "",  ""),
  c("  Cardiovascular",                     col_all["ckm_C"], col_female["ckm_C"], col_male["ckm_C"]),
  c("  Kidney",                             col_all["ckm_K"], col_female["ckm_K"], col_male["ckm_K"]),
  c("  Metabolic",                          col_all["ckm_M"], col_female["ckm_M"], col_male["ckm_M"]),
  c("  >=2 domains",                        col_all["ckm_2plus"], col_female["ckm_2plus"], col_male["ckm_2plus"]),
  c("Distinct baseline diagnoses, median (IQR)", col_all["baseline_dx"], col_female["baseline_dx"], col_male["baseline_dx"]),
  c("Follow-up time, median (IQR), years", col_all["followup"], col_female["followup"], col_male["followup"]),
  c("All-cause deaths, n (%)",              col_all["deaths"], col_female["deaths"], col_male["deaths"]),
  c("Mortality rate per 1000 person-years", col_all["mort_rate"], col_female["mort_rate"], col_male["mort_rate"])
)

# 转 data.table
table1_dt <- rbindlist(lapply(data_rows, function(r) {
  data.table(Characteristic = r[1], All = r[2],
             Female = r[3], Male = r[4])
}))

# Header line as proper column names (with N counts)
n_label <- function(N) sprintf("(n = %s)", format(N, big.mark = ","))
setnames(table1_dt,
         old = c("All", "Female", "Male"),
         new = c(sprintf("All %s", n_label(nrow(cox))),
                 sprintf("Female %s", n_label(n_f)),
                 sprintf("Male %s", n_label(n_m))))

# =========================
# 6) 写入 xlsx
# =========================
out_xlsx <- file.path(OUT_DIR, "Table1_baseline_characteristics.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Table 1")
writeData(wb, "Table 1", table1_dt, startRow = 2, headerStyle = createStyle(
  textDecoration = "bold", border = "TopBottom", borderStyle = "thin",
  halign = "center"
))

# 标题
title_text <- "Table 1. Baseline characteristics of the analytical cohort"
writeData(wb, "Table 1", title_text, startRow = 1, startCol = 1)
addStyle(wb, "Table 1", createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)

# Footnote
foot_row <- nrow(table1_dt) + 3
foot_text <- paste0(
  "Data are mean (SD), median (IQR), or n (%). ",
  "CKM = cardiovascular-kidney-metabolic. ",
  "Cardiovascular, kidney, and metabolic domains are not mutually exclusive; ",
  "the '\u22652 domains' row indicates participants with conditions in two or more domains. ",
  "Mortality rate is calculated as deaths divided by total person-years of follow-up, ",
  "multiplied by 1000."
)
writeData(wb, "Table 1", foot_text,
          startRow = foot_row, startCol = 1)
addStyle(wb, "Table 1", createStyle(fontSize = 9,
                                    textDecoration = "italic",
                                    wrapText = TRUE),
         rows = foot_row, cols = 1)

# 列宽
setColWidths(wb, "Table 1", cols = 1, widths = 50)
setColWidths(wb, "Table 1", cols = 2:4, widths = 25)

# 保存
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n   \u2713 已保存: %s\n", out_xlsx))

# =========================
# 7) Console 输出
# =========================
cat("\n========================================================\n")
cat("Table 1 — preview:\n")
cat("========================================================\n")
print(table1_dt)
cat(sprintf("\n保存路径: %s\n", out_xlsx))