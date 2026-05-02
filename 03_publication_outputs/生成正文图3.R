#!/usr/bin/env Rscript
# ==============================================================================
# Generate Table 3: Hazard ratios for 4-year all-cause mortality
#
# Lancet/BMC Medicine style — 3 nested models in 3 columns
#
# Decision (locked):
#   Paper Table 3 displays 3 models:
#     - Paper "Model 1" = Cox data M1 (Unadjusted)
#     - Paper "Model 2" = Cox data M2 (+ age + sex)
#     - Paper "Model 3" = Cox data M4 (+ age + sex + log baseline diagnoses count)
#   M3 (community burden score adjusted) is OMITTED from Table 3 to avoid
#   collider bias (total_score_unw is derived from dominant community phenotype).
#
# Rows (17): 16 phenotypes (DOM_C1-C16) + MIXED, with C7 as reference (first row)
# Cols (7):  Phenotype | n | Events | Rate per 1000 py | M1 HR | M2 HR | M3 HR
#
# Inputs:
#   - cox_regression_results.xlsx
#       Main_results sheet     (HR data, model labels)
#       Group_sizes_main sheet (n, events, rate per group)
#
# Output:
#   - Table3_cox_hazard_ratios.xlsx
# ==============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# =========================
# Paths
# =========================
COX_XLSX <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results\\cox_regression_results.xlsx"
OUT_DIR  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\manuscript_tables"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# Phenotype labels (与 Table 2, Fig 4 一致)
# =========================
phenotype_labels <- c(
  "DOM_C1"  = "C1 \u2014 Functional & inflammatory GI disorders",
  "DOM_C2"  = "C2 \u2014 Chronic respiratory disease with infections",
  "DOM_C3"  = "C3 \u2014 Dermatologic & venous insufficiency disorders",
  "DOM_C4"  = "C4 \u2014 Benign gynecologic disorders",
  "DOM_C5"  = "C5 \u2014 Age-related ophthalmologic disorders",
  "DOM_C6"  = "C6 \u2014 Cerebrovascular disease with neuropsychiatric sequelae",
  "DOM_C7"  = "C7 \u2014 Upper-airway & orodental inflammation",
  "DOM_C8"  = "C8 \u2014 Hepatobiliary disorders",
  "DOM_C9"  = "C9 \u2014 Advanced solid malignancies with cachexia",
  "DOM_C10" = "C10 \u2014 Benign & malignant breast disorders",
  "DOM_C11" = "C11 \u2014 Urinary stones, BPH & UTI",
  "DOM_C12" = "C12 \u2014 Thyroid dysfunction & neoplasms",
  "DOM_C13" = "C13 \u2014 Lymphoid & myeloid malignancies",
  "DOM_C14" = "C14 \u2014 CKD with anemia, electrolyte & gout",
  "DOM_C15" = "C15 \u2014 T2DM with neuropathy & joint disease",
  "DOM_C16" = "C16 \u2014 CAD, heart failure & arrhythmia",
  "MIXED"   = "MIXED"
)

# Display order: C1-C16 by index, then MIXED, with C7 (Ref.) FIRST
display_order <- c("DOM_C7",
                   paste0("DOM_C", c(1:6, 8:16)),
                   "MIXED")

# =========================
# 1) 加载 Group_sizes_main (n, events, rate)
# =========================
cat("==> Step 1) 加载 Group_sizes_main...\n")
gs <- as.data.table(read.xlsx(COX_XLSX, sheet = "Group_sizes_main"))
cat(sprintf("   行数: %d\n", nrow(gs)))
cat("   列名: ", paste(names(gs), collapse = ", "), "\n", sep = "")

# Rename for clarity
setnames(gs, c("group", "N", "events", "rate_1000py"),
         c("term", "n_grp", "events_grp", "rate_grp"))

# =========================
# 2) 加载 Main_results (HR data)
# =========================
cat("\n==> Step 2) 加载 Main_results...\n")
mr <- as.data.table(read.xlsx(COX_XLSX, sheet = "Main_results"))
cat(sprintf("   行数: %d (期望 64 = 16 communities x 4 models)\n", nrow(mr)))
cat("   列名: ", paste(names(mr), collapse = ", "), "\n", sep = "")

# 仅保留 HR + 95% CI + p_value + model + term
mr_select <- mr[, .(term, HR, HR_lower, HR_upper, p_value, model)]

# =========================
# 3) Filter to 3 paper models: M1 (unadjusted), M2 (+age+sex), M4 (+age+sex+codes)
# =========================
# 注意: model 列实际值如 "Model 1: Unadjusted", "Model 2: + age + sex", etc.
m1_pattern <- "Model 1: Unadjusted"
m2_pattern <- "Model 2: + age + sex"
m4_pattern <- "Model 4: + age + sex + codes"

cat("\n==> Step 3) Filter 3 papers' models...\n")
cat(sprintf("   M1 (paper Model 1): rows = %d\n",
            nrow(mr_select[model == m1_pattern])))
cat(sprintf("   M2 (paper Model 2): rows = %d\n",
            nrow(mr_select[model == m2_pattern])))
cat(sprintf("   M4 (paper Model 3): rows = %d\n",
            nrow(mr_select[model == m4_pattern])))

# =========================
# 4) 格式化 HR (95% CI) for 每个 model
# =========================
fmt_hr_ci <- function(hr, lo, hi) {
  ifelse(is.na(hr) | is.na(lo) | is.na(hi),
         "\u2014",  # em dash for NA
         sprintf("%.2f (%.2f\u2013%.2f)", hr, lo, hi))
}

# Reshape: one row per term, 3 columns for M1/M2/M3 HR
m1_dt <- mr_select[model == m1_pattern,
                   .(term,
                     hr_m1 = fmt_hr_ci(HR, HR_lower, HR_upper))]
m2_dt <- mr_select[model == m2_pattern,
                   .(term,
                     hr_m2 = fmt_hr_ci(HR, HR_lower, HR_upper))]
m3_dt <- mr_select[model == m4_pattern,    # Cox data M4 = paper M3
                   .(term,
                     hr_m3 = fmt_hr_ci(HR, HR_lower, HR_upper))]

# Merge HR data with group sizes
table3_dt <- merge(gs[, .(term, n_grp, events_grp, rate_grp)],
                   m1_dt, by = "term", all.x = TRUE)
table3_dt <- merge(table3_dt, m2_dt, by = "term", all.x = TRUE)
table3_dt <- merge(table3_dt, m3_dt, by = "term", all.x = TRUE)

# =========================
# 5) 排序 + REF 处理
# =========================
table3_dt[, term := factor(term, levels = display_order)]
setorder(table3_dt, term)
table3_dt[, term := as.character(term)]

# C7 是 reference: HR 列全部置为 "1.00 (Reference)"
table3_dt[term == "DOM_C7",
          c("hr_m1", "hr_m2", "hr_m3") := "1.00 (Ref.)"]

# Phenotype label (合并 ID + 名)
table3_dt[, Phenotype := phenotype_labels[term]]

# 数字格式化
table3_dt[, n_str := format(n_grp, big.mark = ",")]
table3_dt[, events_str := format(events_grp, big.mark = ",")]
table3_dt[, rate_str := sprintf("%.2f", rate_grp)]

# =========================
# 6) 组装最终 7 列
# =========================
final_dt <- table3_dt[, .(
  Phenotype = Phenotype,
  `n` = n_str,
  `Events` = events_str,
  `Rate per 1000 person-years` = rate_str,
  `Model 1 HR (95% CI)` = hr_m1,
  `Model 2 HR (95% CI)` = hr_m2,
  `Model 3 HR (95% CI)` = hr_m3
)]

cat("\n==> Step 4) 最终表格预览:\n")
print(final_dt)

# =========================
# 7) 写入 xlsx
# =========================
cat("\n==> Step 5) 写入 xlsx...\n")
out_xlsx <- file.path(OUT_DIR, "Table3_cox_hazard_ratios.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Table 3")

# 标题
title_text <- "Table 3. Hazard ratios for 4-year all-cause mortality by disease community phenotype"
writeData(wb, "Table 3", title_text, startRow = 1, startCol = 1)
addStyle(wb, "Table 3", createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)

# Header + 数据
writeData(wb, "Table 3", final_dt, startRow = 2,
          headerStyle = createStyle(textDecoration = "bold",
                                    border = "TopBottom",
                                    borderStyle = "thin",
                                    halign = "center",
                                    wrapText = TRUE))

# 列对齐: 1 列 left, 其他 center
left_style <- createStyle(halign = "left", valign = "top")
ctr_style  <- createStyle(halign = "center", valign = "top")
n_data_rows <- nrow(final_dt)
addStyle(wb, "Table 3", left_style,
         rows = 3:(2 + n_data_rows),
         cols = 1, gridExpand = TRUE)
addStyle(wb, "Table 3", ctr_style,
         rows = 3:(2 + n_data_rows),
         cols = 2:7, gridExpand = TRUE)

# Footnote
foot_row <- 3 + n_data_rows + 1
foot_lines <- c(
  "Reference: C7 (Upper-airway and orodental inflammation; n = 47,450, events = 1,362).",
  "Model 1 was unadjusted. Model 2 was adjusted for age and sex. Model 3 was Model 2 plus log-transformed total number of baseline diagnoses.",
  "n indicates the number of participants in each phenotype; Events indicates the number of all-cause deaths during follow-up.",
  "MIXED denotes participants without a single dominant community phenotype (no community share \u2265 0.40 at baseline)."
)
for (i in seq_along(foot_lines)) {
  writeData(wb, "Table 3", foot_lines[i],
            startRow = foot_row + i - 1, startCol = 1)
  addStyle(wb, "Table 3",
           createStyle(fontSize = 9, textDecoration = "italic",
                       wrapText = TRUE),
           rows = foot_row + i - 1, cols = 1)
}

# 列宽
setColWidths(wb, "Table 3", cols = 1, widths = 50)
setColWidths(wb, "Table 3", cols = 2, widths = 12)
setColWidths(wb, "Table 3", cols = 3, widths = 12)
setColWidths(wb, "Table 3", cols = 4, widths = 18)
setColWidths(wb, "Table 3", cols = 5:7, widths = 22)

# 保存
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n   \u2713 已保存: %s\n", out_xlsx))

# =========================
# 8) Console 报告
# =========================
cat("\n========================================================\n")
cat("Table 3 \u2014 Decision rationale:\n")
cat("========================================================\n")
cat("\nPaper Model 1 = Cox data M1 (Unadjusted)\n")
cat("Paper Model 2 = Cox data M2 (+ age + sex)\n")
cat("Paper Model 3 = Cox data M4 (+ age + sex + log baseline diagnoses count)\n")
cat("\n*** M3 (community burden score) is OMITTED from Table 3 ***\n")
cat("    Rationale: total_score_unw is derived from dominant community,\n")
cat("    introducing collider/proxy bias when used as a covariate.\n")
cat("    M4 (log baseline diagnoses count) provides cleaner adjustment.\n")
cat("\nPaper Results text describes M2 \u2192 M4 (=paper M3) increment:\n")
cat("    \u0394C = 0.0024 (95%% CI 0.0021\u20130.0027), \"modest but statistically robust\"\n")
cat("\n========================================================\n")
cat("Final Table 3 (17 rows):\n")
cat("========================================================\n")
for (i in seq_len(nrow(final_dt))) {
  cat(sprintf("\n%2d. %s\n", i, final_dt$Phenotype[i]))
  cat(sprintf("    n = %s, Events = %s, Rate = %s /1000py\n",
              final_dt$n[i], final_dt$Events[i],
              final_dt$`Rate per 1000 person-years`[i]))
  cat(sprintf("    M1: %s\n", final_dt$`Model 1 HR (95% CI)`[i]))
  cat(sprintf("    M2: %s\n", final_dt$`Model 2 HR (95% CI)`[i]))
  cat(sprintf("    M3: %s\n", final_dt$`Model 3 HR (95% CI)`[i]))
}

cat(sprintf("\n\u4fdd\u5b58\u8def\u5f84: %s\n", out_xlsx))