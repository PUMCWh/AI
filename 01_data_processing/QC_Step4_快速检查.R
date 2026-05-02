# QC_Step4_快速检查.R
# 读取 Step 4 输出, 打印关键 QC 指标
suppressPackageStartupMessages(library(data.table))

COX_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_dataset_baseline_scores.csv"

dt <- fread(COX_CSV)

cat("\n========== Step 4 QC REPORT ==========\n")
cat(sprintf("Cohort n:                %d\n", nrow(dt)))
cat(sprintf("Events (deaths):         %d (%.2f%%)\n",
            sum(dt$event == 1L), 100 * mean(dt$event == 1L)))
cat(sprintf("Median follow-up:        %.2f years\n", median(dt$time_years)))
cat(sprintf("Person-years total:      %.0f\n", sum(dt$time_years)))
cat(sprintf("Event rate:              %.2f per 1000 py\n",
            1000 * sum(dt$event == 1L) / sum(dt$time_years)))

cat("\n--- Demographics ---\n")
cat(sprintf("Age: mean=%.1f, sd=%.1f, median=%.1f, range=[%.1f, %.1f]\n",
            mean(dt$age_at_index), sd(dt$age_at_index),
            median(dt$age_at_index), min(dt$age_at_index), max(dt$age_at_index)))
cat(sprintf("Male: n=%d (%.2f%%)\n",
            sum(dt$sex_binary == 1L), 100 * mean(dt$sex_binary == 1L)))
cat(sprintf("Female: n=%d (%.2f%%)\n",
            sum(dt$sex_binary == 0L), 100 * mean(dt$sex_binary == 0L)))

cat("\n--- CKM baseline exposure (ITB-safe, 2016-2020 only) ---\n")
ckm_cols <- c("CKM_C_baseline", "CKM_K_baseline", "CKM_M_baseline")
for (cc in ckm_cols) {
  if (cc %in% names(dt)) {
    n <- sum(dt[[cc]] == 1L, na.rm = TRUE)
    cat(sprintf("  %-20s n=%7d (%.2f%%)\n", cc, n, 100 * n / nrow(dt)))
  } else {
    cat(sprintf("  %-20s MISSING!\n", cc))
  }
}

if ("n_ckm_domains_baseline" %in% names(dt)) {
  cat("\n  n_ckm_domains_baseline distribution:\n")
  print(dt[, .N, by = n_ckm_domains_baseline][order(n_ckm_domains_baseline)])
}

cat("\n--- Community scores QC ---\n")
unw_cols <- grep("^comm\\d+_unw$", names(dt), value = TRUE)
cat(sprintf("  Unweighted score columns: %d found\n", length(unw_cols)))
if (length(unw_cols) > 0) {
  score_mat <- as.matrix(dt[, ..unw_cols])
  score_mat[is.na(score_mat)] <- 0
  row_sum <- rowSums(score_mat)
  cat(sprintf("  Patients with any score > 0:  %d (%.2f%%)\n",
              sum(row_sum > 0), 100 * mean(row_sum > 0)))
  cat(sprintf("  Patients with ALL scores = 0: %d (%.2f%%)\n",
              sum(row_sum == 0), 100 * mean(row_sum == 0)))
}

cat("\n--- Role exposure ---\n")
role_cols <- grep("^has_", names(dt), value = TRUE)
for (rc in role_cols) {
  n <- sum(dt[[rc]] == 1L, na.rm = TRUE)
  cat(sprintf("  %-20s n=%7d (%.2f%%)\n", rc, n, 100 * n / nrow(dt)))
}

cat("\n--- Baseline utilization ---\n")
cat(sprintf("  baseline_unique_codes: median=%.0f, IQR=[%.0f, %.0f]\n",
            median(dt$baseline_unique_codes),
            quantile(dt$baseline_unique_codes, 0.25),
            quantile(dt$baseline_unique_codes, 0.75)))

cat("\n======================================\n")