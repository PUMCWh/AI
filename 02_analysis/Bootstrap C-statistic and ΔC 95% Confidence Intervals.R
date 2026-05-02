#!/usr/bin/env Rscript
# ==============================================================================
# Step 12: Bootstrap C-statistic and ΔC 95% Confidence Intervals
#
# Purpose:
#   Produce bootstrap-based 95% confidence intervals for Harrell's C-statistic
#   in Models 1–4 of the main Cox regression, and for the ΔC between successive
#   nested models (M1→M2, M2→M3, M3→M4). These intervals quantify the
#   statistical significance of the incremental prognostic value of each
#   added covariate.
#
# Inputs:
#   - cox_dataset_with_dominant_community.csv (from Step 4.5)
#     Required: time_years, event, dom_community_unw, age_at_index,
#               sex_binary, total_score_unw, baseline_unique_codes
#
# Outputs (written to OUT_DIR):
#   - Step12_bootstrap_C_statistic_results.xlsx
#       Sheet 1: Point_estimates   — original C-stat / ΔC / AIC / BIC
#       Sheet 2: Bootstrap_CI      — 95% CI for each of the above
#       Sheet 3: Delta_pairs       — paired ΔC CI across nested pairs
#       Sheet 4: Run_info          — B, seed, N, events, run time
#
# Method:
#   1) Fit M1-M4 on full cohort; extract point C-statistic and AIC/BIC.
#   2) Draw B stratified bootstrap samples (stratified by event status to
#      preserve event proportion).
#   3) For each bootstrap sample: fit ALL four models, record C-stat.
#   4) Paired ΔC for nested pairs: take bootstrap-sample-level differences
#      (same sample used for both models); percentile CI.
#
# Computational notes:
#   - concordance() is O(n²) worst case but typically fast in the survival pkg
#     (linear-ish). For n = 640K + B=500 + 4 models ≈ 2000 model fits.
#     Expect run time: 30-90 min single-threaded; 10-20 min with 4 workers.
#   - Memory: ~5 GB peak (bootstrap indices + model copies).
#   - Parallel via future.apply for speed.
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(openxlsx)
  library(future.apply)
})

# =========================
# 0) Paths & Parameters
# =========================
COX_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Bootstrap parameters
N_BOOT     <- 500    # number of bootstrap iterations (adjust if too slow; 200 is min acceptable)
SEED       <- 20260423
N_WORKERS  <- 4      # parallel workers (adjust to your machine; set 1 for serial)
REF_GROUP  <- "DOM_C7"   # must match Step 5 reference

# =========================
# 1) Load and prepare Cox dataset
# =========================
cat("==> Step 1) Loading Cox dataset...\n")

if (!file.exists(COX_CSV)) {
  stop("Input not found: ", COX_CSV, "\nRun Steps 4 and 4.5 first.", call. = FALSE)
}

dt <- fread(COX_CSV,
            select = c("pid", "time_years", "event", "age_at_index",
                       "sex_binary", "dom_community_unw",
                       "total_score_unw", "baseline_unique_codes"))

cat("   Loaded:", nrow(dt), "patients\n")

# Drop NONE group (consistent with main Step 5 analysis)
dt <- dt[dom_community_unw != "NONE"]
cat("   After dropping NONE:", nrow(dt), "patients\n")
cat("   Events:", sum(dt$event == 1), sprintf("(%.2f%%)\n",
                                              mean(dt$event == 1) * 100))

# Set reference group for dom_community_unw
dom_levels <- sort(unique(dt$dom_community_unw))
if (!(REF_GROUP %in% dom_levels)) {
  stop("Reference group '", REF_GROUP, "' not found in dom_community_unw.",
       " Available levels: ", paste(dom_levels, collapse = ", "), call. = FALSE)
}
# Put reference first
dom_levels <- c(REF_GROUP, setdiff(dom_levels, REF_GROUP))
dt[, dom_community_unw := factor(dom_community_unw, levels = dom_levels)]

# log-transformed disease count for Model 4
dt[, log_dx_count := log1p(baseline_unique_codes)]

# =========================
# 2) Fit point-estimate models
# =========================
cat("\n==> Step 2) Fitting point-estimate Cox models (M1-M4)...\n")

# Formulas
f_m1 <- Surv(time_years, event) ~ dom_community_unw
f_m2 <- Surv(time_years, event) ~ dom_community_unw + age_at_index + sex_binary
f_m3 <- Surv(time_years, event) ~ dom_community_unw + age_at_index + sex_binary + total_score_unw
f_m4 <- Surv(time_years, event) ~ dom_community_unw + age_at_index + sex_binary + log_dx_count

fit_all <- function(data) {
  list(
    M1 = coxph(f_m1, data = data),
    M2 = coxph(f_m2, data = data),
    M3 = coxph(f_m3, data = data),
    M4 = coxph(f_m4, data = data)
  )
}

pt_fits <- fit_all(dt)

# Extract C-statistic from each
get_C <- function(fit) {
  conc <- survival::concordance(fit)
  as.numeric(conc$concordance)
}

pt_C <- sapply(pt_fits, get_C)
pt_AIC <- sapply(pt_fits, AIC)
pt_BIC <- sapply(pt_fits, BIC)

cat("   Point estimates:\n")
for (m in names(pt_fits)) {
  cat(sprintf("     %s: C = %.4f  |  AIC = %.1f  |  BIC = %.1f\n",
              m, pt_C[m], pt_AIC[m], pt_BIC[m]))
}

# Paired point-estimate ΔC
dC_12 <- pt_C["M2"] - pt_C["M1"]
dC_23 <- pt_C["M3"] - pt_C["M2"]
dC_34 <- pt_C["M4"] - pt_C["M3"]
dC_13 <- pt_C["M3"] - pt_C["M1"]

cat(sprintf("   ΔC (M1→M2) = %.4f\n", dC_12))
cat(sprintf("   ΔC (M2→M3) = %.4f\n", dC_23))
cat(sprintf("   ΔC (M3→M4) = %.4f\n", dC_34))
cat(sprintf("   ΔC (M1→M3) = %.4f\n", dC_13))

# =========================
# 3) Bootstrap setup (stratified by event)
# =========================
cat("\n==> Step 3) Setting up stratified bootstrap (B =", N_BOOT,
    ", workers =", N_WORKERS, ")...\n")

set.seed(SEED)

# Pre-compute row indices for event = 1 vs event = 0 strata
idx_event <- which(dt$event == 1)
idx_ctrl  <- which(dt$event == 0)
n_event <- length(idx_event)
n_ctrl  <- length(idx_ctrl)
cat("   Event stratum:", n_event, " | Control stratum:", n_ctrl, "\n")

# Worker function: one bootstrap iteration
boot_iter <- function(b_seed) {
  set.seed(b_seed)
  samp <- c(
    sample(idx_event, size = n_event, replace = TRUE),
    sample(idx_ctrl,  size = n_ctrl,  replace = TRUE)
  )
  boot_dt <- dt[samp]
  
  # Fit all 4 models; any failure → return NA for that iter
  out <- tryCatch({
    fits <- fit_all(boot_dt)
    c(
      C_M1 = get_C(fits$M1),
      C_M2 = get_C(fits$M2),
      C_M3 = get_C(fits$M3),
      C_M4 = get_C(fits$M4)
    )
  }, error = function(e) {
    warning("Bootstrap iter failed: ", e$message)
    c(C_M1 = NA, C_M2 = NA, C_M3 = NA, C_M4 = NA)
  })
  out
}

boot_seeds <- SEED + seq_len(N_BOOT)

# =========================
# 4) Run bootstrap (parallel)
# =========================
cat("\n==> Step 4) Running", N_BOOT, "bootstrap iterations...\n")
t0 <- Sys.time()

if (N_WORKERS > 1) {
  plan(multisession, workers = N_WORKERS)
  on.exit(plan(sequential), add = TRUE)
  boot_results <- future_sapply(boot_seeds, boot_iter,
                                future.seed = TRUE,
                                simplify = "matrix",
                                future.packages = c("survival", "data.table"),
                                future.globals = list(
                                  fit_all = fit_all, get_C = get_C,
                                  dt = dt, idx_event = idx_event,
                                  idx_ctrl = idx_ctrl,
                                  n_event = n_event, n_ctrl = n_ctrl,
                                  f_m1 = f_m1, f_m2 = f_m2,
                                  f_m3 = f_m3, f_m4 = f_m4
                                ))
} else {
  boot_results <- sapply(boot_seeds, boot_iter, simplify = "matrix")
}

# boot_results is 4 × N_BOOT (rows = C_M1 ... C_M4, columns = iterations)
boot_mat <- t(boot_results)  # now N_BOOT × 4
colnames(boot_mat) <- c("C_M1", "C_M2", "C_M3", "C_M4")
boot_dt_mat <- as.data.table(boot_mat)

n_fail <- sum(is.na(boot_dt_mat$C_M1))
cat(sprintf("   Completed: %d iterations (%d failed)\n",
            sum(!is.na(boot_dt_mat$C_M1)), n_fail))

run_time <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("   Run time: %.1f min\n", as.numeric(run_time)))

# Remove any failed rows
boot_dt_mat <- boot_dt_mat[!is.na(C_M1) & !is.na(C_M2) &
                             !is.na(C_M3) & !is.na(C_M4)]
cat("   Valid bootstrap samples:", nrow(boot_dt_mat), "\n")

# =========================
# 5) Compute percentile CIs
# =========================
cat("\n==> Step 5) Computing 95% percentile CIs...\n")

pct_ci <- function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)

# C-stat CIs
c_ci <- data.table(
  Model       = c("M1", "M2", "M3", "M4"),
  `C-statistic (point)` = round(pt_C, 4),
  `CI lower`  = round(c(pct_ci(boot_dt_mat$C_M1)[1],
                        pct_ci(boot_dt_mat$C_M2)[1],
                        pct_ci(boot_dt_mat$C_M3)[1],
                        pct_ci(boot_dt_mat$C_M4)[1]), 4),
  `CI upper`  = round(c(pct_ci(boot_dt_mat$C_M1)[2],
                        pct_ci(boot_dt_mat$C_M2)[2],
                        pct_ci(boot_dt_mat$C_M3)[2],
                        pct_ci(boot_dt_mat$C_M4)[2]), 4),
  `Boot SD`   = round(c(sd(boot_dt_mat$C_M1, na.rm = TRUE),
                        sd(boot_dt_mat$C_M2, na.rm = TRUE),
                        sd(boot_dt_mat$C_M3, na.rm = TRUE),
                        sd(boot_dt_mat$C_M4, na.rm = TRUE)), 4)
)
c_ci[, `C (95% CI)` := sprintf("%.4f (%.4f\u2013%.4f)",
                               `C-statistic (point)`, `CI lower`, `CI upper`)]

# Paired ΔC CIs (use same-bootstrap-sample differences)
boot_dt_mat[, dC_M1_M2 := C_M2 - C_M1]
boot_dt_mat[, dC_M2_M3 := C_M3 - C_M2]
boot_dt_mat[, dC_M3_M4 := C_M4 - C_M3]
boot_dt_mat[, dC_M1_M3 := C_M3 - C_M1]
boot_dt_mat[, dC_M1_M4 := C_M4 - C_M1]
boot_dt_mat[, dC_M2_M4 := C_M4 - C_M2]

dC_ci <- data.table(
  Comparison           = c("M1 \u2192 M2", "M2 \u2192 M3", "M3 \u2192 M4",
                           "M1 \u2192 M3", "M1 \u2192 M4", "M2 \u2192 M4"),
  `ΔC (point)`         = round(c(dC_12, dC_23, dC_34,
                                 dC_13, pt_C["M4"] - pt_C["M1"],
                                 pt_C["M4"] - pt_C["M2"]), 4),
  `ΔC CI lower`        = round(c(pct_ci(boot_dt_mat$dC_M1_M2)[1],
                                 pct_ci(boot_dt_mat$dC_M2_M3)[1],
                                 pct_ci(boot_dt_mat$dC_M3_M4)[1],
                                 pct_ci(boot_dt_mat$dC_M1_M3)[1],
                                 pct_ci(boot_dt_mat$dC_M1_M4)[1],
                                 pct_ci(boot_dt_mat$dC_M2_M4)[1]), 4),
  `ΔC CI upper`        = round(c(pct_ci(boot_dt_mat$dC_M1_M2)[2],
                                 pct_ci(boot_dt_mat$dC_M2_M3)[2],
                                 pct_ci(boot_dt_mat$dC_M3_M4)[2],
                                 pct_ci(boot_dt_mat$dC_M1_M3)[2],
                                 pct_ci(boot_dt_mat$dC_M1_M4)[2],
                                 pct_ci(boot_dt_mat$dC_M2_M4)[2]), 4)
)
dC_ci[, `Significant at 0.05` := sign(`ΔC CI lower`) == sign(`ΔC CI upper`)]
dC_ci[, `ΔC (95% CI)` := sprintf("%.4f (%.4f\u2013%.4f)",
                                 `ΔC (point)`, `ΔC CI lower`, `ΔC CI upper`)]

cat("   C-statistic CIs:\n")
for (i in seq_len(nrow(c_ci))) {
  cat(sprintf("     %s: %s\n", c_ci$Model[i], c_ci$`C (95% CI)`[i]))
}
cat("   ΔC CIs:\n")
for (i in seq_len(nrow(dC_ci))) {
  cat(sprintf("     %s: %s %s\n", dC_ci$Comparison[i], dC_ci$`ΔC (95% CI)`[i],
              ifelse(dC_ci$`Significant at 0.05`[i], "(sig)", "(ns)")))
}

# =========================
# 6) Write Excel output
# =========================
cat("\n==> Step 6) Writing results to Excel...\n")

# Sheet 1: Point estimates
pt_tbl <- data.table(
  Model        = c("M1", "M2", "M3", "M4"),
  Covariates   = c("dom_community_unw",
                   "+ age + sex",
                   "+ total_score_unw",
                   "+ log1p(baseline_unique_codes)"),
  `C-statistic`= round(pt_C, 4),
  `AIC`        = round(pt_AIC, 1),
  `BIC`        = round(pt_BIC, 1),
  `ΔAIC vs M1` = round(pt_AIC - pt_AIC["M1"], 1),
  `ΔBIC vs M1` = round(pt_BIC - pt_BIC["M1"], 1)
)

# Sheet 4: run info
run_info <- data.table(
  Parameter = c("Bootstrap B", "Random seed base", "Parallel workers",
                "Reference group", "Cohort N", "Events",
                "Valid bootstrap samples", "Run time (minutes)",
                "Stratified by", "CI method"),
  Value     = c(as.character(N_BOOT), as.character(SEED),
                as.character(N_WORKERS), REF_GROUP,
                format(nrow(dt), big.mark = ","),
                format(sum(dt$event == 1), big.mark = ","),
                format(nrow(boot_dt_mat), big.mark = ","),
                sprintf("%.1f", as.numeric(run_time)),
                "event status (event = 1 vs 0)",
                "percentile (2.5% - 97.5%)")
)

wb <- createWorkbook()
hdr_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF",
                         halign = "center", valign = "center",
                         fgFill = "#4472C4", textDecoration = "bold")

add_sheet <- function(name, dt_out) {
  addWorksheet(wb, name)
  writeDataTable(wb, name, dt_out, withFilter = TRUE)
  setColWidths(wb, name, cols = 1:ncol(dt_out), widths = "auto")
  addStyle(wb, name, hdr_style, rows = 1,
           cols = 1:ncol(dt_out), gridExpand = TRUE)
}

add_sheet("Point_estimates", pt_tbl)
add_sheet("Bootstrap_CI",     c_ci)
add_sheet("Delta_pairs",      dC_ci)
add_sheet("Run_info",         run_info)

xlsx_file <- file.path(OUT_DIR, "Step12_bootstrap_C_statistic_results.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", basename(xlsx_file), "\n")

# =========================
# 7) Optional: save raw bootstrap sample for later use
# =========================
rds_file <- file.path(OUT_DIR, "Step12_bootstrap_raw_samples.rds")
saveRDS(list(
  point        = list(C = pt_C, AIC = pt_AIC, BIC = pt_BIC),
  bootstrap    = boot_dt_mat,
  parameters   = list(N_BOOT = N_BOOT, SEED = SEED,
                      REF_GROUP = REF_GROUP, N = nrow(dt),
                      events = sum(dt$event == 1)),
  timestamp    = Sys.time()
), rds_file)
cat("   Saved raw bootstrap samples:", basename(rds_file), "\n")

# =========================
# 8) Console summary for paper
# =========================
cat("\n============================================================\n")
cat("Step 12: Bootstrap C-statistic — complete.\n")
cat("============================================================\n\n")
cat("For direct use in Table 4 of the manuscript:\n\n")
for (i in seq_len(nrow(c_ci))) {
  cat(sprintf("  %s: C = %s\n", c_ci$Model[i], c_ci$`C (95% CI)`[i]))
}
cat("\n")
for (i in seq_len(nrow(dC_ci))) {
  cat(sprintf("  %s: \u0394C = %s%s\n",
              dC_ci$Comparison[i], dC_ci$`ΔC (95% CI)`[i],
              ifelse(dC_ci$`Significant at 0.05`[i],
                     "  (significant)", "")))
}
cat("\nOutputs:\n  ", xlsx_file, "\n")
cat("  ", rds_file, "\n", sep = "")
cat("\n============================================================\n")