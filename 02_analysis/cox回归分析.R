#!/usr/bin/env Rscript
# ==============================================================================
# Step 5: Cox Proportional Hazards Regression — Dominant Community -> Mortality
#
# Main analysis (PRIMARY):
#   Exposure:   dom_community_unw (16 DOM + MIXED, MIXED merged per docx §2.4.3.5)
#   Dataset:    NONE excluded (n=372, <0.1% of cohort)
#   Reference:  AUTO-SELECTED — lowest-mortality DOM community with
#               N >= 15,000 and 35-65% male proportion
#               (at run time this selects DOM_C7 = upper airway & orodental inflammation)
#   Outcome:    all-cause mortality (time_years + event)
#   Landmark:   index = 2021-01-01; max follow-up = 4 years
#
#   Model 1 (unadjusted):  dom_community_unw
#   Model 2 (PRIMARY):     + age_at_index + sex_binary
#   Model 3 (burden ctrl): + total_score_unw  (continuous comorbidity gradient)
#
# Sensitivity analyses:
#   S1: Grace censoring (alternative event/time definitions)
#   S2: TAU = 0.35 (looser dominant-community threshold)
#   S3: TAU = 0.45 (stricter dominant-community threshold)
#   S4: Weighted community scores (dom_community_w)
#   S5: DOM requires baseline_unique_codes >= 2 (restrict rare-code patients)
#   S6: Include NONE (not excluded)
#   S7: Merge NONE into MIXED
#   S8: Weighted + TAU = 0.40
#   S9: Alternative reference = second-lowest-mortality DOM
#
# Outputs:
#   cox_regression_results.xlsx  (Main_results, Sensitivity_*, Model comparison)
# ==============================================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  pkgs <- c("data.table", "survival", "ggplot2", "openxlsx", "scales")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  }
  library(data.table)
  library(survival)
  library(ggplot2)
  library(openxlsx)
  library(scales)
})

options(stringsAsFactors = FALSE)

# =========================
# 0) Paths & Parameters
# =========================
INPUT_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR   <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# *** PRIMARY SETTINGS ***
# Reference group will be auto-selected after loading data (lowest-mortality DOM with n>=10000)
REF_PRIMARY   <- NULL                      # auto-selected below
GROUP_COL     <- "dom_community_unw"       # PRIMARY: merged MIXED (per docx §2.4.3.5)
EXCLUDE_NONE  <- TRUE                      # NONE (n~372) excluded from main analysis

TAU_MAIN <- 0.40
FONT_FAMILY <- "sans"
RUN_CONTINUOUS <- FALSE
MIN_EVENTS_DOM <- 5

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

DOM_PALETTE <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#33A02C", "#FB9A99", "#FDBF6F",
  "#CAB2D6", "#6A3D9A", "#B15928", "#8DD3C7", "#BEBADA"
)
GROUP_COLORS_FIXED <- c(
  "MIXED"      = "#377EB8",
  "NONE"       = "#999999",
  "DOM_RARE"   = "#CCCCCC"
)

# =========================
# 1) Load data + QC
# =========================
cat("==> Step 1) Loading data...\n")
dt_full <- fread(INPUT_CSV)
cat("   Raw rows:", nrow(dt_full), "  Events:", sum(dt_full$event == 1), "\n")

required_cols <- c("time_years", "event", "dom_community_unw", "dom_community_w",
                   "total_score_unw", "age_at_index", "sex_binary", "baseline_unique_codes")
# Only check essential ones (split_w may not exist)
essential <- c("time_years", "event", GROUP_COL, "total_score_unw",
               "age_at_index", "sex_binary", "baseline_unique_codes")
missing_cols <- setdiff(essential, names(dt_full))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

dt_full[, event := as.integer(event == 1L)]
dt_full[, sex_binary := as.integer(sex_binary == 1L)]
dt_full <- dt_full[time_years >= 0 & !is.na(time_years) & !is.na(event)]
n_zero <- sum(dt_full$time_years == 0)
if (n_zero > 0) {
  cat(sprintf("   time_years == 0: %d — replaced with 0.5/365.25\n", n_zero))
  dt_full[time_years == 0, time_years := 0.5 / 365.25]
}

# Grace event
has_grace <- all(c("censored_by_grace", "time_years_grace") %in% names(dt_full))
if (has_grace) {
  dt_full[is.na(censored_by_grace), censored_by_grace := 0L]
  dt_full[, event_grace := fifelse(censored_by_grace == 1L, 0L, as.integer(event))]
}

dt_full[, baseline_codes_log := log(baseline_unique_codes + 1)]

# =========================
# 1b) Merge rare DOM groups (in split column)
# =========================
cat("==> Step 1b) Checking DOM group sizes...\n")

merge_rare_dom <- function(dt_in, col) {
  if (!col %in% names(dt_in)) return(dt_in)
  dom_tab <- dt_in[grepl("^DOM_", get(col)), .(events = sum(event)), by = c(col)]
  setnames(dom_tab, col, "grp")
  rare <- dom_tab[events < MIN_EVENTS_DOM, grp]
  if (length(rare) > 0) {
    n_rare <- dt_in[get(col) %in% rare, .N]
    cat(sprintf("   [%s] Merging %d rare DOM into DOM_RARE: %s (n=%d)\n",
                col, length(rare), paste(rare, collapse = ", "), n_rare))
    dt_in[get(col) %in% rare, (col) := "DOM_RARE"]
  }
  dt_in
}

dt_full <- merge_rare_dom(dt_full, GROUP_COL)
if ("dom_community_unw" %in% names(dt_full))    dt_full <- merge_rare_dom(dt_full, "dom_community_unw")
if ("dom_community_w" %in% names(dt_full))       dt_full <- merge_rare_dom(dt_full, "dom_community_w")
if ("dom_community_w" %in% names(dt_full)) dt_full <- merge_rare_dom(dt_full, "dom_community_w")

# Build color map
all_dom <- sort(unique(grep("^DOM_C", dt_full[[GROUP_COL]], value = TRUE)))
dom_color_map <- setNames(DOM_PALETTE[seq_along(all_dom)], all_dom)
ALL_COLORS <- c(GROUP_COLORS_FIXED, dom_color_map)

# =========================
# 1c) Create main analysis dataset (exclude NONE)
# =========================
n_none <- dt_full[get(GROUP_COL) == "NONE", .N]
n_none_ev <- dt_full[get(GROUP_COL) == "NONE", sum(event)]
cat(sprintf("   NONE group: n=%d, events=%d — %s from main analysis\n",
            n_none, n_none_ev,
            ifelse(EXCLUDE_NONE, "EXCLUDED", "INCLUDED")))

if (EXCLUDE_NONE) {
  dt <- dt_full[get(GROUP_COL) != "NONE"]
} else {
  dt <- copy(dt_full)
}

# =========================
# 1d) Auto-select reference group
# =========================
# Strategy (enhanced): among DOM groups with n >= 15000 AND sex balance (35-65% male),
# pick the one with lowest event rate. Sex balance ensures the reference is interpretable
# in both male and female subgroups. Falls back to n >= 10000 + lowest rate if no candidate.
if (is.null(REF_PRIMARY)) {
  cat("==> Auto-selecting reference group (with sex balance constraint)...\n")
  dom_rates <- dt[grepl("^DOM_", get(GROUP_COL)),
                  .(n = .N, events = sum(event),
                    rate_1000py = 1000 * sum(event) / sum(time_years),
                    pct_male = 100 * mean(sex_binary == 1L)),
                  by = c(GROUP_COL)]
  setnames(dom_rates, GROUP_COL, "grp")
  setorder(dom_rates, rate_1000py)
  cat("   DOM groups ranked by mortality rate:\n")
  for (i in 1:nrow(dom_rates)) {
    eligible_tag <- ""
    if (dom_rates$n[i] >= 15000 &&
        dom_rates$pct_male[i] >= 35 && dom_rates$pct_male[i] <= 65) {
      eligible_tag <- "  [eligible: balanced]"
    } else if (dom_rates$n[i] >= 10000) {
      eligible_tag <- "  [fallback only]"
    }
    cat(sprintf("      %-12s  n=%-7d events=%-6d rate=%.2f/1000py  male=%.1f%%%s\n",
                dom_rates$grp[i], dom_rates$n[i], dom_rates$events[i],
                dom_rates$rate_1000py[i], dom_rates$pct_male[i], eligible_tag))
  }
  # Primary criterion: n>=15000 + 35-65% male + lowest rate
  eligible <- dom_rates[n >= 15000 & pct_male >= 35 & pct_male <= 65]
  if (nrow(eligible) == 0) {
    cat("   WARNING: no sex-balanced candidate; falling back to n>=10000 + lowest rate\n")
    eligible <- dom_rates[n >= 10000]
  }
  if (nrow(eligible) == 0) {
    stop("No DOM group meets reference selection criteria.")
  }
  REF_PRIMARY <- eligible$grp[1]
  cat(sprintf("   ==> Auto-selected reference: %s (n=%d, rate=%.2f/1000py, male=%.1f%%)\n\n",
              REF_PRIMARY, eligible$n[1], eligible$rate_1000py[1], eligible$pct_male[1]))
}

# Verify reference group
if (!REF_PRIMARY %in% unique(dt[[GROUP_COL]])) {
  stop(sprintf("Reference '%s' not found in %s.\nAvailable: %s",
               REF_PRIMARY, GROUP_COL,
               paste(sort(unique(dt[[GROUP_COL]])), collapse = ", ")), call. = FALSE)
}

cat(sprintf("   Main analysis: n=%d, events=%d, ref=%s\n",
            nrow(dt), sum(dt$event), REF_PRIMARY))
cat(sprintf("   Groups in main analysis: %s\n",
            paste(sort(unique(dt[[GROUP_COL]])), collapse = ", ")))

# Print distribution
cat("\n   Group distribution (main analysis dataset):\n")
grp_tab <- dt[, .(N = .N, events = sum(event),
                  rate = round(1000 * sum(event) / sum(time_years), 2)),
              by = c(GROUP_COL)]
setnames(grp_tab, GROUP_COL, "group")
setorder(grp_tab, -rate)
for (i in 1:nrow(grp_tab)) {
  cat(sprintf("      %-15s  n=%-8d events=%-6d rate=%.2f/1000py%s\n",
              grp_tab$group[i], grp_tab$N[i], grp_tab$events[i], grp_tab$rate[i],
              ifelse(grp_tab$group[i] == REF_PRIMARY, "  <-- REF", "")))
}
cat("\n")

# =========================
# 2) Helper: Re-assign dominant community with different TAU
# =========================
reassign_dom <- function(dt_in, score_cols, tau, suffix_out, merge_mixed = FALSE) {
  nums_ord <- as.integer(gsub("comm(\\d+)_(unw|w)", "\\1", score_cols))
  score_cols <- score_cols[order(nums_ord)]
  
  S <- as.matrix(dt_in[, ..score_cols])
  S[is.na(S)] <- 0
  
  comm_nums   <- as.integer(gsub("comm(\\d+)_(unw|w)", "\\1", score_cols))
  comm_labels <- paste0("C", comm_nums)
  
  row_sum    <- rowSums(S)
  row_max    <- apply(S, 1, max)
  row_which  <- max.col(S, ties.method = "first")
  top1_share <- ifelse(row_sum > 0, row_max / row_sum, 0)
  dom_label  <- comm_labels[row_which]
  
  if (merge_mixed) {
    group <- ifelse(
      row_sum == 0, "NONE",
      ifelse(top1_share >= tau, paste0("DOM_", dom_label), "MIXED")
    )
  } else {
    group <- ifelse(
      row_sum == 0, "NONE",
      ifelse(top1_share >= tau, paste0("DOM_", dom_label), "MIXED")
    )
  }
  
  grp_tab_new <- data.table(group = group, event = dt_in$event)
  rare_dom <- grp_tab_new[grepl("^DOM_", group),
                          .(ev = sum(event)), by = group][ev < MIN_EVENTS_DOM, group]
  if (length(rare_dom) > 0) group[group %in% rare_dom] <- "DOM_RARE"
  
  dt_in[, (suffix_out) := group]
  dt_in
}

# =========================
# 3) Table 1
# =========================
cat("==> Step 3) Table 1...\n")

make_table1 <- function(dt_in, group_col) {
  dt_in[, grp := get(group_col)]
  
  calc_stats <- function(d, label) {
    n_ev   <- sum(d$event == 1, na.rm = TRUE)
    py_tot <- sum(d$time_years, na.rm = TRUE)
    fu_all_mean <- mean(d$time_years, na.rm = TRUE)
    fu_all_sd   <- sd(d$time_years, na.rm = TRUE)
    fu_range    <- range(d$time_years, na.rm = TRUE)
    
    d_events <- d[event == 1]
    fu_event_med <- if (nrow(d_events) > 0) {
      sprintf("%.2f (%.2f-%.2f)",
              median(d_events$time_years, na.rm = TRUE),
              quantile(d_events$time_years, 0.25, na.rm = TRUE),
              quantile(d_events$time_years, 0.75, na.rm = TRUE))
    } else "—"
    
    ts_col <- if ("total_score_unw" %in% names(d)) d$total_score_unw else rep(NA, nrow(d))
    
    data.table(
      group             = label,
      n                 = nrow(d),
      age_mean_sd       = sprintf("%.1f (%.1f)", mean(d$age_at_index, na.rm = TRUE),
                                  sd(d$age_at_index, na.rm = TRUE)),
      age_lt18_pct      = sprintf("%.1f", 100 * mean(d$age_at_index < 18, na.rm = TRUE)),
      age_18_39_pct     = sprintf("%.1f", 100 * mean(d$age_at_index >= 18 & d$age_at_index < 40, na.rm = TRUE)),
      age_40_59_pct     = sprintf("%.1f", 100 * mean(d$age_at_index >= 40 & d$age_at_index < 60, na.rm = TRUE)),
      age_60_74_pct     = sprintf("%.1f", 100 * mean(d$age_at_index >= 60 & d$age_at_index < 75, na.rm = TRUE)),
      age_75plus_pct    = sprintf("%.1f", 100 * mean(d$age_at_index >= 75, na.rm = TRUE)),
      male_pct          = sprintf("%.1f", 100 * mean(d$sex_binary == 1, na.rm = TRUE)),
      codes_median_iqr  = sprintf("%.0f (%.0f-%.0f)",
                                  median(d$baseline_unique_codes, na.rm = TRUE),
                                  quantile(d$baseline_unique_codes, 0.25, na.rm = TRUE),
                                  quantile(d$baseline_unique_codes, 0.75, na.rm = TRUE)),
      total_score_median_iqr = sprintf("%.4f (%.4f-%.4f)",
                                       median(ts_col, na.rm = TRUE),
                                       quantile(ts_col, 0.25, na.rm = TRUE),
                                       quantile(ts_col, 0.75, na.rm = TRUE)),
      person_years      = sprintf("%.0f", py_tot),
      fu_mean_sd        = sprintf("%.2f (%.2f)", fu_all_mean, fu_all_sd),
      fu_range          = sprintf("%.3f-%.2f", fu_range[1], fu_range[2]),
      fu_death_median_iqr = fu_event_med,
      n_events          = n_ev,
      event_rate_1000py = sprintf("%.2f", 1000 * n_ev / py_tot)
    )
  }
  
  groups <- sort(unique(dt_in$grp))
  dom_grps <- groups[grepl("^DOM_C", groups)]
  other_grps <- setdiff(groups, dom_grps)
  
  dom_rates <- dt_in[grp %in% dom_grps,
                     .(rate = sum(event) / sum(time_years)), by = grp]
  setorder(dom_rates, -rate)
  
  ordered_grps <- c(intersect(c("NONE", "MIXED", "DOM_RARE"), other_grps),
                    dom_rates$grp)
  
  rows <- list(calc_stats(dt_in, "Overall"))
  for (g in ordered_grps) rows[[length(rows) + 1]] <- calc_stats(dt_in[grp == g], g)
  dt_in[, grp := NULL]
  rbindlist(rows)
}

# Table 1 for main analysis (no NONE)
table1 <- make_table1(dt, GROUP_COL)

# Table 1 for full dataset (including NONE, supplementary)
table1_full <- make_table1(dt_full, GROUP_COL)

# MIXED quartile gradient
cat("==> Step 3b) MIXED quartile mortality gradient...\n")
mixed_dt <- dt[get(GROUP_COL) == "MIXED"]
if (nrow(mixed_dt) > 0 && "total_score_unw" %in% names(mixed_dt)) {
  mixed_dt[, ts_quartile := cut(total_score_unw,
                                breaks = quantile(total_score_unw, probs = c(0, 0.25, 0.5, 0.75, 1),
                                                  na.rm = TRUE),
                                labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                                include.lowest = TRUE)]
  mixed_gradient <- mixed_dt[, .(
    n = .N, events = sum(event),
    person_years = round(sum(time_years)),
    rate_1000py = round(1000 * sum(event) / sum(time_years), 2),
    age_mean = round(mean(age_at_index, na.rm = TRUE), 1),
    male_pct = round(100 * mean(sex_binary == 1, na.rm = TRUE), 1),
    score_median = round(median(total_score_unw, na.rm = TRUE), 4)
  ), by = ts_quartile]
  setorder(mixed_gradient, ts_quartile)
  cat("   MIXED quartile gradient:\n"); print(mixed_gradient)
} else {
  mixed_gradient <- data.table(note = "No MIXED patients")
}
rm(mixed_dt)

# =========================
# 4) Core Cox fitting function
# =========================
fit_cox_models <- function(dt_in, group_col, time_col, event_col,
                           model_label, ref_level = REF_PRIMARY,
                           include_total_score = TRUE) {
  
  dt_in[, grp := factor(get(group_col))]
  
  if (ref_level %in% levels(dt_in$grp)) {
    dt_in[, grp := relevel(grp, ref = ref_level)]
  } else {
    fallback <- dt_in[, .N, by = grp][order(-N)][1, as.character(grp)]
    cat(sprintf("   NOTE [%s]: ref '%s' not found. Using '%s' (n=%d)\n",
                model_label, ref_level, fallback, dt_in[grp == fallback, .N]))
    dt_in[, grp := relevel(grp, ref = fallback)]
    ref_level <- fallback
  }
  
  grp_tab <- dt_in[, .(n = .N, events = sum(get(event_col) == 1)), by = grp]
  small_grps <- grp_tab[events < 20]
  if (nrow(small_grps) > 0) {
    cat(sprintf("   WARNING [%s]: groups <20 events: %s\n",
                model_label,
                paste(sprintf("%s(n=%d,ev=%d)", small_grps$grp, small_grps$n, small_grps$events),
                      collapse = ", ")))
  }
  
  results <- list()
  surv_obj <- Surv(dt_in[[time_col]], dt_in[[event_col]])
  m1 <- NULL; m2 <- NULL; m3 <- NULL; m4 <- NULL
  
  extract_grp <- function(s, model_name, m_obj) {
    coef_names <- rownames(s$conf.int)
    grp_idx <- grep("^grp", coef_names)
    if (length(grp_idx) == 0) return(NULL)
    r <- as.data.table(s$conf.int[grp_idx, , drop = FALSE], keep.rownames = "term")
    setnames(r, c("term", "HR", "neg_HR", "HR_lower", "HR_upper"))
    r[, neg_HR := NULL]
    p <- as.data.table(s$coefficients[grp_idx, , drop = FALSE], keep.rownames = "term")
    r[, p_value := p[["Pr(>|z|)"]]]
    r[, `:=`(model = model_name, analysis = model_label, ref_group = ref_level)]
    r[, term := gsub("^grp", "", term)]
    r[, concordance := round(s$concordance[1], 4)]
    r[, AIC := round(AIC(m_obj), 1)]
    r[, n := s$n]; r[, n_events := s$nevent]
    r
  }
  
  m1 <- tryCatch(coxph(surv_obj ~ grp, data = dt_in), error = function(e) NULL)
  if (!is.null(m1)) results[["m1"]] <- extract_grp(summary(m1), "Model 1: Unadjusted", m1)
  
  m2 <- tryCatch(coxph(surv_obj ~ grp + age_at_index + sex_binary, data = dt_in),
                 error = function(e) NULL)
  if (!is.null(m2)) results[["m2"]] <- extract_grp(summary(m2), "Model 2: + age + sex", m2)
  
  if (include_total_score && "total_score_unw" %in% names(dt_in)) {
    m3 <- tryCatch(coxph(surv_obj ~ grp + age_at_index + sex_binary + total_score_unw,
                         data = dt_in), error = function(e) NULL)
    if (!is.null(m3)) results[["m3"]] <- extract_grp(summary(m3),
                                                     "Model 3: + age + sex + total_score", m3)
  }
  
  m4 <- tryCatch(coxph(surv_obj ~ grp + age_at_index + sex_binary + baseline_codes_log,
                       data = dt_in), error = function(e) NULL)
  if (!is.null(m4)) results[["m4"]] <- extract_grp(summary(m4),
                                                   "Model 4: + age + sex + codes", m4)
  
  dt_in[, grp := NULL]
  list(table = rbindlist(results, fill = TRUE),
       models = list(m1 = m1, m2 = m2, m3 = m3, m4 = m4))
}

# =========================
# 5) Main analysis (merged MIXED, NONE excluded, auto-selected reference)
# =========================
cat(sprintf("==> Step 5) Main Cox regression (%s, ref=%s, NONE excluded)...\n",
            GROUP_COL, REF_PRIMARY))
main_res <- fit_cox_models(dt, GROUP_COL, "time_years", "event",
                           sprintf("Main: split, TAU=0.40, ref=%s, no NONE", REF_PRIMARY),
                           ref_level = REF_PRIMARY)
cat("   Main models fitted.\n")
m2_main <- main_res$table[model == "Model 2: + age + sex"]
setorder(m2_main, -HR)
cat("   Model 2 — top 5 HR:\n")
for (i in 1:min(5, nrow(m2_main))) {
  cat(sprintf("      %-15s HR=%.2f (%.2f-%.2f) p=%s\n",
              m2_main$term[i], m2_main$HR[i], m2_main$HR_lower[i], m2_main$HR_upper[i],
              formatC(m2_main$p_value[i], format = "e", digits = 2)))
}

# =========================
# 6) PH assumption test
# =========================
cat("==> Step 6) Proportional hazards diagnostics...\n")
ph_test <- NULL
if (!is.null(main_res$models$m2)) {
  dt[, grp := factor(get(GROUP_COL))]
  dt[, grp := relevel(grp, ref = REF_PRIMARY)]
  
  ph_raw <- tryCatch(cox.zph(main_res$models$m2, terms = FALSE),
                     error = function(e) { cat("   cox.zph FAILED:", e$message, "\n"); NULL })
  ph_raw_joint <- tryCatch(cox.zph(main_res$models$m2, terms = TRUE),
                           error = function(e) NULL)
  
  if (!is.null(ph_raw)) {
    ph_test <- as.data.table(ph_raw$table, keep.rownames = "term")
    setnames(ph_test, c("term", "chisq", "df", "p"))
    ph_test[, term := gsub("^grp", "", term)]
    
    n_ev <- max(1, main_res$models$m2$nevent)
    ph_test[, rho_abs := sqrt(pmin(1, chisq / n_ev))]
    ph_test[, rho := NA_real_]
    
    y_names <- colnames(ph_raw$y)
    y_names_clean <- gsub("^grp", "", y_names)
    if (!is.null(y_names) && !is.null(ph_raw$x)) {
      for (nm in intersect(ph_test$term, y_names_clean)) {
        j <- match(nm, y_names_clean)
        r_sign <- sign(cor(ph_raw$x, ph_raw$y[, j], use = "complete.obs"))
        ph_test[term == nm, rho := rho_abs * r_sign]
      }
      unmatched <- ph_test[is.na(rho) & term != "GLOBAL", term]
      if (length(unmatched) > 0) {
        term_noglob <- ph_test$term[ph_test$term != "GLOBAL"]
        for (ti in seq_along(term_noglob)) {
          nm <- term_noglob[ti]
          if (nm %in% unmatched && ti <= ncol(ph_raw$y)) {
            r_sign <- sign(cor(ph_raw$x, ph_raw$y[, ti], use = "complete.obs"))
            ph_test[term == nm, rho := rho_abs * r_sign]
          }
        }
      }
    }
    ph_test[term == "GLOBAL", rho := NA_real_]
    ph_test[, rho_abs := NULL]
    ph_test[, abs_rho := abs(rho)]
    ph_test[, interpretation := fcase(
      term == "GLOBAL", "",
      is.na(abs_rho), "",
      abs_rho < 0.05, "Negligible (PH holds)",
      abs_rho < 0.10, "Minor departure",
      abs_rho >= 0.10, "Notable departure",
      default = ""
    )]
    ph_test[, abs_rho := NULL]
    ph_test[, chisq := round(chisq, 2)]; ph_test[, rho := round(rho, 4)]; ph_test[, p := signif(p, 4)]
    
    if (!is.null(ph_raw_joint)) {
      joint_tbl <- as.data.table(ph_raw_joint$table, keep.rownames = "term")
      setnames(joint_tbl, c("term", "chisq", "df", "p"))
      grp_joint <- joint_tbl[term == "grp"]
      if (nrow(grp_joint) > 0) {
        grp_joint[, chisq := round(chisq, 2)]; grp_joint[, p := signif(p, 4)]
        grp_joint[, rho := NA_real_]; grp_joint[, interpretation := "(joint test)"]
        grp_joint[, term := "grp (joint)"]
        ph_test <- rbindlist(list(
          ph_test[term != "GLOBAL"], grp_joint, ph_test[term == "GLOBAL"]
        ), fill = TRUE)
      }
    }
    cat("   PH global p:", ph_test[term == "GLOBAL", p], "\n")
    
    # Schoenfeld residual plots
    cat("   Generating Schoenfeld residual plots...\n")
    tryCatch({
      n_sch <- ncol(ph_raw$y); n_per_page <- 3; MAX_PTS <- 5000
      n_total <- length(ph_raw$x)
      sub_idx <- if (n_total > MAX_PTS) { set.seed(42); sort(sample.int(n_total, MAX_PTS)) } else seq_len(n_total)
      t_sub <- ph_raw$x[sub_idx]; y_sub <- ph_raw$y[sub_idx, , drop = FALSE]
      y_colnames_clean <- gsub("^grp", "", colnames(ph_raw$y))
      pdf(file.path(OUT_DIR, "PH_schoenfeld_residuals.pdf"), width = 10, height = 3.5 * n_per_page, onefile = TRUE)
      for (pg in 1:ceiling(n_sch / n_per_page)) {
        idx_s <- (pg - 1) * n_per_page + 1; idx_e <- min(pg * n_per_page, n_sch)
        par(mfrow = c(idx_e - idx_s + 1, 1), mar = c(4, 5, 3, 1))
        for (i in idx_s:idx_e) {
          if (i > ncol(y_sub)) next
          col_t <- y_colnames_clean[i]
          rho_i <- ph_test[term == col_t, rho][1]; if (is.na(rho_i)) rho_i <- 0
          p_i <- ph_test[term == col_t, p][1]; if (is.na(p_i)) p_i <- 1
          interp <- ph_test[term == col_t, interpretation][1]; if (is.na(interp)) interp <- ""
          plot(t_sub, y_sub[, i], pch = 20, cex = 0.3, col = adjustcolor("grey40", 0.3),
               xlab = "Time", ylab = "Schoenfeld residual",
               main = sprintf("%s (rho=%.4f, p=%s, %s)", col_t, rho_i, formatC(p_i, format = "e", digits = 2), interp),
               cex.main = 0.9, font.main = 1)
          lo <- tryCatch(loess(y_sub[, i] ~ t_sub, span = 0.3), error = function(e) NULL)
          if (!is.null(lo)) { ord <- order(t_sub); pred <- predict(lo, se = TRUE)
          lines(t_sub[ord], pred$fit[ord], col = "red", lwd = 2)
          lines(t_sub[ord], (pred$fit + 1.96 * pred$se.fit)[ord], col = "red", lty = 2, lwd = 0.8)
          lines(t_sub[ord], (pred$fit - 1.96 * pred$se.fit)[ord], col = "red", lty = 2, lwd = 0.8) }
          abline(h = 0, lty = 3)
        }
      }
      dev.off(); cat("   Schoenfeld plots saved.\n")
    }, error = function(e) { cat("   Schoenfeld FAILED:", e$message, "\n"); tryCatch(dev.off(), error = function(e2) NULL) })
    
    # Log-log plots
    cat("   Generating log-log plots...\n")
    tryCatch({
      dt[, grp_ph := factor(get(GROUP_COL))]
      sf_ph <- survfit(Surv(time_years, event) ~ grp_ph, data = dt)
      sdata_ph <- data.table(time = sf_ph$time, surv = sf_ph$surv,
                             stratum_label = rep(names(sf_ph$strata), sf_ph$strata))
      sdata_ph[, stratum_label := gsub("grp_ph=", "", stratum_label)]
      sdata_ph <- sdata_ph[surv > 0 & surv < 1 & time > 0]
      sdata_ph[, log_time := log(time)]; sdata_ph[, log_neg_log_surv := log(-log(surv))]
      
      avail <- intersect(unique(sdata_ph$stratum_label), names(ALL_COLORS))
      sdata_ll <- sdata_ph[stratum_label %in% avail]
      if (nrow(sdata_ll) > 0) {
        p_ll <- ggplot(sdata_ll, aes(x = log_time, y = log_neg_log_surv, color = stratum_label)) +
          geom_line(linewidth = 0.5, alpha = 0.7) + scale_color_manual(values = ALL_COLORS[avail]) +
          labs(title = "Log-log plot", x = "log(time)", y = "log(-log(S(t)))", color = "Group") +
          theme_classic(base_size = 11) + theme(text = element_text(family = FONT_FAMILY),
                                                plot.title = element_text(hjust = 0.5, face = "bold"),
                                                legend.position = "right", legend.text = element_text(size = 7),
                                                legend.key.height = unit(0.35, "cm"))
        ggsave(file.path(OUT_DIR, "PH_loglog_all.pdf"), p_ll, width = 14, height = 8, units = "in", device = PDF_DEVICE)
      }
      dt[, grp_ph := NULL]
    }, error = function(e) { cat("   Log-log FAILED:", e$message, "\n")
      if ("grp_ph" %in% names(dt)) dt[, grp_ph := NULL] })
  }
}
if ("grp" %in% names(dt)) dt[, grp := NULL]

# =========================
# 7) KM curves
# =========================
cat("==> Step 7) KM survival curves...\n")

# Summary KM: DOM(any) vs MIXED (merged)
dt[, grp_simple := fcase(
  grepl("^DOM_", get(GROUP_COL)), "DOM (any)",
  get(GROUP_COL) == "MIXED",       "MIXED",
  default = NA_character_
)]
dt[, grp_simple := factor(grp_simple, levels = c("DOM (any)", "MIXED"))]
dt_km <- dt[!is.na(grp_simple)]
sf <- survfit(Surv(time_years, event) ~ grp_simple, data = dt_km)
sdata <- data.table(time = sf$time, surv = sf$surv,
                    stratum_label = rep(names(sf$strata), sf$strata))
sdata[, stratum_label := gsub("grp_simple=", "", stratum_label)]
km_colors <- c("DOM (any)" = "#4DAF4A", "MIXED" = "#377EB8")
y_min <- max(0.85, min(sdata$surv, na.rm = TRUE) - 0.02)
p_km <- ggplot(sdata, aes(x = time, y = surv, color = stratum_label)) +
  geom_step(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(), limits = c(y_min, 1)) +
  scale_color_manual(name = "Group", values = km_colors) +
  labs(title = "Kaplan-Meier: Dominant community vs Mixed multimorbidity",
       subtitle = sprintf("All-cause mortality | n=%d, events=%d | follow-up 2021-2024",
                          nrow(dt_km), sum(dt_km$event)),
       x = "Follow-up (years)", y = "Survival probability") +
  theme_classic(base_size = 12) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
        legend.position = "bottom")
ggsave(file.path(OUT_DIR, "KM_main_unw.pdf"), p_km, width = 9, height = 7, units = "in", device = PDF_DEVICE)
p_km2 <- p_km + scale_y_continuous(labels = percent_format(), limits = c(0, 1))
ggsave(file.path(OUT_DIR, "KM_main_unw_fullrange.pdf"), p_km2, width = 9, height = 7, units = "in", device = PDF_DEVICE)
rm(dt_km)
dt[, grp_simple := NULL]

# Per-community KM (split panels by mortality rate)
dt[, grp := as.character(get(GROUP_COL))]
dom_groups <- unique(dt$grp[grepl("^DOM_", dt$grp)])
rate_dt <- dt[grp %in% dom_groups,
              .(rate = sum(event) / sum(time_years)), by = grp]
setorder(rate_dt, -rate)
n_half <- ceiling(nrow(rate_dt) / 2)
for (pn in list(
  list(grps = rate_dt[1:n_half, grp], label = "Higher mortality", fname = "KM_per_community_high_risk.pdf"),
  list(grps = rate_dt[(n_half+1):nrow(rate_dt), grp], label = "Lower mortality", fname = "KM_per_community_low_risk.pdf")
)) {
  plot_grps <- unique(c(REF_PRIMARY, pn$grps))
  dt_sub <- dt[grp %in% plot_grps]; dt_sub[, grp := factor(grp, levels = plot_grps)]
  sf_sub <- survfit(Surv(time_years, event) ~ grp, data = dt_sub)
  sdata_sub <- data.table(time = sf_sub$time, surv = sf_sub$surv,
                          stratum_label = rep(names(sf_sub$strata), sf_sub$strata))
  sdata_sub[, stratum_label := gsub("grp=", "", stratum_label)]
  pc <- ALL_COLORS[intersect(plot_grps, names(ALL_COLORS))]
  y_min <- max(0.75, min(sdata_sub$surv, na.rm = TRUE) - 0.03)
  p <- ggplot(sdata_sub, aes(x = time, y = surv, color = stratum_label)) +
    geom_step(linewidth = 0.7) +
    scale_y_continuous(labels = percent_format(), limits = c(y_min, 1)) +
    scale_color_manual(name = "Group", values = pc) +
    labs(title = sprintf("KM — %s communities + %s ref", pn$label, REF_PRIMARY),
         x = "Follow-up (years)", y = "Survival probability") +
    theme_minimal(base_size = 11) +
    theme(text = element_text(family = FONT_FAMILY),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right", legend.text = element_text(size = 8))
  ggsave(file.path(OUT_DIR, pn$fname), p, width = 12, height = 7, units = "in", device = PDF_DEVICE)
}
dt[, grp := NULL]
cat("   KM plots saved.\n")

# =========================
# 8) Forest plot
# =========================
cat("==> Step 8) Forest plots...\n")

plot_forest <- function(res_dt, title_text, filename, dt_source = dt,
                        group_col = GROUP_COL,
                        model_filter = "Model 2: + age + sex") {
  m2_dt <- res_dt[model == model_filter]
  if (nrow(m2_dt) == 0) { cat(sprintf("   [Forest] Skipping %s.\n", filename)); return(invisible(NULL)) }
  
  ref_label <- m2_dt$ref_group[1]
  setorder(m2_dt, -HR)
  
  grp_info <- dt_source[, .(n_grp = .N, events_grp = sum(event)), by = c(group_col)]
  setnames(grp_info, group_col, "term")
  m2_dt <- merge(m2_dt, grp_info, by = "term", all.x = TRUE)
  setorder(m2_dt, -HR)
  
  m2_dt[, HR_CI_fmt := sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper)]
  n_terms <- nrow(m2_dt)
  m2_dt[, ypos := n_terms:1]
  
  x_data_min <- max(0.08, min(m2_dt$HR_lower, na.rm = TRUE) * 0.7)
  x_data_max <- min(50, max(m2_dt$HR_upper, na.rm = TRUE) * 1.4)
  x_sep    <- x_data_max * 1.5
  x_col_hr <- x_data_max * 3.0
  x_canvas <- x_data_max * 5.5
  
  hr_breaks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
  hr_breaks <- hr_breaks[hr_breaks >= x_data_min & hr_breaks <= x_data_max]
  if (!1 %in% hr_breaks) hr_breaks <- sort(c(hr_breaks, 1))
  
  text_size_11 <- 11 / 2.845
  text_size_13 <- 13 / 2.845
  
  p <- ggplot(m2_dt, aes(x = HR, y = ypos)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_segment(aes(x = HR_lower, xend = HR_upper, y = ypos, yend = ypos),
                 linewidth = 0.7, color = "black") +
    geom_point(shape = 15, color = "black", size = 3) +
    geom_vline(xintercept = x_sep, color = "grey75", linewidth = 0.3) +
    annotate("text", x = x_col_hr, y = n_terms + 1.0, label = "HR (95% CI)",
             hjust = 0.5, size = text_size_13, fontface = "bold", family = FONT_FAMILY) +
    geom_text(aes(x = x_col_hr, y = ypos, label = HR_CI_fmt),
              hjust = 0.5, size = text_size_11, family = FONT_FAMILY, color = "black") +
    scale_x_log10(limits = c(x_data_min, x_canvas), breaks = hr_breaks,
                  labels = as.character(hr_breaks)) +
    scale_y_continuous(breaks = m2_dt$ypos, labels = m2_dt$term,
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = title_text,
         subtitle = sprintf("Reference: %s | Adjusted for age, sex", ref_label),
         x = "Hazard Ratio (log scale)", y = "") +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(text = element_text(family = FONT_FAMILY, size = 11),
          plot.title = element_text(hjust = 0, face = "bold", size = 13, margin = margin(b = 2)),
          plot.subtitle = element_text(hjust = 0, color = "grey40", size = 11, margin = margin(b = 10)),
          axis.title.x = element_text(size = 13, margin = margin(t = 8)),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          panel.grid.major.y = element_line(color = "grey94", linewidth = 0.3),
          panel.grid.minor = element_blank(), plot.margin = margin(15, 10, 10, 10))
  
  ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 12, units = "in", device = PDF_DEVICE)
  cat(sprintf("   Forest: %s (%d terms, ref=%s)\n", filename, n_terms, ref_label))
  
  table_dt <- m2_dt[, .(Group = term, Events = events_grp, N = n_grp,
                        HR = sprintf("%.2f", HR),
                        `95% CI` = sprintf("%.2f-%.2f", HR_lower, HR_upper),
                        P = fifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)))]
  ref_row <- data.table(Group = ref_label, Events = NA_integer_, N = NA_integer_,
                        HR = "1.00 (ref)", `95% CI` = "—", P = "—")
  rbind(ref_row, table_dt)
}

forest_table <- plot_forest(main_res$table,
                            sprintf("Forest plot: Main analysis (ref=%s)", REF_PRIMARY),
                            "forest_main.pdf")

# =========================
# 9) Sensitivity analyses
# =========================
cat("==> Step 9) Sensitivity analyses...\n")
sens_results <- list()
unw_cols <- grep("^comm\\d+_unw$", names(dt), value = TRUE)

# --- S1: Grace censoring ---
cat("   S1: Grace censoring...\n")
if (has_grace) {
  dt_s1 <- dt[time_years_grace >= 0 & !is.na(time_years_grace)]
  dt_s1[time_years_grace == 0, time_years_grace := 0.5 / 365.25]
  s1 <- fit_cox_models(dt_s1, GROUP_COL, "time_years_grace", "event_grace",
                       "S1: Grace censoring", ref_level = REF_PRIMARY)
  sens_results[["S1_grace"]] <- s1$table; rm(dt_s1)
}

# --- S2: TAU = 0.35 (split, no NONE) ---
cat("   S2: TAU = 0.35...\n")
dt <- reassign_dom(dt, unw_cols, tau = 0.35, "dom_s2", merge_mixed = TRUE)
dt_s2 <- dt[dom_s2 != "NONE"]
s2 <- fit_cox_models(dt_s2, "dom_s2", "time_years", "event",
                     "S2: TAU=0.35", ref_level = REF_PRIMARY)
sens_results[["S2_tau035"]] <- s2$table; rm(dt_s2)

# --- S3: TAU = 0.45 (split, no NONE) ---
cat("   S3: TAU = 0.45...\n")
dt <- reassign_dom(dt, unw_cols, tau = 0.45, "dom_s3", merge_mixed = TRUE)
dt_s3 <- dt[dom_s3 != "NONE"]
s3 <- fit_cox_models(dt_s3, "dom_s3", "time_years", "event",
                     "S3: TAU=0.45", ref_level = REF_PRIMARY)
sens_results[["S3_tau045"]] <- s3$table; rm(dt_s3)

# --- S4: Weighted scores (split, no NONE) ---
cat("   S4: Weighted scores...\n")
if ("dom_community_w" %in% names(dt)) {
  dt_s4 <- dt[dom_community_w != "NONE"]
  ref_w <- if (REF_PRIMARY %in% unique(dt_s4$dom_community_w)) REF_PRIMARY else {
    dt_s4[grepl("^DOM_", dom_community_w), .N, by = dom_community_w][order(-N)][1, dom_community_w]
  }
  s4 <- fit_cox_models(dt_s4, "dom_community_w", "time_years", "event",
                       sprintf("S4: Weighted (ref=%s)", ref_w), ref_level = ref_w)
  sens_results[["S4_weighted"]] <- s4$table
  plot_forest(s4$table, sprintf("Forest: Weighted (ref=%s)", ref_w),
              "forest_S4_weighted.pdf", dt_source = dt_s4, group_col = "dom_community_w")
  rm(dt_s4)
} else {
  cat("   S4 skipped (dom_community_w not found).\n")
}

# --- S5: DOM requires codes >= 2 ---
cat("   S5: DOM requires codes >= 2...\n")
dt[, dom_s5 := as.character(get(GROUP_COL))]
# Reclassify single-code DOM patients as MIXED (they don't meet strict-DOM criterion)
dt[grepl("^DOM_", dom_s5) & baseline_unique_codes < 2, dom_s5 := "MIXED"]
n_reclass <- dt[as.character(get(GROUP_COL)) != dom_s5, .N]
cat(sprintf("      Reclassified %d patients from DOM to MIXED\n", n_reclass))
dt_s5 <- dt[dom_s5 != "NONE"]
s5 <- fit_cox_models(dt_s5, "dom_s5", "time_years", "event",
                     "S5: DOM requires codes>=2", ref_level = REF_PRIMARY)
sens_results[["S5_strict_dom"]] <- s5$table; rm(dt_s5)

# --- S_includeNONE: Include NONE back (full dataset) ---
cat("   S_includeNONE: Include NONE...\n")
s_incNone <- fit_cox_models(dt_full, GROUP_COL, "time_years", "event",
                            "S_includeNONE: full dataset with NONE",
                            ref_level = REF_PRIMARY)
sens_results[["S_includeNONE"]] <- s_incNone$table
plot_forest(s_incNone$table, "Forest: Include NONE",
            "forest_S_includeNONE.pdf", dt_source = dt_full, group_col = GROUP_COL)

# --- S_mergeNONE: Merge NONE into MIXED ---
cat("   S_mergeNONE: Merge NONE into MIXED...\n")
dt_full[, dom_mergeNone := as.character(get(GROUP_COL))]
dt_full[dom_mergeNone == "NONE", dom_mergeNone := "MIXED"]
s_mergeNone <- fit_cox_models(dt_full, "dom_mergeNone", "time_years", "event",
                              "S_mergeNONE: NONE -> MIXED",
                              ref_level = REF_PRIMARY)
sens_results[["S_mergeNONE"]] <- s_mergeNone$table

sens_all <- rbindlist(sens_results, fill = TRUE)
cat("   All sensitivity analyses complete.\n")

# Sensitivity comparison forest (top/bottom 5 DOM, same ref group only)
cat("   Sensitivity comparison forest...\n")
same_ref_sens <- c("S1_grace", "S2_tau035", "S3_tau045", "S5_strict_dom",
                   "S_includeNONE", "S_mergeNONE")
compare_dt <- rbindlist(c(
  list(main_res$table[model == "Model 2: + age + sex"]),
  lapply(sens_results[intersect(same_ref_sens, names(sens_results))],
         function(x) x[model == "Model 2: + age + sex"])
), fill = TRUE)
if (nrow(compare_dt) > 0) {
  compare_dom <- compare_dt[grepl("^DOM_", term)]
  if (nrow(compare_dom) > 0) {
    main_dom <- main_res$table[model == "Model 2: + age + sex" & grepl("^DOM_", term)]
    setorder(main_dom, -HR)
    key_groups <- unique(c(head(main_dom$term, 5), tail(main_dom$term, 5)))
    compare_key <- compare_dom[term %in% key_groups]
    compare_key[, label_full := paste0(term, " [", analysis, "]")]
    compare_key[, label_full := factor(label_full, levels = rev(unique(label_full)))]
    x_min <- max(0.1, min(compare_key$HR_lower, na.rm = TRUE) * 0.8)
    x_max <- min(30, max(compare_key$HR_upper, na.rm = TRUE) * 1.2)
    p_sens <- ggplot(compare_key, aes(x = HR, y = label_full, shape = analysis)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      geom_segment(aes(x = HR_lower, xend = HR_upper, y = label_full, yend = label_full), linewidth = 0.5) +
      geom_point(size = 2.5) +
      scale_shape_manual(values = c(15, 16, 17, 18, 8, 3, 4, 5, 6, 7, 9, 10)) +
      scale_x_log10(limits = c(x_min, x_max * 1.3)) +
      labs(title = sprintf("Sensitivity: Top/Bottom 5 DOM (ref=%s)", REF_PRIMARY),
           x = "Hazard Ratio (log scale)", y = "", shape = "Analysis") +
      theme_classic(base_size = 10) +
      theme(text = element_text(family = FONT_FAMILY),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom", legend.text = element_text(size = 7),
            axis.text.y = element_text(size = 7), axis.ticks.y = element_blank()) +
      guides(shape = guide_legend(ncol = 2))
    ggsave(file.path(OUT_DIR, "forest_sensitivity_comparison.pdf"), p_sens,
           width = 14, height = max(8, nrow(compare_key) * 0.28 + 3), units = "in", device = PDF_DEVICE)
  }
}

# =========================
# 10) Excel output
# =========================
cat("==> Step 10) Writing Excel...\n")
wb <- createWorkbook()
headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                           valign = "center", fgFill = "#4472C4", textDecoration = "bold")
add_sheet <- function(sn, dt_in) {
  sn <- substr(sn, 1, 31)
  addWorksheet(wb, sn); writeDataTable(wb, sn, dt_in, withFilter = TRUE)
  setColWidths(wb, sn, cols = 1:ncol(dt_in), widths = "auto")
  addStyle(wb, sn, headerStyle, rows = 1, cols = 1:ncol(dt_in), gridExpand = TRUE)
}

add_sheet("Table1_main", table1)
add_sheet("Table1_full_with_NONE", table1_full)
add_sheet("MIXED_quartile_gradient", mixed_gradient)
add_sheet("Main_results", main_res$table)
if (!is.null(ph_test)) add_sheet("PH_test", ph_test)
if (!is.null(forest_table)) add_sheet("Forest_table", forest_table)
add_sheet("Sensitivity_all", sens_all)
for (nm in names(sens_results)) add_sheet(nm, sens_results[[nm]])

grp_size <- dt[, .(N = .N, events = sum(event),
                   rate_1000py = round(1000 * sum(event) / sum(time_years), 2)),
               by = c(GROUP_COL)]
setnames(grp_size, GROUP_COL, "group")
setorder(grp_size, -rate_1000py)
add_sheet("Group_sizes_main", grp_size)

xlsx_file <- file.path(OUT_DIR, "cox_regression_results.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", xlsx_file, "\n")

# =========================
# Cleanup
# =========================
temp_cols <- c("dom_s2", "dom_s3", "dom_s5", "dom_mergeNone",
               "baseline_codes_log", "event_grace", "grp")
dt[, (intersect(temp_cols, names(dt))) := NULL]
dt_full[, (intersect(c("dom_mergeNone", "baseline_codes_log", "event_grace"), names(dt_full))) := NULL]

# =========================
# Summary
# =========================
cat("\n============================================================\n")
cat("Cox regression (dominant community) complete.\n")
cat("============================================================\n")
cat("Output:", OUT_DIR, "\n\n")
cat(sprintf("MAIN ANALYSIS:\n"))
cat(sprintf("  Group column:  %s\n", GROUP_COL))
cat(sprintf("  Groups:        16 DOM + MIXED\n"))
cat(sprintf("  NONE:          EXCLUDED (n=%d, events=%d)\n", n_none, n_none_ev))
cat(sprintf("  Reference:     %s (auto-selected low-mortality, sex-balanced DOM)\n", REF_PRIMARY))
cat(sprintf("  Dataset:       n=%d, events=%d\n", nrow(dt), sum(dt$event)))
cat("\nSensitivity analyses:\n")
cat("  S1:             Grace censoring\n")
cat("  S2:             TAU = 0.35\n")
cat("  S3:             TAU = 0.45\n")
cat("  S4:             Weighted community scores\n")
cat("  S5:             DOM requires codes >= 2\n")
cat("  S_includeNONE:  Include NONE back in model\n")
cat("  S_mergeNONE:    NONE merged into MIXED\n")
cat("\nExcel: cox_regression_results.xlsx\n")
cat("  Table1_main (no NONE) + Table1_full (with NONE)\n")
cat("  MIXED_quartile_gradient\n")
cat("  Main_results, PH_test, Forest_table\n")
cat("  Sensitivity_all + individual sheets\n")
cat("\nPDFs:\n")
cat("  forest_main.pdf                    (main analysis, auto-selected ref)\n")
cat("  forest_S4_weighted.pdf\n")
cat("  forest_S_includeNONE.pdf\n")
cat("  forest_sensitivity_comparison.pdf\n")
cat("  KM_main_unw.pdf / _fullrange.pdf\n")
cat("  KM_per_community_high_risk.pdf / _low_risk.pdf\n")
cat("  PH_schoenfeld_residuals.pdf\n")
cat("  PH_loglog_all.pdf\n")
cat("============================================================\n")