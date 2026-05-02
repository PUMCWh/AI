#!/usr/bin/env Rscript
# ==============================================================================
# Step 6: Subgroup interaction analyses
#   (1) Sex heterogeneity (male vs female)
#   (2) Age heterogeneity (<60 / 60-75 / >=75)
#
# Methods:
#   - Subgroup Cox: independent fits on each subgroup
#       Sex subgroups: adjusted for age (no sex term, since stratified by sex)
#       Age subgroups: adjusted for sex only (no age term, since stratified by age)
#   - Interaction tests on full cohort:
#       Sex: coxph(Surv ~ grp * sex_binary + age) — Wald per-DOM, anova LRT global
#       Age: coxph(Surv ~ grp * age_decade + sex)
#         age_decade = age_at_index / 10 (per-decade scaling for interpretable HRs)
#         Subgroup HRs displayed by clinical cuts <60 / 60-75 / >=75
#   - Reference group: DOM_C7 (sex-balanced low-mortality, matches main analysis)
#
# Inputs:  cox_dataset_with_dominant_community.csv
# Output:  subgroup_interaction_analysis.xlsx
#            sheets: Sex_interaction, Age_interaction
# ==============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table); library(survival); library(openxlsx)
})

INPUT_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR   <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_COL    <- "dom_community_unw"
REF_GROUP    <- "DOM_C7"
EXCLUDE_NONE <- TRUE
AGE_CUTS     <- c(-Inf, 60, 75, Inf)
AGE_LABELS   <- c("<60", "60-75", ">=75")

community_themes <- c(
  "DOM_C1"="Functional & inflammatory GI disorders",
  "DOM_C2"="Chronic respiratory disease with infections",
  "DOM_C3"="Dermatologic & venous insufficiency disorders",
  "DOM_C4"="Benign gynecologic disorders",
  "DOM_C5"="Age-related ophthalmologic disorders",
  "DOM_C6"="Cerebrovascular disease with neuropsychiatric sequelae",
  "DOM_C7"="Upper airway & orodental inflammation",
  "DOM_C8"="Hepatobiliary disorders",
  "DOM_C9"="Advanced solid malignancies with cachexia",
  "DOM_C10"="Benign & malignant breast disorders",
  "DOM_C11"="Urinary stones, BPH & UTI",
  "DOM_C12"="Thyroid dysfunction & neoplasms",
  "DOM_C13"="Lymphoid & myeloid malignancies",
  "DOM_C14"="CKD with anemia, electrolyte & gout",
  "DOM_C15"="T2DM with neuropathy & joint disease",
  "DOM_C16"="CAD, heart failure & arrhythmia",
  "MIXED"="Mixed (no dominant community)"
)

# =========================
# 1) Load
# =========================
cat("==> Loading data...\n")
dt <- fread(INPUT_CSV)
cat(sprintf("   Rows: %d\n", nrow(dt)))
if (EXCLUDE_NONE) {
  dt <- dt[get(GROUP_COL) != "NONE"]
  cat(sprintf("   After excluding NONE: %d\n", nrow(dt)))
}

dt[, grp := factor(get(GROUP_COL))]
dt[, grp := relevel(grp, ref = REF_GROUP)]
dt[, age_decade := age_at_index / 10]
dt[, age_group := cut(age_at_index, breaks = AGE_CUTS, labels = AGE_LABELS, right = FALSE)]

cat(sprintf("   Sex: M=%d (%.1f%%), F=%d (%.1f%%)\n",
            sum(dt$sex_binary == 1L), 100*mean(dt$sex_binary == 1L),
            sum(dt$sex_binary == 0L), 100*mean(dt$sex_binary == 0L)))
cat("   Age groups:\n"); print(dt[, .(N=.N, events=sum(event)), by=age_group][order(age_group)])

# =========================
# 2) Generic subgroup Cox
# =========================
fit_subgroup <- function(dt_sub, label, formula_rhs) {
  dt_sub <- droplevels(dt_sub)
  cat(sprintf("      %-12s N=%d, events=%d\n", label, nrow(dt_sub), sum(dt_sub$event)))
  
  f <- as.formula(paste("Surv(time_years, event) ~", formula_rhs))
  m <- tryCatch(coxph(f, data = dt_sub),
                error = function(e) { cat("      ERR:", e$message, "\n"); NULL })
  if (is.null(m)) return(NULL)
  
  cf <- summary(m)$coefficients
  ci <- confint(m)
  rows <- grep("^grp", rownames(cf))
  res <- data.table(
    term  = sub("^grp", "", rownames(cf)[rows]),
    hr    = exp(cf[rows, "coef"]),
    lower = exp(ci[rows, 1]),
    upper = exp(ci[rows, 2]),
    p     = cf[rows, "Pr(>|z|)"]
  )
  
  stats <- dt_sub[, .(N=.N, Events=sum(event)), by=grp]
  setnames(stats, "grp", "term"); stats[, term := as.character(term)]
  res <- merge(res, stats, by="term", all.x=TRUE)
  
  ref_st <- dt_sub[grp == REF_GROUP, .(N=.N, Events=sum(event))]
  ref_row <- data.table(term=REF_GROUP, hr=1.00, lower=NA_real_, upper=NA_real_,
                        p=NA_real_, N=ref_st$N, Events=ref_st$Events)
  rbindlist(list(ref_row, res), fill=TRUE)
}

# Robust LRT extractor (handles different anova column names)
get_lrt <- function(m_no, m_int) {
  a <- as.data.frame(anova(m_no, m_int))
  chi_col <- grep("Chisq|Chi", names(a), value=TRUE)[1]
  df_col  <- grep("^Df$|Df", names(a), value=TRUE)[1]
  p_col   <- grep("^P|^Pr", names(a), value=TRUE)[1]
  list(chi = a[2, chi_col], df = a[2, df_col], p = a[2, p_col])
}

# =========================
# 3) Sex interaction
# =========================
cat("\n========== SEX INTERACTION ==========\n")
res_male   <- fit_subgroup(dt[sex_binary == 1], "Male",   "grp + age_at_index")
res_female <- fit_subgroup(dt[sex_binary == 0], "Female", "grp + age_at_index")

cat("   Fitting sex interaction model...\n")
m_sex_no  <- coxph(Surv(time_years, event) ~ grp + age_at_index + sex_binary, data=dt)
m_sex_int <- coxph(Surv(time_years, event) ~ grp * sex_binary + age_at_index, data=dt)
lrt_sex <- get_lrt(m_sex_no, m_sex_int)
cat(sprintf("   Global LRT: chi2=%.2f, df=%d, P=%.3g\n", lrt_sex$chi, lrt_sex$df, lrt_sex$p))

cf_si <- summary(m_sex_int)$coefficients
ir <- grep(":sex_binary$", rownames(cf_si))
sex_p_int <- data.table(
  term = sub(":sex_binary$", "", sub("^grp", "", rownames(cf_si)[ir])),
  p_interaction = cf_si[ir, "Pr(>|z|)"]
)

# =========================
# 4) Age interaction (continuous, per-decade)
# =========================
cat("\n========== AGE INTERACTION ==========\n")
res_age <- list()
for (lvl in AGE_LABELS) {
  res_age[[lvl]] <- fit_subgroup(dt[age_group == lvl], paste0("Age ", lvl), "grp + sex_binary")
}

cat("   Fitting age interaction model (per-decade)...\n")
m_age_no  <- coxph(Surv(time_years, event) ~ grp + age_decade + sex_binary, data=dt)
m_age_int <- coxph(Surv(time_years, event) ~ grp * age_decade + sex_binary, data=dt)
lrt_age <- get_lrt(m_age_no, m_age_int)
cat(sprintf("   Global LRT: chi2=%.2f, df=%d, P=%.3g\n", lrt_age$chi, lrt_age$df, lrt_age$p))

cf_ai <- summary(m_age_int)$coefficients
ar <- grep(":age_decade$", rownames(cf_ai))
age_p_int <- data.table(
  term = sub(":age_decade$", "", sub("^grp", "", rownames(cf_ai)[ar])),
  p_interaction = cf_ai[ar, "Pr(>|z|)"]
)

# =========================
# 5) Format helpers
# =========================
fmt_hr <- function(hr, lo, hi, is_ref) {
  fifelse(is_ref | is.na(hr), "1.00 (ref)",
          sprintf("%.2f (%.2f-%.2f)", hr, lo, hi))
}
fmt_p <- function(p) {
  fifelse(is.na(p), "—",
          fifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# =========================
# 6) Build output tables
# =========================
build_sex_table <- function() {
  m <- copy(res_male);   setnames(m, c("hr","lower","upper","p","N","Events"),
                                  c("HR_M","Lo_M","Hi_M","P_M","N_M","Ev_M"))
  f <- copy(res_female); setnames(f, c("hr","lower","upper","p","N","Events"),
                                  c("HR_F","Lo_F","Hi_F","P_F","N_F","Ev_F"))
  tab <- merge(m, f, by="term", all=TRUE)
  tab <- merge(tab, sex_p_int, by="term", all.x=TRUE)
  tab[term == REF_GROUP, p_interaction := NA]
  
  tab[, theme := community_themes[term]]
  tab[is.na(theme), theme := term]
  tab[, term_disp := sub("^DOM_", "", term)]
  is_ref <- tab$term == REF_GROUP
  
  data.table(
    Group = sprintf("%s — %s", tab$term_disp, tab$theme),
    `Male N`            = tab$N_M,
    `Male Events`       = tab$Ev_M,
    `Male HR (95% CI)`  = fmt_hr(tab$HR_M, tab$Lo_M, tab$Hi_M, is_ref),
    `Male P`            = fmt_p(tab$P_M),
    `Female N`          = tab$N_F,
    `Female Events`     = tab$Ev_F,
    `Female HR (95% CI)`= fmt_hr(tab$HR_F, tab$Lo_F, tab$Hi_F, is_ref),
    `Female P`          = fmt_p(tab$P_F),
    P_interaction       = fmt_p(tab$p_interaction)
  )
}

build_age_table <- function() {
  cols <- c("hr","lower","upper","p","N","Events")
  a <- list()
  for (i in 1:3) {
    x <- copy(res_age[[i]])
    setnames(x, cols, paste0(cols, "_", i))
    a[[i]] <- x
  }
  tab <- Reduce(function(x,y) merge(x, y, by="term", all=TRUE), a)
  tab <- merge(tab, age_p_int, by="term", all.x=TRUE)
  tab[term == REF_GROUP, p_interaction := NA]
  
  tab[, theme := community_themes[term]]
  tab[is.na(theme), theme := term]
  tab[, term_disp := sub("^DOM_", "", term)]
  is_ref <- tab$term == REF_GROUP
  
  data.table(
    Group = sprintf("%s — %s", tab$term_disp, tab$theme),
    `<60 N`              = tab$N_1,
    `<60 Events`         = tab$Events_1,
    `<60 HR (95% CI)`    = fmt_hr(tab$hr_1, tab$lower_1, tab$upper_1, is_ref),
    `60-75 N`            = tab$N_2,
    `60-75 Events`       = tab$Events_2,
    `60-75 HR (95% CI)`  = fmt_hr(tab$hr_2, tab$lower_2, tab$upper_2, is_ref),
    `>=75 N`             = tab$N_3,
    `>=75 Events`        = tab$Events_3,
    `>=75 HR (95% CI)`   = fmt_hr(tab$hr_3, tab$lower_3, tab$upper_3, is_ref),
    P_interaction        = fmt_p(tab$p_interaction)
  )
}

sex_table <- build_sex_table()
age_table <- build_age_table()

# Sort: reference at top, then by community number, MIXED at bottom
sort_key <- function(g) {
  is_ref <- grepl("\\(reference\\)|^C7 —", g)
  is_mix <- grepl("^MIXED ", g)
  num <- suppressWarnings(as.integer(sub(".*?C(\\d+).*", "\\1", g)))
  fifelse(is_ref, -1, fifelse(is_mix, 999, num))
}
sex_table <- sex_table[order(sort_key(sex_table$Group))]
age_table <- age_table[order(sort_key(age_table$Group))]

# =========================
# 7) Excel output
# =========================
cat("\n==> Writing Excel...\n")
wb <- createWorkbook()
hdr_style <- createStyle(fontSize=11, fontColour="#FFFFFF", halign="center",
                         valign="center", fgFill="#4472C4", textDecoration="bold")
note_style <- createStyle(fontSize=10, textDecoration="italic", fontColour="#333333")

add_sheet <- function(name, dt_in, lrt_note) {
  addWorksheet(wb, name)
  writeData(wb, name, lrt_note, startRow=1, startCol=1)
  addStyle(wb, name, note_style, rows=1, cols=1)
  writeDataTable(wb, name, dt_in, startRow=3, withFilter=TRUE)
  setColWidths(wb, name, cols=1:ncol(dt_in), widths="auto")
  addStyle(wb, name, hdr_style, rows=3, cols=1:ncol(dt_in), gridExpand=TRUE)
}

add_sheet("Sex_interaction", sex_table,
          sprintf("Global LRT (sex × community): chi2=%.2f, df=%d, P=%.3g  |  Reference: %s  |  Subgroup Cox adjusted for age",
                  lrt_sex$chi, lrt_sex$df, lrt_sex$p, REF_GROUP))
add_sheet("Age_interaction", age_table,
          sprintf("Global LRT (continuous age × community, per-decade): chi2=%.2f, df=%d, P=%.3g  |  Reference: %s  |  Subgroups by clinical cuts (%s), Cox adjusted for sex",
                  lrt_age$chi, lrt_age$df, lrt_age$p, REF_GROUP, paste(AGE_LABELS, collapse=" / ")))

out_path <- file.path(OUT_DIR, "subgroup_interaction_analysis.xlsx")
saveWorkbook(wb, out_path, overwrite=TRUE)
cat(sprintf("   Saved: %s\n", out_path))
cat("==> Done.\n")