#!/usr/bin/env Rscript
# ==============================================================================
# Step 8: Population Attributable Risk (PAR) analysis
#
# Method: Miettinen's formula for adjusted HRs
#   PAR_k = pc_k * (HR_k - 1) / HR_k
#   where pc_k = proportion of cases (deaths) in community k,
#         HR_k = age+sex adjusted HR from Step 5 Model 2
#
# Reference group: DOM_C7 (matches main analysis)
# 95% CIs: computed by substituting HR_lower / HR_upper into the Miettinen formula
#          (Rothman-Greenland approximation, standard in the literature).
#
# Cumulative PAR uses the independence formula:
#   cumPAR = 1 - product(1 - PAR_k)
#
# Subgroups: main (overall), by sex (M/F), by age (<60, 60-75, >=75)
#
# Inputs:
#   cox_regression_results.xlsx  (Step 5 Main_results)
#   subgroup_interaction_analysis.xlsx  (Step 6 Sex_interaction, Age_interaction)
#   cox_dataset_with_dominant_community.csv  (for prevalence/event counts)
#
# Outputs (all in cox_results/):
#   par_analysis.xlsx
#   par_main_barchart.pdf / .png
#   par_sex_barchart.pdf / .png  (side-by-side male/female)
# ==============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(openxlsx); library(patchwork)
})

# =========================
# 0) Paths & parameters
# =========================
XLSX_COX <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results\\cox_regression_results.xlsx"
XLSX_SUB <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results\\subgroup_interaction_analysis.xlsx"
COX_CSV  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REF_GROUP <- "DOM_C7"
MODEL_LBL <- "Model 2: + age + sex"

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# Community theme names
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

tier_colors <- c(
  "Increased risk"     = "#B2182B",
  "Modestly increased" = "#EF8A62",
  "No difference"      = "#999999",
  "Decreased risk"     = "#2166AC",
  "Reference"          = "#000000"
)

# =========================
# 1) Miettinen PAR core
# =========================
miettinen_par <- function(hr, hr_lo, hr_up, pc) {
  # Returns list(par, par_lo, par_up) — matched by POSITION to LHS of := assignment
  par_pt <- pc * (hr - 1) / hr
  par_lo_raw <- pc * (hr_lo - 1) / hr_lo
  par_up_raw <- pc * (hr_up - 1) / hr_up
  # Reorder lo/up if HR<1 flipped them
  out_lo <- pmin(par_lo_raw, par_up_raw)
  out_up <- pmax(par_lo_raw, par_up_raw)
  list(par_pt, out_lo, out_up)  # unnamed: matched positionally
}

# =========================
# 2) Load data
# =========================
cat("==> Loading Step 5 Main_results...\n")
main_res <- as.data.table(read.xlsx(XLSX_COX, sheet = "Main_results"))
m2 <- main_res[model == MODEL_LBL & !is.na(HR) & term != REF_GROUP]
cat(sprintf("   %d non-reference groups in Model 2\n", nrow(m2)))

cat("==> Loading cox dataset for event counts...\n")
dt <- fread(COX_CSV)
dt <- dt[dom_community_unw != "NONE"]

# Total events in overall, male, female, and age groups
total_events <- sum(dt$event == 1L)
total_events_m <- sum(dt$event == 1L & dt$sex_binary == 1L)
total_events_f <- sum(dt$event == 1L & dt$sex_binary == 0L)

dt[, age_group := cut(age_at_index, breaks = c(-Inf, 60, 75, Inf),
                      labels = c("<60", "60-75", ">=75"), right = FALSE)]
total_events_a1 <- sum(dt$event == 1L & dt$age_group == "<60")
total_events_a2 <- sum(dt$event == 1L & dt$age_group == "60-75")
total_events_a3 <- sum(dt$event == 1L & dt$age_group == ">=75")

# Event counts per community (overall)
events_by_grp <- dt[, .(events = sum(event == 1L), N = .N), by = dom_community_unw]
setnames(events_by_grp, "dom_community_unw", "term")

# =========================
# 3) Compute overall PAR
# =========================
cat("\n==> Computing overall PAR (Miettinen formula)...\n")
par_dt <- merge(m2[, .(term, HR, HR_lower, HR_upper)],
                events_by_grp, by = "term", all.x = TRUE)
par_dt[, pc := events / total_events]
par_dt[, c("par", "par_lo", "par_up") := miettinen_par(HR, HR_lower, HR_upper, pc)]
par_dt[, theme := community_themes[term]]
par_dt[, term_display := gsub("^DOM_", "", term)]
par_dt[, label := sprintf("%s — %s", term_display, theme)]

setorder(par_dt, -par)
par_dt[, rank := .I]

# Cumulative PAR (independence assumption)
par_dt[, par_pos := fifelse(is.na(par) | par < 0, 0, par)]
par_dt[, cum_par := 1 - cumprod(1 - par_pos)]

# Risk tier (for coloring, consistent with forest plot)
par_dt[, risk_tier := fcase(
  HR >= 2,   "Increased risk",
  HR >= 1.2, "Modestly increased",
  HR >= 0.8, "No difference",
  default    = "Decreased risk"
)]
par_dt[, risk_tier := factor(risk_tier,
                             levels = c("Increased risk", "Modestly increased",
                                        "No difference", "Decreased risk"))]

cat("   Top 5 PAR contributors:\n")
print(par_dt[1:5, .(term, HR, pc = round(pc, 3), par = round(par, 4),
                    cum_par = round(cum_par, 4))])
cat(sprintf("\n   Cumulative PAR (top 5): %.1f%%\n", 100 * par_dt$cum_par[5]))
cat(sprintf("   Cumulative PAR (top 10): %.1f%%\n", 100 * par_dt$cum_par[10]))
cat(sprintf("   Total PAR (all %d groups): %.1f%%\n",
            nrow(par_dt), 100 * par_dt$cum_par[nrow(par_dt)]))

# =========================
# 4) Subgroup PAR: sex
# =========================
cat("\n==> Computing sex-stratified PAR...\n")
# Parse Step 6 Sex_interaction sheet (row 1 is note, data starts row 3)
sex_raw <- as.data.table(read.xlsx(XLSX_SUB, sheet = "Sex_interaction", startRow = 3))
# Extract numeric HR from "X.XX (X.XX-X.XX)" formatted strings
parse_hr <- function(s) {
  if (is.null(s) || length(s) == 0) {
    return(matrix(numeric(0), ncol = 3,
                  dimnames = list(NULL, c("hr","lo","up"))))
  }
  s <- as.character(s)
  m <- regmatches(s, regexec("([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\)", s))
  vals <- lapply(seq_along(s), function(i) {
    x <- m[[i]]
    if (length(x) < 4) return(c(NA_real_, NA_real_, NA_real_))
    as.numeric(x[2:4])
  })
  mat <- matrix(unlist(vals), ncol = 3, byrow = TRUE)
  colnames(mat) <- c("hr", "lo", "up")
  mat
}

# Fuzzy column lookup: find a column whose name matches a keyword pattern
# Handles openxlsx auto-renaming (spaces/special chars -> dots)
get_col <- function(dt_obj, pattern) {
  cn <- names(dt_obj)
  hit <- grep(pattern, cn, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(hit)) {
    cat(sprintf("   !! Column not found matching '%s'. Available: %s\n",
                pattern, paste(cn, collapse = " | ")))
    return(rep(NA_character_, nrow(dt_obj)))
  }
  dt_obj[[hit]]
}
sex_raw[, term := sub("^([A-Za-z0-9_]+).*$", "\\1", Group)]
sex_raw[, term := fifelse(grepl("^C[0-9]", term), paste0("DOM_", term), term)]

m_hr <- parse_hr(get_col(sex_raw, "Male.*HR"))
f_hr <- parse_hr(get_col(sex_raw, "Female.*HR"))
sex_par <- data.table(
  term = sex_raw$term,
  N_M = as.numeric(get_col(sex_raw, "^Male\\.?N|Male N")),
  events_M = as.numeric(get_col(sex_raw, "Male.*Events")),
  HR_M = m_hr[, 1], HR_M_lo = m_hr[, 2], HR_M_up = m_hr[, 3],
  N_F = as.numeric(get_col(sex_raw, "^Female\\.?N|Female N")),
  events_F = as.numeric(get_col(sex_raw, "Female.*Events")),
  HR_F = f_hr[, 1], HR_F_lo = f_hr[, 2], HR_F_up = f_hr[, 3]
)
# Drop reference row (HR=1, no PAR)
sex_par <- sex_par[term != REF_GROUP]
sex_par[, pc_M := events_M / total_events_m]
sex_par[, pc_F := events_F / total_events_f]

# PAR for male and female separately
sex_par[, c("par_M", "par_M_lo", "par_M_up") :=
          miettinen_par(HR_M, HR_M_lo, HR_M_up, pc_M)]
sex_par[, c("par_F", "par_F_lo", "par_F_up") :=
          miettinen_par(HR_F, HR_F_lo, HR_F_up, pc_F)]
sex_par[, theme := community_themes[term]]
sex_par[, term_display := sub("^DOM_", "", term)]
sex_par[, label := sprintf("%s — %s", term_display, theme)]

# =========================
# 5) Subgroup PAR: age
# =========================
cat("==> Computing age-stratified PAR...\n")
age_raw <- as.data.table(read.xlsx(XLSX_SUB, sheet = "Age_interaction", startRow = 3))
age_raw[, term := sub("^([A-Za-z0-9_]+).*$", "\\1", Group)]
age_raw[, term := fifelse(grepl("^C[0-9]", term), paste0("DOM_", term), term)]

hr_a1 <- parse_hr(get_col(age_raw, "X\\.60.*HR|<60.*HR"))
hr_a2 <- parse_hr(get_col(age_raw, "60.*75.*HR"))
hr_a3 <- parse_hr(get_col(age_raw, "75.*HR|X\\.75.*HR"))

age_par <- data.table(
  term = age_raw$term,
  events_1 = as.numeric(get_col(age_raw, "X\\.60.*Events|<60.*Events")),
  HR_1 = hr_a1[, 1], HR_1_lo = hr_a1[, 2], HR_1_up = hr_a1[, 3],
  events_2 = as.numeric(get_col(age_raw, "60.*75.*Events")),
  HR_2 = hr_a2[, 1], HR_2_lo = hr_a2[, 2], HR_2_up = hr_a2[, 3],
  events_3 = as.numeric(get_col(age_raw, "X\\.75.*Events|>=75.*Events")),
  HR_3 = hr_a3[, 1], HR_3_lo = hr_a3[, 2], HR_3_up = hr_a3[, 3]
)
age_par <- age_par[term != REF_GROUP]
age_par[, pc_1 := events_1 / total_events_a1]
age_par[, pc_2 := events_2 / total_events_a2]
age_par[, pc_3 := events_3 / total_events_a3]
age_par[, c("par_1", "par_1_lo", "par_1_up") := miettinen_par(HR_1, HR_1_lo, HR_1_up, pc_1)]
age_par[, c("par_2", "par_2_lo", "par_2_up") := miettinen_par(HR_2, HR_2_lo, HR_2_up, pc_2)]
age_par[, c("par_3", "par_3_lo", "par_3_up") := miettinen_par(HR_3, HR_3_lo, HR_3_up, pc_3)]
age_par[, theme := community_themes[term]]
age_par[, term_display := sub("^DOM_", "", term)]
age_par[, label := sprintf("%s — %s", term_display, theme)]

# =========================
# 6) Main PAR bar chart
# =========================
cat("\n==> Building main PAR bar chart...\n")

# Drop rows with negative PAR (protective communities — PAR interpretation doesn't apply)
par_plot <- par_dt[par >= 0]
par_plot[, label := factor(label, levels = rev(label))]  # rev so highest is top
par_plot[, par_pct := par * 100]
par_plot[, par_lo_pct := par_lo * 100]
par_plot[, par_up_pct := par_up * 100]
par_plot[, label_text := sprintf("%.1f%%", par_pct)]

p_main <- ggplot(par_plot, aes(x = par_pct, y = label, fill = risk_tier)) +
  geom_col(width = 0.68) +
  geom_errorbarh(aes(xmin = pmax(0, par_lo_pct), xmax = par_up_pct),
                 height = 0.25, color = "grey30", linewidth = 0.4) +
  geom_text(aes(x = par_up_pct, label = label_text),
            hjust = -0.15, size = 3.6, family = "sans",
            color = "grey15") +
  scale_fill_manual(values = tier_colors, name = "Risk tier (HR)", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "Population attributable risk (PAR) for all-cause mortality",
    subtitle = sprintf(
      "Miettinen formula, adjusted for age and sex | Reference: C7 | n = %s, events = %s | Top-5 cumulative PAR = %.1f%%",
      format(nrow(dt), big.mark = ","),
      format(total_events, big.mark = ","),
      100 * par_dt$cum_par[5]),
    x = "Population attributable risk (%)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0,
                                 margin = margin(b = 10)),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

ggsave(file.path(OUT_DIR, "par_main_barchart.pdf"), p_main,
       width = 11, height = 7, units = "in", device = PDF_DEVICE)
ggsave(file.path(OUT_DIR, "par_main_barchart.png"), p_main,
       width = 11, height = 7, units = "in", dpi = 600, bg = "white")
cat("   Saved: par_main_barchart.pdf / .png\n")

# =========================
# 7) Sex-stratified PAR bar chart (side-by-side)
# =========================
cat("==> Building sex-stratified PAR bar chart...\n")

sex_long <- rbind(
  sex_par[, .(term, label, sex = "Male",   par_pct = 100 * par_M,
              par_lo_pct = 100 * par_M_lo, par_up_pct = 100 * par_M_up,
              HR = HR_M)],
  sex_par[, .(term, label, sex = "Female", par_pct = 100 * par_F,
              par_lo_pct = 100 * par_F_lo, par_up_pct = 100 * par_F_up,
              HR = HR_F)]
)
sex_long <- sex_long[!is.na(par_pct) & par_pct > -1]  # drop NA/garbage rows
# Order by average PAR across sexes
ord <- sex_long[, .(avg = mean(par_pct, na.rm = TRUE)), by = label][order(-avg), label]
sex_long[, label := factor(label, levels = rev(ord))]
sex_long[, sex := factor(sex, levels = c("Male", "Female"))]

p_sex <- ggplot(sex_long, aes(x = par_pct, y = label, fill = sex)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbarh(aes(xmin = pmax(0, par_lo_pct), xmax = par_up_pct),
                 height = 0.25, color = "grey30", linewidth = 0.35,
                 position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("Male" = "#00468B", "Female" = "#AD002A"),
                    name = "Sex") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "Sex-stratified PAR for all-cause mortality",
    subtitle = "Miettinen formula, adjusted for age | Reference: C7",
    x = "Population attributable risk (%)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0,
                                 margin = margin(b = 10)),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 6)),
    axis.text = element_text(size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

ggsave(file.path(OUT_DIR, "par_sex_barchart.pdf"), p_sex,
       width = 11, height = 7, units = "in", device = PDF_DEVICE)
ggsave(file.path(OUT_DIR, "par_sex_barchart.png"), p_sex,
       width = 11, height = 7, units = "in", dpi = 600, bg = "white")
cat("   Saved: par_sex_barchart.pdf / .png\n")

# =========================
# 8) Excel output
# =========================
cat("\n==> Writing Excel summary...\n")
fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)
fmt_ci <- function(pt, lo, up) sprintf("%.2f%% (%.2f-%.2f)",
                                       100 * pt, 100 * lo, 100 * up)

# Main PAR table
out_main <- par_dt[, .(
  Rank = rank,
  Community = label,
  N = N,
  Events = events,
  `Exposure prevalence` = sprintf("%.2f%%", 100 * N / nrow(dt)),
  `Cases exposed (pc)`  = sprintf("%.2f%%", 100 * pc),
  `Adjusted HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper),
  `PAR (95% CI)` = fmt_ci(par, par_lo, par_up),
  `PAR %` = sprintf("%.2f", 100 * par),
  `Cumulative PAR %` = sprintf("%.2f", 100 * cum_par)
)]

# Sex PAR table
out_sex <- sex_par[, .(
  Community = label,
  `Male HR` = sprintf("%.2f", HR_M),
  `Male PAR (95% CI)` = fmt_ci(par_M, par_M_lo, par_M_up),
  `Female HR` = sprintf("%.2f", HR_F),
  `Female PAR (95% CI)` = fmt_ci(par_F, par_F_lo, par_F_up)
)]

# Age PAR table
out_age <- age_par[, .(
  Community = label,
  `<60 HR` = sprintf("%.2f", HR_1),
  `<60 PAR (95% CI)` = fmt_ci(par_1, par_1_lo, par_1_up),
  `60-75 HR` = sprintf("%.2f", HR_2),
  `60-75 PAR (95% CI)` = fmt_ci(par_2, par_2_lo, par_2_up),
  `>=75 HR` = sprintf("%.2f", HR_3),
  `>=75 PAR (95% CI)` = fmt_ci(par_3, par_3_lo, par_3_up)
)]

wb <- createWorkbook()
hdr <- createStyle(fontSize=11, fontColour="#FFFFFF", halign="center",
                   valign="center", fgFill="#4472C4", textDecoration="bold")
note_style <- createStyle(fontSize=10, textDecoration="italic", fontColour="#333333")

add_sheet <- function(name, dt_in, note) {
  addWorksheet(wb, name)
  writeData(wb, name, note, startRow = 1, startCol = 1)
  addStyle(wb, name, note_style, rows = 1, cols = 1)
  writeDataTable(wb, name, dt_in, startRow = 3, withFilter = TRUE)
  setColWidths(wb, name, cols = 1:ncol(dt_in), widths = "auto")
  addStyle(wb, name, hdr, rows = 3, cols = 1:ncol(dt_in), gridExpand = TRUE)
}

add_sheet("PAR_main", out_main,
          sprintf("Miettinen formula PAR = pc * (HR-1)/HR | Adjusted HR from Step 5 Model 2 (age + sex) | Reference: %s | Cumulative PAR assumes independence: 1 - prod(1 - PAR_k)",
                  REF_GROUP))
add_sheet("PAR_sex", out_sex,
          "Sex-stratified PAR | Subgroup Cox adjusted for age | Reference: C7 | pc computed within each sex")
add_sheet("PAR_age", out_age,
          "Age-stratified PAR | Subgroup Cox adjusted for sex | Reference: C7 | pc computed within each age group")

saveWorkbook(wb, file.path(OUT_DIR, "par_analysis.xlsx"), overwrite = TRUE)
cat("   Saved: par_analysis.xlsx\n")
cat("\n==> Done.\n")