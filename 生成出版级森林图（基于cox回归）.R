#!/usr/bin/env Rscript
# ==============================================================================
# Step 5b v2: Publication-Quality Forest Plots
#   Main analysis (ref = DOM_C7) + 9 sensitivity analyses
#   Each rendered as Lancet/NEJM-style 4-panel table-forest, single page
#   Output: PDF + PNG per analysis in forest_plots/
# ==============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(openxlsx)
  library(patchwork); library(scales)
})

# Null-coalesce helper
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# =========================
# 0) Paths & parameters
# =========================
XLSX_PATH <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results\\cox_regression_results.xlsx"
COX_CSV   <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR   <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results\\forest_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REF_GROUP_MAIN <- "DOM_C7"
MODEL_LBL      <- "Model 2: + age + sex"

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Community theme names
# =========================
community_themes <- c(
  "DOM_C1"  = "Functional & inflammatory GI disorders",
  "DOM_C2"  = "Chronic respiratory disease with infections",
  "DOM_C3"  = "Dermatologic & venous insufficiency disorders",
  "DOM_C4"  = "Benign gynecologic disorders",
  "DOM_C5"  = "Age-related ophthalmologic disorders",
  "DOM_C6"  = "Cerebrovascular disease with neuropsychiatric sequelae",
  "DOM_C7"  = "Upper airway & orodental inflammation",
  "DOM_C8"  = "Hepatobiliary disorders",
  "DOM_C9"  = "Advanced solid malignancies with cachexia",
  "DOM_C10" = "Benign & malignant breast disorders",
  "DOM_C11" = "Urinary stones, BPH & UTI",
  "DOM_C12" = "Thyroid dysfunction & neoplasms",
  "DOM_C13" = "Lymphoid & myeloid malignancies",
  "DOM_C14" = "CKD with anemia, electrolyte & gout",
  "DOM_C15" = "T2DM with neuropathy & joint disease",
  "DOM_C16" = "CAD, heart failure & arrhythmia",
  "MIXED"      = "Mixed (no dominant community)",
  "MIXED"      = "Mixed (no dominant community)",
  "NONE"       = "None (no community diseases at baseline)",
  "DOM_RARE"   = "Rare DOM (merged)"
)

tier_colors <- c(
  "Increased risk"     = "#B2182B",
  "Modestly increased" = "#EF8A62",
  "No difference"      = "#999999",
  "Decreased risk"     = "#2166AC",
  "Reference"          = "#000000"
)

# =========================
# 2) Group sizes from full cox dataset
# =========================
cat("==> Loading cox dataset...\n")
dt <- fread(COX_CSV)
grp_main <- dt[, .(n = .N, events = sum(event == 1L),
                   rate_1000py = round(1000 * sum(event == 1L) / sum(time_years), 2)),
               by = dom_community_unw]
setnames(grp_main, "dom_community_unw", "term")
cat("   Groups:", nrow(grp_main), "\n")

# =========================
# 3) Forest plot factory
# =========================
make_forest <- function(results_dt, ref_group, plot_title, plot_subtitle,
                        output_basename) {
  cat(sprintf("==> %s (ref=%s)\n", output_basename, ref_group))
  
  m2 <- as.data.table(results_dt)
  m2 <- m2[!is.na(HR) & term != ref_group]
  if (nrow(m2) == 0) { cat("   SKIP (empty)\n"); return(invisible(NULL)) }
  m2[, is_ref := FALSE]
  
  ref_row <- data.table(term = ref_group, HR = 1, HR_lower = NA_real_,
                        HR_upper = NA_real_, p_value = NA_real_, is_ref = TRUE)
  fd <- rbindlist(list(m2[, .(term, HR, HR_lower, HR_upper, p_value, is_ref)],
                       ref_row), fill = TRUE)
  fd[, is_ref := as.logical(is_ref)]
  fd[is.na(is_ref), is_ref := FALSE]
  fd <- merge(fd, grp_main, by = "term", all.x = TRUE)
  
  fd[, theme := community_themes[term]]
  fd[is.na(theme), theme := term]
  fd[, term_display := gsub("^DOM_", "", term)]
  fd[, label := sprintf("%s — %s", term_display, theme)]
  fd[is_ref == TRUE, label := paste0(label, " (reference)")]
  
  setorder(fd, -HR)
  fd[, y_pos := .N:1]
  fd[, label := as.character(label)]
  YL <- c(0.5, max(fd$y_pos) + 1.5)
  
  fd[, n_only_label := ifelse(is.na(n), "—", scales::comma(n))]
  fd[, events_rate_label := ifelse(is.na(events), "—",
                                   sprintf("%s (%.1f)", scales::comma(events), rate_1000py))]
  fd[, hr_label := fifelse(is_ref, "1.00 (reference)",
                           sprintf("%.2f (%.2f, %.2f)", HR, HR_lower, HR_upper))]
  fd[, p_label := fifelse(is_ref, "",
                          fifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)))]
  fd[, risk_tier := fcase(
    is_ref, "Reference",
    HR >= 2,   "Increased risk",
    HR >= 1.2, "Modestly increased",
    HR >= 0.8, "No difference",
    default    = "Decreased risk"
  )]
  fd[, risk_tier := factor(risk_tier,
                           levels = c("Increased risk", "Modestly increased", "No difference",
                                      "Decreased risk", "Reference"))]
  
  bt <- theme_classic(base_size = 12) +
    theme(text = element_text(family = "sans", color = "black"),
          plot.title = element_blank(), plot.margin = margin(5, 5, 5, 5),
          axis.title.y = element_blank(), axis.ticks.y = element_blank())
  
  p_label <- ggplot(fd, aes(y = y_pos)) +
    geom_text(aes(x = 0, label = label, fontface = ifelse(is_ref, "bold", "plain")),
              hjust = 0, size = 3.8, family = "sans") +
    annotate("text", x = 0, y = max(fd$y_pos) + 0.8, label = "Dominant community",
             hjust = 0, size = 3.7, fontface = "bold", family = "sans") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") + labs(x = "") + bt +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          plot.margin = margin(20, 5, 5, 5))
  
  p_counts <- ggplot(fd, aes(y = y_pos)) +
    geom_text(aes(x = 0.30, label = n_only_label),
              hjust = 1, size = 3.6, family = "sans") +
    geom_text(aes(x = 0.55, label = events_rate_label),
              hjust = 0, size = 3.6, family = "sans") +
    annotate("text", x = 0.30, y = max(fd$y_pos) + 0.8, label = "N",
             hjust = 1, size = 3.7, fontface = "bold", family = "sans") +
    annotate("text", x = 0.55, y = max(fd$y_pos) + 0.8, label = "Events (rate/1000py)",
             hjust = 0, size = 3.7, fontface = "bold", family = "sans") +
    scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") + bt +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          plot.margin = margin(20, 5, 5, 5))
  
  xbm <- c(0.5, 1, 2, 4, 8); xbn <- c(0.75, 1.5, 3, 6)
  x_min <- min(0.4, min(fd$HR_lower, na.rm = TRUE) * 0.9)
  x_max <- max(fd$HR_upper, na.rm = TRUE) * 1.1
  
  p_forest <- ggplot(fd, aes(y = y_pos, x = HR)) +
    geom_vline(xintercept = xbn, color = "grey95", linewidth = 0.25) +
    geom_vline(xintercept = xbm, color = "grey90", linewidth = 0.35) +
    geom_vline(xintercept = 1, color = "grey40", linewidth = 0.55, linetype = "dashed") +
    geom_segment(data = fd[is_ref == FALSE],
                 aes(x = HR_lower, xend = HR_upper, y = y_pos, yend = y_pos,
                     color = risk_tier),
                 linewidth = 0.6, inherit.aes = FALSE) +
    geom_point(aes(color = risk_tier), shape = 15, size = 3.5) +
    scale_x_log10(breaks = xbm, minor_breaks = xbn, limits = c(x_min, x_max),
                  labels = function(x) ifelse(x == 1, "1.0", as.character(x))) +
    scale_color_manual(values = tier_colors, name = "Risk vs reference", drop = FALSE) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    labs(x = "HR (log scale, 95% CI)") + bt +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 6)),
          axis.text.x = element_text(size = 12), legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11), legend.key.size = unit(0.5, "cm")) +
    guides(color = guide_legend(nrow = 1))
  
  p_stats <- ggplot(fd, aes(y = y_pos)) +
    geom_text(aes(x = 0, label = hr_label,
                  fontface = ifelse(is_ref, "bold", "plain")),
              hjust = 0, size = 3.6, family = "sans") +
    geom_text(aes(x = 1, label = p_label), hjust = 1, size = 3.6, family = "sans") +
    annotate("text", x = 0, y = max(fd$y_pos) + 0.8, label = "HR (95% CI)",
             hjust = 0, size = 3.7, fontface = "bold", family = "sans") +
    annotate("text", x = 1, y = max(fd$y_pos) + 0.8, label = "P-value",
             hjust = 1, size = 3.7, fontface = "bold", family = "sans") +
    scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") + bt +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          plot.margin = margin(20, 5, 5, 5))
  
  combined <- p_label + p_counts + p_forest + p_stats +
    plot_layout(widths = c(4.0, 2.2, 4.0, 1.8)) +
    plot_annotation(title = plot_title, subtitle = plot_subtitle,
                    theme = theme(
                      plot.title = element_text(size = 16, face = "bold", hjust = 0,
                                                margin = margin(b = 4)),
                      plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0,
                                                   margin = margin(b = 10))))
  
  h <- 0.31 * nrow(fd) + 3.0
  pdf_path <- file.path(OUT_DIR, paste0(output_basename, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(output_basename, ".png"))
  ggsave(pdf_path, combined, width = 13, height = h, units = "in", device = PDF_DEVICE)
  ggsave(png_path, combined, width = 13, height = h, units = "in", dpi = 600, bg = "white")
  cat(sprintf("   Saved (%d rows, h=%.1f in)\n", nrow(fd), h))
  invisible(combined)
}

# =========================
# 4) Main analysis
# =========================
cat("\n========== MAIN ==========\n")
main_res <- as.data.table(read.xlsx(XLSX_PATH, sheet = "Main_results"))
m2_main <- main_res[model == MODEL_LBL]
n_total  <- sum(grp_main$n[grp_main$term != "NONE"], na.rm = TRUE)
ev_total <- sum(grp_main$events[grp_main$term != "NONE"], na.rm = TRUE)

make_forest(
  results_dt    = m2_main,
  ref_group     = REF_GROUP_MAIN,
  plot_title    = "All-cause mortality risk by dominant community group",
  plot_subtitle = sprintf(
    "Cox model adjusted for age and sex | Reference: %s (%s) | n = %s, events = %s | 2021-2024 follow-up",
    gsub("^DOM_", "", REF_GROUP_MAIN), community_themes[REF_GROUP_MAIN],
    scales::comma(n_total), scales::comma(ev_total)),
  output_basename = "forest_main"
)

# =========================
# 5) Sensitivity analyses
# =========================
cat("\n========== SENSITIVITY ==========\n")
sens_configs <- list(
  list(sheet="S1_grace",        ref=REF_GROUP_MAIN, suffix="S1_grace",
       title="Sensitivity S1: grace-period censoring",
       desc="Follow-up censored 365d after last visit if no further records"),
  list(sheet="S2_tau035",       ref=REF_GROUP_MAIN, suffix="S2_tau035",
       title="Sensitivity S2: dominance threshold tau = 0.35",
       desc="More patients reclassified as dominant (vs primary tau=0.40)"),
  list(sheet="S3_tau045",       ref=REF_GROUP_MAIN, suffix="S3_tau045",
       title="Sensitivity S3: dominance threshold tau = 0.45",
       desc="Stricter dominance, more patients reclassified as MIXED"),
  list(sheet="S4_weighted",     ref=REF_GROUP_MAIN, suffix="S4_weighted",
       title="Sensitivity S4: edge-weighted community scores",
       desc="Network-edge-weighted scores instead of unweighted"),
  list(sheet="S5_strict_dom",   ref=REF_GROUP_MAIN, suffix="S5_strict_dom",
       title="Sensitivity S5: dominant requires >=2 baseline codes",
       desc="Excludes single-code dominant classifications"),
  list(sheet="S_includeNONE",   ref=REF_GROUP_MAIN, suffix="S_includeNONE",
       title="Sensitivity: NONE group included",
       desc="Patients with no baseline community diseases retained in model"),
  list(sheet="S_mergeNONE",     ref=REF_GROUP_MAIN, suffix="S_mergeNONE",
       title="Sensitivity: NONE merged into MIXED",
       desc="NONE patients reassigned to MIXED rather than excluded")
)

for (cfg in sens_configs) {
  res <- tryCatch(as.data.table(read.xlsx(XLSX_PATH, sheet = cfg$sheet)),
                  error = function(e) { cat("   SKIP", cfg$sheet, ":", e$message, "\n"); NULL })
  if (is.null(res) || nrow(res) == 0) next
  
  if ("model" %in% names(res)) {
    res_m2 <- res[grepl("Model 2", model)]
    if (nrow(res_m2) == 0) res_m2 <- res
  } else {
    res_m2 <- res
  }
  
  ref_use <- cfg$ref
  if ("ref_group" %in% names(res_m2) && nrow(res_m2) > 0 &&
      !is.na(res_m2$ref_group[1]) && res_m2$ref_group[1] != "") {
    ref_use <- res_m2$ref_group[1]
  }
  
  ref_theme <- community_themes[ref_use] %||% ref_use
  
  make_forest(
    results_dt      = res_m2,
    ref_group       = ref_use,
    plot_title      = cfg$title,
    plot_subtitle   = sprintf("%s | Reference: %s (%s)",
                              cfg$desc, gsub("^DOM_", "", ref_use), ref_theme),
    output_basename = paste0("forest_", cfg$suffix)
  )
}

cat("\n==> All forest plots saved to:\n   ", OUT_DIR, "\n==> Done.\n")