#!/usr/bin/env Rscript
# ==============================================================================
# Step 7 v2: Restricted Cubic Splines (RCS) — dose-response for ALL communities
#   Main figure:    2x3, 6 high-mortality communities (HR > 1.5 in main analysis)
#                   C9, C14, C16  (top row)
#                   C13, C2, C6   (bottom row)
#   Supp figure:    2x5, 10 remaining communities (lower-mortality)
#                   C8, C10, C11, C15, C1  (top row)
#                   C3, C4, C5, C7, C12    (bottom row)
#   MIXED is excluded (composite score has no clean dose-response interpretation)
#
# Method: rms::cph with rcs(score, 4 knots), adjusted for age and sex
# CI ribbon: alpha=0.30, Lancet palette
# ==============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages({
  pkgs <- c("data.table", "ggplot2", "openxlsx", "rms", "patchwork", "scales")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, dependencies = TRUE)
  library(data.table); library(ggplot2); library(openxlsx)
  library(rms); library(patchwork); library(scales)
})

INPUT_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR   <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# Lancet-inspired palette (high contrast, journal-safe)
PAL <- c(
  C9  = "#925E9F",  # purple   - solid tumor
  C14 = "#AD002A",  # dark red - CKD
  C16 = "#00468B",  # dark blue - cardiac
  C13 = "#FDAF91",  # peach    - hematologic malignancy
  C2  = "#1B1919",  # near-black - respiratory
  C6  = "#42B540",  # green    - cerebrovascular
  C8  = "#7B5A47",  # brown    - hepatobiliary
  C10 = "#ED0000",  # bright red - breast
  C11 = "#0099B4",  # teal     - urological
  C15 = "#925E9F",  # purple
  C1  = "#00468B",  # blue
  C3  = "#42B540",  # green
  C4  = "#FDAF91",  # peach
  C5  = "#AD002A",  # dark red
  C7  = "#1B1919",  # black
  C12 = "#7B5A47"   # brown
)

# Community config: HR-ordered for visual narrative
make_cfg <- function(label, score, title, subtitle, color) {
  list(label = label, score = score, title = title,
       subtitle = subtitle, color = color)
}

main_communities <- list(
  make_cfg("C9",  "comm9_unw",  "Solid malignancies with cachexia",
           "C9 — gastrointestinal & lung cancers",       PAL["C9"]),
  make_cfg("C14", "comm14_unw", "Chronic kidney disease",
           "C14 — with anemia, electrolyte & gout",      PAL["C14"]),
  make_cfg("C16", "comm16_unw", "Coronary artery disease",
           "C16 — heart failure & arrhythmia",           PAL["C16"]),
  make_cfg("C13", "comm13_unw", "Hematologic malignancies",
           "C13 — lymphoid & myeloid",                    PAL["C13"]),
  make_cfg("C2",  "comm2_unw",  "Chronic respiratory disease",
           "C2 — with systemic infections",               PAL["C2"]),
  make_cfg("C6",  "comm6_unw",  "Cerebrovascular disease",
           "C6 — neuropsychiatric sequelae",              PAL["C6"])
)
supp_communities <- list(
  make_cfg("C8",  "comm8_unw",  "Hepatobiliary disorders", "C8", PAL["C8"]),
  make_cfg("C10", "comm10_unw", "Breast disorders",        "C10", PAL["C10"]),
  make_cfg("C11", "comm11_unw", "Urinary stones, BPH, UTI", "C11", PAL["C11"]),
  make_cfg("C15", "comm15_unw", "T2DM with neuropathy",    "C15", PAL["C15"]),
  make_cfg("C1",  "comm1_unw",  "GI inflammatory disorders","C1", PAL["C1"]),
  make_cfg("C3",  "comm3_unw",  "Skin & venous disorders", "C3", PAL["C3"]),
  make_cfg("C4",  "comm4_unw",  "Gynecologic disorders",   "C4", PAL["C4"]),
  make_cfg("C5",  "comm5_unw",  "Ophthalmologic disorders","C5", PAL["C5"]),
  make_cfg("C7",  "comm7_unw",  "Upper airway & dental",   "C7", PAL["C7"]),
  make_cfg("C12", "comm12_unw", "Thyroid disorders",       "C12", PAL["C12"])
)

KNOTS  <- 4
N_PRED <- 100

# =========================
# Load data
# =========================
cat("==> Loading data...\n")
dt <- fread(INPUT_CSV)
dt <- dt[dom_community_unw != "NONE"]
cat(sprintf("   Cohort: n=%d, events=%d\n", nrow(dt), sum(dt$event)))
dt[, n_diseases_log := log1p(baseline_unique_codes)]
dd <- datadist(dt); options(datadist = "dd")

# =========================
# RCS fitter
# =========================
fit_rcs <- function(cfg) {
  cat(sprintf("==> Fitting %s (%s)...\n", cfg$label, cfg$score))
  if (sum(dt[[cfg$score]] > 0, na.rm = TRUE) < 100) {
    cat("   SKIP: too few non-zero values\n"); return(NULL)
  }
  
  result <- tryCatch({
    # MINIMAL adjustment: only age + sex.
    # n_diseases_log is intentionally OMITTED — it correlates strongly with the
    # community score itself (more codes in community k -> more total codes),
    # so adjusting for it over-adjusts and washes out the score's natural
    # dose-response signal. The Cox dominant-community analysis (Step 5) already
    # adjusts for total disease count; here we want the unadjusted dose-response.
    # Compute knots from NON-ZERO values to avoid collapse on zero-inflated scores
    nz <- dt[[cfg$score]][dt[[cfg$score]] > 0]
    if (length(unique(nz)) < KNOTS) {
      cat("   SKIP: not enough unique non-zero values for", KNOTS, "knots\n")
      return(NULL)
    }
    # Place knots at quantiles of non-zero distribution + anchor at 0
    # This gives the spline meaningful flexibility across the range where signal exists
    manual_knots <- unique(c(0, quantile(nz, probs = seq(0.1, 0.9, length.out = KNOTS - 1),
                                         na.rm = TRUE)))
    if (length(manual_knots) < KNOTS) {
      manual_knots <- unique(c(0, quantile(nz, probs = seq(0.05, 0.95, length.out = KNOTS),
                                           na.rm = TRUE)))
    }
    cat(sprintf("   manual knots (from non-zero): %s\n",
                paste(round(manual_knots, 4), collapse = ",")))
    
    f <- as.formula(sprintf(
      "Surv(time_years, event) ~ rcs(%s, c(%s)) + age_at_index + sex_binary",
      cfg$score, paste(manual_knots, collapse = ",")))
    
    m <- suppressWarnings(cph(f, data = dt, x = TRUE, y = TRUE, surv = TRUE))
    an <- anova(m)
    
    an_rn <- trimws(rownames(an))
    overall_idx <- which(an_rn == cfg$score)
    if (length(overall_idx) == 0)
      overall_idx <- grep(cfg$score, an_rn, fixed = TRUE)[1]
    nonlin_idx <- which(an_rn == "Nonlinear")
    if (length(nonlin_idx) > 1 && length(overall_idx) > 0) {
      nonlin_idx <- nonlin_idx[nonlin_idx > overall_idx][1]
    } else if (length(nonlin_idx) == 0) {
      nonlin_idx <- NA_integer_
    }
    p_col <- grep("^P$|^Pr", colnames(an), value = TRUE)[1]
    overall_p <- if (length(overall_idx) > 0 && !is.na(overall_idx)) an[overall_idx, p_col] else NA
    nonlin_p  <- if (!is.na(nonlin_idx)) an[nonlin_idx, p_col] else NA
    
    knots <- manual_knots
    
    pred <- as.data.table(Predict(m, name = cfg$score, ref.zero = TRUE,
                                  fun = exp, np = N_PRED))
    # Defensive: Predict's first column may not be named cfg$score in all rms versions
    if (!cfg$score %in% names(pred)) setnames(pred, 1, "score")
    else setnames(pred, cfg$score, "score")
    
    cat(sprintf("   knots: %s | overall P=%.3g | nonlinear P=%.3g\n",
                paste(round(knots, 3), collapse = ","), overall_p, nonlin_p))
    list(cfg = cfg, model = m, pred = pred, knots = knots,
         nonlin_p = nonlin_p, overall_p = overall_p)
  }, error = function(e) {
    cat(sprintf("   ERROR fitting %s: %s\n", cfg$label, conditionMessage(e)))
    NULL
  })
  result
}

# =========================
# Panel builder
# =========================
make_panel <- function(res, y_limits = NULL, panel_letter = NULL, base_size = 12) {
  cfg <- res$cfg
  pred <- res$pred
  
  # Build two separate plotmath labels (left-aligned via two annotate() calls)
  # italic(P)[nonlinearity] forces explicit italic on P; subscripts are auto-italic
  fmt_p_math <- function(p) {
    if (is.na(p)) return("'= NA'")
    if (p < 0.001) "'< 0.001'" else sprintf("'= %.3f'", p)
  }
  p_label_nonlin <- sprintf("italic(P)[nonlinearity] ~ %s", fmt_p_math(res$nonlin_p))
  p_label_overall <- sprintf("italic(P)[overall] ~ %s",     fmt_p_math(res$overall_p))
  ymax <- if (!is.null(y_limits)) y_limits[2] else max(pred$upper, na.rm = TRUE)
  ymin <- if (!is.null(y_limits)) y_limits[1] else min(pred$lower, na.rm = TRUE)
  xmin_d <- min(pred$score, na.rm = TRUE)
  xmax_d <- max(pred$score, na.rm = TRUE)
  
  # Panel title: single-line "C{k} — Full theme name", bold
  title_str <- sprintf("%s — %s", cfg$label, cfg$title)
  
  p <- ggplot(pred, aes(x = score, y = yhat)) +
    geom_hline(yintercept = 1, color = "grey50",
               linewidth = 0.45, linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = cfg$color, alpha = 0.30) +
    geom_line(color = cfg$color, linewidth = 1.2) +
    geom_rug(data = data.table(score = res$knots), aes(x = score),
             inherit.aes = FALSE, length = unit(0.025, "npc"),
             color = cfg$color, linewidth = 0.5, sides = "b") +
    # P-value annotation top-left inside plot area
    annotate("text", x = xmin_d + 0.02 * (xmax_d - xmin_d), y = ymax * 0.78,
             label = p_label_nonlin, parse = TRUE,
             hjust = 0, vjust = 1, size = base_size * 0.32,
             family = "sans", color = "grey25") +
    annotate("text", x = xmin_d + 0.02 * (xmax_d - xmin_d), y = ymax * 0.62,
             label = p_label_overall, parse = TRUE,
             hjust = 0, vjust = 1, size = base_size * 0.32,
             family = "sans", color = "grey25") +
    labs(
      title = title_str,
      x = "Community burden score",
      y = "HR (log scale)"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0,
                                margin = margin(b = 6)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin = margin(10, 12, 8, 12)
    )
  
  # Add panel letter (A, B, C, ...) in upper-LEFT corner using tag
  if (!is.null(panel_letter)) {
    p <- p + labs(tag = panel_letter) +
      theme(plot.tag = element_text(size = base_size + 4, face = "bold",
                                    family = "sans"),
            plot.tag.position = c(0.01, 0.99))
  }
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(trans = "log",
                                breaks = c(0.5, 1, 2, 4, 8, 16),
                                labels = c("0.5","1","2","4","8","16"),
                                limits = y_limits)
  } else {
    p <- p + scale_y_continuous(trans = "log",
                                breaks = c(0.5, 1, 2, 4, 8, 16),
                                labels = c("0.5","1","2","4","8","16"))
  }
  p
}

# =========================
# Fit & plot MAIN figure (2x3, 6 high-HR communities)
# =========================
cat("\n========== MAIN FIGURE (6 high-mortality communities) ==========\n")
res_main <- lapply(main_communities, fit_rcs)
res_main <- res_main[!sapply(res_main, is.null)]

panels_main <- mapply(make_panel, res_main,
                      panel_letter = LETTERS[seq_along(res_main)],
                      MoreArgs = list(y_limits = c(0.5, 20), base_size = 11),
                      SIMPLIFY = FALSE)

main_combined <- wrap_plots(panels_main, nrow = 3, ncol = 2) +
  plot_annotation(
    title = "Dose-response curves for high-mortality community burden and all-cause mortality",
    subtitle = sprintf(
      "Cox model with restricted cubic splines (4 knots), adjusted for age and sex | HR, hazard ratio (vs zero burden) with 95%% CI | n = %s, events = %s",
      scales::comma(nrow(dt)), scales::comma(sum(dt$event))),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0,
                                   margin = margin(b = 10))
    )
  )

ggsave(file.path(OUT_DIR, "rcs_main_6communities.pdf"), main_combined,
       width = 8.27, height = 11.69, units = "in", device = PDF_DEVICE)
ggsave(file.path(OUT_DIR, "rcs_main_6communities.png"), main_combined,
       width = 8.27, height = 11.69, units = "in", dpi = 600, bg = "white")
cat("   Saved: rcs_main_6communities.pdf / .png\n")

# =========================
# Fit & plot SUPPLEMENTARY figure (2x5, 10 lower-HR communities)
# =========================
cat("\n========== SUPPLEMENTARY FIGURE (10 lower-mortality communities) ==========\n")
res_supp <- lapply(supp_communities, fit_rcs)
res_supp <- res_supp[!sapply(res_supp, is.null)]

panels_supp <- mapply(make_panel, res_supp,
                      panel_letter = LETTERS[6 + seq_along(res_supp)],
                      MoreArgs = list(y_limits = c(0.5, 2.5), base_size = 11),
                      SIMPLIFY = FALSE)

supp_combined <- wrap_plots(panels_supp, nrow = 5, ncol = 2) +
  plot_annotation(
    title = "Dose-response curves for lower-mortality communities (Supplementary)",
    subtitle = sprintf(
      "Panels G\u2013P continue from the main figure (panels A\u2013F) | Cox model with restricted cubic splines (4 knots), adjusted for age and sex | HR, hazard ratio (vs zero burden) with 95%% CI | y-axis fixed at HR 0.5\u20132.5 for cross-community comparison | n = %s, events = %s",
      scales::comma(nrow(dt)), scales::comma(sum(dt$event))),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0,
                                   margin = margin(b = 10))
    )
  )

ggsave(file.path(OUT_DIR, "rcs_supplementary_10communities.pdf"), supp_combined,
       width = 8.27, height = 11.69, units = "in", device = PDF_DEVICE)
ggsave(file.path(OUT_DIR, "rcs_supplementary_10communities.png"), supp_combined,
       width = 8.27, height = 11.69, units = "in", dpi = 600, bg = "white")
cat("   Saved: rcs_supplementary_10communities.pdf / .png\n")

# =========================
# Excel summary: all 16 RCS results
# =========================
cat("\n==> Writing Excel summary...\n")
all_res <- c(res_main, res_supp)

summary_dt <- rbindlist(lapply(all_res, function(r) data.table(
  community = r$cfg$label,
  score_column = r$cfg$score,
  group = ifelse(r$cfg$label %in% sapply(main_communities, `[[`, "label"),
                 "main", "supplementary"),
  n_knots = KNOTS,
  knot_positions = paste(round(r$knots, 4), collapse = ", "),
  overall_p = r$overall_p,
  nonlinear_p = r$nonlin_p,
  overall_p_fmt    = ifelse(r$overall_p < 0.001, "<0.001", sprintf("%.3g", r$overall_p)),
  nonlinear_p_fmt  = ifelse(r$nonlin_p  < 0.001, "<0.001", sprintf("%.3g", r$nonlin_p))
)))

wb <- createWorkbook()
hdr <- createStyle(fontSize=11, fontColour="#FFFFFF", halign="center",
                   valign="center", fgFill="#4472C4", textDecoration="bold")
addWorksheet(wb, "RCS_summary_all16")
writeDataTable(wb, "RCS_summary_all16", summary_dt, withFilter = TRUE)
setColWidths(wb, "RCS_summary_all16", cols = 1:ncol(summary_dt), widths = "auto")
addStyle(wb, "RCS_summary_all16", hdr, rows = 1, cols = 1:ncol(summary_dt), gridExpand = TRUE)

for (r in all_res) {
  sn <- paste0("Curve_", r$cfg$label)
  out <- r$pred[, .(score, HR = yhat, HR_lower = lower, HR_upper = upper)]
  addWorksheet(wb, sn)
  writeDataTable(wb, sn, out, withFilter = TRUE)
  setColWidths(wb, sn, cols = 1:ncol(out), widths = "auto")
  addStyle(wb, sn, hdr, rows = 1, cols = 1:ncol(out), gridExpand = TRUE)
}

saveWorkbook(wb, file.path(OUT_DIR, "rcs_all_communities_results.xlsx"), overwrite = TRUE)
cat("   Saved: rcs_all_communities_results.xlsx\n")
cat("\n==> Done.\n")