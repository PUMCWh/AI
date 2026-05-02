#!/usr/bin/env Rscript
# ==============================================================================
# Step 9: ICD-10 Chapter Enrichment Analysis for CKM Disease Communities
#
# Purpose:
#   Validate that the 16 data-driven Leiden communities are not random clusters
#   by testing whether each community is significantly enriched in (or depleted
#   of) diseases from specific ICD-10 chapters. Confirms biological/clinical
#   interpretability of the community structure.
#
# Inputs (from CKM pipeline):
#   - medoid_bridge_metrics_coef0_01_alpha0_05.csv
#     (provides: node → community → ICD10_Chapter mapping for 613 disease nodes)
#
# Outputs (written to OUT_DIR):
#   - FigS2_ICD10_chapter_enrichment_heatmap.pdf
#   - TableS6_ICD10_enrichment_full_results.xlsx
#       Sheet 1: Significant_Enrichment (FDR<0.05, for manuscript)
#       Sheet 2: Full_Results       (all community × chapter tests)
#       Sheet 3: Community_Summary  (one-line summary per community)
#
# Method:
#   For each (community, ICD-10 chapter) cell, build a 2×2 contingency table:
#                    in chapter    not in chapter
#     in community       a              b         | n_community
#     not in comm        c              d         | N - n_community
#                     n_chapter     N - n_chapter  | N
#
#   Run two-sided Fisher's exact test; apply BH FDR correction over all
#   community × chapter combinations; compute observed/expected ratio and
#   log2(O/E) as effect size for heatmap colouring.
#
# Design for CKM paper:
#   - 16 communities (γ = 1.2; n = 640 637 CKM-relevant adults)
#   - Diseases labelled with ICD-10 3-character codes only (no disease names)
#   - ICD-10 chapter codes 1–14 (following CKM paper convention; chapter 15-22
#     typically absent from CKM cohort)
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

# =========================
# 0) Paths & Parameters
# =========================
BRIDGE_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\community_detection\\medoid_bridge_metrics_coef0.01_alpha0.05.csv"
OUT_DIR    <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FDR_THRESH  <- 0.05
FONT_FAMILY <- "Arial"
PDF_W       <- 11   # wider for 16 communities on x-axis
PDF_H       <- 7

# ICD-10 chapter labels (numeric → Roman, following CKM docx convention)
chapter_labels <- c(
  "1"  = "I",
  "2"  = "II",
  "3"  = "III",
  "4"  = "IV",
  "5"  = "V",
  "6"  = "VI",
  "7"  = "VII",
  "8"  = "VIII",
  "9"  = "IX",
  "10" = "X",
  "11" = "XI",
  "12" = "XII",
  "13" = "XIII",
  "14" = "XIV"
)

# Cairo PDF device with fallback
PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Load & prepare data
# =========================
cat("==> Step 1) Loading bridge metrics with community membership and chapter info...\n")

if (!file.exists(BRIDGE_CSV)) {
  stop("Bridge metrics CSV not found: ", BRIDGE_CSV, "\n",
       "This file is produced by Step 3 (Leiden community detection).", call. = FALSE)
}

bridge_dt <- fread(BRIDGE_CSV)
cat("   Loaded:", nrow(bridge_dt), "nodes\n")
cat("   Columns:", paste(names(bridge_dt), collapse = ", "), "\n")

# Keep only nodes assigned to a community (community > 0) and with chapter info
req_cols <- c("node", "community", "ICD10_Chapter")
missing <- setdiff(req_cols, names(bridge_dt))
if (length(missing) > 0) {
  stop("Required columns missing from bridge CSV: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

bridge_dt[, ICD10_Chapter := as.character(ICD10_Chapter)]
nodes_dt <- bridge_dt[community > 0 &
                        !is.na(ICD10_Chapter) &
                        ICD10_Chapter != "" &
                        ICD10_Chapter != "0",
                      .(node, community, ICD10_Chapter)]

cat("   Active nodes (community > 0, with chapter):", nrow(nodes_dt), "\n")
cat("   Communities:", length(unique(nodes_dt$community)), "\n")
cat("   Chapters represented:",
    paste(sort(unique(nodes_dt$ICD10_Chapter)), collapse = ", "), "\n")

# =========================
# 2) Fisher's exact test per (community × chapter) pair
# =========================
cat("\n==> Step 2) Running Fisher's exact tests (community × chapter)...\n")

all_comms    <- sort(unique(nodes_dt$community))
all_chapters <- sort(unique(nodes_dt$ICD10_Chapter))
N_total      <- nrow(nodes_dt)

cat("   Tests to run:", length(all_comms) * length(all_chapters), "\n")
cat("   Background N:", N_total, "nodes\n")

results <- list()
idx <- 0
for (cid in all_comms) {
  n_comm <- nodes_dt[community == cid, .N]
  
  for (ch in all_chapters) {
    # Contingency table cells
    a    <- nodes_dt[community == cid & ICD10_Chapter == ch, .N]
    n_ch <- nodes_dt[ICD10_Chapter == ch, .N]
    b    <- n_comm - a
    c_v  <- n_ch - a
    d    <- (N_total - n_comm) - c_v
    
    mat_2x2 <- matrix(c(a, b, c_v, d), nrow = 2, byrow = TRUE)
    ft <- suppressWarnings(fisher.test(mat_2x2))
    
    or_raw <- as.numeric(ft$estimate)
    if (length(or_raw) == 0) or_raw <- NA_real_
    
    ci_lo <- as.numeric(ft$conf.int[1])
    ci_hi <- as.numeric(ft$conf.int[2])
    if (length(ci_lo) == 0) ci_lo <- NA_real_
    if (length(ci_hi) == 0) ci_hi <- NA_real_
    
    expected      <- n_comm * n_ch / N_total
    obs_exp_ratio <- if (expected > 0) a / expected else NA_real_
    
    idx <- idx + 1
    results[[idx]] <- data.table(
      community        = cid,
      chapter          = ch,
      n_in_comm_ch     = a,
      n_in_comm        = n_comm,
      n_in_chapter     = n_ch,
      n_total          = N_total,
      pct_in_comm      = round(a / n_comm * 100, 1),
      pct_in_network   = round(n_ch / N_total * 100, 1),
      or_raw           = or_raw,
      odds_ratio       = round(or_raw, 3),
      obs_exp_ratio    = round(obs_exp_ratio, 3),
      p_value          = ft$p.value,
      ci_low           = round(ci_lo, 3),
      ci_high          = round(ci_hi, 3)
    )
  }
}

enrich_dt <- rbindlist(results)

# BH FDR correction across ALL community × chapter tests
enrich_dt[, fdr := p.adjust(p_value, method = "BH")]

# Direction (use un-rounded OR to avoid boundary misclassification)
enrich_dt[, direction := fifelse(is.na(or_raw), "unknown",
                                 fifelse(or_raw > 1,   "enriched",
                                         fifelse(or_raw < 1,   "depleted", "neutral")))]
enrich_dt[, significant := fdr < FDR_THRESH]

cat(sprintf("   Total tests: %d\n", nrow(enrich_dt)))
cat(sprintf("   Significant (FDR < %.2f): %d (enriched: %d, depleted: %d)\n",
            FDR_THRESH,
            sum(enrich_dt$significant),
            sum(enrich_dt$significant & enrich_dt$direction == "enriched"),
            sum(enrich_dt$significant & enrich_dt$direction == "depleted")))

# =========================
# 3) Prepare heatmap data
# =========================
cat("\n==> Step 3) Preparing heatmap (log2 O/E fill, significance star overlay)...\n")

# Fill value: log2(O/E); non-significant cells forced to 0 (white)
# O/E = 0 (observed count = 0) at significant FDR → use -max_abs as "fully depleted"
enrich_dt[, log2_oe := fifelse(obs_exp_ratio > 0 & is.finite(obs_exp_ratio),
                               log2(obs_exp_ratio), NA_real_)]
enrich_dt[, log2_oe_plot := fifelse(!significant, 0,
                                    fifelse(!is.na(log2_oe), log2_oe, NA_real_))]

plot_dt <- copy(enrich_dt)

# Community label (C1, C2, ..., C16), ordered by number
plot_dt[, comm_label := factor(paste0("C", community),
                               levels = paste0("C", sort(all_comms)))]

# Chapter labels (Roman + description), ordered I (top) to XIV (bottom)
plot_dt[, chapter_label := chapter_labels[chapter]]
plot_dt[is.na(chapter_label), chapter_label := paste0("Ch.", chapter)]
ch_order <- chapter_labels[as.character(sort(as.integer(all_chapters)))]
plot_dt[, chapter_label := factor(chapter_label, levels = rev(ch_order))]
# rev() so that Roman I appears at TOP of y-axis (natural reading order)

# Significance stars
plot_dt[, sig_label := fifelse(significant & fdr < 0.001, "***",
                               fifelse(significant & fdr < 0.01,  "**",
                                       fifelse(significant & fdr < 0.05,  "*",  "")))]

# In-cell observed count
plot_dt[, cell_text := as.character(n_in_comm_ch)]

# Symmetric colour scale around 0
max_abs <- max(abs(plot_dt$log2_oe_plot), na.rm = TRUE)
max_abs <- max(2, ceiling(max_abs))  # floor at 2 for visual clarity

# Replace NA (significant depletion where count=0) with -max_abs
plot_dt[is.na(log2_oe_plot), log2_oe_plot := -max_abs]

# =========================
# 4) Build heatmap
# =========================
p_heatmap <- ggplot(plot_dt, aes(x = comm_label, y = chapter_label)) +
  geom_tile(aes(fill = log2_oe_plot), color = "grey85", linewidth = 0.3) +
  geom_text(aes(label = cell_text), size = 2.9, color = "grey20",
            family = FONT_FAMILY, vjust = -0.3) +
  geom_text(aes(label = sig_label), size = 3.3, color = "black",
            family = FONT_FAMILY, fontface = "bold", vjust = 1.2) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    name     = expression(log[2](O/E)),
    breaks   = seq(-max_abs, max_abs, by = max(1, floor(max_abs / 3)))
  ) +
  labs(
    title    = "ICD-10 Chapter Enrichment across 16 CKM Disease Communities",
    subtitle = paste0("Fisher's exact test, Benjamini-Hochberg FDR corrected; ",
                      "* FDR<0.05, ** FDR<0.01, *** FDR<0.001; ",
                      "white cells: FDR \u2265 0.05"),
    x = "Disease community",
    y = "ICD-10 chapter"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text              = element_text(family = FONT_FAMILY),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 12, color = "black"),
    plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 9),
    axis.text.x       = element_text(angle = 0,  hjust = 0.5, size = 12, color = "black"),
    axis.text.y       = element_text(size = 11, color = "black"),
    axis.title        = element_text(size = 13, color = "black"),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm"),
    panel.grid        = element_blank()
  )

pdf_file <- file.path(OUT_DIR, "ICD10_chapter_enrichment_heatmap.pdf")
ggsave(pdf_file, p_heatmap, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved:", basename(pdf_file), "\n")

# =========================
# 5) Supplementary Table S6 (multi-sheet Excel)
# =========================
cat("\n==> Step 5) Writing Supplementary Table S6 (Excel with 3 sheets)...\n")

# Sheet 1: significant only (for manuscript)
sig_dt <- enrich_dt[significant == TRUE]
setorder(sig_dt, community, -obs_exp_ratio)

fmt_or_ci <- function(or_val, lo, hi) {
  or_s <- fifelse(is.na(or_val), "\u2014",
                  fifelse(!is.finite(or_val), "\u221e", sprintf("%.2f", or_val)))
  lo_s <- fifelse(is.na(lo), "\u2014",
                  fifelse(!is.finite(lo), "0.00", sprintf("%.2f", lo)))
  hi_s <- fifelse(is.na(hi), "\u2014",
                  fifelse(!is.finite(hi), "\u221e", sprintf("%.2f", hi)))
  fifelse(is.na(or_val), "\u2014",
          paste0(or_s, " (", lo_s, "\u2013", hi_s, ")"))
}

out_sig <- sig_dt[, .(
  Community         = community,
  `ICD-10 chapter`  = chapter_labels[chapter],
  `n observed`      = n_in_comm_ch,
  `Community size`  = n_in_comm,
  `% in community`  = pct_in_comm,
  `% in network`    = pct_in_network,
  `O/E`             = obs_exp_ratio,
  `OR (95% CI)`     = fmt_or_ci(odds_ratio, ci_low, ci_high),
  Direction         = fifelse(direction == "enriched", "Enriched", "Depleted"),
  FDR               = signif(fdr, 3)
)]
setorder(out_sig, Community, -`O/E`)

# Sheet 2: full results
full_dt <- enrich_dt[, .(
  Community         = community,
  `ICD-10 chapter`  = chapter_labels[chapter],
  `n observed`      = n_in_comm_ch,
  `Community size`  = n_in_comm,
  `Chapter size`    = n_in_chapter,
  `Network size`    = n_total,
  `% in community`  = pct_in_comm,
  `% in network`    = pct_in_network,
  `O/E`             = obs_exp_ratio,
  `log2(O/E)`       = round(log2_oe, 3),
  OR                = odds_ratio,
  `CI lower`        = ci_low,
  `CI upper`        = ci_high,
  `P value`         = signif(p_value, 4),
  FDR               = signif(fdr, 4),
  `FDR < 0.05`      = significant
)]
setorder(full_dt, Community, `ICD-10 chapter`)

# Sheet 3: one-line summary per community
enriched_only <- sig_dt[direction == "enriched"]
summary_dt <- enriched_only[, .(
  N_enriched_chapters  = .N,
  `Enriched chapters`  = paste(chapter_labels[chapter], collapse = "; "),
  `Max O/E`            = max(obs_exp_ratio, na.rm = TRUE)
), by = .(Community = community)]
# Add communities with zero significant enrichment
all_comm_tbl <- data.table(Community = all_comms)
summary_dt <- merge(all_comm_tbl, summary_dt, by = "Community", all.x = TRUE)
summary_dt[is.na(N_enriched_chapters), N_enriched_chapters := 0L]
summary_dt[is.na(`Enriched chapters`), `Enriched chapters` := "\u2014"]
summary_dt[is.na(`Max O/E`), `Max O/E` := NA_real_]
setorder(summary_dt, Community)

# Write workbook
wb <- createWorkbook()
hdr_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF",
                         halign = "center", valign = "center",
                         fgFill = "#4472C4", textDecoration = "bold")

addWorksheet(wb, "Significant_Enrichment")
writeDataTable(wb, "Significant_Enrichment", out_sig, withFilter = TRUE)
setColWidths(wb, "Significant_Enrichment", cols = 1:ncol(out_sig), widths = "auto")
addStyle(wb, "Significant_Enrichment", hdr_style, rows = 1,
         cols = 1:ncol(out_sig), gridExpand = TRUE)

# Row colouring: enriched = pink tint, depleted = blue tint
enrich_style  <- createStyle(fgFill = "#FCE4EC")
deplete_style <- createStyle(fgFill = "#E3F2FD")
for (r in seq_len(nrow(out_sig))) {
  sty <- if (out_sig[r, Direction] == "Enriched") enrich_style else deplete_style
  addStyle(wb, "Significant_Enrichment", sty,
           rows = r + 1, cols = 1:ncol(out_sig), gridExpand = TRUE, stack = TRUE)
}

addWorksheet(wb, "Full_Results")
writeDataTable(wb, "Full_Results", full_dt, withFilter = TRUE)
setColWidths(wb, "Full_Results", cols = 1:ncol(full_dt), widths = "auto")
addStyle(wb, "Full_Results", hdr_style, rows = 1,
         cols = 1:ncol(full_dt), gridExpand = TRUE)

addWorksheet(wb, "Community_Summary")
writeDataTable(wb, "Community_Summary", summary_dt, withFilter = TRUE)
setColWidths(wb, "Community_Summary", cols = 1:ncol(summary_dt), widths = "auto")
addStyle(wb, "Community_Summary", hdr_style, rows = 1,
         cols = 1:ncol(summary_dt), gridExpand = TRUE)

xlsx_file <- file.path(OUT_DIR, "TableS6_ICD10_enrichment_full_results.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", basename(xlsx_file), "\n")

# =========================
# 6) Console summary
# =========================
cat("\n============================================================\n")
cat("Step 9: ICD-10 chapter enrichment — complete.\n")
cat("============================================================\n")
cat("\nCommunity-level enrichment summary (enriched chapters only):\n")
for (cid in all_comms) {
  sig_ch <- enrich_dt[community == cid & significant == TRUE & direction == "enriched"]
  n_size <- unique(enrich_dt[community == cid, n_in_comm])
  if (nrow(sig_ch) > 0) {
    setorder(sig_ch, -obs_exp_ratio)
    ch_str <- paste(sprintf("Ch.%s (O/E=%.1f)",
                            sig_ch$chapter, sig_ch$obs_exp_ratio),
                    collapse = ", ")
    cat(sprintf("  C%-2d (n=%3d): %s\n", cid, n_size, ch_str))
  } else {
    cat(sprintf("  C%-2d (n=%3d): no enriched chapter\n", cid, n_size))
  }
}

cat("\nFiles written to:\n  ", OUT_DIR, "\n", sep = "")
cat("  1) ICD10_chapter_enrichment_heatmap.pdf\n")
cat("  2) Table_ICD10_enrichment_full_results.xlsx\n")
cat("     - Sheet 1: Significant_Enrichment (for manuscript)\n")
cat("     - Sheet 2: Full_Results          (for reviewers)\n")
cat("     - Sheet 3: Community_Summary     (one-line per community)\n")
cat("============================================================\n")