#!/usr/bin/env Rscript
# ==============================================================================
# Step 13: Top1-share distribution with τ = 0.40 threshold justification
#
# Purpose:
#   Visualise the empirical distribution of each patient's top-community share
#   (top1_share = largest community score / total score), and show that the
#   τ = 0.40 cut-point separates clearly dominant patterns from diffuse
#   (mixed) multimorbidity. This provides data-driven support for the
#   dominance threshold used in Step 4.5.
#
# Input:
#   - cox_dataset_with_dominant_community.csv  (from Step 4.5)
#     Required columns: dom_share_unw, total_score_unw, dom_community_unw
#
# Output:
#   - FigS6_top1_share_distribution_tau040.pdf
#
# Design:
#   Two-panel layout (portrait A5):
#     Panel A: histogram of top1_share for patients with total_score_unw > 0
#              (i.e., patients with at least one community disease)
#     Panel B: cumulative distribution function (CDF) of top1_share
#
#   Both panels show a vertical dashed line at τ = 0.40; annotated with
#   the proportion of patients classified as DOM vs MIXED.
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# =========================
# 0) Paths & parameters
# =========================
COX_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TAU <- 0.40   # dominance threshold (must match Step 4.5)

FONT_FAMILY <- "Arial"
PDF_W <- 9
PDF_H <- 5

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Load data
# =========================
cat("==> Step 1) Loading Cox dataset with dominant community assignments...\n")

if (!file.exists(COX_CSV)) {
  stop("Input file not found: ", COX_CSV, "\nRun Step 4.5 first.", call. = FALSE)
}

dt <- fread(COX_CSV,
            select = c("dom_share_unw", "total_score_unw", "dom_community_unw"))
cat("   Loaded:", nrow(dt), "patients\n")

# Sanity checks
req_cols <- c("dom_share_unw", "total_score_unw", "dom_community_unw")
missing <- setdiff(req_cols, names(dt))
if (length(missing) > 0) {
  stop("Required columns missing: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Subset to patients WITH community diseases (NONE group excluded; they have no share)
scored <- dt[total_score_unw > 0]
cat("   Patients with at least one community disease (score > 0):", nrow(scored), "\n")
cat("   Excluded NONE group (score = 0):", nrow(dt) - nrow(scored), "\n")

# =========================
# 2) Descriptive statistics
# =========================
cat("\n==> Step 2) Computing descriptive statistics...\n")

n_dom   <- scored[dom_share_unw >= TAU, .N]
n_mixed <- scored[dom_share_unw <  TAU, .N]
n_total <- nrow(scored)

pct_dom   <- n_dom   / n_total * 100
pct_mixed <- n_mixed / n_total * 100

cat(sprintf("   DOM (top1_share >= %.2f):   %d (%.1f%%)\n", TAU, n_dom, pct_dom))
cat(sprintf("   MIXED (top1_share < %.2f):  %d (%.1f%%)\n", TAU, n_mixed, pct_mixed))

# Distribution summary
q <- quantile(scored$dom_share_unw, probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
cat(sprintf("   Distribution (10/25/50/75/90): %.3f / %.3f / %.3f / %.3f / %.3f\n",
            q[1], q[2], q[3], q[4], q[5]))

# =========================
# 3) Panel A: histogram with density overlay
# =========================
cat("\n==> Step 3) Building histogram (Panel A)...\n")

# Annotation positions
y_max_approx <- max(hist(scored$dom_share_unw, breaks = seq(0, 1, 0.02),
                         plot = FALSE)$counts) * 1.10

p_hist <- ggplot(scored, aes(x = dom_share_unw)) +
  # Histogram: fill color distinguishes MIXED (blue, left) vs DOM (red, right)
  geom_histogram(data = scored[dom_share_unw <  TAU],
                 breaks = seq(0, TAU, by = 0.02),
                 fill = "#4575B4", color = "white", linewidth = 0.15,
                 boundary = 0) +
  geom_histogram(data = scored[dom_share_unw >= TAU],
                 breaks = seq(TAU, 1.0001, by = 0.02),
                 fill = "#D73027", color = "white", linewidth = 0.15,
                 boundary = TAU) +
  # Threshold line
  geom_vline(xintercept = TAU, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  # Threshold annotation
  annotate("text", x = TAU, y = y_max_approx * 0.96,
           label = sprintf("\u03C4 = %.2f", TAU),
           hjust = -0.15, vjust = 1,
           size = 4, family = FONT_FAMILY, fontface = "bold") +
  # MIXED / DOM proportions
  annotate("text", x = TAU / 2, y = y_max_approx * 0.80,
           label = sprintf("MIXED\n%d (%.1f%%)", n_mixed, pct_mixed),
           color = "#2C4F7C", size = 4.3, family = FONT_FAMILY,
           fontface = "bold", lineheight = 0.95) +
  annotate("text", x = (TAU + 1) / 2, y = y_max_approx * 0.80,
           label = sprintf("DOM\n%d (%.1f%%)", n_dom, pct_dom),
           color = "#8B1A1F", size = 4.3, family = FONT_FAMILY,
           fontface = "bold", lineheight = 0.95) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "(A) Patient-level top-community share distribution",
    x     = "Top1-share  =  largest community score  /  total score",
    y     = "Number of patients"
  ) +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0, color = "black"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# =========================
# 4) Panel B: empirical CDF
# =========================
cat("==> Step 4) Building cumulative distribution (Panel B)...\n")

p_cdf <- ggplot(scored, aes(x = dom_share_unw)) +
  stat_ecdf(geom = "step", color = "grey15", linewidth = 0.7) +
  geom_vline(xintercept = TAU, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  # CDF at tau (= proportion classified as MIXED)
  geom_hline(yintercept = pct_mixed / 100, linetype = "dotted",
             color = "grey40", linewidth = 0.4) +
  annotate("point", x = TAU, y = pct_mixed / 100,
           color = "black", size = 2.3) +
  annotate("text",  x = TAU + 0.03, y = pct_mixed / 100 - 0.04,
           label = sprintf("F(\u03C4) = %.3f\n(= MIXED fraction)",
                           pct_mixed / 100),
           hjust = 0, size = 3.6, family = FONT_FAMILY, color = "grey20") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "(B) Cumulative distribution",
    x     = "Top1-share",
    y     = "Cumulative proportion of patients"
  ) +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0, color = "black"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# =========================
# 5) Assemble and save
# =========================
cat("\n==> Step 5) Assembling two-panel figure...\n")

p_combined <- p_hist + p_cdf +
  plot_layout(ncol = 2, widths = c(1.1, 1)) +
  plot_annotation(
    caption = sprintf(paste0(
      "Distribution of top-community share for %s patients with at least one ",
      "community-assigned disease. Patients with top1-share \u2265 %.2f are classified ",
      "as having a dominant community phenotype (DOM, red); those below are MIXED (blue)."),
      format(n_total, big.mark = ","), TAU),
    theme = theme(plot.caption = element_text(family = FONT_FAMILY,
                                              hjust = 0, size = 9,
                                              color = "grey25",
                                              margin = margin(t = 8)))
  )

pdf_file <- file.path(OUT_DIR, "FigS6_top1_share_distribution_tau040.pdf")
ggsave(pdf_file, p_combined,
       width = PDF_W, height = PDF_H, units = "in", device = PDF_DEVICE)
cat("   Saved:", basename(pdf_file), "\n")

# =========================
# 6) Console summary
# =========================
cat("\n============================================================\n")
cat("Step 13: Top1-share distribution — complete.\n")
cat("============================================================\n\n")
cat("Key takeaway: how justifiable is \u03C4 = 0.40?\n\n")
cat(sprintf("  Patients with score > 0:      %s\n", format(n_total, big.mark = ",")))
cat(sprintf("  Classified as DOM:            %s (%.1f%%)\n",
            format(n_dom, big.mark = ","), pct_dom))
cat(sprintf("  Classified as MIXED:          %s (%.1f%%)\n",
            format(n_mixed, big.mark = ","), pct_mixed))
cat(sprintf("  Median top1-share:            %.3f\n", q[3]))
cat(sprintf("  IQR top1-share:               %.3f - %.3f\n", q[2], q[4]))
cat("\nInspect the histogram in Panel A:\n")
cat("  - A bimodal shape (dip around \u03C4) supports threshold choice.\n")
cat("  - A smooth unimodal shape suggests \u03C4 is pragmatic, not data-inherent,\n")
cat("    which should be acknowledged in the Methods or Limitations.\n")
cat("\nOutput:\n  ", pdf_file, "\n", sep = "")
cat("============================================================\n")