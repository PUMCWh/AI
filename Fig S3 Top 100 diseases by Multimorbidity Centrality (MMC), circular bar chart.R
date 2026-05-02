#!/usr/bin/env Rscript
# ==============================================================================
# Fig S3: Top 100 diseases by Multimorbidity Centrality (MMC), circular bar chart
#
# Single panel adaptation of user's PhD-style reference plot
# (uploaded as 1777463500351_MMC值绘制参考.R, 3-panel Overall/Female/Male)
#
# Key design points:
#   - Polar bar chart, 100 bars sorted by ICD-10 chapter then MMC (within chapter)
#   - Bar color = ICD-10 chapter (14 chapters, distinct palette)
#   - Disease codes labelled on outer rim, rotated to follow circle
#   - Inner radial axis with MMC scale (0.5, 1, 1.5, 2, 2.5, 3)
#   - Use ICD10_Chapter column DIRECTLY (no need for case_when on code prefix)
#
# Input:
#   - MMC_coef0.01_alpha0.05.csv (DiseaseCode, ICD10_Chapter, Frequency,
#     Prevalence, MMC)
#
# Output:
#   - FigS3_MMC_top100_circular.pdf
# ==============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# =========================
# Paths
# =========================
# MMC 文件在 NETWORK_PCOR_<timestamp> 子文件夹
# 用 list.files 自动找最新的 (避免 hardcode timestamp)
PAPER_ROOT  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL"
pcor_dirs <- list.files(PAPER_ROOT, pattern = "^NETWORK_PCOR_",
                        full.names = TRUE)
if (length(pcor_dirs) == 0) {
  stop("No NETWORK_PCOR_* folder found in ", PAPER_ROOT)
}
PCOR_DIR <- sort(pcor_dirs, decreasing = TRUE)[1]   # latest
MMC_CSV  <- file.path(PCOR_DIR, "MMC_coef0.01_alpha0.05.csv")
cat(sprintf("Using MMC file: %s\n", MMC_CSV))

OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# Parameters
# =========================
TOP_N           <- 100
INNER_RADIUS_PROP <- 1
LABEL_SIZE      <- 3.2
AXIS_LABEL_SIZE <- 4
GRID_COLOR      <- "grey60"
AXIS_COLOR      <- "grey90"
START_PADDING   <- 1
END_PADDING     <- 5

# ICD-10 chapter color palette (与您参考代码一致, 14 chapters)
chapter_cols <- c(
  "1"  = "#B70031", "2"  = "#EC6A3D", "3"  = "#A1C6AD", "4"  = "#0098C5",
  "5"  = "#9A6CAC", "6"  = "#004B9B", "7"  = "#C0C7CB", "8"  = "#AAD9EA",
  "9"  = "#FFE100", "10" = "#E71A10", "11" = "#6EB327", "12" = "#8E006C",
  "13" = "#D2ACC9", "14" = "#ADB5DA"
)

# Chapter labels — Roman numerals only (per user request; chapter names omitted)
chapter_labels <- c(
  "1"  = "I",   "2"  = "II",   "3"  = "III",  "4"  = "IV",
  "5"  = "V",   "6"  = "VI",   "7"  = "VII",  "8"  = "VIII",
  "9"  = "IX",  "10" = "X",    "11" = "XI",   "12" = "XII",
  "13" = "XIII", "14" = "XIV"
)

# =========================
# 1) Load MMC data
# =========================
cat("==> Step 1) Loading MMC data...\n")
mmc_dt <- fread(MMC_CSV)
cat(sprintf("   Total disease nodes: %d\n", nrow(mmc_dt)))
cat(sprintf("   Columns: %s\n", paste(names(mmc_dt), collapse = ", ")))

# =========================
# 2) Top 100 + sort by chapter then MMC
# =========================
cat("\n==> Step 2) Selecting top 100 diseases by MMC...\n")

mmc_dt[, MMC := as.numeric(MMC)]
mmc_dt[, ICD10_Chapter := as.character(ICD10_Chapter)]

mmc_top <- mmc_dt[!is.na(MMC) & MMC >= 0]
setorder(mmc_top, -MMC)
mmc_top <- mmc_top[1:TOP_N]

cat(sprintf("   Top 100 MMC range: %.4f to %.4f\n",
            min(mmc_top$MMC), max(mmc_top$MMC)))

# Sort within: chapter ascending, MMC descending within chapter
# Chapter "0" (unknown) goes to the end
mmc_top[, chapter_factor := factor(ICD10_Chapter,
                                   levels = names(chapter_cols))]
mmc_top[, chapter_sort_key := fifelse(ICD10_Chapter == "0",
                                      99L,
                                      as.integer(ICD10_Chapter))]
setorder(mmc_top, chapter_sort_key, -MMC)
mmc_top[, id := seq_len(.N) - 1L]
mmc_top[, x_pos := id + START_PADDING]

# Add invisible dummy rows for each chapter (ensures all 14 appear in legend
# even if a chapter has no diseases in top 100, e.g. Ch VIII Ear)
# These bars have MMC=0 and are placed at x_pos = -10 (off-canvas), so
# they are not visible in the polar chart but populate the fill scale.
dummy <- data.table(
  DiseaseCode      = paste0("dummy_", names(chapter_cols)),
  ICD10_Chapter    = names(chapter_cols),
  Frequency        = 0L,
  Prevalence       = 0,
  MMC              = 0,
  chapter_factor   = factor(names(chapter_cols),
                            levels = names(chapter_cols)),
  chapter_sort_key = as.integer(names(chapter_cols)),
  id               = -1L,
  x_pos            = -10
)
# Only keep columns present in mmc_top (order match for rbind)
common_cols <- intersect(names(dummy), names(mmc_top))
mmc_top_with_dummy <- rbind(mmc_top[, ..common_cols], dummy[, ..common_cols],
                            fill = TRUE)

# Show chapter distribution (real top 100 only)
chap_dist <- mmc_top[, .N, by = chapter_factor]
cat("\n   Chapter distribution in top 100:\n")
for (i in seq_len(nrow(chap_dist))) {
  ch <- as.character(chap_dist$chapter_factor[i])
  ch_label <- ifelse(ch %in% names(chapter_labels),
                     chapter_labels[ch], "Other")
  cat(sprintf("     Ch %-3s %-22s n = %d\n", ch, ch_label, chap_dist$N[i]))
}

# =========================
# 3) Compute polar coordinates
# =========================
y_max         <- max(mmc_top$MMC, na.rm = TRUE)
inner_radius  <- y_max * INNER_RADIUS_PROP
mmc_breaks    <- c(0.5, 1, 1.5, 2, 2.5)
radial_breaks <- inner_radius + mmc_breaks
plot_x_lim    <- c(0, (TOP_N - 1) + START_PADDING + END_PADDING)
offset        <- y_max * 0.05

# Label data: angle calculation for circular text
label_data <- copy(mmc_top)
label_data[, label_y := inner_radius + MMC + offset]
label_data[, angle   := 90 - (id + 0.5) / .N * 360]
label_data[, hjust   := fifelse(angle > -90 & angle < 90, 0, 1)]
label_data[, angle   := fifelse(angle < -90, angle + 180, angle)]

y_limit_max <- (inner_radius + max(mmc_breaks)) * 1.05

# =========================
# 4) Plot
# =========================
cat("\n==> Step 3) Plotting circular bar chart...\n")

p <- ggplot(mmc_top_with_dummy) +
  # Bars (real ones from mmc_top; dummies have MMC=0 -> invisible)
  geom_rect(aes(xmin = x_pos - 0.5,
                xmax = x_pos + 0.5,
                ymin = inner_radius,
                ymax = inner_radius + MMC,
                fill = chapter_factor),
            colour = "white", linewidth = 0.1) +
  # Inner radial axis line
  annotate("segment", x = 0, xend = 0,
           y = inner_radius,
           yend = inner_radius + max(mmc_breaks),
           colour = AXIS_COLOR, linewidth = 0.6) +
  # Radial axis labels (MMC values: 0.5, 1.0, 1.5, 2.0, 2.5)
  lapply(seq_along(radial_breaks), function(i) {
    annotate("text", x = 0, y = radial_breaks[i],
             label = mmc_breaks[i],
             colour = "black", size = AXIS_LABEL_SIZE,
             hjust = 1.1, vjust = 0.5)
  }) +
  # Disease code labels on outer rim
  geom_text(data = label_data,
            aes(x = x_pos, y = label_y, label = DiseaseCode,
                angle = angle, hjust = hjust),
            size = LABEL_SIZE, family = FONT_FAMILY) +
  # Polar transform
  coord_polar(theta = "x", start = 0, direction = 1, clip = "off") +
  scale_x_continuous(limits = plot_x_lim) +
  scale_y_continuous(limits = c(0, y_limit_max),
                     breaks = radial_breaks,
                     labels = NULL) +
  # Color scale: chapter -> color, with all 14 chapters always shown in legend
  scale_fill_manual(values = chapter_cols,
                    labels = chapter_labels,
                    breaks = names(chapter_cols),    # force 14 entries even if missing in data
                    name   = "ICD-10 chapter",
                    drop   = FALSE,
                    guide  = guide_legend(
                      ncol = 1,
                      keywidth  = unit(0.6, "cm"),
                      keyheight = unit(0.6, "cm"),
                      title.theme = element_text(face = "bold",
                                                 size = 13,
                                                 family = FONT_FAMILY),
                      label.theme = element_text(size = 11,
                                                 family = FONT_FAMILY)
                    )) +
  # Theme: minimal, transparent grids only on radial axis
  theme_minimal(base_size = 11, base_family = FONT_FAMILY) +
  theme(
    axis.text          = element_blank(),
    axis.title         = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = GRID_COLOR, linewidth = 0.3),
    plot.margin        = margin(10, 10, 10, 10),
    legend.position    = "right",
    legend.box.margin  = margin(l = 10)
  )

# =========================
# 5) Save
# =========================
pdf_file <- file.path(OUT_DIR, "FigS3_MMC_top100_circular.pdf")
ggsave(pdf_file, p,
       width = 11, height = 8.5, units = "in",
       device = PDF_DEVICE, dpi = 300)

cat(sprintf("\n   \u2713 Saved: %s\n", pdf_file))

# =========================
# 6) Console report
# =========================
cat("\n=========================================================\n")
cat("Fig S3: Top 100 MMC circular bar chart - DONE\n")
cat("=========================================================\n\n")
cat(sprintf("Top 5 diseases by MMC:\n"))
for (i in 1:5) {
  ch <- mmc_top$ICD10_Chapter[i]
  ch_label <- ifelse(ch %in% names(chapter_labels),
                     chapter_labels[ch], "Other")
  cat(sprintf("  %2d. %-5s MMC = %.4f  (Ch %s: %s)\n",
              i, mmc_top$DiseaseCode[i], mmc_top$MMC[i],
              ch, ch_label))
}
cat(sprintf("\nFile: %s\n", pdf_file))