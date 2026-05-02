#!/usr/bin/env Rscript
# ==============================================================================
# Step 14: STROBE Flow Diagram for Figure 1A (Study Population)
#
# Purpose:
#   Produce a publication-grade STROBE flow diagram showing cohort assembly
#   from raw EHR to analytical cohort of 640 637 CKM-relevant adults.
#
# Output:
#   - Fig1A_STROBE_flow.pdf
#
# PLEASE UPDATE THE NUMBERS BELOW TO MATCH YOUR ACTUAL PIPELINE CHECKPOINTS
# (placeholders marked with <TO VERIFY> — run grep over your Step 1–4 logs to
# find the corresponding exclusion counts).
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

# =========================
# 0) Paths & device
# =========================
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\main_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"
PDF_W <- 7
PDF_H <- 9

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) COHORT NUMBERS — please edit to match actual pipeline
# =========================
#
# Structure: (top) raw EHR → (filter 1) encounters in period →
#            (filter 2) CKM-relevant via ICD-10 → (filter 3) ≥ 1 baseline
#            diagnosis record → (filter 4) death/follow-up linkage →
#            (bottom) analytical cohort of 640 637
#
# You will need to fill in the actual upstream numbers from Step 1–4 logs.

flow <- list(
  # MAIN BOXES (left column)
  list(label = "Yichang regional EHR, 2016–2024\n(N = <TO VERIFY> unique adults)",
       y = 10, fill = "#E8F0F8"),
  
  list(label = sprintf(paste0("Adults with \u2265 1 ICD-10 diagnosis in the CKM ",
                              "domain* during 2016–2024\n(N = <TO VERIFY>)")),
       y = 8,  fill = "#E8F0F8"),
  
  list(label = sprintf(paste0("Adults with \u2265 1 diagnostic record during ",
                              "the baseline window (2016-01-01 to 2020-12-31)\n",
                              "(N = <TO VERIFY>)")),
       y = 6,  fill = "#E8F0F8"),
  
  list(label = sprintf(paste0("Alive at landmark index date\n",
                              "(2021-01-01; N = <TO VERIFY>)")),
       y = 4,  fill = "#E8F0F8"),
  
  list(label = paste0("**Analytical cohort**\n",
                      "N = 640,637\n",
                      "Events = 38,670 (4-y all-cause mortality)\n",
                      "Median follow-up = 4.00 years"),
       y = 2, fill = "#D0E2F0")
)

# EXCLUSION BOXES (right column, between each main box)
exclusions <- list(
  list(label = "Excluded:\nNo ICD-10 diagnosis in C / K / M domains\n(n = <TO VERIFY>)",
       y = 9),
  list(label = "Excluded:\nNo diagnostic activity in baseline window\n(n = <TO VERIFY>)",
       y = 7),
  list(label = "Excluded:\nDied before landmark (2021-01-01)\n(n = <TO VERIFY>)",
       y = 5),
  list(label = "Excluded:\nPatients with no community-assignable diseases\n(NONE group; n = 372)",
       y = 3)
)

# =========================
# 2) Build diagram with grid
# =========================
cat("==> Step 2) Building STROBE flow diagram...\n")

# Open PDF
pdf_file <- file.path(OUT_DIR, "Fig1A_STROBE_flow.pdf")
PDF_DEVICE(pdf_file, width = PDF_W, height = PDF_H, family = FONT_FAMILY)

grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))

# Layout coordinates: main column centered at x_main, exclusion column at x_excl
x_main <- 0.33
x_excl <- 0.78
main_box_w <- 0.58
excl_box_w <- 0.38

# Helper: draw a box with text
draw_box <- function(x_center, y_center, width, label,
                     fill = "#E8F0F8", border = "grey30",
                     fontface = "plain", fontsize = 9) {
  height <- 0.075
  grid.roundrect(x = x_center, y = y_center,
                 width = width, height = height,
                 gp = gpar(fill = fill, col = border, lwd = 0.8),
                 r = unit(2, "mm"))
  # Handle markdown-style **bold** by splitting
  if (grepl("\\*\\*", label)) {
    label <- gsub("\\*\\*", "", label)
    fontface <- "bold"
  }
  grid.text(label, x = x_center, y = y_center,
            gp = gpar(fontface = fontface, fontsize = fontsize,
                      fontfamily = FONT_FAMILY, lineheight = 1.2))
}

# Helper: draw vertical connecting arrow
draw_arrow <- function(x, y_from, y_to, color = "grey40") {
  grid.lines(x = c(x, x), y = c(y_from, y_to),
             gp = gpar(col = color, lwd = 1.2),
             arrow = arrow(length = unit(2.5, "mm"), type = "closed"))
}

# Helper: horizontal exclusion arrow
draw_excl_arrow <- function(x_from, x_to, y, color = "grey40") {
  grid.lines(x = c(x_from, x_to), y = c(y, y),
             gp = gpar(col = color, lwd = 1.2),
             arrow = arrow(length = unit(2.5, "mm"), type = "closed"))
}

# Coordinate system: y 1–11 (normalised to 0–1)
y_to_npc <- function(y) (y - 0.5) / 11

# Draw main boxes
for (item in flow) {
  draw_box(x_main, y_to_npc(item$y),
           width = main_box_w,
           label = item$label,
           fill = item$fill,
           fontsize = 9)
}

# Draw main vertical arrows (between adjacent main boxes)
main_ys <- sapply(flow, function(x) x$y)
for (i in seq_len(length(main_ys) - 1)) {
  y_from <- y_to_npc(main_ys[i]   - 0.6)  # bottom of upper box
  y_to   <- y_to_npc(main_ys[i+1] + 0.6)  # top of lower box
  draw_arrow(x_main, y_from, y_to)
}

# Draw exclusion boxes + horizontal arrows
for (i in seq_along(exclusions)) {
  item <- exclusions[[i]]
  draw_box(x_excl, y_to_npc(item$y),
           width = excl_box_w,
           label = item$label,
           fill = "#FFF5E6",
           fontsize = 8)
  # Arrow from main column to exclusion box
  draw_excl_arrow(x_main + main_box_w/2,
                  x_excl - excl_box_w/2,
                  y_to_npc(item$y))
}

# Title
grid.text("Figure 1A. Study population flow (STROBE)",
          x = 0.5, y = 0.97,
          gp = gpar(fontsize = 12, fontface = "bold",
                    fontfamily = FONT_FAMILY))

# Footnote
grid.text(paste0("*CKM domain diagnoses defined per Supplementary Table S1 ",
                 "(ICD-10 3-character codes)."),
          x = 0.05, y = 0.03, just = "left",
          gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY,
                    col = "grey30"))

popViewport()
dev.off()

cat("   Saved:", pdf_file, "\n")

# =========================
# 3) Export editable data (YAML / CSV) for manual editing
# =========================
cat("\n==> Step 3) Exporting editable flow data (so reviewers/Weihao can edit)...\n")

flow_dt <- data.frame(
  Step       = 1:length(flow),
  Label      = sapply(flow, function(x) gsub("\\n", " | ", gsub("\\*\\*", "", x$label))),
  Y_position = sapply(flow, function(x) x$y),
  Type       = "main"
)

excl_dt <- data.frame(
  Step       = seq_along(exclusions),
  Label      = sapply(exclusions, function(x) gsub("\\n", " | ", x$label)),
  Y_position = sapply(exclusions, function(x) x$y),
  Type       = "exclusion"
)

all_dt <- rbind(flow_dt, excl_dt)
csv_file <- file.path(OUT_DIR, "Fig1A_STROBE_flow_data.csv")
write.csv(all_dt, csv_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("   Saved editable data:", csv_file, "\n")
cat("   (Edit the CSV to update exclusion counts, then re-run this script.)\n")

cat("\n============================================================\n")
cat("Step 14: STROBE flow diagram (Figure 1A) — complete.\n")
cat("============================================================\n\n")
cat("IMPORTANT: the diagram contains several '<TO VERIFY>' placeholders.\n")
cat("Replace them with actual exclusion counts from your Step 1-4 logs:\n\n")
cat("  - Raw EHR adult population (N_total)\n")
cat("  - # with no CKM ICD-10 codes\n")
cat("  - # with no diagnostic activity in 2016-2020 baseline\n")
cat("  - # who died before 2021-01-01 landmark\n")
cat("  - # NONE group (= 372 per current analysis)\n")
cat("\nYou can edit:\n")
cat("  A) Directly in this script (flow/exclusions lists at top), OR\n")
cat("  B) Edit the CSV:", csv_file, "\n")
cat("\n============================================================\n")