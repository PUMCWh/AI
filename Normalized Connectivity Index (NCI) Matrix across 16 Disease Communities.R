#!/usr/bin/env Rscript
# ==============================================================================
# Step 11: Normalized Connectivity Index (NCI) Matrix across 16 Disease Communities
#
# Purpose:
#   Quantify how strongly each pair of disease communities is connected at the
#   network level, after accounting for community size differences (larger
#   communities trivially have more inter-community edges; NCI corrects for this).
#
#   NCI is used to show, for example, whether the cardiovascular and renal
#   communities in CKM are more tightly interlinked than renal and malignancy.
#
# Input:
#   - medoid_edges_with_community_coef0_01_alpha0_05.csv
#     (columns: from, to, weight, p_value, p_adj, from_community, to_community, edge_type)
#
# Outputs (written to OUT_DIR):
#   - FigS5_NCI_community_connectivity_heatmap.pdf   — normalized (main)
#   - FigS5b_raw_total_weight_heatmap.pdf            — raw Σφ, for reference
#   - TableS5b_NCI_pairs_ranked.xlsx                 — top pairs ranked by NCI
#
# Method:
#   For communities i and j, with n_i and n_j members respectively:
#     raw_ij  = Σφ over all inter-community edges (from community i to j)
#     NCI_ij  = raw_ij / (n_i × n_j)
#
#   This normalises by the maximum number of possible inter-community edges,
#   yielding a quantity equivalent to "weighted inter-community edge density".
#   Diagonal (i = j, intra-community) is masked to focus on cross-community links.
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
EDGES_CSV  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\community_detection\\medoid_edges_with_community_coef0.01_alpha0.05.csv"
OUT_DIR    <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"
PDF_W       <- 9
PDF_H       <- 8

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Load edge list and derive community sizes
# =========================
cat("==> Step 1) Loading edges and deriving community membership...\n")

if (!file.exists(EDGES_CSV)) {
  stop("Edges CSV not found: ", EDGES_CSV, call. = FALSE)
}

edges_dt <- fread(EDGES_CSV)
cat("   Loaded:", nrow(edges_dt), "edges\n")
cat("   Columns:", paste(names(edges_dt), collapse = ", "), "\n")

# Validate required columns
req_cols <- c("from", "to", "weight", "from_community", "to_community")
missing <- setdiff(req_cols, names(edges_dt))
if (length(missing) > 0) {
  stop("Required columns missing: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

# Convert community columns to integer
edges_dt[, from_community := as.integer(from_community)]
edges_dt[, to_community   := as.integer(to_community)]

# Derive node → community lookup (two ways: from/to endpoints)
n2c_from <- unique(edges_dt[!is.na(from_community) & from_community > 0,
                            .(node = from, community = from_community)])
n2c_to   <- unique(edges_dt[!is.na(to_community) & to_community > 0,
                            .(node = to,   community = to_community)])
node2comm <- unique(rbind(n2c_from, n2c_to))

if (nrow(node2comm) != uniqueN(node2comm$node)) {
  warning("Some nodes appear in multiple communities. Taking first occurrence.")
  node2comm <- unique(node2comm, by = "node")
}

comm_size <- node2comm[, .(n_nodes = .N), by = community]
setorder(comm_size, community)
cat("   Communities detected:", nrow(comm_size), "\n")
cat("   Size range:", min(comm_size$n_nodes), "-", max(comm_size$n_nodes), "\n")

# =========================
# 2) Compute raw Σφ matrix and NCI matrix
# =========================
cat("\n==> Step 2) Computing raw inter-community weights and NCI...\n")

# Filter to edges with both endpoints in a valid community
edges_valid <- edges_dt[!is.na(from_community) & from_community > 0 &
                          !is.na(to_community)   & to_community   > 0]

# Aggregate Σφ per (from_community, to_community). Because the edge list is
# undirected (upper triangle), we symmetrise the matrix after aggregation.
agg <- edges_valid[, .(sum_phi = sum(weight, na.rm = TRUE)),
                   by = .(from_community, to_community)]

# Build matrix
comm_ids <- sort(unique(c(agg$from_community, agg$to_community)))
comm_labels <- as.character(comm_ids)
n_comms <- length(comm_ids)

raw_mat <- matrix(0, nrow = n_comms, ncol = n_comms,
                  dimnames = list(comm_labels, comm_labels))

for (r in seq_len(nrow(agg))) {
  i <- as.character(agg$from_community[r])
  j <- as.character(agg$to_community[r])
  raw_mat[i, j] <- raw_mat[i, j] + agg$sum_phi[r]
  if (i != j) raw_mat[j, i] <- raw_mat[j, i] + agg$sum_phi[r]
}
# After symmetrising: if edges were already bidirectional, we will have
# double-counted. Edge list format (undirected upper triangle) stores each
# edge once, so this symmetrisation is correct. If your data stores both
# directions, divide by 2 at this point.

# Size vector aligned to matrix dimnames
size_vec <- setNames(comm_size$n_nodes[match(comm_ids, comm_size$community)],
                     comm_labels)

# NCI matrix = raw / (n_i * n_j)
denom_mat <- outer(size_vec, size_vec, "*")
nci_mat   <- raw_mat / denom_mat

# Mask diagonal (intra-community, not the focus here)
diag(nci_mat) <- NA_real_
diag_raw      <- diag(raw_mat)   # keep for reference
diag(raw_mat) <- NA_real_

# Build long data.table versions for plotting
dt_nci <- as.data.table(as.table(nci_mat))
setnames(dt_nci, c("from", "to", "nci"))
dt_nci[, nci := as.numeric(nci)]

dt_raw <- as.data.table(as.table(raw_mat))
setnames(dt_raw, c("from", "to", "raw_sum_phi"))
dt_raw[, raw_sum_phi := as.numeric(raw_sum_phi)]

dt_long <- merge(dt_nci, dt_raw, by = c("from", "to"))

# Off-diagonal stats for console
dt_inter <- dt_long[from != to]
cat("   Off-diagonal NCI range: ",
    sprintf("%.5f to %.5f\n", min(dt_inter$nci, na.rm = TRUE),
            max(dt_inter$nci, na.rm = TRUE)))

# =========================
# 3) Scaling for readable cell labels
# =========================
median_nci <- median(dt_inter[nci > 0, nci], na.rm = TRUE)
cat("   Median nonzero NCI:", sprintf("%.5f", median_nci), "\n")

if (is.na(median_nci) || median_nci < 0.001) {
  scale_factor <- 10000
  scale_label  <- expression("NCI  (\u00D7 10"^4*")")
  scale_suffix <- " (\u00D7 10^4)"
} else if (median_nci < 0.01) {
  scale_factor <- 1000
  scale_label  <- expression("NCI  (\u00D7 10"^3*")")
  scale_suffix <- " (\u00D7 10^3)"
} else {
  scale_factor <- 100
  scale_label  <- expression("NCI  (\u00D7 10"^2*")")
  scale_suffix <- " (\u00D7 10^2)"
}
cat("   Scale factor:", scale_factor, "\n")

# =========================
# 4) Build main NCI heatmap
# =========================
cat("\n==> Step 4) Building NCI heatmap...\n")

dt_plot <- copy(dt_long)
dt_plot[, scaled_nci := nci * scale_factor]
dt_plot[, cell_label := fifelse(!is.na(scaled_nci) & scaled_nci > 0,
                                sprintf("%.1f", scaled_nci), "")]
dt_plot[, from := factor(from, levels = comm_labels)]
dt_plot[, to   := factor(to,   levels = comm_labels)]

# Colour scale upper bound: round up to nearest 5
actual_max <- max(dt_plot$scaled_nci, na.rm = TRUE)
fill_max   <- ceiling(actual_max / 5) * 5
if (!is.finite(fill_max) || fill_max == 0) fill_max <- 5

# Axis label mapping: community ID → "C<id>"
lab_map <- setNames(paste0("C", comm_labels), comm_labels)

# Pretty legend breaks
candidates <- c(1, 2, 5, 10, 20)
step <- candidates[which(fill_max / candidates >= 3 &
                           fill_max / candidates <= 6)[1]]
if (is.na(step)) step <- fill_max / 4
legend_breaks <- seq(0, fill_max, by = step)

p_nci <- ggplot(dt_plot, aes(x = from, y = to)) +
  geom_tile(aes(fill = scaled_nci), color = "grey85", linewidth = 0.3) +
  geom_text(aes(label = cell_label), size = 3, color = "grey15",
            family = FONT_FAMILY) +
  scale_fill_gradient(
    low       = "white",
    high      = "#B2182B",
    limits    = c(0, fill_max),
    breaks    = legend_breaks,
    oob       = scales::squish,
    name      = scale_label,
    na.value  = "grey95"
  ) +
  scale_x_discrete(labels = lab_map) +
  scale_y_discrete(labels = lab_map) +
  labs(
    title    = "Normalized Connectivity Index (NCI) between CKM disease communities",
    subtitle = paste0("NCI = \u03A3\u03C6 / (n_i \u00D7 n_j);  ",
                      "diagonal (intra-community) shown in grey"),
    x = "Community",
    y = "Community"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text              = element_text(family = FONT_FAMILY),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 12, color = "black"),
    plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.text.x       = element_text(size = 11, color = "black"),
    axis.text.y       = element_text(size = 11, color = "black"),
    axis.title        = element_text(size = 12, color = "black"),
    legend.position   = "right",
    legend.key.height = unit(1.3, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.title      = element_text(size = 11, color = "black"),
    legend.text       = element_text(size = 10, color = "black"),
    panel.grid        = element_blank(),
    plot.margin       = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

pdf_file <- file.path(OUT_DIR, "FigS5_NCI_community_connectivity_heatmap.pdf")
ggsave(pdf_file, p_nci, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved:", basename(pdf_file), "\n")

# =========================
# 5) Bonus: raw Σφ heatmap for reference
# =========================
cat("\n==> Step 5) Building raw-weight (\u03A3\u03C6) reference heatmap...\n")

dt_raw_plot <- copy(dt_long)
dt_raw_plot[, from := factor(from, levels = comm_labels)]
dt_raw_plot[, to   := factor(to,   levels = comm_labels)]
dt_raw_plot[, cell_label := fifelse(!is.na(raw_sum_phi) & raw_sum_phi > 0,
                                    sprintf("%.2f", raw_sum_phi), "")]

p_raw <- ggplot(dt_raw_plot, aes(x = from, y = to)) +
  geom_tile(aes(fill = raw_sum_phi), color = "grey85", linewidth = 0.3) +
  geom_text(aes(label = cell_label), size = 2.6, color = "grey15",
            family = FONT_FAMILY) +
  scale_fill_gradient(
    low      = "white",
    high     = "#6A1C2B",
    name     = expression("Total " * Sigma * phi),
    na.value = "grey95"
  ) +
  scale_x_discrete(labels = lab_map) +
  scale_y_discrete(labels = lab_map) +
  labs(
    title    = "Raw inter-community weight (\u03A3\u03C6) — for reference",
    subtitle = "Not size-normalized: larger communities trivially show larger totals",
    x = "Community",
    y = "Community"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text              = element_text(family = FONT_FAMILY),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 12, color = "black"),
    plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.text.x       = element_text(size = 10, color = "black"),
    axis.text.y       = element_text(size = 10, color = "black"),
    axis.title        = element_text(size = 12, color = "black"),
    legend.position   = "right",
    legend.key.height = unit(1.3, "cm"),
    panel.grid        = element_blank()
  ) +
  coord_fixed()

pdf_file_raw <- file.path(OUT_DIR, "FigS5b_raw_total_weight_heatmap.pdf")
ggsave(pdf_file_raw, p_raw, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved:", basename(pdf_file_raw), "\n")

# =========================
# 6) Table S5b: NCI pairs ranked
# =========================
cat("\n==> Step 6) Writing Table S5b (NCI pairs ranked)...\n")

# Upper triangle only (each community pair once)
pairs_dt <- copy(dt_inter)
pairs_dt[, from_i := as.integer(as.character(from))]
pairs_dt[, to_i   := as.integer(as.character(to))]
pairs_dt <- pairs_dt[from_i < to_i]
setorder(pairs_dt, -nci)

pairs_out <- data.table(
  Rank               = seq_len(nrow(pairs_dt)),
  `Community i`      = paste0("C", pairs_dt$from_i),
  `Community j`      = paste0("C", pairs_dt$to_i),
  `n_i`              = size_vec[as.character(pairs_dt$from_i)],
  `n_j`              = size_vec[as.character(pairs_dt$to_i)],
  `Raw Σφ`           = round(pairs_dt$raw_sum_phi, 4),
  NCI                = round(pairs_dt$nci, 6),
  `NCI scaled`       = round(pairs_dt$nci * scale_factor, 2)
)
setnames(pairs_out, "NCI scaled",
         paste0("NCI", scale_suffix))

wb <- createWorkbook()
hdr_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF",
                         halign = "center", valign = "center",
                         fgFill = "#4472C4", textDecoration = "bold")
addWorksheet(wb, "NCI_pairs_ranked")
writeDataTable(wb, "NCI_pairs_ranked", pairs_out, withFilter = TRUE)
setColWidths(wb, "NCI_pairs_ranked", cols = 1:ncol(pairs_out), widths = "auto")
addStyle(wb, "NCI_pairs_ranked", hdr_style, rows = 1,
         cols = 1:ncol(pairs_out), gridExpand = TRUE)

# Community size reference sheet
size_out <- data.table(
  Community = paste0("C", comm_size$community),
  Size      = comm_size$n_nodes,
  `Intra-community Σφ` = round(diag_raw[as.character(comm_size$community)], 4)
)
addWorksheet(wb, "Community_sizes")
writeDataTable(wb, "Community_sizes", size_out, withFilter = TRUE)
setColWidths(wb, "Community_sizes", cols = 1:ncol(size_out), widths = "auto")
addStyle(wb, "Community_sizes", hdr_style, rows = 1,
         cols = 1:ncol(size_out), gridExpand = TRUE)

xlsx_file <- file.path(OUT_DIR, "TableS5b_NCI_pairs_ranked.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", basename(xlsx_file), "\n")

# =========================
# 7) Console summary
# =========================
cat("\n============================================================\n")
cat("Step 11: NCI community connectivity — complete.\n")
cat("============================================================\n")
cat("\nTop 10 most strongly connected community pairs (by NCI):\n")
for (r in seq_len(min(10, nrow(pairs_out)))) {
  cat(sprintf("   %2d. %s \u2013 %s  |  NCI = %.5f (scaled %.2f)  |  raw \u03A3\u03C6 = %.3f\n",
              pairs_out$Rank[r],
              pairs_out$`Community i`[r],
              pairs_out$`Community j`[r],
              pairs_dt$nci[r],
              pairs_dt$nci[r] * scale_factor,
              pairs_dt$raw_sum_phi[r]))
}

cat("\nFiles written to:\n  ", OUT_DIR, "\n", sep = "")
cat("  1) FigS5_NCI_community_connectivity_heatmap.pdf   (main)\n")
cat("  2) FigS5b_raw_total_weight_heatmap.pdf            (reference)\n")
cat("  3) TableS5b_NCI_pairs_ranked.xlsx\n")
cat("============================================================\n")