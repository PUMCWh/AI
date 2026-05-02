#!/usr/bin/env Rscript
# ==============================================================================
# Step 10: Node Role Classification (Guimerà-Amaral framework)
#
# Purpose:
#   Classify each disease node in the CKM comorbidity network into topological
#   roles based on two complementary metrics:
#     - Within-module z-score (z): hub vs non-hub relative to own community
#     - Participation coefficient (P): how connections distribute across
#       communities (0 = all internal, approaches 1 = evenly distributed)
#
#   This identifies "bridge diseases" — conditions that physiologically link
#   multiple disease clusters and may reflect shared pathophysiology or
#   common comorbidity pathways (e.g., dyslipidemia connecting metabolic and
#   cardiovascular clusters).
#
# Inputs:
#   - medoid_bridge_metrics_coef0_01_alpha0_05.csv
#     (columns: node, community, strength, within_module_z_score,
#      participation_coef, top1-3_ext_community, role, ICD10_Chapter, ...)
#
# Outputs (written to OUT_DIR):
#   - FigS4_node_role_scatter_zP.pdf       — main scatter plot
#   - TableS2_top20_bridge_diseases.xlsx   — top-20 by participation coef
#   - TableS2b_all_bridge_nodes.xlsx       — full list P≥0.30 with targets
#
# Role definitions (Guimerà-Amaral):
#   Provincial hub:  z ≥ 2.5  AND P < 0.30   (community-internal core)
#   Connector hub:   z ≥ 2.5  AND P ≥ 0.30   (hub + cross-community bridge)
#   Connector:       z < 2.5  AND 0.30 ≤ P < 0.62   (mid bridge, not hub)
#   Kinless:         z < 2.5  AND P ≥ 0.62   (no clear home community)
#   Peripheral:      z < 2.5  AND P < 0.30   (community-internal, not hub)
#
# Note: disease-name strings are NOT available in CKM dataset (ICD-10 3-char
#       codes only). Labels on scatter = ICD-3 codes. If WHO disease names
#       are needed later, they can be merged via a separate lookup.
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
})

# =========================
# 0) Paths & Parameters
# =========================
BRIDGE_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\community_detection\\medoid_bridge_metrics_coef0.01_alpha0.05.csv"
OUT_DIR    <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Guimerà-Amaral thresholds
Z_THRESH  <- 2.5    # hub threshold
P_THRESH1 <- 0.30   # peripheral | connector boundary
P_THRESH2 <- 0.62   # connector  | kinless    boundary

# Label selection
N_TOP_P   <- 15     # label top-N nodes by participation coefficient (on top of all hubs)
LABEL_MAXNCHAR <- 20  # truncate long ICD codes in labels (defensive)

FONT_FAMILY <- "Arial"
PDF_W       <- 11
PDF_H       <- 7.5

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Load data
# =========================
cat("==> Step 1) Loading bridge metrics...\n")

if (!file.exists(BRIDGE_CSV)) {
  stop("Bridge metrics CSV not found: ", BRIDGE_CSV, call. = FALSE)
}

bridge_dt <- fread(BRIDGE_CSV)
cat("   Loaded:", nrow(bridge_dt), "rows\n")
cat("   Columns:", paste(names(bridge_dt), collapse = ", "), "\n")

# Required columns
req_cols <- c("node", "community", "strength",
              "within_module_z_score", "participation_coef")
missing  <- setdiff(req_cols, names(bridge_dt))
if (length(missing) > 0) {
  stop("Required columns missing: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

# Ensure community is integer and keep only nodes in a community
bridge_dt[, community := as.integer(community)]
active <- bridge_dt[community > 0 & !is.na(participation_coef) &
                      !is.na(within_module_z_score)]
cat("   Active nodes (in a community, with z & P):", nrow(active), "\n")

# =========================
# 2) Classify roles (or validate existing)
# =========================
cat("\n==> Step 2) Classifying node roles (Guimerà-Amaral framework)...\n")

# If `role` column exists in CSV, validate; otherwise classify fresh
active[, role := fcase(
  within_module_z_score >= Z_THRESH & participation_coef >= P_THRESH1, "connector_hub",
  within_module_z_score >= Z_THRESH & participation_coef <  P_THRESH1, "provincial_hub",
  within_module_z_score <  Z_THRESH & participation_coef >= P_THRESH2, "kinless",
  within_module_z_score <  Z_THRESH & participation_coef >= P_THRESH1, "connector",
  within_module_z_score <  Z_THRESH & participation_coef <  P_THRESH1, "peripheral",
  default = "unclassified"
)]

role_labels <- c(
  "provincial_hub" = "Provincial hub",
  "connector_hub"  = "Connector hub",
  "kinless"        = "Kinless",
  "connector"      = "Connector",
  "peripheral"     = "Peripheral",
  "unclassified"   = "Unclassified"
)
role_colors <- c(
  "Provincial hub" = "#E41A1C",   # red
  "Connector hub"  = "#FF7F00",   # orange
  "Kinless"        = "#984EA3",   # purple
  "Connector"      = "#377EB8",   # blue
  "Peripheral"     = "#BDBDBD",   # grey
  "Unclassified"   = "#9E9E9E"
)

active[, role_label := role_labels[role]]
active[is.na(role_label), role_label := "Unclassified"]

# Count per role and drop empty levels for legend cleanliness
role_counts <- active[, .N, by = role_label]
present_roles <- role_counts$role_label
present_roles <- intersect(
  c("Provincial hub", "Connector hub", "Kinless", "Connector",
    "Peripheral", "Unclassified"),
  present_roles
)
active[, role_label := factor(as.character(role_label), levels = present_roles)]
role_colors_use <- role_colors[present_roles]
legend_labels <- setNames(
  paste0(present_roles, " (n=",
         role_counts[match(present_roles, role_label), N], ")"),
  present_roles
)

cat("   Role distribution:\n")
for (rl in present_roles) {
  cat(sprintf("     %-15s : %d\n", rl,
              role_counts[role_label == rl, N]))
}

# =========================
# 3) Build main scatter (Figure S4)
# =========================
cat("\n==> Step 3) Building z-P scatter plot...\n")

# Label nodes: all hubs + top-N by P
hubs_nodes <- active[within_module_z_score >= Z_THRESH, node]
topP_nodes <- active[order(-participation_coef)][seq_len(min(N_TOP_P, .N)), node]
label_nodes <- union(hubs_nodes, topP_nodes)

active[, show_label := fifelse(node %in% label_nodes, as.character(node), NA_character_)]
active[!is.na(show_label) & nchar(show_label) > LABEL_MAXNCHAR,
       show_label := paste0(substr(show_label, 1, LABEL_MAXNCHAR - 1), "\u2026")]

# Axis ranges
z_max <- max(active$within_module_z_score, na.rm = TRUE) * 1.08
z_min <- min(active$within_module_z_score, na.rm = TRUE) - 0.5

p_scatter <- ggplot(active,
                    aes(x = participation_coef, y = within_module_z_score)) +
  # Threshold lines
  geom_hline(yintercept = Z_THRESH,  linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = P_THRESH1, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = P_THRESH2, linetype = "dotted",
             color = "grey70", linewidth = 0.4) +
  # Quadrant labels
  annotate("text", x = P_THRESH1 / 2,            y = z_max * 0.95,
           label = "Provincial\nhub", color = "grey45",
           size = 3.5, fontface = "italic", family = FONT_FAMILY) +
  annotate("text", x = (P_THRESH1 + P_THRESH2) / 2, y = z_max * 0.95,
           label = "Connector\nhub", color = "grey45",
           size = 3.5, fontface = "italic", family = FONT_FAMILY) +
  annotate("text", x = P_THRESH1 / 2,            y = z_min + 0.3,
           label = "Peripheral",  color = "grey45",
           size = 3.5, fontface = "italic", family = FONT_FAMILY) +
  annotate("text", x = (P_THRESH1 + P_THRESH2) / 2, y = z_min + 0.3,
           label = "Connector",  color = "grey45",
           size = 3.5, fontface = "italic", family = FONT_FAMILY) +
  annotate("text", x = (P_THRESH2 + 1) / 2,      y = z_min + 0.3,
           label = "Kinless",    color = "grey45",
           size = 3.5, fontface = "italic", family = FONT_FAMILY) +
  # Peripheral: jitter slightly on x to reduce P=0 overplotting
  geom_point(data = active[role_label == "Peripheral"],
             aes(color = role_label, size = strength),
             alpha = 0.35,
             position = position_jitter(width = 0.006, height = 0, seed = 42)) +
  # Other roles: no jitter
  geom_point(data = active[role_label != "Peripheral"],
             aes(color = role_label, size = strength),
             alpha = 0.75) +
  # Labels with repel
  geom_text_repel(data = active[!is.na(show_label)],
                  aes(label = show_label),
                  size = 3, family = FONT_FAMILY, color = "grey20",
                  segment.color = "grey60", segment.size = 0.3,
                  box.padding = 0.4, point.padding = 0.3,
                  min.segment.length = 0.2,
                  max.overlaps = 25, seed = 42) +
  scale_color_manual(name   = "Node role",
                     values = role_colors_use,
                     labels = legend_labels,
                     drop   = TRUE) +
  scale_size_continuous(name = "Strength",
                        range = c(1.5, 7), breaks = pretty) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(z_min, z_max)) +
  labs(
    x = "Participation coefficient (P)",
    y = "Within-module z-score (z)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 11),
    axis.line        = element_blank(),
    axis.ticks       = element_line(color = "black", linewidth = 0.6),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    legend.position  = "right",
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4)),
         size  = guide_legend(order = 2))

pdf_file <- file.path(OUT_DIR, "FigS4_node_role_scatter_zP.pdf")
ggsave(pdf_file, p_scatter, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved:", basename(pdf_file), "\n")

# =========================
# 4) Table S2: Top 20 bridge diseases (by participation coefficient)
# =========================
cat("\n==> Step 4) Writing Table S2 (Top 20 bridge diseases)...\n")

top20 <- head(active[order(-participation_coef)], 20)

# Helper for "top-N target community (share%)" cells
fmt_target <- function(comm_id, share_total) {
  ifelse(is.na(comm_id), "\u2014",
         sprintf("C%s (%.1f%%)", comm_id, share_total * 100))
}

# Columns available in CSV vary; pull defensively
has_top1 <- all(c("top1_ext_community", "top1_share_total") %in% names(active))
has_top2 <- all(c("top2_ext_community", "top2_share_total") %in% names(active))
has_top3 <- all(c("top3_ext_community", "top3_share_total") %in% names(active))
has_ext  <- "external_strength" %in% names(active)
has_chap <- "ICD10_Chapter" %in% names(active)

tab19 <- data.table(
  Rank              = 1:nrow(top20),
  `ICD-10 code`     = top20$node,
  `ICD-10 chapter`  = if (has_chap) top20$ICD10_Chapter else NA_character_,
  Community         = top20$community,
  Role              = role_labels[top20$role],
  P                 = round(top20$participation_coef, 3),
  z                 = round(top20$within_module_z_score, 2),
  Strength          = round(top20$strength, 3),
  `External strength` = if (has_ext) round(top20$external_strength, 3) else NA_real_
)
if (has_top1) tab19[, `Top 1 target` := fmt_target(top20$top1_ext_community,
                                                   top20$top1_share_total)]
if (has_top2) tab19[, `Top 2 target` := fmt_target(top20$top2_ext_community,
                                                   top20$top2_share_total)]
if (has_top3) tab19[, `Top 3 target` := fmt_target(top20$top3_ext_community,
                                                   top20$top3_share_total)]

# Write Excel
wb20 <- createWorkbook()
hdr_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF",
                         halign = "center", valign = "center",
                         fgFill = "#4472C4", textDecoration = "bold")

addWorksheet(wb20, "Top20_bridges")
writeDataTable(wb20, "Top20_bridges", tab19, withFilter = TRUE)
setColWidths(wb20, "Top20_bridges", cols = 1:ncol(tab19), widths = "auto")
addStyle(wb20, "Top20_bridges", hdr_style, rows = 1,
         cols = 1:ncol(tab19), gridExpand = TRUE)

# Row colour by role
role_fill <- c(
  "Provincial hub" = "#FDDEDE",
  "Connector hub"  = "#FFE8CC",
  "Kinless"        = "#E8D5F5",
  "Connector"      = "#D6EAF8",
  "Peripheral"     = "#F2F2F2"
)
for (r in seq_len(nrow(tab19))) {
  fill_color <- role_fill[tab19[r, Role]]
  if (!is.na(fill_color)) {
    addStyle(wb20, "Top20_bridges", createStyle(fgFill = fill_color),
             rows = r + 1, cols = 1:ncol(tab19),
             gridExpand = TRUE, stack = TRUE)
  }
}

xlsx20 <- file.path(OUT_DIR, "TableS2_top20_bridge_diseases.xlsx")
saveWorkbook(wb20, xlsx20, overwrite = TRUE)
cat("   Saved:", basename(xlsx20), "\n")

# =========================
# 5) Table S2b: Full list of P >= 0.30 (Connector + higher)
# =========================
cat("\n==> Step 5) Writing full bridge node list (P >= 0.30)...\n")

bridge_all <- active[participation_coef >= P_THRESH1][order(-participation_coef)]

tabS2b <- data.table(
  `ICD-10 code`           = bridge_all$node,
  `ICD-10 chapter`        = if (has_chap) bridge_all$ICD10_Chapter else NA_character_,
  Community               = bridge_all$community,
  Role                    = role_labels[bridge_all$role],
  P                       = round(bridge_all$participation_coef, 4),
  z                       = round(bridge_all$within_module_z_score, 3),
  Strength                = round(bridge_all$strength, 3),
  `Within-module strength`= if ("strength_within_module" %in% names(bridge_all))
    round(bridge_all$strength_within_module, 3) else NA_real_,
  `External strength`     = if (has_ext) round(bridge_all$external_strength, 3) else NA_real_
)
if (has_top1) {
  tabS2b[, `Top 1 target community` := fifelse(is.na(bridge_all$top1_ext_community),
                                               "\u2014",
                                               as.character(bridge_all$top1_ext_community))]
  tabS2b[, `Top 1 share (%)`        := fifelse(is.na(bridge_all$top1_share_total),
                                               "\u2014",
                                               sprintf("%.1f", bridge_all$top1_share_total * 100))]
}
if (has_top2) {
  tabS2b[, `Top 2 target community` := fifelse(is.na(bridge_all$top2_ext_community),
                                               "\u2014",
                                               as.character(bridge_all$top2_ext_community))]
  tabS2b[, `Top 2 share (%)`        := fifelse(is.na(bridge_all$top2_share_total),
                                               "\u2014",
                                               sprintf("%.1f", bridge_all$top2_share_total * 100))]
}
if (has_top3) {
  tabS2b[, `Top 3 target community` := fifelse(is.na(bridge_all$top3_ext_community),
                                               "\u2014",
                                               as.character(bridge_all$top3_ext_community))]
  tabS2b[, `Top 3 share (%)`        := fifelse(is.na(bridge_all$top3_share_total),
                                               "\u2014",
                                               sprintf("%.1f", bridge_all$top3_share_total * 100))]
}

cat("   Total bridge nodes (P >= 0.30):", nrow(tabS2b), "\n")

wbS2b <- createWorkbook()
addWorksheet(wbS2b, "Bridge_nodes")
writeDataTable(wbS2b, "Bridge_nodes", tabS2b, withFilter = TRUE)
setColWidths(wbS2b, "Bridge_nodes", cols = 1:ncol(tabS2b), widths = "auto")
addStyle(wbS2b, "Bridge_nodes", hdr_style, rows = 1,
         cols = 1:ncol(tabS2b), gridExpand = TRUE)

# Role colouring
for (r in seq_len(nrow(tabS2b))) {
  fill_color <- role_fill[tabS2b[r, Role]]
  if (!is.na(fill_color)) {
    addStyle(wbS2b, "Bridge_nodes", createStyle(fgFill = fill_color),
             rows = r + 1, cols = 1:ncol(tabS2b),
             gridExpand = TRUE, stack = TRUE)
  }
}

# Role summary sheet
role_summary <- bridge_all[, .(
  N              = .N,
  `Mean P`       = round(mean(participation_coef, na.rm = TRUE), 3),
  `Mean z`       = round(mean(within_module_z_score, na.rm = TRUE), 2),
  `Mean strength`= round(mean(strength, na.rm = TRUE), 2)
), by = .(Role = role_labels[role])]
setorder(role_summary, -`Mean P`)

addWorksheet(wbS2b, "Role_summary")
writeDataTable(wbS2b, "Role_summary", role_summary, withFilter = TRUE)
setColWidths(wbS2b, "Role_summary", cols = 1:ncol(role_summary), widths = "auto")
addStyle(wbS2b, "Role_summary", hdr_style, rows = 1,
         cols = 1:ncol(role_summary), gridExpand = TRUE)

xlsxS2b <- file.path(OUT_DIR, "TableS2b_all_bridge_nodes.xlsx")
saveWorkbook(wbS2b, xlsxS2b, overwrite = TRUE)
cat("   Saved:", basename(xlsxS2b), "\n")

# =========================
# 6) Console summary
# =========================
cat("\n============================================================\n")
cat("Step 10: Node role classification — complete.\n")
cat("============================================================\n")
cat("\nTop 5 bridge diseases (highest participation coefficient):\n")
for (r in seq_len(min(5, nrow(top20)))) {
  row <- top20[r]
  cat(sprintf("   %d. %s  (Community C%d, Ch.%s, role=%s)  P=%.3f  z=%.2f  Strength=%.2f\n",
              r, row$node, row$community,
              if (has_chap) as.character(row$ICD10_Chapter) else "NA",
              row$role,
              row$participation_coef, row$within_module_z_score,
              row$strength))
}

cat("\nFiles written to:\n  ", OUT_DIR, "\n", sep = "")
cat("  1) FigS4_node_role_scatter_zP.pdf\n")
cat("  2) TableS2_top20_bridge_diseases.xlsx\n")
cat("  3) TableS2b_all_bridge_nodes.xlsx  (P >= 0.30)\n")
cat("============================================================\n")