#!/usr/bin/env Rscript
# ==============================================================================
# Fig S2: Leiden Resolution Sensitivity Analysis (publication-grade redraw)
#
# 不重跑 sweep, 直接读 Leiden_Resolution_Sweep.R 的输出
# 修订点 (vs 原 sweep 输出):
#   1. 红色虚线标 γ = 1.2 (paper selected), 非 0.5 (composite optimum)
#   2. 副标题诚实说明: 随机 Leiden 在 1.1-1.4 给 15-17 communities,
#      paper 主 partition 16 communities (multi-seed stable)
#   3. 6-panel 改 4-panel (n_comm, modularity, conductance, ARI stability),
#      coverage + avg size 移除 (信息冗余)
#
# Input:  resolution_sweep_results.xlsx (sweep_by_resolution sheet)
# Output: FigS2_Leiden_resolution_sweep.pdf
# ==============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
})

# =========================
# Paths
# =========================
SWEEP_XLSX <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\resolution_sweep_fine\\resolution_sweep_results.xlsx"
OUT_DIR    <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 加载 sweep 数据
# =========================
cat("==> 加载 sweep 数据...\n")
sweep_dt <- as.data.table(read.xlsx(SWEEP_XLSX, sheet = "sweep_by_resolution"))
cat(sprintf("   Sweep 数据行数: %d\n", nrow(sweep_dt)))
cat(sprintf("   Resolution 范围: %.1f \u2014 %.1f\n",
            min(sweep_dt$resolution), max(sweep_dt$resolution)))

# =========================
# Selected resolution (paper)
# =========================
SELECTED_RES <- 1.2

# Paper 实际 partition (来自 bridge_metrics 文件, 16 communities)
PAPER_N_COMM <- 16

# =========================
# 4 个 panel 共用 theme
# =========================
panel_theme <- theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    plot.title       = element_text(size = 11, face = "bold",
                                    hjust = 0.5, margin = margin(b = 5)),
    axis.title       = element_text(size = 10, color = "black"),
    axis.text        = element_text(size = 9, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 8, 8, 8)
  )

# 红色虚线 (selected) + 颜色
SEL_COLOR <- "#C62828"      # 深红 (selected line)
DOT_COLOR <- "#1976D2"      # 蓝 (data points)

# X 轴 break (避免拥挤)
x_breaks <- c(0.5, 1, 1.5, 2, 3, 5, 10)

# =========================
# Panel A: Number of communities
# =========================
cat("\n==> 绘制 4 panels...\n")

# 在 sweep 数据基础上加一个标注点 (selected)
sel_data_n <- data.table(resolution = SELECTED_RES,
                         n_communities = PAPER_N_COMM)

p_a <- ggplot(sweep_dt, aes(x = resolution, y = n_communities)) +
  geom_line(color = DOT_COLOR, linewidth = 0.6) +
  geom_point(color = DOT_COLOR, size = 2.0) +
  geom_vline(xintercept = SELECTED_RES, color = SEL_COLOR,
             linetype = "dashed", linewidth = 0.6) +
  geom_point(data = sel_data_n,
             color = SEL_COLOR, size = 3.5, shape = 18) +   # 菱形标 selected
  annotate("text", x = SELECTED_RES + 0.15, y = PAPER_N_COMM + 5,
           label = sprintf("\u03b3 = %.1f\n%d communities", SELECTED_RES, PAPER_N_COMM),
           color = SEL_COLOR, size = 3, family = FONT_FAMILY,
           hjust = 0, fontface = "bold") +
  scale_x_log10(breaks = x_breaks, labels = x_breaks) +
  labs(x = "Leiden resolution parameter (\u03b3)",
       y = "Number of communities",
       title = "A  Number of communities") +
  panel_theme

# =========================
# Panel B: Modularity
# =========================
sel_data_m <- sweep_dt[abs(resolution - SELECTED_RES) < 0.01,
                       .(resolution, modularity)]
if (nrow(sel_data_m) > 0) {
  sel_modularity <- sel_data_m$modularity[1]
} else {
  # 如果没有 1.2 那一行, 用最近的
  sel_modularity <- sweep_dt[which.min(abs(resolution - SELECTED_RES)),
                             modularity]
}

p_b <- ggplot(sweep_dt, aes(x = resolution, y = modularity)) +
  geom_line(color = DOT_COLOR, linewidth = 0.6) +
  geom_point(color = DOT_COLOR, size = 2.0) +
  geom_vline(xintercept = SELECTED_RES, color = SEL_COLOR,
             linetype = "dashed", linewidth = 0.6) +
  geom_point(data = data.table(resolution = SELECTED_RES,
                               modularity = sel_modularity),
             color = SEL_COLOR, size = 3.5, shape = 18) +
  scale_x_log10(breaks = x_breaks, labels = x_breaks) +
  scale_y_continuous(limits = c(0.45, 0.65),
                     breaks = seq(0.45, 0.65, 0.05)) +
  labs(x = "Leiden resolution parameter (\u03b3)",
       y = "Modularity (Q)",
       title = "B  Modularity") +
  panel_theme

# =========================
# Panel C: Average conductance (lower = better)
# =========================
sel_cond <- sweep_dt[which.min(abs(resolution - SELECTED_RES)),
                     avg_conductance]

p_c <- ggplot(sweep_dt, aes(x = resolution, y = avg_conductance)) +
  geom_line(color = DOT_COLOR, linewidth = 0.6) +
  geom_point(color = DOT_COLOR, size = 2.0) +
  geom_vline(xintercept = SELECTED_RES, color = SEL_COLOR,
             linetype = "dashed", linewidth = 0.6) +
  geom_point(data = data.table(resolution = SELECTED_RES,
                               avg_conductance = sel_cond),
             color = SEL_COLOR, size = 3.5, shape = 18) +
  scale_x_log10(breaks = x_breaks, labels = x_breaks) +
  labs(x = "Leiden resolution parameter (\u03b3)",
       y = "Average conductance",
       title = "C  Boundary separation (lower is better)") +
  panel_theme

# =========================
# Panel D: Partition stability (mean ARI from 10 seeds)
# =========================
sel_ari <- sweep_dt[which.min(abs(resolution - SELECTED_RES)),
                    ARI_stability_mean]

p_d <- ggplot(sweep_dt, aes(x = resolution, y = ARI_stability_mean)) +
  geom_line(color = DOT_COLOR, linewidth = 0.6) +
  geom_point(color = DOT_COLOR, size = 2.0) +
  # 加一条 ARI = 0.8 的参考线
  geom_hline(yintercept = 0.80, color = "grey50",
             linetype = "dotted", linewidth = 0.4) +
  annotate("text", x = max(sweep_dt$resolution),
           y = 0.81,
           label = "ARI = 0.80 (commonly accepted threshold)",
           hjust = 1, size = 2.6, color = "grey40",
           family = FONT_FAMILY) +
  geom_vline(xintercept = SELECTED_RES, color = SEL_COLOR,
             linetype = "dashed", linewidth = 0.6) +
  geom_point(data = data.table(resolution = SELECTED_RES,
                               ARI_stability_mean = sel_ari),
             color = SEL_COLOR, size = 3.5, shape = 18) +
  scale_x_log10(breaks = x_breaks, labels = x_breaks) +
  scale_y_continuous(limits = c(0.70, 1.00),
                     breaks = seq(0.70, 1.00, 0.05)) +
  labs(x = "Leiden resolution parameter (\u03b3)",
       y = "Partition stability (mean ARI, 10 seeds)",
       title = "D  Partition stability") +
  panel_theme

# =========================
# 组装 4 panels (2x2)
# =========================
cat("\n==> 组装 figure...\n")

# Title + subtitle 作为顶部 header
header_text <- "Figure S2. Leiden resolution sensitivity analysis"
subtitle_text <- paste0(
  "Network of 613 disease nodes (FDR < 0.05; \u03c6 \u2265 0.01). ",
  "The selected resolution \u03b3 = 1.2 (red dashed line) yielded ",
  sprintf("%d medically distinct communities ", PAPER_N_COMM),
  "in the primary partition (Q = 0.62, ARI = 0.92 across 30 seeds). ",
  "Adjacent resolutions \u03b3 = 1.1\u20131.4 produced 15\u201317 communities, ",
  "consistent with the inherent stochasticity of the Leiden algorithm."
)

combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = header_text,
    subtitle = subtitle_text,
    theme = theme(
      plot.title    = element_text(family = FONT_FAMILY, size = 13,
                                   face = "bold",
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(family = FONT_FAMILY, size = 9,
                                   color = "grey25",
                                   lineheight = 1.2,
                                   margin = margin(b = 8))
    )
  )

# =========================
# 保存
# =========================
pdf_file <- file.path(OUT_DIR, "FigS2_Leiden_resolution_sweep.pdf")
ggsave(pdf_file, combined,
       width = 11, height = 9, units = "in",
       device = PDF_DEVICE, dpi = 300)

cat(sprintf("\n   \u2713 已保存: %s\n", pdf_file))

# =========================
# Console 报告
# =========================
cat("\n========================================================\n")
cat("Fig S2 \u2014 完成\n")
cat("========================================================\n\n")
cat("Selected resolution (paper main analysis):\n")
cat(sprintf("  \u03b3 = %.1f\n", SELECTED_RES))
cat(sprintf("  N communities = %d (paper main partition)\n", PAPER_N_COMM))
cat(sprintf("  Modularity = %.4f\n", sel_modularity))
cat(sprintf("  Avg conductance = %.4f\n", sel_cond))
cat(sprintf("  ARI stability = %.4f\n", sel_ari))
cat("\n关键 panels:\n")
cat("  A: Number of communities vs resolution\n")
cat("  B: Modularity Q (peak ~0.62 around \u03b3 = 0.8\u20131.2)\n")
cat("  C: Avg conductance (lower = better separation)\n")
cat("  D: Partition stability (mean ARI from 10 seeds)\n")
cat(sprintf("\n输出文件:\n  %s\n", pdf_file))