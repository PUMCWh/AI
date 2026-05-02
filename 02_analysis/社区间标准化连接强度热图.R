# ==============================================================================
# 社区间标准化连接强度热图
#
# 输入：
#   - consensus_results_all.rds（来自社区发现脚本输出，含medoid分区及社区间矩阵）
#
# 输出：
#   - 社区间连接强度热图.pdf
#
# 方法：
#   原始社区间总权重 = Σφ(社区i与社区j之间所有边的偏相关系数之和)
#   该值受社区规模影响（大社区节点对数多，天然总权重更高），
#   因此以节点对数标准化：
#     标准化连接强度 = Σφ / (n_i × n_j)
#   其中 n_i, n_j 为两个社区的节点数。
#   该指标等价于"加权社区间边密度"，可公平比较不同规模社区间的关联紧密程度。
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# 0) Paths & Params
# =========================
COMM_DIR <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\社区检测"
OUT_DIR  <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\社区检测"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MAIN_TAG <- "coef0.01_alpha0.05"

# 图参数
FONT_FAMILY <- "Arial"
PDF_W <- 10
PDF_H <- 9

# PDF device
PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Load data
# =========================
cat("==> Step 1) Loading data...\n")

consensus_results <- readRDS(file.path(COMM_DIR, "consensus_results_all.rds"))
res <- consensus_results[[MAIN_TAG]]

# medoid分区（社区归属）
medoid_dt <- res$medoid_dt
medoid_dt <- medoid_dt[community > 0]

# 原始社区间总权重矩阵（Σφ）
raw_matrix <- res$medoid_comm_matrix
if (is.null(raw_matrix)) {
  stop("medoid_comm_matrix not found in consensus_results. ",
       "Please re-run 总体-社区发现_update.R first.")
}

cat("   Communities:", length(unique(medoid_dt$community)), "\n")
cat("   Matrix dim:", nrow(raw_matrix), "×", ncol(raw_matrix), "\n")

# =========================
# 2) Compute normalized connectivity
# =========================
cat("==> Step 2) Computing normalized inter-community connectivity...\n")

# 社区规模
comm_size <- medoid_dt[, .N, by = community]
setkey(comm_size, community)
size_vec <- setNames(comm_size$N, as.character(comm_size$community))

# 确保矩阵行列名与社区对齐
if (is.null(rownames(raw_matrix)) || is.null(colnames(raw_matrix))) {
  stop("raw_matrix must have rownames/colnames as community IDs.")
}
if (!all(rownames(raw_matrix) == colnames(raw_matrix))) {
  stop("raw_matrix rownames/colnames are not identical.")
}

# 提取数字社区ID（兼容 'C1' / '1'）
comm_ids <- as.integer(gsub("\\D+", "", rownames(raw_matrix)))
if (anyNA(comm_ids)) stop("Cannot parse community IDs from raw_matrix dimnames.")
ord <- order(comm_ids)
comm_ids <- comm_ids[ord]
comm_labels <- as.character(comm_ids)

# 重排矩阵（按社区编号升序）并统一dimnames为纯数字
raw_matrix <- raw_matrix[ord, ord, drop = FALSE]
rownames(raw_matrix) <- comm_labels
colnames(raw_matrix) <- comm_labels

# 标准化矩阵：norm_ij = raw_ij / (n_i × n_j)（矩阵运算，简洁不易错）
ni <- size_vec[comm_labels]
denom_mat <- outer(ni, ni, "*")
norm_matrix <- raw_matrix / denom_mat

# 对角线设为NA（社区内不展示，避免读者误解）
diag(norm_matrix) <- NA_real_

# 转长格式
dt_long <- as.data.table(as.table(norm_matrix[comm_labels, comm_labels]))
setnames(dt_long, c("from", "to", "norm_weight"))
dt_long[, norm_weight := as.numeric(norm_weight)]

# 同时保留原始总权重（用于格内标注）
dt_raw <- as.data.table(as.table(raw_matrix[comm_labels, comm_labels]))
setnames(dt_raw, c("from", "to", "raw_weight"))
dt_raw[, raw_weight := as.numeric(raw_weight)]
dt_long <- merge(dt_long, dt_raw, by = c("from", "to"))

# 仅保留社区间（去除对角线）用于标注和分析
dt_inter <- dt_long[from != to]

cat("   Normalized matrix computed.\n")
cat("   Off-diagonal range:",
    sprintf("%.4f – %.4f", min(dt_inter$norm_weight), max(dt_inter$norm_weight)), "\n")

# =========================
# 3) Determine scaling and build plot data
# =========================
# 格内文字：标准化值通常很小（0.001量级），需乘以缩放因子使数值易读

# 先看分布决定缩放因子
median_nw <- median(dt_inter[norm_weight > 0, norm_weight])
cat("   Median nonzero normalized weight:", sprintf("%.5f", median_nw), "\n")

# 自动选择缩放：使大多数数值在0.1-10之间
if (median_nw < 0.001) {
  scale_factor <- 10000
  scale_suffix <- expression("NCI" ~ "(×" * 10^4 * ")")
} else if (median_nw < 0.01) {
  scale_factor <- 1000
  scale_suffix <- expression("NCI" ~ "(×" * 10^3 * ")")
} else {
  scale_factor <- 100
  scale_suffix <- expression("NCI" ~ "(×" * 10^2 * ")")
}
# NCI = Normalized Connectivity Index，完整公式见图注

cat("   Scale factor:", scale_factor, "\n")

# 绘图用完整数据（含对角线，对角线 norm_weight=NA → 灰色）
dt_plot <- copy(dt_long)
dt_plot[, scaled_weight := norm_weight * scale_factor]       # 对角线NA → NA * x = NA
dt_plot[, cell_label := fifelse(!is.na(scaled_weight) & scaled_weight > 0,
                                sprintf("%.1f", scaled_weight), "")]
dt_plot[, from := factor(from, levels = comm_labels)]
dt_plot[, to   := factor(to,   levels = comm_labels)]

# =========================
# 4) Create heatmap
# =========================
cat("==> Step 4) Creating heatmap...\n")

# 色阶上限：实际最大值向上取整到最近的5的倍数（美观且覆盖全部数据）
actual_max <- max(dt_plot[!is.na(scaled_weight), scaled_weight])
fill_max <- ceiling(actual_max / 5) * 5  # 向上取整到5的倍数（如17.3→20, 14.2→15）
if (fill_max == 0) fill_max <- 5

# 轴标签映射（命名向量，避免错位）
lab_map <- setNames(paste0("C", comm_labels), comm_labels)

p_heatmap <- ggplot(dt_plot, aes(x = from, y = to)) +
  # 填充
  geom_tile(aes(fill = scaled_weight), color = "grey85", linewidth = 0.3) +
  # 格内数字
  geom_text(aes(label = cell_label),
            size = 3, color = "grey10", family = FONT_FAMILY) +
  # 色阶（顺序色阶：白→红）
  scale_fill_gradient(
    low = "white", high = "#B2182B",
    limits = c(0, fill_max),
    breaks = {
      # 选择步长使得刻度为整数且总数3-6个
      candidates <- c(1, 2, 5, 10, 20)
      step <- candidates[which(fill_max / candidates >= 3 & fill_max / candidates <= 6)[1]]
      if (is.na(step)) step <- fill_max / 4
      seq(0, fill_max, by = step)
    },
    oob = scales::squish,
    name = scale_suffix,
    na.value = "grey95"  # 对角线NA显示为浅灰
  ) +
  # 轴标签
  scale_x_discrete(labels = lab_map) +
  scale_y_discrete(labels = lab_map) +
  labs(x = "Community", y = "Community") +
  # 主题
  theme_minimal(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    axis.text.x      = element_text(size = 12, angle = 0, hjust = 0.5, color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.title       = element_text(size = 14, color = "black"),
    legend.position  = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.title     = element_text(size = 14, color = "black"),
    legend.text      = element_text(size = 12, color = "black"),
    panel.grid       = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

# 保存PDF
pdf_file <- file.path(OUT_DIR, "社区间连接强度热图.pdf")
ggsave(pdf_file, p_heatmap, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved:", pdf_file, "\n")

# =========================
# 5) Console summary: top 10 most connected community pairs
# =========================
cat("\n==> Top 10 most strongly connected community pairs (normalized):\n")
dt_pairs <- copy(dt_inter)
dt_pairs[, from_i := as.integer(as.character(from))]
dt_pairs[, to_i   := as.integer(as.character(to))]
dt_pairs <- dt_pairs[from_i < to_i]
dt_pairs[, c("from_i", "to_i") := NULL]
setorder(dt_pairs, -norm_weight)
top10 <- head(dt_pairs, 10)
for (r in seq_len(nrow(top10))) {
  row <- top10[r]
  cat(sprintf("   C%s – C%s : normalized = %.5f (×%d = %.2f), raw Σφ = %.2f\n",
              row$from, row$to, row$norm_weight, scale_factor,
              row$norm_weight * scale_factor, row$raw_weight))
}

# =========================
# 6) 同时输出原始总权重版本（用于对比/备查）
# =========================
cat("\n==> Bonus: saving raw total weight version for comparison...\n")

dt_plot_raw <- copy(dt_long)
dt_plot_raw[, from := factor(from, levels = comm_labels)]
dt_plot_raw[, to   := factor(to,   levels = comm_labels)]
dt_plot_raw[from == to, raw_weight := NA_real_]  # 对角线NA
dt_plot_raw[, cell_label := fifelse(!is.na(raw_weight) & raw_weight > 0,
                                    sprintf("%.2f", raw_weight), "")]

p_raw <- ggplot(dt_plot_raw, aes(x = from, y = to)) +
  geom_tile(aes(fill = raw_weight), color = "grey85", linewidth = 0.3) +
  geom_text(aes(label = cell_label), size = 2.3, color = "grey10",
            family = FONT_FAMILY) +
  scale_fill_gradient(low = "white", high = "#B70031",
                      name = expression(paste("Total weight (", Sigma, phi, ")")),
                      na.value = "grey95") +
  scale_x_discrete(labels = lab_map) +
  scale_y_discrete(labels = lab_map) +
  labs(x = "Community", y = "Community") +
  theme_minimal(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 8),
    legend.position  = "right",
    legend.key.height = unit(1.5, "cm"),
    panel.grid       = element_blank()
  ) +
  coord_fixed()

ggsave(file.path(OUT_DIR, "社区间连接强度热图_原始总权重.pdf"),
       p_raw, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved: 社区间连接强度热图_原始总权重.pdf (for comparison)\n")

cat("\n========================================\n")
cat("DONE.\n")
cat("  主图 = 标准化版（校正社区规模）\n")
cat("  备查 = 原始总权重版\n")
cat("========================================\n")