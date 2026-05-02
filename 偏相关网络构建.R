# ==============================================================================
# CKM 共病研究 — 步骤2: 偏相关网络构建 (crossprod加速版)
# ==============================================================================
# 算法核心 (来自 共病偏相关网络分析_总体_更新MMC的计算_update.R):
#   1) 从Step 1的稀疏矩阵直接用 crossprod(X) 计算P×P相关矩阵
#   2) 求逆(带ridge正则化)得到偏相关矩阵
#   3) 全程只操作P×P矩阵(~926×926), 不创建N×P密集矩阵(~128万×926)
#   → 内存: ~7MB vs ~9.5GB; 速度: 分钟级 vs 可能卡死
#
# 输入: Step 1 输出的 patient_disease_matrix_strict.rds
# 输出: RDS中间文件 (供Step 3社区发现脚本直接读取) + 网络图 + MMC + Excel
# ==============================================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  pkgs <- c("data.table", "Matrix", "igraph", "ggplot2", "ggraph", "grid", "openxlsx")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  library(data.table); library(Matrix); library(igraph)
  library(ggplot2); library(ggraph); library(grid); library(openxlsx)
})

options(stringsAsFactors = FALSE)
t0 <- Sys.time()

cat("==============================================================================\n")
cat("  CKM 共病研究 — 步骤2: 偏相关网络构建 (crossprod加速版)\n")
cat("==============================================================================\n\n")

# =========================
# 0) 路径与参数
# =========================
STEP1_DIR    <- "C:/Users/HP/Desktop/paper/CKM_FINAL"
OUT_BASE_DIR <- "C:/Users/HP/Desktop/paper/CKM_FINAL"

# 患病率阈值
PREV_THRESHOLD <- 0.0002  # 0.02%

# FDR校正方法
P_ADJ_METHOD <- "BH"

# Ridge正则化 (求逆时数值稳定性)
RIDGE_INIT <- 1e-6
RIDGE_MAX  <- 1e-1

# 阈值组合 (用于多版本网络/MMC/图)
THRESH_GRID <- data.table(
  coef  = c(0.10, 0.05, 0.01, 0.10, 0.05, 0.01),
  alpha = c(0.05, 0.05, 0.05, 0.01, 0.01, 0.01)
)

# 可视化参数 (与已验证的参考脚本完全一致)
PDF_W <- 11.5; PDF_H <- 10
FONT_FAMILY <- "Arial"
SHOW_NODE_LABELS <- TRUE
NODE_LABEL_SIZE  <- 2.2

# 是否保存全量疾病对CSV (P choose 2, 可能很大; 默认关闭)
SAVE_ALL_PAIRS_CSV <- FALSE

# ICD章节颜色
chapter_colors_provided <- c(
  "1" = "#B70031", "2" = "#EC6A3D", "3" = "#A1C6AD", "4" = "#0098C5",
  "5" = "#9A6CAC", "6" = "#004B9B", "7" = "#C0C7CB", "8" = "#AAD9EA",
  "9" = "#FFE100", "10" = "#E71A10", "11" = "#6EB327", "12" = "#8E006C",
  "13" = "#D2ACC9", "14" = "#ADB5DA", "0" = "grey"
)

# ICD章节映射
get_chapter <- function(code) {
  sapply(code, function(c) {
    if      (c >= "A00" & c <= "B99") "1"
    else if (c >= "C00" & c <= "D48") "2"
    else if (c >= "D50" & c <= "D89") "3"
    else if (c >= "E00" & c <= "E90") "4"
    else if (c >= "F00" & c <= "F99") "5"
    else if (c >= "G00" & c <= "G99") "6"
    else if (c >= "H00" & c <= "H59") "7"
    else if (c >= "H60" & c <= "H95") "8"
    else if (c >= "I00" & c <= "I99") "9"
    else if (c >= "J00" & c <= "J99") "10"
    else if (c >= "K00" & c <= "K93") "11"
    else if (c >= "L00" & c <= "L99") "12"
    else if (c >= "M00" & c <= "M99") "13"
    else if (c >= "N00" & c <= "N99") "14"
    else "0"
  })
}

# 输出目录 (带时间戳, 兼容社区发现脚本的自动检测)
ts_tag  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(OUT_BASE_DIR, paste0("NETWORK_PCOR_", ts_tag))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) 读取Step 1输出
# =========================
cat("==> Step 1) 读取Step 1矩阵...\n")

mat_data   <- readRDS(file.path(STEP1_DIR, "patient_disease_matrix_strict.rds"))
X_full     <- mat_data$matrix
dmap_full  <- mat_data$disease_map   # icd_idx, ICD_3
pmap       <- mat_data$patient_map   # pid_idx, 人员ID

n_patients <- nrow(X_full)
n_diseases_raw <- ncol(X_full)
cat("   患者数:", format(n_patients, big.mark = ","), "\n")
cat("   原始疾病数:", n_diseases_raw, "\n")

# 计算频率和患病率
col_freq <- as.integer(Matrix::colSums(X_full))
col_prev <- col_freq / n_patients

# =========================
# 2) 患病率过滤
# =========================
cat("\n==> Step 2) 患病率过滤 (>", PREV_THRESHOLD, ")...\n")

keep_idx <- which(col_prev > PREV_THRESHOLD)
k <- length(keep_idx)
cat("   过滤后疾病数:", k, "\n")
if (k < 2) stop("过滤后疾病数不足2个。")

X <- X_full[, keep_idx]
code_by_id    <- dmap_full$ICD_3[keep_idx]
chapter_by_id <- as.character(get_chapter(code_by_id))
prev_by_id    <- col_prev[keep_idx]
freq_by_id    <- col_freq[keep_idx]

# 构建disease_map (社区发现脚本所需格式)
disease_map <- data.table(
  disease_id    = seq_len(k),
  code          = code_by_id,
  ICD10_Chapter = chapter_by_id,
  freq          = freq_by_id,
  prev          = prev_by_id
)

rm(X_full, dmap_full, col_freq, col_prev); gc()
cat("   稀疏矩阵:", nrow(X), "×", ncol(X),
    "  非零元素:", format(nnzero(X), big.mark = ","), "\n")

# =========================
# 3) 计算Pearson相关矩阵 R (P×P, 通过crossprod)
# =========================
cat("\n==> Step 3) 计算相关矩阵 R (crossprod, 仅P×P)...\n")
t_step <- Sys.time()

n <- n_patients
mu   <- as.numeric(Matrix::colSums(X)) / n
var_s <- (n / (n - 1)) * mu * (1 - mu)
sd_s  <- sqrt(pmax(var_s, 1e-12))

XtX <- as.matrix(Matrix::crossprod(X))
rm(X); gc()

cov_mat <- (XtX - n * tcrossprod(mu, mu)) / (n - 1)
rm(XtX); gc()

R <- cov_mat / tcrossprod(sd_s, sd_s)
rm(cov_mat); gc()

diag(R) <- 1
R[R >  1] <-  1
R[R < -1] <- -1

cat("   相关矩阵维度:", nrow(R), "×", ncol(R), "\n")
cat("   耗时:", round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "秒\n")

# =========================
# 4) 求逆 → 偏相关矩阵 (带ridge正则化)
# =========================
cat("\n==> Step 4) 矩阵求逆 → 偏相关...\n")
t_step <- Sys.time()

ridge <- RIDGE_INIT
Kmat <- NULL

repeat {
  ok <- TRUE
  Kmat <- tryCatch(solve(R + diag(ridge, k)),
                   error = function(e) { ok <<- FALSE; NULL })
  if (ok && all(is.finite(Kmat))) break
  ridge <- ridge * 10
  if (ridge > RIDGE_MAX) stop("矩阵求逆失败, ridge已达上限。")
}
cat("   ridge =", format(ridge, scientific = TRUE), "\n")

dK <- diag(Kmat)
den <- sqrt(tcrossprod(dK, dK))
pcor_mat <- -Kmat / den
diag(pcor_mat) <- 1
rm(Kmat, den, dK, R); gc()

dimnames(pcor_mat) <- list(code_by_id, code_by_id)
cat("   偏相关矩阵完成! 耗时:",
    round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "秒\n")

# =========================
# 5) P值 + FDR校正
# =========================
cat("\n==> Step 5) P值 + FDR校正...\n")

df <- n - k - 2
if (df <= 1) stop("自由度不足 (df <= 1)。")
cat("   df =", df, "\n")

ut <- which(upper.tri(pcor_mat), arr.ind = TRUE)
r_vec <- pcor_mat[ut]
r_vec <- pmin(pmax(r_vec, -0.999999), 0.999999)

t_vec     <- r_vec * sqrt(df / (1 - r_vec^2))
p_vec     <- 2 * pt(-abs(t_vec), df = df)
p_adj_vec <- p.adjust(p_vec, method = P_ADJ_METHOD)

edge_dt <- data.table(
  i       = ut[, 1],
  j       = ut[, 2],
  from    = code_by_id[ut[, 1]],
  to      = code_by_id[ut[, 2]],
  weight  = as.numeric(r_vec),
  p_value = as.numeric(p_vec),
  p_adj   = as.numeric(p_adj_vec)
)
rm(ut, r_vec, t_vec, p_vec, p_adj_vec); gc()

cat("   总疾病对:", format(nrow(edge_dt), big.mark = ","), "\n")
cat("   正偏相关且FDR<0.05:", format(nrow(edge_dt[weight > 0 & p_adj < 0.05]), big.mark = ","), "\n")

if (SAVE_ALL_PAIRS_CSV) {
  fwrite(edge_dt, file.path(OUT_DIR, "partial_correlation_all_pairs_uppertri.csv"))
  cat("   全量疾病对CSV已保存\n")
}

# =========================
# 6) MMC (三个版本)
# =========================
cat("\n==> Step 6) MMC计算...\n")

# 6A) 无阈值版: 全矩阵正偏相关之和
tmp <- pcor_mat; diag(tmp) <- 0
mmc_all <- rowSums(tmp * (tmp > 0), na.rm = TRUE)
rm(tmp); gc()

mmc_all_df <- data.table(
  DiseaseCode = code_by_id, ICD10_Chapter = chapter_by_id,
  Frequency = freq_by_id, Prevalence = prev_by_id,
  MMC = as.numeric(mmc_all)
)
setorder(mmc_all_df, -MMC)

# 6B) coef>=0.01, FDR<0.05
edges_0105 <- edge_dt[p_adj < 0.05 & weight >= 0.01, .(i, j, weight)]
mmc_0105_dt <- rbindlist(list(
  edges_0105[, .(disease_id = i, w = weight)],
  edges_0105[, .(disease_id = j, w = weight)]
))[, .(MMC = sum(w, na.rm = TRUE)), by = disease_id]
mmc_0105 <- rep(0, k)
mmc_0105[mmc_0105_dt$disease_id] <- mmc_0105_dt$MMC

mmc_0105_df <- data.table(
  DiseaseCode = code_by_id, ICD10_Chapter = chapter_by_id,
  Frequency = freq_by_id, Prevalence = prev_by_id,
  MMC = as.numeric(mmc_0105)
)
setorder(mmc_0105_df, -MMC)

# 6C) coef>0, FDR<0.05
edges_0005 <- edge_dt[p_adj < 0.05 & weight > 0, .(i, j, weight)]
mmc_0005_dt <- rbindlist(list(
  edges_0005[, .(disease_id = i, w = weight)],
  edges_0005[, .(disease_id = j, w = weight)]
))[, .(MMC = sum(w, na.rm = TRUE)), by = disease_id]
mmc_0005 <- rep(0, k)
mmc_0005[mmc_0005_dt$disease_id] <- mmc_0005_dt$MMC

mmc_0005_df <- data.table(
  DiseaseCode = code_by_id, ICD10_Chapter = chapter_by_id,
  Frequency = freq_by_id, Prevalence = prev_by_id,
  MMC = as.numeric(mmc_0005)
)
setorder(mmc_0005_df, -MMC)

fwrite(mmc_all_df,  file.path(OUT_DIR, "MMC_no_threshold.csv"))
fwrite(mmc_0105_df, file.path(OUT_DIR, "MMC_coef0.01_alpha0.05.csv"))
fwrite(mmc_0005_df, file.path(OUT_DIR, "MMC_coef0_alpha0.05.csv"))

cat("   Top 10 Hub (无阈值MMC):\n")
print(mmc_all_df[1:10, .(DiseaseCode, MMC = round(MMC, 3), Prevalence = round(Prevalence, 4))])

# =========================
# 7) 网络指标函数
# =========================
calc_metrics <- function(g, coef_thr, alpha_thr, k_total) {
  if (is.null(g) || vcount(g) == 0 || ecount(g) == 0) {
    return(data.table(
      coef = coef_thr, alpha = alpha_thr, nodes_total = k_total,
      nodes = 0L, isolates = 0L, edges = 0L, density = NA_real_,
      components = NA_integer_, giant_nodes = NA_integer_, giant_edges = NA_integer_,
      avg_degree = NA_real_, median_degree = NA_real_, max_degree = NA_real_,
      avg_strength = NA_real_, median_strength = NA_real_, max_strength = NA_real_,
      transitivity_global = NA_real_, assortativity_degree = NA_real_,
      modularity_louvain = NA_real_,
      avg_path_unweighted_giant = NA_real_, diameter_unweighted_giant = NA_real_,
      avg_path_weighted_giant = NA_real_, diameter_weighted_giant = NA_real_,
      degree_centralization = NA_real_, betweenness_centralization = NA_real_,
      edge_w_mean = NA_real_, edge_w_median = NA_real_,
      edge_w_min = NA_real_, edge_w_max = NA_real_, edge_w_sum = NA_real_
    ))
  }
  
  deg <- degree(g); iso <- sum(deg == 0)
  comps <- components(g)
  giant_comp <- which.max(comps$csize)
  vids_giant <- which(comps$membership == giant_comp)
  g_giant <- induced_subgraph(g, vids = vids_giant)
  
  w <- E(g)$weight
  strength_v <- strength(g, weights = w)
  
  modv <- tryCatch({ clu <- cluster_louvain(g, weights = w); modularity(clu) },
                   error = function(e) NA_real_)
  
  apl_u <- NA_real_; dia_u <- NA_real_; apl_w <- NA_real_; dia_w <- NA_real_
  if (vcount(g_giant) >= 2 && ecount(g_giant) >= 1) {
    apl_u <- tryCatch(average.path.length(g_giant, directed = FALSE), error = function(e) NA_real_)
    dia_u <- tryCatch(diameter(g_giant, directed = FALSE), error = function(e) NA_real_)
    ww <- E(g_giant)$weight
    dist_w <- 1 / pmax(ww, 1e-12)
    apl_w <- tryCatch(average.path.length(g_giant, directed = FALSE, weights = dist_w), error = function(e) NA_real_)
    dia_w <- tryCatch(diameter(g_giant, directed = FALSE, weights = dist_w), error = function(e) NA_real_)
  }
  
  deg_cent <- tryCatch(centr_degree(g, mode = "all", normalized = TRUE)$centralization, error = function(e) NA_real_)
  bet_cent <- tryCatch(centr_betw(g, directed = FALSE, normalized = TRUE)$centralization, error = function(e) NA_real_)
  
  data.table(
    coef = coef_thr, alpha = alpha_thr, nodes_total = k_total,
    nodes = vcount(g), isolates = iso, edges = ecount(g),
    density = edge_density(g, loops = FALSE),
    components = comps$no,
    giant_nodes = vcount(g_giant), giant_edges = ecount(g_giant),
    avg_degree = mean(deg), median_degree = as.numeric(median(deg)), max_degree = max(deg),
    avg_strength = mean(strength_v), median_strength = as.numeric(median(strength_v)),
    max_strength = max(strength_v),
    transitivity_global = tryCatch(transitivity(g, type = "global"), error = function(e) NA_real_),
    assortativity_degree = tryCatch(assortativity_degree(g, directed = FALSE), error = function(e) NA_real_),
    modularity_louvain = modv,
    avg_path_unweighted_giant = apl_u, diameter_unweighted_giant = dia_u,
    avg_path_weighted_giant = apl_w, diameter_weighted_giant = dia_w,
    degree_centralization = deg_cent, betweenness_centralization = bet_cent,
    edge_w_mean = mean(w), edge_w_median = as.numeric(median(w)),
    edge_w_min = min(w), edge_w_max = max(w), edge_w_sum = sum(w)
  )
}

# =========================
# 8) 网络绘图函数
# =========================
plot_one_network <- function(g, title_text) {
  base_theme <- theme_void(base_size = 10) +
    theme(
      text = element_text(family = FONT_FAMILY),
      plot.title = element_text(hjust = 0.5, family = FONT_FAMILY),
      legend.title = element_text(size = 14, family = FONT_FAMILY),
      legend.text  = element_text(size = 12, family = FONT_FAMILY),
      legend.key.width  = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.spacing.y  = unit(0.4, "cm")
    )
  
  if (is.null(g) || vcount(g) == 0 || ecount(g) == 0) {
    return(ggplot() + base_theme +
             geom_text(aes(0, 0, label = "No edges under this threshold"),
                       family = FONT_FAMILY, size = 6) +
             labs(title = title_text))
  }
  
  iso_v <- which(degree(g) == 0)
  g_plot <- if (length(iso_v) > 0) delete_vertices(g, iso_v) else g
  
  if (vcount(g_plot) == 0 || ecount(g_plot) == 0) {
    return(ggplot() + base_theme +
             geom_text(aes(0, 0, label = "Only isolates left"),
                       family = FONT_FAMILY, size = 6) +
             labs(title = title_text))
  }
  
  ggraph(g_plot, layout = "fr") +
    geom_edge_link(aes(edge_width = abs(weight)), alpha = 0.9, color = "grey60") +
    geom_node_point(aes(color = chapter, size = prevalence), alpha = 0.85) +
    { if (SHOW_NODE_LABELS) geom_node_text(aes(label = name), family = FONT_FAMILY,
                                           color = "black", size = NODE_LABEL_SIZE) } +
    scale_edge_width(range = c(0.1, 2.5), name = "Partial correlation (\u03C6)",
                     guide = guide_legend(order = 1,
                                          override.aes = list(color = "grey50", linetype = "solid", shape = NA))) +
    scale_color_manual(name = "ICD-10 Chapter", values = chapter_colors_provided,
                       breaks = as.character(1:14),
                       guide = guide_legend(order = 2,
                                            override.aes = list(size = 6, alpha = 1),
                                            keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))) +
    scale_size_continuous(name = "Prevalence", range = c(5, 12), guide = "none") +
    base_theme + labs(title = title_text)
}

# =========================
# 9) 多阈值循环: 网络图 + 边列表 + 指标
# =========================
cat("\n==> Step 9) 多阈值网络构建...\n")

vertex_df_full <- data.frame(
  name = code_by_id,
  chapter = factor(chapter_by_id, levels = c(as.character(1:14), "0")),
  prevalence = prev_by_id,
  stringsAsFactors = FALSE
)

metrics_list <- vector("list", nrow(THRESH_GRID))

for (idx in 1:nrow(THRESH_GRID)) {
  coef_thr  <- THRESH_GRID[idx, coef]
  alpha_thr <- THRESH_GRID[idx, alpha]
  
  edges_sub <- edge_dt[p_adj < alpha_thr & weight >= coef_thr,
                       .(from, to, weight, p_value, p_adj)]
  
  fwrite(edges_sub, file.path(OUT_DIR,
                              sprintf("edge_list_coef_%.2f_alpha_%.2f.csv", coef_thr, alpha_thr)))
  
  g <- graph_from_data_frame(edges_sub[, .(from, to, weight)],
                             directed = FALSE, vertices = vertex_df_full)
  
  metrics_list[[idx]] <- calc_metrics(g, coef_thr, alpha_thr, k_total = k)
  
  title_text <- sprintf("Multimorbidity network (partial correlation) | coef >= %.2f, FDR < %.2f",
                        coef_thr, alpha_thr)
  p <- plot_one_network(g, title_text)
  
  pdf_file <- file.path(OUT_DIR,
                        sprintf("network_coef_%.2f_alpha_%.2f.pdf", coef_thr, alpha_thr))
  ggsave(pdf_file, p, width = PDF_W, height = PDF_H, units = "in", device = grDevices::cairo_pdf)
  cat("   ", basename(pdf_file), "  edges:", nrow(edges_sub), "\n")
}

metrics_dt <- rbindlist(metrics_list, fill = TRUE)

# =========================
# 10) Excel输出
# =========================
cat("\n==> Step 10) Excel输出...\n")

wb <- createWorkbook()
headerStyle <- createStyle(
  fontSize = 11, fontColour = "#FFFFFF", halign = "center", valign = "center",
  fgFill = "#4472C4", textDecoration = "bold"
)

add_sheet <- function(sheet, dt_in) {
  sn <- substr(sheet, 1, 31)
  addWorksheet(wb, sn)
  writeDataTable(wb, sn, dt_in, withFilter = TRUE)
  setColWidths(wb, sn, cols = 1:ncol(dt_in), widths = "auto")
  addStyle(wb, sn, headerStyle, rows = 1, cols = 1:ncol(dt_in), gridExpand = TRUE)
}

add_sheet("network_metrics", metrics_dt)
add_sheet("MMC_no_threshold", mmc_all_df)
add_sheet("MMC_coef0.01_alpha0.05", mmc_0105_df)
add_sheet("MMC_coef0_alpha0.05", mmc_0005_df)

xlsx_file <- file.path(OUT_DIR, "network_metrics_and_MMC.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   已保存:", xlsx_file, "\n")

# =========================
# 11) 保存RDS中间文件 (供社区发现脚本使用)
# =========================
cat("\n==> Step 11) 保存RDS...\n")

saveRDS(disease_map, file.path(OUT_DIR, "disease_map.rds"))
saveRDS(edge_dt,     file.path(OUT_DIR, "edges_all_uppertri_with_fdr.rds"))
saveRDS(pcor_mat,    file.path(OUT_DIR, "pcor_matrix.rds"))
saveRDS(metrics_dt,  file.path(OUT_DIR, "network_metrics.rds"))

cat("   RDS已保存到:", OUT_DIR, "\n")
cat("   → 社区发现脚本的 RDS_DIR 指向此目录即可\n")

# =========================
# 12) 摘要
# =========================
total_time <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

cat("\n==============================================================================\n")
cat("  Step 2 完成  (总用时:", total_time, "分钟)\n")
cat("==============================================================================\n\n")

cat("算法: sparse crossprod → P×P相关 → 求逆(ridge) → 偏相关\n")
cat("  ridge =", format(ridge, scientific = TRUE), "\n")
cat("  疾病数 (患病率>", PREV_THRESHOLD, "):", k, "\n")
cat("  患者数:", format(n, big.mark = ","), "\n\n")

cat("多阈值网络统计:\n")
print(metrics_dt[, .(coef, alpha, nodes, edges, density = round(density, 4),
                     modularity = round(modularity_louvain, 3))])

cat("\nTop 10 Hub (无阈值MMC):\n")
print(mmc_all_df[1:10, .(DiseaseCode, MMC = round(MMC, 3), Prevalence = round(Prevalence, 4))])

cat("\n输出目录:", OUT_DIR, "\n")
cat("  → 下一步: 修改社区发现脚本的 RDS_DIR 指向此目录\n")
cat("  → 运行 Leiden_Resolution_Sweep.R 确定最优分辨率\n")
cat("  → 运行 总体-社区发现_update.R\n")
cat("==============================================================================\n")