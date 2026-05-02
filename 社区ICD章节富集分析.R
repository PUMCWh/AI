
社区ICD章节富集分析完整结果
#
# 输入：
#   - disease_map.rds          （来自总体网络分析输出）
#   - consensus_results_all.rds（来自社区发现脚本输出）
# 输出：
#   - 社区ICD章节富集热图.pdf
#   - 社区ICD章节富集分析.xlsx
#
# 方法：对每个社区×每个ICD-10章节，构建2×2列联表，
#       以Fisher精确检验判断是否显著富集/耗竭，
#       并计算odds ratio。热图填充色为log2(O/E)，星号标注FDR水平。


rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

# =========================
# 0) Paths & Params
# =========================
# --- 请根据实际路径修改 ---
RDS_DIR  <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\总体"
COMM_DIR <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\社区检测"
OUT_DIR  <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\社区检测"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 主分析阈值标签（与社区发现脚本一致）
MAIN_TAG <- "coef0.01_alpha0.05"

# FDR阈值
FDR_THRESH <- 0.05

# 图参数
FONT_FAMILY <- "Arial"
PDF_W <- 12
PDF_H <- 8

# PDF device
PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# ICD-10 章节标签（罗马数字）
chapter_labels <- c(
  "1"  = "I",    "2"  = "II",   "3"  = "III",  "4"  = "IV",
  "5"  = "V",    "6"  = "VI",   "7"  = "VII",  "8"  = "VIII",
  "9"  = "IX",   "10" = "X",    "11" = "XI",   "12" = "XII",
  "13" = "XIII", "14" = "XIV"
)

# =========================
# 1) Load data
# =========================
cat("==> Step 1) Loading data...\n")

# --- disease_map ---
rds_candidates <- list.dirs(RDS_DIR, recursive = FALSE, full.names = TRUE)
rds_candidates <- rds_candidates[grepl("NETWORK_PCOR_", basename(rds_candidates))]
rds_folder <- if (length(rds_candidates) > 0) sort(rds_candidates, decreasing = TRUE)[1] else RDS_DIR

disease_map <- readRDS(file.path(rds_folder, "disease_map.rds"))
disease_map_dt <- as.data.table(disease_map)

# 自动识别code列名（兼容不同数据源）
code_col <- intersect(names(disease_map_dt),
                      c("code", "DiseaseCode", "node", "icd3", "ICD3", "ICD10_3", "ICD10"))[1]
if (is.na(code_col)) stop("Cannot find code column in disease_map.rds. Available: ",
                          paste(names(disease_map_dt), collapse = ", "))
if (code_col != "code") {
  setnames(disease_map_dt, code_col, "code")
  cat("   Renamed column '", code_col, "' -> 'code'\n")
}

# 确保ICD10_Chapter列存在
if (!"ICD10_Chapter" %in% names(disease_map_dt) && "chapter" %in% names(disease_map_dt)) {
  setnames(disease_map_dt, "chapter", "ICD10_Chapter")
}
disease_map_dt[, ICD10_Chapter := as.character(ICD10_Chapter)]
disease_map_dt <- disease_map_dt[!is.na(ICD10_Chapter) & ICD10_Chapter != "" & ICD10_Chapter != "0"]

cat("   Loaded disease_map:", nrow(disease_map_dt), "diseases\n")

# --- medoid partition ---
consensus_results <- readRDS(file.path(COMM_DIR, "consensus_results_all.rds"))
medoid_dt <- consensus_results[[MAIN_TAG]]$medoid_dt
medoid_dt <- medoid_dt[community > 0]  # 排除孤立节点
cat("   Loaded medoid partition:", nrow(medoid_dt), "nodes in",
    length(unique(medoid_dt$community)), "communities\n")

# =========================
# 2) Merge membership + chapter
# =========================
cat("==> Step 2) Merging membership with ICD chapter...\n")

merged <- merge(medoid_dt[, .(node, community)],
                disease_map_dt[, .(code, ICD10_Chapter)],
                by.x = "node", by.y = "code", all.x = TRUE)
# 去除无章节信息的节点
merged <- merged[!is.na(ICD10_Chapter) & ICD10_Chapter != "" & ICD10_Chapter != "0"]

all_comms    <- sort(unique(merged$community))
all_chapters <- sort(unique(merged$ICD10_Chapter), decreasing = FALSE)
N_total      <- nrow(merged)  # 总节点数（有章节信息的）

cat("   Merged:", N_total, "nodes,", length(all_comms), "communities,",
    length(all_chapters), "chapters\n")

# =========================
# 3) Fisher exact test: community × chapter
# =========================
cat("==> Step 3) Running Fisher exact tests...\n")

results <- list()
idx <- 0

for (cid in all_comms) {
  nodes_in_comm <- merged[community == cid, node]
  n_comm <- length(nodes_in_comm)
  
  for (ch in all_chapters) {
    # 2×2 table:
    #              Chapter=ch  Chapter≠ch
    # In community     a           b       | n_comm
    # Not in comm      c           d       | N_total - n_comm
    #                  n_ch       N-n_ch   | N_total
    
    a <- sum(merged[community == cid, ICD10_Chapter == ch])
    b <- n_comm - a
    n_ch <- sum(merged$ICD10_Chapter == ch)
    c_val <- n_ch - a
    d <- (N_total - n_comm) - c_val
    
    # Fisher exact test (two-sided)
    mat_2x2 <- matrix(c(a, b, c_val, d), nrow = 2, byrow = TRUE,
                      dimnames = list(InComm = c("Yes", "No"),
                                      InChapter = c("Yes", "No")))
    ft <- fisher.test(mat_2x2)
    
    # 安全提取OR（保留raw用于方向判定，round仅用于展示）
    or_raw <- suppressWarnings(as.numeric(ft$estimate))
    if (length(or_raw) == 0) or_raw <- NA_real_
    
    ci_lo <- suppressWarnings(as.numeric(ft$conf.int[1]))
    ci_hi <- suppressWarnings(as.numeric(ft$conf.int[2]))
    if (length(ci_lo) == 0) ci_lo <- NA_real_
    if (length(ci_hi) == 0) ci_hi <- NA_real_
    
    # 计算observed/expected比值
    expected <- n_comm * n_ch / N_total
    obs_exp_ratio <- if (expected > 0) a / expected else NA_real_
    
    idx <- idx + 1
    results[[idx]] <- data.table(
      community       = cid,
      chapter         = ch,
      n_in_comm_ch    = a,           # 该社区中属于该章节的疾病数
      n_in_comm       = n_comm,      # 该社区总疾病数
      n_in_chapter    = n_ch,        # 整个网络中该章节总疾病数
      n_total         = N_total,
      pct_in_comm     = round(a / n_comm * 100, 1),  # 社区内该章节占比
      pct_in_network  = round(n_ch / N_total * 100, 1),  # 网络中该章节占比
      or_raw          = or_raw,              # 原始OR，用于方向判定
      odds_ratio      = round(or_raw, 3),    # 展示用OR
      obs_exp_ratio   = round(obs_exp_ratio, 3),
      p_value         = ft$p.value,
      ci_low          = round(ci_lo, 3),
      ci_high         = round(ci_hi, 3)
    )
  }
}

enrich_dt <- rbindlist(results)

# FDR校正（BH法，全部测试一起校正）
enrich_dt[, fdr := p.adjust(p_value, method = "BH")]

# 方向标记（基于原始OR，未round，避免边界误判）
enrich_dt[, direction := fifelse(is.na(or_raw), "unknown",
                                 fifelse(or_raw > 1, "enriched",
                                         fifelse(or_raw < 1, "depleted", "neutral")))]
enrich_dt[, significant := fdr < FDR_THRESH]

cat("   Total tests:", nrow(enrich_dt), "\n")
cat("   Significant (FDR <", FDR_THRESH, "):",
    sum(enrich_dt$significant), "of which",
    sum(enrich_dt$significant & enrich_dt$direction == "enriched"), "enriched,",
    sum(enrich_dt$significant & enrich_dt$direction == "depleted"), "depleted\n")

# 背景口径检查（严格一致，避免图注与实际不符）
if (N_total != 607) {
  stop("Background N_total = ", N_total, " (expected 607). ",
       "Check disease_map chapter mapping / merge key.")
}
cat("   Background N_total:", N_total, "\u2714\n")

# =========================
# 4) Create heatmap (图31)
# =========================
cat("==> Step 4) Creating enrichment heatmap...\n")

# 热图填充：log2(O/E)——展示效应大小，星号展示显著性
# 富集 → 正值（红色），耗竭 → 负值（蓝色）
enrich_dt[, log2_oe := fifelse(obs_exp_ratio > 0 & is.finite(obs_exp_ratio),
                               log2(obs_exp_ratio), NA_real_)]

# 不显著的设为0（白色），仅显著结果着色
# 注意：O/E=0（观察数=0）时 log2_oe=NA，需用 -Inf 标记为最强耗竭
enrich_dt[, log2_oe_plot := fifelse(!significant, 0,
                                    fifelse(!is.na(log2_oe), log2_oe, NA_real_))]
# 暂存NA（O/E=0的显著耗竭），后面用 -max_abs 替换

# 准备绘图数据
plot_dt <- copy(enrich_dt)

# 章节标签
plot_dt[, chapter_label := chapter_labels[chapter]]
plot_dt[is.na(chapter_label), chapter_label := paste0("Ch.", chapter)]

# 社区标签（按编号排序）
plot_dt[, comm_label := factor(paste0("C", community),
                               levels = paste0("C", sort(all_comms)))]

# 章节排序（I在上方，XIV在下方，符合阅读习惯）
ch_order <- chapter_labels[as.character(sort(as.integer(all_chapters)))]
plot_dt[, chapter_label := factor(chapter_label, levels = rev(ch_order))]
# 注：ggplot y轴从下到上，rev()使 I 在最上方

# 显著性星标
plot_dt[, sig_label := fifelse(significant & fdr < 0.001, "***",
                               fifelse(significant & fdr < 0.01,  "**",
                                       fifelse(significant & fdr < 0.05,  "*", "")))]

# 格内数字：观察数（该社区中属于该章节的实际疾病数）
plot_dt[, cell_text := paste0(n_in_comm_ch)]

# 色阶范围对称化
max_abs <- max(abs(plot_dt$log2_oe_plot), na.rm = TRUE)
max_abs <- ceiling(max_abs)
if (max_abs < 2) max_abs <- 2

# O/E=0的显著耗竭 → 用-max_abs（色阶最深蓝）替代NA
plot_dt[is.na(log2_oe_plot), log2_oe_plot := -max_abs]

p_heatmap <- ggplot(plot_dt, aes(x = comm_label, y = chapter_label)) +
  # 底层：填充色
  geom_tile(aes(fill = log2_oe_plot), color = "grey80", linewidth = 0.3) +
  # 格内数字：实际疾病数
  geom_text(aes(label = cell_text), size = 3, color = "grey10",
            family = FONT_FAMILY, vjust = -0.3) +
  # 显著性星标
  geom_text(aes(label = sig_label), size = 3.5, color = "black",
            family = FONT_FAMILY, fontface = "bold", vjust = 1.2) +
  # 双向色阶：红=富集，蓝=耗竭，白=不显著
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = expression(log[2](O/E)),
    breaks = seq(-max_abs, max_abs, by = max(1, floor(max_abs / 3)))
  ) +
  labs(
    title = "ICD-10 Chapter Enrichment by Community",
    subtitle = paste0("Fisher exact test, BH-corrected; ",
                      "* FDR<0.05, ** FDR<0.01, *** FDR<0.001; ",
                      "white = FDR \u2265 0.05"),
    x = "Community",
    y = "ICD-10 Chapter"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13, color = "black"),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.title       = element_text(size = 13, color = "black"),
    legend.position  = "right",
    legend.key.height = unit(1.2, "cm"),
    panel.grid       = element_blank()
  )

ggsave(file.path(OUT_DIR, "图31_社区ICD章节富集热图.pdf"),
       p_heatmap, width = PDF_W, height = PDF_H, units = "in", device = PDF_DEVICE)
cat("   Saved: 图31_社区ICD章节富集热图.pdf\n")

# =========================
# 5) 附表S39：完整富集结果
# =========================
cat("==> Step 5) Writing supplementary table S39...\n")

# 整理输出表（仅保留FDR<0.05的显著结果，列精简）
sig_dt <- enrich_dt[significant == TRUE]
setorder(sig_dt, community, -obs_exp_ratio)

out_dt <- sig_dt[, .(
  Community              = community,
  `ICD-10 Chapter`       = chapter_labels[chapter],
  `n (observed)`         = n_in_comm_ch,
  `% in community`       = pct_in_comm,
  `% in network`         = pct_in_network,
  `O/E`                  = obs_exp_ratio,
  `OR (95% CI)`          = {
    # 格式化OR和CI，处理Inf/0/NA
    or_str <- fifelse(is.na(odds_ratio), "\u2014",
                      fifelse(!is.finite(odds_ratio), "\u221e",
                              sprintf("%.2f", odds_ratio)))
    lo_str <- fifelse(is.na(ci_low), "\u2014",
                      fifelse(!is.finite(ci_low), "0.00",
                              sprintf("%.2f", ci_low)))
    hi_str <- fifelse(is.na(ci_high), "\u2014",
                      fifelse(!is.finite(ci_high), "\u221e",
                              sprintf("%.2f", ci_high)))
    fifelse(is.na(odds_ratio), "\u2014",
            paste0(or_str, " (", lo_str, "\u2013", hi_str, ")"))
  },
  FDR                    = signif(fdr, 3)
)]

cat("   Significant results for S39:", nrow(out_dt), "rows\n")
setorder(out_dt, Community, `ICD-10 Chapter`)

# Excel输出
wb <- createWorkbook()
headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                           valign = "center", fgFill = "#4472C4", textDecoration = "bold")

addWorksheet(wb, "Significant_Enrichment")
writeDataTable(wb, "Significant_Enrichment", out_dt, withFilter = TRUE)
setColWidths(wb, "Significant_Enrichment", cols = 1:ncol(out_dt), widths = "auto")
addStyle(wb, "Significant_Enrichment", headerStyle,
         rows = 1, cols = 1:ncol(out_dt), gridExpand = TRUE)

# 条件格式：O/E>1 富集浅红，O/E<1 耗竭浅蓝
enrichStyle  <- createStyle(fgFill = "#FCE4EC")
depletedStyle <- createStyle(fgFill = "#E3F2FD")
for (r in seq_len(nrow(out_dt))) {
  sty <- if (out_dt[r, `O/E`] > 1) enrichStyle else depletedStyle
  addStyle(wb, "Significant_Enrichment", sty, rows = r + 1,
           cols = 1:ncol(out_dt), gridExpand = TRUE, stack = TRUE)
}

# 第二个sheet：完整308行结果（备查，不纳入正文附表页数）
full_dt <- enrich_dt[, .(
  Community              = community,
  `ICD-10 Chapter`       = chapter_labels[chapter],
  `n (observed)`         = n_in_comm_ch,
  `Community size`       = n_in_comm,
  `Chapter size`         = n_in_chapter,
  `Network size`         = n_total,
  `% in community`       = pct_in_comm,
  `% in network`         = pct_in_network,
  `O/E`                  = obs_exp_ratio,
  `log2(O/E)`            = round(log2_oe, 3),
  OR                     = odds_ratio,
  `CI lower`             = ci_low,
  `CI upper`             = ci_high,
  `P value`              = signif(p_value, 4),
  FDR                    = signif(fdr, 4),
  `FDR < 0.05`           = significant
)]
setorder(full_dt, Community, `ICD-10 Chapter`)

addWorksheet(wb, "Full_Results_308")
writeDataTable(wb, "Full_Results_308", full_dt, withFilter = TRUE)
setColWidths(wb, "Full_Results_308", cols = 1:ncol(full_dt), widths = "auto")
addStyle(wb, "Full_Results_308", headerStyle,
         rows = 1, cols = 1:ncol(full_dt), gridExpand = TRUE)

# 汇总页：每个社区的显著富集章节
sig_summary <- out_dt[`O/E` > 1,
                      .(Enriched_Chapters = paste(`ICD-10 Chapter`, collapse = "; "),
                        N_enriched = .N),
                      by = Community]
setorder(sig_summary, Community)

addWorksheet(wb, "Summary_Enriched")
writeDataTable(wb, "Summary_Enriched", sig_summary, withFilter = TRUE)
setColWidths(wb, "Summary_Enriched", cols = 1:ncol(sig_summary), widths = "auto")
addStyle(wb, "Summary_Enriched", headerStyle,
         rows = 1, cols = 1:ncol(sig_summary), gridExpand = TRUE)

xlsx_file <- file.path(OUT_DIR, "附表S39_社区ICD章节富集分析.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", xlsx_file, "\n")

# =========================
# 6) Console summary
# =========================
cat("\n========================================\n")
cat("DONE. Enrichment analysis complete.\n")
cat("========================================\n\n")

# 简要统计
cat("Community-level summary of significant enrichments:\n")
for (cid in all_comms) {
  sig_ch <- enrich_dt[community == cid & significant == TRUE & direction == "enriched"]
  if (nrow(sig_ch) > 0) {
    ch_str <- paste(sprintf("Ch.%s (OR=%.1f)", sig_ch$chapter, sig_ch$odds_ratio),
                    collapse = ", ")
    cat(sprintf("  C%d (n=%d): %s\n", cid,
                unique(enrich_dt[community == cid, n_in_comm]),
                ch_str))
  } else {
    cat(sprintf("  C%d (n=%d): no significant enrichment\n", cid,
                unique(enrich_dt[community == cid, n_in_comm])))
  }
}

cat("\nOutputs:\n")
cat("  1) 图31_社区ICD章节富集热图.pdf\n")
cat("  2) 附表S39_社区ICD章节富集分析.xlsx\n")
cat("     Sheet 1: Significant_Enrichment (FDR<0.05, for print)\n")
cat("     Sheet 2: Full_Results_308 (all 308 tests, for reviewer)\n")
cat("     Sheet 3: Summary_Enriched (per-community summary)\n")