# ==============================================================================
# 桥接疾病的节点角色分类（z-P散点图）
# 参与系数排名前20的桥接疾病
# 桥接节点完整列表
#
# 输入：
#   - medoid_bridge_metrics_coef0.01_alpha0.05.csv（社区发现脚本输出）
#   - disease_map.rds（疾病名称映射）
#
# 输出：
#   - 节点角色分类zP散点图.pdf
#   - 桥接疾病Top20.xlsx
#   - 桥接节点完整列表.xlsx
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
# 0) Paths & Params
# =========================
RDS_DIR  <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\总体"
COMM_DIR <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\社区检测"
OUT_DIR  <- "C:\\Users\\HP\\Desktop\\PHD\\结果Final\\网络分析\\社区检测"

MAIN_TAG    <- "coef0.01_alpha0.05"
FONT_FAMILY <- "Arial"
PDF_W <- 12
PDF_H <- 8

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# Guimerà–Amaral 阈值
Z_THRESH  <- 2.5   # hub vs non-hub
P_THRESH1 <- 0.30  # peripheral vs connector
P_THRESH2 <- 0.62  # connector vs kinless
P_THRESH3 <- 0.80  # kinless boundary (可选)

# =========================
# 1) Load data
# =========================
cat("==> Step 1) Loading data...\n")

# 桥接指标
bridge_file <- file.path(COMM_DIR, sprintf("medoid_bridge_metrics_%s.csv", MAIN_TAG))
if (!file.exists(bridge_file)) {
  stop("Bridge metrics file not found: ", bridge_file)
}
bridge_dt <- fread(bridge_file)
cat("   Loaded:", nrow(bridge_dt), "nodes\n")
cat("   Bridge CSV columns:", paste(names(bridge_dt), collapse = ", "), "\n")

# 疾病名称映射
rds_candidates <- list.dirs(RDS_DIR, recursive = FALSE, full.names = TRUE)
rds_candidates <- rds_candidates[grepl("NETWORK_PCOR_", basename(rds_candidates))]
rds_folder <- if (length(rds_candidates) > 0) sort(rds_candidates, decreasing = TRUE)[1] else RDS_DIR

disease_map <- as.data.table(readRDS(file.path(rds_folder, "disease_map.rds")))
cat("   disease_map columns:", paste(names(disease_map), collapse = ", "), "\n")

# 自动识别疾病名称列（bridge CSV可能已有，disease_map可能也有）
name_col_bridge <- intersect(names(bridge_dt),
                             c("disease_name", "DiseaseName", "name", "label",
                               "疾病名称", "疾病名称_cn", "name_cn"))[1]
if (!is.na(name_col_bridge)) {
  if (name_col_bridge != "disease_name") setnames(bridge_dt, name_col_bridge, "disease_name")
  cat("   Disease names found in bridge CSV.\n")
} else {
  # 从disease_map补充
  name_col <- intersect(names(disease_map),
                        c("disease_name", "DiseaseName", "name", "label", "description",
                          "疾病名称", "疾病名称_cn", "name_cn", "ICD_name",
                          "disease_label", "Name"))[1]
  code_col <- intersect(names(disease_map),
                        c("code", "DiseaseCode", "node", "icd3", "ICD3"))[1]
  if (is.na(code_col)) stop("Cannot find code column in disease_map")
  
  if (!is.na(name_col)) {
    name_map <- disease_map[, .(code = get(code_col), disease_name = get(name_col))]
    bridge_dt <- merge(bridge_dt, name_map, by.x = "node", by.y = "code", all.x = TRUE)
    cat("   Disease names merged from disease_map.\n")
  } else {
    bridge_dt[, disease_name := node]
    cat("   No disease_name column found, using ICD-3 codes as labels.\n")
  }
}

# 确保ICD10_Chapter存在（bridge CSV通常已有）
if (!"ICD10_Chapter" %in% names(bridge_dt)) {
  ch_col <- intersect(names(disease_map), c("ICD10_Chapter", "chapter"))[1]
  code_col2 <- intersect(names(disease_map), c("code", "DiseaseCode", "node"))[1]
  if (!is.na(ch_col) && !is.na(code_col2)) {
    ch_map <- disease_map[, .(code = get(code_col2), ICD10_Chapter = get(ch_col))]
    bridge_dt <- merge(bridge_dt, ch_map, by.x = "node", by.y = "code", all.x = TRUE)
  }
}

# 仅保留社区内活跃节点
active <- bridge_dt[community > 0 & !is.na(participation_coef)]
cat("   Active nodes in communities:", nrow(active), "\n")

# =========================
# 2) 图33：z-P散点图
# =========================
cat("\n==> Step 2) Creating z-P scatter plot (Figure 33)...\n")

# 节点角色分类（确认已有或重新计算）
if (!"role" %in% names(active)) {
  active[, role := fcase(
    within_module_z_score >= Z_THRESH & participation_coef >= P_THRESH1, "connector_hub",
    within_module_z_score >= Z_THRESH & participation_coef <  P_THRESH1, "provincial_hub",
    within_module_z_score <  Z_THRESH & participation_coef >= P_THRESH2, "kinless",
    within_module_z_score <  Z_THRESH & participation_coef >= P_THRESH1, "connector",
    within_module_z_score <  Z_THRESH & participation_coef <  P_THRESH1, "peripheral",
    default = "unclassified"
  )]
}

# 角色标签和颜色
role_labels <- c(
  "provincial_hub" = "Provincial hub",
  "connector_hub"  = "Connector hub",
  "kinless"        = "Kinless",
  "connector"      = "Connector",
  "peripheral"     = "Peripheral",
  "unclassified"   = "Unclassified"
)
role_colors <- c(
  "Provincial hub" = "#E41A1C",
  "Connector hub"  = "#FF7F00",
  "Kinless"        = "#984EA3",
  "Connector"      = "#377EB8",
  "Peripheral"     = "#BDBDBD",
  "Unclassified"   = "#9E9E9E"
)

active[is.na(role), role := "unclassified"]
active[, role_label := role_labels[role]]
active[is.na(role_label), role_label := "Unclassified"]
active[, role_label := factor(role_label,
                              levels = c("Provincial hub", "Connector hub",
                                         "Kinless", "Connector", "Peripheral",
                                         "Unclassified"))]

# 角色计数（确保全部类都有，即使n=0）
all_roles <- data.table(role_label = factor(
  c("Provincial hub", "Connector hub", "Kinless", "Connector", "Peripheral", "Unclassified"),
  levels = c("Provincial hub", "Connector hub", "Kinless", "Connector", "Peripheral", "Unclassified")))
role_counts <- merge(all_roles, active[, .N, by = role_label], by = "role_label", all.x = TRUE)
role_counts[is.na(N), N := 0L]

# 只保留实际出现的角色（避免图例出现"Unclassified (n=0)"）
keep_levels <- role_counts[N > 0, as.character(role_label)]
active[, role_label := factor(as.character(role_label), levels = keep_levels)]
role_colors_use  <- role_colors[keep_levels]
legend_labels <- setNames(
  paste0(role_counts[N > 0, role_label], " (n=", role_counts[N > 0, N], ")"),
  keep_levels
)

# 标注：hub节点 + 参与系数Top 15
hubs <- active[within_module_z_score >= Z_THRESH]
top_P <- head(active[order(-participation_coef)], 15)
label_nodes <- unique(rbind(hubs[, .(node, disease_name)],
                            top_P[, .(node, disease_name)]))
active[, show_label := fifelse(
  node %in% label_nodes$node,
  fifelse(!is.na(disease_name) & disease_name != node,
          paste0(node, " (", disease_name, ")"),
          node),
  NA_character_
)]
# 过长名称截断
active[!is.na(show_label) & nchar(show_label) > 30,
       show_label := paste0(substr(show_label, 1, 27), "...")]

# 坐标范围
z_max <- max(active$within_module_z_score, na.rm = TRUE) * 1.08
z_min <- min(active$within_module_z_score, na.rm = TRUE) - 0.5

p_scatter <- ggplot(active, aes(x = participation_coef, y = within_module_z_score)) +
  # 分界线
  geom_hline(yintercept = Z_THRESH, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = P_THRESH1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = P_THRESH2, linetype = "dotted", color = "grey70", linewidth = 0.4) +
  # 象限标注
  annotate("text", x = P_THRESH1 / 2, y = z_max * 0.95,
           label = "Provincial\nhub", color = "grey45", size = 4, fontface = "italic",
           family = FONT_FAMILY) +
  annotate("text", x = (P_THRESH1 + P_THRESH2) / 2, y = z_max * 0.95,
           label = "Connector\nhub", color = "grey45", size = 4, fontface = "italic",
           family = FONT_FAMILY) +
  annotate("text", x = P_THRESH1 / 2, y = z_min + 0.3,
           label = "Peripheral", color = "grey45", size = 4, fontface = "italic",
           family = FONT_FAMILY) +
  annotate("text", x = (P_THRESH1 + P_THRESH2) / 2, y = z_min + 0.3,
           label = "Connector", color = "grey45", size = 4, fontface = "italic",
           family = FONT_FAMILY) +
  annotate("text", x = (P_THRESH2 + 1) / 2, y = z_min + 0.3,
           label = "Kinless", color = "grey45", size = 4, fontface = "italic",
           family = FONT_FAMILY) +
  # 散点：Peripheral单独画（微抖避免P=0重叠），其余正常
  geom_point(data = active[role_label == "Peripheral"],
             aes(color = role_label, size = strength), alpha = 0.35,
             position = position_jitter(width = 0.006, height = 0, seed = 42)) +
  geom_point(data = active[role_label != "Peripheral"],
             aes(color = role_label, size = strength), alpha = 0.7) +
  # 标签
  geom_text_repel(data = active[!is.na(show_label)],
                  aes(label = show_label), size = 3, max.overlaps = 25,
                  family = FONT_FAMILY, color = "grey20",
                  segment.color = "grey60", segment.size = 0.3,
                  box.padding = 0.4, point.padding = 0.3,
                  min.segment.length = 0.2, seed = 42) +
  # 色阶
  scale_color_manual(name = "Node role",
                     values = role_colors_use,
                     labels = legend_labels,
                     drop = TRUE) +
  scale_size_continuous(name = "Strength", range = c(1.5, 7),
                        breaks = pretty) +
  # 坐标轴
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(z_min, z_max)) +
  labs(x = "Participation coefficient (P)",
       y = "Within-module z-score (z)") +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12),
    axis.line        = element_blank(),
    axis.ticks       = element_line(color = "black", linewidth = 0.6),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    legend.position  = "right",
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 12),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4)),
         size  = guide_legend(order = 2))

ggsave(file.path(OUT_DIR, "图33_节点角色分类zP散点图.pdf"),
       p_scatter, width = PDF_W, height = PDF_H,
       units = "in", device = PDF_DEVICE)
cat("   Saved: 图33_节点角色分类zP散点图.pdf\n")

# 控制台汇总
cat("\n   Role distribution:\n")
print(active[, .N, by = role_label][order(role_label)])

# =========================
# 3) 表19：参与系数排名前20的桥接疾病
# =========================
cat("\n==> Step 3) Creating Table 19 (Top 20 bridge diseases)...\n")

top20 <- head(active[order(-participation_coef)], 20)

tab19 <- top20[, .(
  Rank                  = 1:20,
  `ICD-3 Code`          = node,
  `Disease Name`        = disease_name,
  `ICD-10 Chapter`      = ICD10_Chapter,
  Community             = community,
  `Role`                = role_labels[role],
  `P`                   = round(participation_coef, 3),
  `z`                   = round(within_module_z_score, 2),
  `Strength`            = round(strength, 2),
  `External strength`   = round(external_strength, 2),
  `Top 1 target`        = fifelse(!is.na(top1_ext_community),
                                  paste0("C", top1_ext_community, " (", round(top1_share_total * 100, 1), "%)"),
                                  "\u2014"),
  `Top 2 target`        = fifelse(!is.na(top2_ext_community),
                                  paste0("C", top2_ext_community, " (", round(top2_share_total * 100, 1), "%)"),
                                  "\u2014"),
  `Top 3 target`        = fifelse(!is.na(top3_ext_community),
                                  paste0("C", top3_ext_community, " (", round(top3_share_total * 100, 1), "%)"),
                                  "\u2014")
)]

# Excel
wb19 <- createWorkbook()
headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                           valign = "center", fgFill = "#4472C4", textDecoration = "bold")
addWorksheet(wb19, "Top20_Bridge")
writeDataTable(wb19, "Top20_Bridge", tab19, withFilter = TRUE)
setColWidths(wb19, "Top20_Bridge", cols = 1:ncol(tab19), widths = "auto")
addStyle(wb19, "Top20_Bridge", headerStyle,
         rows = 1, cols = 1:ncol(tab19), gridExpand = TRUE)

# 角色着色
role_fill <- c("Provincial hub" = "#FDDEDE", "Connector hub" = "#FFE8CC",
               "Kinless" = "#E8D5F5", "Connector" = "#D6EAF8", "Peripheral" = "#F2F2F2")
for (r in 1:nrow(tab19)) {
  fill_color <- role_fill[tab19[r, Role]]
  if (!is.na(fill_color)) {
    addStyle(wb19, "Top20_Bridge", createStyle(fgFill = fill_color),
             rows = r + 1, cols = 1:ncol(tab19), gridExpand = TRUE, stack = TRUE)
  }
}

xlsx19 <- file.path(OUT_DIR, "表19_桥接疾病Top20.xlsx")
saveWorkbook(wb19, xlsx19, overwrite = TRUE)
cat("   Saved:", xlsx19, "\n")

# =========================
# 4) 附表S40：桥接节点完整列表
# =========================
cat("\n==> Step 4) Creating Supplementary Table S40 (full bridge list)...\n")

# 跨社区连接节点：P≥0.30（Connector/Kinless/Hub，Guimerà–Amaral框架阈值）
bridge_all <- active[participation_coef >= P_THRESH1][order(-participation_coef)]

tabS40 <- bridge_all[, .(
  `ICD-3 Code`          = node,
  `Disease Name`        = disease_name,
  `ICD-10 Chapter`      = ICD10_Chapter,
  Community             = community,
  Role                  = role_labels[role],
  `P`                   = round(participation_coef, 4),
  `z`                   = round(within_module_z_score, 3),
  `Strength`            = round(strength, 3),
  `Within-module strength` = round(strength_within_module, 3),
  `External strength`   = round(external_strength, 3),
  `Top 1 target community` = fifelse(is.na(top1_ext_community), "\u2014",
                                     as.character(top1_ext_community)),
  `Top 1 share (%)`     = fifelse(is.na(top1_share_total), "\u2014",
                                  as.character(round(top1_share_total * 100, 1))),
  `Top 2 target community` = fifelse(is.na(top2_ext_community), "\u2014",
                                     as.character(top2_ext_community)),
  `Top 2 share (%)`     = fifelse(is.na(top2_share_total), "\u2014",
                                  as.character(round(top2_share_total * 100, 1))),
  `Top 3 target community` = fifelse(is.na(top3_ext_community), "\u2014",
                                     as.character(top3_ext_community)),
  `Top 3 share (%)`     = fifelse(is.na(top3_share_total), "\u2014",
                                  as.character(round(top3_share_total * 100, 1)))
)]

cat("   Total bridge nodes (P >= 0.30):", nrow(tabS40), "\n")

# Excel
wbS40 <- createWorkbook()
addWorksheet(wbS40, "Bridge_Nodes")
writeDataTable(wbS40, "Bridge_Nodes", tabS40, withFilter = TRUE)
setColWidths(wbS40, "Bridge_Nodes", cols = 1:ncol(tabS40), widths = "auto")
addStyle(wbS40, "Bridge_Nodes", headerStyle,
         rows = 1, cols = 1:ncol(tabS40), gridExpand = TRUE)

# 按角色着色
for (r in 1:nrow(tabS40)) {
  fill_color <- role_fill[tabS40[r, Role]]
  if (!is.na(fill_color)) {
    addStyle(wbS40, "Bridge_Nodes", createStyle(fgFill = fill_color),
             rows = r + 1, cols = 1:ncol(tabS40), gridExpand = TRUE, stack = TRUE)
  }
}

# 角色分布汇总（仅P≥0.30的节点）
role_summary <- bridge_all[, .(
  N = .N,
  `Mean P` = round(mean(participation_coef, na.rm = TRUE), 3),
  `Mean z` = round(mean(within_module_z_score, na.rm = TRUE), 2),
  `Mean strength` = round(mean(strength, na.rm = TRUE), 2)
), by = .(Role = role_labels[role])]
setorder(role_summary, -`Mean P`)

addWorksheet(wbS40, "Role_Summary")
writeDataTable(wbS40, "Role_Summary", role_summary, withFilter = TRUE)
setColWidths(wbS40, "Role_Summary", cols = 1:ncol(role_summary), widths = "auto")
addStyle(wbS40, "Role_Summary", headerStyle,
         rows = 1, cols = 1:ncol(role_summary), gridExpand = TRUE)

xlsxS40 <- file.path(OUT_DIR, "附表S40_桥接节点完整列表.xlsx")
saveWorkbook(wbS40, xlsxS40, overwrite = TRUE)
cat("   Saved:", xlsxS40, "\n")

# =========================
# 5) Console summary
# =========================
cat("\n========================================\n")
cat("DONE.\n")
cat("  图33: z-P散点图 → 图33_节点角色分类zP散点图.pdf\n")
cat("  表19: Top 20桥接疾病 → 表19_桥接疾病Top20.xlsx\n")
cat("  S40:  完整桥接列表 → 附表S40_桥接节点完整列表.xlsx\n")
cat("========================================\n")

cat("\n==> Top 5 bridge diseases:\n")
for (r in 1:min(5, nrow(top20))) {
  row <- top20[r]
  cat(sprintf("   %d. %s (%s) — Community %d, P=%.3f, z=%.2f, role=%s\n",
              r, row$node, row$disease_name, row$community,
              row$participation_coef, row$within_module_z_score, row$role))
}
