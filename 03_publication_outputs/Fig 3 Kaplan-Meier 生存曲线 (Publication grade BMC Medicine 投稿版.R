#!/usr/bin/env Rscript
# ==============================================================================
# Fig 3: Kaplan-Meier 生存曲线 (Publication grade, BMC Medicine 投稿版)
#
# 对应 PhD 图 3-34 的 CKM 版本
#
# 设计 (按 v6 蓝本):
#   - 单 panel, 16 主导社区 + MIXED + 参考组 C7
#   - CKM 核心三层级 (paper 主线 finding): C14 红, C16 橙红, C6 橙 — 粗实线
#   - 其他高死亡 (恶性肿瘤/感染): C9, C13, C2 — 紫色系实线
#   - 低死亡: C12, C4, C10, C5, C15 — 蓝色系实线
#   - 中等死亡 (淡化, 非主线): C1, C3, C8, C11 — 浅灰细线
#   - 黑色虚线: C7 (REFERENCE)
#   - 深灰虚线: MIXED
#   - X 轴: 0–4 年; Y 轴: 累积生存率
#   - 出版级: cairo_pdf, A4 横向, Arial 字体
#
# 输入:
#   - cox_dataset_with_dominant_community.csv (Step 4.5 输出)
#       列: pid, time_years, event, dom_community_unw
#
# 输出:
#   - Fig3_KM_survival_curves.pdf
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(scales)
})

# 检查 ggplot2 版本 (linewidth aesthetic 需要 >= 3.4.0)
GGPLOT_HAS_LINEWIDTH <- utils::packageVersion("ggplot2") >= "3.4.0"
if (!GGPLOT_HAS_LINEWIDTH) {
  cat("提示: ggplot2 版本 < 3.4.0, 将使用统一 linewidth (size aesthetic 已弃用)\n")
}

# =========================
# 0) 路径与参数
# =========================
COX_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community\\cox_dataset_with_dominant_community.csv"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\main_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# 16 社区简短名 (与 Forest plot Fig 4 一致)
community_short_names <- c(
  "DOM_C1"  = "C1 GI inflam.",
  "DOM_C2"  = "C2 Resp. infect.",
  "DOM_C3"  = "C3 Skin/venous",
  "DOM_C4"  = "C4 Gyne",
  "DOM_C5"  = "C5 Ophthalm.",
  "DOM_C6"  = "C6 CVD/neuro",
  "DOM_C7"  = "C7 ENT/oral (Ref.)",
  "DOM_C8"  = "C8 Hepatobil.",
  "DOM_C9"  = "C9 Solid malig.",
  "DOM_C10" = "C10 Breast",
  "DOM_C11" = "C11 Urology",
  "DOM_C12" = "C12 Thyroid",
  "DOM_C13" = "C13 Heme malig.",
  "DOM_C14" = "C14 CKD",
  "DOM_C15" = "C15 T2DM",
  "DOM_C16" = "C16 CAD/HF",
  "MIXED"   = "MIXED"
)

# 颜色策略 — 按 paper 叙事分层配色 (不机械按死亡率排序)
# CKM 核心三层级 (paper 主线 finding): C14 > C16 > C6 — 三原色强对比, 粗线
ckm_core_colors <- c(
  "DOM_C14" = "#4A148C",  # 深紫 — 肾 (Genitourinary 经典色)
  "DOM_C16" = "#B71C1C",  # 深红 — 心 (Cardiac 直觉色)
  "DOM_C6"  = "#0D47A1"   # 深蓝 — 脑 (Nervous/Cerebrovascular)
)
# 其他高死亡率社区 (恶性肿瘤/呼吸感染): 暖色系 (避免与 CKM 紫红撞车)
other_high_colors <- c(
  "DOM_C9"  = "#E64A19",  # 橙红
  "DOM_C13" = "#F57C00",  # 橙
  "DOM_C2"  = "#D81B60"   # 洋红
)
# 低死亡率社区: 蓝色系
low_colors <- c(
  "DOM_C12" = "#1A237E",  # 深蓝
  "DOM_C4"  = "#283593",
  "DOM_C10" = "#3949AB",
  "DOM_C5"  = "#1976D2",
  "DOM_C15" = "#0277BD"
)
# 中等死亡率: 各自不同的中性色, 但饱和度低 (不抢戏 CKM 核心)
mid_colors <- c(
  "DOM_C1"  = "#8D6E63",  # 棕褐
  "DOM_C3"  = "#607D8B",  # 青灰
  "DOM_C8"  = "#455A64",  # 深灰青
  "DOM_C11" = "#A1887F"   # 棕黄
)
mid_groups <- names(mid_colors)

# 特殊
ref_color   <- "#000000"  # C7 黑色虚线
mixed_color <- "#616161"  # MIXED 深灰虚线

# =========================
# 1) 加载数据
# =========================
cat("==> Step 1) 加载 Cox 数据集...\n")

if (!file.exists(COX_CSV)) {
  stop("Input not found: ", COX_CSV, "\nRun Step 4.5 first.", call. = FALSE)
}

dt <- fread(COX_CSV,
            select = c("time_years", "event", "dom_community_unw"))
cat(sprintf("   原始: %d 行\n", nrow(dt)))

# 排除 NONE 组
dt <- dt[dom_community_unw != "NONE"]
cat(sprintf("   排除 NONE 后: %d 行\n", nrow(dt)))

dt[, grp := factor(dom_community_unw, levels = names(community_short_names))]
cat(sprintf("   各组人数:\n"))
gs <- dt[, .(n = .N, events = sum(event),
             rate_1000py = sum(event) / sum(time_years) * 1000), by = grp]
setorder(gs, -rate_1000py)
for (i in seq_len(nrow(gs))) {
  cat(sprintf("     %-13s n=%-7s events=%-6s rate=%.2f /1000py\n",
              gs$grp[i], format(gs$n[i], big.mark=","),
              format(gs$events[i], big.mark=","), gs$rate_1000py[i]))
}

# =========================
# 2) 拟合 KM 曲线
# =========================
cat("\n==> Step 2) 拟合 Kaplan-Meier 曲线...\n")

sf <- survfit(Surv(time_years, event) ~ grp, data = dt)
cat(sprintf("   完成, %d strata\n", length(sf$strata)))

# 提取曲线数据 — 包含起点 (t=0, surv=1)
sdata_list <- list()
strata_names <- names(sf$strata)
strata_lens  <- sf$strata
cum_idx <- 0
for (s_idx in seq_along(strata_names)) {
  n_pts <- strata_lens[s_idx]
  rng <- (cum_idx + 1):(cum_idx + n_pts)
  cum_idx <- cum_idx + n_pts
  
  group_name <- gsub("^grp=", "", strata_names[s_idx])
  
  # 加 t=0 起点
  sdata_list[[s_idx]] <- data.table(
    time  = c(0, sf$time[rng]),
    surv  = c(1, sf$surv[rng]),
    n_risk = c(sum(dt$grp == group_name), sf$n.risk[rng]),
    group = group_name
  )
}
sdata <- rbindlist(sdata_list)

# 4 年截断 (有的曲线可能有 > 4 年的尾)
sdata <- sdata[time <= 4.05]

cat(sprintf("   KM 数据点数: %d\n", nrow(sdata)))

# =========================
# 3) 配置颜色 + 线型映射
# =========================
all_colors <- c(ckm_core_colors, other_high_colors, low_colors,
                mid_colors,
                "DOM_C7" = ref_color,
                "MIXED" = mixed_color)

# 线型映射: REF + MIXED 是虚线, 其他实线
all_linetypes <- setNames(rep("solid", length(all_colors)), names(all_colors))
all_linetypes["DOM_C7"] <- "dashed"
all_linetypes["MIXED"]  <- "dashed"

# =========================
# 4) 主 KM panel
# =========================
cat("\n==> Step 3) 绘制主 KM panel...\n")

# 重新排序 group factor (上方暖色 → 中间灰 → 下方冷色 → REF/MIXED)
group_order <- c(names(ckm_core_colors), names(other_high_colors),
                 names(low_colors), mid_groups, "DOM_C7", "MIXED")
sdata[, group := factor(group, levels = group_order)]

# 计算 y 轴下限 (避免顶部空太多)
y_min <- max(0.80, min(sdata$surv, na.rm = TRUE) - 0.02)

# 颜色与线型
color_vec    <- all_colors[group_order]
linetype_vec <- all_linetypes[group_order]

# Legend label: 仅 ID (C1, C2, ..., MIXED) — 不带名称和死亡率, 保持图例简洁
label_vec <- sapply(group_order, function(g) {
  if (g == "MIXED") "MIXED"
  else if (g == "DOM_C7") "C7 (Ref.)"
  else sub("DOM_", "", g)   # DOM_C9 -> C9
})

p_km <- ggplot(sdata, aes(x = time, y = surv, color = group,
                          linetype = group, group = group)) +
  # 中等死亡率社区 (mid_groups): 浅灰细线 — 先画在底层
  geom_step(data = sdata[group %in% mid_groups],
            linewidth = 0.40) +
  # 低死亡率社区 + REF + MIXED: 中等粗细
  geom_step(data = sdata[group %in% c(names(low_colors),
                                      "DOM_C7", "MIXED")],
            linewidth = 0.75) +
  # 其他高死亡率: 粗
  geom_step(data = sdata[group %in% names(other_high_colors)],
            linewidth = 0.85) +
  # CKM 核心三层级 (paper 主线 finding): 最粗, 画在最上层
  geom_step(data = sdata[group %in% names(ckm_core_colors)],
            linewidth = 1.05) +
  scale_color_manual(values = color_vec, labels = label_vec, name = NULL) +
  scale_linetype_manual(values = linetype_vec, labels = label_vec, name = NULL) +
  scale_x_continuous(breaks = 0:4, limits = c(0, 4),
                     expand = expansion(mult = c(0.005, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(y_min, 1.001),
                     expand = c(0, 0)) +
  labs(x = "Follow-up time (years)",
       y = "Cumulative survival probability") +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    axis.title       = element_text(size = 12, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.text        = element_text(size = 10, color = "black"),
    legend.key.height  = unit(0.5, "cm"),
    legend.key.width   = unit(1.2, "cm"),
    plot.margin        = margin(10, 10, 5, 10)
  ) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE,
                              override.aes = list(linewidth = 1.0)),
         linetype = guide_legend(ncol = 1, byrow = TRUE))

# =========================
# 5) 保存
# =========================
cat("==> Step 4) 保存图形...\n")

pdf_file <- file.path(OUT_DIR, "Fig3_KM_survival_curves.pdf")
ggsave(pdf_file, p_km,
       width = 11, height = 7, units = "in",
       device = PDF_DEVICE, dpi = 300)

cat(sprintf("\n   ✓ 已保存: %s\n", basename(pdf_file)))

# =========================
# 8) Console 报告
# =========================
cat("\n============================================================\n")
cat("Fig 3: KM 生存曲线 — 完成\n")
cat("============================================================\n\n")
cat("分组配色策略 (按 paper 叙事):\n")
cat("  CKM 核心三层级 (粗线): C14 (renal) > C16 (cardiac) > C6 (cerebrovascular)\n")
cat("  其他高死亡率: C9, C13, C2 (紫色系)\n")
cat("  低死亡率: C12, C4, C10, C5, C15 (蓝色系)\n")
cat("  中等死亡率 (淡化非主线): C1, C3, C8, C11 (浅灰细线)\n")
cat("  REFERENCE: C7 (黑色虚线)\n")
cat("  MIXED:     深灰虚线\n\n")
cat("各组 4 年死亡率排名:\n")
for (i in seq_len(nrow(gs))) {
  cat(sprintf("   %2d. %-15s rate = %5.2f /1000py  (n = %s, events = %s)\n",
              i, community_short_names[as.character(gs$grp[i])],
              gs$rate_1000py[i],
              format(gs$n[i], big.mark = ","),
              format(gs$events[i], big.mark = ",")))
}
cat("\n输出文件:\n  ", pdf_file, "\n", sep = "")
cat("============================================================\n")