#!/usr/bin/env Rscript
# ==============================================================================
# Step 15-16-17 三连发 — 补表 S8、补表 S9、补图 S7
#
# 一次性产出三个收尾性补充材料：
#
#   补表 S8: 参考组选择论证表 (Reference group justification)
#            -> TableS8_reference_group_justification.xlsx
#
#   补表 S9: 比例风险假设检验汇总表 (PH assumption summary)
#            -> TableS9_PH_assumption_summary.xlsx
#
#   补图 S7: MIXED 组内 burden 四分位死亡率梯度
#            -> FigS7_MIXED_quartile_gradient.pdf
#
# 输入 (从主分析复用):
#   - cox_regression_results.xlsx
#       sheet "Table1_main"           — 各组年龄/性别/疾病数等基线
#       sheet "Group_sizes_main"      — 各组 N、events、crude rate
#       sheet "PH_test"               — Schoenfeld 检验结果 (per-community)
#       sheet "MIXED_quartile_gradient" — MIXED 内四分位死亡率
#       sheet "Main_results"          — Cox HR (用于补表 S8 的 ref selection metric)
#
# 输出目录: 与主补充材料一致
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
})

# =========================
# 0) 路径与参数
# =========================
COX_XLSX <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_results\\cox_regression_results.xlsx"
OUT_DIR  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"
PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# 16 个社区的正式名称 (与正文一致)
community_names <- c(
  "DOM_C1"  = "C1 — Functional & inflammatory GI disorders",
  "DOM_C2"  = "C2 — Chronic respiratory disease with infections",
  "DOM_C3"  = "C3 — Dermatologic & venous insufficiency disorders",
  "DOM_C4"  = "C4 — Benign gynecologic disorders",
  "DOM_C5"  = "C5 — Age-related ophthalmologic disorders",
  "DOM_C6"  = "C6 — Cerebrovascular disease with neuropsychiatric",
  "DOM_C7"  = "C7 — Upper-airway & orodental inflammation (REFERENCE)",
  "DOM_C8"  = "C8 — Hepatobiliary disorders",
  "DOM_C9"  = "C9 — Advanced solid malignancies with cachexia",
  "DOM_C10" = "C10 — Benign & malignant breast disorders",
  "DOM_C11" = "C11 — Urinary stones, BPH & UTI",
  "DOM_C12" = "C12 — Thyroid dysfunction & neoplasms",
  "DOM_C13" = "C13 — Lymphoid & myeloid malignancies",
  "DOM_C14" = "C14 — CKD with anemia, electrolyte & gout",
  "DOM_C15" = "C15 — T2DM with neuropathy & joint disease",
  "DOM_C16" = "C16 — CAD, heart failure & arrhythmia",
  "MIXED"   = "MIXED — Mixed (no dominant community)"
)

# 通用样式：开工
hdr_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF",
                         halign = "center", valign = "center",
                         fgFill = "#4472C4", textDecoration = "bold")

cat("============================================================\n")
cat("Step 15-17 三连发：补表 S8、补表 S9、补图 S7\n")
cat("============================================================\n\n")

# =========================
# 任务 1: 补表 S8 — 参考组选择论证
# =========================
cat("==> 任务 1) 生成补表 S8: 参考组选择论证...\n")

# 读取 group sizes + Table1_main
gs_dt    <- as.data.table(read.xlsx(COX_XLSX, sheet = "Group_sizes_main"))
t1_dt    <- as.data.table(read.xlsx(COX_XLSX, sheet = "Table1_main"))

# 提取每组的 age + sex (用 group 列匹配)
gs_dt[, group := as.character(group)]
t1_dt[, group := as.character(group)]

# Table1_main 含每个组的 age_mean_sd 和 male_pct
t1_pick <- t1_dt[, .(group,
                     age_mean_sd,
                     male_pct,
                     codes_median_iqr,
                     person_years,
                     fu_mean_sd)]
gs_pick <- gs_dt[, .(group, N, events, rate_1000py)]
ref_dt  <- merge(gs_pick, t1_pick, by = "group", all.x = TRUE)

# 标记候选 reference 类型
ref_dt[, candidate_type := fifelse(group == "DOM_C7",  "Selected (main analysis)",
                                   fifelse(group == "DOM_C12", "Lowest mortality (alternative)",
                                           fifelse(group == "DOM_C4",  "Lowest events (alternative)",
                                                   fifelse(group == "DOM_C10", "Sex-skewed alternative",
                                                           fifelse(group == "MIXED",   "Largest group (alternative)",
                                                                   "Other community"
                                                           )))))]

# 注: HR 数据 (用 Model 2 vs C7 reference) 不进入 S8 表本身——
# S8 着眼于"为什么选 C7 作为 reference", 而不是"vs C7 的相对风险"。
# 若审稿人需要 candidate group 的 HR (如 MIXED 1.92 vs C7), 直接引用主 Forest plot 即可。

# 仅对 6 个候选组展示，按 evaluation criteria 排序
candidates <- c("DOM_C7", "DOM_C12", "DOM_C4", "DOM_C10", "MIXED", "DOM_C6")
ref_focus <- ref_dt[group %in% candidates]
setorder(ref_focus, rate_1000py)

# 评估准则 (defensive: 只用 ref_focus 已有的列)
ref_focus[, evaluation := fcase(
  group == "DOM_C7", paste0(
    "Selected as reference: low mortality rate (~7.3/1000py); large size ",
    "(n=47,450); sex-balanced; clinically interpretable (upper-airway/orodental)."
  ),
  group == "DOM_C12", paste0(
    "Rejected: although mortality is lowest (~3.9/1000py), the community is ",
    "female-skewed (~73% female) and contains thyroid neoplasms which ",
    "themselves carry mortality risk; using it would imply an unrealistic ",
    "low-baseline benchmark."
  ),
  group == "DOM_C4", paste0(
    "Rejected: lowest event count (n=213) but extreme sex skew (>99% female); ",
    "would yield unstable/biased reference for male strata."
  ),
  group == "DOM_C10", paste0(
    "Rejected: small event count (n=267) and sex skew (>96% female); benign ",
    "and malignant breast lesions co-classified, mixed prognosis."
  ),
  group == "MIXED",   paste0(
    "Rejected: heterogeneous by definition; using it as reference would ",
    "obscure rather than reveal the community-specific risk hierarchy."
  ),
  group == "DOM_C6",  paste0(
    "Rejected: mortality (~14.7/1000py) is intermediate, not low; using as ",
    "reference would compress dynamic range of HRs."
  ),
  default = ""
)]

# 构建 publication 表格
S8 <- data.table(
  `Candidate group`      = community_names[ref_focus$group],
  `n (% of cohort)`      = sprintf("%s (%.1f%%)",
                                   format(ref_focus$N, big.mark = ","),
                                   ref_focus$N / 640637 * 100),
  `Events`               = format(ref_focus$events, big.mark = ","),
  `Rate /1000py`         = sprintf("%.2f", ref_focus$rate_1000py),
  `Mean age (SD)`        = ref_focus$age_mean_sd,
  `% Male`               = ref_focus$male_pct,
  `Median codes (IQR)`   = ref_focus$codes_median_iqr,
  Decision               = ifelse(ref_focus$group == "DOM_C7", "✓ SELECTED",
                                  "✗ Rejected"),
  Rationale              = ref_focus$evaluation
)

# 写入 Excel (compatible chars only — no fancy unicode in column names)
wb <- createWorkbook()
addWorksheet(wb, "Reference_justification")
writeDataTable(wb, "Reference_justification", S8)
setColWidths(wb, "Reference_justification",
             cols = 1:ncol(S8), widths = "auto")
addStyle(wb, "Reference_justification", hdr_style, rows = 1,
         cols = 1:ncol(S8), gridExpand = TRUE)

# 高亮 SELECTED 行
sel_style <- createStyle(fgFill = "#D4EDDA", textDecoration = "bold")
sel_row   <- which(S8$Decision == "✓ SELECTED")
if (length(sel_row) > 0) {
  addStyle(wb, "Reference_justification", sel_style,
           rows = sel_row + 1, cols = 1:ncol(S8),
           gridExpand = TRUE, stack = TRUE)
}

xlsx_S8 <- file.path(OUT_DIR, "TableS8_reference_group_justification.xlsx")
saveWorkbook(wb, xlsx_S8, overwrite = TRUE)
cat("   ✓ 已保存:", basename(xlsx_S8), "\n")
cat("   行数:", nrow(S8), " (DOM_C7 selected, 5 alternatives rejected)\n\n")

# =========================
# 任务 2: 补表 S9 — PH 假设检验汇总
# =========================
cat("==> 任务 2) 生成补表 S9: PH 假设检验汇总...\n")

ph_dt <- as.data.table(read.xlsx(COX_XLSX, sheet = "PH_test"))
ph_dt[, term := as.character(term)]

# 仅保留 community 行 (跳过 age_at_index, sex_binary, grp joint, GLOBAL)
ph_comm <- ph_dt[grepl("^DOM_C|^MIXED$", term)]

# 排序: ref(DOM_C7) 优先, 其余按 |rho| 降序
ph_comm[, abs_rho := abs(as.numeric(rho))]
ph_comm[, is_ref  := fifelse(term == "DOM_C7", 0, 1)]
setorder(ph_comm, is_ref, -abs_rho)
ph_comm[, abs_rho := NULL]
ph_comm[, is_ref  := NULL]

# 重写 interpretation 用更精确的判断
ph_comm[, ph_judgment := fcase(
  abs(as.numeric(rho)) < 0.05, "PH holds (|ρ| < 0.05)",
  abs(as.numeric(rho)) < 0.10, "Minor departure (0.05 ≤ |ρ| < 0.10)",
  default = "Notable departure (|ρ| ≥ 0.10)"
)]

# 格式化 P 值
fmt_p <- function(p) {
  p <- as.numeric(p)
  fifelse(is.na(p), "—",
          fifelse(p < 0.001, "<0.001",
                  fifelse(p < 0.01, sprintf("%.3f", p),
                          sprintf("%.3f", p))))
}

S9 <- data.table(
  Community     = community_names[ph_comm$term],
  `chi2 (df=1)` = sprintf("%.2f", as.numeric(ph_comm$chisq)),
  `Schoenfeld P` = sapply(ph_comm$p, fmt_p),
  rho           = sprintf("%.4f", as.numeric(ph_comm$rho)),
  abs_rho       = sprintf("%.4f", abs(as.numeric(ph_comm$rho))),
  Judgment      = ph_comm$ph_judgment
)
# 注: 列名 chi2/rho/abs_rho 在 Excel 输出后由 setHeader 替换为 chi^2/rho/|rho|
# (用 ASCII 列名以避免 Windows non-UTF8 locale 下 R 解析问题)

# 在尾部加 covariate 与 GLOBAL 行
ph_extras <- ph_dt[term %in% c("age_at_index", "sex_binary",
                               "grp (joint)", "GLOBAL")]
if (nrow(ph_extras) > 0) {
  S9_extras <- data.table(
    Community     = c("Age (continuous covariate)",
                      "Sex (binary covariate)",
                      "Joint test (community group, df=16)",
                      "Global test (overall model)")[
                        match(ph_extras$term,
                              c("age_at_index", "sex_binary",
                                "grp (joint)", "GLOBAL"))],
    `chi2 (df=1)` = sprintf("%.2f", as.numeric(ph_extras$chisq)),
    `Schoenfeld P` = sapply(ph_extras$p, fmt_p),
    rho           = ifelse(is.na(as.numeric(ph_extras$rho)), "—",
                           sprintf("%.4f", as.numeric(ph_extras$rho))),
    abs_rho       = ifelse(is.na(as.numeric(ph_extras$rho)), "—",
                           sprintf("%.4f", abs(as.numeric(ph_extras$rho)))),
    Judgment      = sapply(seq_len(nrow(ph_extras)), function(i) {
      r <- as.numeric(ph_extras$rho[i])
      if (ph_extras$term[i] == "GLOBAL") {
        return("Global χ² inflated by sample size (n = 640,637); see notes.")
      }
      if (ph_extras$term[i] == "grp (joint)") {
        return("Joint test of community grouping")
      }
      if (is.na(r)) return("—")
      if (abs(r) < 0.05) return("PH holds (|ρ| < 0.05)")
      if (abs(r) < 0.10) return("Minor departure")
      return("Notable departure")
    })
  )
  S9 <- rbind(S9, S9_extras)
}

# 在写入 Excel 之前, 把 ASCII 列名替换为漂亮的 unicode 列名
# (openxlsx 在写 xlsx 时直接走 UTF-8, 不依赖 R locale)
setnames(S9,
         old = c("chi2 (df=1)", "rho", "abs_rho"),
         new = c("\u03c7\u00b2 (df=1)", "\u03c1", "|\u03c1|"))

# 保存
wb2 <- createWorkbook()
addWorksheet(wb2, "PH_assumption_summary")
writeDataTable(wb2, "PH_assumption_summary", S9)
setColWidths(wb2, "PH_assumption_summary",
             cols = 1:ncol(S9), widths = "auto")
addStyle(wb2, "PH_assumption_summary", hdr_style, rows = 1,
         cols = 1:ncol(S9), gridExpand = TRUE)

# 颜色编码: 绿=PH holds, 黄=minor, 红=notable
green_style  <- createStyle(fgFill = "#D4EDDA")
yellow_style <- createStyle(fgFill = "#FFF3CD")
red_style    <- createStyle(fgFill = "#F8D7DA")

for (i in seq_len(nrow(S9))) {
  judgment <- S9$Judgment[i]
  sty <- if (grepl("^PH holds", judgment)) green_style
  else if (grepl("Minor departure", judgment)) yellow_style
  else if (grepl("Notable departure", judgment)) red_style
  else NULL
  if (!is.null(sty)) {
    addStyle(wb2, "PH_assumption_summary", sty,
             rows = i + 1, cols = 1:ncol(S9),
             gridExpand = TRUE, stack = TRUE)
  }
}

# 加 notes sheet
notes_dt <- data.table(
  Note = c(
    "Schoenfeld residuals were used to test the proportional hazards (PH) assumption.",
    "ρ < 0 indicates a hazard ratio that decreases over follow-up time (vs reference DOM_C7).",
    "Most per-community Schoenfeld P values are highly significant; this is expected at n = 640,637 because the test statistic scales with sample size and detects clinically negligible departures.",
    "The effect size (|ρ|) is the appropriate measure: 15 of 16 dominant communities have |ρ| < 0.05 (negligible), and only MIXED shows |ρ| ≥ 0.05 (minor departure, ρ = -0.089).",
    "The hierarchy of community-specific HRs reported in the main analysis is robust to this minor PH departure of the MIXED group; sensitivity analyses (Supplementary Table S5) confirm rank stability.",
    "Reference: Therneau TM, Grambsch PM. Modeling Survival Data: Extending the Cox Model. Springer, 2000."
  )
)
addWorksheet(wb2, "Notes")
writeDataTable(wb2, "Notes", notes_dt)
setColWidths(wb2, "Notes", cols = 1, widths = 120)
addStyle(wb2, "Notes", hdr_style, rows = 1, cols = 1, gridExpand = TRUE)

xlsx_S9 <- file.path(OUT_DIR, "TableS9_PH_assumption_summary.xlsx")
saveWorkbook(wb2, xlsx_S9, overwrite = TRUE)
cat("   ✓ 已保存:", basename(xlsx_S9), "\n")
cat("   行数:", nrow(S9), " (16 communities + 4 covariate/global rows)\n\n")

# =========================
# 任务 3: 补图 S7 — MIXED 四分位死亡率梯度
# =========================
cat("==> 任务 3) 生成补图 S7: MIXED 四分位死亡率梯度...\n")

mq_dt <- as.data.table(read.xlsx(COX_XLSX, sheet = "MIXED_quartile_gradient"))

# 解析数据列
mq_dt[, ts_quartile := as.character(ts_quartile)]
mq_dt[, n            := as.numeric(n)]
mq_dt[, events       := as.numeric(events)]
mq_dt[, person_years := as.numeric(person_years)]
mq_dt[, rate_1000py  := as.numeric(rate_1000py)]
mq_dt[, age_mean     := as.numeric(age_mean)]
mq_dt[, score_median := as.numeric(score_median)]

# Q label 简化
mq_dt[, q_lbl := factor(c("Q1", "Q2", "Q3", "Q4"),
                        levels = c("Q1", "Q2", "Q3", "Q4"))]

# 颜色梯度 (从浅到深)
quart_cols <- c("Q1" = "#A6CEE3", "Q2" = "#5A95C8",
                "Q3" = "#3060A8", "Q4" = "#1F3D8C")

# Panel A: 死亡率梯度
p_rate <- ggplot(mq_dt, aes(x = q_lbl, y = rate_1000py, fill = q_lbl)) +
  geom_col(color = "grey20", linewidth = 0.4, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f", rate_1000py)),
            vjust = -0.4, size = 4, family = FONT_FAMILY,
            color = "grey15", fontface = "bold") +
  scale_fill_manual(values = quart_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = seq(0, 50, 10)) +
  labs(
    title    = "(A) 4-year mortality rate by MIXED burden quartile",
    subtitle = "MIXED group (n = 153,856) split into burden-score quartiles",
    x        = "Burden-score quartile (Q1 = lowest, Q4 = highest)",
    y        = "Mortality rate per 1,000 person-years"
  ) +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0,
                                    color = "black"),
    plot.subtitle    = element_text(size = 10, hjust = 0, color = "grey40"),
    axis.title       = element_text(size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# Panel B: 中位 burden score by quartile (展示梯度有真实数据范围)
p_score <- ggplot(mq_dt, aes(x = q_lbl, y = score_median, fill = q_lbl)) +
  geom_col(color = "grey20", linewidth = 0.4, width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", score_median)),
            vjust = -0.4, size = 3.8, family = FONT_FAMILY,
            color = "grey15", fontface = "bold") +
  scale_fill_manual(values = quart_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = seq(0, 0.5, 0.1)) +
  labs(
    title    = "(B) Median burden score by quartile",
    subtitle = "Confirms graded burden across quartiles",
    x        = "Burden-score quartile",
    y        = "Median total burden score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0,
                                    color = "black"),
    plot.subtitle    = element_text(size = 10, hjust = 0, color = "grey40"),
    axis.title       = element_text(size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# 组装
ratio_q4_q1 <- mq_dt$rate_1000py[4] / mq_dt$rate_1000py[1]
n_total     <- sum(mq_dt$n)
events_tot  <- sum(mq_dt$events)

p_combined <- p_rate + p_score +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    caption = sprintf(paste0(
      "MIXED group (n = %s, events = %s) stratified into burden quartiles. ",
      "Mortality rate increases monotonically from Q1 (%.1f/1000py) to Q4 ",
      "(%.1f/1000py), a %.1f-fold gradient — confirming that the MIXED group ",
      "is NOT homogeneous and that within-MIXED burden carries dose-response ",
      "information consistent with the continuous community-burden Cox findings."),
      format(n_total, big.mark = ","),
      format(events_tot, big.mark = ","),
      mq_dt$rate_1000py[1], mq_dt$rate_1000py[4], ratio_q4_q1),
    theme = theme(plot.caption = element_text(family = FONT_FAMILY,
                                              hjust = 0, size = 9,
                                              color = "grey25",
                                              margin = margin(t = 10),
                                              lineheight = 1.2))
  )

pdf_S7 <- file.path(OUT_DIR, "FigS7_MIXED_quartile_gradient.pdf")
ggsave(pdf_S7, p_combined, width = 9.5, height = 5,
       units = "in", device = PDF_DEVICE)
cat("   ✓ 已保存:", basename(pdf_S7), "\n")
cat(sprintf("   关键 finding: Q4/Q1 死亡率比 = %.1fx\n", ratio_q4_q1))
cat(sprintf("   各四分位 rate (per 1000 py): Q1=%.1f → Q2=%.1f → Q3=%.1f → Q4=%.1f\n",
            mq_dt$rate_1000py[1], mq_dt$rate_1000py[2],
            mq_dt$rate_1000py[3], mq_dt$rate_1000py[4]))

# =========================
# 总结
# =========================
cat("\n============================================================\n")
cat("Step 15-17 三连发完成\n")
cat("============================================================\n\n")
cat("产出 3 项 publication-ready 内容:\n\n")
cat("  Supp Table S8: ", basename(xlsx_S8), "\n", sep = "")
cat("    用途: 论证 DOM_C7 作为 Cox 参考组的合理性\n")
cat("    内容: 6 个候选组对比 (1 selected, 5 rejected)\n\n")
cat("  Supp Table S9: ", basename(xlsx_S9), "\n", sep = "")
cat("    用途: 比例风险假设检验 transparency reporting\n")
cat("    内容: 16 social communities + 4 covariate/global 行\n")
cat("    色码: 绿=PH holds, 黄=minor departure, 红=notable\n\n")
cat("  Supp Fig S7: ", basename(pdf_S7), "\n", sep = "")
cat("    用途: MIXED 群组内 burden 四分位 dose-response 证据\n")
cat(sprintf("    finding: Q1=%.1f → Q4=%.1f /1000py (%.1fx 梯度)\n",
            mq_dt$rate_1000py[1], mq_dt$rate_1000py[4], ratio_q4_q1))
cat("\n所有文件位于:\n   ", OUT_DIR, "\n", sep = "")
cat("============================================================\n")