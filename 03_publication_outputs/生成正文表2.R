#!/usr/bin/env Rscript
# ==============================================================================
# Generate Table 2: Disease community composition and topological features
#
# Lancet/BMC Medicine style — strict terminology, ICD codes only (no disease names)
#
# Columns (6):
#   1. Phenotype                       — e.g. "C14 — CKD with anemia, electrolyte & gout"
#   2. No. of diseases                 — count of ICD-10 nodes in community
#   3. Most prevalent diseases         — top 3 by baseline prevalence
#   4. Hub diseases                    — top 3 by within-community node strength
#   5. Connector diseases              — top 3 by participation coefficient
#   6. Most enriched ICD-10 chapter    — chapter with highest O/E among FDR<0.05 enriched
#
# Inputs:
#   - medoid_bridge_metrics_coef0.01_alpha0.05.csv (node, community, strength_within_module,
#                                                   participation_coef, prev)
#   - TableS6_ICD10_enrichment_full_results.xlsx (Significant_Enrichment sheet)
#
# Output:
#   - Table2_community_composition.xlsx
# ==============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# =========================
# Paths
# =========================
BRIDGE_CSV  <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\community_detection\\medoid_bridge_metrics_coef0.01_alpha0.05.csv"
ENRICH_XLSX <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\supplementary_figures\\TableS6_ICD10_enrichment_full_results.xlsx"
OUT_DIR     <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\manuscript_tables"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# Phenotype labels (与 v6 plan + 其他正文图一致)
# =========================
phenotype_labels <- c(
  "C1"  = "C1 \u2014 Functional & inflammatory GI disorders",
  "C2"  = "C2 \u2014 Chronic respiratory disease with infections",
  "C3"  = "C3 \u2014 Dermatologic & venous insufficiency disorders",
  "C4"  = "C4 \u2014 Benign gynecologic disorders",
  "C5"  = "C5 \u2014 Age-related ophthalmologic disorders",
  "C6"  = "C6 \u2014 Cerebrovascular disease with neuropsychiatric sequelae",
  "C7"  = "C7 \u2014 Upper-airway & orodental inflammation",
  "C8"  = "C8 \u2014 Hepatobiliary disorders",
  "C9"  = "C9 \u2014 Advanced solid malignancies with cachexia",
  "C10" = "C10 \u2014 Benign & malignant breast disorders",
  "C11" = "C11 \u2014 Urinary stones, BPH & UTI",
  "C12" = "C12 \u2014 Thyroid dysfunction & neoplasms",
  "C13" = "C13 \u2014 Lymphoid & myeloid malignancies",
  "C14" = "C14 \u2014 CKD with anemia, electrolyte & gout",
  "C15" = "C15 \u2014 T2DM with neuropathy & joint disease",
  "C16" = "C16 \u2014 CAD, heart failure & arrhythmia"
)

# ICD-10 chapter Roman numeral -> short name
chapter_short <- c(
  "I"    = "I Infectious",
  "II"   = "II Neoplasms",
  "III"  = "III Blood",
  "IV"   = "IV Endocrine/Metabolic",
  "V"    = "V Mental",
  "VI"   = "VI Nervous",
  "VII"  = "VII Eye",
  "VIII" = "VIII Ear",
  "IX"   = "IX Circulatory",
  "X"    = "X Respiratory",
  "XI"   = "XI Digestive",
  "XII"  = "XII Skin",
  "XIII" = "XIII Musculoskeletal",
  "XIV"  = "XIV Genitourinary"
)

# =========================
# 1) 加载 bridge metrics
# =========================
cat("==> Step 1) 加载节点拓扑指标...\n")
bridge_dt <- fread(BRIDGE_CSV)
bridge_dt[, community := as.integer(community)]
bridge_dt[, prev := as.numeric(prev)]
bridge_dt[, strength_within_module := as.numeric(strength_within_module)]
bridge_dt[, participation_coef := as.numeric(participation_coef)]
bridge_dt <- bridge_dt[community >= 1 & community <= 16]
cat(sprintf("   节点数: %d, 社区数: %d\n",
            nrow(bridge_dt), uniqueN(bridge_dt$community)))

# =========================
# 2) 加载 enrichment data
# =========================
cat("\n==> Step 2) 加载 ICD 富集结果...\n")
enrich_dt <- as.data.table(read.xlsx(ENRICH_XLSX, sheet = "Significant_Enrichment"))
cat(sprintf("   显著富集行数: %d\n", nrow(enrich_dt)))
# 列名 sanitize: 应该有 Community, ICD-10.chapter, O/E, Direction
# 仅保留 Enriched (不要 Depleted)
oe_col_name <- grep("^O", names(enrich_dt), value = TRUE)[1]  # "O/E" or "O.E"
chap_col    <- grep("ICD", names(enrich_dt), value = TRUE)[1]
dir_col     <- grep("Direction", names(enrich_dt), value = TRUE)[1]

setnames(enrich_dt, c(chap_col, oe_col_name, dir_col),
         c("chapter", "OE", "direction"))
enrich_dt[, OE := as.numeric(OE)]
enrich_dt <- enrich_dt[direction == "Enriched"]
cat(sprintf("   仅 Enriched (排除 Depleted): %d 行\n", nrow(enrich_dt)))

# =========================
# 3) 每个社区取 top 3 (3 metrics × 16 communities)
# =========================
cat("\n==> Step 3) 计算各社区 top 3 metrics...\n")

# 辅助函数: 格式化 top 3 列表
fmt_top3_prev <- function(dt) {
  # 按 prev 降序, 取前 3
  top <- dt[order(-prev)][1:min(3, .N)]
  paste(sprintf("%s (%.1f%%)", top$node, 100 * top$prev), collapse = "; ")
}

fmt_top3_hub <- function(dt) {
  top <- dt[order(-strength_within_module)][1:min(3, .N)]
  paste(sprintf("%s (%.2f)", top$node, top$strength_within_module),
        collapse = "; ")
}

fmt_top3_conn <- function(dt) {
  top <- dt[order(-participation_coef)][1:min(3, .N)]
  paste(sprintf("%s (%.2f)", top$node, top$participation_coef),
        collapse = "; ")
}

# 每个社区最显著富集 chapter (max O/E)
get_top_chapter <- function(cid) {
  sub <- enrich_dt[Community == cid]
  if (nrow(sub) == 0) return("\u2014")  # em dash
  top <- sub[order(-OE)][1]
  short_name <- chapter_short[as.character(top$chapter)]
  if (is.na(short_name)) short_name <- as.character(top$chapter)
  sprintf("%s (O/E %.2f)", short_name, top$OE)
}

# =========================
# 4) 组装 Table 2
# =========================
cat("\n==> Step 4) 组装 Table 2...\n")

table2_list <- list()
for (cid in 1:16) {
  sub <- bridge_dt[community == cid]
  if (nrow(sub) == 0) next
  
  table2_list[[cid]] <- data.table(
    Phenotype = phenotype_labels[paste0("C", cid)],
    `No. of diseases` = nrow(sub),
    `Most prevalent diseases` = fmt_top3_prev(sub),
    `Hub diseases` = fmt_top3_hub(sub),
    `Connector diseases` = fmt_top3_conn(sub),
    `Most enriched ICD-10 chapter` = get_top_chapter(cid)
  )
}
table2_dt <- rbindlist(table2_list)

cat("   预览第一行:\n")
print(table2_dt[1])

# =========================
# 5) 写入 xlsx
# =========================
cat("\n==> Step 5) 写入 xlsx...\n")

out_xlsx <- file.path(OUT_DIR, "Table2_community_composition.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Table 2")

# 标题行
title_text <- "Table 2. Disease community composition and topological features"
writeData(wb, "Table 2", title_text, startRow = 1, startCol = 1)
addStyle(wb, "Table 2", createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)

# Header + 数据
writeData(wb, "Table 2", table2_dt, startRow = 2,
          headerStyle = createStyle(textDecoration = "bold",
                                    border = "TopBottom",
                                    borderStyle = "thin",
                                    halign = "center",
                                    wrapText = TRUE))

# 数据行 wrapText
n_data_rows <- nrow(table2_dt)
addStyle(wb, "Table 2", createStyle(wrapText = TRUE, valign = "top"),
         rows = 3:(2 + n_data_rows),
         cols = 1:ncol(table2_dt), gridExpand = TRUE)

# Footnote
foot_row <- 3 + n_data_rows + 1
foot_lines <- c(
  "Diseases are reported using ICD-10 three-character codes; corresponding disease names are provided in Supplementary Table S2.",
  "Most prevalent diseases: three diseases with the highest baseline period prevalence within each phenotype; values in parentheses indicate prevalence (%).",
  "Hub diseases: three diseases with the highest within-community node strength, reflecting central position within the phenotype; values in parentheses indicate within-community strength.",
  "Connector diseases: three diseases with the highest participation coefficient, reflecting cross-phenotype connectivity; values in parentheses indicate participation coefficient.",
  "Most enriched ICD-10 chapter: chapter with the highest observed-to-expected ratio (Fisher exact test, BH-adjusted FDR < 0.05). All 16 phenotypes were significantly enriched in at least one chapter; full enrichment results are in Supplementary Table S6."
)
for (i in seq_along(foot_lines)) {
  writeData(wb, "Table 2", foot_lines[i],
            startRow = foot_row + i - 1, startCol = 1)
  addStyle(wb, "Table 2", createStyle(fontSize = 9,
                                      textDecoration = "italic",
                                      wrapText = TRUE),
           rows = foot_row + i - 1, cols = 1)
}

# 列宽
setColWidths(wb, "Table 2", cols = 1, widths = 50)
setColWidths(wb, "Table 2", cols = 2, widths = 14)
setColWidths(wb, "Table 2", cols = 3, widths = 35)
setColWidths(wb, "Table 2", cols = 4, widths = 35)
setColWidths(wb, "Table 2", cols = 5, widths = 35)
setColWidths(wb, "Table 2", cols = 6, widths = 30)

# 行高 (data rows)
setRowHeights(wb, "Table 2",
              rows = 3:(2 + n_data_rows), heights = 35)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n   \u2713 已保存: %s\n", out_xlsx))

# =========================
# 6) Console 输出
# =========================
cat("\n========================================================\n")
cat("Table 2 — preview (16 rows):\n")
cat("========================================================\n")
for (i in 1:nrow(table2_dt)) {
  row <- table2_dt[i]
  cat(sprintf("\n[%s]\n", row$Phenotype))
  cat(sprintf("  Diseases: %d\n", row$`No. of diseases`))
  cat(sprintf("  Most prevalent:    %s\n", row$`Most prevalent diseases`))
  cat(sprintf("  Hub:               %s\n", row$`Hub diseases`))
  cat(sprintf("  Connector:         %s\n", row$`Connector diseases`))
  cat(sprintf("  Most enriched:     %s\n", row$`Most enriched ICD-10 chapter`))
}
cat(sprintf("\n保存路径: %s\n", out_xlsx))