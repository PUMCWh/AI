#!/usr/bin/env Rscript
# ==============================================================================
# Master script: Generate Supplementary Tables S2 - S9 in publication format
#
# Produces 7 xlsx files in OUT_DIR, all Lancet/BMC Medicine appendix style:
#   TableS2_full_MMC_ranking.xlsx
#   TableS3_top_bridge_diseases.xlsx
#   TableS5_sex_stratified_HR.xlsx
#   TableS6_age_stratified_HR.xlsx
#   TableS7_ICD10_enrichment_full.xlsx
#   TableS8_sensitivity_analyses.xlsx
#   TableS9_PH_assumption_test.xlsx
#
# Inputs (auto-detected):
#   - NETWORK_PCOR_<timestamp>/MMC_coef0.01_alpha0.05.csv
#   - community_detection/medoid_bridge_metrics_coef0.01_alpha0.05.csv
#   - subgroup_interaction_analysis.xlsx (Sex_interaction + Age_interaction)
#   - cox_results/cox_regression_results.xlsx (Sensitivity_all + PH_test)
#   - supplementary_figures/TableS6_ICD10_enrichment_full_results.xlsx
#
# Output: C:\Users\HP\Desktop\paper\CKM_FINAL\manuscript_tables\
# ==============================================================================

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# =========================
# Paths (auto-detect when possible)
# =========================
PAPER_ROOT <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL"
OUT_DIR    <- file.path(PAPER_ROOT, "manuscript_tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Auto-find PCOR network folder (latest timestamp)
pcor_dirs <- list.files(PAPER_ROOT, pattern = "^NETWORK_PCOR_",
                        full.names = TRUE)
if (length(pcor_dirs) == 0) {
  stop("No NETWORK_PCOR_* folder found in ", PAPER_ROOT)
}
PCOR_DIR <- sort(pcor_dirs, decreasing = TRUE)[1]   # latest by timestamp

MMC_CSV     <- file.path(PCOR_DIR,
                         "MMC_coef0.01_alpha0.05.csv")
BRIDGE_CSV  <- file.path(PAPER_ROOT, "community_detection",
                         "medoid_bridge_metrics_coef0.01_alpha0.05.csv")
SUBGRP_XLSX <- file.path(PAPER_ROOT, "cox_analysis", "cox_results",
                         "subgroup_interaction_analysis.xlsx")
COX_XLSX    <- file.path(PAPER_ROOT, "cox_analysis", "cox_results",
                         "cox_regression_results.xlsx")
ENRICH_XLSX <- file.path(PAPER_ROOT, "supplementary_figures",
                         "TableS6_ICD10_enrichment_full_results.xlsx")

# Fallback paths if primary location not found
if (!file.exists(SUBGRP_XLSX)) {
  for (alt in c(file.path(PAPER_ROOT, "cox_analysis",
                          "subgroup_interaction_analysis.xlsx"),
                file.path(PAPER_ROOT,
                          "subgroup_interaction_analysis.xlsx"))) {
    if (file.exists(alt)) { SUBGRP_XLSX <- alt; break }
  }
}

cat("Input paths:\n")
cat(sprintf("  MMC:        %s  (%s)\n", basename(MMC_CSV),
            ifelse(file.exists(MMC_CSV), "OK", "MISSING")))
cat(sprintf("  Bridge:     %s  (%s)\n", basename(BRIDGE_CSV),
            ifelse(file.exists(BRIDGE_CSV), "OK", "MISSING")))
cat(sprintf("  Subgroup:   %s  (%s)\n", basename(SUBGRP_XLSX),
            ifelse(file.exists(SUBGRP_XLSX), "OK", "MISSING")))
cat(sprintf("  Cox xlsx:   %s  (%s)\n", basename(COX_XLSX),
            ifelse(file.exists(COX_XLSX), "OK", "MISSING")))
cat(sprintf("  Enrichment: %s  (%s)\n", basename(ENRICH_XLSX),
            ifelse(file.exists(ENRICH_XLSX), "OK", "MISSING")))

# =========================
# Phenotype labels (consistent with main tables)
# =========================
phenotype_labels <- c(
  "DOM_C1"  = "C1 \u2014 Functional & inflammatory GI disorders",
  "DOM_C2"  = "C2 \u2014 Chronic respiratory disease with infections",
  "DOM_C3"  = "C3 \u2014 Dermatologic & venous insufficiency disorders",
  "DOM_C4"  = "C4 \u2014 Benign gynecologic disorders",
  "DOM_C5"  = "C5 \u2014 Age-related ophthalmologic disorders",
  "DOM_C6"  = "C6 \u2014 Cerebrovascular disease with neuropsychiatric sequelae",
  "DOM_C7"  = "C7 \u2014 Upper-airway & orodental inflammation",
  "DOM_C8"  = "C8 \u2014 Hepatobiliary disorders",
  "DOM_C9"  = "C9 \u2014 Advanced solid malignancies with cachexia",
  "DOM_C10" = "C10 \u2014 Benign & malignant breast disorders",
  "DOM_C11" = "C11 \u2014 Urinary stones, BPH & UTI",
  "DOM_C12" = "C12 \u2014 Thyroid dysfunction & neoplasms",
  "DOM_C13" = "C13 \u2014 Lymphoid & myeloid malignancies",
  "DOM_C14" = "C14 \u2014 CKD with anemia, electrolyte & gout",
  "DOM_C15" = "C15 \u2014 T2DM with neuropathy & joint disease",
  "DOM_C16" = "C16 \u2014 CAD, heart failure & arrhythmia",
  "MIXED"   = "MIXED"
)
display_order <- c("DOM_C7",
                   paste0("DOM_C", c(1:6, 8:16)),
                   "MIXED")

# Roman numeral chapter labels
roman_chapters <- c(
  "1" = "I", "2" = "II", "3" = "III", "4" = "IV", "5" = "V",
  "6" = "VI", "7" = "VII", "8" = "VIII", "9" = "IX", "10" = "X",
  "11" = "XI", "12" = "XII", "13" = "XIII", "14" = "XIV"
)

# =========================
# Helper: write publication-style xlsx
# =========================
write_pub_xlsx <- function(dt, file_path, sheet_name,
                           title_text, footnotes,
                           col_widths = NULL,
                           left_cols = 1) {
  wb <- createWorkbook()
  addWorksheet(wb, sheet_name)
  
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1)
  addStyle(wb, sheet_name,
           createStyle(textDecoration = "bold", fontSize = 12),
           rows = 1, cols = 1)
  
  writeData(wb, sheet_name, dt, startRow = 2,
            headerStyle = createStyle(textDecoration = "bold",
                                      border = "TopBottom",
                                      borderStyle = "thin",
                                      halign = "center",
                                      wrapText = TRUE))
  
  n_rows <- nrow(dt)
  n_cols <- ncol(dt)
  if (n_rows > 0) {
    # Left-align specified left cols
    addStyle(wb, sheet_name,
             createStyle(halign = "left", valign = "top"),
             rows = 3:(2 + n_rows),
             cols = seq_len(left_cols),
             gridExpand = TRUE)
    # Center-align rest
    if (n_cols > left_cols) {
      addStyle(wb, sheet_name,
               createStyle(halign = "center", valign = "top"),
               rows = 3:(2 + n_rows),
               cols = (left_cols + 1):n_cols,
               gridExpand = TRUE)
    }
  }
  
  # Footnotes
  foot_row <- 3 + n_rows + 1
  for (i in seq_along(footnotes)) {
    writeData(wb, sheet_name, footnotes[i],
              startRow = foot_row + i - 1, startCol = 1)
    addStyle(wb, sheet_name,
             createStyle(fontSize = 9, textDecoration = "italic",
                         wrapText = TRUE),
             rows = foot_row + i - 1, cols = 1)
  }
  
  # Column widths
  if (!is.null(col_widths)) {
    for (i in seq_along(col_widths)) {
      setColWidths(wb, sheet_name, cols = i, widths = col_widths[i])
    }
  }
  
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

# Helper: format HR (95% CI) — accepts already-formatted strings or numeric tuples
fmt_hr <- function(hr, lo, hi, digits = 2) {
  ifelse(is.na(hr) | is.na(lo) | is.na(hi),
         "\u2014",
         sprintf(paste0("%.", digits, "f (%.", digits, "f\u2013%.", digits, "f)"),
                 hr, lo, hi))
}

# Helper: format P value in Lancet/NEJM journal style
# Rules:
#   NA / em-dash    -> em-dash
#   already-formatted strings like "<0.001", "0.045" -> pass through (preserve)
#   p >= 0.10       -> 2 decimals (e.g. "0.12")
#   0.001 <= p < 0.10 -> 3 decimals (e.g. "0.043")
#   p < 0.001       -> scientific with unicode superscripts (e.g. "9.3 \u00d7 10\u207b\u00b9\u2077")
fmt_p <- function(p) {
  out <- character(length(p))
  for (i in seq_along(p)) {
    pi <- p[i]
    
    # Already-formatted string handling
    if (is.character(pi)) {
      pi_chr <- pi
      # NA-like values
      if (is.na(pi_chr) || pi_chr == "" || pi_chr == "NA" ||
          pi_chr == "\u2014" || pi_chr == "-") {
        out[i] <- "\u2014"
        next
      }
      # Pre-formatted with < or > prefix: pass through unchanged
      if (grepl("^[<>]", pi_chr)) {
        out[i] <- pi_chr
        next
      }
      # Otherwise try numeric conversion
      pi_num <- suppressWarnings(as.numeric(pi_chr))
    } else {
      pi_num <- pi
    }
    
    if (is.na(pi_num)) {
      out[i] <- "\u2014"
    } else if (pi_num >= 0.10) {
      out[i] <- sprintf("%.2f", pi_num)
    } else if (pi_num >= 0.001) {
      out[i] <- sprintf("%.3f", pi_num)
    } else if (pi_num > 0) {
      # Pretty scientific notation: "9.3 \u00d7 10\u207b\u00b9\u2077"
      s <- formatC(pi_num, format = "e", digits = 1)
      m <- regmatches(s, regexec("([0-9.]+)e([+-]?[0-9]+)", s))[[1]]
      if (length(m) >= 3) {
        base_str <- m[2]
        exp_int  <- as.integer(m[3])
        sup_digits <- c("\u2070", "\u00b9", "\u00b2", "\u00b3", "\u2074",
                        "\u2075", "\u2076", "\u2077", "\u2078", "\u2079")
        sign_str <- if (exp_int < 0) "\u207b" else ""
        digits_str <- ""
        for (d in strsplit(as.character(abs(exp_int)), "")[[1]]) {
          digits_str <- paste0(digits_str, sup_digits[as.integer(d) + 1])
        }
        out[i] <- paste0(base_str, " \u00d7 10", sign_str, digits_str)
      } else {
        out[i] <- sprintf("%.1e", pi_num)   # fallback
      }
    } else {
      # p == 0 (machine underflow)
      out[i] <- "<1 \u00d7 10\u207b\u00b3\u2070\u2070"
    }
  }
  out
}

# Helper: build standardized abbreviations footnote (alphabetical, semicolon-separated)
build_abbrev_footnote <- function(abbrev_map) {
  # abbrev_map: named character vector; names = abbreviations, values = full forms
  sorted_keys <- sort(names(abbrev_map))
  parts <- vapply(sorted_keys,
                  function(k) sprintf("%s = %s", k, abbrev_map[[k]]),
                  character(1))
  paste0("Abbreviations: ", paste(parts, collapse = "; "), ".")
}

# =========================
# Table S2: Top 50 diseases by MMC + mean partial correlation
# =========================
cat("\n========================================================\n")
cat("Table S2: Top 50 diseases by MMC (with mean partial correlation)\n")
cat("========================================================\n")

# Edges file path (used to compute mean partial correlation per node)
EDGES_CSV <- file.path(PAPER_ROOT, "community_detection",
                       "medoid_edges_with_community_coef0.01_alpha0.05.csv")

if (!file.exists(MMC_CSV)) {
  warning("MMC_CSV not found, skipping Table S2")
} else if (!file.exists(EDGES_CSV)) {
  warning("EDGES_CSV not found, skipping Table S2")
} else {
  mmc <- fread(MMC_CSV)
  mmc[, ICD10_Chapter := as.character(ICD10_Chapter)]
  
  # Compute mean partial correlation per disease from edges
  cat("   Computing mean partial correlation per disease...\n")
  edges <- fread(EDGES_CSV)
  # Each disease appears as 'from' or 'to' in undirected edges; combine
  edges_long <- rbind(
    edges[, .(node = from, weight)],
    edges[, .(node = to,   weight)]
  )
  mean_phi <- edges_long[, .(mean_phi = mean(weight, na.rm = TRUE),
                             n_edges = .N), by = node]
  cat(sprintf("   Computed mean phi for %d nodes\n", nrow(mean_phi)))
  
  # Merge MMC with mean phi
  setorder(mmc, -MMC)
  mmc[, Rank := seq_len(.N)]
  mmc <- merge(mmc, mean_phi, by.x = "DiseaseCode", by.y = "node",
               all.x = TRUE)
  setorder(mmc, Rank)
  
  # Top 50
  TOP_N_S2 <- 50
  mmc_top <- mmc[Rank <= TOP_N_S2]
  mmc_top[, Chapter_Roman := roman_chapters[ICD10_Chapter]]
  mmc_top[is.na(Chapter_Roman), Chapter_Roman := "Other"]
  mmc_top[, MMC_str  := sprintf("%.4f", MMC)]
  mmc_top[, Prev_str := sprintf("%.3f%%", Prevalence * 100)]
  mmc_top[, mean_phi_str := ifelse(is.na(mean_phi), "\u2014",
                                   sprintf("%.4f", mean_phi))]
  mmc_top[, n_edges_str := ifelse(is.na(n_edges), "\u2014",
                                  format(as.integer(n_edges),
                                         big.mark = ","))]
  
  s2_dt <- mmc_top[, .(
    Rank,
    `ICD-10 code` = DiseaseCode,
    Chapter = Chapter_Roman,
    `Frequency (n)` = format(Frequency, big.mark = ","),
    Prevalence = Prev_str,
    MMC = MMC_str,
    `No. of edges` = n_edges_str,
    mean_phi_col = mean_phi_str
  )]
  # Set unicode column name AFTER table built (avoid backtick + \u parser issue)
  setnames(s2_dt, "mean_phi_col", "Mean \u03c6")
  
  cat(sprintf("   Rows: %d (top %d)\n", nrow(s2_dt), TOP_N_S2))
  
  write_pub_xlsx(
    dt = s2_dt,
    file_path = file.path(OUT_DIR, "TableS2_top50_MMC_ranking.xlsx"),
    sheet_name = "Table S2",
    title_text = "Table S2. Top 50 diseases ranked by multimorbidity centrality",
    footnotes = c(
      paste0("Diseases are ranked in descending order of multimorbidity centrality (MMC), ",
             "a weighted composite measure of within-network connectivity that quantifies each disease's structural role in the comorbidity network. ",
             "Higher MMC values indicate greater integration into the network."),
      paste0("Mean \u03c6 is the arithmetic mean of partial correlation coefficients across all significant edges incident to that disease ",
             "(thresholds: \u03c6 \u2265 0.01, Benjamini-Hochberg FDR < 0.05). ",
             "The partial correlation network was constructed from 1.28 million disease-instance records across 640,637 baseline-period adults."),
      paste0("The full ranking of all 613 diseases is available as a supplementary data file."),
      build_abbrev_footnote(c(
        "FDR" = "false discovery rate",
        "ICD-10" = "International Classification of Diseases, 10th revision",
        "MMC" = "multimorbidity centrality",
        "\u03c6" = "partial correlation coefficient (network edge weight)"
      ))
    ),
    col_widths = c(8, 14, 10, 16, 14, 12, 14, 14),
    left_cols = 2
  )
  cat(sprintf("   \u2713 Saved: TableS2_top50_MMC_ranking.xlsx\n"))
}

# =========================
# Table S3: Top 100 bridge diseases (by participation coefficient)
# =========================
cat("\n========================================================\n")
cat("Table S3: Top 100 bridge diseases\n")
cat("========================================================\n")

if (!file.exists(BRIDGE_CSV)) {
  warning("BRIDGE_CSV not found, skipping Table S3")
} else {
  brg <- fread(BRIDGE_CSV)
  setorder(brg, -participation_coef)
  TOP_N_S3 <- 100
  brg_top <- brg[seq_len(min(TOP_N_S3, .N))]
  brg_top[, Rank := seq_len(.N)]
  brg_top[, ICD10_Chapter := as.character(ICD10_Chapter)]
  brg_top[, Chapter_Roman := roman_chapters[ICD10_Chapter]]
  brg_top[is.na(Chapter_Roman), Chapter_Roman := "Other"]
  brg_top[, PC_str := sprintf("%.3f", participation_coef)]
  brg_top[, Z_str  := sprintf("%.3f", within_module_z_score)]
  brg_top[, Strength_str := sprintf("%.2f", strength)]
  
  s3_dt <- brg_top[, .(
    Rank,
    `ICD-10 code` = node,
    Chapter = Chapter_Roman,
    `Primary community` = paste0("C", community),
    `Participation coefficient` = PC_str,
    `Within-module z-score` = Z_str,
    `Total strength` = Strength_str,
    Role = role
  )]
  
  cat(sprintf("   Rows: %d (top %d)\n", nrow(s3_dt), TOP_N_S3))
  
  write_pub_xlsx(
    dt = s3_dt,
    file_path = file.path(OUT_DIR, "TableS3_top100_bridge_diseases.xlsx"),
    sheet_name = "Table S3",
    title_text = "Table S3. Top 100 bridge diseases ranked by participation coefficient",
    footnotes = c(
      paste0("Bridge diseases connect multiple communities and are essential to the functional integration of the multimorbidity network. ",
             "The participation coefficient (PC, Guimer\u00e0 and Amaral, 2005) measures the diversity of a node's connections across communities; ",
             "PC ranges from 0 (all connections within one community) to 1 (connections evenly distributed across all communities)."),
      paste0("The within-module z-score quantifies a node's connectivity within its primary community ",
             "relative to other nodes in the same community. ",
             "Combined with the participation coefficient, it places each disease into one of seven Guimer\u00e0\u2013Amaral roles."),
      paste0("Roles in this dataset: R7 = connector hub (z \u2265 2.5, PC > 0.62); R6 = connector non-hub (z < 2.5, 0.30 < PC \u2264 0.75); ",
             "R3 = peripheral non-hub (z < 2.5, 0.05 < PC \u2264 0.30); R2 = ultra-peripheral (z < 2.5, PC \u2264 0.05). ",
             "Diseases are ranked by PC in descending order."),
      build_abbrev_footnote(c(
        "ICD-10" = "International Classification of Diseases, 10th revision",
        "PC" = "participation coefficient",
        "z-score" = "within-module z-score of total edge weight"
      ))
    ),
    col_widths = c(8, 14, 10, 18, 24, 22, 14, 18),
    left_cols = 2
  )
  cat(sprintf("   \u2713 Saved: TableS3_top100_bridge_diseases.xlsx\n"))
}

# =========================
# Table S5: Sex-stratified HR
# =========================
cat("\n========================================================\n")
cat("Table S5: Sex-stratified HR (16 communities x 2 sex)\n")
cat("========================================================\n")

if (!file.exists(SUBGRP_XLSX)) {
  warning("SUBGRP_XLSX not found, skipping Table S5")
} else {
  # Sheet has metadata in row 1, headers in row 2; skip to row 2 as header
  s5_raw <- as.data.table(read.xlsx(SUBGRP_XLSX,
                                    sheet = "Sex_interaction",
                                    startRow = 2))
  cat(sprintf("   Rows: %d  Cols: %d\n", nrow(s5_raw), ncol(s5_raw)))
  cat(sprintf("   Cols: %s\n", paste(names(s5_raw), collapse = ", ")))
  
  # Group column contains full phenotype label (e.g. "C7 — Upper airway & orodental inflammation")
  # NOT "DOM_C7". So extract community ID first, then map to standardized label.
  s5_raw[, comm_id := regmatches(Group,
                                 regexpr("^(C\\d+|MIXED)", Group))]
  # Convert to DOM_C* key (or keep MIXED)
  s5_raw[, dom_key := fifelse(comm_id == "MIXED", "MIXED",
                              paste0("DOM_", comm_id))]
  s5_raw[, Phenotype := phenotype_labels[dom_key]]
  s5_raw[is.na(Phenotype), Phenotype := Group]   # fallback to original
  s5_raw[, term_factor := factor(dom_key, levels = display_order)]
  setorder(s5_raw, term_factor)
  cat(sprintf("   Communities found: %s\n",
              paste(s5_raw$comm_id, collapse = ", ")))
  
  # Build display table
  # Identify columns with format-friendly names (read.xlsx sanitizes)
  male_n_col      <- grep("Male.N$|Male\\.N$|^Male\\.N", names(s5_raw),
                          value = TRUE)[1]
  male_e_col      <- grep("Male.Events", names(s5_raw), value = TRUE)[1]
  male_hr_col     <- grep("Male.HR", names(s5_raw), value = TRUE)[1]
  male_p_col      <- grep("^Male.P$|Male\\.P$", names(s5_raw), value = TRUE)[1]
  female_n_col    <- grep("Female.N", names(s5_raw), value = TRUE)[1]
  female_e_col    <- grep("Female.Events", names(s5_raw), value = TRUE)[1]
  female_hr_col   <- grep("Female.HR", names(s5_raw), value = TRUE)[1]
  female_p_col    <- grep("^Female.P$|Female\\.P$", names(s5_raw),
                          value = TRUE)[1]
  pinter_col      <- grep("interaction", names(s5_raw), value = TRUE,
                          ignore.case = TRUE)[1]
  
  cat(sprintf("   Male HR col: %s, Female HR col: %s, P_int: %s\n",
              male_hr_col, female_hr_col, pinter_col))
  
  # Helper for safe integer formatting (handles NA)
  fmt_int <- function(x) {
    ifelse(is.na(x), "\u2014",
           format(as.integer(x), big.mark = ","))
  }
  
  s5_dt <- s5_raw[, .(
    Phenotype = Phenotype,
    `Male n` = fmt_int(get(male_n_col)),
    `Male events` = fmt_int(get(male_e_col)),
    `Male HR (95% CI)` = ifelse(is.na(get(male_hr_col)) | get(male_hr_col) == "",
                                "\u2014", get(male_hr_col)),
    `Female n` = fmt_int(get(female_n_col)),
    `Female events` = fmt_int(get(female_e_col)),
    `Female HR (95% CI)` = ifelse(is.na(get(female_hr_col)) | get(female_hr_col) == "",
                                  "\u2014", get(female_hr_col)),
    `P-interaction` = fmt_p(get(pinter_col))
  )]
  s5_dt[Phenotype == phenotype_labels["DOM_C7"],
        c("Male HR (95% CI)", "Female HR (95% CI)") := "1.00 (reference)"]
  
  write_pub_xlsx(
    dt = s5_dt,
    file_path = file.path(OUT_DIR, "TableS5_sex_stratified_HR.xlsx"),
    sheet_name = "Table S5",
    title_text = "Table S5. Sex-stratified hazard ratios for 4-year all-cause mortality by disease community phenotype",
    footnotes = c(
      paste0("Data are HR (95% CI) versus C7 (Upper-airway & orodental inflammation, reference community). ",
             "Subgroup Cox proportional hazards models adjusted for age."),
      paste0("Global likelihood-ratio test for sex \u00d7 community interaction: \u03c7\u00b2 = 113.58, df = 16, ",
             "P = 9.3 \u00d7 10\u207b\u00b9\u2077."),
      paste0("MIXED denotes patients without a single dominant community phenotype ",
             "(no community share \u2265 0.40 at baseline)."),
      build_abbrev_footnote(c(
        "CI" = "confidence interval",
        "HR" = "hazard ratio",
        "LRT" = "likelihood-ratio test",
        "MIXED" = "patients without a single dominant community phenotype"
      ))
    ),
    col_widths = c(50, 12, 12, 22, 12, 12, 22, 16),
    left_cols = 1
  )
  cat(sprintf("   \u2713 Saved: TableS5_sex_stratified_HR.xlsx\n"))
}

# =========================
# Table S6: Age-stratified HR
# =========================
cat("\n========================================================\n")
cat("Table S6: Age-stratified HR (16 communities x 3 age groups)\n")
cat("========================================================\n")

if (!file.exists(SUBGRP_XLSX)) {
  warning("SUBGRP_XLSX not found, skipping Table S6")
} else {
  s6_raw <- as.data.table(read.xlsx(SUBGRP_XLSX,
                                    sheet = "Age_interaction",
                                    startRow = 2))
  cat(sprintf("   Rows: %d  Cols: %d\n", nrow(s6_raw), ncol(s6_raw)))
  cat(sprintf("   Cols: %s\n", paste(names(s6_raw), collapse = ", ")))
  
  # Group column contains full phenotype label (same as S5). Extract community ID.
  s6_raw[, comm_id := regmatches(Group,
                                 regexpr("^(C\\d+|MIXED)", Group))]
  s6_raw[, dom_key := fifelse(comm_id == "MIXED", "MIXED",
                              paste0("DOM_", comm_id))]
  s6_raw[, Phenotype := phenotype_labels[dom_key]]
  s6_raw[is.na(Phenotype), Phenotype := Group]
  s6_raw[, term_factor := factor(dom_key, levels = display_order)]
  setorder(s6_raw, term_factor)
  
  # Try multiple naming patterns for age groups
  young_n_col   <- grep("60.N|^X..60.N|<60.N",   names(s6_raw), value = TRUE)[1]
  young_e_col   <- grep("60.Events|<60.Events",  names(s6_raw), value = TRUE)[1]
  young_hr_col  <- grep("60.HR|<60.HR",          names(s6_raw), value = TRUE)[1]
  mid_n_col     <- grep("60.75.N",               names(s6_raw), value = TRUE)[1]
  mid_e_col     <- grep("60.75.Events",          names(s6_raw), value = TRUE)[1]
  mid_hr_col    <- grep("60.75.HR",              names(s6_raw), value = TRUE)[1]
  old_n_col     <- grep("75.N$|>=75.N",          names(s6_raw), value = TRUE)[1]
  old_e_col     <- grep("75.Events|>=75.Events", names(s6_raw), value = TRUE)[1]
  old_hr_col    <- grep("75.HR|>=75.HR",         names(s6_raw), value = TRUE)[1]
  pinter_col    <- grep("interaction",           names(s6_raw),
                        value = TRUE, ignore.case = TRUE)[1]
  
  cat(sprintf("   <60 HR col: %s, 60-75 HR: %s, >=75 HR: %s\n",
              young_hr_col, mid_hr_col, old_hr_col))
  
  # Helper for safe integer formatting + n/events string
  fmt_int_safe <- function(x) {
    ifelse(is.na(x), "\u2014",
           format(as.integer(x), big.mark = ","))
  }
  fmt_ne <- function(n_col, e_col) {
    sprintf("%s / %s", fmt_int_safe(n_col), fmt_int_safe(e_col))
  }
  fmt_hr_safe <- function(x) {
    ifelse(is.na(x) | x == "", "\u2014", x)
  }
  
  s6_dt <- s6_raw[, .(
    Phenotype = Phenotype,
    `<60: n / events` = fmt_ne(get(young_n_col), get(young_e_col)),
    `<60 HR (95% CI)` = fmt_hr_safe(get(young_hr_col)),
    `60-75: n / events` = fmt_ne(get(mid_n_col), get(mid_e_col)),
    `60-75 HR (95% CI)` = fmt_hr_safe(get(mid_hr_col)),
    `>=75: n / events` = fmt_ne(get(old_n_col), get(old_e_col)),
    `>=75 HR (95% CI)` = fmt_hr_safe(get(old_hr_col)),
    `P-interaction` = fmt_p(get(pinter_col))
  )]
  
  s6_dt[Phenotype == phenotype_labels["DOM_C7"],
        c("<60 HR (95% CI)", "60-75 HR (95% CI)", ">=75 HR (95% CI)") := "1.00 (reference)"]
  
  write_pub_xlsx(
    dt = s6_dt,
    file_path = file.path(OUT_DIR, "TableS6_age_stratified_HR.xlsx"),
    sheet_name = "Table S6",
    title_text = "Table S6. Age-stratified hazard ratios for 4-year all-cause mortality by disease community phenotype",
    footnotes = c(
      paste0("Data are HR (95% CI) versus C7 (Upper-airway & orodental inflammation, reference community). ",
             "Subgroup Cox models adjusted for sex. Age groups defined by clinical cuts: <60, 60\u201374, and \u226575 years."),
      paste0("Global likelihood-ratio test for continuous age \u00d7 community interaction (per decade): ",
             "\u03c7\u00b2 = 850.01, df = 16, P = 1.3 \u00d7 10\u207b\u00b9\u2077\u2070."),
      build_abbrev_footnote(c(
        "CI" = "confidence interval",
        "HR" = "hazard ratio",
        "LRT" = "likelihood-ratio test",
        "MIXED" = "patients without a single dominant community phenotype"
      ))
    ),
    col_widths = c(50, 18, 22, 18, 22, 18, 22, 16),
    left_cols = 1
  )
  cat(sprintf("   \u2713 Saved: TableS6_age_stratified_HR.xlsx\n"))
}

# =========================
# Table S7: ICD-10 chapter enrichment
# =========================
cat("\n========================================================\n")
cat("Table S7: ICD-10 chapter enrichment (full results)\n")
cat("========================================================\n")

if (!file.exists(ENRICH_XLSX)) {
  warning("ENRICH_XLSX not found, skipping Table S7")
} else {
  # Try common sheet names
  enrich_sheets <- getSheetNames(ENRICH_XLSX)
  cat(sprintf("   Available sheets: %s\n", paste(enrich_sheets, collapse = ", ")))
  
  # Prefer "Significant_Enrichment" or similar
  sig_sheet <- enrich_sheets[grepl("[Ss]ignif", enrich_sheets)][1]
  if (is.na(sig_sheet)) sig_sheet <- enrich_sheets[1]
  cat(sprintf("   Using sheet: %s\n", sig_sheet))
  
  s7_raw <- as.data.table(read.xlsx(ENRICH_XLSX, sheet = sig_sheet))
  cat(sprintf("   Rows: %d  Cols: %s\n", nrow(s7_raw),
              paste(names(s7_raw), collapse = ", ")))
  
  # Save as-is (already publication-ready usually); just add styling wrapper
  # First standardize column types
  setnames(s7_raw,
           names(s7_raw),
           gsub("\\.", " ", names(s7_raw)))
  
  write_pub_xlsx(
    dt = s7_raw,
    file_path = file.path(OUT_DIR, "TableS7_ICD10_enrichment_full.xlsx"),
    sheet_name = "Table S7",
    title_text = "Table S7. Full results of ICD-10 chapter enrichment analysis across 16 disease communities",
    footnotes = c(
      paste0("Chapter enrichment was assessed using a hypergeometric test comparing the observed proportion of ICD-10 chapter membership ",
             "within each community against the chapter distribution in the full network background of 613 disease nodes. ",
             "P values were adjusted for multiple testing using the Benjamini-Hochberg false discovery rate procedure across all community-by-chapter combinations."),
      paste0("Higher observed-over-expected ratios indicate stronger over-representation of a chapter within the community. ",
             "Communities are referenced by their numeric label (C1\u2013C16) as defined in the main analysis."),
      build_abbrev_footnote(c(
        "FDR" = "false discovery rate (Benjamini-Hochberg)",
        "ICD-10" = "International Classification of Diseases, 10th revision",
        "O/E" = "observed-to-expected ratio"
      ))
    ),
    col_widths = NULL,
    left_cols = 2
  )
  cat(sprintf("   \u2713 Saved: TableS7_ICD10_enrichment_full.xlsx\n"))
}

# =========================
# Table S8: 7 sensitivity analyses
# =========================
cat("\n========================================================\n")
cat("Table S8: 7 sensitivity analyses (Model 2 HR pivoted)\n")
cat("========================================================\n")

s8_raw <- as.data.table(read.xlsx(COX_XLSX, sheet = "Sensitivity_all"))
cat(sprintf("   Rows: %d (expected ~452)\n", nrow(s8_raw)))

# Filter to Model 2 only (most informative; full 4 models would be too wide)
s8_m2 <- s8_raw[grepl("Model 2", model)]
cat(sprintf("   Model 2 rows: %d\n", nrow(s8_m2)))

# Format HR (95% CI) per row
s8_m2[, HR_str := fmt_hr(as.numeric(HR),
                         as.numeric(HR_lower),
                         as.numeric(HR_upper), 2)]

# Pivot: rows = term, columns = analysis
analyses_order <- c("S1: Grace censoring", "S2: tau=0.35", "S3: tau=0.45",
                    "S4: Weighted score", "S5: Strict dominant",
                    "S6: Include NONE", "S7: Merge NONE")
# Match user's actual analysis labels (substr-based)
unique_analyses <- unique(s8_m2$analysis)
cat(sprintf("   Unique sensitivity analyses: %s\n",
            paste(unique_analyses, collapse = " | ")))

s8_wide <- dcast(s8_m2, term ~ analysis, value.var = "HR_str", fill = "\u2014")
cat(sprintf("   Wide table dimensions: %d rows x %d cols\n",
            nrow(s8_wide), ncol(s8_wide)))

# CRITICAL: dcast does NOT create a DOM_C7 row because it's the REF
# (Cox doesn't output a coefficient for the reference level). We must add it.
analysis_cols <- setdiff(names(s8_wide), "term")
c7_row <- as.list(c("DOM_C7", rep("1.00 (reference)", length(analysis_cols))))
names(c7_row) <- c("term", analysis_cols)
s8_wide <- rbind(s8_wide, as.data.table(c7_row))
cat(sprintf("   Added DOM_C7 row. New rows: %d (expected 17 = 16 + REF)\n",
            nrow(s8_wide)))

s8_wide[, term_factor := factor(term, levels = display_order)]
setorder(s8_wide, term_factor)
s8_wide[, term_factor := NULL]

# Replace term with phenotype label
s8_wide[, Phenotype := phenotype_labels[term]]
s8_wide[is.na(Phenotype), Phenotype := term]
setcolorder(s8_wide, c("Phenotype", setdiff(names(s8_wide), c("term", "Phenotype"))))
s8_wide[, term := NULL]

# Clean sensitivity column headers: "S1: Grace censoring" -> "Grace censoring"
# "S2: tau=0.35" -> "\u03c4 = 0.35"
clean_sens_name <- function(x) {
  out <- gsub("^S\\d+:\\s*", "", x)              # drop "S1: " etc
  out <- gsub("^S_", "", out)                     # drop "S_" prefix
  out <- gsub("tau\\s*=\\s*", "\u03c4 = ", out)  # tau -> \u03c4
  out <- gsub("includeNONE", "include unclassified", out)
  out <- gsub("mergeNONE", "merge unclassified into MIXED", out)
  out
}
new_names <- sapply(names(s8_wide), clean_sens_name, USE.NAMES = FALSE)
setnames(s8_wide, names(s8_wide), new_names)

write_pub_xlsx(
  dt = s8_wide,
  file_path = file.path(OUT_DIR, "TableS8_sensitivity_analyses.xlsx"),
  sheet_name = "Table S8",
  title_text = "Table S8. Hazard ratios for 4-year all-cause mortality across seven sensitivity analyses (Model 2: adjusted for age and sex)",
  footnotes = c(
    paste0("Data are HR (95% CI) for community phenotype against the reference (C7: Upper-airway & orodental inflammation). ",
           "Each column represents a separate sensitivity analysis testing robustness to alternative analytical decisions:"),
    paste0("\u2022 Grace censoring: 90-day grace period applied to baseline incident events. ",
           "\u2022 \u03c4 = 0.35 / \u03c4 = 0.45: lower / higher community-share threshold for dominance assignment. ",
           "\u2022 Weighted score: community burden score weighted by per-disease mortality risk. ",
           "\u2022 Strict dominant: only patients with single-community dominance \u2265 0.50 retained. ",
           "\u2022 Include unclassified: previously excluded NONE group analysed as a separate stratum. ",
           "\u2022 Merge unclassified into MIXED: NONE group reclassified into MIXED."),
    paste0("Spearman rank correlation \u03c1 of point estimates across all seven sensitivity analyses ranged from 0.95 to 0.99, ",
           "indicating that the primary findings are robust to these analytical choices."),
    build_abbrev_footnote(c(
      "CI" = "confidence interval",
      "HR" = "hazard ratio",
      "MIXED" = "patients without a single dominant community phenotype",
      "NONE" = "patients without any qualifying community membership at baseline"
    ))
  ),
  col_widths = c(50, rep(22, ncol(s8_wide) - 1)),
  left_cols = 1
)
cat(sprintf("   \u2713 Saved: TableS8_sensitivity_analyses.xlsx\n"))

# =========================
# Table S9: PH assumption test
# =========================
cat("\n========================================================\n")
cat("Table S9: PH assumption test\n")
cat("========================================================\n")

s9_raw <- as.data.table(read.xlsx(COX_XLSX, sheet = "PH_test"))
cat(sprintf("   Rows: %d  Cols: %s\n", nrow(s9_raw),
            paste(names(s9_raw), collapse = ", ")))
cat(sprintf("   Unique terms: %s\n",
            paste(unique(s9_raw$term), collapse = " | ")))

# PH_test contains: 15 communities (DOM_C1-C16 except C7=REF) + MIXED
# + age_at_index, sex_binary, grp (joint), GLOBAL = 20 rows total.
# My display_order has DOM_C7 first which is NOT in PH_test - this is correct,
# we don't show a row for the reference. But factor levels must include
# the special covariates / global rows.

# Build proper factor levels including special terms in fixed order
ph_special_order <- c("age_at_index", "sex_binary", "grp (joint)", "GLOBAL")
ph_full_levels <- c(setdiff(display_order, "DOM_C7"),  # 16 community terms
                    ph_special_order)
s9_raw[, term_factor := factor(term, levels = ph_full_levels)]
setorder(s9_raw, term_factor)
s9_raw[, term_factor := NULL]

# Build Phenotype labels for special rows manually
special_labels <- c(
  "age_at_index" = "Age at index date (continuous)",
  "sex_binary"   = "Sex (male vs female)",
  "grp (joint)"  = "Joint test across all community terms",
  "GLOBAL"       = "Global test (all model covariates)"
)
s9_raw[, Phenotype := phenotype_labels[term]]
s9_raw[is.na(Phenotype) & term %in% names(special_labels),
       Phenotype := special_labels[term]]
s9_raw[is.na(Phenotype), Phenotype := term]   # final fallback

# Format columns (handle NA / missing values gracefully)
s9_raw[, chisq_str := ifelse(is.na(chisq), "\u2014",
                             sprintf("%.2f", as.numeric(chisq)))]
s9_raw[, df_str    := ifelse(is.na(df), "\u2014", as.character(df))]
s9_raw[, p_num     := suppressWarnings(as.numeric(p))]
s9_raw[, p_str     := fmt_p(p_num)]
s9_raw[, rho_num   := suppressWarnings(as.numeric(rho))]
s9_raw[, rho_str   := ifelse(is.na(rho_num), "\u2014",
                             sprintf("%.4f", rho_num))]
# interpretation may have NA for joint/global rows; replace
s9_raw[is.na(interpretation) | interpretation == "",
       interpretation := "Joint test (no individual interpretation)"]

s9_dt <- s9_raw[, .(
  Term = Phenotype,
  chisq_col = chisq_str,
  df = df_str,
  p_value_col = p_str,
  rho_col = rho_str,
  Interpretation = interpretation
)]
# Set unicode column names after table built (avoid backtick + \u issues)
setnames(s9_dt,
         c("chisq_col", "p_value_col", "rho_col"),
         c("\u03c7\u00b2", "P value", "\u03c1 (correlation)"))

write_pub_xlsx(
  dt = s9_dt,
  file_path = file.path(OUT_DIR, "TableS9_PH_assumption_test.xlsx"),
  sheet_name = "Table S9",
  title_text = "Table S9. Test of proportional hazards assumption based on Schoenfeld residuals",
  footnotes = c(
    paste0("The Schoenfeld residual test (Grambsch-Therneau) was computed for every covariate ",
           "in the primary Cox proportional hazards model (Model 2: community phenotype + age + sex). ",
           "A small P value suggests the proportional hazards assumption may not hold for that covariate."),
    paste0("In a large cohort (n = 640,637), even minor deviations from proportional hazards may yield small P values; ",
           "interpretation should therefore consider the magnitude of \u03c1 (the correlation between scaled Schoenfeld residuals and log time). ",
           "An absolute value of \u03c1 below 0.05 is generally considered negligible."),
    paste0("The Joint test row evaluates non-proportionality across all 16 community terms simultaneously; ",
           "the Global test row evaluates non-proportionality across all model covariates jointly."),
    build_abbrev_footnote(c(
      "df" = "degrees of freedom",
      "MIXED" = "patients without a single dominant community phenotype",
      "PH" = "proportional hazards",
      "\u03c1" = "Pearson correlation between scaled Schoenfeld residuals and log time",
      "\u03c7\u00b2" = "chi-squared test statistic"
    ))
  ),
  col_widths = c(50, 12, 8, 18, 22, 30),
  left_cols = 1
)
cat(sprintf("   \u2713 Saved: TableS9_PH_assumption_test.xlsx\n"))

# =========================
# Final summary
# =========================
cat("\n=========================================================\n")
cat("MASTER SCRIPT COMPLETE\n")
cat("=========================================================\n\n")
cat(sprintf("Output directory: %s\n\n", OUT_DIR))
out_files <- list.files(OUT_DIR, pattern = "^TableS[2-9]_",
                        full.names = FALSE)
cat(sprintf("Generated %d files:\n", length(out_files)))
for (f in sort(out_files)) cat(sprintf("  - %s\n", f))