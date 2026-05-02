#!/usr/bin/env Rscript
# ==============================================================================
# Community Score Mapping to Individuals
# Produces:
#   1) cox_dataset_baseline_scores.csv   — baseline-period cohort for Cox regression
#   2) profile_dataset_all_scores.csv    — full-period profile for all patients
#   3) QC_and_summaries.xlsx             — quality-control checks
#   4) community_id_reference.csv        — community ID reference table
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  pkgs <- c("data.table", "duckdb", "DBI", "openxlsx")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
    }
  }
  library(data.table)
  library(duckdb)
  library(DBI)
  library(openxlsx)
})

options(stringsAsFactors = FALSE)

# =========================
# 0) Paths & Parameters
# =========================
INPUT_PARQUET <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\analysis_CKM_strict.parquet"
DEATH_PARQUET <- "C:\\Users\\HP\\Desktop\\yichang_data\\死亡数据\\pumc_death_merged_all.parquet"
NETWORK_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\NETWORK_PCOR_20260410_221305"
COMMUNITY_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\community_detection"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Study design: landmark analysis (per docx §2.4.3)
# 5-year baseline exposure window + 4-year follow-up
BASELINE_START <- as.Date("2016-01-01")
BASELINE_END   <- as.Date("2020-12-31")
INDEX_DATE     <- as.Date("2021-01-01")
STUDY_END      <- as.Date("2024-12-31")
PROFILE_START  <- as.Date("2016-01-01")
PROFILE_END    <- as.Date("2024-12-31")

# Threshold (medoid partition) and edge filter for strength
THRESH_TAG    <- "coef0.01_alpha0.05"
EDGE_COEF_MIN <- 0.01
EDGE_ALPHA_MAX <- 0.05

# Optional grace period for sensitivity censoring (days)
GRACE_DAYS <- 365

# =========================
# 1) Locate Community & Network Outputs
# =========================
cat("==> Step 1) Loading community and network outputs...\n")

consensus_rds <- file.path(COMMUNITY_DIR, "consensus_results_all.rds")
if (!file.exists(consensus_rds)) {
  stop("Missing consensus_results_all.rds in COMMUNITY_DIR.", call. = FALSE)
}
consensus_results <- readRDS(consensus_rds)
if (!THRESH_TAG %in% names(consensus_results)) {
  stop("Threshold tag not found in consensus_results_all.rds.", call. = FALSE)
}
medoid_dt <- as.data.table(consensus_results[[THRESH_TAG]]$medoid_dt)
if (is.null(medoid_dt) || nrow(medoid_dt) == 0) {
  stop("medoid_dt missing for the threshold tag.", call. = FALSE)
}

# Standardise column names
if (!("node" %in% names(medoid_dt))) {
  cand <- intersect(names(medoid_dt), c("code", "icd3", "ICD_3", "DiseaseCode"))
  if (length(cand) == 1) setnames(medoid_dt, cand, "node")
}
if (!("community" %in% names(medoid_dt))) {
  cand <- intersect(names(medoid_dt), c("comm", "cluster", "partition", "module"))
  if (length(cand) == 1) setnames(medoid_dt, cand, "community")
}
if (!all(c("node", "community") %in% names(medoid_dt))) {
  stop("medoid_dt must include node/community columns.", call. = FALSE)
}

medoid_dt[, node := toupper(substr(trimws(as.character(node)), 1, 3))]
medoid_dt[, community := as.integer(community)]
medoid_dt <- medoid_dt[!is.na(node) & node != "" & !is.na(community) & community > 0]
medoid_dt <- unique(medoid_dt, by = "node")

# ---- NO renumbering: keep original community IDs from Leiden/medoid ----
cat("QC: community IDs retained from medoid partition:",
    paste(sort(unique(medoid_dt$community)), collapse = ", "), "\n")
cat("QC: number of communities =", length(unique(medoid_dt$community)), "\n")

# Output community reference table for cross-script consistency
comm_ref <- medoid_dt[, .(n_diseases = .N), by = community]
setorder(comm_ref, community)
fwrite(comm_ref, file.path(OUT_DIR, "community_id_reference.csv"))
cat("   Saved community_id_reference.csv\n")

# Load edges
edges_rds <- file.path(NETWORK_DIR, "edges_all_uppertri_with_fdr.rds")
if (!file.exists(edges_rds)) {
  stop("Missing edges_all_uppertri_with_fdr.rds in NETWORK_DIR.", call. = FALSE)
}
edge_dt <- readRDS(edges_rds)

# Load bridge metrics (for node role exposure variables)
bridge_csv <- file.path(COMMUNITY_DIR,
                        sprintf("medoid_bridge_metrics_%s.csv", THRESH_TAG))
if (file.exists(bridge_csv)) {
  bridge_dt <- fread(bridge_csv)
  cat("   Loaded bridge metrics:", nrow(bridge_dt), "nodes\n")
} else {
  # Fallback: try RDS
  bridge_rds <- file.path(COMMUNITY_DIR, "bridge_metrics_all.rds")
  if (file.exists(bridge_rds)) {
    bridge_all <- readRDS(bridge_rds)
    bridge_dt <- bridge_all[threshold == THRESH_TAG]
  } else {
    warning("No bridge metrics found. Role exposure variables will not be computed.")
    bridge_dt <- NULL
  }
}

# Build role lookup: disease code → topological role
if (!is.null(bridge_dt)) {
  role_lookup <- bridge_dt[!is.na(role) & role != "isolate",
                           .(code = toupper(substr(trimws(as.character(node)), 1, 3)),
                             role = as.character(role))]
  role_lookup <- unique(role_lookup, by = "code")
  cat("   Role distribution:\n")
  print(role_lookup[, .N, by = role][order(-N)])
} else {
  role_lookup <- NULL
}

# =========================
# 2) Build Node Strengths (network weights)
# =========================
cat("==> Step 2) Computing node strengths...\n")

edge_f <- edge_dt[p_adj < EDGE_ALPHA_MAX & weight >= EDGE_COEF_MIN]
edge_f <- edge_f[weight > 0]

strength_dt <- rbindlist(list(
  edge_f[, .(code = from, w = weight)],
  edge_f[, .(code = to,   w = weight)]
))[, .(strength = sum(w, na.rm = TRUE)), by = code]

comm_map <- merge(medoid_dt, strength_dt, by.x = "node", by.y = "code", all.x = TRUE)
comm_map[is.na(strength), strength := 0]
comm_map <- comm_map[community > 0]
comm_map_db <- comm_map[, .(
  code      = toupper(substr(trimws(as.character(node)), 1, 3)),
  community = as.integer(community),
  strength  = as.numeric(strength)
)]
comm_map_db <- comm_map_db[!is.na(code) & code != "" & !is.na(community) & community > 0]

comm_totals <- comm_map_db[, .(
  comm_size      = .N,
  strength_total = sum(strength, na.rm = TRUE)
), by = community]

# =========================
# 3) DuckDB: load diagnosis data
# =========================
cat("==> Step 3) DuckDB: loading diagnosis data...\n")

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
DBI::dbExecute(con, "PRAGMA threads=16;")
DBI::dbExecute(con, "PRAGMA enable_progress_bar=true;")

p_in <- normalizePath(INPUT_PARQUET, winslash = "/", mustWork = TRUE)
p_in <- gsub("'", "''", p_in)

DBI::dbExecute(con, sprintf("
CREATE OR REPLACE TEMP VIEW diag AS
SELECT
  upper(CAST(pid AS VARCHAR)) AS pid,
  upper(substr(trim(CAST(icd3 AS VARCHAR)), 1, 3)) AS code,
  CAST(dx_date AS DATE) AS dx_date,
  CAST(sex AS VARCHAR) AS sex,
  CAST(birth_date AS DATE) AS birth_date,
  CAST(CKM AS INTEGER) AS CKM,
  CAST(CKM_C AS INTEGER) AS CKM_C,
  CAST(CKM_K AS INTEGER) AS CKM_K,
  CAST(CKM_M AS INTEGER) AS CKM_M,
  CAST(n_ckm_domains AS INTEGER) AS n_ckm_domains
FROM read_parquet('%s')
WHERE pid IS NOT NULL
  AND icd3 IS NOT NULL
  AND dx_date IS NOT NULL;
", p_in))

# Patient demographics (全期，仅用于年龄/性别)
patient_demo <- as.data.table(DBI::dbGetQuery(con, "
SELECT
  pid,
  any_value(sex) AS sex,
  min(birth_date) AS birth_date
FROM diag
GROUP BY pid;
"))

# Baseline-period CKM domain exposure (避免不朽时间偏倚)
# 仅统计 2016-2020 基线期内发生的 CKM 诊断
cat("==> Step 3b) Computing baseline-period CKM domain exposure (ITB-safe)...\n")
ckm_baseline <- as.data.table(DBI::dbGetQuery(con, sprintf("
SELECT
  pid,
  CAST(max(CKM_C) AS INTEGER) AS CKM_C_baseline,
  CAST(max(CKM_K) AS INTEGER) AS CKM_K_baseline,
  CAST(max(CKM_M) AS INTEGER) AS CKM_M_baseline
FROM diag
WHERE dx_date BETWEEN DATE '%s' AND DATE '%s'
GROUP BY pid;
", BASELINE_START, BASELINE_END)))
ckm_baseline[, n_ckm_domains_baseline :=
               CKM_C_baseline + CKM_K_baseline + CKM_M_baseline]
cat(sprintf("   Baseline-period CKM: n=%d patients with at least one diagnosis\n",
            nrow(ckm_baseline)))
cat(sprintf("     CKM_C_baseline: %d (%.1f%%)\n",
            sum(ckm_baseline$CKM_C_baseline == 1),
            100 * mean(ckm_baseline$CKM_C_baseline == 1)))
cat(sprintf("     CKM_K_baseline: %d (%.1f%%)\n",
            sum(ckm_baseline$CKM_K_baseline == 1),
            100 * mean(ckm_baseline$CKM_K_baseline == 1)))
cat(sprintf("     CKM_M_baseline: %d (%.1f%%)\n",
            sum(ckm_baseline$CKM_M_baseline == 1),
            100 * mean(ckm_baseline$CKM_M_baseline == 1)))

# =========================
# 4) Helper to build community scores
# =========================
build_scores <- function(start_date, end_date, label) {
  cat(sprintf("==> Building scores for %s (%s to %s)...\n", label, start_date, end_date))
  
  dbWriteTable(con, "comm_map", comm_map_db, overwrite = TRUE, temporary = TRUE)
  
  sql_scores <- sprintf("
SELECT
  d.pid,
  cm.community,
  COUNT(DISTINCT d.code) AS n_codes,
  SUM(cm.strength) AS strength_sum
FROM (
  SELECT DISTINCT pid, code
  FROM diag
  WHERE dx_date BETWEEN DATE '%s' AND DATE '%s'
) d
JOIN comm_map cm
  ON d.code = cm.code
GROUP BY d.pid, cm.community;
", start_date, end_date)
  
  score_long <- as.data.table(DBI::dbGetQuery(con, sql_scores))
  if (nrow(score_long) == 0) {
    return(data.table(pid = character(), community = integer(),
                      n_codes = integer(), strength_sum = numeric(),
                      score_unw = numeric(), score_w = numeric()))
  }
  
  score_long <- merge(score_long, comm_totals, by = "community", all.x = TRUE)
  score_long[, score_unw := n_codes / comm_size]
  score_long[, score_w   := fifelse(strength_total > 0,
                                    strength_sum / strength_total, 0)]
  score_long[, .(pid, community, n_codes, strength_sum, score_unw, score_w)]
}

# Helper: widen scores and rename columns
widen_scores <- function(score_long_dt) {
  wide <- dcast(score_long_dt, pid ~ community,
                value.var = c("n_codes", "strength_sum", "score_unw", "score_w"),
                fill = 0)
  cols <- setdiff(names(wide), "pid")
  # Rename: n_codes_3 -> comm3_n, strength_sum_3 -> comm3_str,
  #         score_unw_3 -> comm3_unw, score_w_3 -> comm3_w
  setnames(wide, cols, gsub("^n_codes_(\\d+)$",      "comm\\1_n",   cols))
  setnames(wide, names(wide), gsub("^strength_sum_(\\d+)$", "comm\\1_str", names(wide)))
  setnames(wide, names(wide), gsub("^score_unw_(\\d+)$",    "comm\\1_unw", names(wide)))
  setnames(wide, names(wide), gsub("^score_w_(\\d+)$",      "comm\\1_w",   names(wide)))
  
  num_cols <- setdiff(names(wide), "pid")
  wide[, (num_cols) := lapply(.SD, function(x) round(x, 6)), .SDcols = num_cols]
  wide
}

# =========================
# 5) Baseline scores (for Cox)
# =========================
baseline_scores <- build_scores(BASELINE_START, BASELINE_END, "baseline")
scores_wide <- widen_scores(baseline_scores)

baseline_pids <- as.data.table(DBI::dbGetQuery(con, sprintf("
SELECT DISTINCT pid
FROM diag
WHERE dx_date BETWEEN DATE '%s' AND DATE '%s';
", BASELINE_START, BASELINE_END)))

cat("==> Step 5.5) Computing baseline burden/utilization covariates...\n")
baseline_cov <- as.data.table(DBI::dbGetQuery(con, sprintf("
SELECT
  pid,
  COUNT(*)              AS baseline_rows,
  COUNT(DISTINCT dx_date) AS baseline_dx_days,
  COUNT(DISTINCT code)    AS baseline_unique_codes
FROM diag
WHERE dx_date BETWEEN DATE '%s' AND DATE '%s'
GROUP BY pid;
", BASELINE_START, BASELINE_END)))
baseline_cov <- merge(baseline_pids, baseline_cov, by = "pid", all.x = TRUE)
baseline_cov[is.na(baseline_rows),         baseline_rows := 0]
baseline_cov[is.na(baseline_dx_days),      baseline_dx_days := 0]
baseline_cov[is.na(baseline_unique_codes), baseline_unique_codes := 0]

cox_dt <- merge(baseline_pids, patient_demo, by = "pid", all.x = TRUE)
cox_dt <- merge(cox_dt, scores_wide, by = "pid", all.x = TRUE)
cox_dt <- merge(cox_dt,
                baseline_cov[, .(pid, baseline_rows, baseline_dx_days, baseline_unique_codes)],
                by = "pid", all.x = TRUE)

# Merge baseline-period CKM domain exposure (ITB-safe)
cox_dt <- merge(cox_dt, ckm_baseline, by = "pid", all.x = TRUE)
for (cc in c("CKM_C_baseline", "CKM_K_baseline", "CKM_M_baseline",
             "n_ckm_domains_baseline")) {
  cox_dt[is.na(get(cc)), (cc) := 0L]
}

# Fill NAs in score + covariate columns
score_cols <- setdiff(names(cox_dt), c("pid", "sex", "birth_date"))
cox_dt[, (score_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = score_cols]
cox_dt[is.na(baseline_rows),         baseline_rows := 0]
cox_dt[is.na(baseline_dx_days),      baseline_dx_days := 0]
cox_dt[is.na(baseline_unique_codes), baseline_unique_codes := 0]

# Sex binary
cox_dt[, sex_binary := fifelse(grepl("男|male|^m$|^1$", sex, ignore.case = TRUE), 1L, 0L)]

# =========================
# 5.7) Role exposure variables
# =========================
# For each patient: does their baseline disease portfolio include
# a connector_hub / provincial_hub / kinless / connector disease?
if (!is.null(role_lookup)) {
  cat("==> Step 5.7) Computing role exposure variables...\n")
  
  # Get each patient's baseline diseases
  baseline_diseases <- as.data.table(DBI::dbGetQuery(con, sprintf("
  SELECT DISTINCT pid, code
  FROM diag
  WHERE dx_date BETWEEN DATE '%s' AND DATE '%s';
  ", BASELINE_START, BASELINE_END)))
  
  # Merge with role lookup
  baseline_roles <- merge(baseline_diseases, role_lookup, by = "code", all.x = TRUE)
  
  # Patient-level role exposure (binary: has at least one disease of this role?)
  role_exposure <- baseline_roles[!is.na(role), .(
    has_connector_hub  = as.integer(any(role == "connector_hub")),
    has_provincial_hub = as.integer(any(role == "provincial_hub")),
    has_kinless        = as.integer(any(role == "kinless")),
    has_connector      = as.integer(any(role == "connector")),
    has_any_hub        = as.integer(any(role %in% c("connector_hub", "provincial_hub"))),
    n_connector_hub_diseases  = sum(role == "connector_hub"),
    n_hub_diseases            = sum(role %in% c("connector_hub", "provincial_hub")),
    n_kinless_diseases        = sum(role == "kinless")
  ), by = pid]
  
  cox_dt <- merge(cox_dt, role_exposure, by = "pid", all.x = TRUE)
  
  # Fill NAs (patients with no role-assigned diseases)
  role_cols <- c("has_connector_hub", "has_provincial_hub", "has_kinless",
                 "has_connector", "has_any_hub",
                 "n_connector_hub_diseases", "n_hub_diseases", "n_kinless_diseases")
  for (rc in role_cols) {
    cox_dt[is.na(get(rc)), (rc) := 0L]
  }
  
  cat("   Role exposure:\n")
  cat("     has_connector_hub:", sum(cox_dt$has_connector_hub), "/", nrow(cox_dt), "\n")
  cat("     has_provincial_hub:", sum(cox_dt$has_provincial_hub), "/", nrow(cox_dt), "\n")
  cat("     has_kinless:", sum(cox_dt$has_kinless), "/", nrow(cox_dt), "\n")
  cat("     has_any_hub:", sum(cox_dt$has_any_hub), "/", nrow(cox_dt), "\n")
}

# =========================
# 6) Death data + follow-up
# =========================
cat("==> Step 6) Loading death data and computing follow-up...\n")

death_path <- normalizePath(DEATH_PARQUET, winslash = "/", mustWork = TRUE)
death_path <- gsub("'", "''", death_path)
death_dt <- as.data.table(DBI::dbGetQuery(con, sprintf("
SELECT
  upper(id) AS pid,
  CAST(deathdate AS DATE) AS death_date
FROM read_parquet('%s')
WHERE id IS NOT NULL;
", death_path)))

death_dt <- death_dt[!is.na(death_date)]
death_dt <- death_dt[, .(death_date = min(death_date)), by = pid]

cox_dt <- merge(cox_dt, death_dt, by = "pid", all.x = TRUE)
cox_dt <- cox_dt[is.na(death_date) | death_date >= INDEX_DATE]

cox_dt[, age_at_index := as.numeric(difftime(INDEX_DATE, birth_date, units = "days")) / 365.25]
cox_dt[, event := fifelse(!is.na(death_date) &
                            death_date >= INDEX_DATE &
                            death_date <= STUDY_END, 1L, 0L)]
cox_dt[, stop_date := fifelse(!is.na(death_date) & death_date <= STUDY_END,
                              death_date, STUDY_END)]
cox_dt[, time_days  := as.numeric(difftime(stop_date, INDEX_DATE, units = "days"))]
cox_dt[, time_days  := fifelse(time_days <= 0, 0.5, time_days)]
cox_dt[, time_years := time_days / 365.25]

# Sensitivity censoring (grace after last diagnosis)
last_dx <- as.data.table(DBI::dbGetQuery(con, sprintf("
SELECT pid, max(dx_date) AS last_dx_date
FROM diag
WHERE dx_date <= DATE '%s'
GROUP BY pid;
", STUDY_END)))

cox_dt <- merge(cox_dt, last_dx, by = "pid", all.x = TRUE)
cox_dt <- cox_dt[(last_dx_date + GRACE_DAYS) >= INDEX_DATE]
cox_dt[, grace_censor := as.Date(pmin(as.numeric(STUDY_END),
                                      as.numeric(last_dx_date + GRACE_DAYS),
                                      na.rm = TRUE),
                                 origin = "1970-01-01")]
cox_dt[, stop_date_grace := fifelse(event == 1L, stop_date, grace_censor)]
cox_dt[, time_days_grace  := as.numeric(difftime(stop_date_grace, INDEX_DATE, units = "days"))]
cox_dt[, time_days_grace  := fifelse(time_days_grace <= 0, 0.5, time_days_grace)]
cox_dt[, time_years_grace := time_days_grace / 365.25]
cox_dt[, censored_by_grace := as.integer(event == 0L & stop_date_grace < STUDY_END)]
cox_dt[, lost_like := as.integer(event == 0L & (last_dx_date + GRACE_DAYS) < STUDY_END)]

# QC: all-zero score check
score_mat_cols <- grep("^comm\\d+_(unw|w)$", names(cox_dt), value = TRUE)
if (length(score_mat_cols) > 0) {
  n_all_zero <- cox_dt[, sum(rowSums(.SD) == 0), .SDcols = score_mat_cols]
  cat("QC: patients with ALL community scores = 0 :", n_all_zero, " / ", nrow(cox_dt), "\n")
}

# Z-standardisation
score_cols_unw <- grep("^comm\\d+_unw$", names(cox_dt), value = TRUE)
score_cols_w   <- grep("^comm\\d+_w$",   names(cox_dt), value = TRUE)
z_transform <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}
for (cc in c(score_cols_unw, score_cols_w)) {
  cox_dt[, (paste0(cc, "_z")) := z_transform(get(cc))]
}

# =========================
# 7) Profile scores (all patients, full period)
# =========================
profile_scores <- build_scores(PROFILE_START, PROFILE_END, "profile")
profile_wide   <- widen_scores(profile_scores)

# Full-period utilisation covariates (NOT baseline — matches profile time window)
cat("==> Step 7.5) Computing full-period covariates for profile...\n")
profile_cov <- as.data.table(DBI::dbGetQuery(con, sprintf("
SELECT
  pid,
  COUNT(*)              AS total_rows,
  COUNT(DISTINCT dx_date) AS total_dx_days,
  COUNT(DISTINCT code)    AS total_unique_codes
FROM diag
WHERE dx_date BETWEEN DATE '%s' AND DATE '%s'
GROUP BY pid;
", PROFILE_START, PROFILE_END)))

profile_dt <- merge(patient_demo, profile_wide, by = "pid", all.x = TRUE)
profile_dt <- merge(profile_dt, profile_cov, by = "pid", all.x = TRUE)
profile_dt[is.na(total_rows),         total_rows := 0]
profile_dt[is.na(total_dx_days),      total_dx_days := 0]
profile_dt[is.na(total_unique_codes), total_unique_codes := 0]
profile_dt[, age_at_index := as.numeric(difftime(INDEX_DATE, birth_date, units = "days")) / 365.25]
profile_dt[, sex_binary := fifelse(grepl("男|male|^m$|^1$", sex, ignore.case = TRUE), 1L, 0L)]

# =========================
# 8) Output
# =========================
cat("==> Step 8) Writing outputs...\n")

fwrite(cox_dt,     file.path(OUT_DIR, "cox_dataset_baseline_scores.csv"))
fwrite(profile_dt, file.path(OUT_DIR, "profile_dataset_all_scores.csv"))

wb <- createWorkbook()
qc_dt <- data.table(
  baseline_start      = as.character(BASELINE_START),
  baseline_end        = as.character(BASELINE_END),
  index_date          = as.character(INDEX_DATE),
  study_end           = as.character(STUDY_END),
  grace_days          = GRACE_DAYS,
  n_cohort            = nrow(cox_dt),
  n_event             = sum(cox_dt$event == 1L, na.rm = TRUE),
  n_censored_by_grace = sum(cox_dt$censored_by_grace == 1L, na.rm = TRUE),
  n_lost_like         = sum(cox_dt$lost_like == 1L, na.rm = TRUE),
  n_neg_fu            = sum(cox_dt$time_days < 0, na.rm = TRUE),
  n_neg_fu_grace      = sum(cox_dt$time_days_grace < 0, na.rm = TRUE),
  n_missing_birth     = sum(is.na(cox_dt$birth_date)),
  n_missing_sex       = sum(is.na(cox_dt$sex) | cox_dt$sex == ""),
  n_CKM_C_baseline    = sum(cox_dt$CKM_C_baseline == 1L, na.rm = TRUE),
  n_CKM_K_baseline    = sum(cox_dt$CKM_K_baseline == 1L, na.rm = TRUE),
  n_CKM_M_baseline    = sum(cox_dt$CKM_M_baseline == 1L, na.rm = TRUE),
  n_male              = sum(cox_dt$sex_binary == 1L, na.rm = TRUE),
  male_pct            = round(100 * mean(cox_dt$sex_binary == 1L, na.rm = TRUE), 2),
  n_profile           = nrow(profile_dt)
)
addWorksheet(wb, "QC_summary")
writeDataTable(wb, "QC_summary", qc_dt)

addWorksheet(wb, "cox_preview_50k")
writeDataTable(wb, "cox_preview_50k", cox_dt[1:min(.N, 50000)])

prof_sum <- profile_dt[, .(
  n          = .N,
  mean_age   = mean(age_at_index, na.rm = TRUE),
  male_pct   = mean(sex == "男性", na.rm = TRUE),
  missing_birth = sum(is.na(birth_date)),
  missing_sex   = sum(is.na(sex) | sex == "")
)]
addWorksheet(wb, "profile_summary")
writeDataTable(wb, "profile_summary", prof_sum)

# Community reference
addWorksheet(wb, "community_reference")
writeDataTable(wb, "community_reference", comm_ref)

saveWorkbook(wb, file.path(OUT_DIR, "QC_and_summaries.xlsx"), overwrite = TRUE)

DBI::dbDisconnect(con, shutdown = TRUE)

cat("\n========================================\n")
cat("QC: negative follow-up (time_days) = ", sum(cox_dt$time_days < 0, na.rm = TRUE), "\n")
cat("QC: negative follow-up grace       = ", sum(cox_dt$time_days_grace < 0, na.rm = TRUE), "\n")
cat("QC: deaths in follow-up            = ", sum(cox_dt$event == 1L, na.rm = TRUE), "\n")
cat("QC: cohort n (Cox)                 = ", nrow(cox_dt), "\n")
cat("QC: profile n                      = ", nrow(profile_dt), "\n")
cat("========================================\n")
cat("DONE. Outputs saved to:\n", OUT_DIR, "\n", sep = "")