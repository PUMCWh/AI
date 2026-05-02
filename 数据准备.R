# =============================================================
# CKM 共病研究 — 步骤1: 数据准备 v3 (DuckDB加速版)
# =============================================================
# 与v2的区别: 使用DuckDB直接查询parquet, 不加载全部数据到R内存
# 6400万行数据预计5分钟内完成 (v2需要1小时+)
#
# 编码定义依据:
#   AHA 2023 Presidential Advisory (Circulation 148:1606-1635)
#   JACC Advances 2025 CKM-Attributable Mortality Study
#
# 相比原稿(Supplementary Table S1)的编码变更:
#   [心血管域] 删除I42(心肌病)和I49(其他心律失常) — AHA Stage 4未列举
#              脑血管病扩展为I60-I69全范围 (含I68)
#   [肾脏域]   保留N17(AKI) — JACC 2025纳入, AKI→CKD转化证据充分
#   [代谢域]   E88.81→E16.803 (中国版ICD-10 GB/T 14396)
#              删除R73 (上游排除第18章)
#              删除E10/E11冗余子码(E10.00等), 保留E10-E14全3位码范围
#              保留E10-E14全范围(含T1DM): 中国EHR中E14使用极频繁
#
# 两层分离设计:
#   层1 (入组): 至少1次CKM相关诊断 → 进入CKM分析队列
#   层2 (网络): 所有疾病统一严格规则 → 构建患者-疾病矩阵
# =============================================================

rm(list = ls()); gc()

library(data.table)
library(arrow)
library(DBI)
library(duckdb)

setDTthreads(0)
t0 <- Sys.time()

cat("=============================================================\n")
cat("  CKM 共病研究 — 数据准备 v3 (DuckDB加速版)\n")
cat("=============================================================\n\n")

# =============================================
# 0) 参数
# =============================================
INPUT_PARQUET <- "C:/Users/HP/Desktop/multimorbidity_test/data_for_analysis/统一数据集Final/数据Final/筛选诊断时间_pre2025/post_clean_neoplasm/unified_pre2025_grace_FINAL.parquet"
OUT_DIR       <- "C:/Users/HP/Desktop/paper/CKM_FINAL"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_OP_VISITS   <- 2L
MIN_OP_GAP_DAYS <- 30L
MIN_IP_VISITS   <- 1L
MIN_AGE         <- 18

STUDY_START   <- as.Date("2016-01-01")
STUDY_END     <- as.Date("2024-12-31")
DISCOVERY_END <- as.Date("2021-12-31")

# =============================================
# 1) CKM编码定义
# =============================================
cat("▶ 1) 构建CKM编码表...\n")
# 编码依据: AHA 2023 Presidential Advisory (Circulation 148:1606-1635)
# 参考校验: JACC Advances 2025 CKM-Attributable Mortality Study

# --- 心血管域 (Cardiovascular) ---
# AHA Stage 4 clinical CVD: CHD, HF, Stroke, PAD, AF
# + Stage 2 hypertension (已发表CKM-ICD研究均归入心血管域)
keep3_cardio <- c(
  "I10", "I11", "I13", "I15",           # 高血压 (AHA Stage 2, 但文献惯例归CV)
  sprintf("I%02d", 20:25),              # I20-I25 缺血性心脏病 (CHD)
  "I48",                                # 房颤 (AF) — AHA仅列AF, 不含I49其他心律失常
  "I50",                                # 心力衰竭 (HF)
  sprintf("I%02d", 60:69),              # I60-I69 脑血管病全范围 (Stroke + sequelae)
  sprintf("I%02d", 70:79)               # I70-I79 外周血管病 (PAD)
)
# 注: I42(心肌病)和I49(其他心律失常)已删除 — AHA Stage 4未列举, JACC 2025未纳入
# 注: I12(高血压性肾病)归入肾脏域, 不在此处

# --- 肾脏域 (Kidney) ---
# AHA: moderate-to-high risk CKD; JACC 2025 使用 N17-N19
keep3_renal <- c(
  "N17",                                # 急性肾损伤 (JACC 2025纳入; AKI→CKD转化有充分证据)
  "N18",                                # 慢性肾脏病
  "N19",                                # 未特指肾衰竭
  "N03", "N04", "N05", "N08",           # 肾炎/肾病综合征 (CKD重要病因)
  "I12",                                # 高血压性肾病
  "N25"                                 # 肾小管功能障碍
)

# --- 代谢域 (Metabolic) ---
# AHA Stage 2: T2DM, obesity, hypertriglyceridemia, metabolic syndrome
# E10-E14保留全范围: 中国EHR中E14(未特指糖尿病)使用极频繁, 删除会丢失大量患者
# Methods中注明: sensitivity analysis using E11 only
keep3_meta <- c(sprintf("E%02d", 10:14), "E66", "E78")

# 4位码 (中国版ICD-10 GB/T 14396)
keep4_meta <- c("E16.803")  # 代谢综合征

exclude_prefix_3 <- c("I23")
exclude_exact_4  <- c("I46.9")

all_ckm_core_3 <- unique(c(keep3_cardio, keep3_renal, keep3_meta))

cat("  心血管域:", length(keep3_cardio), "个3位码\n")
cat("  肾脏域:  ", length(keep3_renal), "个3位码\n")
cat("  代谢域:  ", length(keep3_meta), "个3位码 + E16.803\n")

# =============================================
# 2) 连接DuckDB, 直接读parquet
# =============================================
cat("\n▶ 2) 连接DuckDB...\n")
t_start <- Sys.time()

con <- dbConnect(duckdb(), dbdir = ":memory:")
dbExecute(con, "PRAGMA threads=16;")
dbExecute(con, "PRAGMA memory_limit='8GB';")

p_in <- normalizePath(INPUT_PARQUET, winslash = "/", mustWork = TRUE)
p_in <- gsub("'", "''", p_in)

# 创建视图 (不加载数据到内存)
dbExecute(con, sprintf("
CREATE OR REPLACE VIEW diag AS
SELECT
  CAST(\"人员ID\" AS VARCHAR) AS pid,
  upper(substr(trim(CAST(\"ICD_3\" AS VARCHAR)), 1, 3)) AS icd3,
  trim(CAST(\"疾病编码_update\" AS VARCHAR)) AS icd_full,
  CAST(\"诊断日期\" AS DATE) AS dx_date,
  CAST(\"诊断场景\" AS VARCHAR) AS dx_scene,
  CAST(\"age_at_diagnose\" AS DOUBLE) AS age_dx,
  CAST(\"性别\" AS VARCHAR) AS sex,
  CAST(\"出生日期\" AS DATE) AS birth_date
FROM read_parquet('%s')
WHERE \"人员ID\" IS NOT NULL
  AND \"ICD_3\" IS NOT NULL
  AND \"诊断日期\" IS NOT NULL
  AND CAST(\"age_at_diagnose\" AS DOUBLE) >= %d;
", p_in, MIN_AGE))

# 排除码
dbExecute(con, sprintf("
CREATE OR REPLACE VIEW diag_clean AS
SELECT * FROM diag
WHERE icd3 NOT IN (%s)
  AND NOT starts_with(icd_full, 'I46.9');
", paste0("'", exclude_prefix_3, "'", collapse = ",")))

n_total <- dbGetQuery(con, "SELECT COUNT(*) AS n FROM diag_clean")$n
n_patients <- dbGetQuery(con, "SELECT COUNT(DISTINCT pid) AS n FROM diag_clean")$n
cat("  清洗后记录:", format(n_total, big.mark = ","), "\n")
cat("  唯一患者:  ", format(n_patients, big.mark = ","), "\n")
cat("  耗时:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "秒\n")

# =============================================
# 3) 层1: CKM入组
# =============================================
cat("\n▶ 3) 层1: CKM队列入组...\n")
t_step <- Sys.time()

# 写入CKM编码表到DuckDB
ckm_codes_dt <- data.table(
  icd3 = all_ckm_core_3,
  domain = fifelse(all_ckm_core_3 %in% keep3_cardio, "C",
                   fifelse(all_ckm_core_3 %in% keep3_renal, "K", "M"))
)
# 处理同时属于多个域的编码 (如I12同时在cardio和renal)
ckm_codes_dt <- rbindlist(list(
  data.table(icd3 = keep3_cardio, domain = "C"),
  data.table(icd3 = keep3_renal,  domain = "K"),
  data.table(icd3 = keep3_meta,   domain = "M")
))
dbWriteTable(con, "ckm_codes", ckm_codes_dt, overwrite = TRUE)

# 3位码入组
ckm_patients_3 <- dbGetQuery(con, "
SELECT DISTINCT pid FROM diag_clean
WHERE icd3 IN (SELECT DISTINCT icd3 FROM ckm_codes);
")

# 4位码入组 (E16.803 代谢综合征)
ckm_patients_4 <- dbGetQuery(con, sprintf("
SELECT DISTINCT pid FROM diag_clean
WHERE starts_with(icd_full, '%s');
", keep4_meta[1]))

ckm_pids <- unique(c(ckm_patients_3$pid, ckm_patients_4$pid))
cat("  CKM队列患者数:", format(length(ckm_pids), big.mark = ","), "\n")

# 写入CKM患者ID到DuckDB
dbWriteTable(con, "ckm_pids", data.table(pid = ckm_pids), overwrite = TRUE)

# CKM域标记 (每个患者有哪些域)
ckm_domain <- as.data.table(dbGetQuery(con, "
SELECT
  d.pid,
  MAX(CASE WHEN c.domain = 'C' THEN 1 ELSE 0 END) AS CKM_C,
  MAX(CASE WHEN c.domain = 'K' THEN 1 ELSE 0 END) AS CKM_K,
  MAX(CASE WHEN c.domain = 'M' THEN 1 ELSE 0 END) AS CKM_M
FROM diag_clean d
JOIN ckm_codes c ON d.icd3 = c.icd3
WHERE d.pid IN (SELECT pid FROM ckm_pids)
GROUP BY d.pid;
"))

# 补充E16.803的代谢域
e16_pids <- dbGetQuery(con, sprintf("
SELECT DISTINCT pid FROM diag_clean
WHERE pid IN (SELECT pid FROM ckm_pids)
  AND starts_with(icd_full, '%s');
", keep4_meta[1]))
if (nrow(e16_pids) > 0) {
  ckm_domain[pid %in% e16_pids$pid, CKM_M := 1L]
  cat("  E16.803(代谢综合征)补充患者数:", nrow(e16_pids), "\n")
}

ckm_domain[, n_ckm_domains := CKM_C + CKM_K + CKM_M]
ckm_map <- c("100" = 1L, "010" = 2L, "001" = 3L,
             "110" = 4L, "101" = 5L, "011" = 6L, "111" = 7L)
ckm_domain[, CKM := ckm_map[paste0(CKM_C, CKM_K, CKM_M)]]
ckm_domain <- ckm_domain[!is.na(CKM)]

cat("  域=1:", format(sum(ckm_domain$n_ckm_domains == 1), big.mark = ","),
    " 域=2:", format(sum(ckm_domain$n_ckm_domains == 2), big.mark = ","),
    " 域=3:", format(sum(ckm_domain$n_ckm_domains == 3), big.mark = ","), "\n")
cat("  耗时:", round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "秒\n")

# =============================================
# 4) 层2: 统一严格疾病定义 (全部在DuckDB中完成)
# =============================================
cat("\n▶ 4) 层2: 统一严格疾病定义 (DuckDB)...\n")
t_step <- Sys.time()
cat("  规则: 门诊≥", MIN_OP_VISITS, "次(间隔≥", MIN_OP_GAP_DAYS,
    "天) 或 住院≥", MIN_IP_VISITS, "次\n\n")

# 4.1 住院路径
cat("  4.1 住院路径...\n")
ip_qualified <- as.data.table(dbGetQuery(con, sprintf("
SELECT pid, icd3
FROM diag_clean
WHERE pid IN (SELECT pid FROM ckm_pids)
  AND dx_scene != '门诊'
GROUP BY pid, icd3
HAVING COUNT(*) >= %d;
", MIN_IP_VISITS)))
cat("    住院合格:", format(nrow(ip_qualified), big.mark = ","), "对\n")

# 4.2 门诊路径 (核心优化: SQL完成间隔计算)
cat("  4.2 门诊路径...\n")
op_qualified <- as.data.table(dbGetQuery(con, sprintf("
SELECT pid, icd3
FROM (
  SELECT
    pid,
    icd3,
    COUNT(DISTINCT dx_date) AS n_dates,
    datediff('day', MIN(dx_date), MAX(dx_date)) AS date_span_days
  FROM diag_clean
  WHERE pid IN (SELECT pid FROM ckm_pids)
    AND dx_scene = '门诊'
  GROUP BY pid, icd3
) sub
WHERE n_dates >= %d AND date_span_days >= %d;
", MIN_OP_VISITS, MIN_OP_GAP_DAYS)))
cat("    门诊合格:", format(nrow(op_qualified), big.mark = ","), "对\n")

# 4.3 合并
all_pairs_strict <- unique(rbind(ip_qualified, op_qualified))
cat("  严格定义合格:", format(nrow(all_pairs_strict), big.mark = ","), "对\n")

# 4.4 宽松定义
n_relaxed <- dbGetQuery(con, "
SELECT COUNT(*) AS n FROM (
  SELECT DISTINCT pid, icd3 FROM diag_clean
  WHERE pid IN (SELECT pid FROM ckm_pids)
);
")$n
cat("  宽松定义合格:", format(n_relaxed, big.mark = ","), "对\n")

reduction <- round(100 * (1 - nrow(all_pairs_strict) / n_relaxed), 1)
cat("  过滤减少比例:", reduction, "%\n")
cat("  层2耗时:", round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "秒\n")

# =============================================
# 4.5) 宽松定义数据 (敏感性分析)
# =============================================
cat("\n▶ 4.5) 获取宽松定义数据...\n")
all_pairs_relaxed <- as.data.table(dbGetQuery(con, "
SELECT DISTINCT pid, icd3 FROM diag_clean
WHERE pid IN (SELECT pid FROM ckm_pids);
"))

# =============================================
# 5) 患者层面汇总
# =============================================
cat("\n▶ 5) 患者层面汇总...\n")

ps_strict <- all_pairs_strict[, .(n_conditions_strict = uniqueN(icd3)), by = pid]
ps_relaxed <- all_pairs_relaxed[, .(n_conditions_relaxed = uniqueN(icd3)), by = pid]

# 统一列名 (与DuckDB输出的pid对齐)
setnames(ckm_domain, "pid", "pid")
ckm_patients <- merge(ckm_domain, ps_strict, by = "pid", all.x = TRUE)
ckm_patients <- merge(ckm_patients, ps_relaxed, by = "pid", all.x = TRUE)
ckm_patients[is.na(n_conditions_strict), n_conditions_strict := 0L]
ckm_patients[is.na(n_conditions_relaxed), n_conditions_relaxed := 0L]

cat("  CKM患者总数:", format(nrow(ckm_patients), big.mark = ","), "\n")
cat("  严格定义平均疾病数:", round(mean(ckm_patients$n_conditions_strict), 2), "\n")
cat("  宽松定义平均疾病数:", round(mean(ckm_patients$n_conditions_relaxed), 2), "\n")

# =============================================
# 6) 时间分段
# =============================================
cat("\n▶ 6) 时间分段...\n")

first_ckm <- as.data.table(dbGetQuery(con, sprintf("
SELECT pid, MIN(dx_date) AS first_ckm_date
FROM diag_clean
WHERE pid IN (SELECT pid FROM ckm_pids)
  AND (icd3 IN (SELECT DISTINCT icd3 FROM ckm_codes)
       OR starts_with(icd_full, '%s'))
GROUP BY pid;
", keep4_meta[1])))

ckm_patients <- merge(ckm_patients, first_ckm, by = "pid", all.x = TRUE)
ckm_patients[, cohort := fifelse(
  first_ckm_date <= DISCOVERY_END, "discovery", "validation"
)]

n_disc <- sum(ckm_patients$cohort == "discovery", na.rm = TRUE)
n_val  <- sum(ckm_patients$cohort == "validation", na.rm = TRUE)
cat("  发现集(2016-2021):", format(n_disc, big.mark = ","), "\n")
cat("  验证集(2022-2024):", format(n_val, big.mark = ","), "\n")

# =============================================
# 7) 疾病定义对比
# =============================================
cat("\n▶ 7) 疾病定义对比...\n")

freq_strict  <- all_pairs_strict[, .(freq_strict = uniqueN(pid)), by = icd3]
freq_relaxed <- all_pairs_relaxed[, .(freq_relaxed = uniqueN(pid)), by = icd3]

comparison <- merge(freq_relaxed, freq_strict, by = "icd3", all = TRUE)
setnames(comparison, "icd3", "ICD_3")  # 统一列名, 兼容后续脚本
comparison[is.na(freq_strict), freq_strict := 0L]
comparison[is.na(freq_relaxed), freq_relaxed := 0L]
comparison[, reduction_pct := round(100 * (1 - freq_strict / freq_relaxed), 1)]
comparison[freq_relaxed == 0, reduction_pct := NA_real_]
comparison[, is_ckm_core := ICD_3 %in% all_ckm_core_3]
comparison[, n_ckm_patients := nrow(ckm_patients)]
comparison[, prev_strict  := round(freq_strict / nrow(ckm_patients), 6)]
comparison[, prev_relaxed := round(freq_relaxed / nrow(ckm_patients), 6)]
setorder(comparison, -freq_relaxed)

# 重点关注
for (code in c("K29", "N17", "I10", "E11", "E78")) {
  row <- comparison[ICD_3 == code]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: 宽松=%s → 严格=%s (减少%s%%)\n",
                code,
                format(row$freq_relaxed, big.mark = ","),
                format(row$freq_strict, big.mark = ","),
                row$reduction_pct))
  }
}

fwrite(comparison, file.path(OUT_DIR, "disease_definition_comparison.csv"))

# =============================================
# 8) 构建患者-疾病二元矩阵
# =============================================
cat("\n▶ 8) 构建患者-疾病二元矩阵...\n")

build_matrix <- function(pair_dt, id_col = "pid", code_col = "icd3") {
  pair_sub <- unique(pair_dt[, .SD, .SDcols = c(id_col, code_col)])
  setnames(pair_sub, c("pid", "icd3"))
  pair_sub[, pid_idx := as.integer(factor(pid))]
  pair_sub[, icd_idx := as.integer(factor(icd3))]
  
  # 列名兼容下游脚本 (总CKM网络分析.R 期望 ICD_3; 社区映射回个体.R 期望 人员ID)
  disease_map <- unique(pair_sub[, .(icd_idx, ICD_3 = icd3)])
  setorder(disease_map, icd_idx)
  patient_map <- unique(pair_sub[, .(pid_idx, 人员ID = pid)])
  setorder(patient_map, pid_idx)
  
  mat <- Matrix::sparseMatrix(
    i = pair_sub$pid_idx, j = pair_sub$icd_idx,
    dims = c(max(pair_sub$pid_idx), max(pair_sub$icd_idx)), x = 1
  )
  list(matrix = mat, disease_map = disease_map, patient_map = patient_map)
}

mat_strict  <- build_matrix(all_pairs_strict)
mat_relaxed <- build_matrix(all_pairs_relaxed)

cat("  严格:", nrow(mat_strict$matrix), "×", ncol(mat_strict$matrix), "\n")
cat("  宽松:", nrow(mat_relaxed$matrix), "×", ncol(mat_relaxed$matrix), "\n")

# =============================================
# 9) 保存
# =============================================
cat("\n▶ 9) 保存...\n")
t_save <- Sys.time()

# 患者表
# 统一列名为后续脚本兼容 (人员ID)
ckm_out <- copy(ckm_patients)
setnames(ckm_out, "pid", "人员ID")
fwrite(ckm_out, file.path(OUT_DIR, "ckm_patients_strict.csv"))

ckm_out_relaxed <- copy(ckm_out)
ckm_out_relaxed[, n_conditions_strict := NULL]
setnames(ckm_out_relaxed, "n_conditions_relaxed", "n_conditions")
fwrite(ckm_out_relaxed, file.path(OUT_DIR, "ckm_patients_relaxed.csv"))

# 矩阵
saveRDS(mat_strict,  file.path(OUT_DIR, "patient_disease_matrix_strict.rds"))
saveRDS(mat_relaxed, file.path(OUT_DIR, "patient_disease_matrix_relaxed.rds"))

# 记录级数据 (只导出CKM患者的记录, 用DuckDB高效完成)
cat("  导出parquet (严格)...\n")
dbWriteTable(con, "ckm_info", ckm_patients, overwrite = TRUE)

# 构建输出路径 (避免normalizePath对不存在文件的问题)
out_strict_pq <- gsub("\\\\", "/", file.path(OUT_DIR, "analysis_CKM_strict.parquet"))
out_relaxed_pq <- gsub("\\\\", "/", file.path(OUT_DIR, "analysis_CKM_relaxed.parquet"))

dbExecute(con, sprintf("
COPY (
  SELECT d.*, c.CKM, c.CKM_C, c.CKM_K, c.CKM_M,
         c.n_ckm_domains, c.n_conditions_strict, c.cohort
  FROM diag_clean d
  JOIN ckm_info c ON d.pid = c.pid
) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD);
", gsub("'", "''", out_strict_pq)))

cat("  导出parquet (宽松)...\n")
dbExecute(con, sprintf("
COPY (
  SELECT d.*, c.CKM, c.CKM_C, c.CKM_K, c.CKM_M, c.n_ckm_domains
  FROM diag_clean d
  JOIN ckm_info c ON d.pid = c.pid
) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD);
", gsub("'", "''", out_relaxed_pq)))

# CKM编码参考表
ckm_codebook <- rbindlist(list(
  data.table(domain = "Cardiovascular", ICD_3 = keep3_cardio, level = "3-char"),
  data.table(domain = "Kidney",         ICD_3 = keep3_renal,  level = "3-char"),
  data.table(domain = "Metabolic",      ICD_3 = keep3_meta,   level = "3-char"),
  data.table(domain = "Metabolic",      ICD_3 = keep4_meta,   level = "4/5-char")
))
fwrite(ckm_codebook, file.path(OUT_DIR, "ckm_codebook.csv"))

# 队列摘要
summary_dt <- data.table(
  item = c(
    "数据时间范围", "清洗后总记录", "总唯一患者(≥18岁)",
    "CKM队列患者", "心血管域(C)", "肾脏域(K)", "代谢域(M)",
    "域=1", "域=2", "域=3",
    "严格定义(患者,疾病)对", "宽松定义(患者,疾病)对", "过滤减少(%)",
    "严格定义平均疾病数", "宽松定义平均疾病数",
    "严格矩阵疾病维度", "宽松矩阵疾病维度",
    "发现集(2016-2021)", "验证集(2022-2024)",
    "编码变更: R73删除", "编码变更: I42+I49删除", "编码变更: N17保留",
    "编码变更: I60-I69全范围", "编码变更: E88.81→E16.803"
  ),
  value = c(
    paste(STUDY_START, "至", STUDY_END),
    format(n_total, big.mark = ","),
    format(n_patients, big.mark = ","),
    format(nrow(ckm_patients), big.mark = ","),
    format(sum(ckm_patients$CKM_C), big.mark = ","),
    format(sum(ckm_patients$CKM_K), big.mark = ","),
    format(sum(ckm_patients$CKM_M), big.mark = ","),
    format(sum(ckm_patients$n_ckm_domains == 1), big.mark = ","),
    format(sum(ckm_patients$n_ckm_domains == 2), big.mark = ","),
    format(sum(ckm_patients$n_ckm_domains == 3), big.mark = ","),
    format(nrow(all_pairs_strict), big.mark = ","),
    format(n_relaxed, big.mark = ","),
    as.character(reduction),
    as.character(round(mean(ckm_patients$n_conditions_strict), 2)),
    as.character(round(mean(ckm_patients$n_conditions_relaxed), 2)),
    as.character(ncol(mat_strict$matrix)),
    as.character(ncol(mat_relaxed$matrix)),
    format(n_disc, big.mark = ","),
    format(n_val, big.mark = ","),
    "上游排除第18章", "AHA未列举,JACC未纳入", "JACC纳入,AKI→CKD证据充分",
    "含I68,与JACC一致", "中国版ICD-10 (GB/T 14396)"
  )
)
fwrite(summary_dt, file.path(OUT_DIR, "cohort_summary.csv"))

cat("  保存耗时:", round(as.numeric(difftime(Sys.time(), t_save, units = "secs")), 1), "秒\n")

# =============================================
# 10) 断开DuckDB
# =============================================
dbDisconnect(con, shutdown = TRUE)

# =============================================
# 11) 摘要
# =============================================
total_time <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

cat("\n=============================================================\n")
cat("  数据准备 v3 完成  (总用时:", total_time, "分钟)\n")
cat("=============================================================\n\n")

cat("设计:\n")
cat("  层1(入组): ≥1次CKM诊断 → 进入队列\n")
cat("  层2(网络): 所有疾病统一严格规则\n")
cat("  加速: DuckDB直接查询parquet, 无需全量加载到R内存\n\n")

cat("编码表 (中国版ICD-10, 依据AHA 2023 + JACC 2025):\n")
cat("  [删除] R73 (第18章已被上游排除)\n")
cat("  [删除] I42/心肌病, I49/其他心律失常 (AHA Stage 4未列举)\n")
cat("  [保留] N17/AKI (JACC 2025纳入; AKI→CKD证据充分)\n")
cat("  [扩展] I60-I69全范围 (与JACC一致)\n")
cat("  [修正] E88.81→E16.803/代谢综合征 (中国版ICD-10)\n")
cat("  [保留] E10-E14全范围 (中国EHR中E14使用频繁; 敏感性分析用E11)\n\n")

cat("输出目录:", OUT_DIR, "\n")
cat("  analysis_CKM_strict.parquet\n")
cat("  analysis_CKM_relaxed.parquet\n")
cat("  patient_disease_matrix_strict.rds\n")
cat("  patient_disease_matrix_relaxed.rds\n")
cat("  ckm_patients_strict.csv\n")
cat("  ckm_patients_relaxed.csv\n")
cat("  disease_definition_comparison.csv\n")
cat("  ckm_codebook.csv\n")
cat("  cohort_summary.csv\n\n")

print(summary_dt)

cat("\n下一步: 将 patient_disease_matrix_strict.rds 输入网络分析脚本\n")
cat("=============================================================\n")