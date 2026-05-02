#!/usr/bin/env Rscript

# Regenerated analysis + preliminary plotting script (standalone)
# Does NOT source legacy scripts.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
  library(survival)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

CFG <- list(
  input_csv = Sys.getenv("CKM_COX_DATA", unset = "data/cox_dataset_with_dominant_community.csv"),
  out_dir = Sys.getenv("CKM_ANALYSIS_OUT", unset = "outputs/02_analysis_regenerated"),
  tau = as.numeric(Sys.getenv("CKM_TAU", unset = "0.40")),
  ref_group = Sys.getenv("CKM_REF", unset = "DOM_C7")
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

message("[1/5] Loading data: ", CFG$input_csv)
dt <- fread(CFG$input_csv)
stopifnot(all(c("time_years", "event", "dom_community_unw", "age_at_index", "sex_binary") %in% names(dt)))

# ---- Main Cox models ----
message("[2/5] Fitting Cox models")
dt <- dt[dom_community_unw != "NONE"]
dt[, dom_community_unw := factor(dom_community_unw)]
if (!(CFG$ref_group %in% levels(dt$dom_community_unw))) {
  CFG$ref_group <- levels(dt$dom_community_unw)[1]
}
dt[, dom_community_unw := relevel(dom_community_unw, ref = CFG$ref_group)]

m1 <- coxph(Surv(time_years, event) ~ dom_community_unw, data = dt)
m2 <- coxph(Surv(time_years, event) ~ dom_community_unw + age_at_index + sex_binary, data = dt)

extract_hr <- function(fit, model_name) {
  s <- summary(fit)
  tb <- as.data.table(s$coefficients, keep.rownames = "term")
  ci <- as.data.table(s$conf.int, keep.rownames = "term")
  out <- merge(tb, ci[, .(term, `lower .95`, `upper .95`, `exp(coef)`)], by = "term")
  out[, model := model_name]
  out[, .(model, term, hr = `exp(coef)`, lower95 = `lower .95`, upper95 = `upper .95`, p = `Pr(>|z|)`)]
}

cox_res <- rbindlist(list(extract_hr(m1, "M1_unadjusted"), extract_hr(m2, "M2_age_sex")), fill = TRUE)
fwrite(cox_res, file.path(CFG$out_dir, "cox_results_regenerated.csv"))

# ---- Subgroup summary ----
message("[3/5] Subgroup summaries")
dt[, age_group := cut(age_at_index, c(-Inf, 60, 75, Inf), labels = c("<60", "60-75", ">=75"), right = FALSE)]
sub_sum <- dt[, .(n = .N, deaths = sum(event == 1), py = sum(time_years), rate_per_1000py = 1000 * sum(event == 1) / sum(time_years)), by = .(sex_binary, age_group, dom_community_unw)]
fwrite(sub_sum, file.path(CFG$out_dir, "subgroup_summary_regenerated.csv"))

# ---- Preliminary plot: top1 share threshold check ----
message("[4/5] Plotting top1-share distribution")
if (all(c("dom_share_unw", "total_score_unw") %in% names(dt))) {
  pp <- ggplot(dt[total_score_unw > 0], aes(x = dom_share_unw)) +
    geom_histogram(binwidth = 0.02, fill = "#4575b4", color = "white") +
    geom_vline(xintercept = CFG$tau, linetype = "dashed", linewidth = 0.8, color = "#d73027") +
    labs(x = "Top1 share", y = "Count", title = sprintf("Top1-share distribution (tau=%.2f)", CFG$tau)) +
    theme_bw(base_size = 12)
  ggsave(file.path(CFG$out_dir, "prelim_top1_share_tau040_regenerated.pdf"), pp, width = 7, height = 5)
}

# ---- Export workbook ----
message("[5/5] Writing workbook")
wb <- createWorkbook()
addWorksheet(wb, "cox")
writeData(wb, "cox", cox_res)
addWorksheet(wb, "subgroup")
writeData(wb, "subgroup", sub_sum)
saveWorkbook(wb, file.path(CFG$out_dir, "analysis_summary_regenerated.xlsx"), overwrite = TRUE)

message("Done. Outputs in: ", CFG$out_dir)
