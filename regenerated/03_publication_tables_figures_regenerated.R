#!/usr/bin/env Rscript

# Regenerated publication tables/figures script (standalone)
# Does NOT source legacy scripts.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
  library(survival)
})

CFG <- list(
  input_csv = Sys.getenv("CKM_COX_DATA", unset = "data/cox_dataset_with_dominant_community.csv"),
  out_dir = Sys.getenv("CKM_PUB_OUT", unset = "outputs/03_publication_regenerated"),
  ref_group = Sys.getenv("CKM_REF", unset = "DOM_C7")
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
dt <- fread(CFG$input_csv)

stopifnot(all(c("time_years", "event", "dom_community_unw", "age_at_index", "sex_binary") %in% names(dt)))

dt <- dt[dom_community_unw != "NONE"]
dt[, dom_community_unw := factor(dom_community_unw)]
if (CFG$ref_group %in% levels(dt$dom_community_unw)) dt[, dom_community_unw := relevel(dom_community_unw, ref = CFG$ref_group)]

# --- Table 1 (baseline characteristics) ---
t1 <- dt[, .(
  N = .N,
  Age_mean = mean(age_at_index, na.rm = TRUE),
  Age_sd = sd(age_at_index, na.rm = TRUE),
  Deaths = sum(event == 1, na.rm = TRUE),
  Followup_median = median(time_years, na.rm = TRUE)
), by = .(Sex = fifelse(sex_binary == 1, "Male", "Female"))]
fwrite(t1, file.path(CFG$out_dir, "Table1_baseline_regenerated.csv"))

# --- Table 3-like HR table ---
m1 <- coxph(Surv(time_years, event) ~ dom_community_unw, data = dt)
m2 <- coxph(Surv(time_years, event) ~ dom_community_unw + age_at_index + sex_binary, data = dt)

mk <- function(fit, lbl) {
  s <- summary(fit)
  x <- as.data.table(s$conf.int, keep.rownames = "term")
  p <- as.data.table(s$coefficients, keep.rownames = "term")[, .(term, p = `Pr(>|z|)`)]
  y <- merge(x, p, by = "term")
  y[, model := lbl]
  y[, .(model, term, HR = `exp(coef)`, lower95 = `lower .95`, upper95 = `upper .95`, p)]
}

hr <- rbindlist(list(mk(m1, "Model1"), mk(m2, "Model2")), fill = TRUE)
fwrite(hr, file.path(CFG$out_dir, "Table3_HR_regenerated.csv"))

# --- KM figure (publication-style basic) ---
fit <- survfit(Surv(time_years, event) ~ dom_community_unw, data = dt)
ss <- summary(fit)
km <- data.table(time = ss$time, surv = ss$surv, strata = gsub("dom_community_unw=", "", ss$strata))

p <- ggplot(km, aes(x = time, y = surv, color = strata, group = strata)) +
  geom_step(linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 1)) +
  labs(x = "Years", y = "Survival probability", color = "Phenotype", title = "Kaplan-Meier survival by dominant community") +
  theme_bw(base_size = 11)

ggsave(file.path(CFG$out_dir, "Fig3_KM_regenerated.pdf"), p, width = 11, height = 8)

# --- Export workbook ---
wb <- createWorkbook()
addWorksheet(wb, "Table1")
writeData(wb, "Table1", t1)
addWorksheet(wb, "Table3_HR")
writeData(wb, "Table3_HR", hr)
saveWorkbook(wb, file.path(CFG$out_dir, "publication_tables_regenerated.xlsx"), overwrite = TRUE)

message("Done. Outputs in: ", CFG$out_dir)
