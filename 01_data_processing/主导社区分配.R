#!/usr/bin/env Rscript
# ==============================================================================
# Step 4.5: Assign Dominant Community to Each Patient
#
# Input:
#   cox_dataset_baseline_scores.csv     (from Step 4)
#   profile_dataset_all_scores.csv      (from Step 4, optional)
#
# Output (written to <OUT_DIR>):
#   cox_dataset_with_dominant_community.csv
#   profile_dataset_with_dominant_community.csv  (if profile CSV exists)
#
# New columns added to each file (per score type: unw and w):
#   dom_community_{unw|w}  : "NONE" / "MIXED" / "DOM_C{k}"
#   dom_share_{unw|w}      : top-community share = top1 score / total score
#   total_score_{unw|w}    : sum of all community scores (UNROUNDED, full precision)
#
# Classification rules (TAU = 0.40):
#   total_score == 0        -> NONE      (no community diseases at baseline)
#   top1_share >= TAU       -> DOM_C{k}  (one community clearly dominates)
#   top1_share <  TAU       -> MIXED     (no single dominant community)
#
# Notes:
#   - This is the "merged MIXED" design used as the PRIMARY analysis (docx §2.4.3.5).
#   - MIXED_HIGH/LOW split variants explored in earlier drafts are NOT produced here.
#     Downstream Cox scripts (Step 5/5b/6/7/8) use the merged MIXED column only.
#
# Implementation notes:
#   - Score columns are sorted by community number (reproducible tie-breaking)
#   - total_score_{unw|w} saved UNROUNDED (full precision for Cox covariates)
#   - Defensive check: stops with clear error if no score columns detected
# ==============================================================================

rm(list = ls()); gc()
library(data.table)

# =========================
# Paths & threshold
# =========================
COX_CSV     <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\cox_dataset_baseline_scores.csv"
PROFILE_CSV <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\profile_dataset_all_scores.csv"
OUT_DIR     <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\cox_analysis\\dominant_community"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TAU <- 0.40   # dominance threshold (docx §2.4.3.5)

# =========================
# Helpers
# =========================
sort_score_cols <- function(cols) {
  # Sort "comm12_unw" etc. by community number for deterministic tie-breaking
  nums <- as.integer(gsub("comm(\\d+)_(unw|w)", "\\1", cols))
  cols[order(nums)]
}

assign_dominant <- function(dt, score_cols, suffix) {
  if (length(score_cols) == 0) {
    stop(sprintf("No score columns found for suffix '%s'. Expected columns like comm1_%s, comm2_%s, ...",
                 suffix, suffix, suffix), call. = FALSE)
  }
  
  score_cols <- sort_score_cols(score_cols)
  
  comm_nums <- as.integer(gsub("comm(\\d+)_(unw|w)", "\\1", score_cols))
  cat(sprintf("   [%s] Using %d score columns: comm%d..comm%d\n",
              suffix, length(score_cols), min(comm_nums), max(comm_nums)))
  
  S <- as.matrix(dt[, ..score_cols])
  S[is.na(S)] <- 0
  
  comm_labels <- paste0("C", comm_nums)
  
  row_sum    <- rowSums(S)
  row_max    <- apply(S, 1, max)
  row_which  <- max.col(S, ties.method = "first")
  top1_share <- ifelse(row_sum > 0, row_max / row_sum, 0)
  dom_label  <- comm_labels[row_which]
  
  group <- ifelse(
    row_sum == 0, "NONE",
    ifelse(top1_share >= TAU, paste0("DOM_", dom_label), "MIXED")
  )
  
  n_total <- length(group)
  n_none  <- sum(group == "NONE")
  n_dom   <- sum(grepl("^DOM_", group))
  n_mixed <- sum(group == "MIXED")
  cat(sprintf("   [%s] NONE=%d (%.2f%%) | DOM=%d (%.2f%%) | MIXED=%d (%.2f%%)\n",
              suffix,
              n_none,  100 * n_none  / n_total,
              n_dom,   100 * n_dom   / n_total,
              n_mixed, 100 * n_mixed / n_total))
  
  dom_tab <- sort(table(group[grepl("^DOM_", group)]), decreasing = TRUE)
  if (length(dom_tab) > 0) {
    cat(sprintf("   [%s] DOM breakdown: %s\n", suffix,
                paste(names(dom_tab), dom_tab, sep = "=", collapse = ", ")))
  }
  
  col_dom   <- paste0("dom_community_", suffix)
  col_share <- paste0("dom_share_", suffix)
  col_total <- paste0("total_score_", suffix)
  
  dt[, (col_dom)   := group]
  dt[, (col_share) := round(top1_share, 4)]
  dt[, (col_total) := row_sum]   # UNROUNDED
  
  dt
}

# =========================
# Process COX dataset
# =========================
cat("==> Processing cox_dataset_baseline_scores.csv...\n")
cox_dt <- fread(COX_CSV)
cat("   Rows:", nrow(cox_dt), "\n")

unw_cols <- grep("^comm\\d+_unw$", names(cox_dt), value = TRUE)
w_cols   <- grep("^comm\\d+_w$",   names(cox_dt), value = TRUE)
cat("   Detected", length(unw_cols), "unw cols,", length(w_cols), "w cols\n\n")

cox_dt <- assign_dominant(cox_dt, unw_cols, "unw")
cox_dt <- assign_dominant(cox_dt, w_cols,   "w")

cox_out <- file.path(OUT_DIR, "cox_dataset_with_dominant_community.csv")
fwrite(cox_dt, cox_out)
cat("\n   Saved:", cox_out, "\n")

# =========================
# Process PROFILE dataset (optional)
# =========================
if (!file.exists(PROFILE_CSV)) {
  cat("\n==> PROFILE_CSV not found, skipping profile dataset processing.\n")
  cat("    Expected:", PROFILE_CSV, "\n")
} else {
  cat("\n==> Processing profile_dataset_all_scores.csv...\n")
  prof_dt <- fread(PROFILE_CSV)
  cat("   Rows:", nrow(prof_dt), "\n")
  
  unw_cols_p <- grep("^comm\\d+_unw$", names(prof_dt), value = TRUE)
  w_cols_p   <- grep("^comm\\d+_w$",   names(prof_dt), value = TRUE)
  cat("   Detected", length(unw_cols_p), "unw cols,", length(w_cols_p), "w cols\n\n")
  
  prof_dt <- assign_dominant(prof_dt, unw_cols_p, "unw")
  prof_dt <- assign_dominant(prof_dt, w_cols_p,   "w")
  
  prof_out <- file.path(OUT_DIR, "profile_dataset_with_dominant_community.csv")
  fwrite(prof_dt, prof_out)
  cat("\n   Saved:", prof_out, "\n")
}

# =========================
# Summary
# =========================
cat("\n============================================\n")
cat("Step 4.5: Dominant Community Assignment complete.\n")
cat("============================================\n")
cat("\nColumns added (per file, per score type unw/w):\n")
cat("  dom_community_{unw|w}  : NONE / MIXED / DOM_C{k}\n")
cat("  dom_share_{unw|w}      : top1 score / total score\n")
cat("  total_score_{unw|w}    : sum of all scores (UNROUNDED)\n")
cat(sprintf("\nClassification rules (TAU = %.2f):\n", TAU))
cat("  total_score == 0        -> NONE\n")
cat("  top1_share >= TAU       -> DOM_C{k}\n")
cat("  top1_share <  TAU       -> MIXED\n")
cat("\nFiles saved:\n")
cat("  ", cox_out, "\n")
if (exists("prof_out")) cat("  ", prof_out, "\n")
cat("\nDownstream expected columns (Cox Step 5):\n")
cat("  dom_community_unw (main exposure)\n")
cat("  dom_community_w   (weighted sensitivity, S4)\n")
cat("  total_score_unw   (continuous covariate, Model 3)\n")
cat("============================================\n")