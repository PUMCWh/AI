#!/usr/bin/env Rscript
# ==============================================================================
# Leiden Resolution Sweep
# Purpose: Test multiple resolution values to find optimal community granularity
# for a ~500-node multimorbidity network.
#
# Outputs:
#   - Comparison table (Excel) with quality metrics per resolution
#   - Panel plot (6 key metrics vs resolution)
#   - Composite quality score plot
#   - Community size distribution boxplot
#   - Modularity vs n_communities elbow plot
#
# Usage: Adjust RDS_DIR / OUT_DIR, then source() or Rscript this file.
# ==============================================================================
rm(list = ls()); gc()

suppressPackageStartupMessages({
  pkgs <- c("data.table", "Matrix", "igraph", "ggplot2", "openxlsx", "scales")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  }
  library(data.table); library(Matrix); library(igraph)
  library(ggplot2); library(openxlsx); library(scales)
})

options(stringsAsFactors = FALSE)

# =========================
# 0) Paths & Params  *** MODIFY THESE ***
# =========================
RDS_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\resolution_sweep_fine"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Edge filter — use your primary threshold
COEF_THR  <- 0.01
ALPHA_THR <- 0.05

# Resolution grid: coarse full range + fine 1.0-2.0 supplement
RESOLUTION_GRID <- sort(unique(c(
  0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0,  # coarse
  1.1, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9                                           # fine supplement
)))

# Leiden params
LEIDEN_N_ITER <- 100
LEIDEN_SEED   <- 42
N_SEEDS_STABILITY <- 10   # per-resolution stability check (快速版, 正式分析用30)

FONT_FAMILY <- "sans"   # avoid Arial font issue on some Windows systems

# PDF device
PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) grDevices::pdf)

# =========================
# 1) Load network
# =========================
cat("==> Step 1) Loading network...\n")
rds_candidates <- list.dirs(RDS_DIR, recursive = FALSE, full.names = TRUE)
rds_candidates <- rds_candidates[grepl("NETWORK_PCOR_", basename(rds_candidates))]
if (length(rds_candidates) == 0) stop("No NETWORK_PCOR_* folder found in RDS_DIR")
rds_folder <- sort(rds_candidates, decreasing = TRUE)[1]
cat("   RDS folder:", rds_folder, "\n")

disease_map <- readRDS(file.path(rds_folder, "disease_map.rds"))
edge_dt     <- readRDS(file.path(rds_folder, "edges_all_uppertri_with_fdr.rds"))

code_by_id    <- disease_map$code
chapter_by_id <- as.character(disease_map$ICD10_Chapter)
chapter_by_id[is.na(chapter_by_id) | chapter_by_id == ""] <- "0"
prev_by_id    <- disease_map$prev; prev_by_id[is.na(prev_by_id)] <- 0
k <- length(code_by_id)
cat("   Total disease nodes: k =", k, "\n")

# Build graph
edges_sub <- edge_dt[p_adj < ALPHA_THR & weight >= COEF_THR, .(from, to, weight)]
vertex_df <- data.frame(name = code_by_id, chapter = chapter_by_id,
                        prevalence = prev_by_id, stringsAsFactors = FALSE)
g_full <- graph_from_data_frame(edges_sub, directed = FALSE, vertices = vertex_df)
n_connected <- sum(degree(g_full) > 0)
cat("   Edges:", ecount(g_full), "  Connected nodes:", n_connected,
    "  Isolates:", k - n_connected, "\n")

# =========================
# 2) Leiden function (identical to your main script)
# =========================
run_leiden <- function(g, resolution = 1.0, n_iter = 100, seed = 42) {
  set.seed(seed)
  w <- E(g)$weight
  g_pos <- if (any(w <= 0)) delete_edges(g, which(w <= 0)) else g
  iso_v <- which(degree(g_pos) == 0)
  node_names_all <- V(g_pos)$name
  g_sub <- if (length(iso_v) > 0) delete_vertices(g_pos, iso_v) else g_pos
  if (vcount(g_sub) == 0 || ecount(g_sub) == 0) {
    return(data.table(node = node_names_all, community = 0L))
  }
  clu <- cluster_leiden(g_sub, objective_function = "modularity",
                        resolution = resolution, n_iterations = n_iter,
                        weights = E(g_sub)$weight)
  result <- data.table(node = V(g_sub)$name, community = as.integer(membership(clu)))
  if (length(iso_v) > 0) {
    result <- rbind(result, data.table(node = node_names_all[iso_v], community = 0L))
  }
  result
}

# =========================
# 3) Quality metrics
# =========================

# Coverage: fraction of intra-community edge weight / total edge weight
calc_coverage <- function(g, comm_dt) {
  w <- E(g)$weight
  g_pos <- if (any(w <= 0)) delete_edges(g, which(w <= 0)) else g
  if (ecount(g_pos) == 0) return(NA_real_)
  
  comm_map <- setNames(comm_dt$community, comm_dt$node)
  el <- as.data.table(as_data_frame(g_pos, what = "edges"))
  el[, from_c := comm_map[from]]
  el[, to_c   := comm_map[to]]
  
  total_w <- sum(el$weight)
  intra_w <- sum(el$weight[el$from_c == el$to_c & el$from_c > 0 & el$to_c > 0])
  intra_w / total_w
}

# Average conductance (lower = better community separation)
calc_avg_conductance <- function(g, comm_dt) {
  w <- E(g)$weight
  g_pos <- if (any(w <= 0)) delete_edges(g, which(w <= 0)) else g
  if (ecount(g_pos) == 0) return(NA_real_)
  
  comm_map <- setNames(comm_dt$community, comm_dt$node)
  el <- as.data.table(as_data_frame(g_pos, what = "edges"))
  el[, from_c := comm_map[from]]
  el[, to_c   := comm_map[to]]
  
  comms <- sort(unique(comm_dt$community[comm_dt$community > 0]))
  if (length(comms) < 2) return(NA_real_)
  
  cond_vals <- numeric(length(comms))
  for (i in seq_along(comms)) {
    cid <- comms[i]
    nodes_in <- comm_dt$node[comm_dt$community == cid]
    
    # Edges touching this community (at least one endpoint inside)
    touching <- el[from %in% nodes_in | to %in% nodes_in]
    if (nrow(touching) == 0) { cond_vals[i] <- 0; next }
    
    # Volume: total weight of edges touching community
    vol_c <- sum(touching$weight)
    
    # Cut: weight of edges crossing community boundary
    cut_c <- sum(touching$weight[
      (touching$from %in% nodes_in & !touching$to %in% nodes_in) |
        (!touching$from %in% nodes_in & touching$to %in% nodes_in)
    ])
    
    total_vol <- sum(el$weight)
    denom <- min(vol_c, total_vol - vol_c)
    cond_vals[i] <- if (denom > 0) cut_c / denom else 0
  }
  mean(cond_vals)
}

# Mean internal edge density across communities
calc_mean_internal_density <- function(g, comm_dt) {
  comms <- sort(unique(comm_dt$community[comm_dt$community > 0]))
  if (length(comms) < 2) return(NA_real_)
  
  w <- E(g)$weight
  g_pos <- if (any(w <= 0)) delete_edges(g, which(w <= 0)) else g
  
  int_dens <- numeric(length(comms))
  for (i in seq_along(comms)) {
    nodes_in <- comm_dt$node[comm_dt$community == comms[i]]
    if (length(nodes_in) < 2) { int_dens[i] <- 0; next }
    g_sub <- induced_subgraph(g_pos, vids = which(V(g_pos)$name %in% nodes_in))
    max_edges <- choose(length(nodes_in), 2)
    int_dens[i] <- if (max_edges > 0) ecount(g_sub) / max_edges else 0
  }
  mean(int_dens, na.rm = TRUE)
}

# =========================
# 4) Sweep loop
# =========================
cat("\n==> Step 4) Running resolution sweep...\n")
cat(sprintf("   Testing %d resolution values: %s\n",
            length(RESOLUTION_GRID), paste(RESOLUTION_GRID, collapse = ", ")))

sweep_results <- list()
comm_size_details <- list()

t0 <- Sys.time()

for (ri in seq_along(RESOLUTION_GRID)) {
  res <- RESOLUTION_GRID[ri]
  cat(sprintf("\n   [%2d/%d] resolution = %5.1f ... ",
              ri, length(RESOLUTION_GRID), res))
  
  # Primary partition (seed = LEIDEN_SEED)
  comm_dt <- run_leiden(g_full, resolution = res, n_iter = LEIDEN_N_ITER, seed = LEIDEN_SEED)
  
  active <- comm_dt[community > 0]
  n_comm <- length(unique(active$community))
  n_iso  <- sum(comm_dt$community == 0)
  
  # Community sizes
  comm_sizes <- active[, .N, by = community]
  
  # Modularity
  m_vec <- setNames(comm_dt$community, comm_dt$node)
  mod_val <- tryCatch(
    modularity(g_full, membership = m_vec[V(g_full)$name], weights = E(g_full)$weight),
    error = function(e) NA_real_)
  
  # Quality metrics
  coverage <- calc_coverage(g_full, comm_dt)
  avg_cond <- calc_avg_conductance(g_full, comm_dt)
  int_dens <- calc_mean_internal_density(g_full, comm_dt)
  
  # Singleton & tiny communities
  n_singletons <- sum(comm_sizes$N == 1)
  n_tiny       <- sum(comm_sizes$N <= 3)
  
  # Multi-seed stability (quick, N_SEEDS_STABILITY seeds)
  ari_vals  <- numeric(N_SEEDS_STABILITY)
  ncomm_vals <- integer(N_SEEDS_STABILITY)
  for (s in 1:N_SEEDS_STABILITY) {
    test_dt <- run_leiden(g_full, resolution = res, n_iter = LEIDEN_N_ITER, seed = s)
    ncomm_vals[s] <- length(unique(test_dt[community > 0, community]))
    merged <- merge(comm_dt[community > 0], test_dt[community > 0],
                    by = "node", suffixes = c("_ref", "_test"))
    ari_vals[s] <- if (nrow(merged) > 0) {
      compare(merged$community_ref, merged$community_test, method = "adjusted.rand")
    } else NA_real_
  }
  
  row <- data.table(
    resolution          = res,
    n_communities       = n_comm,
    modularity          = round(mod_val, 4),
    coverage            = round(coverage, 4),
    avg_conductance     = round(avg_cond, 4),
    internal_density    = round(int_dens, 4),
    isolates            = n_iso,
    singleton_comms     = n_singletons,
    tiny_comms_le3      = n_tiny,
    avg_comm_size       = round(mean(comm_sizes$N), 1),
    median_comm_size    = median(comm_sizes$N),
    min_comm_size       = min(comm_sizes$N),
    max_comm_size       = max(comm_sizes$N),
    sd_comm_size        = round(sd(comm_sizes$N), 1),
    ARI_stability_mean  = round(mean(ari_vals, na.rm = TRUE), 4),
    ARI_stability_sd    = round(sd(ari_vals, na.rm = TRUE), 4),
    ARI_stability_min   = round(min(ari_vals, na.rm = TRUE), 4),
    n_comms_seeds_range = paste0(min(ncomm_vals), "-", max(ncomm_vals))
  )
  
  sweep_results[[ri]] <- row
  comm_size_details[[ri]] <- comm_sizes[, .(resolution = res, community, n_diseases = N)]
  
  cat(sprintf("%2d comms, mod=%.3f, cov=%.3f, cond=%.3f, ARI=%.3f",
              n_comm, mod_val, coverage, avg_cond, mean(ari_vals, na.rm = TRUE)))
}

sweep_dt <- rbindlist(sweep_results)
sizes_dt <- rbindlist(comm_size_details)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat(sprintf("\n\n==> Sweep complete. Time: %.1f minutes.\n\n", elapsed))

# Print summary table
cat("===== SWEEP RESULTS =====\n")
print(sweep_dt[, .(resolution, n_communities, modularity, coverage,
                   avg_conductance, internal_density,
                   ARI_stability_mean, avg_comm_size, singleton_comms)])

# =========================
# 5) Composite quality score
# =========================
# Normalize each metric to [0, 1], then weighted average.
# Higher is better:  modularity, coverage, internal_density, ARI_stability_mean
# Lower is better:   avg_conductance, singleton_comms  (inverted)

normalize_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] - rng[1] < 1e-12) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

sweep_dt[, norm_modularity   := normalize_01(modularity)]
sweep_dt[, norm_coverage     := normalize_01(coverage)]
sweep_dt[, norm_int_density  := normalize_01(internal_density)]
sweep_dt[, norm_stability    := normalize_01(ARI_stability_mean)]
sweep_dt[, norm_conductance  := 1 - normalize_01(avg_conductance)]   # inverted
sweep_dt[, norm_singletons   := 1 - normalize_01(singleton_comms)]   # inverted

# Composite (equal weight — adjust if you prefer different emphasis)
sweep_dt[, quality_score := round(
  (norm_modularity + norm_coverage + norm_int_density +
     norm_stability + norm_conductance + norm_singletons) / 6, 4)]

# Rank
sweep_dt[, rank := frank(-quality_score, ties.method = "min")]
setorder(sweep_dt, resolution)

best_res <- sweep_dt[which.max(quality_score), resolution]

cat(sprintf("\n==> BEST resolution by composite score: %.1f\n", best_res))
cat(sprintf("   Communities: %d | Modularity: %.3f | Avg size: %.1f | ARI: %.3f\n",
            sweep_dt[resolution == best_res, n_communities],
            sweep_dt[resolution == best_res, modularity],
            sweep_dt[resolution == best_res, avg_comm_size],
            sweep_dt[resolution == best_res, ARI_stability_mean]))

# =========================
# 6) Plots
# =========================
cat("\n==> Step 6) Generating plots...\n")

# --- 6a) Multi-panel: 6 key metrics vs resolution ---
plot_dt <- melt(sweep_dt,
                id.vars = "resolution",
                measure.vars = c("n_communities", "modularity", "coverage",
                                 "avg_conductance", "ARI_stability_mean", "avg_comm_size"),
                variable.name = "metric", value.name = "value")

metric_labels <- c(
  n_communities       = "Number of communities",
  modularity          = "Modularity",
  coverage            = "Coverage (intra-community weight fraction)",
  avg_conductance     = "Avg conductance (lower = better separation)",
  ARI_stability_mean  = "Partition stability (mean ARI, 10 seeds)",
  avg_comm_size       = "Avg community size (diseases)"
)
plot_dt[, metric_label := factor(metric_labels[as.character(metric)],
                                 levels = metric_labels)]

p_panel <- ggplot(plot_dt, aes(x = resolution, y = value)) +
  geom_line(linewidth = 0.8, color = "#4472C4") +
  geom_point(size = 2.5, color = "#4472C4") +
  geom_vline(xintercept = best_res, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = RESOLUTION_GRID[seq(1, length(RESOLUTION_GRID), by = 1)]) +
  labs(title = "Leiden resolution sweep: key metrics",
       subtitle = sprintf("Red dashed line = best composite (res = %.1f) | %d nodes, coef >= %.2f, FDR < %.2f",
                          best_res, k, COEF_THR, ALPHA_THR),
       x = "Resolution parameter", y = "") +
  theme_minimal(base_size = 11) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(OUT_DIR, "resolution_sweep_panel.pdf"),
       p_panel, width = 13, height = 10, units = "in", device = PDF_DEVICE)
cat("   Panel plot saved.\n")

# --- 6b) Composite quality score ---
p_quality <- ggplot(sweep_dt, aes(x = resolution, y = quality_score)) +
  geom_line(linewidth = 1, color = "#ED7D31") +
  geom_point(size = 3, color = "#ED7D31") +
  geom_point(data = sweep_dt[resolution == best_res],
             size = 6, color = "red", shape = 18) +
  geom_text(data = sweep_dt[resolution == best_res],
            aes(label = sprintf("Best: res=%.1f\n%d communities", resolution, n_communities)),
            vjust = -1.5, hjust = 0.5, size = 3.5, color = "red", fontface = "bold",
            family = FONT_FAMILY) +
  scale_x_continuous(breaks = RESOLUTION_GRID) +
  labs(title = "Composite quality score by resolution",
       subtitle = "Equal weight: modularity + coverage + internal density + stability - conductance - singletons",
       x = "Resolution parameter", y = "Composite quality score [0-1]") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "resolution_quality_score.pdf"),
       p_quality, width = 10, height = 6, units = "in", device = PDF_DEVICE)
cat("   Quality score plot saved.\n")

# --- 6c) Community size distribution boxplot per resolution ---
sizes_dt[, resolution_f := factor(sprintf("%.1f", resolution),
                                  levels = sprintf("%.1f", RESOLUTION_GRID))]
p_sizes <- ggplot(sizes_dt, aes(x = resolution_f, y = n_diseases)) +
  geom_boxplot(fill = "#4472C4", alpha = 0.4, outlier.size = 1.5) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1, color = "#4472C4") +
  geom_hline(yintercept = c(5, 10, 20, 50), linetype = "dotted", color = "grey60") +
  annotate("text", x = 0.5, y = c(5, 10, 20, 50),
           label = c("5", "10", "20", "50"),
           hjust = -0.2, color = "grey50", size = 3, family = FONT_FAMILY) +
  labs(title = "Community size distribution by resolution",
       subtitle = sprintf("%d disease nodes | coef >= %.2f, FDR < %.2f", k, COEF_THR, ALPHA_THR),
       x = "Resolution", y = "Community size (number of diseases)") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

ggsave(file.path(OUT_DIR, "community_size_distribution.pdf"),
       p_sizes, width = 13, height = 7, units = "in", device = PDF_DEVICE)
cat("   Size distribution plot saved.\n")

# --- 6d) Modularity vs n_communities (elbow plot) ---
p_elbow <- ggplot(sweep_dt, aes(x = n_communities, y = modularity)) +
  geom_line(linewidth = 0.8, color = "grey50") +
  geom_point(aes(color = resolution), size = 4) +
  geom_text(aes(label = sprintf("r=%.1f", resolution)),
            vjust = -1.2, size = 2.8, family = FONT_FAMILY) +
  scale_color_viridis_c(name = "Resolution", option = "C") +
  labs(title = "Modularity vs number of communities",
       subtitle = "Look for the elbow: modularity plateau before dropping off",
       x = "Number of communities", y = "Modularity") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

ggsave(file.path(OUT_DIR, "modularity_vs_n_communities.pdf"),
       p_elbow, width = 10, height = 7, units = "in", device = PDF_DEVICE)
cat("   Elbow plot saved.\n")

# --- 6e) Stability vs resolution ---
p_stab <- ggplot(sweep_dt, aes(x = resolution, y = ARI_stability_mean)) +
  geom_ribbon(aes(ymin = ARI_stability_mean - ARI_stability_sd,
                  ymax = pmin(ARI_stability_mean + ARI_stability_sd, 1)),
              fill = "#4472C4", alpha = 0.2) +
  geom_line(linewidth = 0.8, color = "#4472C4") +
  geom_point(size = 2.5, color = "#4472C4") +
  geom_point(aes(y = ARI_stability_min), size = 1.5, color = "red", shape = 4) +
  geom_vline(xintercept = best_res, linetype = "dashed", color = "red", linewidth = 0.4) +
  scale_x_continuous(breaks = RESOLUTION_GRID) +
  labs(title = "Partition stability across resolutions",
       subtitle = sprintf("Mean ARI ± SD vs reference seed (10 seeds) | red × = min ARI"),
       x = "Resolution", y = "ARI (adjusted Rand index)") +
  ylim(0, 1.05) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "stability_vs_resolution.pdf"),
       p_stab, width = 10, height = 6, units = "in", device = PDF_DEVICE)
cat("   Stability plot saved.\n")

# =========================
# 7) Excel output
# =========================
cat("\n==> Step 7) Writing Excel...\n")
wb <- createWorkbook()
headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                           valign = "center", fgFill = "#4472C4", textDecoration = "bold")

add_sheet <- function(sn, dt_in) {
  sn <- substr(sn, 1, 31)
  addWorksheet(wb, sn)
  writeDataTable(wb, sn, dt_in, withFilter = TRUE)
  setColWidths(wb, sn, cols = 1:ncol(dt_in), widths = "auto")
  addStyle(wb, sn, headerStyle, rows = 1, cols = 1:ncol(dt_in), gridExpand = TRUE)
}

# Main comparison (sorted by quality)
sweep_ranked <- copy(sweep_dt)
setorder(sweep_ranked, -quality_score)
add_sheet("sweep_ranked", sweep_ranked)

# Same data sorted by resolution
sweep_by_res <- copy(sweep_dt)
setorder(sweep_by_res, resolution)
add_sheet("sweep_by_resolution", sweep_by_res)

# Community sizes per resolution
add_sheet("community_sizes", sizes_dt)

# Recommendation
rec_dt <- data.table(
  item = c("Network info",
           "Total disease nodes (k)",
           "Edges (coef>=0.01, FDR<0.05)",
           "Connected nodes",
           "",
           "Best resolution (composite score)",
           "Communities at best resolution",
           "Modularity at best resolution",
           "Avg community size at best resolution",
           "ARI stability at best resolution",
           "",
           "Current resolution (your main script)",
           "Communities at current resolution",
           "",
           "Interpretation guide",
           "Modularity",
           "Coverage",
           "Avg conductance",
           "Internal density",
           "ARI stability",
           "Singletons",
           "",
           "Next steps"),
  value = c("",
            as.character(k),
            as.character(ecount(g_full)),
            as.character(n_connected),
            "",
            as.character(best_res),
            as.character(sweep_dt[resolution == best_res, n_communities]),
            as.character(sweep_dt[resolution == best_res, modularity]),
            as.character(sweep_dt[resolution == best_res, avg_comm_size]),
            as.character(sweep_dt[resolution == best_res, ARI_stability_mean]),
            "",
            "1.0",
            as.character(sweep_dt[resolution == 1.0, n_communities]),
            "",
            "",
            "Higher = better defined communities (typically 0.3-0.7)",
            "Higher = more edge weight captured within communities",
            "Lower = better separation between communities",
            "Higher = denser internal connections",
            "Higher = more reproducible across random seeds (>0.8 ideal)",
            "Fewer = better (singleton communities are often noise)",
            "",
            paste0("1) Review plots. 2) Choose resolution balancing quality + interpretability. ",
                   "3) Set LEIDEN_RESOLUTION <- X in main script. 4) Re-run full analysis."))
)
add_sheet("recommendation", rec_dt)

xlsx_file <- file.path(OUT_DIR, "resolution_sweep_results.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", xlsx_file, "\n")

# =========================
# Final Summary
# =========================
cat("\n============================================================\n")
cat("  RESOLUTION SWEEP COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", OUT_DIR, "\n\n")
cat("Files:\n")
cat("  resolution_sweep_panel.pdf          6 metrics vs resolution\n")
cat("  resolution_quality_score.pdf        Composite score curve\n")
cat("  community_size_distribution.pdf     Size boxplots\n")
cat("  modularity_vs_n_communities.pdf     Elbow plot\n")
cat("  stability_vs_resolution.pdf         ARI stability ribbon\n")
cat("  resolution_sweep_results.xlsx       Full comparison table\n")
cat("\nTop 5 resolutions by composite quality:\n")
top5 <- copy(sweep_dt)
setorder(top5, -quality_score)
print(top5[1:min(5, .N), .(rank, resolution, n_communities,
                           modularity, coverage, ARI_stability_mean,
                           avg_comm_size, quality_score)])
cat(sprintf("\n==> RECOMMENDATION: Set LEIDEN_RESOLUTION <- %.1f\n", best_res))
cat("    Then re-run your main community detection script.\n")
cat("============================================================\n")