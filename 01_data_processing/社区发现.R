# ==============================================================================
# Leiden Community Detection on Multimorbidity Partial-Correlation Network
# Enhanced version with:
#   - Connection profiles (per-node cross-community strength shares)
#   - Community-community connectivity matrix + heatmap
#   - Enhanced community summary (Top10 by prevalence & by within-strength)
#   - Multi-seed stability (ARI/NMI across 30 seeds)
#   - Bridge node metrics (Guimerà & Amaral node roles)
# ==============================================================================
rm(list = ls()); gc()

suppressPackageStartupMessages({
  pkgs <- c("data.table", "Matrix", "igraph", "ggplot2", "ggraph", "grid",
            "openxlsx", "RColorBrewer")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  }
  library(data.table); library(Matrix); library(igraph)
  library(ggplot2); library(ggraph); library(grid)
  library(openxlsx); library(RColorBrewer)
})

options(stringsAsFactors = FALSE)

# =========================
# 0) Paths & Params
# =========================
RDS_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL"
OUT_DIR <- "C:\\Users\\HP\\Desktop\\paper\\CKM_FINAL\\community_detection"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

THRESH_GRID <- data.table(
  coef  = c(0.01),
  alpha = c(0.05)
)
# 注: 其他阈值(0.05/0.10)可作为敏感性分析, 需要时取消下面的注释
# THRESH_GRID <- data.table(
#   coef  = c(0.01, 0.05, 0.10),
#   alpha = c(0.05, 0.05, 0.05)
# )

LEIDEN_RESOLUTION <- 1.2
LEIDEN_N_ITER     <- 100
LEIDEN_SEED       <- 42
N_SEEDS_STABILITY <- 30   # multi-seed stability test

PDF_W <- 14; PDF_H <- 12
FONT_FAMILY <- "Arial"
NODE_LABEL_SIZE <- 2.2
SHOW_NODE_LABELS <- TRUE

chapter_colors_provided <- c(
  "1" = "#B70031", "2" = "#EC6A3D", "3" = "#A1C6AD", "4" = "#0098C5",
  "5" = "#9A6CAC", "6" = "#004B9B", "7" = "#C0C7CB", "8" = "#AAD9EA",
  "9" = "#FFE100", "10" = "#E71A10", "11" = "#6EB327", "12" = "#8E006C",
  "13" = "#D2ACC9", "14" = "#ADB5DA", "0" = "grey"
)

# PDF device: prefer cairo_pdf (better Unicode/font support), fall back to pdf
PDF_DEVICE <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 1, height = 1); dev.off(); file.remove(tmp)
  grDevices::cairo_pdf
}, error = function(e) {
  cat("   Note: cairo_pdf unavailable, using standard pdf device.\n")
  grDevices::pdf
})

# =========================
# 1) Load RDS intermediates
# =========================
cat("==> Step 1) Loading RDS intermediates...\n")
rds_candidates <- list.dirs(RDS_DIR, recursive = FALSE, full.names = TRUE)
rds_candidates <- rds_candidates[grepl("NETWORK_PCOR_", basename(rds_candidates))]
if (length(rds_candidates) > 0) {
  rds_folder <- sort(rds_candidates, decreasing = TRUE)[1]
  cat("   Auto-detected RDS folder:", rds_folder, "\n")
} else {
  rds_folder <- RDS_DIR
}

disease_map <- readRDS(file.path(rds_folder, "disease_map.rds"))
edge_dt     <- readRDS(file.path(rds_folder, "edges_all_uppertri_with_fdr.rds"))
pcor_mat    <- readRDS(file.path(rds_folder, "pcor_matrix.rds"))

code_by_id    <- disease_map$code
chapter_by_id <- as.character(disease_map$ICD10_Chapter)
chapter_by_id[is.na(chapter_by_id) | chapter_by_id == ""] <- "0"
prev_by_id    <- disease_map$prev
prev_by_id[is.na(prev_by_id)] <- 0
freq_by_id    <- disease_map$freq
freq_by_id[is.na(freq_by_id)] <- 0
k             <- length(code_by_id)
cat("   Loaded:", k, "diseases,", nrow(edge_dt), "pairs\n")

# =========================
# 2) Leiden community detection function
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
# 3) Community summary (enhanced: Top10 by prevalence & by within-strength)
# =========================
summarize_communities <- function(comm_dt, g, disease_map_dt, bridge_dt = NULL) {
  merged <- merge(comm_dt, disease_map_dt,
                  by.x = "node", by.y = "code", all.x = TRUE)
  
  comm_summary <- merged[community > 0, .(
    n_diseases       = .N,
    diseases         = paste(sort(node), collapse = "; "),
    chapters         = paste(sort(unique(ICD10_Chapter)), collapse = ", "),
    dominant_chapter = { tab <- table(ICD10_Chapter); names(tab)[which.max(tab)] },
    mean_prevalence  = mean(prev, na.rm = TRUE),
    max_prevalence   = max(prev, na.rm = TRUE),
    total_freq       = sum(freq, na.rm = TRUE)
  ), by = community]
  setorder(comm_summary, community)
  
  # Top10 by prevalence
  top10_prev <- merged[community > 0, {
    dt_ord <- .SD[order(-prev)]
    .(top10_by_prevalence = paste(head(dt_ord$node, 10), collapse = "; "))
  }, by = community]
  
  comm_summary <- merge(comm_summary, top10_prev, by = "community", all.x = TRUE)
  
  # Top10 by within-community strength (if bridge_dt available)
  if (!is.null(bridge_dt)) {
    br_s <- bridge_dt[community > 0, .(node, community, strength_within_module)]
    top10_str <- br_s[, {
      dt_ord <- .SD[order(-strength_within_module)]
      .(top10_by_within_strength = paste(head(dt_ord$node, 10), collapse = "; "))
    }, by = community]
    comm_summary <- merge(comm_summary, top10_str, by = "community", all.x = TRUE)
  }
  
  # Internal density & mean weight
  for (cid in comm_summary$community) {
    nodes_in <- comm_dt[community == cid, node]
    if (length(nodes_in) < 2) {
      comm_summary[community == cid, `:=`(internal_density = NA_real_,
                                          internal_mean_weight = NA_real_)]
      next
    }
    g_sub <- induced_subgraph(g, vids = which(V(g)$name %in% nodes_in))
    d <- if (ecount(g_sub) == 0) 0 else edge_density(g_sub, loops = FALSE)
    mw <- if (ecount(g_sub) == 0) NA_real_ else mean(E(g_sub)$weight, na.rm = TRUE)
    comm_summary[community == cid, `:=`(internal_density = d,
                                        internal_mean_weight = mw)]
  }
  
  comm_summary
}

# =========================
# 4a) Bridge metrics + connection profile
# =========================
calc_bridge_metrics <- function(g, comm_dt) {
  
  w <- E(g)$weight
  g_pos <- if (any(w <= 0)) delete_edges(g, which(w <= 0)) else g
  
  comm_map <- setNames(comm_dt$community, comm_dt$node)
  membership_vec <- comm_map[V(g_pos)$name]
  membership_vec[is.na(membership_vec)] <- 0L
  
  node_names <- V(g_pos)$name
  n_v <- length(node_names)
  adj <- as_adjacency_matrix(g_pos, attr = "weight", sparse = TRUE)
  
  # Strength
  s_i <- strength(g_pos, weights = E(g_pos)$weight)
  
  # Within-module strength
  all_comms <- sort(unique(membership_vec[membership_vec > 0]))
  s_within <- numeric(n_v)
  for (v in 1:n_v) {
    my_comm <- membership_vec[v]
    if (my_comm == 0) next
    nb <- which(adj[v, ] > 0)
    s_within[v] <- sum(adj[v, nb[membership_vec[nb] == my_comm]])
  }
  
  # Within-module z-score
  z_score <- numeric(n_v)
  for (cid in all_comms) {
    idx_c <- which(membership_vec == cid)
    sw_c <- s_within[idx_c]
    mu_c <- mean(sw_c); sd_c <- sd(sw_c)
    if (is.na(sd_c) || sd_c < 1e-12) { z_score[idx_c] <- 0 }
    else { z_score[idx_c] <- (sw_c - mu_c) / sd_c }
  }
  
  # Participation coefficient + connection profile per node
  P_coef <- numeric(n_v)
  # Store per-node strength to each community
  s_to_comm <- matrix(0, nrow = n_v, ncol = length(all_comms))
  colnames(s_to_comm) <- as.character(all_comms)
  
  for (v in 1:n_v) {
    if (s_i[v] < 1e-12 || membership_vec[v] == 0) { P_coef[v] <- NA_real_; next }
    nb <- which(adj[v, ] > 0)
    if (length(nb) == 0) { P_coef[v] <- 0; next }
    
    sum_sq <- 0
    for (ci in seq_along(all_comms)) {
      cid <- all_comms[ci]
      in_c <- nb[membership_vec[nb] == cid]
      s_ic <- sum(adj[v, in_c])
      s_to_comm[v, ci] <- s_ic
      sum_sq <- sum_sq + (s_ic / s_i[v])^2
    }
    P_coef[v] <- 1 - sum_sq
  }
  
  # External strength
  external_strength <- numeric(n_v)
  for (v in 1:n_v) {
    if (membership_vec[v] == 0) { external_strength[v] <- NA_real_; next }
    own_idx <- which(all_comms == membership_vec[v])
    external_strength[v] <- s_i[v] - s_to_comm[v, own_idx]
  }
  
  # Top 3 EXTERNAL target communities per node (exclude own community)
  # share_total    = strength to target community / total strength
  # share_external = strength to target community / external strength only
  top1_comm <- rep(NA_integer_, n_v); top1_share_total <- rep(NA_real_, n_v); top1_share_external <- rep(NA_real_, n_v)
  top2_comm <- rep(NA_integer_, n_v); top2_share_total <- rep(NA_real_, n_v); top2_share_external <- rep(NA_real_, n_v)
  top3_comm <- rep(NA_integer_, n_v); top3_share_total <- rep(NA_real_, n_v); top3_share_external <- rep(NA_real_, n_v)
  
  for (v in 1:n_v) {
    my_comm <- membership_vec[v]
    if (s_i[v] < 1e-12 || my_comm == 0) next
    
    own_idx <- which(all_comms == my_comm)
    ext_str <- s_i[v] - s_to_comm[v, own_idx]
    
    # KEY FIX: if no external strength, keep all top* as NA (no fake "targets")
    if (!is.finite(ext_str) || ext_str <= 1e-12) next
    
    shares_t <- s_to_comm[v, ] / s_i[v]
    shares_e <- s_to_comm[v, ] / ext_str
    
    # Mask own community so it's never selected as a "target"
    shares_t[own_idx] <- -Inf
    shares_e[own_idx] <- -Inf
    # Mask zero-contribution communities (avoid fake "targets")
    shares_t[s_to_comm[v, ] <= 1e-12] <- -Inf
    shares_e[s_to_comm[v, ] <= 1e-12] <- -Inf
    
    ord <- order(-shares_t)
    if (length(ord) >= 1 && is.finite(shares_t[ord[1]])) {
      top1_comm[v] <- all_comms[ord[1]]
      top1_share_total[v] <- shares_t[ord[1]]
      top1_share_external[v] <- shares_e[ord[1]]
    }
    if (length(ord) >= 2 && is.finite(shares_t[ord[2]])) {
      top2_comm[v] <- all_comms[ord[2]]
      top2_share_total[v] <- shares_t[ord[2]]
      top2_share_external[v] <- shares_e[ord[2]]
    }
    if (length(ord) >= 3 && is.finite(shares_t[ord[3]])) {
      top3_comm[v] <- all_comms[ord[3]]
      top3_share_total[v] <- shares_t[ord[3]]
      top3_share_external[v] <- shares_e[ord[3]]
    }
  }
  
  result <- data.table(
    node                   = node_names,
    community              = as.integer(membership_vec),
    strength               = as.numeric(s_i),
    strength_within_module = as.numeric(s_within),
    external_strength      = as.numeric(external_strength),
    within_module_z_score  = as.numeric(z_score),
    participation_coef     = as.numeric(P_coef),
    top1_ext_community     = top1_comm,
    top1_share_total       = round(top1_share_total, 4),
    top1_share_external    = round(top1_share_external, 4),
    top2_ext_community     = top2_comm,
    top2_share_total       = round(top2_share_total, 4),
    top2_share_external    = round(top2_share_external, 4),
    top3_ext_community     = top3_comm,
    top3_share_total       = round(top3_share_total, 4),
    top3_share_external    = round(top3_share_external, 4)
  )
  
  # Node role classification
  result[, role := fcase(
    community == 0,                                            "isolate",
    within_module_z_score >= 2.5 & participation_coef >= 0.30, "connector_hub",
    within_module_z_score >= 2.5 & participation_coef <  0.30, "provincial_hub",
    within_module_z_score <  2.5 & participation_coef >= 0.62, "kinless",
    within_module_z_score <  2.5 & participation_coef >= 0.30, "connector",
    within_module_z_score <  2.5 & participation_coef <  0.30, "peripheral",
    default = "unclassified"
  )]
  
  # Return both the node table and the full s_to_comm matrix (for community matrix)
  list(node_dt = result, s_to_comm = s_to_comm, all_comms = all_comms,
       membership_vec = membership_vec)
}

# =========================
# 4b) Community-community connectivity matrix
# =========================
calc_community_matrix <- function(g, comm_dt) {
  w <- E(g)$weight
  g_pos <- if (any(w <= 0)) delete_edges(g, which(w <= 0)) else g
  
  comm_map <- setNames(comm_dt$community, comm_dt$node)
  el <- as.data.table(as_data_frame(g_pos, what = "edges"))
  el[, from_c := comm_map[from]]
  el[, to_c   := comm_map[to]]
  
  # Only consider edges between active (non-isolate) nodes
  el <- el[from_c > 0 & to_c > 0]
  
  # Use ALL communities from comm_dt (not just those appearing in edges)
  all_comms <- sort(unique(comm_dt$community[comm_dt$community > 0]))
  n_c <- length(all_comms)
  comm_labels <- as.character(all_comms)
  mat <- matrix(0, n_c, n_c, dimnames = list(comm_labels, comm_labels))
  
  # Vectorised aggregation: canonical pair (min, max) then sum
  if (nrow(el) > 0) {
    el[, ca := pmin(from_c, to_c)]
    el[, cb := pmax(from_c, to_c)]
    pair_sum <- el[, .(total = sum(weight)), by = .(ca, cb)]
    for (r in 1:nrow(pair_sum)) {
      a <- as.character(pair_sum$ca[r])
      b <- as.character(pair_sum$cb[r])
      mat[a, b] <- mat[a, b] + pair_sum$total[r]
      if (a != b) mat[b, a] <- mat[b, a] + pair_sum$total[r]
    }
  }
  
  # Top contributing bridge nodes per inter-community pair
  inter <- el[from_c != to_c]
  bridge_contrib <- NULL
  if (nrow(inter) > 0) {
    inter[, pair := paste(pmin(from_c, to_c), pmax(from_c, to_c), sep = "-")]
    contrib_from <- inter[, .(node = from, pair, weight)]
    contrib_to   <- inter[, .(node = to,   pair, weight)]
    contrib_all  <- rbind(contrib_from, contrib_to)
    contrib_agg  <- contrib_all[, .(total_weight = sum(weight)), by = .(pair, node)]
    bridge_contrib <- contrib_agg[, {
      dt_ord <- .SD[order(-total_weight)]
      head(dt_ord, 5)
    }, by = pair]
    bridge_contrib[, rank := seq_len(.N), by = pair]
  }
  
  list(matrix = mat, bridge_contrib = bridge_contrib)
}

# =========================
# 4c) Plot: community hull network (same as before)
# =========================
plot_community_hull <- function(g, comm_dt, title_text) {
  base_theme <- theme_void(base_size = 10) +
    theme(
      text = element_text(family = FONT_FAMILY),
      plot.title = element_text(hjust = 0.5, family = FONT_FAMILY, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, family = FONT_FAMILY, size = 10,
                                   color = "grey40", margin = margin(b = 8)),
      legend.title = element_text(size = 11, family = FONT_FAMILY),
      legend.text  = element_text(size = 10, family = FONT_FAMILY),
      legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.5, "cm"),
      legend.spacing.y = unit(0.35, "cm"), plot.margin = margin(10, 10, 10, 10)
    )
  
  if (is.null(g) || vcount(g) == 0 || ecount(g) == 0) {
    return(ggplot() + base_theme +
             geom_text(aes(0, 0, label = "No edges under this threshold"),
                       family = FONT_FAMILY, size = 6) + labs(title = title_text))
  }
  
  iso_v <- which(degree(g) == 0)
  g_plot <- if (length(iso_v) > 0) delete_vertices(g, iso_v) else g
  if (vcount(g_plot) == 0 || ecount(g_plot) == 0) {
    return(ggplot() + base_theme +
             geom_text(aes(0, 0, label = "Only isolates (no edges)"),
                       family = FONT_FAMILY, size = 6) + labs(title = title_text))
  }
  
  comm_map <- setNames(comm_dt$community, comm_dt$node)
  V(g_plot)$community <- as.integer(comm_map[V(g_plot)$name])
  V(g_plot)$community[is.na(V(g_plot)$community)] <- 0L
  
  set.seed(LEIDEN_SEED)
  lay <- create_layout(g_plot, layout = "fr")
  
  hull_data <- data.table(x = lay$x, y = lay$y, community = lay$community)
  hull_data <- hull_data[community > 0]
  
  expand_hull <- function(x, y, factor = 0.08) {
    cx <- mean(x); cy <- mean(y)
    list(x = cx + (x - cx) * (1 + factor), y = cy + (y - cy) * (1 + factor))
  }
  
  hull_polys <- hull_data[, {
    if (.N >= 3) {
      idx <- chull(x, y); idx <- c(idx, idx[1])
      ex <- expand_hull(x[idx], y[idx]); .(hx = ex$x, hy = ex$y)
    } else if (.N == 2) {
      dx <- diff(x); dy <- diff(y)
      len <- max(sqrt(dx^2 + dy^2), 1e-12)
      nx <- -dy/len*0.18; ny <- dx/len*0.18
      .(hx = c(x[1]-dx*0.1+nx, x[2]+dx*0.1+nx, x[2]+dx*0.1-nx, x[1]-dx*0.1-nx, x[1]-dx*0.1+nx),
        hy = c(y[1]-dy*0.1+ny, y[2]+dy*0.1+ny, y[2]+dy*0.1-ny, y[1]-dy*0.1-ny, y[1]-dy*0.1+ny))
    } else NULL
  }, by = community]
  
  label_data <- hull_data[, .(lx = mean(x), ly = max(y)), by = community]
  label_data[, label := paste0("C", community)]
  y_range <- diff(range(lay$y, na.rm = TRUE))
  label_data[, ly := ly + y_range * 0.035]
  
  n_comm <- length(unique(hull_data$community))
  subtitle_text <- sprintf("%d communities | %d nodes | %d edges",
                           n_comm, vcount(g_plot), ecount(g_plot))
  
  p <- ggraph(lay)
  if (nrow(hull_polys) > 0) {
    p <- p + geom_polygon(data = hull_polys, aes(x = hx, y = hy, group = community),
                          fill = NA, color = "grey50", linewidth = 0.75, linetype = "dashed")
  }
  if (nrow(label_data) > 0) {
    p <- p + geom_text(data = label_data, aes(x = lx, y = ly, label = label),
                       family = FONT_FAMILY, fontface = "bold", color = "grey35", size = 3.8)
  }
  p <- p +
    geom_edge_link(aes(edge_width = abs(weight)), alpha = 0.6, color = "grey65") +
    geom_node_point(aes(color = chapter, size = prevalence), alpha = 0.85)
  if (SHOW_NODE_LABELS) {
    p <- p + geom_node_text(aes(label = name), family = FONT_FAMILY,
                            color = "black", size = NODE_LABEL_SIZE)
  }
  p + scale_edge_width(range = c(0.1, 2.5),
                       name = expression("Partial correlation (" * varphi * ")"),
                       guide = guide_legend(order = 1, override.aes = list(color = "grey50", linetype = "solid", shape = NA))) +
    scale_color_manual(name = "ICD-10 Chapter", values = chapter_colors_provided, breaks = as.character(1:14),
                       guide = guide_legend(order = 2, override.aes = list(size = 6, alpha = 1),
                                            keywidth = unit(0.5,"cm"), keyheight = unit(0.5,"cm"))) +
    scale_size_continuous(name = "Prevalence", range = c(5, 12), guide = "none") +
    base_theme + labs(title = title_text, subtitle = subtitle_text)
}

# =========================
# 4d) Plot: individual community subnetwork
# =========================
plot_single_community <- function(g, comm_dt, community_id, coef_thr, alpha_thr) {
  nodes_in <- comm_dt[community == community_id, node]
  g_sub <- induced_subgraph(g, vids = which(V(g)$name %in% nodes_in))
  iso_v <- which(degree(g_sub) == 0)
  if (length(iso_v) > 0 && length(iso_v) < vcount(g_sub)) g_sub <- delete_vertices(g_sub, iso_v)
  
  base_theme <- theme_void(base_size = 10) +
    theme(text = element_text(family = FONT_FAMILY),
          plot.title = element_text(hjust = 0.5, family = FONT_FAMILY, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, family = FONT_FAMILY, size = 10, color = "grey40"),
          legend.title = element_text(size = 11, family = FONT_FAMILY),
          legend.text = element_text(size = 10, family = FONT_FAMILY),
          plot.margin = margin(10, 10, 10, 10))
  
  if (vcount(g_sub) == 0 || ecount(g_sub) == 0) {
    return(ggplot() + base_theme +
             geom_text(aes(0, 0, label = "Single node or no internal edges"),
                       family = FONT_FAMILY, size = 5) + labs(title = sprintf("Class %d", community_id)))
  }
  
  set.seed(LEIDEN_SEED)
  ggraph(g_sub, layout = "fr") +
    geom_edge_link(aes(edge_width = abs(weight)), alpha = 0.7, color = "grey60") +
    geom_node_point(aes(color = chapter, size = prevalence), alpha = 0.85) +
    geom_node_text(aes(label = name), family = FONT_FAMILY, color = "black", size = NODE_LABEL_SIZE) +
    scale_edge_width(range = c(0.3, 3.0), limits = c(0, 0.5), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
                     name = expression("Partial correlation (" * varphi * ")"),
                     guide = guide_legend(order = 1, override.aes = list(color = "grey50", linetype = "solid", shape = NA))) +
    scale_color_manual(name = "ICD-10 Chapter", values = chapter_colors_provided,
                       breaks = as.character(1:14), drop = TRUE,
                       guide = guide_legend(order = 2, override.aes = list(size = 6, alpha = 1),
                                            keywidth = unit(0.5,"cm"), keyheight = unit(0.5,"cm"))) +
    scale_size_continuous(name = "Prevalence", range = c(6, 14), guide = "none") +
    base_theme +
    labs(title = sprintf("Class %d  |  coef >= %.2f, FDR < %.2f", community_id, coef_thr, alpha_thr),
         subtitle = sprintf("%d diseases | %d edges", vcount(g_sub), ecount(g_sub)))
}

# =========================
# 4e) Community-community heatmap
# =========================
plot_community_heatmap <- function(comm_matrix, title_text) {
  mat <- comm_matrix
  # Sort community IDs numerically
  comm_ids_num <- sort(as.integer(rownames(mat)))
  comm_ids <- as.character(comm_ids_num)
  
  # Full symmetric matrix (not just upper triangle)
  dt_long <- as.data.table(as.table(mat[comm_ids, comm_ids]))
  setnames(dt_long, c("from_community", "to_community", "weight"))
  dt_long[, weight := as.numeric(weight)]
  
  # Off-diagonal, positive weight only
  dt_inter <- dt_long[from_community != to_community & weight > 0]
  
  # Force numeric factor ordering
  dt_inter[, from_community := factor(from_community, levels = comm_ids)]
  dt_inter[, to_community   := factor(to_community,   levels = comm_ids)]
  
  show_labels <- (length(comm_ids) <= 25)
  
  ggplot(dt_inter, aes(x = from_community, y = to_community, fill = weight)) +
    geom_tile(color = "white", linewidth = 0.5) +
    { if (show_labels) geom_text(aes(label = round(weight, 2)), size = 3, family = FONT_FAMILY) } +
    scale_fill_gradient(low = "white", high = "#B70031", name = "Total\nweight") +
    labs(title = title_text,
         x = "Community", y = "Community") +
    theme_minimal(base_size = 11) +
    theme(text = element_text(family = FONT_FAMILY),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 0),
          panel.grid = element_blank()) +
    coord_fixed()
}

# =========================
# 5) Main loop
# =========================
cat("==> Step 5) Running Leiden community detection...\n")

vertex_df_full <- data.frame(
  name = code_by_id,
  chapter = factor(chapter_by_id, levels = c(as.character(1:14), "0")),
  prevalence = prev_by_id, stringsAsFactors = FALSE
)
disease_map_dt <- data.table(code = code_by_id, ICD10_Chapter = chapter_by_id,
                             freq = freq_by_id, prev = prev_by_id)

all_membership <- list(); all_comm_summary <- list()
all_metrics <- list(); all_bridge <- list()
all_comm_matrix <- list(); all_bridge_contrib <- list()

for (idx in 1:nrow(THRESH_GRID)) {
  coef_thr  <- THRESH_GRID[idx, coef]
  alpha_thr <- THRESH_GRID[idx, alpha]
  tag <- sprintf("coef%.2f_alpha%.2f", coef_thr, alpha_thr)
  cat(sprintf("\n--- Threshold: coef >= %.2f, FDR < %.2f ---\n", coef_thr, alpha_thr))
  
  edges_sub <- edge_dt[p_adj < alpha_thr & weight >= coef_thr, .(from, to, weight, p_value, p_adj)]
  cat("   Edges:", nrow(edges_sub), "\n")
  
  g <- graph_from_data_frame(edges_sub[, .(from, to, weight)], directed = FALSE, vertices = vertex_df_full)
  
  # Leiden
  comm_dt <- run_leiden(g, resolution = LEIDEN_RESOLUTION, n_iter = LEIDEN_N_ITER, seed = LEIDEN_SEED)
  n_comm <- length(unique(comm_dt[community > 0, community]))
  n_isolates <- sum(comm_dt$community == 0)
  cat("   Communities:", n_comm, "  Isolates:", n_isolates, "\n")
  
  m_vec <- setNames(comm_dt$community, comm_dt$node)
  mod_val <- tryCatch(modularity(g, membership = m_vec[V(g)$name], weights = E(g)$weight),
                      error = function(e) NA_real_)
  cat("   Modularity:", round(mod_val, 4), "\n")
  
  # ---- Bridge metrics + connection profiles ----
  cat("   Computing bridge metrics & connection profiles...\n")
  bridge_res <- calc_bridge_metrics(g, comm_dt)
  bridge_dt <- bridge_res$node_dt
  bridge_dt <- merge(bridge_dt, disease_map_dt, by.x = "node", by.y = "code", all.x = TRUE)
  setorder(bridge_dt, -participation_coef)
  bridge_dt[, threshold := tag]
  all_bridge[[tag]] <- copy(bridge_dt)
  fwrite(bridge_dt, file.path(OUT_DIR, sprintf("bridge_metrics_%s.csv", tag)))
  
  # ---- Community summary (enhanced with Top10) ----
  comm_summary <- summarize_communities(comm_dt, g, disease_map_dt, bridge_dt)
  
  # Network metrics
  deg <- degree(g); comps <- components(g)
  giant_size <- comps$csize[which.max(comps$csize)]
  
  metrics_row <- data.table(
    threshold = tag, coef_threshold = coef_thr, alpha_threshold = alpha_thr,
    total_diseases = k, nodes_with_edges = sum(deg > 0), isolates = sum(deg == 0),
    edges = ecount(g), density = round(edge_density(g, loops = FALSE), 6),
    n_components = comps$no, giant_component_nodes = giant_size,
    giant_component_pct = round(100 * giant_size / vcount(g), 2),
    n_communities_leiden = n_comm, modularity_leiden = round(mod_val, 4),
    avg_degree = round(mean(deg), 2), median_degree = median(deg), max_degree = max(deg),
    avg_community_size = round(mean(comm_summary$n_diseases), 2),
    median_community_size = median(comm_summary$n_diseases),
    min_community_size = min(comm_summary$n_diseases),
    max_community_size = max(comm_summary$n_diseases),
    resolution_param = LEIDEN_RESOLUTION
  )
  
  all_membership[[tag]] <- copy(comm_dt)[, threshold := tag]
  all_comm_summary[[tag]] <- copy(comm_summary)[, threshold := tag]
  all_metrics[[idx]] <- metrics_row
  
  # ---- Community-community connectivity matrix ----
  cat("   Computing community-community connectivity matrix...\n")
  cm_res <- calc_community_matrix(g, comm_dt)
  all_comm_matrix[[tag]] <- cm_res$matrix
  if (!is.null(cm_res$bridge_contrib)) {
    all_bridge_contrib[[tag]] <- copy(cm_res$bridge_contrib)[, threshold := tag]
  }
  
  # Heatmap
  p_heat <- tryCatch(
    plot_community_heatmap(cm_res$matrix,
                           sprintf("Inter-community connectivity | coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr)),
    error = function(e) { cat("   Heatmap failed:", e$message, "\n"); NULL }
  )
  if (!is.null(p_heat)) {
    n_c_plot <- nrow(cm_res$matrix)
    heat_size <- max(6, n_c_plot * 0.6 + 2)
    ggsave(file.path(OUT_DIR, sprintf("community_connectivity_heatmap_%s.pdf", tag)),
           p_heat, width = heat_size, height = heat_size, units = "in", device = PDF_DEVICE)
    cat("   Heatmap saved.\n")
  }
  
  # ---- Top 30 bridge bar chart ----
  bridge_active <- bridge_dt[community > 0 & !is.na(participation_coef)]
  if (nrow(bridge_active) >= 1) {
    top_n <- min(30, nrow(bridge_active))
    top_bridge <- head(bridge_active[order(-participation_coef)], top_n)
    top_bridge[, node := factor(node, levels = rev(node))]
    
    p_bridge <- ggplot(top_bridge, aes(x = node, y = participation_coef, fill = ICD10_Chapter)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = round(participation_coef, 3)), hjust = -0.1, size = 3, family = FONT_FAMILY) +
      coord_flip() +
      scale_fill_manual(name = "ICD-10 Chapter", values = chapter_colors_provided, breaks = as.character(1:14)) +
      labs(title = sprintf("Top %d bridge diseases (participation coefficient)", top_n),
           subtitle = sprintf("coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr),
           x = "Disease (ICD-3 code)", y = "Participation coefficient") +
      theme_minimal(base_size = 11) +
      theme(text = element_text(family = FONT_FAMILY),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            axis.text.y = element_text(size = 8)) +
      ylim(0, min(1, max(top_bridge$participation_coef, na.rm = TRUE) * 1.15))
    
    ggsave(file.path(OUT_DIR, sprintf("top_bridge_diseases_%s.pdf", tag)),
           p_bridge, width = 10, height = max(6, top_n * 0.28),
           units = "in", device = PDF_DEVICE)
  }
  
  # ---- z vs P scatter ----
  if (nrow(bridge_active) >= 2) {
    z_max <- max(bridge_active$within_module_z_score, na.rm = TRUE)
    p_role <- ggplot(bridge_active, aes(x = participation_coef, y = within_module_z_score)) +
      geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = 0.30, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = 0.62, linetype = "dotted", color = "grey70", linewidth = 0.4) +
      annotate("text", x = 0.15, y = z_max*0.95, label = "Provincial\nhub", color = "grey45", size = 3, family = FONT_FAMILY) +
      annotate("text", x = 0.50, y = z_max*0.95, label = "Connector\nhub", color = "grey45", size = 3, family = FONT_FAMILY) +
      annotate("text", x = 0.15, y = -0.5, label = "Peripheral", color = "grey45", size = 3, family = FONT_FAMILY) +
      annotate("text", x = 0.50, y = -0.5, label = "Connector", color = "grey45", size = 3, family = FONT_FAMILY) +
      geom_point(aes(color = ICD10_Chapter, size = strength), alpha = 0.75) +
      scale_color_manual(name = "ICD-10 Chapter", values = chapter_colors_provided, breaks = as.character(1:14)) +
      scale_size_continuous(name = "Strength", range = c(2, 8)) +
      labs(title = "Node role classification (Guimer\u00e0 & Amaral)",
           subtitle = sprintf("coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr),
           x = "Participation coefficient (P)", y = "Within-module z-score (z)") +
      theme_minimal(base_size = 11) +
      theme(text = element_text(family = FONT_FAMILY),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
    
    ggsave(file.path(OUT_DIR, sprintf("node_role_scatter_%s.pdf", tag)),
           p_role, width = 10, height = 8, units = "in", device = PDF_DEVICE)
  }
  
  # ---- Global network PDF ----
  p_net <- plot_community_hull(g, comm_dt,
                               sprintf("Multimorbidity network (Leiden) | coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr))
  ggsave(file.path(OUT_DIR, sprintf("community_network_%s.pdf", tag)),
         p_net, width = PDF_W, height = PDF_H, units = "in", device = PDF_DEVICE)
  
  # ---- Edge CSV with community ----
  edges_comm <- copy(edges_sub)
  m_map <- setNames(comm_dt$community, comm_dt$node)
  edges_comm[, from_community := m_map[from]]
  edges_comm[, to_community   := m_map[to]]
  edges_comm[, edge_type := fifelse(from_community == to_community, "intra-community", "inter-community")]
  fwrite(edges_comm, file.path(OUT_DIR, sprintf("edges_with_community_%s.csv", tag)))
  
  # ---- Community subnetwork plots ----
  comm_ids <- sort(unique(comm_dt[community > 0, community]))
  if (length(comm_ids) > 0) {
    sub_dir <- file.path(OUT_DIR, tag)
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    cat("   Plotting", length(comm_ids), "community subnetworks...\n")
    for (cid in comm_ids) {
      n_nodes_c <- sum(comm_dt$community == cid)
      p_sub <- tryCatch(plot_single_community(g, comm_dt, cid, coef_thr, alpha_thr),
                        error = function(e) NULL)
      if (!is.null(p_sub)) {
        ggsave(file.path(sub_dir, sprintf("Class_%02d_%d_diseases.pdf", cid, n_nodes_c)),
               p_sub, width = 10, height = 11, units = "in", device = PDF_DEVICE)
      }
    }
  }
  cat("   Threshold", tag, "complete.\n")
}

# =========================
# 6) Cross-threshold comparison + Multi-seed stability
# =========================
cat("\n==> Step 6) Cross-threshold & multi-seed stability...\n")
metrics_all      <- rbindlist(all_metrics, fill = TRUE)
membership_all   <- rbindlist(all_membership, fill = TRUE)
comm_summary_all <- rbindlist(all_comm_summary, fill = TRUE)

# 6a) Pairwise ARI / NMI across thresholds
tags <- unique(membership_all$threshold)
ari_results <- list()
if (length(tags) >= 2) {
  for (i in 1:(length(tags) - 1)) {
    for (j in (i + 1):length(tags)) {
      dt_i <- membership_all[threshold == tags[i], .(node, ci = community)]
      dt_j <- membership_all[threshold == tags[j], .(node, cj = community)]
      merged <- merge(dt_i, dt_j, by = "node")[ci > 0 & cj > 0]
      if (nrow(merged) > 0) {
        ari <- compare(merged$ci, merged$cj, method = "adjusted.rand")
        nmi <- compare(merged$ci, merged$cj, method = "nmi")
      } else { ari <- NA_real_; nmi <- NA_real_ }
      ari_results[[length(ari_results) + 1]] <- data.table(
        threshold_1 = tags[i], threshold_2 = tags[j],
        n_shared_active_nodes = nrow(merged),
        adjusted_rand_index = round(ari, 4), normalized_mutual_info = round(nmi, 4))
    }
  }
}
ari_dt <- if (length(ari_results) > 0) rbindlist(ari_results) else data.table(note = "Only 1 threshold")

# 6b) Multi-seed stability per threshold (seeds 1..N_SEEDS_STABILITY)
#     ALSO stores all seed memberships for consensus partition
cat("   Multi-seed stability (", N_SEEDS_STABILITY, " seeds per threshold)...\n")
seed_stability <- list()
seed_memberships <- list()   # store actual partitions for consensus

for (idx in 1:nrow(THRESH_GRID)) {
  coef_thr  <- THRESH_GRID[idx, coef]
  alpha_thr <- THRESH_GRID[idx, alpha]
  tag <- sprintf("coef%.2f_alpha%.2f", coef_thr, alpha_thr)
  cat("     ", tag, ": ")
  
  edges_sub <- edge_dt[p_adj < alpha_thr & weight >= coef_thr, .(from, to, weight)]
  g <- graph_from_data_frame(edges_sub, directed = FALSE, vertices = vertex_df_full)
  
  # Reference partition (seed = LEIDEN_SEED)
  ref_comm <- run_leiden(g, LEIDEN_RESOLUTION, LEIDEN_N_ITER, seed = LEIDEN_SEED)
  
  ari_vals <- numeric(N_SEEDS_STABILITY)
  nmi_vals <- numeric(N_SEEDS_STABILITY)
  n_comm_vals <- integer(N_SEEDS_STABILITY)
  seed_mem_list <- vector("list", N_SEEDS_STABILITY)
  
  for (s in 1:N_SEEDS_STABILITY) {
    test_comm <- run_leiden(g, LEIDEN_RESOLUTION, LEIDEN_N_ITER, seed = s)
    n_comm_vals[s] <- length(unique(test_comm[community > 0, community]))
    seed_mem_list[[s]] <- test_comm   # store for consensus
    
    merged <- merge(ref_comm[community > 0], test_comm[community > 0], by = "node",
                    suffixes = c("_ref", "_test"))
    if (nrow(merged) > 0) {
      ari_vals[s] <- compare(merged$community_ref, merged$community_test, method = "adjusted.rand")
      nmi_vals[s] <- compare(merged$community_ref, merged$community_test, method = "nmi")
    } else {
      ari_vals[s] <- NA_real_; nmi_vals[s] <- NA_real_
    }
  }
  
  seed_stability[[tag]] <- data.table(
    threshold = tag, seed = 1:N_SEEDS_STABILITY,
    n_communities = n_comm_vals,
    ARI_vs_reference = round(ari_vals, 4),
    NMI_vs_reference = round(nmi_vals, 4)
  )
  seed_memberships[[tag]] <- list(
    g = g, ref_comm = ref_comm, seed_comms = seed_mem_list,
    ari_vals = ari_vals, nmi_vals = nmi_vals,
    coef_thr = coef_thr, alpha_thr = alpha_thr
  )
  cat("ARI mean=", round(mean(ari_vals, na.rm=TRUE), 3),
      " NMI mean=", round(mean(nmi_vals, na.rm=TRUE), 3), "\n")
}

seed_stab_all <- rbindlist(seed_stability, fill = TRUE)

seed_stab_summary <- seed_stab_all[, .(
  ARI_mean   = round(mean(ARI_vs_reference, na.rm=TRUE), 4),
  ARI_sd     = round(sd(ARI_vs_reference, na.rm=TRUE), 4),
  ARI_min    = round(min(ARI_vs_reference, na.rm=TRUE), 4),
  NMI_mean   = round(mean(NMI_vs_reference, na.rm=TRUE), 4),
  NMI_sd     = round(sd(NMI_vs_reference, na.rm=TRUE), 4),
  NMI_min    = round(min(NMI_vs_reference, na.rm=TRUE), 4),
  n_communities_median = median(n_communities),
  n_communities_range  = paste0(min(n_communities), "-", max(n_communities))
), by = threshold]

# Seed stability box plot
p_stab <- ggplot(seed_stab_all, aes(x = threshold, y = ARI_vs_reference)) +
  geom_boxplot(fill = "#4472C4", alpha = 0.5, width = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  labs(title = sprintf("Leiden partition stability (%d random seeds)", N_SEEDS_STABILITY),
       subtitle = "Adjusted Rand Index vs reference partition (seed=42)",
       x = "Threshold", y = "ARI") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = FONT_FAMILY),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40")) +
  ylim(0, 1.05)
ggsave(file.path(OUT_DIR, "seed_stability_ARI_boxplot.pdf"), p_stab,
       width = 8, height = 6, units = "in", device = PDF_DEVICE)
cat("   Seed stability boxplot saved.\n")

# =========================
# 6c) Consensus Partition + Medoid Seed (per threshold)
# =========================
cat("\n==> Step 6c) Consensus partition & medoid seed...\n")

CONSENSUS_TAU <- 0.5   # standard threshold for co-assignment matrix

consensus_results <- list()

for (tag in names(seed_memberships)) {
  sm <- seed_memberships[[tag]]
  g <- sm$g
  coef_thr <- sm$coef_thr
  alpha_thr <- sm$alpha_thr
  seed_comms <- sm$seed_comms
  n_seeds <- length(seed_comms)
  
  cat(sprintf("   %s: ", tag))
  
  # --- Get active nodes: active (community>0) in ≥50% of seeds ---
  all_nodes <- seed_comms[[1]]$node
  min_active_seeds <- ceiling(0.5 * n_seeds)
  
  active_count <- integer(length(all_nodes))
  names(active_count) <- all_nodes
  for (s in 1:n_seeds) {
    sc <- seed_comms[[s]]
    act <- sc$node[sc$community > 0]
    active_count[act] <- active_count[act] + 1L
  }
  active_nodes <- names(active_count)[active_count >= min_active_seeds]
  n_active <- length(active_nodes)
  
  if (n_active < 2) {
    cat("too few active nodes, skipping.\n")
    next
  }
  cat(n_active, "active nodes (>=50% seeds). ")
  
  # --- Build co-assignment matrix (N_active x N_active) ---
  # C[i,j] = proportion of seeds where node i and j are in the same community
  node_idx <- setNames(seq_along(active_nodes), active_nodes)
  coassign <- matrix(0, n_active, n_active)
  
  for (s in 1:n_seeds) {
    sc <- seed_comms[[s]]
    comm_full <- rep(0L, n_active)
    names(comm_full) <- active_nodes
    sc_active <- sc[node %in% active_nodes]
    matched <- intersect(sc_active$node, active_nodes)
    comm_full[matched] <- setNames(sc_active$community, sc_active$node)[matched]
    
    # Pairwise co-assignment: same community AND both > 0
    for (c_id in unique(comm_full[comm_full > 0])) {
      members <- which(comm_full == c_id)
      if (length(members) >= 2) {
        pairs <- combn(members, 2)
        coassign[pairs[1, ], pairs[2, ]] <- coassign[pairs[1, ], pairs[2, ]] + 1
        coassign[pairs[2, ], pairs[1, ]] <- coassign[pairs[2, ], pairs[1, ]] + 1
      }
    }
  }
  coassign <- coassign / n_seeds
  diag(coassign) <- 1
  
  # --- Threshold at tau -> WEIGHTED consensus graph (not binary) ---
  co_w <- coassign
  co_w[co_w < CONSENSUS_TAU] <- 0
  diag(co_w) <- 0
  
  g_consensus <- graph_from_adjacency_matrix(co_w, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
  V(g_consensus)$name <- active_nodes
  
  # Leiden on weighted consensus graph
  if (ecount(g_consensus) > 0) {
    set.seed(LEIDEN_SEED)
    clu_cons <- cluster_leiden(g_consensus, objective_function = "modularity",
                               resolution = LEIDEN_RESOLUTION, n_iterations = LEIDEN_N_ITER,
                               weights = E(g_consensus)$weight)
    cons_membership <- as.integer(membership(clu_cons))
  } else {
    cons_membership <- rep(1L, n_active)
  }
  
  cons_dt <- data.table(node = active_nodes, community = cons_membership)
  # Add isolates back as community 0
  iso_nodes <- setdiff(all_nodes, active_nodes)
  if (length(iso_nodes) > 0) {
    cons_dt <- rbind(cons_dt, data.table(node = iso_nodes, community = 0L))
  }
  
  n_cons_comm <- length(unique(cons_dt[community > 0, community]))
  cat(n_cons_comm, "consensus communities. ")
  
  # --- Medoid seed: highest mean ARI to all other seeds (on consistent node set) ---
  # Helper: build full membership vector on active_nodes (missing -> 0)
  make_full_mem <- function(dt, nodes) {
    v <- integer(length(nodes)); names(v) <- nodes
    m <- setNames(dt$community, dt$node)
    hit <- intersect(names(m), nodes)
    v[hit] <- m[hit]
    v
  }
  
  ari_matrix <- matrix(NA_real_, n_seeds, n_seeds)
  for (i in 1:n_seeds) {
    vi <- make_full_mem(seed_comms[[i]], active_nodes)
    for (j in i:n_seeds) {
      if (i == j) { ari_matrix[i, j] <- 1; next }
      vj <- make_full_mem(seed_comms[[j]], active_nodes)
      a <- compare(vi, vj, method = "adjusted.rand")
      ari_matrix[i, j] <- a; ari_matrix[j, i] <- a
    }
  }
  mean_ari_per_seed <- rowMeans(ari_matrix, na.rm = TRUE)
  medoid_seed_idx <- which.max(mean_ari_per_seed)
  medoid_dt <- copy(seed_comms[[medoid_seed_idx]])
  
  cat("Medoid seed =", medoid_seed_idx, "\n")
  
  # --- Bridge metrics for consensus partition ---
  cons_bridge_res <- calc_bridge_metrics(g, cons_dt)
  cons_bridge_dt <- cons_bridge_res$node_dt
  cons_bridge_dt <- merge(cons_bridge_dt, disease_map_dt, by.x = "node", by.y = "code", all.x = TRUE)
  
  # --- Bridge metrics for medoid partition ---
  med_bridge_res <- calc_bridge_metrics(g, medoid_dt)
  med_bridge_dt <- med_bridge_res$node_dt
  med_bridge_dt <- merge(med_bridge_dt, disease_map_dt, by.x = "node", by.y = "code", all.x = TRUE)
  
  # --- Consistency check 1: Top 30 bridge overlap (consensus vs medoid) ---
  cons_top30 <- head(cons_bridge_dt[community > 0][order(-participation_coef)], 30)$node
  med_top30  <- head(med_bridge_dt[community > 0][order(-participation_coef)], 30)$node
  overlap_30 <- length(intersect(cons_top30, med_top30))
  jaccard_30 <- length(intersect(cons_top30, med_top30)) / length(union(cons_top30, med_top30))
  
  # --- Consistency check 2: ARI between consensus and medoid ---
  merged_cm <- merge(cons_dt[community > 0], medoid_dt[community > 0], by = "node",
                     suffixes = c("_cons", "_med"))
  ari_cons_med <- if (nrow(merged_cm) > 0) {
    compare(merged_cm$community_cons, merged_cm$community_med, method = "adjusted.rand")
  } else NA_real_
  nmi_cons_med <- if (nrow(merged_cm) > 0) {
    compare(merged_cm$community_cons, merged_cm$community_med, method = "nmi")
  } else NA_real_
  
  # --- Consistency check 3: mean ARI of consensus vs all seeds ---
  ari_cons_seeds <- numeric(n_seeds)
  for (s in 1:n_seeds) {
    ms <- merge(cons_dt[community > 0], seed_comms[[s]][community > 0], by = "node",
                suffixes = c("_cons", "_s"))
    ari_cons_seeds[s] <- if (nrow(ms) > 0) {
      compare(ms$community_cons, ms$community_s, method = "adjusted.rand")
    } else NA_real_
  }
  
  # --- Community summary for consensus ---
  cons_comm_summary <- summarize_communities(cons_dt, g, disease_map_dt, cons_bridge_dt)
  
  # --- Community-community matrix for consensus ---
  cons_cm <- calc_community_matrix(g, cons_dt)
  
  # --- Store results ---
  consensus_results[[tag]] <- list(
    consensus_dt        = cons_dt,
    consensus_bridge    = cons_bridge_dt,
    consensus_summary   = cons_comm_summary,
    consensus_comm_matrix = cons_cm$matrix,
    medoid_seed_idx     = medoid_seed_idx,
    medoid_dt           = medoid_dt,
    medoid_bridge       = med_bridge_dt,
    n_consensus_communities = n_cons_comm,
    n_medoid_communities = length(unique(medoid_dt[community > 0, community])),
    consistency = data.table(
      threshold           = tag,
      tau                 = CONSENSUS_TAU,
      n_consensus_communities = n_cons_comm,
      n_medoid_communities = length(unique(medoid_dt[community > 0, community])),
      medoid_seed         = medoid_seed_idx,
      top30_bridge_overlap = overlap_30,
      top30_bridge_jaccard = round(jaccard_30, 4),
      ARI_consensus_vs_medoid = round(ari_cons_med, 4),
      NMI_consensus_vs_medoid = round(nmi_cons_med, 4),
      ARI_consensus_vs_seeds_mean = round(mean(ari_cons_seeds, na.rm=TRUE), 4),
      ARI_consensus_vs_seeds_min  = round(min(ari_cons_seeds, na.rm=TRUE), 4)
    )
  )
  
  # --- Consensus network plot ---
  p_cons <- plot_community_hull(g, cons_dt,
                                sprintf("Consensus partition (tau=%.2f) | coef >= %.2f, FDR < %.2f",
                                        CONSENSUS_TAU, coef_thr, alpha_thr))
  ggsave(file.path(OUT_DIR, sprintf("consensus_network_%s.pdf", tag)),
         p_cons, width = PDF_W, height = PDF_H, units = "in", device = PDF_DEVICE)
  
  # --- Consensus community subnetwork plots ---
  cons_comm_ids <- sort(unique(cons_dt[community > 0, community]))
  if (length(cons_comm_ids) > 0) {
    cons_sub_dir <- file.path(OUT_DIR, paste0("consensus_", tag))
    dir.create(cons_sub_dir, recursive = TRUE, showWarnings = FALSE)
    for (cid in cons_comm_ids) {
      n_c <- sum(cons_dt$community == cid)
      p_sub <- tryCatch(plot_single_community(g, cons_dt, cid, coef_thr, alpha_thr),
                        error = function(e) NULL)
      if (!is.null(p_sub)) {
        ggsave(file.path(cons_sub_dir, sprintf("Class_%02d_%d_diseases.pdf", cid, n_c)),
               p_sub, width = 10, height = 11, units = "in", device = PDF_DEVICE)
      }
    }
  }
  
  # --- Consensus bridge bar chart ---
  cons_active <- cons_bridge_dt[community > 0 & !is.na(participation_coef)]
  if (nrow(cons_active) >= 1) {
    top_n <- min(30, nrow(cons_active))
    top_cb <- head(cons_active[order(-participation_coef)], top_n)
    top_cb[, node := factor(node, levels = rev(node))]
    p_cb <- ggplot(top_cb, aes(x = node, y = participation_coef, fill = ICD10_Chapter)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = round(participation_coef, 3)), hjust = -0.1, size = 3, family = FONT_FAMILY) +
      coord_flip() +
      scale_fill_manual(name = "ICD-10 Chapter", values = chapter_colors_provided, breaks = as.character(1:14)) +
      labs(title = sprintf("Top %d bridge diseases (consensus, tau=%.2f)", top_n, CONSENSUS_TAU),
           subtitle = sprintf("coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr),
           x = "Disease", y = "Participation coefficient") +
      theme_minimal(base_size = 11) +
      theme(text = element_text(family = FONT_FAMILY),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            axis.text.y = element_text(size = 8)) +
      ylim(0, min(1, max(top_cb$participation_coef, na.rm = TRUE) * 1.15))
    ggsave(file.path(OUT_DIR, sprintf("consensus_top_bridge_%s.pdf", tag)),
           p_cb, width = 10, height = max(6, top_n * 0.28),
           units = "in", device = PDF_DEVICE)
  }
  
  # --- Consensus heatmap ---
  p_heat_cons <- tryCatch(
    plot_community_heatmap(cons_cm$matrix,
                           sprintf("Consensus inter-community connectivity | %s", tag)),
    error = function(e) NULL)
  if (!is.null(p_heat_cons)) {
    n_c_plot <- nrow(cons_cm$matrix)
    ggsave(file.path(OUT_DIR, sprintf("consensus_heatmap_%s.pdf", tag)),
           p_heat_cons, width = max(6, n_c_plot*0.6+2), height = max(6, n_c_plot*0.6+2),
           units = "in", device = PDF_DEVICE)
  }
  
  cat("   ", tag, "consensus outputs saved.\n")
  
  # =================================================================
  # --- MEDOID SEED: full output suite (PRIMARY partition) ---
  # =================================================================
  cat("   ", tag, "generating medoid seed outputs (primary)...\n")
  
  # Community summary for medoid
  med_comm_summary <- summarize_communities(medoid_dt, g, disease_map_dt, med_bridge_dt)
  
  # Community-community matrix for medoid
  med_cm <- calc_community_matrix(g, medoid_dt)
  
  n_med_comm <- length(unique(medoid_dt[community > 0, community]))
  
  # Store medoid summary + matrix in consensus_results for Excel
  consensus_results[[tag]]$medoid_summary    <- med_comm_summary
  consensus_results[[tag]]$medoid_comm_matrix <- med_cm$matrix
  if (!is.null(med_cm$bridge_contrib)) {
    consensus_results[[tag]]$medoid_bridge_contrib <- med_cm$bridge_contrib
  }
  
  # --- Medoid network plot ---
  p_med <- plot_community_hull(g, medoid_dt,
                               sprintf("Medoid partition (seed=%d) | coef >= %.2f, FDR < %.2f",
                                       medoid_seed_idx, coef_thr, alpha_thr))
  ggsave(file.path(OUT_DIR, sprintf("medoid_network_%s.pdf", tag)),
         p_med, width = PDF_W, height = PDF_H, units = "in", device = PDF_DEVICE)
  
  # --- Medoid community subnetwork plots ---
  med_comm_ids <- sort(unique(medoid_dt[community > 0, community]))
  if (length(med_comm_ids) > 0) {
    med_sub_dir <- file.path(OUT_DIR, paste0("medoid_", tag))
    dir.create(med_sub_dir, recursive = TRUE, showWarnings = FALSE)
    for (cid in med_comm_ids) {
      n_c <- sum(medoid_dt$community == cid)
      p_sub <- tryCatch(plot_single_community(g, medoid_dt, cid, coef_thr, alpha_thr),
                        error = function(e) NULL)
      if (!is.null(p_sub)) {
        ggsave(file.path(med_sub_dir, sprintf("Class_%02d_%d_diseases.pdf", cid, n_c)),
               p_sub, width = 10, height = 11, units = "in", device = PDF_DEVICE)
      }
    }
  }
  
  # --- Medoid bridge bar chart ---
  med_active <- med_bridge_dt[community > 0 & !is.na(participation_coef)]
  if (nrow(med_active) >= 1) {
    top_n <- min(30, nrow(med_active))
    top_mb <- head(med_active[order(-participation_coef)], top_n)
    top_mb[, node := factor(node, levels = rev(node))]
    p_mb <- ggplot(top_mb, aes(x = node, y = participation_coef, fill = ICD10_Chapter)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = round(participation_coef, 3)), hjust = -0.1, size = 3, family = FONT_FAMILY) +
      coord_flip() +
      scale_fill_manual(name = "ICD-10 Chapter", values = chapter_colors_provided, breaks = as.character(1:14)) +
      labs(title = sprintf("Top %d bridge diseases (medoid seed=%d)", top_n, medoid_seed_idx),
           subtitle = sprintf("coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr),
           x = "Disease", y = "Participation coefficient") +
      theme_minimal(base_size = 11) +
      theme(text = element_text(family = FONT_FAMILY),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            axis.text.y = element_text(size = 8)) +
      ylim(0, min(1, max(top_mb$participation_coef, na.rm = TRUE) * 1.15))
    ggsave(file.path(OUT_DIR, sprintf("medoid_top_bridge_%s.pdf", tag)),
           p_mb, width = 10, height = max(6, top_n * 0.28),
           units = "in", device = PDF_DEVICE)
  }
  
  # --- Medoid z vs P scatter (node role) ---
  if (nrow(med_active) >= 2) {
    z_max_m <- max(med_active$within_module_z_score, na.rm = TRUE)
    p_role_m <- ggplot(med_active, aes(x = participation_coef, y = within_module_z_score)) +
      geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = 0.30, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = 0.62, linetype = "dotted", color = "grey70", linewidth = 0.4) +
      annotate("text", x = 0.15, y = z_max_m*0.95, label = "Provincial\nhub", color = "grey45", size = 3, family = FONT_FAMILY) +
      annotate("text", x = 0.50, y = z_max_m*0.95, label = "Connector\nhub", color = "grey45", size = 3, family = FONT_FAMILY) +
      annotate("text", x = 0.15, y = -0.5, label = "Peripheral", color = "grey45", size = 3, family = FONT_FAMILY) +
      annotate("text", x = 0.50, y = -0.5, label = "Connector", color = "grey45", size = 3, family = FONT_FAMILY) +
      geom_point(aes(color = ICD10_Chapter, size = strength), alpha = 0.75) +
      scale_color_manual(name = "ICD-10 Chapter", values = chapter_colors_provided, breaks = as.character(1:14)) +
      scale_size_continuous(name = "Strength", range = c(2, 8)) +
      labs(title = sprintf("Node role classification (medoid seed=%d)", medoid_seed_idx),
           subtitle = sprintf("coef >= %.2f, FDR < %.2f", coef_thr, alpha_thr),
           x = "Participation coefficient (P)", y = "Within-module z-score (z)") +
      theme_minimal(base_size = 11) +
      theme(text = element_text(family = FONT_FAMILY),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
    ggsave(file.path(OUT_DIR, sprintf("medoid_node_role_scatter_%s.pdf", tag)),
           p_role_m, width = 10, height = 8, units = "in", device = PDF_DEVICE)
  }
  
  # --- Medoid heatmap ---
  p_heat_med <- tryCatch(
    plot_community_heatmap(med_cm$matrix,
                           sprintf("Medoid inter-community connectivity (seed=%d) | %s", medoid_seed_idx, tag)),
    error = function(e) NULL)
  if (!is.null(p_heat_med)) {
    n_c_plot <- nrow(med_cm$matrix)
    ggsave(file.path(OUT_DIR, sprintf("medoid_heatmap_%s.pdf", tag)),
           p_heat_med, width = max(6, n_c_plot*0.6+2), height = max(6, n_c_plot*0.6+2),
           units = "in", device = PDF_DEVICE)
  }
  
  # --- Medoid edge CSV ---
  edges_sub_tag <- edge_dt[p_adj < alpha_thr & weight >= coef_thr, .(from, to, weight, p_value, p_adj)]
  med_map <- setNames(medoid_dt$community, medoid_dt$node)
  edges_sub_tag[, from_community := med_map[from]]
  edges_sub_tag[, to_community   := med_map[to]]
  edges_sub_tag[, edge_type := fifelse(from_community == to_community, "intra-community", "inter-community")]
  fwrite(edges_sub_tag, file.path(OUT_DIR, sprintf("medoid_edges_with_community_%s.csv", tag)))
  
  # --- Medoid bridge CSV ---
  fwrite(med_bridge_dt, file.path(OUT_DIR, sprintf("medoid_bridge_metrics_%s.csv", tag)))
  
  cat("   ", tag, "medoid outputs saved.\n")
}

# Collect consistency summary
consistency_all <- rbindlist(lapply(consensus_results, function(x) x$consistency), fill = TRUE)

# =========================
# 7) Styled Excel
# =========================
cat("==> Step 7) Writing styled Excel...\n")
xlsx_file <- file.path(OUT_DIR, "Leiden_community_detection_results.xlsx")
wb <- createWorkbook()

headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                           valign = "center", fgFill = "#4472C4", textDecoration = "bold")
add_sheet <- function(sheet_name, dt_in) {
  sn <- substr(sheet_name, 1, 31)
  addWorksheet(wb, sn)
  writeDataTable(wb, sn, dt_in, withFilter = TRUE)
  setColWidths(wb, sn, cols = 1:ncol(dt_in), widths = "auto")
  addStyle(wb, sn, headerStyle, rows = 1, cols = 1:ncol(dt_in), gridExpand = TRUE)
}

add_sheet("metrics_comparison", metrics_all)
add_sheet("cross_threshold_ARI", ari_dt)
add_sheet("seed_stability_summary", seed_stab_summary)
add_sheet("seed_stability_detail", seed_stab_all)

for (tag in names(all_comm_summary)) add_sheet(paste0("summary_", tag), all_comm_summary[[tag]])
for (tag in names(all_membership)) {
  mem <- merge(copy(all_membership[[tag]]), disease_map_dt, by.x = "node", by.y = "code", all.x = TRUE)
  setorder(mem, community, node)
  add_sheet(paste0("members_", tag), mem)
}

bridge_all <- rbindlist(all_bridge, fill = TRUE)
for (tag in names(all_bridge)) {
  add_sheet(paste0("bridge_", tag), all_bridge[[tag]][order(-participation_coef)])
}

# Community-community connectivity matrices as sheets
for (tag in names(all_comm_matrix)) {
  mat <- all_comm_matrix[[tag]]
  mat_dt <- as.data.table(mat, keep.rownames = "community")
  add_sheet(paste0("comm_matrix_", tag), mat_dt)
}

# Inter-community bridge contributors
if (length(all_bridge_contrib) > 0) {
  bc_all <- rbindlist(all_bridge_contrib, fill = TRUE)
  add_sheet("inter_comm_bridges", bc_all)
}

# --- Consensus partition sheets ---
if (exists("consistency_all") && nrow(consistency_all) > 0) {
  add_sheet("consensus_consistency", consistency_all)
}

for (tag in names(consensus_results)) {
  cr <- consensus_results[[tag]]
  
  # Consensus membership
  cons_mem <- merge(copy(cr$consensus_dt), disease_map_dt,
                    by.x = "node", by.y = "code", all.x = TRUE)
  setorder(cons_mem, community, node)
  add_sheet(paste0("cons_members_", tag), cons_mem)
  
  # Consensus community summary
  add_sheet(paste0("cons_summary_", tag), cr$consensus_summary)
  
  # Consensus bridge metrics
  cons_br <- copy(cr$consensus_bridge)
  setorder(cons_br, -participation_coef)
  add_sheet(paste0("cons_bridge_", tag), cons_br)
  
  # Consensus community-community matrix
  cmat <- cr$consensus_comm_matrix
  cmat_dt <- as.data.table(cmat, keep.rownames = "community")
  add_sheet(paste0("cons_cmatrix_", tag), cmat_dt)
  
  # Medoid membership
  med_mem <- merge(copy(cr$medoid_dt), disease_map_dt,
                   by.x = "node", by.y = "code", all.x = TRUE)
  setorder(med_mem, community, node)
  add_sheet(paste0("medoid_members_", tag), med_mem)
  
  # Medoid community summary
  if (!is.null(cr$medoid_summary)) {
    add_sheet(paste0("medoid_summary_", tag), cr$medoid_summary)
  }
  
  # Medoid bridge metrics
  med_br <- copy(cr$medoid_bridge)
  setorder(med_br, -participation_coef)
  add_sheet(paste0("medoid_bridge_", tag), med_br)
  
  # Medoid community-community matrix
  if (!is.null(cr$medoid_comm_matrix)) {
    mmat <- cr$medoid_comm_matrix
    mmat_dt <- as.data.table(mmat, keep.rownames = "community")
    add_sheet(paste0("med_cmatrix_", tag), mmat_dt)
  }
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("   Saved:", xlsx_file, "\n")

# =========================
# 8) Summary bar charts
# =========================
cat("==> Step 8) Summary bar charts...\n")

p_bar <- ggplot(metrics_all, aes(x = threshold, y = n_communities_leiden)) +
  geom_col(fill = "#4472C4", width = 0.6) +
  geom_text(aes(label = n_communities_leiden), vjust = -0.5, size = 5, family = FONT_FAMILY) +
  labs(title = "Number of Leiden communities by threshold",
       subtitle = sprintf("Resolution = %.1f", LEIDEN_RESOLUTION), x = "Threshold", y = "Communities") +
  theme_minimal(base_size = 13) +
  theme(text = element_text(family = FONT_FAMILY), plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"), axis.text.x = element_text(angle = 15, hjust = 1)) +
  ylim(0, max(metrics_all$n_communities_leiden, na.rm = TRUE) * 1.25)
ggsave(file.path(OUT_DIR, "communities_count_comparison.pdf"), p_bar,
       width = 8, height = 6, units = "in", device = PDF_DEVICE)

p_mod <- ggplot(metrics_all, aes(x = threshold, y = modularity_leiden)) +
  geom_col(fill = "#ED7D31", width = 0.6) +
  geom_text(aes(label = modularity_leiden), vjust = -0.5, size = 5, family = FONT_FAMILY) +
  labs(title = "Modularity (Leiden) by threshold", x = "Threshold", y = "Modularity") +
  theme_minimal(base_size = 13) +
  theme(text = element_text(family = FONT_FAMILY), plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1)) + ylim(0, 1)
ggsave(file.path(OUT_DIR, "modularity_comparison.pdf"), p_mod,
       width = 8, height = 6, units = "in", device = PDF_DEVICE)

# =========================
# 9) RDS intermediates
# =========================
saveRDS(membership_all,      file.path(OUT_DIR, "leiden_membership_all.rds"))
saveRDS(comm_summary_all,    file.path(OUT_DIR, "leiden_community_summary_all.rds"))
saveRDS(metrics_all,         file.path(OUT_DIR, "leiden_metrics_comparison.rds"))
saveRDS(bridge_all,          file.path(OUT_DIR, "bridge_metrics_all.rds"))
saveRDS(seed_stab_all,       file.path(OUT_DIR, "seed_stability_all.rds"))
saveRDS(consensus_results,   file.path(OUT_DIR, "consensus_results_all.rds"))

cat("\n========================================\n")
cat("DONE. All results saved to:\n", OUT_DIR, "\n")
cat("========================================\n")
cat("\nOutputs:\n")
cat("  === PRIMARY (medoid seed partition) ===\n")
cat("  1)  medoid_network_*.pdf               (global network, hull borders)\n")
cat("  2)  medoid_*/ subfolders               (community subnetwork PDFs)\n")
cat("  3)  medoid_top_bridge_*.pdf            (Top 30 bridge diseases)\n")
cat("  4)  medoid_node_role_scatter_*.pdf     (z vs P, Guimera & Amaral)\n")
cat("  5)  medoid_heatmap_*.pdf               (community-community connectivity)\n")
cat("  6)  medoid_edges_with_community_*.csv\n")
cat("  7)  medoid_bridge_metrics_*.csv\n")
cat("  === SUPPLEMENTARY (single-seed=42 + consensus) ===\n")
cat("  8)  community_network_*.pdf            (seed=42 Leiden)\n")
cat("  9)  consensus_network_*.pdf            (consensus tau=0.5)\n")
cat("  10) coef*_alpha*/ subfolders           (seed=42 community subnetworks)\n")
cat("  11) consensus_*/ subfolders            (consensus community subnetworks)\n")
cat("  12) top_bridge_diseases_*.pdf / consensus_top_bridge_*.pdf\n")
cat("  13) node_role_scatter_*.pdf            (seed=42)\n")
cat("  14) community_connectivity_heatmap_*.pdf / consensus_heatmap_*.pdf\n")
cat("  === STABILITY & COMPARISON ===\n")
cat("  15) seed_stability_ARI_boxplot.pdf\n")
cat("  16) communities_count_comparison.pdf / modularity_comparison.pdf\n")
cat("  === EXCEL ===\n")
cat("  17) Leiden_community_detection_results.xlsx\n")
cat("      Sheets: metrics_comparison, cross_threshold_ARI\n")
cat("      Sheets: seed_stability_summary / detail, consensus_consistency\n")
cat("      Sheets: medoid_members_*, medoid_summary_*, medoid_bridge_*, med_cmatrix_*\n")
cat("      Sheets: summary_*, members_*, bridge_*, comm_matrix_*  (seed=42)\n")
cat("      Sheets: cons_members_*, cons_summary_*, cons_bridge_*, cons_cmatrix_*\n")
cat("      Sheets: inter_comm_bridges\n")
cat("  === RDS ===\n")
cat("  18) RDS intermediates (incl. consensus_results_all.rds)\n")