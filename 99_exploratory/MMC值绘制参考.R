###############################################################################
# 极坐标 MMC 柱状图函数（最终修正版）
#   • 将图例字体样式移至 guide_legend() 内部，确保其生效
###############################################################################

# ─────────────────────────────────────────────────────────────
# 0. 依赖
# ─────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)

# ─────────────────────────────────────────────────────────────
# 1. 绘图函数 (最终修正版)
# ─────────────────────────────────────────────────────────────
plot_circular <- function(mmc_file,
                          plot_title,
                          show_legend        = TRUE,
                          title_vjust        = -2,
                          title_size         = 20,
                          title_margin_b     = -4,
                          legend_title_size  = 14,
                          legend_text_size   = 12,
                          legend_key_width_cm= 0.6,
                          legend_key_height_cm= 0.6) {
  
  # 1.1 读数据 ─────
  if (file.exists(mmc_file)) {
    mmc_df <- read_csv(mmc_file, show_col_types = FALSE)
  } else {
    set.seed(123)
    mmc_df <- data.frame(
      DiseaseCode = paste0(sample(LETTERS, 150, TRUE), sample(0:99, 150, TRUE)),
      MMC = abs(rnorm(150, 5, 4)) + 1
    )
  }
  
  # 1.2 参数 ─────
  top_n_diseases        <- 100
  inner_radius_proportion <- 1
  label_size            <- 3
  axis_label_size       <- 3.9
  grid_line_color       <- "grey60"
  axis_line_color       <- "grey90"
  axis_text_color       <- "black"
  start_padding         <- 1
  end_padding           <- 5
  
  chapter_cols <- c(
    "1" = "#B70031", "2"="#EC6A3D", "3"="#A1C6AD", "4"="#0098C5",
    "5" = "#9A6CAC", "6"="#004B9B", "7"="#C0C7CB", "8"="#AAD9EA",
    "9" = "#FFE100", "10"="#E71A10", "11"="#6EB327", "12"="#8E006C",
    "13" = "#D2ACC9", "14"="#ADB5DA"
  )
  
  # 1.3 预处理 ─────
  mmc_processed <- mmc_df |>
    mutate(MMC = as.numeric(MMC)) |>
    filter(!is.na(MMC) & MMC >= 0) |>
    arrange(desc(MMC)) |>
    slice_head(n = top_n_diseases) |>
    mutate(
      DiseaseCode = as.character(DiseaseCode),
      chapter = case_when(
        DiseaseCode >= "A00" & DiseaseCode <= "B99" ~ "1",
        DiseaseCode >= "C00" & DiseaseCode <= "D48" ~ "2",
        DiseaseCode >= "D50" & DiseaseCode <= "D89" ~ "3",
        DiseaseCode >= "E00" & DiseaseCode <= "E90" ~ "4",
        DiseaseCode >= "F00" & DiseaseCode <= "F99" ~ "5",
        DiseaseCode >= "G00" & DiseaseCode <= "G99" ~ "6",
        DiseaseCode >= "H00" & DiseaseCode <= "H59" ~ "7",
        DiseaseCode >= "H60" & DiseaseCode <= "H95" ~ "8",
        DiseaseCode >= "I00" & DiseaseCode <= "I99" ~ "9",
        DiseaseCode >= "J00" & DiseaseCode <= "J99" ~ "10",
        DiseaseCode >= "K00" & DiseaseCode <= "K93" ~ "11",
        DiseaseCode >= "L00" & DiseaseCode <= "L99" ~ "12",
        DiseaseCode >= "M00" & DiseaseCode <= "M99" ~ "13",
        DiseaseCode >= "N00" & DiseaseCode <= "N99" ~ "14",
        TRUE ~ "0"
      ),
      chapter_factor   = factor(chapter, levels = names(chapter_cols)),
      chapter_sort_key = if_else(chapter == "0", 99, as.numeric(chapter))
    ) |>
    arrange(chapter_sort_key, desc(MMC)) |>
    mutate(id = row_number() - 1, x_pos = id + start_padding)
  
  # 1.4 计算坐标 ─────
  y_max         <- max(mmc_processed$MMC, na.rm = TRUE)
  inner_radius  <- y_max * inner_radius_proportion
  mmc_breaks    <- c(0.5, 1, 1.5, 2, 2.5, 3)
  radial_breaks <- inner_radius + mmc_breaks
  plot_x_limits <- c(0, (top_n_diseases - 1) + start_padding + end_padding)
  offset        <- y_max * 0.05
  label_data    <- mmc_processed |>
    mutate(
      label_y = inner_radius + MMC + offset,
      angle   = 90 - (id + 0.5) / n() * 360,
      hjust   = ifelse(angle > -90 & angle < 90, 0, 1),
      angle   = ifelse(angle < -90, angle + 180, angle)
    )
  y_limit_max   <- (inner_radius + 3) * 1.05
  
  # 1.5 绘图 ─────
  p <- ggplot(mmc_processed) +
    geom_rect(aes(xmin = x_pos - 0.5, xmax = x_pos + 0.5,
                  ymin = inner_radius, ymax = inner_radius + MMC,
                  fill = chapter_factor),
              colour = "white", linewidth = 0.1) +
    annotate("segment", x = 0, xend = 0, y = inner_radius, yend = inner_radius + 3,
             colour = axis_line_color, linewidth = 0.6) +
    lapply(seq_along(radial_breaks), function(i) {
      annotate("text", x = 0, y = radial_breaks[i], label = mmc_breaks[i],
               colour = axis_text_color, size = axis_label_size, hjust = 1.1, vjust = 0.5)
    }) +
    geom_text(data = label_data, aes(x = x_pos, y = label_y, label = DiseaseCode, angle = angle, hjust = hjust),
              size = label_size) +
    coord_polar(theta = "x", start = 0, direction = 1, clip = "off") +
    scale_x_continuous(limits = plot_x_limits) +
    scale_y_continuous(limits = c(0, y_limit_max), breaks = radial_breaks, labels = NULL) +
    ggtitle(plot_title) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text          = element_blank(),
      axis.title         = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = grid_line_color, linewidth = 0.3),
      plot.title.position = "plot",
      plot.title         = element_text(
        face = "bold", hjust = 0, vjust = title_vjust,
        size = title_size, margin = margin(b = title_margin_b)),
      # 移除这里的图例字体设置，因为它们将被 guide_legend 覆盖
      # legend.title       = element_text(...),
      # legend.text        = element_text(...),
      plot.margin        = margin(0, 0, 0, 0)
    )
  
  # 1.6 图例控制 (核心修正) ─────
  if (show_legend) {
    p <- p + scale_fill_manual(
      values = chapter_cols,
      name   = "ICD-10 Chapter",
      drop   = FALSE,
      guide  = guide_legend(
        # 色块（Key）尺寸控制
        keywidth  = unit(legend_key_width_cm, 'cm'),
        keyheight = unit(legend_key_height_cm, 'cm'),
        
        # 字体样式控制（高优先级） <--- 核心修正
        title.theme = element_text(face = "bold", size = legend_title_size),
        label.theme = element_text(size = legend_text_size)
      )
    )
  } else {
    p <- p + scale_fill_manual(values = chapter_cols, guide = "none")
  }
  
  return(p)
}


# ─────────────────────────────────────────────────────────────
# 2. CSV 路径
# ─────────────────────────────────────────────────────────────
csv_overall <- 
csv_female  <- 
csv_male    <- 

# ─────────────────────────────────────────────────────────────
# 3. 生成三张子图 (现在所有图例参数都将正确生效)
# ─────────────────────────────────────────────────────────────
p_overall <- plot_circular(
  csv_overall, "Overall",
  title_vjust = -6,
  legend_title_size  = 16, # <-- 现在这个参数会生效
  legend_text_size   = 14, # <-- 现在这个参数会生效
  legend_key_width_cm= 0.8,  # <-- 这个参数之前就生效
  legend_key_height_cm= 0.8   # <-- 这个参数之前就生效
)

p_female  <- plot_circular(
  csv_female,  "Female",
  show_legend = FALSE, 
  title_vjust = -4
)

p_male    <- plot_circular(
  csv_male,    "Male",
  show_legend = FALSE,
  title_vjust = -4
)

# ─────────────────────────────────────────────────────────────
# 4. 组合
# ─────────────────────────────────────────────────────────────
legend_obj <- get_legend(p_overall)
p_overall_clean <- p_overall + theme(legend.position = "none")

left_column  <- plot_grid(p_overall_clean, p_male, ncol = 1, rel_heights = c(1, 1))
right_column <- plot_grid(p_female, legend_obj, ncol = 1, rel_heights = c(1, 1))
combined_plot <- plot_grid(left_column, right_column, ncol = 2, rel_widths = c(1, 0.9))

print(combined_plot)