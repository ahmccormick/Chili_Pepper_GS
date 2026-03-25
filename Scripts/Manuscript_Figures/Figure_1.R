###################################################################################################
# FIGURE 1 - GWAS Circos + Manhattan bubble plot
# Uses new BLUE-adjusted GWAS tables:
#   - TableS5: significant SNPs in >= 2 models (per trait)
#   - TableS7: multi-trait SNPs >= 3 traits and >= 3 models
###################################################################################################

rm(list = ls())
library(tidyverse)
library(circlize)
library(ggrepel)

#############################################
# 0) Shared setup
#############################################

chr_lengths <- c(333, 177, 290, 249, 255, 253, 266, 174, 278, 210, 275, 260)
chr_names   <- paste0("chr", 1:12)
chr_offsets <- cumsum(c(0, chr_lengths[-12] + 15))
chr_mids    <- chr_offsets + chr_lengths / 2

cat_colors <- c(
  "Canopy spectral"  = "#4E79A7",
  "Leaf area"        = "#59A14F",
  "Yield components" = "#E15759",
  "Growth"           = "#F28E2B",
  "Physiology"       = "#76B7B2",
  "Phenology"        = "#B07AA1",
  "Pollen"           = "#A0522D",
  "Yield"            = "#E8C245"
)

cat_levels_order <- c(
  "Canopy spectral", "Leaf area", "Yield components",
  "Growth", "Physiology", "Phenology", "Pollen", "Yield"
)

outdir <- "~/R/WorldVeg_Capsicum/STI_GWAS/GWAS_tables/"

#############################################
# 1) Load new tables
#############################################

# S5: per trait per SNP, significant in >= 2 models
s5 <- read.csv(file.path(outdir, "TableS5_significant_SNPs_2plus_models.csv"),
               stringsAsFactors = FALSE)

# S7: multi-trait SNPs >= 3 traits AND >= 3 models
s7 <- read.csv(file.path(outdir, "TableS7_multitrait_SNPs_3plus_models.csv"),
               stringsAsFactors = FALSE) %>%
  select(-any_of(c("X", "X.1")))  # remove blank columns if present

message("S5 rows: ", nrow(s5))
message("S7 rows: ", nrow(s7))

#############################################
# 2) Prepare S5 for circos outer tracks
#    (one row per trait-SNP with n_models)
#############################################

df_a <- s5 %>%
  mutate(Category = factor(Category, levels = cat_levels_order))

#############################################
# 3) Prepare S7 for inner track
#############################################

df_b <- s7 %>%
  mutate(
    dominant_category = factor(dominant_category, levels = cat_levels_order),
    x_pos = chr_offsets[Chr] + Pos / 1e6
  )

message("S7 SNPs for inner track: ", nrow(df_b))
print(df_b %>% count(dominant_category))

#############################################
# 4) prep_track helper
#############################################

prep_track <- function(df, chr_col = "Chr", pos_col = "Pos") {
  df %>%
    mutate(
      chr   = paste0("chr", .data[[chr_col]]),
      start = .data[[pos_col]],
      end   = .data[[pos_col]] + 1
    )
}

df_a_track <- prep_track(df_a)
df_b_track <- prep_track(df_b)

#############################################
# 5) FIGURE 1A — Circos plot
#############################################

pdf(file.path(outdir, "Figure1A_GWAS_circos_newdata.pdf"), width = 10, height = 10)

circos.clear()
circos.par(
  "track.height" = 0.08,
  "cell.padding" = c(0, 0, 0, 0),
  "gap.degree"   = 2,
  "start.degree" = 90,
  "clock.wise"   = TRUE
)

cytoband <- data.frame(
  V1 = chr_names,
  V2 = 0,
  V3 = chr_lengths * 1e6,
  V4 = chr_names,
  V5 = "gpos50"
)

circos.initializeWithIdeogram(
  cytoband,
  chromosome.index = chr_names,
  plotType         = c("axis", "labels"),
  axis.labels.cex  = 0.8,
  labels.cex       = 1.4
)

# outer tracks — one per category from S5
cats <- levels(df_a$Category)

for (cat in cats) {
  df_cat <- df_a_track %>% filter(Category == cat)
  col    <- cat_colors[cat]
  
  circos.track(
    ylim         = c(0, 4),
    track.height = 0.07,
    bg.border    = "grey90",
    bg.col       = "grey99",
    panel.fun = function(region, value, ...) {
      chr_now <- CELL_META$sector.index
      df_chr  <- df_cat %>% filter(chr == chr_now)
      if (nrow(df_chr) == 0) return()
      circos.points(
        x   = df_chr$start,
        y   = df_chr$n_models,
        pch = 16,
        cex = df_chr$n_models * 0.3,
        col = adjustcolor(col, alpha.f = 0.7)
      )
    }
  )
}

# innermost track — S7 multi-trait SNPs
circos.track(
  ylim         = c(0, 12),
  track.height = 0.12,
  bg.border    = "grey85",
  bg.col       = "white",
  panel.fun = function(region, value, ...) {
    chr_now <- CELL_META$sector.index
    df_chr  <- df_b_track %>% filter(chr == chr_now)
    if (nrow(df_chr) == 0) return()
    
    # all SNPs as circles
    circos.points(
      x   = df_chr$start,
      y   = df_chr$n_traits,
      pch = 21,
      cex = 0.3 + (df_chr$n_traits / 12) * 0.9,
      col = "black",
      bg  = cat_colors[as.character(df_chr$dominant_category)]
    )
    
    # additional black ring for 4-model hits
    df_4model <- df_chr %>% filter(n_models_total == 4)
    if (nrow(df_4model) > 0) {
      circos.points(
        x   = df_4model$start,
        y   = df_4model$n_traits,
        pch = 21,
        cex = 0.3 + (df_4model$n_traits / 12) * 0.9,
        col = "black",
        bg  = NA
      )
    }
  }
)

legend(
  "bottomright",
  legend = names(cat_colors),
  fill   = cat_colors,
  border = NA,
  cex    = 0.7,
  title  = "Trait category",
  bty    = "n"
)

circos.clear()
dev.off()
message("Circos saved.")

#############################################
#############################################
# 6) FIGURE 1B — Manhattan bubble plot
#############################################

# install.packages(c("scatterpie", "ggrepel", "ggforce", "tidyverse", "scales"))
library(scatterpie)
library(tidyverse)
library(ggrepel)
library(scales)
library(ggforce)

gwas_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12),
    axis.text.y        = element_text(size = 10, face = "bold"),
    axis.title.y       = element_blank(),
    axis.title.x       = element_text(size = 13),
    legend.position    = "bottom",
    legend.box         = "vertical",
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    panel.background   = element_rect(fill = "white", color = NA)
  )

# y-axis ordering
single_cats <- cat_levels_order[
  cat_levels_order %in% unique(df_b$Categories[df_b$n_categories == 1])
]
multi_cats <- sort(unique(df_b$Categories[df_b$n_categories > 1]))
cats_present <- c(single_cats, multi_cats)

y_levels <- rev(cats_present)
y_lookup <- setNames(seq_along(y_levels), y_levels)

df_b_plot <- df_b %>%
  mutate(
    Categories = factor(Categories, levels = cats_present),
    y_num      = y_lookup[as.character(Categories)]
  )

# ------------------------------------------
# normalize x within each chromosome
# ------------------------------------------
df_b_plot <- df_b_plot %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = FALSE, convert = TRUE) %>%
  mutate(
    CHR = as.integer(CHR),
    BP  = as.numeric(BP)
  )

chr_lengths <- df_b_plot %>%
  group_by(CHR) %>%
  summarise(chr_max = max(BP, na.rm = TRUE), .groups = "drop")

df_b_plot <- df_b_plot %>%
  left_join(chr_lengths, by = "CHR") %>%
  mutate(
    rel_pos = if_else(chr_max > 0, BP / chr_max, 0.5),
    x_plot  = CHR - 1 + rel_pos
  )

# alternating row bands
row_bands <- tibble(
  ynum = seq_along(y_levels),
  fill_band = rep(c("white", "grey98"), length.out = length(y_levels))
)

# split categories per SNP into columns for scatterpie
all_cats <- names(cat_colors)

pie_data <- df_b_plot %>%
  mutate(Categories = as.character(Categories)) %>%
  rowwise() %>%
  mutate(
    cats_list = list(trimws(strsplit(Categories, ",")[[1]]))
  ) %>%
  ungroup()

for (cat in all_cats) {
  pie_data[[cat]] <- mapply(
    function(cats_list, n_cat) {
      if (cat %in% cats_list) 1 / n_cat else 0
    },
    pie_data$cats_list,
    pie_data$n_categories
  )
}

# sensible radius now that x-axis is normalized
pie_data <- pie_data %>%
  mutate(
    radius = scales::rescale(n_traits, to = c(0.10, 0.22))
  )

# ------------------------------------------
# spread out crowded local clusters within each chromosome/category row
# ------------------------------------------
pie_data <- pie_data %>%
  group_by(CHR, y_num) %>%
  arrange(x_plot, .by_group = TRUE) %>%
  mutate(
    gap = c(Inf, diff(x_plot)),
    cluster = cumsum(gap > 0.18)
  ) %>%
  group_by(CHR, y_num, cluster) %>%
  mutate(
    x_offset = if (n() == 1) 0 else seq(-0.15, 0.15, length.out = n()),
    x_plot_jitter = x_plot + x_offset
  ) %>%
  ungroup() %>%
  select(-gap, -cluster, -x_offset)

ring_data <- pie_data %>%
  filter(n_models_total == 4)

# dummy data only to generate a size legend below the plot
size_legend_data <- tibble(
  x_plot_jitter = NA_real_,
  y_num = NA_real_,
  n_traits = c(3, 5, 7, 9, 11)
)

p_b <- ggplot() +
  geom_rect(
    data = row_bands %>% filter(fill_band == "grey98"),
    aes(
      ymin = ynum - 0.5, ymax = ynum + 0.5,
      xmin = -Inf, xmax = Inf
    ),
    fill = "grey98",
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = 1:11,
    color = "grey85",
    linewidth = 0.3
  ) +
  geom_scatterpie(
    data = pie_data,
    aes(x = x_plot_jitter, y = y_num, r = radius),
    cols = all_cats,
    alpha = 0.9,
    color = "grey35",
    linewidth = 0.3
  ) +
  ggforce::geom_circle(
    data = ring_data,
    aes(x0 = x_plot_jitter, y0 = y_num, r = radius),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8,
    fill = NA,
    show.legend = FALSE
  ) +
  # invisible layer to create a bottom size legend
  geom_point(
    data = size_legend_data,
    aes(x = x_plot_jitter, y = y_num, size = n_traits),
    alpha = 0,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = cat_colors,
    name   = "Trait category"
  ) +
  scale_size_continuous(
    range = c(3, 8),
    breaks = c(3, 5, 7, 9, 11),
    name = "Traits significant"
  ) +
  geom_text_repel(
    data = pie_data,
    aes(x = x_plot_jitter, y = y_num, label = SNP),
    size = 3,
    color = "grey30",
    box.padding = 0.35,
    point.padding = 0.25,
    force = 2,
    min.segment.length = 0.2,
    max.overlaps = Inf
  ) +
  scale_x_continuous(
    breaks = seq(0.5, 11.5, by = 1),
    labels = 1:12,
    limits = c(-0.15, 12.25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(y_levels),
    labels = y_levels,
    limits = c(0.5, length(y_levels) + 0.5),
    expand = c(0, 0)
  ) +
  coord_fixed(ratio = 0.9) +
  gwas_theme +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 21, size = 5)
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(alpha = 1, shape = 21, fill = "grey80", color = "grey35")
    )
  ) +
  labs(x = "Chromosome")

print(p_b)

ggsave(
  file.path(outdir, "Figure1B_GWAS_manhattan_pie_newdata.pdf"),
  p_b,
  width = 38,
  height = 8
)