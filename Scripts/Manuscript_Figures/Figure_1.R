#############
# FIGURE 1A
#############

#install.packages("circlize")
library(circlize)
library(readxl)
library(tidyverse)

# load data
df_a <- read_xlsx("~/R/World_Veg_Project/filtered_data/STI_GWAS_summary_updated/table_S5_348_SNPs.xlsx") %>%
  mutate(Category = factor(Category, levels = c(
    "Canopy spectral", "Leaf area", "Yield components",
    "Growth", "Physiology", "Phenology", "Pollen", "Yield")))

df_b <- read_xlsx("~/R/World_Veg_Project/filtered_data/STI_GWAS_summary_updated/table_s7_16_SNPs.xlsx")

df_special <- read_xlsx("~/R/World_Veg_Project/filtered_data/STI_GWAS_summary_updated/table_S5_348_SNPs.xlsx") %>%
  filter(Category %in% c("Yield", "Pollen")) %>%
  distinct(SNP, Chr, Pos, Category, Trait, n_models, Models) %>%
  mutate(`Max # of GWAS models` = n_models)

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

# chromosome lengths
chr_lengths <- c(333, 177, 290, 249, 255, 253, 266, 174, 278, 210, 275, 260) #Zhangshugang ref https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030867735.1/
chr_names   <- paste0("chr", 1:12)

# prepare genome reference
cytoband <- data.frame(
  V1 = chr_names,
  V2 = 0,
  V3 = chr_lengths * 1e6,
  V4 = chr_names,
  V5 = "gpos50"
)

# prepare SNP data for each category track
prep_track <- function(df, chr_col = "Chr", pos_col = "Pos") {
  df %>%
    mutate(
      chr   = paste0("chr", .data[[chr_col]]),
      start = .data[[pos_col]],
      end   = .data[[pos_col]] + 1
    )
}

df_a_track <- prep_track(df_a)
df_b_track <- prep_track(df_b, pos_col = "Position")
df_sp_track <- prep_track(df_special)

# open PDF
pdf("Figure1_GWAS_circos.pdf", width = 10, height = 10)

circos.clear()
circos.par(
  "track.height"    = 0.08,
  "cell.padding"    = c(0, 0, 0, 0),
  "gap.degree"      = 2,
  "start.degree"    = 90,
  "clock.wise"      = TRUE
)

# initialize with chromosome ideogram
circos.initializeWithIdeogram(
  cytoband,
  chromosome.index = chr_names,
  plotType         = c("axis", "labels"),
  axis.labels.cex  = 0.8,
  labels.cex       = 1.4
)

# one track per category for Panel A data
cats <- levels(df_a$Category)

for (cat in cats) {
  df_cat <- df_a_track %>% filter(Category == cat)
  col    <- cat_colors[cat]
  
  circos.track(
    ylim = c(0, 4),
    track.height = 0.07,
    bg.border    = "grey90",
    bg.col       = "grey99",
    panel.fun = function(region, value, ...) {
      chr_now <- CELL_META$sector.index
      df_chr  <- df_cat %>% filter(chr == chr_now)
      if (nrow(df_chr) == 0) return()
      circos.points(
        x    = df_chr$start,
        y    = df_chr$n_models,
        pch  = 16,
        cex  = df_chr$n_models * 0.3,
        col  = adjustcolor(col, alpha.f = 0.7)
      )
    }
  )
}

# innermost track — 16 gold SNPs + pollen/yield diamonds
circos.track(
  ylim = c(0, 12),
  track.height = 0.12,
  bg.border    = "grey85",
  bg.col       = "white",
  panel.fun = function(region, value, ...) {
    chr_now <- CELL_META$sector.index
    
    # gold SNPs
    df_chr <- df_b_track %>% filter(chr == chr_now)
    if (nrow(df_chr) > 0) {
      circos.points(
        x   = df_chr$Position,
        y   = df_chr$Traits_with_SNP,
        pch = 21,
        cex = df_chr$Traits_with_SNP * 0.15,
        col = "black",
        bg  = cat_colors[as.character(df_chr$Category)]
      )
    }
    
    # pollen/yield diamonds
    df_sp <- df_sp_track %>% filter(chr == chr_now)
    if (nrow(df_sp) > 0) {
      circos.points(
        x   = df_sp$Pos,
        y   = rep(2, nrow(df_sp)),
        pch = 23,
        cex = 1.2,
        col = "black",
        bg  = cat_colors[as.character(df_sp$Category)]
      )
    }
  }
)

# add category labels on the right
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
getwd()

####################################################################
# Inner track with chromosome labels and Mb axis
# Uses your existing cytoband, chr_names, df_b_track,
# df_sp_track and cat_colors objects
####################################################################

pdf("Figure1_GWAS_inner_track.pdf", width = 10, height = 10)

circos.clear()

circos.par(
  cell.padding = c(0, 0, 0, 0),
  gap.degree = 2,
  start.degree = 90,
  clock.wise = TRUE
)

# Chromosome labels and Mb numbers, matching the original
circos.initializeWithIdeogram(
  cytoband,
  chromosome.index = chr_names,
  plotType = c("axis", "labels"),
  axis.labels.cex = 0.8,
  labels.cex = 1.4
)

# Single SNP track
circos.trackPlotRegion(
  ylim = c(0, 12),
  track.height = 0.25,
  bg.border = "grey70",
  bg.col = "grey99",
  panel.fun = function(x, y) {
    
    chr_now <- CELL_META$sector.index
    
    # Multi-trait SNPs: circles
    df_chr <- df_b_track %>%
      dplyr::filter(chr == chr_now)
    
    if (nrow(df_chr) > 0) {
      circos.points(
        x = df_chr$start,
        y = df_chr$Traits_with_SNP,
        pch = 21,
        cex = df_chr$Traits_with_SNP * 0.15,
        col = "black",
        bg = unname(
          cat_colors[as.character(df_chr$Category)]
        )
      )
    }
    
    # Pollen/yield SNPs: diamonds
    df_sp <- df_sp_track %>%
      dplyr::filter(chr == chr_now)
    
    if (nrow(df_sp) > 0) {
      circos.points(
        x = df_sp$start,
        y = rep(2, nrow(df_sp)),
        pch = 23,
        cex = 1.2,
        col = "black",
        bg = unname(
          cat_colors[as.character(df_sp$Category)]
        )
      )
    }
  }
)

# Category legend in the empty centre
legend(
  "center",
  legend = names(cat_colors),
  fill = unname(cat_colors),
  border = NA,
  cex = 0.85,
  title = "Trait category",
  bty = "n"
)

circos.clear()
dev.off()


#############
# FIGURE 1B
#############
library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)

# shared setup
chr_lengths <- c(333, 177, 290, 249, 255, 253, 266, 174, 278, 210, 275, 260) #Zhangshugang ref https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030867735.1/
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

gwas_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    axis.text.x         = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y         = element_text(size = 10, face = "bold"),
    axis.title.y        = element_blank(),
    legend.position     = "bottom",
    legend.title        = element_text(size = 9, face = "bold"),
    legend.text         = element_text(size = 8),
    plot.title          = element_blank(),
    plot.subtitle       = element_blank(),
    panel.border        = element_rect(color = "grey70", linewidth = 0.5),
    panel.background    = element_rect(fill = "white")
  )

cat_levels_order <- c(
  "Canopy spectral", "Leaf area", "Yield components",
  "Growth", "Physiology", "Phenology", "Pollen", "Yield"
)

row_bands <- tibble(
  Category = factor(cat_levels_order, levels = cat_levels_order),
  ymin = seq_along(cat_levels_order) - 0.5,
  ymax = seq_along(cat_levels_order) + 0.5,
  fill = rep(c("white", "grey97"), length.out = length(cat_levels_order))
)
############################################
# Panel B — 16 gold SNPs + yield/pollen
############################################

df_b <- read_xlsx("~/R/World_Veg_Project/filtered_data/STI_GWAS_summary_updated/table_s7_16_SNPs.xlsx") %>%
  mutate(
    x_pos    = chr_offsets[Chr] + Position / 1e6,
    Category = factor(Category, levels = cat_levels_order)
  )

df_special <- read_xlsx("~/R/World_Veg_Project/filtered_data/STI_GWAS_summary_updated/table_S5_348_SNPs.xlsx") %>%
  filter(Category %in% c("Yield", "Pollen")) %>%
  distinct(SNP, Chr, Pos, Category, Trait, n_models, Models) %>%
  mutate(
    x_pos                  = chr_offsets[Chr] + Pos / 1e6,
    Traits_with_SNP        = 1,
    Category               = factor(Category, levels = cat_levels_order),
    `Max # of GWAS models` = n_models
  )


p_b <- ggplot(df_b, aes(x = x_pos, y = Category,
                        size = Traits_with_SNP,
                        fill = Category, color = Category)) +
  geom_rect(data = row_bands,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = fill),
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c(setNames(cat_colors, names(cat_colors)),
                               "white" = "white", "grey97" = "grey97"),
                    guide = "none") +
  geom_vline(xintercept = chr_offsets[-1] - 7.5, color = "grey85", linewidth = 0.3) +
  geom_point(shape = 21, alpha = 0.85, stroke = 0.6) +
  geom_point(data = df_b %>% filter(`Max # of GWAS models` == 4),
             shape = 21, fill = NA, color = "black", stroke = 1.5, show.legend = FALSE) +
  geom_point(data = df_special,
             aes(x = x_pos, y = Category, fill = Category, color = Category),
             shape = 23, size = 4, alpha = 0.9, stroke = 0.8, show.legend = FALSE) +
  geom_point(data = df_special %>% filter(`Max # of GWAS models` == 4),
             aes(x = x_pos, y = Category),
             shape = 23, fill = NA, color = "black", size = 4, stroke = 1.5, show.legend = FALSE) +
  geom_text_repel(aes(label = SNP), size = 3.5, color = "grey30",
                  box.padding = 0.4, min.segment.length = 0.2, max.overlaps = 20) +
  geom_text_repel(data = df_special, aes(x = x_pos, y = Category, label = SNP),
                  size = 3.5, color = "grey30",
                  box.padding = 0.4, min.segment.length = 0.2, max.overlaps = 20) +
  scale_color_manual(values = cat_colors, guide = "none") +
  scale_size_area(max_size = 16, breaks = c(3, 5, 7, 9, 11), name = "Traits significant") +
  scale_x_continuous(breaks = chr_mids, labels = 1:12, expand = c(0.02, 0)) +
  scale_y_discrete(limits = rev) +
  gwas_theme +
  theme(
    base_size      = 14,
    axis.text.x    = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y    = element_text(size = 12, face = "bold"),
    axis.title.x   = element_text(size = 13),
    legend.title   = element_text(size = 11, face = "bold"),
    legend.text    = element_text(size = 10)
  ) +
  guides(size = guide_legend(override.aes = list(fill = "grey50", color = "grey50", shape = 21))) +
  annotate("text", x = Inf, y = -Inf,
           label = "◆ Yield / Pollen (multi-model, single traits category)",
           hjust = 1.05, vjust = -0.5, size = 3.5, color = "grey30") +
  labs(x = "Chromosome")
p_b
#ggsave("Figure1B_GWAS.pdf", p_b, width = 20, height = 7)







