####################################################################
# Figure 3: Quality GEBV distributions
# Core n = 397; global n = 7,843
# Z-score scaling across both collections combined, within each trait
####################################################################

library(tidyverse)

####################################################################
# 1. Paths
####################################################################

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/quality"

core_folder <- file.path(folder, "core")
global_folder <- file.path(folder, "global")

####################################################################
# 2. Load quality GEBVs
####################################################################

load_quality <- function(path, collection) {
  
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::select(-any_of(c("Set", "Group"))) %>%
    tidyr::pivot_longer(
      cols = -Sample_ID,
      names_to = "Trait",
      values_to = "GEBV"
    ) %>%
    dplyr::mutate(Collection = collection) %>%
    dplyr::filter(is.finite(GEBV))
}


gebv_q_long <- bind_rows(
  load_quality(
    file.path(core_folder, "GEBVs_Quality_all_397.csv"),
    "Core"
  ),
  load_quality(
    file.path(global_folder, "GEBVs_Quality_all_7843.csv"),
    "Global"
  )
)

####################################################################
# 3. Scale within each trait across both collections together
####################################################################

gebv_q_scaled <- gebv_q_long %>%
  group_by(Trait) %>%
  mutate(
    GEBV_scaled = if (sd(GEBV) > 0) {
      (GEBV - mean(GEBV)) / sd(GEBV)
    } else {
      NA_real_
    }
  ) %>%
  ungroup() %>%
  filter(is.finite(GEBV_scaled))

####################################################################
# 4. Read rrBLUP prediction accuracies
####################################################################

load_pa <- function(path, collection) {
  
  readRDS(path)$xval.result %>%
    transmute(
      Trait = as.character(trait),
      Collection = collection,
      PA = as.numeric(as.character(r.mean))
    ) %>%
    filter(!is.na(Trait), nzchar(Trait))
}

pa_quality <- bind_rows(
  load_pa(
    file.path(core_folder, "xval_rrblup_23traits_kfold_10.RDS"),
    "Core"
  ),
  load_pa(
    file.path(global_folder, "xval_rrblup_23traits_kfold_10.RDS"),
    "Global"
  )
)

stopifnot(
  !anyDuplicated(pa_quality[c("Trait", "Collection")])
)

####################################################################
# 5. Create facet labels and order by core PA
####################################################################

format_pa <- function(x) {
  ifelse(is.finite(x), sprintf("%.2f", x), "NA")
}

pa_wide <- pa_quality %>%
  pivot_wider(
    names_from = Collection,
    values_from = PA,
    names_prefix = "PA_"
  ) %>%
  mutate(
    facet_label = paste0(
      gsub("_", " ", Trait),
      "\n(PA core = ", format_pa(PA_Core),
      ", global = ", format_pa(PA_Global), ")"
    )
  ) %>%
  arrange(desc(PA_Core))

missing_pa <- setdiff(
  unique(gebv_q_scaled$Trait),
  pa_wide$Trait
)

if (length(missing_pa) > 0) {
  stop(
    "Traits missing from PA results: ",
    paste(missing_pa, collapse = ", ")
  )
}

gebv_q_plot <- gebv_q_scaled %>%
  dplyr::left_join(
    pa_wide %>% dplyr::select(Trait, facet_label),
    by = "Trait"
  ) %>%
  dplyr::mutate(
    facet_label = factor(
      facet_label,
      levels = pa_wide$facet_label
    ),
    Collection = factor(
      Collection,
      levels = c("Global", "Core")
    )
  )

####################################################################
# 6. Keep the same trait exclusions as your original figure
####################################################################

traits_to_remove <- c(
  "Fruit_predominant_shape_oblate",
  "Fruit_external_color_b",
  "Fruit_external_color_L",
  "Immature_fruit_external_color_L",
  "Fruit_external_color_a",
  "External_immature_fruit_color_green",
  "Immature_fruit_external_color_a",
  "Immature_fruit_external_color_b",
  "Fruit_load"
)

plot_data <- gebv_q_plot %>%
  filter(!Trait %in% traits_to_remove)

cat("Traits plotted:", n_distinct(plot_data$Trait), "\n")

####################################################################
####################################################################
# 7. Calculate means and quartiles
####################################################################

quality_summary <- plot_data %>%
  group_by(Trait, facet_label, Collection) %>%
  summarise(
    N = n(),
    Mean = mean(GEBV_scaled),
    Q25 = unname(quantile(GEBV_scaled, 0.25)),
    Median = median(GEBV_scaled),
    Q75 = unname(quantile(GEBV_scaled, 0.75)),
    .groups = "drop"
  )

quality_lines <- quality_summary %>%
  pivot_longer(
    cols = c(Mean, Q25, Median, Q75),
    names_to = "Statistic",
    values_to = "Position"
  ) %>%
  mutate(
    Line_type = case_when(
      Statistic == "Mean" ~ "Mean",
      Statistic == "Median" ~ "Median",
      TRUE ~ "25th / 75th percentile"
    ),
    Line_type = factor(
      Line_type,
      levels = c(
        "Mean",
        "Median",
        "25th / 75th percentile"
      )
    )
  )

####################################################################
# 8. Plot with population-specific means and quartiles
####################################################################

p_quality <- ggplot(
  plot_data,
  aes(
    x = GEBV_scaled,
    fill = Collection,
    color = Collection
  )
) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  
  geom_vline(
    data = quality_lines,
    aes(
      xintercept = Position,
      color = Collection,
      linetype = Line_type
    ),
    linewidth = 0.5,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  
  scale_fill_manual(
    values = c(
      Global = "grey60",
      Core = "#B07AA1"
    ),
    breaks = c("Global", "Core"),
    labels = c("Global (n = 7,843)", "Core (n = 397)"),
    name = NULL
  ) +
  
  scale_color_manual(
    values = c(
      Global = "grey40",
      Core = "#7B4F7B"
    ),
    breaks = c("Global", "Core"),
    labels = c("Global (n = 7,843)", "Core (n = 397)"),
    name = NULL
  ) +
  
  scale_linetype_manual(
    values = c(
      "Mean" = "solid",
      "Median" = "dashed",
      "25th / 75th percentile" = "dotted"
    ),
    name = NULL
  ) +
  
  coord_cartesian(xlim = c(-5, 5)) +
  
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 4
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = "none",
    linetype = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(color = "black")
    )
  ) +
  
  labs(
    x = "Scaled GEBV (z-score)",
    y = "Density",
    title = "Quality trait GEBV distributions — core vs global",
    subtitle = "Scaled within trait across both collections combined"
  )

print(p_quality)

####################################################################
# 9. Save at original dimensions
####################################################################

ggsave(
  filename = file.path(
    folder,
    "Figure3_Quality_GEBV_distributions.pdf"
  ),
  plot = p_quality,
  width = 14,
  height = 10,
  units = "in"
)

readr::write_csv(
  quality_summary,
  file.path(
    folder,
    "Figure3_Quality_population_mean_quantiles.csv"
  )
)


