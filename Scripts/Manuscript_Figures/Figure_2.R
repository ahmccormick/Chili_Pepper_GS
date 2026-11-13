####################################################################
# Figure 2: STI GEBV distributions by trait category
# Core n = 397; global n = 7,843
####################################################################

library(tidyverse)

####################################################################
# 1. Paths
####################################################################

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/STI"

core_folder <- file.path(folder, "core")
global_folder <- file.path(folder, "global")

core_pa_file <- file.path(
  core_folder,
  "xval_rrblup_STI_core397_kfold10.RDS"
)

global_pa_file <- file.path(
  global_folder,
  "xval_rrblup_STI_global7843_kfold10.RDS"
)

####################################################################
# 2. Category colours
####################################################################

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

####################################################################
# 3. Category lookup using new trait names
####################################################################

category_traits <- list(
  "Canopy spectral" = c(
    "huedaily", "huemean", "huebin0", "huebin1", "huebin5",
    "ndvimean", "ndvibin0", "ndvibin1", "ndvibin2",
    "ndvibin3", "ndvibin4", "ndvibin5",
    "green", "greenb0", "greenb1", "greenb2", "greenb3",
    "npcimean", "psrimean"
  ),
  "Leaf area" = c(
    "leafareadaily", "leafareadailyplant", "leafareaindex",
    "leafareaprojecteddaily", "leafareaprojectedplant",
    "leafareafinalplot", "leafareafinalplant",
    "leafareaindexfinal", "leafareafinalprojectedplot",
    "leafareafinalprojectedplant"
  ),
  "Yield components" = c(
    "fruitno", "fruitlength", "fruitwidth",
    "fruitshapeindex", "fruitweight"
  ),
  "Growth" = c(
    "growthrateplot", "growthrateplant",
    "height", "heightmax", "heightfinal", "heighmaxfinal",
    "biomassfinalplot", "biomassfinalplant"
  ),
  "Physiology" = c(
    "leaftemp", "leafminushobotemp",
    "leafangledaily", "leafanglemean", "leafanglenoon",
    "leafanglemidnight", "leafanglemidnightminusnoon",
    "leafanglebright", "leafanglelow", "leafanglelowminusbright",
    "leafanglecold", "leafanglehot", "leafanglecoldminushot",
    "leafinclinationmean", "leafinclinationhot",
    "leafinclinationcold", "leafinclinationbright",
    "leafinclinationdark", "leafinclinationnoon",
    "leafinclinationmidnight",
    "lightpendaily", "lightpenmeanplot",
    "lightpenmeanplant", "lightpenetrationfinal"
  ),
  "Phenology" = c(
    "dasanthesis", "dasmaturity", "datanthesis", "datmaturity"
  ),
  "Pollen" = c(
    "pollenconcentration", "pollenactivity"
  ),
  "Yield" = "yield"
)

cat_lookup <- tibble(
  Category = rep(names(category_traits), lengths(category_traits)),
  Trait_key = unlist(category_traits, use.names = FALSE)
)

stopifnot(!anyDuplicated(cat_lookup$Trait_key))

####################################################################
# 4. Read GEBVs and standardize trait names for matching
####################################################################

load_sti <- function(path, collection) {
  
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::select(-tidyselect::any_of(c("Set", "Group"))) %>%
    tidyr::pivot_longer(
      cols = -Sample_ID,
      names_to = "Trait",
      values_to = "GEBV"
    ) %>%
    dplyr::mutate(
      Trait_key = tolower(Trait),
      Collection = collection
    )
}

gebv_long <- dplyr::bind_rows(
  load_sti(
    file.path(core_folder, "GEBVs_STI_all_397.csv"),
    "Core"
  ),
  load_sti(
    file.path(global_folder, "GEBVs_STI_all_7843.csv"),
    "Global"
  )
)

unmapped_traits <- setdiff(
  unique(gebv_long$Trait_key),
  cat_lookup$Trait_key
)

if (length(unmapped_traits) > 0) {
  stop(
    "Traits missing from category lookup: ",
    paste(unmapped_traits, collapse = ", ")
  )
}

gebv_long <- gebv_long %>%
  dplyr::left_join(cat_lookup, by = "Trait_key")

####################################################################
# 5. Scale within trait AND collection, as in your previous code
####################################################################

gebv_scaled <- gebv_long %>%
  dplyr::filter(is.finite(GEBV)) %>%
  dplyr::group_by(Collection, Trait_key) %>%
  dplyr::mutate(
    GEBV_scaled = if (n() >= 2 && sd(GEBV) > 0) {
      (GEBV - mean(GEBV)) / sd(GEBV)
    } else {
      NA_real_
    }
  ) %>%
  dplyr::ungroup()

# Average standardized traits for each accession/category
gebv_combined <- gebv_scaled %>%
  dplyr::filter(is.finite(GEBV_scaled)) %>%
  dplyr::group_by(Collection, Sample_ID, Category) %>%
  dplyr::summarise(
    mean_scaled_GEBV = mean(GEBV_scaled),
    .groups = "drop"
  )

####################################################################
# 6. Load PA results, with a placeholder for missing core results
####################################################################

load_pa <- function(path, collection) {
  
  if (!file.exists(path)) {
    message(collection, " rrBLUP PA file not yet available.")
    
    return(tibble(
      Trait_key = character(),
      Collection = character(),
      PA = numeric()
    ))
  }
  
  readRDS(path)$xval.result %>%
    dplyr::transmute(
      Trait_key = tolower(as.character(trait)),
      Collection = collection,
      PA = as.numeric(as.character(r.mean))
    ) %>%
    dplyr::filter(
      !is.na(Trait_key),
      nzchar(Trait_key)
    )
}

pa_all <- dplyr::bind_rows(
  load_pa(core_pa_file, "Core"),
  load_pa(global_pa_file, "Global")
)

stopifnot(
  !anyDuplicated(pa_all[c("Trait_key", "Collection")])
)

# Category counts based on traits present in both GEBV files
shared_traits <- gebv_long %>%
  dplyr::distinct(Collection, Trait_key, Category) %>%
  dplyr::count(Trait_key, Category, name = "N_collections") %>%
  dplyr::filter(N_collections == 2)

cat_n_traits <- shared_traits %>%
  dplyr::count(Category, name = "n_traits")

pa_cat <- pa_all %>%
  dplyr::inner_join(
    shared_traits %>% dplyr::select(Trait_key, Category),
    by = "Trait_key"
  ) %>%
  dplyr::group_by(Category, Collection) %>%
  dplyr::summarise(
    n_PA = sum(is.finite(PA)),
    mean_PA = if (any(is.finite(PA))) {
      mean(PA[is.finite(PA)])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

####################################################################
# 7. Build facet labels
####################################################################

category_summary <- tidyr::expand_grid(
  Category = names(cat_colors),
  Collection = c("Core", "Global")
) %>%
  dplyr::left_join(cat_n_traits, by = "Category") %>%
  dplyr::left_join(
    pa_cat,
    by = c("Category", "Collection")
  ) %>%
  dplyr::mutate(
    PA_text = dplyr::case_when(
      Collection == "Core" & !file.exists(core_pa_file) ~ "pending",
      is.finite(mean_PA) ~ sprintf("%.2f", mean_PA),
      TRUE ~ "NA"
    )
  )

category_labels <- category_summary %>%
  dplyr::select(Category, n_traits, Collection, PA_text) %>%
  tidyr::pivot_wider(
    names_from = Collection,
    values_from = PA_text
  ) %>%
  dplyr::mutate(
    Category_label = paste0(
      Category,
      "\n(n = ", n_traits, " traits; mean PA core = ",
      Core, ", global = ", Global, ")"
    )
  ) %>%
  dplyr::arrange(match(Category, names(cat_colors)))

gebv_plot <- gebv_combined %>%
  dplyr::left_join(
    category_labels %>%
      dplyr::select(Category, Category_label),
    by = "Category"
  ) %>%
  dplyr::mutate(
    Category_label = factor(
      Category_label,
      levels = category_labels$Category_label
    )
  )




####################################################################
# 8. Calculate population means and quartiles for each category
####################################################################

distribution_summary <- gebv_plot %>%
  dplyr::filter(is.finite(mean_scaled_GEBV)) %>%
  dplyr::group_by(Category, Category_label, Collection) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Mean = mean(mean_scaled_GEBV),
    Q25 = unname(quantile(mean_scaled_GEBV, 0.25)),
    Median = median(mean_scaled_GEBV),
    Q75 = unname(quantile(mean_scaled_GEBV, 0.75)),
    .groups = "drop"
  )

# Prepare reference lines.
# Core lines use the category colour; global lines use grey.
reference_lines <- distribution_summary %>%
  tidyr::pivot_longer(
    cols = c(Mean, Q25, Median, Q75),
    names_to = "Statistic",
    values_to = "Position"
  ) %>%
  dplyr::mutate(
    Colour_group = ifelse(
      Collection == "Core",
      as.character(Category),
      "Global"
    ),
    Line_type = dplyr::case_when(
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
# 9. Plot density curves with means and quartiles
####################################################################

p_sti <- ggplot(gebv_plot, aes(x = mean_scaled_GEBV)) +
  
  # Global density
  geom_density(
    data = dplyr::filter(gebv_plot, Collection == "Global"),
    aes(fill = "Global", color = "Global"),
    alpha = 0.5,
    linewidth = 0.5
  ) +
  
  # Core density
  geom_density(
    data = dplyr::filter(gebv_plot, Collection == "Core"),
    aes(fill = Category, color = Category),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  # Population-specific means and quartiles
  geom_vline(
    data = reference_lines,
    aes(
      xintercept = Position,
      color = Colour_group,
      linetype = Line_type
    ),
    linewidth = 0.55,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  
  coord_cartesian(xlim = c(-3, 3)) +
  
  scale_fill_manual(
    values = c(cat_colors, Global = "grey70"),
    breaks = c(names(cat_colors), "Global"),
    labels = c(
      paste0(names(cat_colors), " — Core"),
      "Global (n = 7,843)"
    ),
    name = NULL
  ) +
  
  scale_color_manual(
    values = c(cat_colors, Global = "grey40"),
    breaks = c(names(cat_colors), "Global"),
    labels = c(
      paste0(names(cat_colors), " — Core"),
      "Global (n = 7,843)"
    ),
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
  
  facet_wrap(
    ~ Category_label,
    scales = "free_y",
    ncol = 2
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 8)
  ) +
  
  guides(
    fill = guide_legend(
      order = 1,
      ncol = 3
    ),
    color = "none",
    linetype = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(color = "black")
    )
  ) +
  
  labs(
    x = "Mean of trait-standardized GEBVs",
    y = "Density",
    subtitle = "Core (n = 397) vs Global (n = 7,843)"
  )

print(p_sti)

#####################################################################
# 10. Save at original dimensions
####################################################################

ggsave(
  filename = file.path(folder, "Figure2_STI_core_vs_global.pdf"),
  plot = p_sti,
  width = 14,
  height = 8,
  units = "in"
)

ggsave(
  filename = file.path(folder, "Figure2_STI_core_vs_global.png"),
  plot = p_sti,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300
)

readr::write_csv(
  category_summary,
  file.path(folder, "Figure2_STI_category_PA_summary.csv")
)

readr::write_csv(
  distribution_summary,
  file.path(folder, "Figure2_STI_population_mean_quantiles.csv")
)
