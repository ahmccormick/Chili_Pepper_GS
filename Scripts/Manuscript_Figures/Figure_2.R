###################################################################################################
# FIGURE 2 - Core vs Global GEBV density + category-level PA summary
# Panel A: Density of mean scaled GEBVs by category
# Panel B: Mean rrBLUP prediction accuracy (PA) by category
###################################################################################################

library(tidyverse)

#############################################
# 0) File paths
#############################################

core_file   <- "~/R/WorldVeg_Capsicum/STI_GS/GEBVs_STI_73traits_BLUEadjusted_ALL_n423.csv"
global_file <- "~/R/WorldVeg_Capsicum/STI_GS_10k/GEBVs_STI_73traits_BLUEadjusted_ALL_10k.csv"

core_pa_file   <- "~/R/WorldVeg_Capsicum/STI_GS/PA_summary_73traits_n423.csv"
global_pa_file <- "~/R/WorldVeg_Capsicum/STI_GS_10k/PA_summary_all_traits_10k.csv"

out_base <- "~/R/WorldVeg_Capsicum/Figures/"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

#############################################
# 1) Shared color palette
#############################################

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

#############################################
# 2) Category lookup
#############################################

cat_lookup_sti <- tribble(
  ~Trait,                        ~Category,
  "Huedaily",                    "Canopy spectral",
  "Huemean",                     "Canopy spectral",
  "Huebin0",                     "Canopy spectral",
  "Huebin1",                     "Canopy spectral",
  "Huebin5",                     "Canopy spectral",
  "NDVImean",                    "Canopy spectral",
  "NDVIbin0",                    "Canopy spectral",
  "NDVIbin1",                    "Canopy spectral",
  "NDVIbin2",                    "Canopy spectral",
  "NDVIbin3",                    "Canopy spectral",
  "NDVIbin4",                    "Canopy spectral",
  "NDVIbin5",                    "Canopy spectral",
  "Green",                       "Canopy spectral",
  "Greenb0",                     "Canopy spectral",
  "Greenb1",                     "Canopy spectral",
  "Greenb2",                     "Canopy spectral",
  "Greenb3",                     "Canopy spectral",
  "NPCImean",                    "Canopy spectral",
  "PSRImean",                    "Canopy spectral",
  
  "Growthrateplot",              "Growth",
  "Growthrateplant",             "Growth",
  "Height",                      "Growth",
  "Heightmax",                   "Growth",
  "Biomassfinalplot",            "Growth",
  "Biomassfinalplant",           "Growth",
  "Heightfinal",                 "Growth",
  "Heighmaxfinal",               "Growth",
  
  "Leafareadaily",               "Leaf area",
  "Leafareadailyplant",          "Leaf area",
  "Leafareaindex",               "Leaf area",
  "Leafareaprojecteddaily",      "Leaf area",
  "Leafareaprojectedplant",      "Leaf area",
  "Leafareafinalplot",           "Leaf area",
  "Leafareafinalplant",          "Leaf area",
  "Leafareaindexfinal",          "Leaf area",
  "Leafareafinalprojectedplot",  "Leaf area",
  "Leafareafinalprojectedplant", "Leaf area",
  
  "DASanthesis",                 "Phenology",
  "DASmaturity",                 "Phenology",
  "DATanthesis",                 "Phenology",
  "DATmaturity",                 "Phenology",
  
  "Leaftemp",                    "Physiology",
  "Leafminushobotemp",           "Physiology",
  "Leafangledaily",              "Physiology",
  "Leafanglemean",               "Physiology",
  "Leafanglenoon",               "Physiology",
  "Leafanglemidnight",           "Physiology",
  "Leafanglemidnightminusnoon",  "Physiology",
  "Leafanglebright",             "Physiology",
  "Leafanglecold",               "Physiology",
  "Leafanglecoldminushot",       "Physiology",
  "LeafangleHot",                "Physiology",
  "Leafanglelow",                "Physiology",
  "Leafanglelowminusbright",     "Physiology",
  "Leafinclinationbright",       "Physiology",
  "Leafinclinationcold",         "Physiology",
  "Leafinclinationdark",         "Physiology",
  "Leafinclinationhot",          "Physiology",
  "Leafinclinationmean",         "Physiology",
  "Leafinclinationmidnight",     "Physiology",
  "Leafinclinationnoon",         "Physiology",
  "Lightpendaily",               "Physiology",
  "Lightpenetrationfinal",       "Physiology",
  "Lightpenmeanplant",           "Physiology",
  "Lightpenmeanplot",            "Physiology",
  
  "Pollenactivity",              "Pollen",
  "Pollenconcentration",         "Pollen",
  
  "Yield",                       "Yield",
  
  "Fruitno",                     "Yield components",
  "Fruitlength",                 "Yield components",
  "Fruitshapeindex",             "Yield components",
  "Fruitweight",                 "Yield components",
  "Fruitwidth",                  "Yield components"
)

#############################################
# 3) Number of traits per category
#############################################

cat_n_traits <- cat_lookup_sti %>%
  group_by(Category) %>%
  summarise(n_traits = n_distinct(Trait), .groups = "drop")

#############################################
# 4) Helper functions
#############################################

clean_trait <- function(x) {
  x %>%
    tolower() %>%
    gsub("\\.", "", .) %>%
    gsub("_", "", .)
}

process_gebv_sti <- function(df, collection_name) {
  df %>%
    pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
    mutate(
      Trait = gsub("^GEBV_", "", Trait_raw),
      Trait_match = clean_trait(Trait)
    ) %>%
    left_join(
      cat_lookup_sti %>%
        mutate(Trait_match = clean_trait(Trait)) %>%
        select(Trait_match, Category),
      by = "Trait_match"
    ) %>%
    filter(!is.na(Category)) %>%
    group_by(Trait) %>%
    mutate(GEBV_scaled = as.numeric(scale(GEBV))) %>%
    ungroup() %>%
    group_by(Line, Category) %>%
    summarise(
      mean_scaled_GEBV = mean(GEBV_scaled, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Collection = collection_name)
}

make_labels_with_pa <- function(df, pa_summary) {
  
  base_labels <- cat_lookup_sti %>%
    group_by(Category) %>%
    summarise(n_traits = n_distinct(Trait), .groups = "drop") %>%
    left_join(pa_summary, by = "Category") %>%
    mutate(
      Category_label = paste0(
        Category,
        "\n(n = ", n_traits, " traits)",
        "\nCore PA = ", Core_PA_label,
        ", Global PA = ", Global_PA_label
      )
    )
  
  label_levels <- base_labels %>%
    mutate(Category = factor(Category, levels = names(cat_colors))) %>%
    arrange(Category) %>%
    pull(Category_label)
  
  df %>%
    left_join(base_labels %>% select(Category, Category_label), by = "Category") %>%
    mutate(Category_label = factor(Category_label, levels = label_levels))
}

density_theme <- theme_bw(base_size = 11) +
  theme(
    strip.text         = element_text(face = "bold", size = 8, lineheight = 1.0),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9)
  )

#############################################
# 5) Read rrBLUP PA tables
#############################################

pa_core <- readr::read_csv(core_pa_file, show_col_types = FALSE) %>%
  filter(!is.na(trait), trait != "") %>%
  transmute(
    Trait_raw = trait,
    PA_core   = as.numeric(r.mean_rrBLUP),
    SD_core   = as.numeric(r.sd_rrBLUP)
  )

pa_global <- readr::read_csv(global_pa_file, show_col_types = FALSE) %>%
  filter(!is.na(trait), trait != "") %>%
  transmute(
    Trait_raw  = trait,
    PA_global  = as.numeric(r.mean_rrBLUP),
    SD_global  = as.numeric(r.sd_rrBLUP)
  )

#############################################
# 6) Harmonize PA trait names to lookup table
#############################################

cat_lookup_clean <- cat_lookup_sti %>%
  mutate(Trait_match = clean_trait(Trait)) %>%
  mutate(
    Trait_match = case_when(
      Trait_match == "heightmaxfinal" ~ "heighmaxfinal",
      TRUE ~ Trait_match
    )
  )

pa_core_cat <- pa_core %>%
  mutate(Trait_match = clean_trait(Trait_raw)) %>%
  left_join(cat_lookup_clean %>% select(Trait_match, Category), by = "Trait_match")

pa_global_cat <- pa_global %>%
  mutate(Trait_match = clean_trait(Trait_raw)) %>%
  left_join(cat_lookup_clean %>% select(Trait_match, Category), by = "Trait_match")

message("Unmatched core PA traits:")
print(pa_core_cat %>% filter(is.na(Category)))

message("Unmatched global PA traits:")
print(pa_global_cat %>% filter(is.na(Category)))

#############################################
# 7) Category-level PA summary
#############################################

pa_category_core <- pa_core_cat %>%
  filter(!is.na(Category)) %>%
  group_by(Category) %>%
  summarise(
    mean_PA_core = mean(PA_core, na.rm = TRUE),
    sd_PA_core   = sd(PA_core, na.rm = TRUE),
    n_traits_core = n(),
    .groups = "drop"
  )

pa_category_global <- pa_global_cat %>%
  filter(!is.na(Category)) %>%
  group_by(Category) %>%
  summarise(
    mean_PA_global = mean(PA_global, na.rm = TRUE),
    sd_PA_global   = sd(PA_global, na.rm = TRUE),
    n_traits_global = n(),
    .groups = "drop"
  )

pa_category_both <- full_join(pa_category_core, pa_category_global, by = "Category") %>%
  mutate(
    Core_PA_label   = sprintf("%.2f", mean_PA_core),
    Global_PA_label = sprintf("%.2f", mean_PA_global)
  )

#############################################
# 8) Read GEBV files
#############################################

gebv_core <- read.csv(core_file, stringsAsFactors = FALSE) %>%
  rename(Line = 1) %>%
  select(-any_of("Set"))

gebv_global <- read.csv(global_file, stringsAsFactors = FALSE) %>%
  rename(Line = 1) %>%
  select(-any_of("Set"))

message("Core lines: ", nrow(gebv_core))
message("Global lines: ", nrow(gebv_global))


#############################################
# 9) Process GEBVs and combine - scaled TOGETHER #############################################
#############################################

# combine both collections first, then scale
gebv_combined <- bind_rows(
  gebv_core   %>% mutate(Collection = "Core (n = 423)"),
  gebv_global %>% mutate(Collection = "Global (n = 10,026)")
) %>%
  pivot_longer(-c(Line, Collection), names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(
    Trait       = gsub("^GEBV_", "", Trait_raw),
    Trait_match = clean_trait(Trait)
  ) %>%
  left_join(
    cat_lookup_sti %>%
      mutate(Trait_match = clean_trait(Trait)) %>%
      select(Trait_match, Category),
    by = "Trait_match"
  ) %>%
  filter(!is.na(Category)) %>%
  group_by(Trait) %>%                          # scale across BOTH collections together
  mutate(GEBV_scaled = as.numeric(scale(GEBV))) %>%
  ungroup() %>%
  group_by(Line, Category, Collection) %>%
  summarise(
    mean_scaled_GEBV = mean(GEBV_scaled, na.rm = TRUE),
    .groups = "drop"
  )

print(gebv_combined %>% count(Collection))

gebv_combined_labeled <- make_labels_with_pa(
  gebv_combined,
  pa_summary = pa_category_both %>% select(Category, Core_PA_label, Global_PA_label)
)


#############################################
# 10) Panel A - Density plot
#############################################

p_sti <- gebv_combined_labeled %>%
  mutate(
    Category   = factor(Category, levels = names(cat_colors)),
    Collection = factor(Collection, levels = c("Core (n = 423)", "Global (n = 10,026)"))
  ) %>%
  ggplot(aes(x = mean_scaled_GEBV)) +
  geom_density(
    data = . %>% filter(Collection == "Global (n = 10,026)"),
    aes(fill = "Global (n = 10,026)", color = "Global (n = 10,026)"),
    alpha = 0.35, linewidth = 0.6
  ) +
  geom_density(
    data = . %>% filter(Collection == "Core (n = 423)"),
    aes(fill = Category, color = Category),
    alpha = 0.45, linewidth = 0.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  coord_cartesian(xlim = c(-3, 3)) +
  scale_fill_manual(
    values = c(cat_colors, "Global (n = 10,026)" = "grey70"),
    breaks = c("Global (n = 10,026)"),
    name   = NULL
  ) +
  scale_color_manual(
    values = c(cat_colors, "Global (n = 10,026)" = "grey50"),
    breaks = c("Global (n = 10,024)"),
    name   = NULL
  ) +
  facet_wrap(~ Category_label, scales = "free_y", ncol = 2) +
  density_theme +
  labs(
    x = "Mean scaled GEBV (z-score)",
    y = "Density",
    title = NULL,
    subtitle = NULL
  )

print(p_sti)

out_base <- "~/R/WorldVeg_Capsicum/Figures/"
ggsave(
  file.path(out_base, "Figure2_GEBV_core_vs_global_BLUEadjusted.pdf"),
  p_sti,
  width  = 10,
  height = 6
)


######################################################################################################################################################################################################
############## curious about singles
###################################################################################################
# Supplemental plot: all 73 traits, core vs global distributions with rrBLUP PA
###################################################################################################

# build per-trait PA labels from already-loaded pa_core and pa_global
pa_for_labels <- full_join(
  pa_core %>% transmute(Trait_match = clean_trait(Trait_raw), PA_core),
  pa_global %>% transmute(Trait_match = clean_trait(Trait_raw), PA_global),
  by = "Trait_match"
) %>%
  mutate(
    PA_core_label   = ifelse(is.na(PA_core),   "NA", sprintf("%.2f", PA_core)),
    PA_global_label = ifelse(is.na(PA_global), "NA", sprintf("%.2f", PA_global))
  )

# long format for both collections
gebv_core_long <- gebv_core %>%
  pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(
    Trait      = gsub("^GEBV_", "", Trait_raw),
    Collection = "Core (n = 423)"
  )

gebv_global_long <- gebv_global %>%
  pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(
    Trait      = gsub("^GEBV_", "", Trait_raw),
    Collection = "Global (n = 10,026)"
  )

# scale together
gebv_trait_long <- bind_rows(gebv_core_long, gebv_global_long) %>%
  group_by(Trait) %>%
  mutate(GEBV_scaled = as.numeric(scale(GEBV))) %>%
  ungroup()

#this is the second one
# scale separately
#gebv_trait_long <- bind_rows(gebv_core_long, gebv_global_long) %>%
  #group_by(Trait, Collection) %>%        # scale separately
 # mutate(GEBV_scaled = as.numeric(scale(GEBV))) %>%
  #ungroup()


# add category and PA labels
trait_category_lookup <- cat_lookup_sti %>%
  mutate(Trait_match = clean_trait(Trait))

gebv_trait_long <- gebv_trait_long %>%
  mutate(Trait_match = clean_trait(Trait)) %>%
  left_join(
    trait_category_lookup %>% select(Trait_match, Category),
    by = "Trait_match"
  ) %>%
  left_join(pa_for_labels, by = "Trait_match") %>%
  mutate(
    facet_label = ifelse(
      is.na(Category),
      Trait,
      paste0(Trait, "\n(", Category, ")",
             "\nCore PA = ", PA_core_label,
             ", Global PA = ", PA_global_label)
    )
  )

# keep trait order from lookup
trait_order <- cat_lookup_sti %>%
  pull(Trait) %>%
  unique()

gebv_trait_long <- gebv_trait_long %>%
  mutate(
    Trait = factor(Trait, levels = trait_order),
    facet_label = factor(
      facet_label,
      levels = gebv_trait_long %>%
        distinct(Trait, facet_label) %>%
        arrange(Trait) %>%
        pull(facet_label)
    ),
    Collection = factor(Collection,
                        levels = c("Global (n = 10,026)", "Core (n = 423)"))
  )

p_all_traits <- ggplot(gebv_trait_long,
                       aes(x = GEBV_scaled, fill = Collection, color = Collection)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.25) +
  scale_fill_manual(values = c("Global (n = 10,026)" = "grey70",
                               "Core (n = 423)"      = "#4E79A7")) +
  scale_color_manual(values = c("Global (n = 10,026)" = "grey40",
                                "Core (n = 423)"      = "#2F5D8A")) +
  coord_cartesian(xlim = c(-5, 5)) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 9) +
  theme(
    strip.text         = element_text(face = "bold", size = 6, lineheight = 0.9),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 8)
  ) +
  labs(
    x = "Scaled GEBV (z-score)",
    y = "Density",
    title = NULL,
    subtitle = NULL
  )

print(p_all_traits)




######################################################################################################################################################################################################
# 16 overlap traits
# Core vs Global for overlapping high-PA traits only (individual trait facets)
# scaled together
###################################################################################################

library(tidyverse)

#############################################
# 0) File paths
#############################################

core_file      <- "~/R/WorldVeg_Capsicum/STI_GS/GEBVs_STI_73traits_BLUEadjusted_ALL_n423.csv"
global_file    <- "~/R/WorldVeg_Capsicum/STI_GS_10k/GEBVs_STI_73traits_BLUEadjusted_ALL_10k.csv"
core_pa_file   <- "~/R/WorldVeg_Capsicum/STI_GS/PA_summary_73traits_n423.csv"
global_pa_file <- "~/R/WorldVeg_Capsicum/STI_GS_10k/PA_summary_all_traits_10k.csv"
out_base       <- "~/R/WorldVeg_Capsicum/Figures/"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

#############################################
# 1) Load GEBV files
#############################################

gebv_core <- read.csv(core_file, stringsAsFactors = FALSE) %>%
  rename(Line = 1) %>%
  select(-any_of("Set"))

gebv_global <- read.csv(global_file, stringsAsFactors = FALSE) %>%
  rename(Line = 1) %>%
  select(-any_of("Set"))

#############################################
# 2) Load PA files and extract rrBLUP PA
#############################################

clean_trait <- function(x) {
  x %>% tolower() %>% gsub("\\.", "", .) %>% gsub("_", "", .)
}

pa_core <- readr::read_csv(core_pa_file, show_col_types = FALSE) %>%
  filter(!is.na(trait), trait != "") %>%
  transmute(
    Trait      = clean_trait(trait),
    PA_core    = as.numeric(r.mean_rrBLUP)
  )

pa_global <- readr::read_csv(global_pa_file, show_col_types = FALSE) %>%
  filter(!is.na(trait), trait != "") %>%
  transmute(
    Trait      = clean_trait(trait),
    PA_global  = as.numeric(r.mean_rrBLUP)
  )

pa_both <- full_join(pa_core, pa_global, by = "Trait")

#############################################
# 3) Define shared high-PA traits
#############################################

shared_traits <- c(
  "biomassfinalplant",
  "fruitshapeindex",
  "fruitweight",
  "fruitwidth",
  "growthrateplant",
  "leafareafinalprojectedplant",
  "leafareaprojectedplant",
  "ndvibin5",
  "ndvimean",
  "yield"
)


shared_traits <- c(
  "biomassfinalplant",
  "fruitshapeindex",
  "fruitweight",
  "fruitwidth",
  "growthrateplant",
  "leafareafinalprojectedplant",
  "leafareaprojectedplant",
  "ndvibin5",
  "ndvimean",
  "yield",
  "dasanthesis",
  "dasmaturity",
  "datanthesis",
  "datmaturity"
)
#############################################
# 4) Trait full name labels
#############################################

trait_info <- tribble(
  ~Trait,                        ~Full_label,
  "biomassfinalplant",           "Biomass final per plant",
  "fruitshapeindex",             "Fruit shape index",
  "fruitweight",                 "Average fruit weight",
  "fruitwidth",                  "Average fruit width",
  "growthrateplant",             "Digital biomass daily change per plant",
  "leafareafinalprojectedplant", "Leaf area projected final per plant",
  "leafareaprojectedplant",      "Leaf area (projected) daily change per plant",
  "ndvibin5",                    "NDVI bin 5",
  "ndvimean",                    "NDVI average",
  "yield",                       "Fruit yield"
)

trait_info <- tribble(
  ~Trait,                        ~Full_label,
  "biomassfinalplant",           "Biomass final per plant",
  "fruitshapeindex",             "Fruit shape index",
  "fruitweight",                 "Average fruit weight",
  "fruitwidth",                  "Average fruit width",
  "growthrateplant",             "Digital biomass daily change per plant",
  "leafareafinalprojectedplant", "Leaf area projected final per plant",
  "leafareaprojectedplant",      "Leaf area (projected) daily change per plant",
  "ndvibin5",                    "NDVI bin 5",
  "ndvimean",                    "NDVI average",
  "yield",                       "Fruit yield",
  "dasanthesis",                 "Days after sowing to anthesis",
  "dasmaturity",                 "Days after sowing to maturity",
  "datanthesis",                 "Days after transplant to anthesis",
  "datmaturity",                 "Days after transplant to maturity"
)
#############################################
# 5) Pivot GEBVs to long format
#############################################

gebv_core_long <- gebv_core %>%
  pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(Trait = clean_trait(gsub("^GEBV_", "", Trait_raw)),
         Collection = "Core (n = 423)")

gebv_global_long <- gebv_global %>%
  pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(Trait = clean_trait(gsub("^GEBV_", "", Trait_raw)),
         Collection = "Global (n = 10,026)")

gebv_shared <- bind_rows(gebv_core_long, gebv_global_long) %>%
  filter(Trait %in% shared_traits) %>%
  group_by(Trait) %>%
  mutate(GEBV_scaled = as.numeric(scale(GEBV))) %>%
  ungroup()

#############################################
# 6) Add PA labels and facet labels
#############################################

pa_shared <- pa_both %>%
  filter(Trait %in% shared_traits) %>%
  mutate(
    PA_core_label   = sprintf("%.2f", PA_core),
    PA_global_label = sprintf("%.2f", PA_global)
  )

gebv_shared <- gebv_shared %>%
  left_join(trait_info, by = "Trait") %>%
  left_join(pa_shared %>% select(Trait, PA_core_label, PA_global_label), by = "Trait") %>%
  mutate(
    facet_label = paste0(
      Full_label,
      "\n(", Trait, ")",
      "\nCore PA = ", PA_core_label,
      ", Global PA = ", PA_global_label
    )
  )

# set factor levels for trait order
facet_levels <- gebv_shared %>%
  distinct(Trait, facet_label) %>%
  mutate(Trait = factor(Trait, levels = shared_traits)) %>%
  arrange(Trait) %>%
  pull(facet_label)

gebv_shared <- gebv_shared %>%
  mutate(
    Trait       = factor(Trait, levels = shared_traits),
    facet_label = factor(facet_label, levels = facet_levels),
    Collection  = factor(Collection, levels = c("Global (n = 10,026)", "Core (n = 423)"))
  )

#############################################
# 7) Plot
#############################################

p_shared_traits <- ggplot(gebv_shared, aes(x = GEBV_scaled)) +
  geom_density(
    data = subset(gebv_shared, Collection == "Global (n = 10,026)"),
    aes(fill = "Global (n = 10,026)", color = "Global (n = 10,026)"),
    alpha = 0.35, linewidth = 0.6
  ) +
  geom_density(
    data = subset(gebv_shared, Collection == "Core (n = 423)"),
    aes(fill = "Core (n = 423)", color = "Core (n = 423)"),
    alpha = 0.35, linewidth = 0.6
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  coord_cartesian(xlim = c(-5, 5)) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c("Global (n = 10,026)" = "grey70",
               "Core (n = 423)"      = "#4E79A7"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Global (n = 10,026)" = "grey45",
               "Core (n = 423)"      = "#2F5D8A"),
    name = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(face = "bold", size = 7, lineheight = 0.95),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9)
  ) +
  labs(
    x = "Scaled GEBV (z-score)",
    y = "Density",
    title = NULL,
    subtitle = NULL
  )

print(p_shared_traits)
