###########################################
#23 quality
library(tidyverse)

# read quality GEBVs
gebv_q_core   <- read.csv("~/R/World_Veg_Project/filtered_data/GEBVs_7_categories/GEBVs_quality_23trait_n423.csv",
                          stringsAsFactors = FALSE)
gebv_q_global <- read.csv("~/R/World_Veg_Project/filtered_data/GEBVs_7_categories/GEBVs_quality_23trait_n10026.csv",
                          stringsAsFactors = FALSE)

# check column names
colnames(gebv_q_core)

# read PA files
pa_423  <- read.csv("~/R/World_Veg_Project/filtered_data/GEBVs_7_categories/n423_PAs_96traits.csv",  stringsAsFactors = FALSE)
pa_10k  <- read.csv("~/R/World_Veg_Project/filtered_data/GEBVs_7_categories/n10k_PAs_96traits.csv",  stringsAsFactors = FALSE)

# filter to just the 23 quality traits — keep Title_Case names
quality_traits_423 <- pa_423 %>%
  filter(grepl("^[A-Z]", trait)) %>%   # quality traits start with capital letter
  select(Trait = trait, PA = r.mean)

quality_traits_10k <- pa_10k %>%
  filter(grepl("^[A-Z]", trait)) %>%
  select(Trait = trait, PA = r.mean)

print(quality_traits_423)
print(quality_traits_10k)


colnames(gebv_q_core)


#plot
library(tidyverse)

# pivot core to long
gebv_q_core_long <- gebv_q_core %>%
  select(-Group) %>%
  pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(Trait = gsub("^GEBV_", "", Trait_raw),
         Collection = "Core (n = 423)")

# pivot global to long
gebv_q_global_long <- gebv_q_global %>%
  select(-any_of("Group")) %>%
  pivot_longer(-Line, names_to = "Trait_raw", values_to = "GEBV") %>%
  mutate(Trait = gsub("^GEBV_", "", Trait_raw),
         Collection = "Global (n = 10,026)")

# combine
gebv_q_long <- bind_rows(gebv_q_core_long, gebv_q_global_long)

# scale within each trait AND collection separately
#gebv_q_scaled <- gebv_q_long %>%
  #group_by(Trait, Collection) %>%
  #mutate(GEBV_scaled = scale(GEBV)[,1]) %>%
  #ungroup()


# second - scale together
gebv_q_scaled <- gebv_q_long %>%
  group_by(Trait) %>%                # remove Collection from group_by
  mutate(GEBV_scaled = scale(GEBV)[,1]) %>%
  ungroup()


# join PA for both collections — keep only 23 quality traits
# filter to just Title_Case traits (not DAS/DAT)
quality_traits_423_clean <- quality_traits_423 %>%
  filter(grepl("^[A-Z][a-z]", Trait)) %>%
  mutate(Collection = "Core (n = 423)")

quality_traits_10k_clean <- quality_traits_10k %>%
  filter(grepl("^[A-Z][a-z]", Trait)) %>%
  mutate(Collection = "Global (n = 10,026)")

pa_quality <- bind_rows(quality_traits_423_clean, quality_traits_10k_clean)

# get PA label per trait per collection for facet labels
pa_label <- pa_quality %>%
  mutate(PA_label = round(PA, 2))

# join PA onto scaled data for labeling
gebv_q_plot <- gebv_q_scaled %>%
  left_join(pa_label %>% select(Trait, Collection, PA_label), by = c("Trait", "Collection"))

# create facet label with both PAs shown
pa_wide <- pa_label %>%
  select(Trait, Collection, PA_label) %>%
  pivot_wider(names_from = Collection, values_from = PA_label) %>%
  rename(PA_core = `Core (n = 423)`, PA_global = `Global (n = 10,026)`) %>%
  mutate(facet_label = paste0(gsub("_", " ", Trait),
                              "\n(PA core = ", PA_core,
                              ", global = ", PA_global, ")"))

# join facet labels
gebv_q_plot <- gebv_q_plot %>%
  left_join(pa_wide %>% select(Trait, facet_label), by = "Trait") %>%
  mutate(
    facet_label = factor(facet_label,
                         levels = pa_wide %>%
                           arrange(desc(PA_core)) %>%
                           pull(facet_label)),
    Collection = factor(Collection,
                        levels = c("Global (n = 10,026)", "Core (n = 423)"))
  )

# plot
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
  # add more here as needed
)

p_quality <- gebv_q_plot %>%
  filter(!Trait %in% traits_to_remove) %>%   # add this line
  ggplot(aes(x = GEBV_scaled, fill = Collection, color = Collection)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values  = c("Global (n = 10,026)" = "grey60",
                                "Core (n = 423)"      = "#B07AA1")) +
  scale_color_manual(values = c("Global (n = 10,026)" = "grey40",
                                "Core (n = 423)"      = "#7B4F7B")) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(face = "bold", size = 7),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 9)
  ) +
  labs(x        = "Scaled GEBV (z-score)",
       y        = "Density",
       title    = "Quality trait GEBV distributions — core vs global collection",
       subtitle = "Scaled within trait and collection")

p_quality

#longer X
p_quality <- gebv_q_plot %>%
  filter(!Trait %in% traits_to_remove) %>%
  ggplot(aes(x = GEBV_scaled, fill = Collection, color = Collection)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values  = c("Global (n = 10,026)" = "grey60",
                                "Core (n = 423)"      = "#B07AA1")) +
  scale_color_manual(values = c("Global (n = 10,026)" = "grey40",
                                "Core (n = 423)"      = "#7B4F7B")) +
  coord_cartesian(xlim = c(-5, 5)) +    # add this line
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(face = "bold", size = 7),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 9)
  ) +
  labs(x        = "Scaled GEBV (z-score)",
       y        = "Density",
       title    = "",
       subtitle = "")

p_quality


ggsave(
  filename = "~/R/WorldVeg_Capsicum/Figures/Figure3_quality_density.pdf",
  plot = p_quality,
  width = 10,
  height = 6,
  units = "in"
)

#########################
# all 23
p_quality <- gebv_q_plot %>%
  ggplot(aes(x = GEBV_scaled, fill = Collection, color = Collection)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values  = c("Global (n = 10,026)" = "grey60",
                                "Core (n = 423)"      = "#B07AA1")) +
  scale_color_manual(values = c("Global (n = 10,026)" = "grey40",
                                "Core (n = 423)"      = "#7B4F7B")) +
  coord_cartesian(xlim = c(-5, 5)) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(face = "bold", size = 7),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 9)
  ) +
  labs(x        = "Scaled GEBV (z-score)",
       y        = "Density",
       title    = "",
       subtitle = "")

p_quality
