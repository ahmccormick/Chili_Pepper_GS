########################################
#Individual Marker Effect Sizes - YIELD - ctrl timepoint
########################################
library(rrBLUP)
library(dplyr)
library(tidyr)

# 1. Use already-loaded genotype matrix
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)


# 2. Load phenotype (yield)
pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)

# 3. Filter both to common lines
common_ids <- intersect(rownames(geno_df), rownames(pheno))
geno_df <- geno_df[common_ids, ]
pheno <- pheno[common_ids, ]

# 4. Load training genotypes
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# 5. Subset training data
geno_train <- geno_df[rownames(geno_df) %in% train_genotypes, ]
pheno_train <- pheno[rownames(pheno) %in% train_genotypes, ]

# 6. Run rrBLUP
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"
yield_values <- pheno_train$yield
solve_yield <- mixed.solve(y = yield_values, Z = geno_train)
u.hat_yield <- solve_yield$u  # SNP effects

# 7. Apply to all individuals
geno_all <- as.matrix(geno_df)
mode(geno_all) <- "numeric"
individual_marker_effects <- geno_all * matrix(u.hat_yield, nrow = nrow(geno_all), ncol = ncol(geno_all), byrow = TRUE)

# 8. Reshape + export
individual_marker_effects_df <- as.data.frame(individual_marker_effects)
individual_marker_effects_df$Individual <- rownames(geno_all)

individual_marker_effects_long <- individual_marker_effects_df %>%
  pivot_longer(cols = -Individual, names_to = "SNP", values_to = "Effect")

write.csv(individual_marker_effects_long,
          "~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP.csv",
          row.names = FALSE)


########################################
#Individual Marker Effect Sizes - YIELD - heat1 timepoint
########################################
library(rrBLUP)
library(dplyr)
library(tidyr)

# 1. Use already-loaded genotype matrix
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
#geno_df <- geno_df_imputed  # already in memory and imputed

# 2. Load phenotype (yield)
#pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)


# 3. Filter both to common lines
common_ids <- intersect(rownames(geno_df), rownames(pheno))
geno_df <- geno_df[common_ids, ]
pheno <- pheno[common_ids, ]

# 4. Load training genotypes
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# 5. Subset training data
geno_train <- geno_df[rownames(geno_df) %in% train_genotypes, ]
pheno_train <- pheno[rownames(pheno) %in% train_genotypes, ]

# 6. Run rrBLUP
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"
yield_values <- pheno_train$yield
solve_yield <- mixed.solve(y = yield_values, Z = geno_train)
u.hat_yield <- solve_yield$u  # SNP effects

# 7. Apply to all individuals
geno_all <- as.matrix(geno_df)
mode(geno_all) <- "numeric"
individual_marker_effects <- geno_all * matrix(u.hat_yield, nrow = nrow(geno_all), ncol = ncol(geno_all), byrow = TRUE)

# 8. Reshape + export
individual_marker_effects_df <- as.data.frame(individual_marker_effects)
individual_marker_effects_df$Individual <- rownames(geno_all)

individual_marker_effects_long <- individual_marker_effects_df %>%
  pivot_longer(cols = -Individual, names_to = "SNP", values_to = "Effect")

write.csv(individual_marker_effects_long,
          "~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP_heat1.csv",
          row.names = FALSE)


########################################
#Individual Marker Effect Sizes - YIELD - heat2 timepoint
########################################
library(rrBLUP)
library(dplyr)
library(tidyr)

# 1. Use already-loaded genotype matrix
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
#geno_df <- geno_df_imputed  # already in memory and imputed

# 2. Load phenotype (yield)
#pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)


# 3. Filter both to common lines
common_ids <- intersect(rownames(geno_df), rownames(pheno))
geno_df <- geno_df[common_ids, ]
pheno <- pheno[common_ids, ]

# 4. Load training genotypes
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# 5. Subset training data
geno_train <- geno_df[rownames(geno_df) %in% train_genotypes, ]
pheno_train <- pheno[rownames(pheno) %in% train_genotypes, ]

# 6. Run rrBLUP
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"
yield_values <- pheno_train$yield
solve_yield <- mixed.solve(y = yield_values, Z = geno_train)
u.hat_yield <- solve_yield$u  # SNP effects

# 7. Apply to all individuals
geno_all <- as.matrix(geno_df)
mode(geno_all) <- "numeric"
individual_marker_effects <- geno_all * matrix(u.hat_yield, nrow = nrow(geno_all), ncol = ncol(geno_all), byrow = TRUE)

# 8. Reshape + export
individual_marker_effects_df <- as.data.frame(individual_marker_effects)
individual_marker_effects_df$Individual <- rownames(geno_all)

individual_marker_effects_long <- individual_marker_effects_df %>%
  pivot_longer(cols = -Individual, names_to = "SNP", values_to = "Effect")

write.csv(individual_marker_effects_long,
          "~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP_heat2.csv",
          row.names = FALSE)


###############################################################################
#Hot 28
###############################################################################
# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

library(readr)

control <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP.csv",
                    col_types = cols(
                      Individual = col_character(),
                      Effect = col_double(),
                      SNP = col_character()  # <- critical fix
                    ))

heat1 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP_heat1.csv",
                  col_types = cols(
                    Individual = col_character(),
                    Effect = col_double(),
                    SNP = col_character()
                  ))

# Load and fix sample ID column
hot_28  <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/Extreme_V_SampleIDs.csv")
colnames(hot_28)[1] <- "Sample_ID"

# Subset to the 28 heat-loving lines
selected_ids <- hot_28$Sample_ID
control_28 <- control %>% filter(Individual %in% selected_ids)
heat1_28   <- heat1   %>% filter(Individual %in% selected_ids)

# Merge control and heat1 on Individual and SNP
merged <- control_28 %>%
  rename(Effect_control = Effect) %>%
  inner_join(heat1_28 %>% rename(Effect_heat1 = Effect), by = c("Individual", "SNP"))


# Step 4: Create wide format for paired t-tests
wide <- merged %>%
  pivot_wider(names_from = Individual, values_from = c(Effect_control, Effect_heat1))

# Step 5: Calculate t-tests per SNP (with Bonferroni correction)
test_results <- merged %>%
  group_by(SNP) %>%
  summarize(
    t_stat = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$statistic
    } else { NA },
    p_value = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$p.value
    } else { NA },
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "bonferroni"),
    significant = adjusted_p_value < 0.001
  )


# Step 6: Extract significant SNPs and plot
significant_snps <- merged %>%
  filter(SNP %in% test_results$SNP[test_results$significant]) %>%
  group_by(SNP) %>%
  summarize(
    Mean_Effect_Heat1 = mean(Effect_heat1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Chromosome = gsub(":.*", "", SNP),
    Position = as.numeric(gsub(".*:", "", SNP)),
    Position_Mb = Position / 1e6
  ) %>%
  arrange(as.numeric(Chromosome), Position) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = as.character(1:12)),
    Genomic_Position = row_number()
  )

significant_snps <- significant_snps %>%
  filter(abs(Mean_Effect_Heat1) > 0.001)

write_csv(significant_snps, "~/R/World_Veg_Collab_Pepper/outputs/significant_snps_hot28.csv")

# Step 7: Plot
plot_faceted <- ggplot(significant_snps, aes(x = Position / 1e6, y = Mean_Effect_Heat1)) +
  geom_point(color = "red4", size = 1.2, alpha = 0.7) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  labs(
    title = "Significant SNP Marker Effects by Chromosome",
    x = "Genomic Position (Mb)",
    y = "Mean Effect Size (Heat1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

#print(plot_faceted)

# Optional save
ggsave("~/R/World_Veg_Collab_Pepper/FIGS/Heat_Loving_SNPs_Heat1_vs_Control_315snp.pdf",
       plot = plot_faceted, width = 10, height = 6)

###############################################################################
#Cold 57
###############################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Step 1: Load full marker effect datasets
control <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP.csv",
                    col_types = cols(
                      Individual = col_character(),
                      Effect = col_double(),
                      SNP = col_character()
                    ))

heat1 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP_heat1.csv",
                  col_types = cols(
                    Individual = col_character(),
                    Effect = col_double(),
                    SNP = col_character()
                  ))

# Step 2: Load cold-tolerant line IDs
cold_57 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/Extreme_InvertedV_SampleIDs.csv")
colnames(cold_57)[1] <- "Sample_ID"
selected_ids <- cold_57$Sample_ID

# Step 3: Subset genotype effect data to only cold-tolerant lines
control_57 <- control %>% filter(Individual %in% selected_ids)
heat1_57   <- heat1   %>% filter(Individual %in% selected_ids)

# Step 4: Merge datasets by SNP and Individual
merged <- control_57 %>%
  rename(Effect_control = Effect) %>%
  inner_join(heat1_57 %>% rename(Effect_heat1 = Effect), by = c("Individual", "SNP"))

# Step 5: (Optional) Wide format if needed
wide <- merged %>%
  pivot_wider(names_from = Individual, values_from = c(Effect_control, Effect_heat1))

# Step 6: Paired t-tests for each SNP
test_results <- merged %>%
  group_by(SNP) %>%
  summarize(
    t_stat = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$statistic
    } else { NA },
    p_value = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$p.value
    } else { NA },
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "bonferroni"),
    significant = adjusted_p_value < 0.001
  )

# Step 7: Summarize significant SNPs
significant_snps <- merged %>%
  filter(SNP %in% test_results$SNP[test_results$significant]) %>%
  group_by(SNP) %>%
  summarize(
    Mean_Effect_Heat1 = mean(Effect_heat1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Chromosome = gsub(":.*", "", SNP),
    Position = as.numeric(gsub(".*:", "", SNP)),
    Position_Mb = Position / 1e6
  ) %>%
  arrange(as.numeric(Chromosome), Position) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = as.character(1:12)),
    Genomic_Position = row_number()
  )

# Optional: Filter by effect size threshold
significant_snps <- significant_snps %>%
  filter(abs(Mean_Effect_Heat1) > 0.001)

write_csv(significant_snps, "~/R/World_Veg_Collab_Pepper/outputs/significant_snps_cold57.csv")

# Step 8: Plot faceted by chromosome
plot_faceted <- ggplot(significant_snps, aes(x = Position_Mb, y = Mean_Effect_Heat1)) +
  geom_point(color = "blue4", size = 1.2, alpha = 0.7) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  labs(
    title = "Significant SNP Marker Effects by Chromosome (Cold-Loving Lines)",
    x = "Genomic Position (Mb)",
    y = "Mean Effect Size (Heat1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Step 9: Save plot
ggsave("~/R/World_Veg_Collab_Pepper/FIGS/Cold_Loving_SNPs_Heat1_vs_Control_57lines.pdf",
       plot = plot_faceted, width = 10, height = 6)



