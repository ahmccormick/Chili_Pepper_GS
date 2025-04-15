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

