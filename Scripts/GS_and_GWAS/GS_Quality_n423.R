##########################################
#McLeod Phenotype data - 23 traits (for n=329 lines)
##########################################
#Downloaded pheno data from 
#https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/S7JVEM
pheno_data <- read.csv("~/R/World_Veg_Collab_Pepper/McLeod_pheno_GS/G2P-SOL_Pepper_CC_Pheno_data_2023.csv", sep = ";")

head(pheno_data)
library(dplyr)

per_line_pheno <- pheno_data %>%
  group_by(G2PSOL_code) %>%  # Group by your line name
  summarize(across(where(is.numeric), mean, na.rm = TRUE))  # Average numeric columns

# View the result
head(per_line_pheno)

library(dplyr)

# Step 1: Keep only lines with no NA values
per_line_pheno_complete <- per_line_pheno %>%
  filter(if_all(where(is.numeric), ~ !is.na(.)))

per_line_pheno_complete <- per_line_pheno_complete %>%
  rename(
    `Axis_length` = AL,
    `Brix` = BX,
    `Fruit_external_color_a` = FCa,
    `Fruit_external_color_b` = FCb,
    `Fruit_external_color_L` = FCL,
    `Fruit_maximum_length` = FLe,
    `Flowering_time` = FlT,
    `Fruit_shape_index` = FShI,
    `Fruit_weight` = FWe,
    `Fruit_maximum_width` = FWi,
    `Immature_fruit_external_color_a` = IFCa,
    `Immature_fruit_external_color_b` = IFCb,
    `Immature_fruit_external_color_L` = IFCL,
    `Locule_number` = LN,
    `Total_plant_height` = PHe,
    `Pericarp_thickness` = PTh,
    `Total_fruit_number` = TFN,
    `Total_fruit_weight` = TFWe,
    `Fruit_fasciation` = FF,
    `Fruit_load` = FLo,
    `Fruit_pungency` = FP,
    `Fruit_predominant_shape_oblate` = FShO,
    `External_immature_fruit_color_green` = IFCG
  )

# Step 3: View your cleaned dataset
View(per_line_pheno_complete)

#write.csv(per_line_pheno_complete, "~/R/World_Veg_Collab_Pepper/McLeod_pheno_GS/McLeod_per_line_pheno_data_n329.csv", row.names = FALSE)

##########################################

#cut out columns 
#write.csv(per_line_pheno_complete, "~/R/World_Veg_Collab_Pepper/McLeod_pheno_GS/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = FALSE)

############################################################################################################
# Step 2 _GS
############################################################################################################
install.packages("rrBLUP")
library(rrBLUP)
library(dplyr)

##################
# GENOTYPE IMPORT (IMPUTED)
##################
geno_df_imputed <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv",
                            row.names = 1, check.names = FALSE)

##################
# PHENOTYPE IMPORT
##################
pheno_McLeod <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)

training_set_lines <-read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", header = TRUE)

# Get vector of training line names from genotype file
training_line_names <- training_set_lines[[1]]

# Subset phenotype data to match those training lines
pheno_train <- pheno_McLeod[rownames(pheno_McLeod) %in% training_line_names, ] #gives 132 overlap or 31%

##################
# TRAINING SET FILTERING
##################
# Train = lines that have phenotype data
train_genotypes <- rownames(pheno_train)

# Test = all remaining lines without phenotype data
test_genotypes <- setdiff(rownames(geno_df_imputed), train_genotypes)


geno_train <- geno_df_imputed[train_genotypes, ]
geno_test  <- geno_df_imputed[test_genotypes, ]

# Make sure matrices are numeric
geno_train <- as.matrix(sapply(geno_train, as.numeric))
rownames(geno_train) <- train_genotypes

geno_test  <- as.matrix(sapply(geno_test, as.numeric))
rownames(geno_test)  <- test_genotypes



run_rrBLUP_all_traits <- function(geno_train, geno_test, pheno_train, output_prefix) {
  pheno_train <- pheno_train[, colnames(pheno_train) != "X"]  # Remove 'X' column if present
  
  traits <- colnames(pheno_train)
  GEBVs_train_list <- list()
  GEBVs_test_list <- list()
  
  for (trait in traits) {
    cat("Training rrBLUP model for", trait, "...\n")
    
    y <- pheno_train[[trait]]
    y[is.na(y)] <- mean(y, na.rm = TRUE)
    
    if (var(y, na.rm = TRUE) == 0) {
      cat("Skipping", trait, "due to zero variance.\n")
      next
    }
    
    rrblup_model <- mixed.solve(y = y, Z = geno_train)
    
    GEBVs_train <- geno_train %*% rrblup_model$u
    GEBVs_test  <- geno_test %*% rrblup_model$u
    
    GEBVs_train_list[[trait]] <- as.data.frame(GEBVs_train)
    GEBVs_test_list[[trait]]  <- as.data.frame(GEBVs_test)
    
    colnames(GEBVs_train_list[[trait]]) <- paste0("GEBV_", trait, "_Train")
    colnames(GEBVs_test_list[[trait]])  <- paste0("GEBV_", trait, "_Test")
  }
  
  GEBVs_train_final <- do.call(cbind, GEBVs_train_list)
  GEBVs_test_final  <- do.call(cbind, GEBVs_test_list)
  
  rownames(GEBVs_train_final) <- rownames(geno_train)
  rownames(GEBVs_test_final)  <- rownames(geno_test)
  
  write.csv(GEBVs_train_final, paste0(output_prefix, "_Train.csv"), row.names = TRUE)
  write.csv(GEBVs_test_final,  paste0(output_prefix, "_Test.csv"),  row.names = TRUE)
  
  cat("GEBV predictions saved for all traits.\n")
}



run_rrBLUP_all_traits(geno_train, geno_test, pheno_train, "/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_23traits_n423")
################

# Load both files
gebv_train <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_23traits_n423_Train.csv", row.names = 1)
gebv_test  <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_23traits_n423_Test.csv", row.names = 1)

# Remove _Train and _Test from column names
colnames(gebv_train) <- gsub("_Train$", "", colnames(gebv_train))
colnames(gebv_test)  <- gsub("_Test$", "", colnames(gebv_test))

# Add a column indicating Set
gebv_train$Set <- "Train"
gebv_test$Set  <- "Test"

# Make sure Set is first column (optional)
gebv_train <- gebv_train[, c("Set", setdiff(names(gebv_train), "Set"))]
gebv_test  <- gebv_test[, c("Set", setdiff(names(gebv_test), "Set"))]

# Combine
gebv_all <- rbind(gebv_train, gebv_test)

# Save
write.csv(gebv_all, "/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv", row.names = TRUE)


################
#PLOTTING 23X GEBV HEATMAP
################
# Load your combined 23-trait GEBVs
gebv_23traits <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv", row.names = 1)


library(pheatmap)

# Drop non-numeric columns if needed
gebv_mat <- gebv_23traits %>% select(where(is.numeric))

# Optionally scale the data (row-wise)
pheatmap(gebv_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 5,
         main = "Heatmap of GEBVs across Traits and Lines")



####SCALE BY GEBV RANGES TO COMPARE
# Standardize the matrix (z-score by column)
scaled_mat <- scale(gebv_mat)

pheatmap(scaled_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Z-scored GEBVs across Traits and Lines")


#########################################################################################################
#Step 3
#########################################################################################################
#Model Prediction Accuracies - ran on HPC
#########################################################################################################

####################
#Method one: rrBLUP
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
#install.packages("hibayes")
library(ggplot2)


####################
#STEP ONE
####################
# Load genotype matrix (imputed, same as before)
geno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load phenotype matrix for 23 traits
pheno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)

# Load training set line names
train_genotypes <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", header = TRUE)[[1]]

# Filter genotype matrix to match phenotype
common_genotypes <- intersect(rownames(geno_data), rownames(pheno_data))
geno_filtered <- geno_data[common_genotypes, ]
pheno_filtered <- pheno_data[common_genotypes, ]

geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

# Filter training genotypes to only those present in the data
train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes_filtered, ]
geno_test  <- geno_matrix[setdiff(rownames(geno_matrix), train_genotypes_filtered), ]

# Subset phenotype matrix to match training set
pheno_train <- pheno_filtered[train_genotypes_filtered, ]
pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

####################
#STEP TWO
####################

pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

####################
#STEP THREE
####################
# Load custom cross-validation function
#source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_kfold_functions.R")
source("/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

# Drop Genotype column before cross-validation
pheno_train_numeric <- pheno_train %>% select(-Genotype)

# Run cross-validation
xval_k10_rrblup <- k.xval(
  g.in = geno_train,
  y.in = pheno_train_numeric,
  y.trainset = pheno_train_numeric,
  k.fold = 10,
  reps = 50
)
# Save output
saveRDS(xval_k10_rrblup, "/Users/annamccormick/R/World_Veg_Project/data/outputs/xval_rrblup_23traits_kfold_10.RDS")


####################
#Method Two: Gaussian Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

# Load genotype and phenotype data
geno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
pheno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)
train_genotypes <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", header = TRUE)[[1]]

# Filter and format
common_genotypes <- intersect(rownames(geno_data), rownames(pheno_data))
geno_matrix <- as.matrix(geno_data[common_genotypes, ])
mode(geno_matrix) <- "numeric"
pheno_filtered <- pheno_data[common_genotypes, ]

train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))
geno_train <- geno_matrix[train_genotypes_filtered, ]
pheno_train <- pheno_filtered[train_genotypes_filtered, ]

# Add 'Genotype' column for kinship models
pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

# Also drop 'Genotype' for matrix-only methods
pheno_train_numeric <- pheno_train %>% select(-Genotype)

####################
source("/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

# Compute Gaussian distance matrix (Euclidean distance on A.mat)
K <- A.mat(geno_train)
k_dist <- dist(K)

# Run k-fold CV
xval_k10_GAUSS <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = pheno_train,
  y.trainset = pheno_train,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)

saveRDS(xval_k10_GAUSS, "/Users/annamccormick/R/World_Veg_Project/data/outputs/xval_GAUSS_kfold_10.RDS")


####################
# Method Three: Exponential Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
# STEP ONE: Load & Filter Data
####################

# Load genotype matrix
geno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv", 
                      row.names = 1, check.names = FALSE)

# Load phenotype matrix (control condition)
control_pheno_check <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", 
                                row.names = 1)

# Load training genotypes
train_genotypes <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", 
                            header = TRUE)[[1]]

# Filter for genotypes in both genotype and phenotype datasets
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))
geno_matrix <- as.matrix(geno_data[common_genotypes, ])
mode(geno_matrix) <- "numeric"
control_pheno_filtered <- control_pheno_check[common_genotypes, ]

# Filter training set
train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))
geno_train <- geno_matrix[train_genotypes_filtered, ]
control_train <- control_pheno_filtered[train_genotypes_filtered, ]

# Add 'Genotype' column for kin.blup compatibility
control_train$Genotype <- rownames(control_train)
control_train <- control_train %>% relocate(Genotype)

####################
# STEP TWO: Compute Exponential Kernel
####################
source("/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

# Compute exponential kernel matrix
K.Exp <- Kernel_computation(X = geno_train, name = "exponential", degree = NULL, nL = NULL)
exp_dist <- dist(K.Exp)
exp_dist_mat <- as.matrix(exp_dist)

# Set row/col names explicitly
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)

# Sanity checks
stopifnot(all(rownames(exp_dist_mat) == rownames(geno_train)))
stopifnot(all(rownames(exp_dist_mat) == control_train$Genotype))

####################
# STEP THREE: Run Cross-Validation
####################
xval_k10_EXP <- k.xval.EXP(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)

# Save output
saveRDS(xval_k10_EXP, "/Users/annamccormick/R/World_Veg_Project/data/outputs/xval_EXP_kfold_10.RDS")




###################################################################################################
#Plot prediction accuracies
###################################################################################################

library(ggplot2)
library(dplyr)
library(readr)

# Load cross-validation results (these should contain the $xval.result list)
rrblup_results <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_rrblup_23traits_kfold_10.RDS")$xval.result
gauss_results  <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_GAUSS_kfold_10.RDS")$xval.result
exp_results    <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_EXP_kfold_10.RDS")$xval.result

# Add model labels
rrblup_results$model <- "rrBLUP"
gauss_results$model  <- "Gaussian"
exp_results$model    <- "Exponential"

# Combine into one dataframe
all_models <- bind_rows(rrblup_results, gauss_results, exp_results)

# Ensure correct data types
all_models <- all_models %>%
  mutate(
    r.mean = as.numeric(r.mean),
    r.sd = as.numeric(r.sd),
    trait = as.factor(trait),
    model = factor(model, levels = c("rrBLUP", "Gaussian", "Exponential"))
  )

# Plot
ggplot(all_models, aes(x = model, y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
                width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(~trait, scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color = "red", linewidth = 1.2, linetype = "longdash") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 10)
  ) +
  labs(
    title = "Prediction Accuracy (r) for 23 Traits",
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  ylim(0, 1)

#export results as table
library(tidyr)

pa_table <- all_models %>%
  select(trait, model, r.mean, r.sd) %>%
  arrange(trait, model) %>%
  pivot_wider(
    names_from = model,
    values_from = c(r.mean, r.sd)
  )

readr::write_csv(
  pa_table,
  "~/R/World_Veg_Project/data/outputs/PA_23trait_core423/PA_summary_23traits.csv"
)

###################################################################################################

