
######################Outline for document#########################
# Step 1: 
  #1A - Pre-processing - conversion of .vcf file to genotype matrix
  #1B - Convert genotype matrix to rrBLUP format
  #1C - Transpose genotype matrix for downstream use
  #1D - Genotype import & cleaning/imputation 
####################################################################
# Step 2:
  # Genotype import
  # Phenotype import
  # Set training and test populations
  # Run GS
####################################################################
# Step 3:
# Prediction accuracy using
 # rrBLUP
 # Gaussian Kernel
 # Exponential Kernel
####################################################################



####################################################################
#1A FROM VCF TO GENOTYPE MATRIX - rows=individuals  columns =SNPs, in 0/1/2 format
####################################################################
library(vcfR)
library(adegenet)
vcf <- read.vcfR("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_filtered_and_LD_pruned.vcf.gz")
genotype_matrix <- extract.gt(vcf, element = "GT")

dim(genotype_matrix)  # Should return (SNPs × Samples)
write.csv(genotype_matrix, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_raw.csv", row.names = TRUE)

genotype_matrix_raw <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_raw.csv", 
                                row.names = 1, check.names = FALSE)

# Check dimensions
dim(genotype_matrix_raw)  # Should match (SNPs × Samples)

# Preview 
head(genotype_matrix_raw[, 1:5])  #first 5 samples

####################################################################
#1B - Convert genotype matrix to rrBLUP format
####################################################################
# Function to convert genotypes to numeric format for rrBLUP
convert_geno <- function(x) {
  ifelse(x == "0/0", 0, 
         ifelse(x == "0/1" | x == "1/0", 1, 
                ifelse(x == "1/1", 2, NA)))
}

# Apply conversion to the entire matrix
geno_numeric <- apply(genotype_matrix_raw, 2, convert_geno)

# Convert to dataframe
geno_df <- as.data.frame(geno_numeric)

# Check dimensions after conversion
dim(geno_df)  #(340734 × 423)

# Save the converted genotype matrix
write.csv(geno_df, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format.csv", row.names = TRUE)

####################################################################
#1C - Transpose genotype matrix for downstream use
####################################################################
library(rrBLUP)

# genotype matrix in rrBLUP format
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format.csv", 
                    row.names = 1, check.names = FALSE)

# Check dimensions to confirm SNPs x Samples
dim(geno_df)

# Preview first few values
head(geno_df[, 1:5])

geno_df <- t(geno_df)  # Transpose so genotypes (samples) become rows
geno_df <- as.data.frame(geno_df)  # Convert back to a dataframe
head(geno_df[, 1:5])
cat("Corrected genotype matrix dimensions:", dim(geno_df), "\n")


write.csv(geno_df, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_transposed.csv", row.names = TRUE)


##################
#1D - GENOTYPE IMPORT & CLEANING/IMPUTATION
##################

geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_transposed.csv", 
                    row.names = 1, check.names = FALSE)

# Define imputation function (replace missing values with the most common allele)
impute_geno <- function(x) {
  x[is.na(x)] <- as.numeric(names(which.max(table(x, useNA = "no"))))  # Replace NA with mode
  return(x)
}

# Apply the imputation function to each SNP column
cat("Imputing missing values in genotype data...\n")
geno_df_imputed <- apply(geno_df, 2, impute_geno)

# Convert back to dataframe
geno_df_imputed <- as.data.frame(geno_df_imputed)

# Confirm missing values after imputation
num_missing_after <- sum(is.na(geno_df_imputed))
cat("Missing values after imputation:", num_missing_after, "\n") # Should be 0

# Save the cleaned genotype matrix
write.csv(geno_df_imputed, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = TRUE)


############################################################################################################
# Step 2
############################################################################################################
library(rrBLUP)
library(dplyr)

##################
# GENOTYPE IMPORT (IMPUTED)
##################
geno_df_imputed <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv",
                            row.names = 1, check.names = FALSE)

##################
# PHENOTYPE IMPORT
##################

control_pheno_check <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)

##################
#TRAINING
##################
set.seed(123)  # For reproducibility

# Identify lines with complete phenotype data across ALL conditions
complete_lines <- rownames(control_pheno_check)[
  rowSums(is.na(control_pheno_check)) == 0 &
    rowSums(is.na(heat_stress_1_check)) == 0 &
    rowSums(is.na(heat_stress_2_check)) == 0
]

# Randomly sample 150 genotypes only from the fully complete lines
train_genotypes <- sample(complete_lines, 150)

##################
#TEST
##################
# Define test set as remaining lines (all genotyped lines except those in training)
test_genotypes <- setdiff(rownames(geno_df_imputed), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_df_imputed[train_genotypes, ]
geno_test <- geno_df_imputed[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

# Check missing values again
cat("Missing values in control_train:", sum(is.na(control_train)), "\n")
cat("Missing values in heat_stress_1_train:", sum(is.na(heat_stress_1_train)), "\n")
cat("Missing values in heat_stress_2_train:", sum(is.na(heat_stress_2_train)), "\n")



#convert genotype matrix to numeric for rrBLUP compatibility
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"

geno_test <- as.matrix(geno_test)
mode(geno_test) <- "numeric"

############################################################################################################
library(rrBLUP)

# Function to train rrBLUP and predict GEBVs for multiple traits
run_rrBLUP_all_traits <- function(geno_train, geno_test, pheno_train, output_prefix) {
  
  traits <- colnames(pheno_train)  # Get trait names
  GEBVs_train_list <- list()
  GEBVs_test_list <- list()
  
  for (trait in traits) {
    cat("Training rrBLUP model for", trait, "...\n")
    
    y <- pheno_train[[trait]]  # Extract trait values
    y[is.na(y)] <- mean(y, na.rm = TRUE)  # Impute missing values (just in case)
    
    # Skip if variance is zero (to prevent errors)
    if (var(y, na.rm = TRUE) == 0) {
      cat("Skipping", trait, "due to zero variance.\n")
      next
    }
    
    # Train rrBLUP model
    rrblup_model <- mixed.solve(y = y, Z = geno_train)
    
    # Predict GEBVs for training and test sets
    GEBVs_train <- geno_train %*% rrblup_model$u
    GEBVs_test <- geno_test %*% rrblup_model$u
    
    # Convert to data frames and store
    GEBVs_train_list[[trait]] <- as.data.frame(GEBVs_train)
    GEBVs_test_list[[trait]] <- as.data.frame(GEBVs_test)
    
    # Rename columns
    colnames(GEBVs_train_list[[trait]]) <- paste0("GEBV_", trait, "_Train")
    colnames(GEBVs_test_list[[trait]]) <- paste0("GEBV_", trait, "_Test")
  }
  
  # Combine GEBVs for all traits
  GEBVs_train_final <- do.call(cbind, GEBVs_train_list)
  GEBVs_test_final <- do.call(cbind, GEBVs_test_list)
  
  # Add row names
  rownames(GEBVs_train_final) <- rownames(geno_train)
  rownames(GEBVs_test_final) <- rownames(geno_test)
  
  # Save results
  write.csv(GEBVs_train_final, paste0(output_prefix, "_Train.csv"), row.names = TRUE)
  write.csv(GEBVs_test_final, paste0(output_prefix, "_Test.csv"), row.names = TRUE)
  
  cat("GEBV predictions saved for all traits.\n")
}

# Run rrBLUP for all traits in all conditions
run_rrBLUP_all_traits(geno_train, geno_test, control_train, 
                      "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control")
run_rrBLUP_all_traits(geno_train, geno_test, heat_stress_1_train, 
                      "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1")
run_rrBLUP_all_traits(geno_train, geno_test, heat_stress_2_train, 
                      "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2")



# Save the 150 training genotypes to a CSV file
write.csv(data.frame(Training_Genotypes = train_genotypes), 
          "~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv", 
          row.names = FALSE)

#########################################################################################################
#Step 3
#########################################################################################################
#Model Prediction Accuracies
#########################################################################################################

####################
#Method one: rrBLUP
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
#STEP ONE
####################

# Load Genotype Data (in rrBLUP format)
#geno_data <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Environmental and Phenotypic Data
control_pheno_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/heat_stress_2_filtered_286.csv", row.names = 1)

# Identify common genotypes
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))

# Filter the genotype matrix to match phenotype data
geno_filtered <- geno_data[common_genotypes, ]

# Ensure genotype data is in numeric format
geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
#STEP TWO
####################
# Load the previously saved training genotypes
#train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
train_genotypes <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# Define the test set as the remaining genotypes
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes, ]
geno_test <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]


####################
#STEP THREE
####################

# Load custom k-fold cross-validation function
source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
xval_k10_rrblup <- k.xval(g.in = geno_train, y.in = control_train, y.trainset = control_train, k.fold = 10, reps = 50)
#saveRDS(xval_k10_rrblup, "~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")
saveRDS(xval_k10_rrblup, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")

####################
#Method Two: Gaussian Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
#STEP ONE
####################

# Load Genotype Data (in rrBLUP format)
#geno_data <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Environmental and Phenotypic Data
control_pheno_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/heat_stress_2_filtered_286.csv", row.names = 1)

# Identify common genotypes
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))

# Filter the genotype matrix to match phenotype data
geno_filtered <- geno_data[common_genotypes, ]

# Ensure genotype data is in numeric format
geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
#STEP TWO
####################
# Load the previously saved training genotypes
#train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
train_genotypes <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# Define the test set as the remaining genotypes
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes, ]
geno_test <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

##############################################################################################################################################################
#add genotype column names to control_train as kin.blup needs it
control_train$Genotype <- rownames(control_train)

control_train <- control_train %>%
  relocate(Genotype)

print(colnames(control_train)[1])  # This should be "Genotype"

####################
#STEP THREE
####################

# Load custom k-fold cross-validation function
source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
#xval_k10_rrblup <- k.xval(g.in = geno_train, y.in = control_train, y.trainset = control_train, k.fold = 10, reps = 50)
#saveRDS(xval_k10_rrblup, "~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")
#saveRDS(xval_k10_rrblup, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")

# Compute relationship matrix using Gaussian Kernel
K <- A.mat(geno_train)
k_dist <- dist(K)

# Run k-fold cross-validation for Gaussian Kernel
xval_k10_GAUSS <- k.xval.GAUSS(g.in = geno_train, y.in = control_train, y.trainset = control_train, k_dist = k_dist, k.fold = 10, reps = 50)
#saveRDS(xval_k10_GAUSS, "~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_kfold_10.RData")
saveRDS(xval_k10_GAUSS, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/xval_GAUSS_kfold_10.RData")


####################
#Method Three: Exponential Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
#STEP ONE
####################

# Load Genotype Data (in rrBLUP format)
#geno_data <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Environmental and Phenotypic Data
control_pheno_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/heat_stress_2_filtered_286.csv", row.names = 1)

# Identify common genotypes
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))

# Filter the genotype matrix to match phenotype data
geno_filtered <- geno_data[common_genotypes, ]

# Ensure genotype data is in numeric format
geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
#STEP TWO
####################
# Load the previously saved training genotypes
#train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
train_genotypes <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# Define the test set as the remaining genotypes
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes, ]
geno_test <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

##############################################################################################################################################################
#add genotype column names to control_train as kin.blup needs it
control_train$Genotype <- rownames(control_train)

control_train <- control_train %>%
  relocate(Genotype)

print(colnames(control_train)[1])  # This should be "Genotype"

####################
#STEP THREE
####################

# Load custom k-fold cross-validation function
source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
#xval_k10_rrblup <- k.xval(g.in = geno_train, y.in = control_train, y.trainset = control_train, k.fold = 10, reps = 50)
#saveRDS(xval_k10_rrblup, "~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")
#saveRDS(xval_k10_rrblup, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")

# Compute relationship matrix using Gaussian Kernel
K <- A.mat(geno_train)
k_dist <- dist(K)

# Run k-fold cross-validation for Gaussian Kernel
#xval_k10_GAUSS <- k.xval.GAUSS(g.in = geno_train, y.in = control_train, y.trainset = control_train, k_dist = k_dist, k.fold = 10, reps = 50)
#saveRDS(xval_k10_GAUSS, "~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_kfold_10.RData")
#saveRDS(xval_k10_GAUSS, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_GAUSS_kfold_10.RData")

# Compute relationship matrix using Exponential Kernel
K.Exp <- Kernel_computation(X = geno_train, name = "exponential", degree = NULL, nL = NULL)
exp_dist <- dist(K.Exp)

##############################################################################################################################################################

rownames(geno_train) <- as.character(trimws(rownames(geno_train)))
control_train$Genotype <- as.character(trimws(control_train$Genotype))

print(setdiff(rownames(geno_train), control_train$Genotype))  # Should return character(0)
print(setdiff(control_train$Genotype, rownames(geno_train)))  # Should return character(0)

exp_dist_mat <- as.matrix(exp_dist)  # Convert from dist object
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)  # Assign correct row/col names

print(setdiff(rownames(exp_dist_mat), rownames(geno_train)))  # Should return character(0)
print(setdiff(rownames(exp_dist_mat), control_train$Genotype))  # Should return character(0)

##############################################################################################################################################################
# Run k-fold cross-validation for Exponential Kernel
xval_k10_EXP <- k.xval.EXP(g.in = geno_train, y.in = control_train, y.trainset = control_train, k_dist = exp_dist_mat, k.fold = 10, reps = 50)
#saveRDS(xval_k10_EXP, "~/R/World_Veg_Collab_Pepper/outputs/xval_EXP_kfold_10.RData")
saveRDS(xval_k10_EXP, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/xval_EXP_kfold_10.RData")

#########################################################################################################

####################
#Plotting Prediction Accuracies
####################
# Load cross-validation results
rrblup_kfold10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")$xval.result
gauss_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_kfold_10.RData")$xval.result
EXP_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_EXP_kfold_10.RData")$xval.result

#rrblup_kfold10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")$xval.result
#gauss_kfold_10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_GAUSS_kfold_10.RData")$xval.result
#EXP_kfold_10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_EXP_kfold_10.RData")$xval.result
#bayescpi_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_BayesCpi_kfold_10.RData")$xval.result

# Assign model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
#bayescpi_kfold_10$model <- "BayesCpi"

# Combine all results
#all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10, bayescpi_kfold_10)
all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

# Convert standard deviation to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

# Select traits to analyze
all_gs_traits <- all_models %>% filter(trait %in% unique(all_models$trait))  # Automatically selects all available traits

####################
#STEP FIVE
####################
#Plot prediction accuracy
ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +  # Adjust the number of rows as needed
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +  # Baseline for interpretation
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits in Chili Dataset",
       x = "Model",
       y = "Prediction Accuracy (r)")


############################
#1st model rrBLUP PAs
############################
# Convert columns to numeric (fixing character/factor issues)
all_gs_traits$r.mean <- as.numeric(all_gs_traits$r.mean)
all_gs_traits$r.sd <- as.numeric(all_gs_traits$r.sd)

# Check if there are any NA values after conversion
summary(all_gs_traits$r.mean)
summary(all_gs_traits$r.sd)

# Remove rows where r.mean or r.sd is NA (caused by conversion errors)
all_gs_traits <- all_gs_traits %>% filter(!is.na(r.mean), !is.na(r.sd))

ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits in Chili Dataset",
       x = "Model",
       y = "Prediction Accuracy (r)")

# Filter traits with mean prediction accuracy above 0.5
all_gs_traits_filtered <- all_gs_traits %>% filter(r.mean > 0.5)

# Generate the plot
ggplot(all_gs_traits_filtered, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for Traits Above 0.5",
       x = "Model",
       y = "Prediction Accuracy (r)")

#########################################################################################################



