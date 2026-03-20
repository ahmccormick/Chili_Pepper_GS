#############################################
# STI Genomic Selection - rrBLUP, Gaussian, Exponential Kernel
# Phenotype: pre-computed STI from BLUE-adjusted data
# GEBV predictions on all 423 lines
#############################################

rm(list = ls())

library(rrBLUP)
library(dplyr)
library(readr)
library(ggplot2)

#############################################
# 0) File paths -- EDIT THESE
#############################################

geno_file       <- "~/R/World_Veg_Project/filtered_data/genotype_matrix_imputed.csv"
sti_file        <- "~/R/World_Veg_Project/anna_BLUES/anna_BLUES_rebuild/STI_from_BLUEs_all_73_traits.csv"
training_file   <- "~/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv"
xval_functions  <- "~/R/World_Veg_Project/data/inputs/xval_kfold_functions.R"

# OUTPUT FOLDER
out_base <- "~/R/WorldVeg_Capsicum/STI_GS/"

#############################################
# 1) Load and clean STI phenotype data
#############################################

pheno_data <- read.csv(sti_file,
                       na.strings = c("NA", "", "NaN", " "),
                       check.names = FALSE) %>%
  rename(Genotype = accession) %>%
  mutate(Genotype = trimws(Genotype)) %>%
  filter(!is.na(Genotype), Genotype != "") %>%
  filter(if_any(where(is.numeric), ~ !is.na(.)))

message("Phenotype rows after cleanup: ", nrow(pheno_data))
message("Unique genotypes: ", n_distinct(pheno_data$Genotype))
message("Traits: ", ncol(pheno_data) - 1)

# Save cleaned phenotype file
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
write_csv(pheno_data,
          file.path(out_base, "cleaned_pheno_STI_73traits_BLUEadjusted.csv"))

#############################################
# 2) Load genotype matrix
#############################################

geno_data <- read.csv(geno_file, row.names = 1, check.names = FALSE)
geno_matrix_full <- as.matrix(geno_data)
mode(geno_matrix_full) <- "numeric"
message("Genotype matrix: ", nrow(geno_matrix_full), " samples x ", ncol(geno_matrix_full), " markers")

#############################################
# 3) Load training set and align data
#############################################

norm_id <- function(x) toupper(trimws(as.character(x)))

training_set   <- read.csv(training_file, header = TRUE)
training_names <- norm_id(training_set[[1]])

# normalize IDs
pheno_data$Genotype      <- norm_id(pheno_data$Genotype)
rownames(geno_matrix_full) <- norm_id(rownames(geno_matrix_full))

# common genotypes across geno and pheno (for training set alignment)
common_genotypes <- intersect(rownames(geno_matrix_full), pheno_data$Genotype)
pheno_filtered   <- pheno_data %>% filter(Genotype %in% common_genotypes)
rownames(pheno_filtered) <- pheno_filtered$Genotype

message("Common genotypes (geno + pheno): ", length(common_genotypes))

# training set -- must be in both geno and pheno
train_genotypes <- intersect(training_names, common_genotypes)

# test set -- ALL remaining lines in full genotype matrix (includes lines without phenotype)
test_genotypes  <- setdiff(rownames(geno_matrix_full), train_genotypes)

geno_train <- geno_matrix_full[train_genotypes, ]
geno_test  <- geno_matrix_full[test_genotypes,  ]

pheno_train         <- pheno_filtered[train_genotypes, ]
pheno_train_numeric <- pheno_train %>% select(-Genotype)

message("Training N: ",          length(train_genotypes))
message("Test N (all remaining): ", length(test_genotypes))
# Training = 150, Test = ~273 (all 423 minus 150)

#############################################
# 4) GEBV predictions -- rrBLUP
#############################################

trait_names <- names(pheno_train_numeric)
message("Running GEBV predictions for ", length(trait_names), " traits...")

GEBVs_train_list <- list()
GEBVs_test_list  <- list()

for (trait in trait_names) {
  
  y <- pheno_train_numeric[[trait]]
  if (all(is.na(y))) { message("Skipping ", trait, " - all NA"); next }
  y[is.na(y)] <- mean(y, na.rm = TRUE)
  if (var(y, na.rm = TRUE) == 0) { message("Skipping ", trait, " - zero variance"); next }
  
  fit <- rrBLUP::mixed.solve(y = y, Z = geno_train)
  
  GEBVs_train_list[[trait]] <- setNames(
    as.data.frame(geno_train %*% fit$u), paste0("GEBV_", trait, "_Train"))
  GEBVs_test_list[[trait]]  <- setNames(
    as.data.frame(geno_test  %*% fit$u), paste0("GEBV_", trait, "_Test"))
}

GEBVs_train_final <- do.call(cbind, GEBVs_train_list)
GEBVs_test_final  <- do.call(cbind, GEBVs_test_list)
rownames(GEBVs_train_final) <- train_genotypes
rownames(GEBVs_test_final)  <- test_genotypes

# OUTPUT: GEBV files
write.csv(GEBVs_train_final,
          file.path(out_base, "GEBVs_STI_73traits_BLUEadjusted_Train.csv"),
          row.names = TRUE)
write.csv(GEBVs_test_final,
          file.path(out_base, "GEBVs_STI_73traits_BLUEadjusted_Test.csv"),
          row.names = TRUE)
message("GEBV train and test files saved.")

# Combine train + test
gebv_train_clean <- GEBVs_train_final
gebv_test_clean  <- GEBVs_test_final
colnames(gebv_train_clean) <- gsub("_Train$", "", colnames(gebv_train_clean))
colnames(gebv_test_clean)  <- gsub("_Test$",  "", colnames(gebv_test_clean))
gebv_train_clean$Set <- "Train"
gebv_test_clean$Set  <- "Test"
gebv_all <- rbind(
  gebv_train_clean[, c("Set", setdiff(names(gebv_train_clean), "Set"))],
  gebv_test_clean[ , c("Set", setdiff(names(gebv_test_clean),  "Set"))]
)

# OUTPUT: combined GEBV file (all 423)
write.csv(gebv_all,
          file.path(out_base, "GEBVs_STI_73traits_BLUEadjusted_ALL_n423.csv"),
          row.names = TRUE)
message("Combined GEBV file (all 423) saved.")




####################################################################################################################################################################
# PAs (on HPC)
####################################################################################################################################################################
rm(list = ls())

library(rrBLUP)
library(dplyr)
library(readr)
library(ggplot2)

#############################################
# 0) File paths
#############################################

geno_file       <- "~/R/World_Veg_Project/filtered_data/genotype_matrix_imputed.csv"
sti_file        <- "~/R/World_Veg_Project/anna_BLUES/anna_BLUES_rebuild/STI_from_BLUEs_all_73_traits.csv"
training_file   <- "~/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv"
xval_functions  <- "~/R/World_Veg_Project/data/inputs/xval_kfold_functions.R"
out_base        <- "~/R/WorldVeg_Capsicum/STI_GS/"

#############################################
# 1) Load phenotype
#############################################

pheno_data <- read.csv(sti_file,
                       na.strings = c("NA", "", "NaN", " "),
                       check.names = FALSE) %>%
  rename(Genotype = accession) %>%
  mutate(Genotype = trimws(Genotype)) %>%
  filter(!is.na(Genotype), Genotype != "") %>%
  filter(if_any(where(is.numeric), ~ !is.na(.)))

#############################################
# 2) Load genotype matrix
#############################################

geno_data <- read.csv(geno_file, row.names = 1, check.names = FALSE)
geno_matrix_full <- as.matrix(geno_data)
mode(geno_matrix_full) <- "numeric"

#############################################
# 3) Align and split training/test
#############################################

norm_id <- function(x) toupper(trimws(as.character(x)))

training_set   <- read.csv(training_file, header = TRUE)
training_names <- norm_id(training_set[[1]])

pheno_data$Genotype        <- norm_id(pheno_data$Genotype)
rownames(geno_matrix_full) <- norm_id(rownames(geno_matrix_full))

common_genotypes <- intersect(rownames(geno_matrix_full), pheno_data$Genotype)
pheno_filtered   <- pheno_data %>% filter(Genotype %in% common_genotypes)
rownames(pheno_filtered) <- pheno_filtered$Genotype

train_genotypes     <- intersect(training_names, common_genotypes)
geno_train          <- geno_matrix_full[train_genotypes, ]
pheno_train         <- pheno_filtered[train_genotypes, ]
pheno_train_numeric <- pheno_train %>% select(-Genotype)

message("Training N: ", length(train_genotypes))
# Should be 150 
#############################################
# 5) Cross-validation -- rrBLUP
# Note: CV runs on training set only (150 lines)
#############################################

source(xval_functions)

message("Running k-fold cross-validation (rrBLUP)...")
xval_rrblup <- k.xval(
  g.in       = geno_train,
  y.in       = pheno_train_numeric,
  y.trainset = pheno_train_numeric,
  k.fold     = 10,
  reps       = 50
)

# OUTPUT: rrBLUP cross-validation RDS
saveRDS(xval_rrblup,
        file.path(out_base, "xval_rrblup_STI_73traits_BLUEadjusted_kfold10.RDS"))

#############################################
# 6) Cross-validation -- Gaussian Kernel
#############################################

message("Computing Gaussian kernel...")
K      <- A.mat(geno_train)
k_dist <- dist(K)

message("Running k-fold cross-validation (Gaussian)...")
xval_gauss <- k.xval.GAUSS(
  g.in       = geno_train,
  y.in       = pheno_train,
  y.trainset = pheno_train,
  k_dist     = k_dist,
  k.fold     = 10,
  reps       = 50
)

# OUTPUT: Gaussian cross-validation RDS
saveRDS(xval_gauss,
        file.path(out_base, "xval_GAUSS_STI_73traits_BLUEadjusted_kfold10.RDS"))
message("Gaussian cross-validation saved.")

#############################################
# 7) Cross-validation -- Exponential Kernel
#############################################

message("Computing Exponential kernel...")
K_exp        <- Kernel_computation(X = geno_train, name = "exponential",
                                   degree = NULL, nL = NULL)
exp_dist_mat <- as.matrix(dist(K_exp))
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)

stopifnot(all(rownames(exp_dist_mat) == rownames(geno_train)))
stopifnot(all(rownames(exp_dist_mat) == pheno_train$Genotype))

message("Running k-fold cross-validation (Exponential)...")
xval_exp <- k.xval.EXP(
  g.in       = geno_train,
  y.in       = pheno_train,
  y.trainset = pheno_train,
  k_dist     = exp_dist_mat,
  k.fold     = 10,
  reps       = 50
)

# Add this after the xval_exp call in section 7
saveRDS(xval_exp,
        file.path(out_base, "xval_EXP_STI_73traits_BLUEadjusted_kfold10.RDS"))
message("Exponential cross-validation saved.")



#######################################################################################################################################
# Post HPC
#######################################################################################################################################
# 8) Plot prediction accuracies -- all traits
#############################################
# OUTPUT FOLDER
out_base <- "~/R/WorldVeg_Capsicum/STI_GS/"

rrblup_results <- readRDS(file.path(out_base, "xval_rrblup_STI_73traits_BLUEadjusted_kfold10.RDS"))$xval.result
gauss_results  <- readRDS(file.path(out_base, "xval_GAUSS_STI_73traits_BLUEadjusted_kfold10.RDS"))$xval.result
exp_results    <- readRDS(file.path(out_base, "xval_EXP_STI_73traits_BLUEadjusted_kfold10.RDS"))$xval.result

rrblup_results$model <- "rrBLUP"
gauss_results$model  <- "Gaussian"
exp_results$model    <- "Exponential"

all_models <- bind_rows(rrblup_results, gauss_results, exp_results) %>%
  mutate(
    r.mean = as.numeric(r.mean),
    r.sd   = as.numeric(r.sd),
    trait  = as.factor(trait),
    model  = factor(model, levels = c("rrBLUP", "Gaussian", "Exponential"))
  )

p_all <- ggplot(all_models, aes(x = model, y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
                width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(~trait, scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color = "red", linewidth = 1.2, linetype = "longdash") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    strip.text  = element_text(size = 10)
  ) +
  labs(title = "Prediction Accuracy (r) - 73 STI Traits (BLUE-adjusted)",
       x = "Model", y = "Prediction Accuracy (r)") +
  ylim(0, 1)

p_all

ggsave(
  filename = "~/R/WorldVeg_Capsicum/STI_GS/PA_all_traits_n423.pdf",
  plot = p_all,
  width = 24,
  height = 12,
  units = "in"
)


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
  "~/R/WorldVeg_Capsicum/STI_GS/PA_summary_73traits.csv"
)


#############################################
# 9) Plot -- r.mean > 0.5 in all models
#############################################

all_models_filtered <- all_models %>%
  group_by(trait) %>%
  filter(all(r.mean > 0.5)) %>%
  ungroup()

p_filtered <- ggplot(all_models_filtered, aes(x = model, y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
                width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(~trait, scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color = "red", linewidth = 1.2, linetype = "longdash") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    strip.text  = element_text(size = 10)
  ) +
  labs(title = "Prediction Accuracy (r) - Traits Where All Models > 0.5",
       x = "Model", y = "Prediction Accuracy (r)") +
  ylim(0, 1)

p_filtered

#list 
traits_above_05 <- all_models %>%
  group_by(trait) %>%
  filter(all(r.mean > 0.5)) %>%
  ungroup() %>%
  distinct(trait) %>%
  pull(trait)

traits_above_05

write.csv(
  data.frame(trait = traits_above_05),
  "~/R/WorldVeg_Capsicum/STI_GS/traits_above_0.5_all_models.csv",
  row.names = FALSE
)
