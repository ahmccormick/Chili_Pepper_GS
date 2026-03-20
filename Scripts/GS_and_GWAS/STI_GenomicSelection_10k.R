#############################################
# STI Genomic Selection - 10k Global Collection
# Model: rrBLUP
# Phenotype: pre-computed STI from BLUE-adjusted data (core 423)
# Predicts GEBVs for all 10k lines
#############################################

rm(list = ls())

library(rrBLUP)
library(dplyr)
library(readr)

#############################################
# 0) File paths -- EDIT THESE
#############################################

geno_file      <- "~/R/World_Veg_Project/STI_GS_10k/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv"
sti_file       <- "~/R/World_Veg_Project/anna_BLUES/anna_BLUES_rebuild/STI_from_BLUEs_all_73_traits.csv"
training_file  <- "~/R/World_Veg_Project/STI_GS_10k/overlapping_sample_ids_282.csv"

# OUTPUT FOLDER
out_base <- "~/R/WorldVeg_Capsicum/STI_GS_10k/"

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
message("Traits: ", ncol(pheno_data) - 1)

# Save cleaned phenotype file
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
write_csv(pheno_data,
          file.path(out_base, "cleaned_pheno_STI_73traits_BLUEadjusted_10k.csv"))

#############################################
# 2) Load genotype matrix (10k)
#############################################

geno_df_imputed <- read.csv(geno_file, row.names = 1, check.names = FALSE)
geno_matrix_full <- as.matrix(geno_df_imputed)
mode(geno_matrix_full) <- "numeric"
rownames(geno_matrix_full) <- trimws(rownames(geno_matrix_full))

message("Genotype matrix: ", nrow(geno_matrix_full), " samples x ",
        ncol(geno_matrix_full), " markers")

#############################################
# 3) Define training and test sets
#############################################

norm_id <- function(x) toupper(trimws(as.character(x)))

# normalize IDs
pheno_data$Genotype      <- norm_id(pheno_data$Genotype)
rownames(geno_matrix_full) <- norm_id(rownames(geno_matrix_full))

# load training IDs
training_ids <- norm_id(read.csv(training_file, header = TRUE)[[1]])

# training = overlap between training file, genotype matrix, and phenotype data
train_genotypes <- intersect(training_ids,
                             intersect(rownames(geno_matrix_full),
                                       pheno_data$Genotype))

# test = ALL remaining lines in full 10k genotype matrix
test_genotypes  <- setdiff(rownames(geno_matrix_full), train_genotypes)

message("Training N: ", length(train_genotypes))
message("Test N (all remaining 10k): ", length(test_genotypes))

# save training set used
write.csv(data.frame(Genotype = train_genotypes),
          file.path(out_base, "training_set_BLUEadjusted.csv"),
          row.names = FALSE)

# subset genotype matrices
geno_train <- geno_matrix_full[train_genotypes, ]
geno_test  <- geno_matrix_full[test_genotypes,  ]

# subset phenotype to training lines only
rownames(pheno_data) <- pheno_data$Genotype
pheno_train          <- pheno_data[train_genotypes, ] %>% select(-Genotype)

message("Phenotype training matrix: ", nrow(pheno_train), " lines x ",
        ncol(pheno_train), " traits")

#############################################
# 4) Run rrBLUP and predict GEBVs
#############################################

trait_names      <- colnames(pheno_train)
GEBVs_train_list <- list()
GEBVs_test_list  <- list()

for (trait in trait_names) {
  
  y <- pheno_train[[trait]]
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

# OUTPUT: separate train and test GEBV files
write.csv(GEBVs_train_final,
          file.path(out_base, "GEBVs_STI_73traits_BLUEadjusted_10k_Train.csv"),
          row.names = TRUE)
write.csv(GEBVs_test_final,
          file.path(out_base, "GEBVs_STI_73traits_BLUEadjusted_10k_Test.csv"),
          row.names = TRUE)

#############################################
# 5) Combine train + test into single file
#############################################

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

# OUTPUT: combined GEBV file (all 10k lines)
write.csv(gebv_all,
          file.path(out_base, "GEBVs_STI_73traits_BLUEadjusted_ALL_10k.csv"),
          row.names = TRUE)




####################################################################################################################################################################
####################################################################################################################################################################
# PAS
####################################################################################################################################################################
#local file paths
geno_file      <- "~/R/World_Veg_Project/STI_GS_10k/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv"
training_file  <- "~/R/World_Veg_Project/STI_GS_10k/overlapping_sample_ids_282.csv"
xval_functions <- "/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R"
####################################################################################################################################################################
# STI PAs - 10k dataset
# Phenotype: BLUE-adjusted STI (73 traits)
# Training: overlapping_sample_ids_282
####################################################################################################################################################################
rm(list = ls())

library(rrBLUP)
library(dplyr)
library(readr)
library(ggplot2)

# HPC graphics fix
options(bitmapType = "cairo")

#############################################
# 0) File paths -- EDIT THESE
#############################################

geno_file      <- "/home/ahmccorm/kantar_koastore/anna/WorldVeg_Capsicum/STI_GS_10k/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv"
sti_file       <- "/home/ahmccorm/kantar_koastore/anna/WorldVeg_Capsicum/STI_GS_10k/STI_from_BLUEs_all_73_traits.csv"
training_file  <- "/home/ahmccorm/kantar_koastore/anna/WorldVeg_Capsicum/STI_GS_10k/overlapping_sample_ids_282.csv"
xval_functions <- "/home/ahmccorm/kantar_koastore/anna/WorldVeg_Capsicum/STI_GS_10k/xval_kfold_functions.R"
out_base       <- "/home/ahmccorm/kantar_koastore/anna/WorldVeg_Capsicum/STI_GS_10k/"

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

message("Phenotype rows: ", nrow(pheno_data))
message("Traits: ", ncol(pheno_data) - 1)

#############################################
# 2) Load genotype matrix
#############################################

geno_data <- read.csv(geno_file, row.names = 1, check.names = FALSE)
geno_matrix <- as.matrix(geno_data)
mode(geno_matrix) <- "numeric"
rownames(geno_matrix) <- trimws(rownames(geno_matrix))

message("Genotype matrix: ", nrow(geno_matrix), " samples x ", ncol(geno_matrix), " markers")

#############################################
# 3) Align and define training set
#############################################

norm_id <- function(x) toupper(trimws(as.character(x)))

training_ids <- norm_id(read.csv(training_file, header = TRUE)[[1]])

pheno_data$Genotype  <- norm_id(pheno_data$Genotype)
rownames(geno_matrix) <- norm_id(rownames(geno_matrix))

# common genotypes across geno and pheno
common_genotypes <- intersect(rownames(geno_matrix), pheno_data$Genotype)
geno_filtered    <- geno_matrix[common_genotypes, ]
pheno_filtered   <- pheno_data %>% filter(Genotype %in% common_genotypes)
rownames(pheno_filtered) <- pheno_filtered$Genotype

# training = overlap of training IDs with common genotypes
train_genotypes <- intersect(training_ids, common_genotypes)

geno_train          <- geno_filtered[train_genotypes, ]
pheno_train         <- pheno_filtered[train_genotypes, ]
pheno_train_numeric <- pheno_train %>% select(-Genotype)

message("Training N: ", length(train_genotypes))

#############################################
# 4) Cross-validation -- rrBLUP
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

saveRDS(xval_rrblup,
        file.path(out_base, "xval_rrblup_STI_73traits_BLUEadjusted_10k_kfold10.RDS"))
message("rrBLUP cross-validation saved.")

#############################################
# 5) Cross-validation -- Gaussian Kernel
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

saveRDS(xval_gauss,
        file.path(out_base, "xval_GAUSS_STI_73traits_BLUEadjusted_10k_kfold10.RDS"))
message("Gaussian cross-validation saved.")

#############################################
# 6) Cross-validation -- Exponential Kernel
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

saveRDS(xval_exp,
        file.path(out_base, "xval_EXP_STI_73traits_BLUEadjusted_10k_kfold10.RDS"))
message("Exponential cross-validation saved.")

message("\nAll done. RDS files saved to: ", out_base)

##################################################################################################################################################################
###################################################################################################
#Plot prediction accuracies
###################################################################################################

library(ggplot2)
library(dplyr)
library(readr)

# Load cross-validation results (these should contain the $xval.result list)
rrblup_results <- readRDS("~/R/WorldVeg_Capsicum/STI_GS_10k/xval_rrblup_STI_73traits_BLUEadjusted_10k_kfold10.RDS")$xval.result
gauss_results  <- readRDS("~/R/WorldVeg_Capsicum/STI_GS_10k/xval_GAUSS_STI_73traits_BLUEadjusted_10k_kfold10.RDS")$xval.result
exp_results    <- readRDS("~/R/WorldVeg_Capsicum/STI_GS_10k/xval_EXP_STI_73traits_BLUEadjusted_10k_kfold10.RDS")$xval.result

# Add model labels
rrblup_results$model <- "rrBLUP"
gauss_results$model  <- "Gaussian"
exp_results$model    <- "Exponential"

# Combine into one dataframe
all_models <- bind_rows(rrblup_results, gauss_results, exp_results)
#all_models <- bind_rows(gauss_results, exp_results)

# Ensure correct data types
all_models <- all_models %>%
  mutate(
    r.mean = as.numeric(r.mean),
    r.sd = as.numeric(r.sd),
    trait = as.factor(trait),
    model = factor(model, levels = c("rrBLUP", "Gaussian", "Exponential"))
  )

# Plot
p<- ggplot(all_models, aes(x = model, y = r.mean, color = model)) +
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
    title = "Prediction Accuracy (r) for 73 Traits",
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  ylim(0, 1)

p

ggsave(
  filename = "~/R/WorldVeg_Capsicum/STI_GS_10k/PA_all_traits.pdf",
  plot = p,
  width = 24,
  height = 12,
  units = "in"
)

###################################################################################################
# Keep only traits where ALL models have r.mean > 0.5
all_models_filtered <- all_models %>%
  group_by(trait) %>%
  filter(all(r.mean > 0.5)) %>%
  ungroup()

# Plot only the filtered traits
ggplot(all_models_filtered, aes(x = model, y = r.mean, color = model)) +
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
    title = "Prediction Accuracy (r) for Traits Where All Models > 0.5",
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  ylim(0, 1)

pa_table <- all_models %>%
  select(trait, model, r.mean, r.sd) %>%
  arrange(trait, model) %>%
  tidyr::pivot_wider(
    names_from = model,
    values_from = c(r.mean, r.sd)
  )


# save PAs as table
pa_table <- all_models %>%
  select(trait, model, r.mean, r.sd) %>%
  arrange(trait, model) %>%
  tidyr::pivot_wider(
    names_from = model,
    values_from = c(r.mean, r.sd)
  )


readr::write_csv(
  pa_table,
  "~/R/WorldVeg_Capsicum/STI_GS_10k/PA_summary_all_traits.csv"
)


#list
traits_above_05 <- all_models_filtered %>%
  distinct(trait) %>%
  pull(trait)

traits_above_05

readr::write_csv(
  data.frame(trait = traits_above_05),
  "~/R/WorldVeg_Capsicum/STI_GS_10k/traits_above_0.5_all_models.csv"
)

