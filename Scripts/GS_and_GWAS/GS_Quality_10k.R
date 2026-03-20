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
library(rrBLUP)
library(dplyr)

##################
# GENOTYPE IMPORT (IMPUTED)
##################
geno_df_imputed <- read.csv(
  "~/R/World_Veg_Collab_Pepper/10250_GS/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv",
  row.names = 1,
  check.names = FALSE
)

##################
# PHENOTYPE IMPORT
##################
pheno_McLeod       <- read.csv("~/R/World_Veg_Collab_Pepper/McLeod_pheno_GS/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)


##################
# TRAINING SET FILTERING
##################
# Train = lines that have phenotype data
train_genotypes <- rownames(pheno_McLeod)

# Test = all remaining lines without phenotype data
test_genotypes <- setdiff(rownames(geno_df_imputed), train_genotypes)


geno_train <- geno_df_imputed[train_genotypes, ]
geno_test  <- geno_df_imputed[test_genotypes, ]

# Make sure matrices are numeric
geno_train <- as.matrix(sapply(geno_train, as.numeric))
rownames(geno_train) <- train_genotypes

geno_test  <- as.matrix(sapply(geno_test, as.numeric))
rownames(geno_test)  <- test_genotypes


# Phenotype data for training
pheno_train <- pheno_McLeod


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



run_rrBLUP_all_traits(geno_train, geno_test, pheno_train, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_McLeod_23traits_10k")

# Load both files
gebv_train <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_McLeod_23traits_10k_Train.csv", row.names = 1)
gebv_test  <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_McLeod_23traits_10k_Test.csv", row.names = 1)

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
write.csv(gebv_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_McLeod_23traits_10k_ALL.csv", row.names = TRUE)
###################################################################################################

###################################################################################################
###################################################################################################

# Load your combined 23-trait GEBVs
gebv_23traits <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_McLeod_23traits_10k_ALL.csv", row.names = 1)

# Bring Line ID back as a column
gebv_23traits$Line <- rownames(gebv_23traits)

library(readxl)
library(dplyr)

tripodi <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/10038_Tripodi_2021.csv")

tripodi_excel <- read_excel(
  "~/R/World_Veg_Collab_Pepper/10250_GS/pnas.2104315118.sd01.xlsx",
  sheet = "Dataset S1",
  skip = 1
)

# Country to Region map
country_to_region <- tibble::tibble(
  Country = c(
    "Japan", "China", "Taiwan", "Korea", "Mongolia"
  ),
  Region = "East Asia"
) %>%
  add_row(Country = c(
    "India", "Sri Lanka", "Nepal", "Bangladesh", "Bhutan", "Vietnam",
    "Thailand", "Philippines", "Indonesia", "Malaysia", "Malasya", "Myanmar", "Cambodia", "Laos", "Maldives"
  ), Region = "South/South East Asia") %>%
  add_row(Country = c(
    "Uzbekistan", "Tajikistan", "Kyrgyzstan", "AfGhananistan"
  ), Region = "Central Asia") %>%
  add_row(Country = c(
    "Hungary", "Czech rep", "Czech rep;Slovakia", "Slovakia", "Austria",
    "Germany", "Netherlands", "France", "Belgium", "Switzerland", "Denmark", "United Kingdom", "Portugal", "Spain", "Italy"
  ), Region = "Western/Central Europe") %>%
  add_row(Country = c(
    "Romania", "Russia", "Ukraine", "Serbia", "Bulgaria", "Croatia", "Slovenia", "Moldova", "Albania", "Yugoslavia", "Sembria (Montenegro)"
  ), Region = "Eastern Europe") %>%
  add_row(Country = c(
    "Ethiopia", "Congo", "Ivory Coast", "Togo", "Sudan", "Nigeria", "Tanzania",
    "South Africa", "Burkina Faso", "Zambia", "Zimbabwe", "Ghana", "Cameroon", "Niger",
    "Sierra Leone", "Uganda", "Angola", "Benin"
  ), Region = "Africa") %>%
  add_row(Country = c(
    "Israel", "Turkey", "Syria", "Iraq", "Iran", "Yemen"
  ), Region = "Middle East") %>%
  add_row(Country = c(
    "Mexico", "Guatemala", "Costa Rica", "Nicaragua", "El salvador", "Honduras", "Panama",
    "Puerto Rico", "Puerto rico", "Trinidad", "Tobago", "Saint Lucia", "Haiti", "Grenada",
    "Bahamas", "Cabo Verde", "Northern Mariana Islands", "Virgin Islands"
  ), Region = "Mesoamerica") %>%
  add_row(Country = c(
    "United States", "Canada"
  ), Region = "North America") %>%
  add_row(Country = c(
    "Perù", "Argentina", "Colombia", "Bolivia", "Brazil", "Ecuador", "Paraguay",
    "Chile", "Venezuela", "Uruguay", "Guyana", "Suriname Rep"
  ), Region = "South America") %>%
  add_row(Country = c(
    "Australia", "The Republic of Fiji", "Vanatu", "Solomon Islands", "Papua N Guinea"
  ), Region = "Oceania")


tripodi_geo <- tripodi %>%
  select(Sample_ID_G2P_Sol_database, material_source_geographic_location) %>%
  left_join(country_to_region, by = c("material_source_geographic_location" = "Country"))


# Make sure IDs are character
tripodi_geo$Sample_ID_G2P_Sol_database <- as.character(tripodi_geo$Sample_ID_G2P_Sol_database)
gebv_23traits$Line <- as.character(gebv_23traits$Line)

# Fix Line column: remove final "0" if there are 4 trailing numbers + extra 0
gebv_23traits$Line <- gsub("0$", "", gebv_23traits$Line)


# Merge
gebv_23traits_geo <- left_join(
  gebv_23traits,
  tripodi_geo[, c("Sample_ID_G2P_Sol_database", "Region")],
  by = c("Line" = "Sample_ID_G2P_Sol_database")
)


# Remove NAs
gebv_23traits_geo_clean <- gebv_23traits_geo %>%
  filter(!is.na(Region))

# Save if you want
write.csv(gebv_23traits_geo_clean, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_McLeod_23traits_10k_ALL_fixed_with_Region.csv", row.names = FALSE)



######################
#Plots
#########################
library(dplyr)
library(tidyr)
library(ggplot2)

# Summarize: mean GEBV for each Region and each Trait
gebv_summary <- gebv_23traits_geo_clean %>%
  group_by(Region) %>%
  summarize(across(starts_with("GEBV_"), mean, na.rm = TRUE))

# Make it longer format (Region, Trait, Mean_GEBV)
gebv_long <- gebv_summary %>%
  pivot_longer(
    cols = starts_with("GEBV_"),
    names_to = "Trait",
    values_to = "Mean_GEBV"
  )
############
#GEBVS by region for 23 traits
ggplot(gebv_long, aes(x = Trait, y = Region, fill = Mean_GEBV)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(fill = "Mean GEBV")

############
#by line
gebv_lines_long <- gebv_23traits_geo_clean %>%
  select(Line, starts_with("GEBV_")) %>%   # Keep Line and GEBVs
  pivot_longer(
    cols = starts_with("GEBV_"),
    names_to = "Trait",
    values_to = "GEBV"
  )

# No summarizing here — keep individual lines
gebv_long_lines <- gebv_23traits_geo_clean %>%
  pivot_longer(
    cols = starts_with("GEBV_"),
    names_to = "Trait",
    values_to = "GEBV"
  )
ggplot(gebv_long_lines, aes(x = Region, y = GEBV, fill = Region)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(title = "GEBV distributions by Region for each Trait",
       x = "Region",
       y = "GEBV")

############
#boxplot
gebv_long_lines <- gebv_23traits_geo_clean %>%
  pivot_longer(
    cols = starts_with("GEBV_"),
    names_to = "Trait",
    values_to = "GEBV"
  )

library(ggplot2)

ggplot(gebv_long_lines, aes(x = Trait, y = GEBV)) +
  geom_boxplot(outlier.size = 0.5, fill = "lightblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title.x = element_blank()
  ) +
  labs(
    title = "Distribution of GEBVs Across All Lines",
    y = "GEBV Value"
  )

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#Prediction accuracies for 10k for 12 quality traits

#########################################################################################################
# Genomic Selection Prediction Accuracy - 23 Traits (10K Dataset)
#########################################################################################################

library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
# STEP ONE: Load Data
####################

# Genotype data (imputed, rrBLUP format)
geno_data <- read.csv(
  "~/R/World_Veg_Project/data/inputs/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv",
  row.names = 1,
  check.names = FALSE
)

# Phenotype data (McLeod, 23 traits)
pheno_data <- read.csv(
  "~/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv",
  row.names = 1
)

# Training genotype names
train_genotypes <- read.csv(
  "~/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv",
  header = TRUE
)[[1]]

####################
# STEP TWO: Filter and Format
####################

common_genotypes <- intersect(rownames(geno_data), rownames(pheno_data))
geno_filtered <- geno_data[common_genotypes, ]
pheno_filtered <- pheno_data[common_genotypes, ]

geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))
geno_train <- geno_matrix[train_genotypes_filtered, ]
pheno_train <- pheno_filtered[train_genotypes_filtered, ]
pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

####################
# STEP THREE: Run rrBLUP CV
####################

source("~/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

pheno_train_numeric <- pheno_train %>% select(-Genotype)

xval_k10_rrblup <- k.xval(
  g.in = geno_train,
  y.in = pheno_train_numeric,
  y.trainset = pheno_train_numeric,
  k.fold = 10,
  reps = 50
)

saveRDS(xval_k10_rrblup, "~/R/World_Veg_Project/data/outputs/xval_rrblup_23traits_kfold_10_10k.RDS")

####################
# STEP FOUR: Gaussian Kernel
####################

K <- A.mat(geno_train)
k_dist <- dist(K)

xval_k10_GAUSS <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = pheno_train,
  y.trainset = pheno_train,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)

saveRDS(xval_k10_GAUSS, "~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_23traits_kfold_10_10k.RDS")

####################
# STEP FIVE: Exponential Kernel
####################

K.Exp <- Kernel_computation(X = geno_train, name = "exponential", degree = 1)

exp_dist <- dist(K.Exp)
exp_dist_mat <- as.matrix(exp_dist)
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)

stopifnot(all(rownames(exp_dist_mat) == rownames(geno_train)))
stopifnot(all(rownames(exp_dist_mat) == pheno_train$Genotype))

xval_k10_EXP <- k.xval.EXP(
  g.in = geno_train,
  y.in = pheno_train,
  y.trainset = pheno_train,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)

saveRDS(xval_k10_EXP, "~/R/World_Veg_Collab_Pepper/outputs/xval_EXP_23traits_kfold_10_10k.RDS")

###################################################################################################
#Plot prediction accuracies
###################################################################################################

library(ggplot2)
library(dplyr)
library(readr)

# Load cross-validation results (these should contain the $xval.result list)
rrblup_results <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_10k/xval_rrblup_23traits_kfold_10_10k.RDS")$xval.result
gauss_results  <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_10k/xval_GAUSS_23traits_kfold_10_10k.RDS")$xval.result
exp_results    <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_10k/xval_EXP_23traits_kfold_10_10k.RDS")$xval.result

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


###################################################################################################