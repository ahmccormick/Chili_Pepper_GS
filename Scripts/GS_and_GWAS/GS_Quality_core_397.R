# Supplemental Figure S20
#####################################################################################################################
# Quality traits for core collection (HPC)

library(vcfR)
library(rrBLUP)

####################################################################
# 1. Paths
####################################################################

folder <- "/home/ahmccorm/scratch/world_veg/quality/core_n397_quality"

vcf_file <- "Pepper_397_annuum_LD_pruned_beagle.vcf.gz"

pheno_file <- "McLeod_per_line_pheno_data_n329_23traits.csv"

output_folder <- file.path(folder, "outputs")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

####################################################################
# 2. Read BEAGLE genotypes
####################################################################

vcf <- read.vcfR(file.path(folder, vcf_file))
gt <- extract.gt(vcf, element = "GT")

# Accept phased and unphased calls
gt[] <- gsub("|", "/", gt, fixed = TRUE)

geno_numeric <- matrix(
  NA_real_,
  nrow = nrow(gt),
  ncol = ncol(gt),
  dimnames = dimnames(gt)
)

geno_numeric[which(gt == "0/0")] <- 0
geno_numeric[which(gt == "0/1" | gt == "1/0")] <- 1
geno_numeric[which(gt == "1/1")] <- 2

geno <- t(geno_numeric)

stopifnot(
  !anyNA(geno),
  !anyDuplicated(rownames(geno))
)

cat("Genotyped accessions:", nrow(geno), "\n")
cat("SNPs:", ncol(geno), "\n")

rm(vcf, gt, geno_numeric)
gc()

####################################################################
# 3. Read quality phenotypes and match accession IDs
####################################################################

pheno_data <- read.csv(
  file.path(folder, pheno_file),
  row.names = 1,
  check.names = FALSE
)

pheno_data <- pheno_data[
  , colnames(pheno_data) != "X",
  drop = FALSE
]

stopifnot(
  !anyDuplicated(rownames(pheno_data)),
  all(vapply(pheno_data, is.numeric, logical(1))),
  ncol(pheno_data) == 23
)

train_ids <- intersect(rownames(geno), rownames(pheno_data))
unmatched_ids <- setdiff(rownames(pheno_data), rownames(geno))

cat("Phenotype accessions:", nrow(pheno_data), "\n")
cat("Matching training accessions:", length(train_ids), "\n")
cat("Unmatched phenotype accessions:", length(unmatched_ids), "\n")

write.csv(
  data.frame(Sample_ID = unmatched_ids),
  file.path(output_folder, "Unmatched_quality_phenotype_IDs.csv"),
  row.names = FALSE
)

pheno_train <- pheno_data[train_ids, , drop = FALSE]
geno_train <- geno[train_ids, , drop = FALSE]

stopifnot(
  identical(rownames(geno_train), rownames(pheno_train))
)

# This prepared n329 file is expected to contain complete phenotypes.
# Stop rather than silently impute missing phenotype values.
if (any(!is.finite(as.matrix(pheno_train)))) {
  stop("Quality phenotype file contains missing/non-finite values.")
}

if (nrow(pheno_train) < 10) {
  stop("Too few matching accessions for 10-fold CV.")
}

variable_traits <- vapply(
  pheno_train,
  function(y) var(y) > 0,
  logical(1)
)

if (any(!variable_traits)) {
  message(
    "Skipping zero-variance traits: ",
    paste(names(variable_traits)[!variable_traits], collapse = ", ")
  )
  pheno_train <- pheno_train[, variable_traits, drop = FALSE]
}

if (ncol(pheno_train) == 0) {
  stop("No variable phenotype traits remain.")
}

write.csv(
  data.frame(Sample_ID = train_ids),
  file.path(output_folder, "Quality_training_IDs.csv"),
  row.names = FALSE
)

####################################################################
# 4. Final GS models: GEBVs for all VCF accessions
####################################################################

gebv <- matrix(
  NA_real_,
  nrow = nrow(geno),
  ncol = ncol(pheno_train),
  dimnames = list(rownames(geno), colnames(pheno_train))
)

# Remove SNPs invariant in the training population
keep <- colSums(
  sweep(geno_train, 2, geno_train[1, ], FUN = "!=")
) > 0

if (!any(keep)) {
  stop("No variable SNPs in the training population.")
}

marker_means <- colMeans(geno_train[, keep, drop = FALSE])

Z_train <- sweep(
  geno_train[, keep, drop = FALSE],
  2, marker_means, FUN = "-"
)

for (trait in colnames(pheno_train)) {
  
  cat("Fitting final GS model:", trait, "\n")
  
  model <- mixed.solve(
    y = pheno_train[[trait]],
    Z = Z_train
  )
  
  gebv[, trait] <- as.numeric(
    geno[, keep, drop = FALSE] %*% model$u
  ) - sum(marker_means * as.numeric(model$u))
}

write.csv(
  data.frame(
    Sample_ID = rownames(geno),
    Set = ifelse(rownames(geno) %in% train_ids, "Train", "Test"),
    gebv,
    row.names = NULL,
    check.names = FALSE
  ),
  file.path(
    output_folder,
    paste0("GEBVs_Quality_all_", nrow(geno), ".csv")
  ),
  row.names = FALSE
)

rm(Z_train, gebv, model)
gc()

####################################################################
# 5. Prediction accuracies: original custom functions
####################################################################

source(file.path(folder, "xval_kfold_functions.R"))

# rrBLUP receives numeric phenotype columns only
set.seed(123)

xval_rrblup <- k.xval(
  g.in = geno_train,
  y.in = pheno_train,
  y.trainset = pheno_train,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_rrblup,
  file.path(output_folder, "xval_rrblup_23traits_kfold_10.RDS")
)

####################################################################
# 6. Gaussian: original construction
####################################################################

pheno_kernel <- data.frame(
  Genotype = rownames(pheno_train),
  pheno_train,
  row.names = rownames(pheno_train),
  check.names = FALSE
)

# A.mat requires -1/0/1 coding
K <- A.mat(geno_train - 1)
k_dist <- dist(K)

set.seed(123)

xval_gauss <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = pheno_kernel,
  y.trainset = pheno_kernel,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_gauss,
  file.path(output_folder, "xval_GAUSS_23traits_kfold_10.RDS")
)

####################################################################
# 7. Exponential: original implementation
####################################################################

K.Exp <- Kernel_computation(
  X = geno_train,
  name = "exponential",
  degree = NULL,
  nL = NULL
)

exp_dist_mat <- as.matrix(dist(K.Exp))

rownames(exp_dist_mat) <- colnames(exp_dist_mat) <-
  rownames(geno_train)

set.seed(123)

xval_exp <- k.xval.EXP(
  g.in = geno_train,
  y.in = pheno_kernel,
  y.trainset = pheno_kernel,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_exp,
  file.path(output_folder, "xval_EXP_23traits_kfold_10.RDS")
)

cat("\nQuality GS and PA analyses completed.\n")



###############
# PLOT
###############
library(ggplot2)
library(dplyr)

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/quality/core"

# Load prediction accuracy results
rrblup_results <- readRDS(
  file.path(folder, "xval_rrblup_23traits_kfold_10.RDS")
)$xval.result

gauss_results <- readRDS(
  file.path(folder, "xval_GAUSS_23traits_kfold_10.RDS")
)$xval.result

exp_results <- readRDS(
  file.path(folder, "xval_EXP_23traits_kfold_10.RDS")
)$xval.result

# Add model labels
rrblup_results$model <- "rrBLUP"
gauss_results$model <- "Gaussian"
exp_results$model <- "Exponential"

# Combine and format results
all_models_core <- bind_rows(
  rrblup_results, gauss_results, exp_results
) %>%
  mutate(
    r.mean = as.numeric(as.character(r.mean)),
    r.sd = as.numeric(as.character(r.sd)),
    trait = as.factor(trait),
    model = factor(
      model,
      levels = c("rrBLUP", "Gaussian", "Exponential")
    )
  ) %>%
  filter(
    !is.na(trait),
    is.finite(r.mean),
    is.finite(r.sd)
  )

# Plot
p_core <- ggplot(
  all_models_core,
  aes(x = model, y = r.mean, color = model)
) +
  theme_bw() +
  geom_errorbar(
    aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
    width = 0.3,
    position = position_dodge(0.2)
  ) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(~ trait, scales = "free_x", nrow = 5) +
  geom_hline(
    yintercept = 0.5,
    color = "red",
    linewidth = 1.2,
    linetype = "longdash"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 10
    ),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 10)
  ) +
  labs(
    title = "Prediction Accuracy (r) for Quality Traits — Core Collection",
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  coord_cartesian(ylim = c(0, 1))

print(p_core)

ggsave(
  filename = file.path(folder, "prediction_accuracy_23traits_core.pdf"),
  plot = p_core,
  width = 20,
  height = 12,
  units = "in"
)

#####################################################################################################################
