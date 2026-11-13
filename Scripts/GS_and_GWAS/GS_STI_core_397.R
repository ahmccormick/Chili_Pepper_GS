# Supplemental Figure S28
#####################################################################################################################
# Core STI: BEAGLE genotypes → GS → prediction accuracies
# All matching STI-phenotyped accessions
# Original CV functions: 10-fold × 50 repetitions

library(vcfR)
library(rrBLUP)
library(dplyr)
library(tidyr)
library(ggplot2)

####################################################################
# 1. Paths
####################################################################

folder <- "/home/ahmccorm/scratch/world_veg/STI/core"

output_folder <- file.path(folder, "outputs")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(
  folder, "Pepper_397_annuum_LD_pruned_beagle.vcf.gz"
)

sti_file <- file.path(folder, "STI_from_BLUEs_all_73_traits.csv")
xval_functions <- file.path(folder, "xval_kfold_functions.R")

stopifnot(
  file.exists(vcf_file),
  file.exists(sti_file),
  file.exists(xval_functions)
)

####################################################################
# 2. Read BEAGLE genotypes
####################################################################

vcf <- read.vcfR(vcf_file)

gt <- extract.gt(vcf, element = "GT")
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

# Match accession IDs consistently
norm_id <- function(x) {
  toupper(trimws(as.character(x)))
}

rownames(geno) <- norm_id(rownames(geno))

stopifnot(
  !anyNA(geno),
  !anyDuplicated(rownames(geno))
)

cat("Genotyped accessions:", nrow(geno), "\n")
cat("Markers:", ncol(geno), "\n")

rm(vcf, gt, geno_numeric)
gc()

####################################################################
# 3. Read STI phenotypes and align by accession
####################################################################

pheno <- read.csv(
  sti_file,
  check.names = FALSE,
  na.strings = c("NA", "", "NaN", " ")
)

stopifnot("accession" %in% names(pheno))

pheno$accession <- norm_id(pheno$accession)

stopifnot(
  !anyNA(pheno$accession),
  all(nzchar(pheno$accession)),
  !anyDuplicated(pheno$accession)
)

rownames(pheno) <- pheno$accession
pheno$accession <- NULL

stopifnot(
  ncol(pheno) == 73,
  all(vapply(pheno, is.numeric, logical(1)))
)

# Treat non-finite STI values as missing
pheno[] <- lapply(pheno, function(y) {
  y[!is.finite(y)] <- NA_real_
  y
})

matched_ids <- intersect(rownames(geno), rownames(pheno))
unmatched_ids <- setdiff(rownames(pheno), rownames(geno))

cat("STI phenotype records:", nrow(pheno), "\n")
cat("Matching phenotype records:", length(matched_ids), "\n")
cat("Unmatched phenotype records:", length(unmatched_ids), "\n")

write.csv(
  data.frame(Sample_ID = unmatched_ids),
  file.path(output_folder, "Unmatched_STI_accessions.csv"),
  row.names = FALSE
)

pheno_train <- pheno[matched_ids, , drop = FALSE]
geno_train <- geno[matched_ids, , drop = FALSE]

stopifnot(
  identical(rownames(geno_train), rownames(pheno_train))
)

write.csv(
  data.frame(Sample_ID = matched_ids),
  file.path(output_folder, "STI_matching_accessions.csv"),
  row.names = FALSE
)

####################################################################
# 4. Final rrBLUP GS models
# Use observed STI values per trait; predict all VCF accessions
####################################################################

traits <- names(pheno_train)

gebv <- matrix(
  NA_real_,
  nrow = nrow(geno),
  ncol = length(traits),
  dimnames = list(rownames(geno), traits)
)

summary_list <- list()
training_list <- list()

for (trait in traits) {
  
  observed <- is.finite(pheno_train[[trait]])
  train_ids <- matched_ids[observed]
  y <- pheno_train[[trait]][observed]
  
  training_list[[trait]] <- data.frame(
    Trait = rep(trait, length(train_ids)),
    Sample_ID = train_ids
  )
  
  status <- "Fitted"
  
  if (length(y) < 3) {
    
    status <- "Skipped: fewer than 3 observations"
    
  } else if (var(y) == 0) {
    
    status <- "Skipped: zero phenotype variance"
    
  } else {
    
    cat("Fitting final STI GS model:", trait,
        "- training n =", length(y), "\n")
    
    Z <- geno[train_ids, , drop = FALSE]
    
    # Remove markers invariant in this trait's training set
    keep <- colSums(
      sweep(Z, 2, Z[1, ], FUN = "!=")
    ) > 0
    
    if (!any(keep)) {
      
      status <- "Skipped: no variable markers"
      
    } else {
      
      marker_means <- colMeans(Z[, keep, drop = FALSE])
      
      Z <- sweep(
        Z[, keep, drop = FALSE],
        2, marker_means, FUN = "-"
      )
      
      model <- mixed.solve(y = y, Z = Z)
      
      # Centered genetic contribution for every accession
      gebv[, trait] <- as.numeric(
        geno[, keep, drop = FALSE] %*% model$u
      ) - sum(marker_means * as.numeric(model$u))
      
      rm(model)
    }
    
    rm(Z)
  }
  
  summary_list[[trait]] <- data.frame(
    Trait = trait,
    N_training = length(y),
    N_without_observed_STI = nrow(geno) - length(y),
    Status = status
  )
}

training_summary <- bind_rows(summary_list)

write.csv(
  data.frame(
    Sample_ID = rownames(geno),
    gebv,
    row.names = NULL,
    check.names = FALSE
  ),
  file.path(
    output_folder,
    paste0("GEBVs_STI_all_", nrow(geno), ".csv")
  ),
  row.names = FALSE
)

write.csv(
  training_summary,
  file.path(output_folder, "STI_training_summary.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(training_list),
  file.path(output_folder, "STI_training_IDs_by_trait.csv"),
  row.names = FALSE
)

rm(gebv)
gc()

####################################################################
# 5. Prepare phenotypes for original CV functions
####################################################################

# Exclude traits with too few observations for 10-fold CV,
# zero variance, or an unsuccessful final GS fit.
cv_traits <- training_summary %>%
  filter(Status == "Fitted", N_training >= 10) %>%
  pull(Trait)

if (length(cv_traits) == 0) {
  stop("No traits suitable for cross-validation.")
}

pheno_cv <- pheno_train[, cv_traits, drop = FALSE]

# Preserve NA values: original functions exclude them per trait
pheno_kernel <- data.frame(
  Genotype = rownames(pheno_cv),
  pheno_cv,
  row.names = rownames(pheno_cv),
  check.names = FALSE
)

source(xval_functions)

####################################################################
# 6. rrBLUP prediction accuracies
####################################################################

set.seed(123)

cat("\nStarting rrBLUP STI cross-validation...\n")

xval_rrblup <- k.xval(
  g.in = geno_train,
  y.in = pheno_cv,
  y.trainset = pheno_cv,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_rrblup,
  file.path(output_folder, "xval_rrblup_STI_core397_kfold10.RDS")
)

####################################################################
# 7. Gaussian: original construction
####################################################################

# A.mat requires -1/0/1 coding
K <- A.mat(geno_train - 1)
k_dist <- dist(K)

set.seed(123)

cat("\nStarting Gaussian STI cross-validation...\n")

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
  file.path(output_folder, "xval_GAUSS_STI_core397_kfold10.RDS")
)

####################################################################
# 8. Exponential: original implementation
####################################################################

K_exp <- Kernel_computation(
  X = geno_train,
  name = "exponential",
  degree = NULL,
  nL = NULL
)

exp_dist_mat <- as.matrix(dist(K_exp))

rownames(exp_dist_mat) <- colnames(exp_dist_mat) <-
  rownames(geno_train)

stopifnot(
  identical(rownames(exp_dist_mat), pheno_kernel$Genotype)
)

set.seed(123)

cat("\nStarting Exponential STI cross-validation...\n")

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
  file.path(output_folder, "xval_EXP_STI_core397_kfold10.RDS")
)

####################################################################
# 9. Combine PA results
####################################################################

rrblup_results <- xval_rrblup$xval.result
gauss_results <- xval_gauss$xval.result
exp_results <- xval_exp$xval.result

rrblup_results$model <- "rrBLUP"
gauss_results$model <- "Gaussian"
exp_results$model <- "Exponential"

all_models <- bind_rows(
  rrblup_results,
  gauss_results,
  exp_results
) %>%
  mutate(
    trait = as.character(trait),
    r.mean = as.numeric(as.character(r.mean)),
    r.sd = as.numeric(as.character(r.sd)),
    model = factor(
      model,
      levels = c("rrBLUP", "Gaussian", "Exponential")
    )
  ) %>%
  filter(
    !is.na(trait),
    nzchar(trait),
    is.finite(r.mean),
    is.finite(r.sd)
  )

pa_table <- all_models %>%
  select(trait, model, r.mean, r.sd) %>%
  pivot_wider(
    names_from = model,
    values_from = c(r.mean, r.sd)
  )

write.csv(
  pa_table,
  file.path(output_folder, "STI_PA_summary_core397.csv"),
  row.names = FALSE
)

####################################################################
# 10. Plot prediction accuracies
####################################################################

p <- ggplot(
  all_models,
  aes(x = model, y = r.mean, color = model)
) +
  theme_bw() +
  geom_errorbar(
    aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
    width = 0.3
  ) +
  geom_point(size = 3) +
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
    title = paste0(
      "Core STI Prediction Accuracies — ",
      n_distinct(all_models$trait), " Traits"
    ),
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  coord_cartesian(ylim = c(0, 1))

ggsave(
  filename = file.path(output_folder, "STI_PA_all_traits_core397.pdf"),
  plot = p,
  width = 24,
  height = 12,
  units = "in"
)

# Require all three models to pass the threshold
high_traits <- all_models %>%
  group_by(trait) %>%
  summarise(
    passes = n_distinct(model) == 3 && all(r.mean > 0.5),
    .groups = "drop"
  ) %>%
  filter(passes) %>%
  pull(trait)

cat("\nTraits with PA > 0.5 for all three models:",
    length(high_traits), "\n")

write.csv(
  data.frame(Trait = high_traits),
  file.path(output_folder, "STI_traits_all_models_above_0.5.csv"),
  row.names = FALSE
)

cat("\nSTI GS and PA analyses completed.\n")

#####################################################################################################################
