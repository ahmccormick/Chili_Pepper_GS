####################################################################################################
# PVE (percent variance explained) for Table S7 SNPs
# Method: incremental r² from regression of STI ~ PC1-5 + SNP, partialling out population structure
####################################################################################################

library(tidyverse)
library(data.table)

#############################################
# 0) Paths
#############################################

geno_path  <- "~/R/WorldVeg_Capsicum/STI_GWAS/genotype_matrix_imputed.csv"
blues_path <- "~/R/WorldVeg_Capsicum/STI_GWAS/STI_from_BLUEs_all_73_traits.csv"
s7_path    <- "~/R/WorldVeg_Capsicum/STI_GWAS/GWAS_tables/TableS7_multitrait_SNPs_3plus_models_annotated.csv"
pc_path    <- "~/R/WorldVeg_Capsicum/STI_GWAS/MLM_PC_K/MLM_PC_K/green/GAPIT.Genotype.PCA.csv"
outdir     <- "~/R/WorldVeg_Capsicum/STI_GWAS/GWAS_tables/"

#############################################
# 1) Load S7 and get SNP x trait pairs
#############################################

s7 <- read.csv(s7_path, stringsAsFactors = FALSE) %>%
  select(-any_of(c("X", "X.1")))

snp_trait_pairs <- s7 %>%
  select(SNP, Traits) %>%
  mutate(Trait = strsplit(Traits, ", ")) %>%
  unnest(Trait) %>%
  select(SNP, Trait) %>%
  distinct()

message("SNP x trait pairs: ", nrow(snp_trait_pairs))

#############################################
# 2) Load genotype matrix, subset to S7 SNPs
#############################################

message("Loading genotype matrix...")
geno <- fread(geno_path, data.table = FALSE)
rownames(geno) <- geno[[1]]
geno <- geno[, -1]

s7_snps       <- unique(s7$SNP)
snps_found    <- intersect(s7_snps, colnames(geno))
snps_missing  <- setdiff(s7_snps, colnames(geno))

if (length(snps_missing) > 0) {
  message("WARNING: ", length(snps_missing), " SNPs not found in genotype matrix:")
  print(snps_missing)
}

geno_s7 <- geno[, snps_found, drop = FALSE]
message("Genotype subset: ", nrow(geno_s7), " taxa x ", ncol(geno_s7), " SNPs")
rm(geno); gc()

#############################################
# 3) Load STI values and PCs
#############################################

blues <- read.csv(blues_path, stringsAsFactors = FALSE, row.names = 1)
message("STI values: ", nrow(blues), " taxa x ", ncol(blues), " traits")

pcs <- read.csv(pc_path, stringsAsFactors = FALSE) %>%
  select(taxa, PC1, PC2, PC3, PC4, PC5) %>%
  column_to_rownames("taxa")
message("PCs loaded: ", nrow(pcs), " taxa")

# check traits
traits_missing <- setdiff(unique(snp_trait_pairs$Trait), colnames(blues))
if (length(traits_missing) > 0) {
  message("WARNING: traits not in STI file: ", paste(traits_missing, collapse = ", "))
  snp_trait_pairs <- snp_trait_pairs %>% filter(Trait %in% colnames(blues))
}

#############################################
# 4) Align taxa across all three matrices
#############################################

common_taxa <- intersect(rownames(geno_s7), rownames(blues))
common_taxa <- intersect(common_taxa, rownames(pcs))
message("Taxa in common: ", length(common_taxa))

geno_s7 <- geno_s7[common_taxa, , drop = FALSE]
blues   <- blues[common_taxa, , drop = FALSE]
pcs     <- pcs[common_taxa, , drop = FALSE]

#############################################
# 5) PC-adjusted PVE function
#############################################

compute_pve_adjusted <- function(snp_id, trait_id, geno_mat, blues_mat, pc_mat) {
  
  if (!snp_id %in% colnames(geno_mat)) return(
    tibble(SNP = snp_id, Trait = trait_id, n = NA_integer_,
           MAF = NA_real_, r2_null = NA_real_, r2_full = NA_real_,
           r2_snp = NA_real_, PVE_pct = NA_real_,
           note = "SNP not in genotype matrix"))
  
  if (!trait_id %in% colnames(blues_mat)) return(
    tibble(SNP = snp_id, Trait = trait_id, n = NA_integer_,
           MAF = NA_real_, r2_null = NA_real_, r2_full = NA_real_,
           r2_snp = NA_real_, PVE_pct = NA_real_,
           note = "Trait not in STI file"))
  
  g  <- geno_mat[, snp_id]
  y  <- blues_mat[, trait_id]
  pc <- as.data.frame(pc_mat)
  
  keep <- !is.na(g) & !is.na(y) & complete.cases(pc)
  g  <- g[keep]; y <- y[keep]; pc <- pc[keep, , drop = FALSE]
  n  <- sum(keep)
  
  if (n < 20) return(
    tibble(SNP = snp_id, Trait = trait_id, n = n,
           MAF = NA_real_, r2_null = NA_real_, r2_full = NA_real_,
           r2_snp = NA_real_, PVE_pct = NA_real_,
           note = "insufficient data (n < 20)"))
  
  maf      <- min(mean(g) / 2, 1 - mean(g) / 2)
  fit_null <- lm(y ~ ., data = pc)
  fit_full <- lm(y ~ g + ., data = cbind(g = g, pc))
  r2_null  <- summary(fit_null)$r.squared
  r2_full  <- summary(fit_full)$r.squared
  r2_snp   <- max(0, r2_full - r2_null)
  
  tibble(
    SNP     = snp_id, Trait   = trait_id, n = n,
    MAF     = round(maf, 4),
    r2_null = round(r2_null, 4),
    r2_full = round(r2_full, 4),
    r2_snp  = round(r2_snp, 4),
    PVE_pct = round(r2_snp * 100, 2),
    note    = "ok"
  )
}

#############################################
# 6) Run
#############################################

message("Computing PC-adjusted PVE for ", nrow(snp_trait_pairs), " SNP x trait pairs...")

pve_results_adj <- map2_dfr(
  snp_trait_pairs$SNP,
  snp_trait_pairs$Trait,
  ~ compute_pve_adjusted(.x, .y, geno_s7, blues, pcs)
)

message("Done. Notes breakdown:")
print(table(pve_results_adj$note))

#############################################
# 7) Summarise per SNP
#############################################

pve_per_snp_adj <- pve_results_adj %>%
  filter(note == "ok") %>%
  group_by(SNP) %>%
  summarise(
    n_traits_computed = n(),
    mean_PVE_pct  = round(mean(PVE_pct), 2),
    max_PVE_pct   = round(max(PVE_pct), 2),
    min_PVE_pct   = round(min(PVE_pct), 2),
    mean_r2_null  = round(mean(r2_null), 2),
    best_trait    = Trait[which.max(PVE_pct)],
    .groups = "drop"
  )

#############################################
# 8) Join to S7 and write
#############################################

s7_final <- s7 %>%
  left_join(pve_per_snp_adj, by = "SNP")

write.csv(pve_results_adj,
          file.path(outdir, "TableS7_PVE_adjusted_per_SNP_trait.csv"),
          row.names = FALSE)

write.csv(s7_final,
          file.path(outdir, "TableS7_annotated_with_PVE_adjusted.csv"),
          row.names = FALSE)




