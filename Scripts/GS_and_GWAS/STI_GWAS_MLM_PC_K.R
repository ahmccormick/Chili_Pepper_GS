#############################################
# STI GWAS in Capsicum annuum using GAPIT
# Model: MLM (PC + K)
# Phenotype: pre-computed STI from BLUE-adjusted data
#############################################

rm(list = ls())

library(readxl)
library(dplyr)

source("http://zzlab.net/GAPIT/gapit_functions.txt")

#############################################
# 0) File paths -- EDIT THESE
#############################################

geno_file <- "~/R/World_Veg_Project/filtered_data/genotype_matrix_imputed.csv"
meta_file <- "~/R/World_Veg_Project/filtered_data/423.xlsx"
sti_file  <- "~/R/World_Veg_Project/anna_BLUES/anna_BLUES_rebuild/STI_from_BLUEs_all_73_traits.csv"

# OUTPUT FOLDER
out_base  <- "~/R/WorldVeg_Capsicum/STI_GWAS/MLM_PC_K/"

# Analysis settings
id_col    <- "accession"
pca_total <- 5
min_n     <- 30

#############################################
# 1) Load genotype matrix
#############################################

geno_df_imputed <- read.csv(geno_file, row.names = 1, check.names = FALSE)
rownames(geno_df_imputed) <- trimws(rownames(geno_df_imputed))
message("Genotype matrix loaded: ", nrow(geno_df_imputed), " samples x ", ncol(geno_df_imputed), " markers")

#############################################
# 2) Restrict to C. annuum
#############################################

metadata <- read_excel(meta_file) %>%
  rename(Sample_ID = `Sample ID`) %>%
  mutate(Sample_ID = trimws(as.character(Sample_ID)),
         Species   = trimws(as.character(Species)))

annuum_ids         <- metadata$Sample_ID[metadata$Species == "Capsicum annuum"]
annuum_ids_in_geno <- intersect(annuum_ids, rownames(geno_df_imputed))
geno_df_imputed_annuum <- geno_df_imputed[annuum_ids_in_geno, , drop = FALSE]
message("C. annuum genotype matrix: ", nrow(geno_df_imputed_annuum), " samples x ", ncol(geno_df_imputed_annuum), " markers")

#############################################
# 3) Build SNP map
#############################################

snp_names <- colnames(geno_df_imputed_annuum)
GM_annuum <- data.frame(
  SNP        = snp_names,
  Chromosome = suppressWarnings(as.numeric(sub(":.*", "", snp_names))),
  Position   = suppressWarnings(as.numeric(sub(".*:", "", snp_names))),
  stringsAsFactors = FALSE
)
GM_annuum <- GM_annuum[!is.na(GM_annuum$Chromosome) & !is.na(GM_annuum$Position), ]
geno_df_imputed_annuum <- geno_df_imputed_annuum[, GM_annuum$SNP, drop = FALSE]
annuum_geno_ids <- rownames(geno_df_imputed_annuum)
message("Markers retained after map parsing: ", ncol(geno_df_imputed_annuum))

#############################################
# 4) Load pre-computed STI matrix
#############################################

sti_wide <- read.csv(sti_file, stringsAsFactors = FALSE)
sti_wide[[id_col]] <- trimws(as.character(sti_wide[[id_col]]))

traits_common <- sort(setdiff(colnames(sti_wide), id_col))
message("Traits in STI file: ", length(traits_common))
print(traits_common)

#############################################
# 5) Helper: clean trait names for folders
#############################################

clean_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^[:alnum:]_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

#############################################
# 6) Main GWAS loop
#############################################

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
gwas_log <- list()

for (tr in traits_common) {
  
  message("\n============================================")
  message("Processing trait: ", tr)
  
  # extract STI values for this trait
  sti_df <- sti_wide %>%
    select(accession = all_of(id_col), STI = all_of(tr)) %>%
    filter(!is.na(STI)) %>%
    filter(accession %in% annuum_geno_ids)
  
  n_sti <- nrow(sti_df)
  message("  N after filter = ", n_sti)
  
  if (n_sti < min_n) {
    message("  Skipping trait (N < ", min_n, ")")
    gwas_log[[length(gwas_log) + 1]] <- data.frame(
      trait  = tr,
      N      = n_sti,
      status = "skipped_lowN",
      stringsAsFactors = FALSE
    )
    next
  }
  
  common_ids <- intersect(annuum_geno_ids, sti_df$accession)
  
  GD <- data.frame(Taxa = common_ids,
                   geno_df_imputed_annuum[common_ids, , drop = FALSE],
                   check.names = FALSE)
  
  trait_label <- paste0("STI_", clean_name(tr))
  Y <- data.frame(Taxa = sti_df$accession, STI = sti_df$STI, check.names = FALSE)
  colnames(Y)[2] <- trait_label
  Y <- Y[Y$Taxa %in% common_ids, ]
  
  # TRAIT-LEVEL OUTPUT SUBFOLDER (inside out_base above)
  trait_outdir <- file.path(out_base, clean_name(tr))
  dir.create(trait_outdir, recursive = TRUE, showWarnings = FALSE)
  
  old_wd <- getwd()
  setwd(trait_outdir)
  
  res <- tryCatch({
    GAPIT(
      Y         = Y,
      GD        = GD,
      GM        = GM_annuum,
      PCA.total = pca_total,
      model     = "MLM"
    )
    data.frame(trait = tr, N = n_sti, status = "ran", stringsAsFactors = FALSE)
  }, error = function(e) {
    message("    ERROR: ", e$message)
    data.frame(trait = tr, N = n_sti, status = paste0("error: ", e$message), stringsAsFactors = FALSE)
  })
  
  setwd(old_wd)
  gwas_log[[length(gwas_log) + 1]] <- res
}

#############################################
# 7) Save run log
#############################################

gwas_log_df <- bind_rows(gwas_log)

# LOG FILE -- saved inside out_base
write.csv(gwas_log_df,
          file = file.path(out_base, "STI_GWAS_MLM_BLUEadjusted_log.csv"),
          row.names = FALSE)
