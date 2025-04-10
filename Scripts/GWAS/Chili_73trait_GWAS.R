# Load genotype data 
geno_df_imputed <- read.csv(
  "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv",
  row.names = 1, check.names = FALSE
)

# Function to run GAPIT per environment and trait
run_gwas_for_environment <- function(env_name, pheno_file, outdir) {
  # Load phenotype data
  pheno_data <- read.csv(pheno_file, row.names = 1)
  
  # Loop over all traits in phenotype file
  for (trait in colnames(pheno_data)) {
    message("Running GWAS for ", trait, " in ", env_name)
    
    # Filter to samples in both geno and pheno
    common_ids <- intersect(rownames(geno_df_imputed), rownames(pheno_data))
    geno <- geno_df_imputed[common_ids, ]
    pheno <- pheno_data[common_ids, , drop = FALSE]  # Keep as data.frame
    
    # Format for GAPIT
    geno_gwas <- data.frame(Taxa = rownames(geno), geno, check.names = FALSE)
    rownames(geno_gwas) <- NULL
    pheno_gwas <- data.frame(Taxa = rownames(pheno), Trait = pheno[[trait]])
    colnames(pheno_gwas)[2] <- trait  # Name second column as trait
    
    # Create genetic map
    snp_names <- colnames(geno_gwas)[-1]
    GM <- data.frame(
      SNP = snp_names,
      Chromosome = as.numeric(sub(":.*", "", snp_names)),
      Position = as.numeric(sub(".*:", "", snp_names)),
      stringsAsFactors = FALSE
    )
    
    # Set working directory to save results
    trait_dir <- file.path(outdir, paste0("GAPIT_", env_name, "_", trait))
    dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
    setwd(trait_dir)
    
    # Load GAPIT functions if not already loaded
    source("http://zzlab.net/GAPIT/gapit_functions.txt")
    
    # Run GAPIT
    GAPIT(
      Y = pheno_gwas,
      GD = geno_gwas,
      GM = GM,
      PCA.total = 5,
      model = "FarmCPU"
    )
  }
}


# Run for Control
run_gwas_for_environment(
  env_name = "Control",
  pheno_file = "~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv",
  outdir = "~/R/World_Veg_Collab_Pepper/GAPIT_GWAS/GAPIT_GWAS_all_control/"
)

# Run for Heat Stress 1
run_gwas_for_environment(
  env_name = "Heat1",
  pheno_file = "~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv",
  outdir = "~/R/World_Veg_Collab_Pepper/GAPIT_GWAS/GAPIT_GWAS_all_heat1/"
)

# Run for Heat Stress 2
run_gwas_for_environment(
  env_name = "Heat2",
  pheno_file = "~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv",
  outdir = "~/R/World_Veg_Collab_Pepper/GAPIT_GWAS/GAPIT_GWAS_all_heat2/"
)
