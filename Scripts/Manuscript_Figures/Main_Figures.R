###########################################
# FIGURE 1
###########################################

###################################
# GWAS (Capsicum annuum only) using BLUEs per environment (GAPIT FarmCPU) 
# BLUEs were generated using the supplemental script Anna_BLUEs.R
###################################
library(readxl)
library(dplyr)

###################################
# Load genotype matrix
###################################
geno_df_imputed <- read.csv(
  "~/R/World_Veg_Project/filtered_data/genotype_matrix_imputed.csv",
  row.names = 1,
  check.names = FALSE
)

###################################
# Load metadata and filter to C. annuum
###################################

library(readxl)
library(dplyr)

metadata <- readxl::read_excel(
  "~/R/World_Veg_Project/filtered_data/423.xlsx"
)

metadata <- dplyr::rename(
  metadata,
  Sample_ID = `Sample ID`
)

metadata$Sample_ID <- trimws(as.character(metadata$Sample_ID))

annuum_ids <- metadata$Sample_ID[metadata$Species == "Capsicum annuum"]

# Subset genotype matrix to annuum-only
annuum_ids_in_geno <- intersect(annuum_ids, rownames(geno_df_imputed))
geno_df_imputed_annuum <- geno_df_imputed[annuum_ids_in_geno, , drop = FALSE]

message("Genotypes (all):   ", nrow(geno_df_imputed))
message("Genotypes (annuum):", nrow(geno_df_imputed_annuum), " x ", ncol(geno_df_imputed_annuum))

###################################
# Load BLUEs tables (Control / Heat1 / Heat2)
###################################
blues_ctrl <- read.csv(
  "~/R/World_Veg_Project/filtered_data/BLUEs/BLUEs_Control.csv",
  stringsAsFactors = FALSE
)
blues_h1 <- read.csv(
  "~/R/World_Veg_Project/filtered_data/BLUEs/BLUEs_Heat_Stress_1.csv",
  stringsAsFactors = FALSE
)
blues_h2 <- read.csv(
  "~/R/World_Veg_Project/filtered_data/BLUEs/BLUEs_Heat_Stress_2.csv",
  stringsAsFactors = FALSE
)

# Standardize accession IDs
blues_ctrl$accession <- trimws(as.character(blues_ctrl$accession))
blues_h1$accession   <- trimws(as.character(blues_h1$accession))
blues_h2$accession   <- trimws(as.character(blues_h2$accession))

###################################
# Quick overlap checks 
###################################
geno_ids_annuum <- rownames(geno_df_imputed_annuum)

message("Overlap Control (max): ", length(intersect(blues_ctrl$accession, geno_ids_annuum)))
message("Overlap Heat1   (max): ", length(intersect(blues_h1$accession,   geno_ids_annuum)))
message("Overlap Heat2   (max): ", length(intersect(blues_h2$accession,   geno_ids_annuum)))

###################################
# GAPIT GWAS runner for BLUEs
###################################
run_gwas_from_blues <- function(env_name, blues_df, outdir, geno_mat,
                                id_col = "accession",
                                drop_cols = c("plantno."),   # columns to exclude from GWAS if present
                                pca_total = 5,
                                model = "FarmCPU",
                                min_n = 30) {
  
  stopifnot(id_col %in% colnames(blues_df))
  
  # Ensure ID column is clean + set rownames
  blues_df[[id_col]] <- trimws(as.character(blues_df[[id_col]]))
  rownames(blues_df) <- blues_df[[id_col]]
  
  # Trait columns = all except ID and any dropped columns
  trait_cols <- setdiff(colnames(blues_df), c(id_col, drop_cols))
  
  # Load GAPIT once per environment run
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  
  for (trait in trait_cols) {
    message("Running GWAS for ", trait, " in ", env_name)
    
    # Pull phenotype vector for this trait
    y <- blues_df[, trait, drop = FALSE]
    
    # Coerce to numeric if needed (CSV sometimes loads as character)
    if (!is.numeric(y[[trait]])) {
      suppressWarnings(y[[trait]] <- as.numeric(y[[trait]]))
    }
    
    # Drop missing phenotypes for this trait
    y <- y[!is.na(y[[trait]]), , drop = FALSE]
    
    # Intersect genotype and phenotype IDs
    common_ids <- intersect(rownames(geno_mat), rownames(y))
    message("  N common (geno + BLUE) = ", length(common_ids))
    
    if (length(common_ids) < min_n) {
      message("  Skipping ", trait, " (N < ", min_n, " after filtering).")
      next
    }
    
    geno  <- geno_mat[common_ids, , drop = FALSE]
    pheno <- y[common_ids, , drop = FALSE]
    
    # GAPIT input formats
    geno_gwas <- data.frame(Taxa = rownames(geno), geno, check.names = FALSE)
    rownames(geno_gwas) <- NULL
    
    pheno_gwas <- data.frame(Taxa = rownames(pheno), Trait = pheno[[trait]])
    colnames(pheno_gwas)[2] <- trait
    
    # Genetic map from SNP names like "Chr:Pos"
    snp_names <- colnames(geno_gwas)[-1]
    GM <- data.frame(
      SNP        = snp_names,
      Chromosome = as.numeric(sub(":.*", "", snp_names)),
      Position   = as.numeric(sub(".*:", "", snp_names)),
      stringsAsFactors = FALSE
    )
    
    # Output directory per trait
    trait_dir <- file.path(outdir, paste0("GAPIT_", env_name, "_", trait))
    dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
    setwd(trait_dir)
    
    # Run GAPIT
    GAPIT(
      Y = pheno_gwas,
      GD = geno_gwas,
      GM = GM,
      PCA.total = pca_total,
      model = model
    )
  }
}

get_N_one <- function(blues_df, trait, geno_mat = geno_df_imputed_annuum) {
  y <- blues_df[, c("accession", trait)]
  y$accession <- trimws(as.character(y$accession))
  y <- y[!is.na(y[[trait]]), ]
  length(intersect(rownames(geno_mat), y$accession))
}

get_N_one(blues_ctrl, "yield")
get_N_one(blues_h1,  "yield")
get_N_one(blues_h2,  "yield")


###################################
# Run GWAS for each environment 
###################################

# Control BLUEs
run_gwas_from_blues(
  env_name = "Control_annuum_BLUE",
  blues_df = blues_ctrl,
  outdir   = "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_control/",
  geno_mat = geno_df_imputed_annuum
)

# Heat Stress 1 BLUEs
run_gwas_from_blues(
  env_name = "Heat1_annuum_BLUE",
  blues_df = blues_h1,
  outdir   = "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_heat1/",
  geno_mat = geno_df_imputed_annuum
)

# Heat Stress 2 BLUEs
run_gwas_from_blues(
  env_name = "Heat2_annuum_BLUE",
  blues_df = blues_h2,
  outdir   = "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_heat2/",
  geno_mat = geno_df_imputed_annuum
)


#############################################
# BLUE-based STI for ALL traits + GWAS
#############################################

library(dplyr)

# 1) Identify common trait columns across the three BLUE tables
id_col <- "accession"
drop_cols <- c("plantno.")  # non-trait columns to exclude if present

traits_common <- Reduce(intersect, list(
  setdiff(colnames(blues_ctrl), c(id_col, drop_cols)),
  setdiff(colnames(blues_h1),   c(id_col, drop_cols)),
  setdiff(colnames(blues_h2),   c(id_col, drop_cols))
))

# Optional: keep only numeric-ish traits (can still be read as character)
# We'll coerce to numeric inside the loop anyway.
message("Common traits across environments: ", length(traits_common))
print(traits_common)

# 2) Function to compute STI from BLUEs for one trait
make_sti_from_blues_one_trait <- function(trait) {
  
  dC  <- blues_ctrl %>% dplyr::select(dplyr::all_of(id_col), C  = dplyr::all_of(trait))
  dS1 <- blues_h1   %>% dplyr::select(dplyr::all_of(id_col), S1 = dplyr::all_of(trait))
  dS2 <- blues_h2   %>% dplyr::select(dplyr::all_of(id_col), S2 = dplyr::all_of(trait))
  
  d <- dC %>%
    dplyr::inner_join(dS1, by = id_col) %>%
    dplyr::inner_join(dS2, by = id_col) %>%
    dplyr::mutate(
      C  = suppressWarnings(as.numeric(C)),
      S1 = suppressWarnings(as.numeric(S1)),
      S2 = suppressWarnings(as.numeric(S2))
    ) %>%
    dplyr::filter(!is.na(C), !is.na(S1), !is.na(S2))
  
  Xc_bar <- mean(d$C, na.rm = TRUE)
  
  d <- d %>%
    dplyr::mutate(
      STI = sqrt(S1 * S2 * C) / (Xc_bar^2)
    ) %>%
    dplyr::filter(is.finite(STI)) %>%
    dplyr::select(accession, STI)
  
  d
}


# 3) Loop over traits -> compute STI -> GWAS
#    Output folders: one subfolder per trait
out_base <- "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_STI_alltraits/"

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Annnum-only genotype IDs (optional but recommended)
annuum_geno_ids <- rownames(geno_df_imputed_annuum)

# keep a log of sample sizes and whether trait ran
sti_log <- list()

for (tr in traits_common) {
  message("\n==============================")
  message("STI + GWAS for trait: ", tr)
  
  sti_df <- make_sti_from_blues_one_trait(tr)
  
  # restrict to annuum genotypes (so GWAS is annuum-only)
  sti_df <- sti_df %>% filter(accession %in% annuum_geno_ids)
  
  n_sti <- nrow(sti_df)
  message("  N STI (after merges + annuum filter) = ", n_sti)
  
  # Skip if too few samples
  if (n_sti < 30) {
    message("  Skipping ", tr, " (N < 30).")
    sti_log[[tr]] <- data.frame(trait = tr, N = n_sti, status = "skipped_lowN")
    next
  }
  
  # Run GWAS on STI (treat STI as the phenotype trait)
  # run_gwas_from_blues expects a df with accession + trait columns
  sti_for_gwas <- sti_df %>% rename(!!paste0("STI_", tr) := STI)
  
  trait_outdir <- file.path(out_base, tr)
  dir.create(trait_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # If something errors for a trait (e.g., all NA after numeric coercion), keep going
  tryCatch({
    run_gwas_from_blues(
      env_name = paste0("STI_annuum_BLUE_", tr),
      blues_df = sti_for_gwas,
      outdir   = trait_outdir,
      geno_mat = geno_df_imputed_annuum,
      drop_cols = character(0)
    )
    sti_log[[tr]] <- data.frame(trait = tr, N = n_sti, status = "ran")
  }, error = function(e) {
    message("  ERROR for trait ", tr, ": ", e$message)
    sti_log[[tr]] <- data.frame(trait = tr, N = n_sti, status = paste0("error: ", e$message))
  })
}

###########################################
# FIGURE: Circular GWAS plots (annuum BLUEs)
# Control vs Heat1 vs Heat2 vs STI (Yield)
###########################################

library(dplyr)
library(purrr)
library(CMplot)

trait <- "yield"

# --- Existing per-environment trait dirs (as you had) ---
ctrl_trait_dir <- "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_control/GAPIT_Control_annuum_BLUE_yield"
h1_trait_dir   <- "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_heat1/GAPIT_Heat1_annuum_BLUE_yield"
h2_trait_dir   <- "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_heat2/GAPIT_Heat2_annuum_BLUE_yield"

# --- NEW: STI yield trait dir ---
# Your STI pipeline used: out_base/<trait>/  and env_name = "STI_annuum_BLUE_<trait>"
# so the trait directory typically becomes:
#   out_base/yield/GAPIT_STI_annuum_BLUE_yield
sti_trait_dir  <- "~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_STI_alltraits/yield/GAPIT_STI_annuum_BLUE_yield_STI_yield"


find_gapit_gwas_file <- function(trait_dir, trait = "yield") {
  # Look for the FarmCPU GWAS results CSV produced by GAPIT
  # Try common patterns first
  candidates <- c(
    list.files(trait_dir, pattern = "GAPIT\\.Association\\.GWAS_Results\\.FarmCPU.*\\.csv$", full.names = TRUE),
    list.files(trait_dir, pattern = "GAPIT\\.Association\\.GWAS_Results.*\\.csv$", full.names = TRUE),
    list.files(trait_dir, pattern = "GAPIT\\.GWAS\\.Results.*\\.csv$", full.names = TRUE)
  ) %>% unlist() %>% unique()
  
  if (length(candidates) == 0) {
    stop("No GAPIT GWAS results CSV found in: ", trait_dir)
  }
  
  # Prefer files that mention the trait name (case-insensitive)
  trait_hits <- candidates[grepl(trait, basename(candidates), ignore.case = TRUE)]
  if (length(trait_hits) > 0) return(trait_hits[1])
  
  # Otherwise just take the first candidate
  candidates[1]
}


# Find GAPIT GWAS results files
ctrl_file <- find_gapit_gwas_file(ctrl_trait_dir, trait = trait)
h1_file   <- find_gapit_gwas_file(h1_trait_dir,   trait = trait)
h2_file   <- find_gapit_gwas_file(h2_trait_dir,   trait = trait)
sti_file  <- find_gapit_gwas_file(sti_trait_dir,  trait = trait)  # <- NEW

message("Control GWAS file: ", ctrl_file)
message("Heat1   GWAS file: ", h1_file)
message("Heat2   GWAS file: ", h2_file)
message("STI     GWAS file: ", sti_file)

clean_gapit_for_cmplot <- function(path, p_colname) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  
  # Standardize possible column names
  # GAPIT outputs often include: SNP, Chr, Pos, P.value
  # Sometimes: Chromosome, Position, P.value
  col_map <- list(
    SNP        = c("SNP", "snp", "Marker", "marker"),
    Chr        = c("Chr", "chr", "Chromosome", "chromosome"),
    Pos        = c("Pos", "pos", "Position", "position"),
    P.value    = c("P.value", "P.Value", "p.value", "P", "p", "Pvalue", "pvalue")
  )
  
  pick_col <- function(df, options) {
    hit <- options[options %in% colnames(df)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  c_snp <- pick_col(df, col_map$SNP)
  c_chr <- pick_col(df, col_map$Chr)
  c_pos <- pick_col(df, col_map$Pos)
  c_p   <- pick_col(df, col_map$P.value)
  
  if (any(is.na(c(c_snp, c_chr, c_pos, c_p)))) {
    stop(
      "Missing required columns in: ", path, "\n",
      "Found columns: ", paste(colnames(df), collapse = ", "), "\n",
      "Need something like SNP/Chr/Pos/P.value (or Chromosome/Position)."
    )
  }
  
  out <- df %>%
    transmute(
      SNP        = .data[[c_snp]],
      Chromosome = readr::parse_number(as.character(.data[[c_chr]])),
      Position   = as.numeric(.data[[c_pos]]),
      P.value    = as.numeric(.data[[c_p]])
    ) %>%
    mutate(
      P.value = ifelse(is.na(P.value) | P.value <= 0, 1e-300, P.value)
    ) %>%
    transmute(SNP, Chromosome, Position, !!p_colname := P.value)
  
  out
}

# Read + standardize for CMplot
control <- clean_gapit_for_cmplot(ctrl_file, "P_control")
heat1   <- clean_gapit_for_cmplot(h1_file,   "P_heat1")
heat2   <- clean_gapit_for_cmplot(h2_file,   "P_heat2")
sti     <- clean_gapit_for_cmplot(sti_file,  "P_sti")   # <- NEW

# Merge all (full joins keep SNPs that exist in any run)
gwas_wide4 <- purrr::reduce(
  list(control, heat1, heat2, sti),
  dplyr::full_join,
  by = c("SNP","Chromosome","Position")
) %>%
  mutate(Chromosome = as.numeric(Chromosome))

# IMPORTANT:
# To make STI the OUTER ring, put P_sti first (CMplot uses p-columns order for rings)
gwas_wide4 <- purrr::reduce(
  list(control, heat1, heat2, sti),
  dplyr::full_join,
  by = c("SNP","Chromosome","Position")
) %>%
  dplyr::mutate(Chromosome = as.numeric(Chromosome)) %>%
  dplyr::select(SNP, Chromosome, Position, P_sti, P_control, P_heat1, P_heat2)


# Circos plot (4 rings: STI outer -> Control -> Heat1 -> Heat2 inner)
CMplot(
  gwas_wide4,
  type = "p",
  plot.type = "c",
  r = 0.30,
  col = c("grey30", "grey60"),
  chr.labels = paste("Chr", 1:12),
  threshold = c(1e-6, 1e-4),
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  signal.line = 1,
  signal.col = c("red", "blue"),
  chr.den.col = c("blue", "white", "red"),
  bin.size = 1e6,
  outward = FALSE,
  file = "pdf",
  file.name = "GWAS_CompositeRings_AnnuumBLUE_STI_Control_Heat1_Heat2_Yield",
  dpi = 300,
  file.output = TRUE,
  width = 15,
  height = 15
)
#Plots are stored in: /Users/annamccormick/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_STI_alltraits/lightpenetrationfinal/GAPIT_STI_annuum_BLUE_lightpenetrationfinal_STI_lightpenetrationfinal 

##linear plot
library(CMplot)

CMplot(
  gwas_wide4,
  type = "p",
  plot.type = "m",
  multracks = TRUE,
  LOG10 = TRUE,
  cex = 0.4,          # base point size
  signal.cex = 0.5,   # significant SNP size (default is ~1.2–1.5)
  threshold = c(1e-6, 1e-4),
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  file = "pdf",
  file.name = "Manhattan_Multitrack_STI_Control_Heat1_Heat2_Yield_dots",
  width = 20,
  height = 6
)

##table of significant SNPs per condition
library(dplyr)
library(tidyr)

thr_strict <- 1e-6

sig_all_strict <- gwas_wide4 %>%
  tidyr::pivot_longer(
    cols = c(P_sti, P_control, P_heat1, P_heat2),
    names_to = "Environment",
    values_to = "P"
  ) %>%
  dplyr::filter(!is.na(P), P < thr_strict) %>%
  dplyr::mutate(
    Environment = dplyr::recode(
      Environment,
      P_sti     = "STI",
      P_control = "Control",
      P_heat1   = "Heat1",
      P_heat2   = "Heat2"
    ),
    negLog10P = -log10(P)
  ) %>%
  dplyr::arrange(Chromosome, Position, P)

sig_all_strict
write.csv(sig_all_strict, "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/Significant_SNPs_strict_all_conditions.csv", row.names = FALSE)


##################################################################################################
# Table S5 annotations
##################################################################################################
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readr)

gff_path <- "/Users/annamccormick/R/World_Veg_Project/Reference_sequence_cornell/Zhangshugang.gff3"

gff <- import(gff_path)

# What chromosomes/seqlevels does the GFF use?
seqlevels(gff) %>% head()

# Usually you want genes as the main annotation target
genes <- gff[gff$type == "gene"]

# Sometimes gene IDs are stored in different columns; inspect:
mcols(genes) %>% names()

# Common fields: ID, Name, gene, gene_name, Parent, product, etc.
head(mcols(genes)[, intersect(names(mcols(genes)), c("ID","Name","gene","gene_name","product","locus_tag","Parent"))])


read_gapit_hits <- function(gapit_csv, p_thresh = 1e-6) {
  df <- read.csv(gapit_csv, stringsAsFactors = FALSE)
  
  # Defensive column mapping
  cn <- colnames(df)
  SNPcol <- intersect(c("SNP","Marker","snp"), cn)[1]
  Chrcol <- intersect(c("Chr","Chromosome","chr"), cn)[1]
  Poscol <- intersect(c("Pos","Position","pos"), cn)[1]
  Pcol   <- intersect(c("P.value","P.Value","p.value","P","p"), cn)[1]
  
  stopifnot(!is.na(SNPcol), !is.na(Chrcol), !is.na(Poscol), !is.na(Pcol))
  
  df %>%
    transmute(
      SNP = .data[[SNPcol]],
      Chr = as.character(.data[[Chrcol]]),
      Pos = as.integer(.data[[Poscol]]),
      P   = as.numeric(.data[[Pcol]])
    ) %>%
    filter(!is.na(P), P < p_thresh)
}


gff_chrs <- seqlevels(genes)
head(gff_chrs)


to_gff_chr <- function(x) {
  x <- trimws(as.character(x))
  
  # strip any leading "chr"/"Chr"/"CHR"
  x2 <- gsub("^(?i)chr", "", x, perl = TRUE)
  
  # try to interpret as integer chromosome number
  num <- suppressWarnings(as.integer(x2))
  
  # if numeric, convert to Chr01..Chr12; otherwise keep as-is
  out <- ifelse(!is.na(num), sprintf("Chr%02d", num), x)
  
  out
}

hits_to_gr <- function(hits_df) {
  hits_df <- hits_df %>%
    mutate(Chr = to_gff_chr(Chr))
  
  GRanges(
    seqnames = hits_df$Chr,
    ranges   = IRanges(start = hits_df$Pos, end = hits_df$Pos),
    SNP      = hits_df$SNP,
    P        = hits_df$P
  )
}


annotate_snps_with_genes <- function(snp_gr, genes_gr) {
  # Keep only SNPs that are on chromosomes present in genes
  snp_gr <- keepSeqlevels(snp_gr, intersect(seqlevels(snp_gr), seqlevels(genes_gr)), pruning.mode="coarse")
  
  # Overlaps: SNP falls inside a gene
  ov <- findOverlaps(snp_gr, genes_gr, ignore.strand = TRUE)
  
  in_gene <- tibble(
    SNP = mcols(snp_gr)$SNP[queryHits(ov)],
    P   = mcols(snp_gr)$P[queryHits(ov)],
    Chr = as.character(seqnames(snp_gr)[queryHits(ov)]),
    Pos = start(snp_gr)[queryHits(ov)],
    feature = "gene",
    gene_id   = as.character(mcols(genes_gr)$ID[subjectHits(ov)]),
    gene_name = as.character(if ("Name" %in% names(mcols(genes_gr))) mcols(genes_gr)$Name[subjectHits(ov)] else NA),
    gene_start = start(genes_gr)[subjectHits(ov)],
    gene_end   = end(genes_gr)[subjectHits(ov)]
  )
  
  # Nearest gene for every SNP (including in-gene; distance 0 if inside)
  nearest_idx <- nearest(snp_gr, genes_gr, ignore.strand = TRUE)
  nearest_tbl <- tibble(
    SNP = mcols(snp_gr)$SNP,
    P   = mcols(snp_gr)$P,
    Chr = as.character(seqnames(snp_gr)),
    Pos = start(snp_gr),
    nearest_gene_id   = as.character(mcols(genes_gr)$ID[nearest_idx]),
    nearest_gene_name = as.character(if ("Name" %in% names(mcols(genes_gr))) mcols(genes_gr)$Name[nearest_idx] else NA),
    nearest_gene_start = start(genes_gr)[nearest_idx],
    nearest_gene_end   = end(genes_gr)[nearest_idx],
    dist_to_nearest_bp = distanceToNearest(snp_gr, genes_gr[nearest_idx], ignore.strand=TRUE)@elementMetadata$distance
  )
  
  list(in_gene = in_gene, nearest = nearest_tbl)
}

p_strict <- 1e-6

ctrl_hits <- read_gapit_hits("~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_control/GAPIT_Control_annuum_BLUE_yield/GAPIT.Association.GWAS_Results.FarmCPU.yield(NYC).csv", p_strict)
h1_hits   <- read_gapit_hits("~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_heat1/GAPIT_Heat1_annuum_BLUE_yield/GAPIT.Association.GWAS_Results.FarmCPU.yield(NYC).csv", p_strict)
h2_hits   <- read_gapit_hits("~/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_heat2/GAPIT_Heat2_annuum_BLUE_yield/GAPIT.Association.GWAS_Results.FarmCPU.yield(NYC).csv", p_strict)



ctrl_gr <- hits_to_gr(ctrl_hits)
h1_gr   <- hits_to_gr(h1_hits)
h2_gr   <- hits_to_gr(h2_hits)

ctrl_ann <- annotate_snps_with_genes(ctrl_gr, genes)
h1_ann   <- annotate_snps_with_genes(h1_gr, genes)
h2_ann   <- annotate_snps_with_genes(h2_gr, genes)

write.csv(ctrl_ann$in_gene,  "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/yield_control_sig_snps_in_genes.csv", row.names=FALSE)
write.csv(ctrl_ann$nearest,  "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/yield_control_sig_snps_nearest_gene.csv", row.names=FALSE)

write.csv(h1_ann$in_gene,    "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/yield_heat1_sig_snps_in_genes.csv", row.names=FALSE)
write.csv(h1_ann$nearest,    "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/yield_heat1_sig_snps_nearest_gene.csv", row.names=FALSE)

write.csv(h2_ann$in_gene,    "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/yield_heat2_sig_snps_in_genes.csv", row.names=FALSE)
write.csv(h2_ann$nearest,    "~/R/World_Veg_Project/filtered_data/GWAS_HITS_ANNOTATIONS/yield_heat2_sig_snps_nearest_gene.csv", row.names=FALSE)


###################################
#Supplemental Tables S6 and S7
###################################
library(tidyverse)
library(stringr)

out_base <- "/Users/annamccormick/R/World_Veg_Project/filtered_data/GAPIT_GWAS_annuum_BLUE_STI_alltraits/"
exclude_traits <- c("heightmax", "leafanglecoldminushot", "leafangledaily", "leafminushobotemp")
logp_thresh <- 6
p_thresh <- 10^(-logp_thresh)  # 1e-6


# 1) Find ALL CSVs recursively (because outputs are inside trait folders)
all_csv <- list.files(out_base, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)



# 2) Helper: get trait from path (first folder under out_base)
trait_from_path <- function(path, out_base) {
  rel <- sub(paste0("^", normalizePath(out_base), "/?"), "", normalizePath(path))
  strsplit(rel, "/")[[1]][1]
}


# 3) Helper: detect p-value column
detect_p_col <- function(df) {
  p_col <- grep("P\\.value$|P_value$|P\\.value\\.$|P\\.Value$|P\\.value|P_value", names(df),
                ignore.case = TRUE, value = TRUE)[1]
  if (is.na(p_col)) return(NA_character_)
  p_col
}


# 4) Make a simple trait->one-file table (pick ONE CSV per trait)
#    This mimics your old "one CSV per trait" workflow.
#    Preference: a CSV containing "Association" or "GWAS_Results" in the filename,
#    else any CSV with "GAPIT", else the first one.

file_tbl <- tibble(path = all_csv) %>%
  mutate(Trait = map_chr(path, ~trait_from_path(.x, out_base))) %>%
  filter(!(tolower(Trait) %in% tolower(exclude_traits))) %>%
  group_by(Trait) %>%
  summarise(
    gwas_csv = {
      files <- path
      assoc <- files[grepl("Association|GWAS_Results|GWAS.Results", basename(files), ignore.case = TRUE)]
      gapit <- files[grepl("GAPIT", basename(files), ignore.case = TRUE)]
      if (length(assoc) > 0) assoc[1]
      else if (length(gapit) > 0) gapit[1]
      else files[1]
    },
    .groups = "drop"
  )

cat("Traits detected (after exclusion):", nrow(file_tbl), "\n")

###################################
# TABLE S6: count significant SNPs per trait at -log10(p) >= 6
###################################
count_sig_snps <- function(file, trait) {
  df <- suppressWarnings(read.csv(file))
  p_col <- detect_p_col(df)
  if (is.na(p_col) || !("SNP" %in% names(df))) {
    return(tibble(Trait = trait, Significant_SNPs = NA_integer_))
  }
  
  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  sig_n <- sum(pvals <= p_thresh, na.rm = TRUE)
  
  tibble(Trait = trait, Significant_SNPs = sig_n)
}

table_s6_all <- pmap_dfr(
  list(file_tbl$gwas_csv, file_tbl$Trait),
  ~count_sig_snps(..1, ..2)
) %>%
  arrange(desc(Significant_SNPs), Trait)

# “publication” version: only traits with >=1 sig SNP
table_s6 <- table_s6_all %>%
  filter(!is.na(Significant_SNPs), Significant_SNPs > 0)

cat("Traits with >=1 significant SNP:", nrow(table_s6), "/", nrow(file_tbl), "\n")

write.csv(table_s6,
          file.path(out_base, "STI_BLUE_Table_S6_sigSNPcounts_perTrait_logp6.csv"),
          row.names = FALSE)

# optional QA: include zeros too
write.csv(table_s6_all,
          file.path(out_base, "STI_BLUE_Table_S6_ALLtraits_includingZeros.csv"),
          row.names = FALSE)

###################################
# TABLE S7: SNPs significant in >=2 traits (shared)
###################################
extract_sig_snps <- function(file, trait) {
  df <- suppressWarnings(read.csv(file))
  p_col <- detect_p_col(df)
  if (is.na(p_col) || !("SNP" %in% names(df))) return(tibble())
  
  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  
  tibble(SNP = as.character(df$SNP), pval = pvals) %>%
    filter(!is.na(SNP), !is.na(pval), pval <= p_thresh) %>%
    distinct(SNP) %>%
    mutate(Trait = trait)
}

sig_long <- pmap_dfr(
  list(file_tbl$gwas_csv, file_tbl$Trait),
  ~extract_sig_snps(..1, ..2)
)

table_s7 <- sig_long %>%
  group_by(SNP) %>%
  summarise(
    Traits_with_SNP = n_distinct(Trait),
    Traits = paste(sort(unique(Trait)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(Traits_with_SNP >= 2) %>%
  arrange(desc(Traits_with_SNP), SNP)

write.csv(table_s7,
          file.path(out_base, "STI_BLUE_Table_S7_sharedSigSNPs_acrossTraits_logp6.csv"),
          row.names = FALSE)


###########################################
# FIGURE 2A
###########################################
library(tidyverse)
library(ggpubr)
library(janitor)

# Load and clean
df <- readxl::read_excel("/Users/annamccormick/R/World_Veg_Project/data/inputs/Phenospec_climatedata_timepoints3X.xlsx") %>%
  clean_names()

# Make sure day is a Date
df$day <- as.Date(df$day)

# Define correct trait names
traits_of_interest <- c("temperature_c",
                        "relative_humidity_percent",
                        "solar_radiation_daily_mj_m2",
                        "precipitation_mm",
                        "wind_speed_m_s",
                        "wind_direction")

# Convert traits to numeric
df <- df %>%
  mutate(across(all_of(traits_of_interest), as.numeric))

# Pivot to long format
df_long <- df %>%
  pivot_longer(cols = all_of(traits_of_interest),
               names_to = "trait",
               values_to = "value")

# Plot all traits by timepoint
ggboxplot(df_long, x = "timepoint", y = "value",
          fill = "timepoint", add = "jitter") +
  facet_wrap(~ trait, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons = list(
                       c("Control", "HeatStress1"),
                       c("Control", "HeatStress2"),
                       c("HeatStress1", "HeatStress2")),
                     p.adjust.method = "bonferroni") +
  labs(title = "Climate Traits Across Timepoints",
       y = "Value",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1.5, "lines"),       # ← Adjust this value as needed
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## temp and solar rad
library(tidyverse)
library(ggpubr)
library(janitor)

# Load and clean
df <- readxl::read_excel("/Users/annamccormick/R/World_Veg_Project/data/inputs/Phenospec_climatedata_timepoints3X.xlsx") %>%
  clean_names()

# Make sure day is a Date
df$day <- as.Date(df$day)

# Traits of interest
traits_of_interest <- c("temperature_c", "solar_radiation_daily_mj_m2")

# Convert traits to numeric
df <- df %>%
  mutate(across(all_of(traits_of_interest), as.numeric))

# Pivot to long format (only those 2 traits)
df_long <- df %>%
  pivot_longer(cols = all_of(traits_of_interest),
               names_to = "trait",
               values_to = "value")

# Plot
ggboxplot(df_long, x = "timepoint", y = "value",
          fill = "timepoint", add = "jitter") +
  facet_wrap(~ trait, scales = "free_y", ncol = 2) +
  stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons = list(
                       c("Control", "HeatStress1"),
                       c("Control", "HeatStress2"),
                       c("HeatStress1", "HeatStress2")),
                     p.adjust.method = "bonferroni") +
  labs(title = "",
       y = "Value",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Lable clean-up

# Nice labels for traits
trait_labels <- c(
  "temperature_c" = "Temperature (°C)",
  "solar_radiation_daily_mj_m2" = "Solar radiation (MJ/m²)"
)

# Plot
ggboxplot(df_long, x = "timepoint", y = "value",
          fill = "timepoint", add = "jitter") +
  facet_wrap(~ trait, scales = "free_y", ncol = 2,
             labeller = labeller(trait = trait_labels)) +
  stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons = list(
                       c("Control", "HeatStress1"),
                       c("Control", "HeatStress2"),
                       c("HeatStress1", "HeatStress2")),
                     p.adjust.method = "bonferroni") +
  labs(title = "",
       y = "Value",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggboxplot(df_long, x = "timepoint", y = "value",
          fill = "timepoint", add = "jitter") +
  facet_wrap(~ trait, scales = "free_y", ncol = 2,
             labeller = labeller(trait = trait_labels)) +
  stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons = list(
                       c("Control", "HeatStress1"),
                       c("Control", "HeatStress2"),
                       c("HeatStress1", "HeatStress2")),
                     p.adjust.method = "bonferroni") +
  labs(title = "",
       y = "Value",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c(
    "Control"     = "#1f77b4",
    "HeatStress1" = "#ff7f0e",
    "HeatStress2" = "red"
  ))

#fix plot labels
# Custom labels
time_labels <- c(
  "Control"     = "Control",
  "HeatStress1" = "Heat stress-1",
  "HeatStress2" = "Heat stress-2"
)

# Relabel factor levels in the data
df_long <- df_long %>%
  mutate(timepoint = recode(timepoint, !!!time_labels))

# Plot
ggboxplot(df_long, x = "timepoint", y = "value",
          fill = "timepoint", add = "jitter") +
  facet_wrap(~ trait, scales = "free_y", ncol = 2,
             labeller = labeller(trait = trait_labels)) +
  stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons = list(
                       c("Control", "Heat stress-1"),
                       c("Control", "Heat stress-2"),
                       c("Heat stress-1", "Heat stress-2")),
                     p.adjust.method = "bonferroni") +
  labs(title = "",
       y = "Value",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c(
    "Control"        = "#1f77b4",
    "Heat stress-1"  = "#ff7f0e",
    "Heat stress-2"  = "red"
  ))


###########################################
# FIGURE 2B
###########################################
################################################################################################
# Supplemental Figure 17
################################################################################################
# Core
rrblup_core_control <- readRDS("~/R/World_Veg_Project/data/outputs/core423/xval_rrblup_kfold_10.RData")$xval.result %>%
  mutate(model = "rrBLUP", set = "Core", condition = "Control")

rrblup_core_heat1 <- readRDS("~/R/World_Veg_Project/data/outputs/core423/xval_rrblup_kfold_10_heat1_model_1.RData")$xval.result %>%
  mutate(model = "rrBLUP", set = "Core", condition = "Heat1")

rrblup_core_heat2 <- readRDS("~/R/World_Veg_Project/data/outputs/core423/xval_rrblup_kfold_10_heat2_model_1.RData")$xval.result %>%
  mutate(model = "rrBLUP", set = "Core", condition = "Heat2")



# Global datasets
rrblup_global_control <- readRDS("~/R/World_Veg_Project/data/outputs/10250_GS/xval_rrblup_kfold_control_10k.RDS")$xval.result %>%
  mutate(model = "rrBLUP", set = "Global", condition = "Control")

rrblup_global_heat1 <- readRDS("~/R/World_Veg_Project/data/outputs/10250_GS/xval_rrblup_kfold_heat1_10k.RDS")$xval.result %>%
  mutate(model = "rrBLUP", set = "Global", condition = "Heat1")

rrblup_global_heat2 <- readRDS("~/R/World_Veg_Project/data/outputs/10250_GS/xval_rrblup_kfold_heat2_10k.RDS")$xval.result %>%
  mutate(model = "rrBLUP", set = "Global", condition = "Heat2")


library(dplyr)

# Combine all six result data frames into one
all_pa <- bind_rows(
  rrblup_core_control,
  rrblup_core_heat1,
  rrblup_core_heat2,
  rrblup_global_control,
  rrblup_global_heat1,
  rrblup_global_heat2
)

# Filter for traits with r.mean > 0.5
high_pa <- all_pa %>% filter(r.mean > 0.5)

# Count how many conditions each trait passed per set
trait_condition_count <- high_pa %>%
  group_by(set, trait) %>%
  summarise(n_conditions = n_distinct(condition), .groups = "drop") %>%
  filter(n_conditions == 3)

# Extract traits that passed in all 3 conditions for both sets
core_traits <- trait_condition_count %>% filter(set == "Core") %>% pull(trait)
global_traits <- trait_condition_count %>% filter(set == "Global") %>% pull(trait)

common_traits <- intersect(core_traits, global_traits)
length(common_traits)
print(common_traits)

#####################
#plotting the 13 common traits

common_traits <- c(
  "biomassfinalplant", "biomassfinalplot", "fruitlength", "fruitno", 
  "fruitshapeindex", "fruitweight", "fruitwidth", "leafangledaily", 
  "leafareadaily", "leafareadailyplant", "leafareaprojecteddaily", 
  "leafareaprojectedplant", "yield"
)

common_traits_prefixed <- paste0("GEBV_", common_traits)


library(tidyverse)

# Load Global (10k) GEBVs
gebv_global_control <- read_csv("~/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_all_73_Control_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Control", Group = "Global")

gebv_global_heat1 <- read_csv("~/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_all_73_Heat_Stress_1_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Heat Stress 1", Group = "Global")

gebv_global_heat2 <- read_csv("~/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_all_73_Heat_Stress_2_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Heat Stress 2", Group = "Global")

# Load Core (423) GEBVs
gebv_core_control <- read_csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_all_n423_73trait_control.csv") %>%
  mutate(Condition = "Control", Group = "Core")

gebv_core_heat1 <- read_csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_all_n423_73trait_Heatstress_1.csv") %>%
  mutate(Condition = "Heat Stress 1", Group = "Core")

gebv_core_heat2 <- read_csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_all_n423_73trait_Heatstress_2.csv") %>%
  mutate(Condition = "Heat Stress 2", Group = "Core")

# Combine all datasets
gebv_all <- bind_rows(
  gebv_global_control,
  gebv_global_heat1,
  gebv_global_heat2,
  gebv_core_control,
  gebv_core_heat1,
  gebv_core_heat2
)

# Pivot to long format
gebv_all_long <- gebv_all %>%
  pivot_longer(cols = all_of(common_traits_prefixed), names_to = "Trait", values_to = "GEBV")


# Filter for traits of interest (already narrowed)
gebv_filtered <- gebv_all_long %>%
  mutate(Trait = str_remove(Trait, "GEBV_"))

ggplot(gebv_filtered, aes(x = GEBV, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Trait + Condition, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Common High-Accuracy Traits: GEBV Distributions by Condition and Group",
    x = "GEBV",
    y = "Density",
    fill = "Group"
  )

################################################################################################
# Supplemental Figure 18
################################################################################################

common_traits <- c(
  "biomassfinalplant", "biomassfinalplot", "fruitlength", "fruitno", 
  "fruitshapeindex", "fruitweight", "fruitwidth", "leafangledaily", 
  "leafareadaily", "leafareadailyplant", "leafareaprojecteddaily", 
  "leafareaprojectedplant", "yield"
)

common_traits_prefixed <- paste0("GEBV_", common_traits)

#############
#trait distributions within core only with median line
library(tidyverse)
library(dplyr)
library(ggplot2)

# Load Global (10k) GEBVs
gebv_global_control <- read_csv("~/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_all_73_Control_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Control", Group = "Global")

gebv_global_heat1 <- read_csv("~/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_all_73_Heat_Stress_1_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Heat Stress 1", Group = "Global")

gebv_global_heat2 <- read_csv("~/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_all_73_Heat_Stress_2_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Heat Stress 2", Group = "Global")

# Load Core (423) GEBVs
gebv_core_control <- read_csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_all_n423_73trait_control.csv") %>%
  mutate(Condition = "Control", Group = "Core")

gebv_core_heat1 <- read_csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_all_n423_73trait_Heatstress_1.csv") %>%
  mutate(Condition = "Heat Stress 1", Group = "Core")

gebv_core_heat2 <- read_csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_all_n423_73trait_Heatstress_2.csv") %>%
  mutate(Condition = "Heat Stress 2", Group = "Core")

# Combine all datasets
gebv_all <- bind_rows(
  gebv_global_control,
  gebv_global_heat1,
  gebv_global_heat2,
  gebv_core_control,
  gebv_core_heat1,
  gebv_core_heat2
)

# Pivot to long format
gebv_all_long <- gebv_all %>%
  pivot_longer(cols = all_of(common_traits_prefixed), names_to = "Trait", values_to = "GEBV")


# 1. Filter for core lines only
gebv_core_only <- gebv_all_long %>%
  filter(Group == "Core") %>%
  mutate(Trait = str_remove(Trait, "GEBV_"))  # Optional cleanup

medians_df <- gebv_core_only %>%
  group_by(Trait, Condition) %>%
  summarize(median_gebv = median(GEBV, na.rm = TRUE), .groups = "drop")

ggplot(gebv_core_only, aes(x = GEBV, fill = Condition)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Trait, scales = "free", ncol = 4) +
  geom_vline(
    data = medians_df,
    aes(xintercept = median_gebv, color = Condition),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#1f77b4",       # blue
      "Heat Stress 1" = "#ff7f0e", # orange
      "Heat Stress 2" = "red"      # red
    )
  ) +
  scale_color_manual(
    values = c(
      "Control" = "#1f77b4",
      "Heat Stress 1" = "#ff7f0e",
      "Heat Stress 2" = "red"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "GEBV",
    y = "Density",
    fill = "Condition"
  )


#fix labels
# Define your raw and pretty labels
common_traits <- c(
  "biomassfinalplant", "biomassfinalplot", "fruitlength", "fruitno", 
  "fruitshapeindex", "fruitweight", "fruitwidth", "leafangledaily", 
  "leafareadaily", "leafareadailyplant", "leafareaprojecteddaily", 
  "leafareaprojectedplant", "yield"
)

pretty_labels <- c(
  "biomassfinalplant"      = "biomass final plant",
  "biomassfinalplot"       = "biomass final plot",
  "fruitlength"            = "fruit length",
  "fruitno"                = "fruit no",
  "fruitshapeindex"        = "fruit shape index",
  "fruitweight"            = "fruit weight",
  "fruitwidth"             = "fruit width",
  "leafangledaily"         = "leaf angle daily",
  "leafareadaily"          = "leaf area daily",
  "leafareadailyplant"     = "leaf area daily plant",
  "leafareaprojecteddaily" = "leaf area projected daily",
  "leafareaprojectedplant" = "leaf area projected plant",
  "yield"                  = "yield"
)

# Use labeller in facet_wrap
p<- ggplot(gebv_core_only, aes(x = GEBV, fill = Condition)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Trait, scales = "free", ncol = 4,
             labeller = labeller(Trait = pretty_labels)) +   # <-- HERE
  geom_vline(
    data = medians_df,
    aes(xintercept = median_gebv, color = Condition),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#1f77b4",
      "Heat Stress 1" = "#ff7f0e",
      "Heat Stress 2" = "red"
    )
  ) +
  scale_color_manual(
    values = c(
      "Control" = "#1f77b4",
      "Heat Stress 1" = "#ff7f0e",
      "Heat Stress 2" = "red"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "GEBV",
    y = "Density",
    fill = "Condition"
  )

p

ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/F2B_GEBV_core_density.pdf",
  plot     = p,
  width    = 15,     # adjust as you like
  height   = 8,     # adjust as you like
  units    = "in"
)

################################################################################################
# Figure 3A - heat loving lines
################################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Import GEBV files
GEBVs_control <- read.csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_Fruit_Yield_Control_Merged.csv", row.names = 1)
GEBVs_heat1   <- read.csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_Fruit_Yield_Heat_Stress_1_Merged.csv", row.names = 1)
GEBVs_heat2   <- read.csv("~/R/World_Veg_Project/data/outputs/core423/GEBVs_Fruit_Yield_Heat_Stress_2_Merged.csv", row.names = 1)

# Merge into one dataframe
GEBVs_all <- data.frame(
  Sample_ID   = rownames(GEBVs_control),
  GEBV_Heat1  = GEBVs_heat1[, 1],
  GEBV_Control= GEBVs_control[, 1],
  GEBV_Heat2  = GEBVs_heat2[, 1]
)

# Rank genotypes (lower GEBV = better rank) 
GEBVs_all <- GEBVs_all %>%
  mutate(
    Rank_Heat1   = rank(-GEBV_Heat1,   ties.method = "average"),
    Rank_Control = rank(-GEBV_Control, ties.method = "average"),
    Rank_Heat2   = rank(-GEBV_Heat2,   ties.method = "average")
  )

# Reshape for plotting 
GEBVs_long <- GEBVs_all %>%
  select(Sample_ID, Rank_Heat1, Rank_Control, Rank_Heat2) %>%
  pivot_longer(cols = c("Rank_Heat1", "Rank_Control", "Rank_Heat2"),
               names_to = "Condition", values_to = "Rank")

# Ensure condition order
GEBVs_long$Condition <- factor(GEBVs_long$Condition, 
                               levels = c("Rank_Heat1", "Rank_Control", "Rank_Heat2"))

# Identify extreme V-shaped genotypes 
GEBVs_all <- GEBVs_all %>%
  mutate(Total_V_Intensity = abs(Rank_Control - Rank_Heat1) + abs(Rank_Control - Rank_Heat2))

threshold <- quantile(GEBVs_all$Total_V_Intensity, 0.80)  # top 20%

extreme_V_genotypes <- GEBVs_all %>%
  filter(Rank_Heat1 < Rank_Control & Rank_Heat2 < Rank_Control & Total_V_Intensity >= threshold)

# Prepare plotting data
v_genotype_plot_data <- GEBVs_long %>%
  filter(Sample_ID %in% extreme_V_genotypes$Sample_ID)

# Final Plot 
ggplot(v_genotype_plot_data, 
       aes(x = Condition, y = Rank, group = Sample_ID, color = Sample_ID)) +
  geom_line(alpha = 0.8, linewidth = 1.1) +
  geom_point(size = 2) +
  scale_y_reverse() +  # Lower rank = better
  scale_color_manual(values = rep("firebrick3", length(unique(v_genotype_plot_data$Sample_ID)))) +
  labs(
    title = "                                                    Heat tolerant", 
    x = "Condition", 
    y = "Rank"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

#write.csv(extreme_V_genotypes$Sample_ID, 
#"~/R/World_Veg_Collab_Pepper/outputs/Extreme_V_SampleIDs.csv", 
#row.names = FALSE)


#bigger axis
p<-ggplot(v_genotype_plot_data, 
       aes(x = Condition, y = Rank, group = Sample_ID, color = Sample_ID)) +
  geom_line(alpha = 0.8, linewidth = 1.1) +
  geom_point(size = 2) +
  scale_y_reverse() +
  scale_color_manual(values = rep("firebrick3", length(unique(v_genotype_plot_data$Sample_ID)))) +
  labs(
    title = "                                                    Heat tolerant", 
    x = "Condition", 
    y = "Rank"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),        # tick labels bigger
    axis.title.x = element_text(size = 18),     # X-axis title bigger
    axis.title.y = element_text(size = 18),     # Y-axis title bigger
    plot.title = element_text(size = 16, face = "bold")
  )

p

ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/heat_tolerant_28.pdf",
  plot = p,
  width = 10,
  height = 7,
  units = "in"
)


################################################################################################
# Figure 3C (cold-loving)
################################################################################################
#continue from above
######################
library(dplyr)
library(ggplot2)
library(tidyr)

# Compute rank change intensity for inverted V-shaped patterns
GEBVs_all <- GEBVs_all %>%
  mutate(
    Total_InvertedV_Intensity = abs(Rank_Control - Rank_Heat1) + abs(Rank_Control - Rank_Heat2)
  )

# Set threshold (top 20% of Total_InvertedV_Intensity)
threshold_invertedV <- quantile(GEBVs_all$Total_InvertedV_Intensity, 0.80)

# Extract genotypes showing the strongest Inverted V-shape
extreme_invertedV_genotypes <- GEBVs_all %>%
  filter(Rank_Heat1 > Rank_Control & Rank_Heat2 > Rank_Control & Total_InvertedV_Intensity >= threshold_invertedV)

# Save to CSV (optional)
# write.csv(extreme_invertedV_genotypes, "~/R/World_Veg_Collab_Pepper/outputs/Extreme_InvertedV_Genotypes.csv", row.names = FALSE)


# Prepare data for plotting
GEBVs_long$Condition <- factor(GEBVs_long$Condition, levels = c("Rank_Heat1", "Rank_Control", "Rank_Heat2"))

# Plot Inverted V-shaped genotypes
ggplot(GEBVs_long %>% filter(Sample_ID %in% extreme_invertedV_genotypes$Sample_ID), 
       aes(x = Condition, y = Rank, group = Sample_ID, color = Sample_ID)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  scale_y_reverse() +  # Lower rank is better
  labs(title = "", 
       x = "Condition", y = "Rank (Lower is Better)") +
  theme(legend.position = "none")

# Prepare data for plotting
invertedV_plot_data <- GEBVs_long %>%
  filter(Sample_ID %in% extreme_invertedV_genotypes$Sample_ID)

# Plot
ggplot(invertedV_plot_data, 
       aes(x = Condition, y = Rank, group = Sample_ID, color = Sample_ID)) +
  geom_line(alpha = 0.8, linewidth = 1.1) +
  geom_point(size = 2) +
  scale_y_reverse() +  # Lower rank is better
  scale_color_manual(values = rep("royalblue4", length(unique(invertedV_plot_data$Sample_ID)))) +
  labs(
    title = "                                                     Heat sensitive", 
    x = "Condition", 
    y = "Rank "
  ) +
  theme_minimal() +  # Clean white background
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

#bigger axis
p<-ggplot(invertedV_plot_data, 
       aes(x = Condition, y = Rank, group = Sample_ID, color = Sample_ID)) +
  geom_line(alpha = 0.8, linewidth = 1.1) +
  geom_point(size = 2) +
  scale_y_reverse() +  
  scale_color_manual(values = rep("royalblue4", length(unique(invertedV_plot_data$Sample_ID)))) +
  labs(
    title = "                                                     Heat sensitive", 
    x = "Condition", 
    y = "Rank "
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),        # tick labels
    axis.title.x = element_text(size = 18),     # X-axis label bigger
    axis.title.y = element_text(size = 18),     # Y-axis label bigger
    plot.title = element_text(size = 16, face = "bold")
  )

ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/cold_tolerant_57.pdf",
  plot = p,
  width = 10,
  height = 7,
  units = "in"
)

#write.csv(extreme_V_genotypes$Sample_ID, 
         # "~/R/World_Veg_Collab_Pepper/outputs/Extreme_V_SampleIDs.csv", 
         # row.names = FALSE)

#write.csv(extreme_invertedV_genotypes$Sample_ID, 
         # "~/R/World_Veg_Collab_Pepper/outputs/Extreme_InvertedV_SampleIDs.csv", 
        # row.names = FALSE)

########################################################################################################################
# Required preface for the marker effect changes for Figure 2B and 2D
########################################################################################################################
#Individual Marker Effect Sizes - YIELD - ctrl timepoint
########################################
library(rrBLUP)
library(dplyr)
library(tidyr)

# 1. Use already-loaded genotype matrix
#geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_df <- geno_df_imputed  # already in memory and imputed

# 2. Load phenotype (yield)
#pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
pheno <- control_pheno_check  # already loaded

# 3. Filter both to common lines
common_ids <- intersect(rownames(geno_df), rownames(pheno))
geno_df <- geno_df[common_ids, ]
pheno <- pheno[common_ids, ]

# 4. Load training genotypes
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# 5. Subset training data
geno_train <- geno_df[rownames(geno_df) %in% train_genotypes, ]
pheno_train <- pheno[rownames(pheno) %in% train_genotypes, ]

# 6. Run rrBLUP
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"
yield_values <- pheno_train$yield
solve_yield <- mixed.solve(y = yield_values, Z = geno_train)
u.hat_yield <- solve_yield$u  # SNP effects

# 7. Apply to all individuals
geno_all <- as.matrix(geno_df)
mode(geno_all) <- "numeric"
individual_marker_effects <- geno_all * matrix(u.hat_yield, nrow = nrow(geno_all), ncol = ncol(geno_all), byrow = TRUE)

# 8. Reshape + export
individual_marker_effects_df <- as.data.frame(individual_marker_effects)
individual_marker_effects_df$Individual <- rownames(geno_all)

individual_marker_effects_long <- individual_marker_effects_df %>%
  pivot_longer(cols = -Individual, names_to = "SNP", values_to = "Effect")

write.csv(individual_marker_effects_long,
          "~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP.csv",
          row.names = FALSE)

########################################
#Individual Marker Effect Sizes - YIELD - heat1 timepoint
########################################
library(rrBLUP)
library(dplyr)
library(tidyr)

# 1. Use already-loaded genotype matrix
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
#geno_df <- geno_df_imputed  # already in memory and imputed

# 2. Load phenotype (yield)
#pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)


# 3. Filter both to common lines
common_ids <- intersect(rownames(geno_df), rownames(pheno))
geno_df <- geno_df[common_ids, ]
pheno <- pheno[common_ids, ]

# 4. Load training genotypes
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# 5. Subset training data
geno_train <- geno_df[rownames(geno_df) %in% train_genotypes, ]
pheno_train <- pheno[rownames(pheno) %in% train_genotypes, ]

# 6. Run rrBLUP
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"
yield_values <- pheno_train$yield
solve_yield <- mixed.solve(y = yield_values, Z = geno_train)
u.hat_yield <- solve_yield$u  # SNP effects

# 7. Apply to all individuals
geno_all <- as.matrix(geno_df)
mode(geno_all) <- "numeric"
individual_marker_effects <- geno_all * matrix(u.hat_yield, nrow = nrow(geno_all), ncol = ncol(geno_all), byrow = TRUE)

# 8. Reshape + export
individual_marker_effects_df <- as.data.frame(individual_marker_effects)
individual_marker_effects_df$Individual <- rownames(geno_all)

individual_marker_effects_long <- individual_marker_effects_df %>%
  pivot_longer(cols = -Individual, names_to = "SNP", values_to = "Effect")

write.csv(individual_marker_effects_long,
          "~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP_heat1.csv",
          row.names = FALSE)


########################################
#Individual Marker Effect Sizes - YIELD - heat2 timepoint
########################################
library(rrBLUP)
library(dplyr)
library(tidyr)

# 1. Use already-loaded genotype matrix
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
#geno_df <- geno_df_imputed  # already in memory and imputed

# 2. Load phenotype (yield)
#pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)


# 3. Filter both to common lines
common_ids <- intersect(rownames(geno_df), rownames(pheno))
geno_df <- geno_df[common_ids, ]
pheno <- pheno[common_ids, ]

# 4. Load training genotypes
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# 5. Subset training data
geno_train <- geno_df[rownames(geno_df) %in% train_genotypes, ]
pheno_train <- pheno[rownames(pheno) %in% train_genotypes, ]

# 6. Run rrBLUP
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"
yield_values <- pheno_train$yield
solve_yield <- mixed.solve(y = yield_values, Z = geno_train)
u.hat_yield <- solve_yield$u  # SNP effects

# 7. Apply to all individuals
geno_all <- as.matrix(geno_df)
mode(geno_all) <- "numeric"
individual_marker_effects <- geno_all * matrix(u.hat_yield, nrow = nrow(geno_all), ncol = ncol(geno_all), byrow = TRUE)

# 8. Reshape + export
individual_marker_effects_df <- as.data.frame(individual_marker_effects)
individual_marker_effects_df$Individual <- rownames(geno_all)

individual_marker_effects_long <- individual_marker_effects_df %>%
  pivot_longer(cols = -Individual, names_to = "SNP", values_to = "Effect")

write.csv(individual_marker_effects_long,
          "~/R/World_Veg_Collab_Pepper/outputs/individual_marker_effects_yield_per_SNP_heat2.csv",
          row.names = FALSE)

################################################################################################
# Figure 3B (Hot 28)
################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

control <- read_csv("~/R/World_Veg_Project/data/outputs/individual_marker_effects_yield_per_SNP.csv",
                    col_types = cols(
                      Individual = col_character(),
                      Effect = col_double(),
                      SNP = col_character()  # <- critical fix
                    ))

heat1 <- read_csv("~/R/World_Veg_Project/data/outputs/individual_marker_effects_yield_per_SNP_heat1.csv",
                  col_types = cols(
                    Individual = col_character(),
                    Effect = col_double(),
                    SNP = col_character()
                  ))

# Load and fix sample ID column
hot_28  <- read_csv("~/R/World_Veg_Project/data/outputs/Extreme_V_SampleIDs.csv")
colnames(hot_28)[1] <- "Sample_ID"

# Subset to the 28 heat-loving lines
selected_ids <- hot_28$Sample_ID
control_28 <- control %>% filter(Individual %in% selected_ids)
heat1_28   <- heat1   %>% filter(Individual %in% selected_ids)

# Merge control and heat1 on Individual and SNP
merged <- control_28 %>%
  dplyr::rename(Effect_control = Effect) %>%
  inner_join(
    heat1_28 %>% dplyr::rename(Effect_heat1 = Effect),
    by = c("Individual", "SNP")
  )


# Step 4: Create wide format for paired t-tests
wide <- merged %>%
  pivot_wider(names_from = Individual, values_from = c(Effect_control, Effect_heat1))

# Step 5: Calculate t-tests per SNP (with Bonferroni correction)
test_results <- merged %>%
  group_by(SNP) %>%
  summarize(
    t_stat = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$statistic
    } else { NA },
    p_value = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$p.value
    } else { NA },
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "bonferroni"),
    significant = adjusted_p_value < 0.001
  )


# Step 6: Extract significant SNPs and plot
significant_snps <- merged %>%
  filter(SNP %in% test_results$SNP[test_results$significant]) %>%
  group_by(SNP) %>%
  summarize(
    Mean_Effect_Heat1 = mean(Effect_heat1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Chromosome = gsub(":.*", "", SNP),
    Position = as.numeric(gsub(".*:", "", SNP)),
    Position_Mb = Position / 1e6
  ) %>%
  arrange(as.numeric(Chromosome), Position) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = as.character(1:12)),
    Genomic_Position = row_number()
  )

significant_snps <- significant_snps %>%
  filter(abs(Mean_Effect_Heat1) > 0.001)

#write_csv(significant_snps, "~/R/World_Veg_Collab_Pepper/outputs/significant_snps_hot28.csv")

# Step 7: Plot
plot_faceted <- ggplot(significant_snps, aes(x = Position / 1e6, y = Mean_Effect_Heat1)) +
  geom_point(color = "red4", size = 1.2, alpha = 0.7) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  labs(
    title = "Significant SNP Marker Effects by Chromosome",
    x = "Genomic Position (Mb)",
    y = "Mean Effect Size (Heat1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(plot_faceted)

# Optional save
#ggsave("~/R/World_Veg_Collab_Pepper/FIGS/Heat_Loving_SNPs_Heat1_vs_Control_315snp.pdf",
#plot = plot_faceted, width = 10, height = 6)

################################################################################################
# Figure 3D (Cold 57)
################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Step 1: Load full marker effect datasets
control <- read_csv("~/R/World_Veg_Project/data/outputs/individual_marker_effects_yield_per_SNP.csv",
                    col_types = cols(
                      Individual = col_character(),
                      Effect = col_double(),
                      SNP = col_character()
                    ))

heat1 <- read_csv("~/R/World_Veg_Project/data/outputs/individual_marker_effects_yield_per_SNP_heat1.csv",
                  col_types = cols(
                    Individual = col_character(),
                    Effect = col_double(),
                    SNP = col_character()
                  ))

# Step 2: Load cold-tolerant line IDs
cold_57 <- read_csv("~/R/World_Veg_Project/data/outputs/Extreme_InvertedV_SampleIDs.csv")
colnames(cold_57)[1] <- "Sample_ID"
selected_ids <- cold_57$Sample_ID

# Step 3: Subset genotype effect data to only cold-tolerant lines
control_57 <- control %>% filter(Individual %in% selected_ids)
heat1_57   <- heat1   %>% filter(Individual %in% selected_ids)

# Step 4: Merge datasets by SNP and Individual
merged <- control_57 %>%
  dplyr::rename(Effect_control = Effect) %>%
  inner_join(
    heat1_57 %>% dplyr::rename(Effect_heat1 = Effect),
    by = c("Individual", "SNP")
  )

# Step 5: (Optional) Wide format if needed
wide <- merged %>%
  pivot_wider(names_from = Individual, values_from = c(Effect_control, Effect_heat1))

# Step 6: Paired t-tests for each SNP
test_results <- merged %>%
  group_by(SNP) %>%
  summarize(
    t_stat = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$statistic
    } else { NA },
    p_value = if (length(unique(Effect_control)) > 1 && length(unique(Effect_heat1)) > 1) {
      t.test(Effect_control, Effect_heat1, paired = TRUE)$p.value
    } else { NA },
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "bonferroni"),
    significant = adjusted_p_value < 0.001
  )

# Step 7: Summarize significant SNPs
significant_snps <- merged %>%
  filter(SNP %in% test_results$SNP[test_results$significant]) %>%
  group_by(SNP) %>%
  summarize(
    Mean_Effect_Heat1 = mean(Effect_heat1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Chromosome = gsub(":.*", "", SNP),
    Position = as.numeric(gsub(".*:", "", SNP)),
    Position_Mb = Position / 1e6
  ) %>%
  arrange(as.numeric(Chromosome), Position) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = as.character(1:12)),
    Genomic_Position = row_number()
  )

# Optional: Filter by effect size threshold
significant_snps <- significant_snps %>%
  filter(abs(Mean_Effect_Heat1) > 0.001)

#write_csv(significant_snps, "~/R/World_Veg_Collab_Pepper/outputs/significant_snps_cold57.csv")

# Step 8: Plot faceted by chromosome
plot_faceted <- ggplot(significant_snps, aes(x = Position_Mb, y = Mean_Effect_Heat1)) +
  geom_point(color = "blue4", size = 1.2, alpha = 0.7) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  labs(
    title = "Significant SNP Marker Effects by Chromosome (Cold-Loving Lines)",
    x = "Genomic Position (Mb)",
    y = "Mean Effect Size (Heat1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
plot_faceted

# Step 9: Save plot
#ggsave("~/R/World_Veg_Collab_Pepper/FIGS/Cold_Loving_SNPs_Heat1_vs_Control_57lines.pdf",
#plot = plot_faceted, width = 10, height = 6)

################################################################################################
# Figure 4A – Climate vs Quality GEBV correlations
################################################################################################
library(tidyverse)
library(pheatmap)
library(ggfortify)

# datasets 
gebv_control <- read_csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_all_n423_73trait_control.csv", show_col_types = FALSE)
gebv_heat1   <- read_csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_all_n423_73trait_Heatstress_1.csv", show_col_types = FALSE)
gebv_heat2   <- read_csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_all_n423_73trait_Heatstress_2.csv", show_col_types = FALSE)

gebv_quality <- read_csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv",
                         show_col_types = FALSE)
colnames(gebv_quality)[1] <- "Line"

# Standardize IDs
gebv_control <- gebv_control %>% mutate(Line = trimws(as.character(Line)))
gebv_heat1   <- gebv_heat1   %>% mutate(Line = trimws(as.character(Line)))
gebv_heat2   <- gebv_heat2   %>% mutate(Line = trimws(as.character(Line)))
gebv_quality <- gebv_quality %>% mutate(Line = trimws(as.character(Line)))


gebv_quality_long <- gebv_quality %>%
  select(-Set) %>%
  pivot_longer(-Line, names_to = "QualityTrait", values_to = "QualityGEBV") %>%
  mutate(QualityGEBV = as.numeric(QualityGEBV))   # <-- ensure numeric



# Helper for correlations ---
get_correlations <- function(climate_df, quality_long) {
  climate_long <- climate_df %>%
    pivot_longer(-Line, names_to = "Trait", values_to = "GEBV")
  
  climate_avg <- climate_long %>%
    group_by(Line) %>%
    summarise(MeanClimateGEBV = mean(GEBV, na.rm = TRUE), .groups = "drop")
  
  merged <- inner_join(climate_avg, quality_long, by = "Line")
  message("Merged rows: ", nrow(merged))  # debug check
  
  if (nrow(merged) == 0) {
    stop("Join failed: no overlapping Line IDs between climate and quality datasets.")
  }
  
  cor_df <- merged %>%
    group_by(QualityTrait) %>%
    summarise(Correlation = cor(MeanClimateGEBV, QualityGEBV, use = "complete.obs"),
              .groups = "drop")
  
  return(cor_df)
}

# Compute correlations 
cor_control <- get_correlations(gebv_control, gebv_quality_long) %>%
  dplyr::rename_with(~ "Control", "Correlation")

cor_heat1 <- get_correlations(gebv_heat1, gebv_quality_long) %>%
  dplyr::rename_with(~ "Heat1", "Correlation")

cor_heat2 <- get_correlations(gebv_heat2, gebv_quality_long) %>%
  dplyr::rename_with(~ "Heat2", "Correlation")



# Merge correlation tables 
cor_matrix <- cor_control %>%
  inner_join(cor_heat1, by = "QualityTrait") %>%
  inner_join(cor_heat2, by = "QualityTrait")

pca_input <- cor_matrix %>%
  column_to_rownames("QualityTrait") %>%
  scale()

# PCA plot 
pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)
autoplot(pca_res, loadings = TRUE, loadings.label = TRUE) +
  theme_minimal() +
  labs(title = "PCA of Quality Trait Correlation Dynamics",
       subtitle = "Trait grouping across Control, Heat Stress 1, Heat Stress 2")

# Heatmap 
pheatmap(pca_input,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Clustering of Quality Traits by Climate Correlation Pattern",
         color = colorRampPalette(c("blue", "white", "red"))(100))

# fix plot labels

# --- Clean labels before heatmap ---
rownames(pca_input) <- rownames(pca_input) %>%
  sub("^GEBV_", "", .) %>%   # remove prefix
  gsub("_", " ", .)          # replace underscores with spaces

# --- Heatmap with cleaned labels ---
pheatmap(pca_input,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "",
         color = colorRampPalette(c("blue", "white", "red"))(100))


################################################################################################
# Figure 4B
################################################################################################
#all 23 traits 
library(tidyverse)
library(dplyr)
library(readr)




# Core (423) quality GEBVs
gebv_quality_core <- read_csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv") %>%
  rename(Line = `...1`) %>%
  select(-Set) %>%
  mutate(Group = "Core")

# Global (10k) quality GEBVs
gebv_quality_global <- read_csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_McLeod_23traits_10k_ALL.csv") %>%
  rename(Line = `...1`) %>%
  select(-Set) %>%
  mutate(Group = "Global")

# Combine and pivot to long format in one step
gebv_quality_all_long <- bind_rows(gebv_quality_core, gebv_quality_global) %>%
  pivot_longer(cols = -c(Line, Group), names_to = "Trait", values_to = "GEBV")

# Plot GEBV distributions
ggplot(gebv_quality_all_long, aes(x = GEBV, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Trait, scales = "free", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(title = "Quality Trait GEBV Distributions: Core (423) vs Global (10k)",
       x = "GEBV",
       y = "Density",
       fill = "Group")
################################################################################
# above 0.5 PA for quality across core and global plots
library(dplyr)
library(readr)

# Load only rrBLUP results
rrblup_core <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_rrblup_23traits_kfold_10.RDS")$xval.result %>%
  mutate(Source = "Core")

rrblup_global <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_10k/xval_rrblup_23traits_kfold_10_10k.RDS")$xval.result %>%
  mutate(Source = "Global")

# Combine and keep only rrBLUP
rrblup_all <- bind_rows(rrblup_core, rrblup_global)

# Ensure numeric and factor types
rrblup_all <- rrblup_all %>%
  mutate(
    r.mean = as.numeric(r.mean),
    trait = as.factor(trait)
  )

# Filter for traits with r.mean > 0.5
rrblup_above_0.5 <- rrblup_all %>%
  filter(r.mean > 0.5)

# Find traits that are above 0.5 in BOTH Core and Global
common_high_accuracy_traits <- rrblup_above_0.5 %>%
  group_by(trait) %>%
  summarize(n_sources = n_distinct(Source)) %>%
  filter(n_sources == 2) %>%
  pull(trait)

# View the traits
print(common_high_accuracy_traits)
#################################################################################
library(tidyverse)

# Define the 16 high-accuracy traits
common_high_accuracy_traits <- c(
  "GEBV_Brix", "GEBV_Axis_length", "GEBV_Immature_fruit_external_color_L",
  "GEBV_External_immature_fruit_color_green", "GEBV_Fruit_external_color_a",
  "GEBV_Fruit_fasciation", "GEBV_Fruit_maximum_length", "GEBV_Fruit_maximum_width",
  "GEBV_Fruit_pungency", "GEBV_Fruit_shape_index", "GEBV_Fruit_weight",
  "GEBV_Locule_number", "GEBV_Total_plant_height", "GEBV_Pericarp_thickness",
  "GEBV_Total_fruit_weight", "GEBV_Total_fruit_number"
)

# --- Load Core (423) ---
gebv_quality_core <- read_csv(
  "/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv",
  show_col_types = FALSE
)
colnames(gebv_quality_core)[1] <- "Line"   # fix first col name
gebv_quality_core <- gebv_quality_core %>%
  select(-Set) %>%
  mutate(Group = "Core")

# --- Load Global (10k) ---
gebv_quality_global <- read_csv(
  "/Users/annamccormick/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_McLeod_23traits_10k_ALL.csv",
  show_col_types = FALSE
)
colnames(gebv_quality_global)[1] <- "Line"  # fix first col name
gebv_quality_global <- gebv_quality_global %>%
  select(-Set) %>%
  mutate(Group = "Global")

# --- Combine & pivot ---
gebv_quality_all_long <- bind_rows(gebv_quality_core, gebv_quality_global) %>%
  pivot_longer(cols = -c(Line, Group), names_to = "Trait", values_to = "GEBV")

# --- Filter for high-accuracy traits ---
gebv_quality_filtered <- gebv_quality_all_long %>%
  filter(Trait %in% common_high_accuracy_traits)

# Density plot ---
ggplot(gebv_quality_filtered, aes(x = GEBV, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Trait, scales = "free", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(
    title = "High-Accuracy Traits: Core vs Global GEBV Distributions",
    x = "GEBV",
    y = "Density",
    fill = "Group"
  )

# plot
ggplot(gebv_quality_filtered, aes(x = GEBV, fill = Group)) +
  geom_density(position = "identity", alpha = 0.4, bins = 30) +
  facet_wrap(~ Trait, scales = "free", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(
    title = "High-Accuracy Traits: Core vs Global GEBV Distributions",
    x = "GEBV",
    y = "Density",
    fill = "Group"
  )


#fix plot labels
# Create nicer labels: drop "GEBV_" and replace "_" with spaces
label_vec <- setNames(
  gsub("_", " ", sub("^GEBV_", "", common_high_accuracy_traits)),
  common_high_accuracy_traits
)

# Density plot with clean facet labels
ggplot(gebv_quality_filtered, aes(x = GEBV, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Trait, scales = "free", ncol = 4,
             labeller = labeller(Trait = label_vec)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "GEBV",
    y = "Density",
    fill = "Group"
  )

############
#Fig 4B again


################################################################################
# Libraries
################################################################################
library(tidyverse)
library(dplyr)
library(readr)

################################################################################
# Safe import function for GEBV CSV files
################################################################################
load_quality_csv <- function(path, group_name) {
  df <- read_csv(path, show_col_types = FALSE)
  
  # Rename FIRST column to "Line", regardless of its weird name
  names(df)[1] <- "Line"
  
  df %>%
    select(-Set) %>%
    mutate(Group = group_name)
}

################################################################################
# Load Core (423) & Global (10k) datasets
################################################################################

gebv_quality_core <- load_quality_csv(
  "/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv",
  "Core"
)

gebv_quality_global <- load_quality_csv(
  "/Users/annamccormick/R/World_Veg_Project/data/outputs/10250_GS/GEBVs_McLeod_23traits_10k_ALL.csv",
  "Global"
)

################################################################################
# Combine and pivot long
################################################################################

gebv_quality_all_long <- bind_rows(gebv_quality_core, gebv_quality_global) %>%
  pivot_longer(cols = -c(Line, Group),
               names_to = "Trait",
               values_to = "GEBV")

################################################################################
# Load rrBLUP prediction accuracy for Core + Global
################################################################################

rrblup_core <- readRDS(
  "/Users/annamccormick/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_rrblup_23traits_kfold_10.RDS"
)$xval.result %>%
  mutate(Source = "Core")

rrblup_global <- readRDS(
  "/Users/annamccormick/R/World_Veg_Project/data/outputs/PA_23trait_10k/xval_rrblup_23traits_kfold_10_10k.RDS"
)$xval.result %>%
  mutate(Source = "Global")

rrblup_all <- bind_rows(rrblup_core, rrblup_global) %>%
  mutate(
    r.mean = as.numeric(r.mean),
    trait  = as.character(trait)
  )

################################################################################
# Identify traits with > 0.5 accuracy in BOTH datasets
################################################################################

common_high_accuracy_traits <- rrblup_all %>%
  filter(r.mean > 0.5) %>%
  group_by(trait) %>%
  summarise(n_sources = n_distinct(Source)) %>%
  filter(n_sources == 2) %>%
  pull(trait)

print("Traits with >0.5 PA in both Core and Global:")
print(common_high_accuracy_traits)

################################################################################
# FIX — Add GEBV_ prefix to match column names in your GEBV tables
################################################################################

common_high_accuracy_traits_prefixed <- paste0("GEBV_", common_high_accuracy_traits)

################################################################################
# Filter long GEBV data for these traits
################################################################################

gebv_quality_filtered <- gebv_quality_all_long %>%
  filter(Trait %in% common_high_accuracy_traits_prefixed)

################################################################################
# Clean facet labels
################################################################################

label_vec <- setNames(
  gsub("_", " ", sub("^GEBV_", "", common_high_accuracy_traits_prefixed)),
  common_high_accuracy_traits_prefixed
)

################################################################################
# Plot
################################################################################

p <- ggplot(gebv_quality_filtered, aes(x = GEBV, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(
    ~ Trait,
    scales = "free",
    ncol = 4,
    labeller = labeller(Trait = label_vec)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "GEBV",
    y = "Density",
    fill = "Group"
  )

print(p)

################################################################################
# Save as PDF
################################################################################

ggsave(
  "/Users/annamccormick/R/World_Veg_Project/PANNELS/4B.pdf",
  plot = p,
  width = 20,
  height = 10,
  units = "in"
)




####### FIN










