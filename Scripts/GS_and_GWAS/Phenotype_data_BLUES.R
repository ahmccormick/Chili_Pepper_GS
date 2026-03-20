######################################################
#LAYOUT
#Filter raw data to 286 intersect IDs → three filtered raw files
#Within-environment BLUEs → BLUEs_Control.csv, BLUEs_HS1.csv, BLUEs_HS2.csv
#STI from those three BLUE files → STI_from_BLUEs_all_73_traits.csv

######################################################
# Raw Phenotype data subsetting to the 286 lines common across all environmental timepoints. 
######################################################

library(dplyr)
library(readr)

raw_ctrl <- "~/R/World_Veg_Project/anna_BLUES/Control_Pheno_data.csv"
raw_hs1  <- "~/R/World_Veg_Project/anna_BLUES/Heat_stress_1_data.csv"
raw_hs2  <- "~/R/World_Veg_Project/anna_BLUES/Heat_stress_2_data.csv"

f286_ctrl <- "~/R/World_Veg_Project/anna_BLUES/control_pheno_filtered_286.csv"
f286_hs1  <- "~/R/World_Veg_Project/anna_BLUES/heat_stress_1_filtered_286.csv"
f286_hs2  <- "~/R/World_Veg_Project/anna_BLUES/heat_stress_2_filtered_286.csv"

outdir <- "anna_BLUES_rebuild"
#dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read the 286 file and extract accession IDs 
get_ids_from_286 <- function(path_286) {
  x <- read.csv(path_286, stringsAsFactors = FALSE, check.names = FALSE)
  # accession IDs are in the first column (often named "" or "Unnamed: 0")
  ids <- trimws(as.character(x[[1]]))
  ids[!is.na(ids) & ids != ""]
}

ids_ctrl <- get_ids_from_286(f286_ctrl)
ids_hs1  <- get_ids_from_286(f286_hs1)
ids_hs2  <- get_ids_from_286(f286_hs2)

# usually these should be (almost) identical; use union or intersect as you prefer
ids_union <- sort(unique(c(ids_ctrl, ids_hs1, ids_hs2)))
ids_intersect <- Reduce(intersect, list(ids_ctrl, ids_hs1, ids_hs2))

cat("IDs in ctrl_286:", length(ids_ctrl), "\n")
cat("IDs in hs1_286 :", length(ids_hs1), "\n")
cat("IDs in hs2_286 :", length(ids_hs2), "\n")
cat("Union IDs      :", length(ids_union), "\n")
cat("Intersect IDs  :", length(ids_intersect), "\n\n")

# ---- read raw RCBD files ----
ctrl_raw <- read.csv(raw_ctrl, stringsAsFactors = FALSE, check.names = FALSE)
hs1_raw  <- read.csv(raw_hs1,  stringsAsFactors = FALSE, check.names = FALSE)
hs2_raw  <- read.csv(raw_hs2,  stringsAsFactors = FALSE, check.names = FALSE)

# standardize id and rep column names
library(dplyr)

standardize_raw <- function(df, env_label) {
  
  # trim any accidental spaces in column names
  names(df) <- trimws(names(df))
  
  # --- ID column: g2pname -> accession ---
  if (!("g2pname" %in% names(df))) {
    stop("No g2pname column found. Columns are: ", paste(names(df), collapse = ", "))
  }
  df$g2pname <- trimws(as.character(df$g2pname))
  names(df)[names(df) == "g2pname"] <- "accession"
  
  # --- rep column: accept Rep or rep (case-insensitive) ---
  rep_idx <- which(tolower(names(df)) == "rep")
  if (length(rep_idx) != 1) {
    stop("Could not uniquely find rep/Rep column. Columns are: ", paste(names(df), collapse = ", "))
  }
  names(df)[rep_idx] <- "rep"
  df$rep <- trimws(as.character(df$rep))
  
  # --- environment label ---
  df$environment <- env_label
  
  df
}

# re-read raw files (fresh)
ctrl_raw <- read.csv(raw_ctrl, stringsAsFactors = FALSE, check.names = FALSE)
hs1_raw  <- read.csv(raw_hs1,  stringsAsFactors = FALSE, check.names = FALSE)
hs2_raw  <- read.csv(raw_hs2,  stringsAsFactors = FALSE, check.names = FALSE)

# standardize
ctrl_raw <- standardize_raw(ctrl_raw, "Control")
hs1_raw  <- standardize_raw(hs1_raw,  "HS1")
hs2_raw  <- standardize_raw(hs2_raw,  "HS2")

# sanity check
names(ctrl_raw)[1:5]
head(ctrl_raw[, c("accession","rep","environment")])


#FILTER raw data to the SAME accession set as the 286 files 

ids_keep <- ids_intersect

ctrl_filt <- ctrl_raw %>% dplyr::filter(accession %in% ids_keep)
hs1_filt  <- hs1_raw  %>% dplyr::filter(accession %in% ids_keep)
hs2_filt  <- hs2_raw  %>% dplyr::filter(accession %in% ids_keep)

cat("Filtered raw rows (Control):", nrow(ctrl_filt), "\n")
cat("Filtered raw rows (HS1)    :", nrow(hs1_filt), "\n")
cat("Filtered raw rows (HS2)    :", nrow(hs2_filt), "\n\n")

# quick check: how many reps (rows) per accession?
check_reps <- function(df, label) {
  tmp <- df %>%
    dplyr::count(accession, name = "n_rows") %>%
    dplyr::count(n_rows, name = "n_accessions") %>%
    dplyr::arrange(n_rows)
  cat("Rep-count distribution for", label, ":\n")
  print(tmp)
  cat("\n")
}

check_reps(ctrl_filt, "Control")
check_reps(hs1_filt,  "HS1")
check_reps(hs2_filt,  "HS2")

# write filtered raw RCBD files (ready for BLUEs)
write.csv(
  ctrl_filt,
  "~/R/World_Veg_Project/anna_BLUES/Control_raw_filtered_to_286IDs.csv",
  row.names = FALSE
)

write.csv(
  hs1_filt,
  "~/R/World_Veg_Project/anna_BLUES/HS1_raw_filtered_to_286IDs.csv",
  row.names = FALSE
)

write.csv(
  hs2_filt,
  "~/R/World_Veg_Project/anna_BLUES/HS2_raw_filtered_to_286IDs.csv",
  row.names = FALSE
)


################################################################################
# FULL PIPELINE: BLUEs from 3 raw environment phenotype files (RCBD)
################################################################################
library(lme4)
library(emmeans)
library(dplyr)
library(purrr)
library(readr)

# inputs
file_ctrl <- "~/R/World_Veg_Project/anna_BLUES/Control_raw_filtered_to_286IDs.csv"
file_hs1  <- "~/R/World_Veg_Project/anna_BLUES/HS1_raw_filtered_to_286IDs.csv"
file_hs2  <- "~/R/World_Veg_Project/anna_BLUES/HS2_raw_filtered_to_286IDs.csv"

outdir <- "~/R/World_Veg_Project/anna_BLUES/anna_BLUES_rebuild"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

accession_col <- "accession"
rep_col <- "rep"
extra_non_trait_cols <- c("genotype", "plantno.")
exclude_traits <- character(0)

env_ctrl <- "Control"
env_hs1  <- "HS1"
env_hs2  <- "HS2"


# Read + standardize

read_one_env <- function(path, env_label) {
  df <- read.csv(path, na.strings = c("NA", ".", "", " ", "NaN"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))
  
  df[[accession_col]] <- factor(trimws(df[[accession_col]]))
  df[[rep_col]]       <- factor(trimws(df[[rep_col]]))
  df$environment      <- factor(env_label)
  
  non_trait <- c(accession_col, rep_col, "environment", extra_non_trait_cols)
  trait_cols <- setdiff(names(df), non_trait)
  
  df <- df %>%
    mutate(across(all_of(trait_cols), ~ parse_number(as.character(.x))))
  
  df
}

ctrl <- read_one_env(file_ctrl, env_ctrl)
hs1  <- read_one_env(file_hs1, env_hs1)
hs2  <- read_one_env(file_hs2, env_hs2)

dat <- bind_rows(ctrl, hs1, hs2)


# Traits
traits_common <- Reduce(intersect, list(
  setdiff(names(ctrl), c(accession_col, rep_col, "environment", extra_non_trait_cols)),
  setdiff(names(hs1),  c(accession_col, rep_col, "environment", extra_non_trait_cols)),
  setdiff(names(hs2),  c(accession_col, rep_col, "environment", extra_non_trait_cols))
))

traits_common <- setdiff(traits_common, exclude_traits)
cat("Traits:", length(traits_common), "\n")  # should be 73

######################################################
# Within-environment BLUEs (one file per env)
######################################################
calc_blues_within_env <- function(trait, env_label) {
  d <- dat %>%
    filter(environment == env_label) %>%
    filter(!is.na(.data[[trait]]))
  
  # Need at least 2 accessions to estimate accession effects
  if (n_distinct(d[[accession_col]]) < 2) return(NULL)
  
  m <- lmer(
    as.formula(paste0(trait, " ~ ", accession_col, " + (1|", rep_col, ")")),
    data = d, REML = TRUE
  )
  
  emmeans(m, accession_col) %>%
    as.data.frame() %>%
    transmute(
      accession = .data[[accession_col]],
      !!trait := emmean
    )
}

write_env_blues <- function(env_label, out_file) {
  cat("Computing within-env BLUEs for", env_label, "...\n")
  
  blues_env <- traits_common %>%
    purrr::map(~ calc_blues_within_env(.x, env_label)) %>%
    purrr::compact() %>%
    purrr::reduce(dplyr::left_join, by = accession_col)
  
  write.csv(blues_env, file.path(outdir, out_file), row.names = FALSE)
  cat("Wrote:", file.path(outdir, out_file), "\n")
  blues_env
}


blues_ctrl <- write_env_blues("Control", "BLUEs_Control.csv")
blues_hs1  <- write_env_blues("HS1",     "BLUEs_HS1.csv")
blues_hs2  <- write_env_blues("HS2",     "BLUEs_HS2.csv")



############################################
# STI made from BLUE phenotype data
############################################
library(tidyverse)

# ── file paths ────────────────────────────────────────────────────────────────
blue_ctrl_file <- "~/R/World_Veg_Project/anna_BLUEs/anna_BLUES_rebuild/BLUEs_Control.csv"
blue_h1_file   <- "~/R/World_Veg_Project/anna_BLUEs/anna_BLUES_rebuild/BLUEs_HS1.csv"
blue_h2_file   <- "~/R/World_Veg_Project/anna_BLUEs/anna_BLUES_rebuild/BLUEs_HS2.csv"
outfile        <- "~/R/World_Veg_Project/anna_BLUEs/anna_BLUES_rebuild/STI_from_BLUEs_all_73_traits.csv"

id_col      <- "accession"
drop_cols   <- c("plantno.")

# ── load BLUEs ────────────────────────────────────────────────────────────────
blues_ctrl <- read.csv(blue_ctrl_file, stringsAsFactors = FALSE) %>%
  mutate(across(all_of(id_col), ~ trimws(as.character(.))))
blues_h1   <- read.csv(blue_h1_file,   stringsAsFactors = FALSE) %>%
  mutate(across(all_of(id_col), ~ trimws(as.character(.))))
blues_h2   <- read.csv(blue_h2_file,   stringsAsFactors = FALSE) %>%
  mutate(across(all_of(id_col), ~ trimws(as.character(.))))

# ── identify common traits ────────────────────────────────────────────────────
traits_common <- Reduce(intersect, list(
  setdiff(colnames(blues_ctrl), c(id_col, drop_cols)),
  setdiff(colnames(blues_h1),   c(id_col, drop_cols)),
  setdiff(colnames(blues_h2),   c(id_col, drop_cols))
))
traits_common <- sort(traits_common)
message("Traits to process: ", length(traits_common))

# ── STI function (same as your GWAS code) ────────────────────────────────────
make_sti_from_blues_one_trait <- function(trait, blues_ctrl, blues_h1, blues_h2,
                                          id_col = "accession") {
  dC  <- blues_ctrl %>% select(all_of(id_col), C  = all_of(trait))
  dS1 <- blues_h1   %>% select(all_of(id_col), S1 = all_of(trait))
  dS2 <- blues_h2   %>% select(all_of(id_col), S2 = all_of(trait))
  
  d <- dC %>%
    inner_join(dS1, by = id_col) %>%
    inner_join(dS2, by = id_col) %>%
    mutate(across(c(C, S1, S2), ~ suppressWarnings(as.numeric(.)))) %>%
    filter(!is.na(C), !is.na(S1), !is.na(S2))
  
  if (nrow(d) == 0) return(NULL)
  
  Xc_bar <- mean(d$C, na.rm = TRUE)
  
  d %>%
    mutate(STI = (sqrt(S1 * S2) * C) / (Xc_bar^2)) %>%
    filter(is.finite(STI)) %>%
    select(accession = all_of(id_col), STI) %>%
    rename(!!trait := STI)
}

# ── loop over traits and collect ──────────────────────────────────────────────
sti_list <- list()

for (tr in traits_common) {
  message("Computing STI for: ", tr)
  result <- make_sti_from_blues_one_trait(tr, blues_ctrl, blues_h1, blues_h2, id_col)
  if (!is.null(result)) sti_list[[tr]] <- result
}

# ── join all traits into wide matrix ─────────────────────────────────────────
sti_wide <- reduce(sti_list, full_join, by = "accession")

message("STI matrix: ", nrow(sti_wide), " lines x ", ncol(sti_wide) - 1, " traits")
head(sti_wide[, 1:5])

# ── save ──────────────────────────────────────────────────────────────────────
write.csv(sti_wide, outfile, row.names = FALSE)
###################    STI_from_BLUEs_all_73_traits.csv final file ################### 


