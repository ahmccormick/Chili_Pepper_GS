####################################################################################################
# POST GWAS TABLES - BLUE-adjusted STI GWAS results
# Table S5: All significant SNPs (FDR < 0.05) per trait per model + counts
# Table S6: Pleiotropic SNPs significant in >= 3 models within a trait
#           (for pollen: >= 2 models; for yield: >= 1 model since only 2/1 traits)
# Table S7: Multi-trait SNPs significant across >= 2 traits AND >= 3 models
#           (with special handling for pollen and yield categories)
####################################################################################################

library(tidyverse)

#############################################
# 0) File paths -- EDIT THESE
#############################################

model_dirs <- c(
  BLINK   = "~/R/WorldVeg_Capsicum/STI_GWAS/BLINK/BLINK",
  MLMM    = "~/R/WorldVeg_Capsicum/STI_GWAS/MLMM/MLMM",
  MLM     = "~/R/WorldVeg_Capsicum/STI_GWAS/MLM_PC_K/MLM_PC_K",
  FarmCPU = "~/R/WorldVeg_Capsicum/STI_GWAS/FarmCPU/FarmCPU"
)

outdir <- "~/R/WorldVeg_Capsicum/STI_GWAS/GWAS_tables/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#############################################
# 1) Category lookup
#############################################

trait_categories <- tribble(
  ~Trait,                        ~Category,
  "biomassfinalplant",           "Growth",
  "biomassfinalplot",            "Growth",
  "DASanthesis",                 "Phenology",
  "DASmaturity",                 "Phenology",
  "DATanthesis",                 "Phenology",
  "DATmaturity",                 "Phenology",
  "fruitlength",                 "Yield components",
  "fruitno",                     "Yield components",
  "fruitshapeindex",             "Yield components",
  "fruitweight",                 "Yield components",
  "fruitwidth",                  "Yield components",
  "green",                       "Canopy spectral",
  "greenb0",                     "Canopy spectral",
  "greenb1",                     "Canopy spectral",
  "greenb2",                     "Canopy spectral",
  "greenb3",                     "Canopy spectral",
  "growthrateplant",             "Growth",
  "growthrateplot",              "Growth",
  "heighmaxfinal",               "Growth",
  "height",                      "Growth",
  "heightfinal",                 "Growth",
  "heightmax",                   "Growth",
  "huebin0",                     "Canopy spectral",
  "huebin1",                     "Canopy spectral",
  "huebin5",                     "Canopy spectral",
  "huedaily",                    "Canopy spectral",
  "huemean",                     "Canopy spectral",
  "leafanglebright",             "Physiology",
  "leafanglecold",               "Physiology",
  "leafanglecoldminushot",       "Physiology",
  "leafangledaily",              "Physiology",
  "leafangleHot",                "Physiology",
  "leafanglelow",                "Physiology",
  "leafanglelowminusbright",     "Physiology",
  "leafanglemean",               "Physiology",
  "leafanglemidnight",           "Physiology",
  "leafanglemidnightminusnoon",  "Physiology",
  "leafanglenoon",               "Physiology",
  "leafareadaily",               "Leaf area",
  "leafareadailyplant",          "Leaf area",
  "leafareafinalplant",          "Leaf area",
  "leafareafinalplot",           "Leaf area",
  "leafareafinalprojectedplant", "Leaf area",
  "leafareafinalprojectedplot",  "Leaf area",
  "leafareaindex",               "Leaf area",
  "leafareaindexfinal",          "Leaf area",
  "leafareaprojecteddaily",      "Leaf area",
  "leafareaprojectedplant",      "Leaf area",
  "leafinclinationbright",       "Physiology",
  "leafinclinationcold",         "Physiology",
  "leafinclinationdark",         "Physiology",
  "leafinclinationhot",          "Physiology",
  "leafinclinationmean",         "Physiology",
  "leafinclinationmidnight",     "Physiology",
  "leafinclinationnoon",         "Physiology",
  "leafminushobotemp",           "Physiology",
  "leaftemp",                    "Physiology",
  "lightpendaily",               "Physiology",
  "lightpenetrationfinal",       "Physiology",
  "lightpenmeanplant",           "Physiology",
  "lightpenmeanplot",            "Physiology",
  "ndvibin0",                    "Canopy spectral",
  "ndvibin1",                    "Canopy spectral",
  "ndvibin2",                    "Canopy spectral",
  "ndvibin3",                    "Canopy spectral",
  "ndvibin4",                    "Canopy spectral",
  "ndvibin5",                    "Canopy spectral",
  "ndvimean",                    "Canopy spectral",
  "npcimean",                    "Canopy spectral",
  "pollenactivity",              "Pollen",
  "pollenconcentration",         "Pollen",
  "psrimean",                    "Canopy spectral",
  "yield",                       "Yield"
)

# categories with limited traits -- need lower model thresholds
limited_categories <- c("Pollen", "Yield")  # 2 and 1 trait respectively

#############################################
# 2) Find GWAS result files
#############################################

find_gwas_files <- function(base_dir, model_name) {
  files <- list.files(
    base_dir,
    pattern = paste0("^GAPIT\\.Association\\.GWAS_Results\\.", model_name, ".*\\(NYC\\)\\.csv$"),
    recursive = TRUE,
    full.names = TRUE
  )
  tibble(
    Model = model_name,
    File  = files,
    Trait = basename(dirname(files))
  )
}

file_table <- bind_rows(
  find_gwas_files(model_dirs["BLINK"],   "BLINK"),
  find_gwas_files(model_dirs["MLMM"],    "MLMM"),
  find_gwas_files(model_dirs["MLM"],     "MLM"),
  find_gwas_files(model_dirs["FarmCPU"], "FarmCPU")
)

message("Total GWAS files found: ", nrow(file_table))
print(file_table %>% count(Model))

#############################################
# 3) Read and process GWAS files
#############################################

process_gwas_file <- function(file, trait, model) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  p_col   <- grep("^P.value$|^P_value$|^P\\.value$", colnames(df), value = TRUE)[1]
  snp_col <- grep("^SNP$", colnames(df), value = TRUE)[1]
  chr_col <- grep("^Chr$|^Chromosome$", colnames(df), value = TRUE)[1]
  pos_col <- grep("^Pos$|^Position$", colnames(df), value = TRUE)[1]
  
  if (is.na(p_col) || is.na(snp_col)) {
    warning("Skipping (missing SNP or P.value): ", file)
    return(NULL)
  }
  
  out <- tibble(
    Trait   = trait,
    Model   = model,
    SNP     = df[[snp_col]],
    P.value = as.numeric(df[[p_col]])
  )
  
  if (!is.na(chr_col)) out$Chr <- as.numeric(df[[chr_col]])
  if (!is.na(pos_col)) out$Pos <- as.numeric(df[[pos_col]])
  
  out %>%
    filter(!is.na(SNP), !is.na(P.value)) %>%
    mutate(
      FDR  = p.adjust(P.value, method = "BH"),
      logP = -log10(P.value)
    )
}

#############################################
# 4) Loop through files
#############################################

sig_list   <- list()
count_list <- list()

for (i in seq_len(nrow(file_table))) {
  message("Processing ", i, "/", nrow(file_table), ": ",
          file_table$Model[i], " | ", file_table$Trait[i])
  
  gwas <- process_gwas_file(
    file  = file_table$File[i],
    trait = file_table$Trait[i],
    model = file_table$Model[i]
  )
  
  if (is.null(gwas)) next
  
  count_list[[length(count_list) + 1]] <- tibble(
    Trait                  = file_table$Trait[i],
    Model                  = file_table$Model[i],
    Significant_SNPs_FDR05 = sum(gwas$FDR < 0.05, na.rm = TRUE)
  )
  
  sig_snps <- gwas %>% filter(FDR < 0.05)
  if (nrow(sig_snps) > 0) sig_list[[length(sig_list) + 1]] <- sig_snps
  
  rm(gwas, sig_snps)
  gc()
}

sig_counts <- bind_rows(count_list) %>% arrange(Trait, Model)
sig_all    <- bind_rows(sig_list)
message("Total significant SNPs: ", nrow(sig_all))

#############################################
#############################################
# 5) Add category to sig_all
#############################################

sig_all <- sig_all %>%
  left_join(trait_categories, by = "Trait")

# check for unmatched traits
unmatched <- sig_all %>% filter(is.na(Category)) %>% distinct(Trait)
if (nrow(unmatched) > 0) {
  message("WARNING: Traits with no category match:")
  print(unmatched)
}

# S5c: significant SNPs passing FDR in >= 2 models
sig_2plus_models <- sig_all %>%
  distinct(Trait, Model, SNP, .keep_all = TRUE) %>%
  group_by(Trait, SNP) %>%
  summarise(
    n_models  = n_distinct(Model),
    Models    = paste(sort(unique(Model)), collapse = ", "),
    Category  = first(Category),
    Chr       = first(Chr),
    Pos       = first(Pos),
    min_P     = min(P.value, na.rm = TRUE),
    min_FDR   = min(FDR, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  filter(n_models >= 2) %>%
  relocate(Category, .after = Trait) %>%
  arrange(Category, Trait, desc(n_models), Chr, Pos)

message("SNPs significant in >= 2 models: ", nrow(sig_2plus_models))
message("Breakdown by n_models:")
print(sig_2plus_models %>% count(n_models))

write.csv(
  sig_2plus_models,
  file.path(outdir, "TableS5_significant_SNPs_2plus_models.csv"),
  row.names = FALSE
)

#############################################
# 7) TABLE S6: Pleiotropic SNPs
#    >= 3 models within a trait
#    Exception: Pollen >= 2 models, Yield >= 1 model
#############################################

# define model threshold per category
get_model_threshold <- function(category) {
  case_when(
    category == "Yield"  ~ 1L,   # only 1 trait so any model counts
    category == "Pollen" ~ 2L,   # only 2 traits so lower threshold
    TRUE                 ~ 3L    # all other categories
  )
}
#############################################
# TABLE S6: SNPs significant in >= 3 traits
# (any model, just needs to pass FDR < 0.05)
# Special handling: Pollen >= 2 traits (only 2 exist)
# Yield included (only 1 trait)
#############################################

table_s6 <- sig_all %>%
  distinct(Trait, SNP, .keep_all = TRUE) %>%
  group_by(SNP) %>%
  summarise(
    Traits_with_SNP   = n_distinct(Trait),
    Traits            = paste(sort(unique(Trait)), collapse = ", "),
    n_categories      = n_distinct(Category),
    Categories        = paste(sort(unique(Category)), collapse = ", "),
    dominant_category = names(sort(table(Category), decreasing = TRUE))[1],
    Chr               = first(Chr),
    Pos               = first(Pos),
    min_P             = min(P.value, na.rm = TRUE),
    min_FDR           = min(FDR, na.rm = TRUE),
    best_trait        = Trait[which.min(P.value)],
    .groups           = "drop"
  ) %>%
  filter(
    Traits_with_SNP >= 3 |                                           # standard: >= 3 traits
    (dominant_category == "Pollen" & Traits_with_SNP >= 2) |        # pollen: >= 2 traits
    (dominant_category == "Yield"  & Traits_with_SNP >= 1)          # yield: any significance
  ) %>%
  arrange(desc(Traits_with_SNP), Chr, Pos)

message("Table S6 SNPs: ", nrow(table_s6))
message("Breakdown by traits per SNP:")
print(table_s6 %>% count(Traits_with_SNP) %>% arrange(desc(Traits_with_SNP)))
message("Breakdown by dominant category:")
print(table_s6 %>% count(dominant_category))

write.csv(
  table_s6,
  file.path(outdir, "TableS6_multitrait_SNPs_3plus_traits.csv"),
  row.names = FALSE
)



#############################################
# 8) TABLE S7: Multi-trait SNPs
#    Significant across >= 2 traits AND >= 3 models total
#    Exception: Pollen >= 2 traits not applicable (only 2 traits),
#    include if significant in both pollen traits
#    Yield: single trait so excluded from multi-trait by definition
#############################################

table_s7 <- sig_all %>%
  distinct(Trait, Model, SNP, .keep_all = TRUE) %>%
  group_by(SNP) %>%
  summarise(
    n_traits          = n_distinct(Trait),
    n_models_total    = n_distinct(Model),
    Models            = paste(sort(unique(Model)), collapse = ", "),
    Traits            = paste(sort(unique(Trait)), collapse = ", "),
    n_categories      = n_distinct(Category),
    Categories        = paste(sort(unique(Category)), collapse = ", "),
    dominant_category = names(sort(table(Category), decreasing = TRUE))[1],
    Chr               = first(Chr),
    Pos               = first(Pos),
    min_P             = min(P.value, na.rm = TRUE),
    min_FDR           = min(FDR, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  filter(
    (n_traits >= 3 & n_models_total >= 3) |                    # standard
      (n_traits >= 2 & dominant_category == "Pollen") |          # pollen exception
      (dominant_category == "Yield" & n_models_total >= 3)       # yield exception: single trait but >= 3 models
  ) %>%
  arrange(desc(n_traits), desc(n_models_total), Chr, Pos)

write.csv(
  table_s7,
  file.path(outdir, "TableS7_multitrait_SNPs_3plus_models.csv"),
  row.names = FALSE
)

message("Table S7 saved: ", nrow(table_s7), " SNPs")






