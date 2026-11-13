########################################################################################################################
# Supplemental Figure S1
########################################################################################################################
#Figure S1. Raw climate data for the control, heat stress-1 and heat stress-2 timepoints 
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


########################################################################################################################
# Supplemental Figure S2
########################################################################################################################
#Figure S2. Population structure for the core (n=397) capsicum annuum only collection
#core population structure/training plots
#SNPRelate
#########
library(SNPRelate)
getwd()
setwd("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/")

#vcf.fn <- "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_filtered_and_LD_pruned.vcf.gz"
vcf.fn <- "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Pepper_397_annuum_LD_pruned_beagle.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "pepper.gds", method="copy.num.of.ref")
snpgdsSummary("pepper.gds")
genofile <- snpgdsOpen("pepper.gds")
set.seed(1234)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
#write.csv(tab, "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC.csv")
write.csv(tab, "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Pepper_PC.csv")


library(ggplot2)
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC.csv")
pca <- read.csv("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Pepper_PC.csv")

# Basic PCA plot: EV1 vs. EV2
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2)) +
  geom_point(size = 3, alpha = 0.8) + # Points for samples
  theme_minimal() +                   # Clean theme
  labs(
    title = "Pepper PCA 1",
    x = "EV1 (2.57 %)",
    y = "EV2 (1.74%)"
  )

# Display the plot
print(pca_plot)

###############################################################################
# METADATA

# Load necessary libraries
library(readr)    # For reading CSV files
library(readxl)   # For reading Excel files
library(dplyr)    # For data manipulation

# Load the Pepper_PC.csv file
#pepper_pc <- read_csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Pepper_PC.csv")
pepper_pc <- read_csv("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Pepper_PC.csv")

# Load the Global Capsicum Core Collection Excel file (first sheet)
#capsicum_core <- read_excel("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Global capsicum core collection.xlsx")
capsicum_core <- read_excel("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Global capsicum core collection.xlsx")


# Print column names to check structure
colnames(pepper_pc)
colnames(capsicum_core)

# Perform an inner join to match based on the sample.id column
matched_data <- pepper_pc %>%
  inner_join(capsicum_core, by = c("sample.id" = "id")) 

# Save the matched data to a new CSV file
#write_csv(matched_data, "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv")
#write_csv(matched_data, "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv") 
write_csv(matched_data, "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Matched_Pepper_Data_397.csv") 

################################################################################
#missing three

library(dplyr)
library(readr)

# Load data
pepper_pc <- read_csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Pepper_PC.csv")
capsicum_core <- read_excel("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Global capsicum core collection.xlsx")

# Perform an inner join
matched_data <- pepper_pc %>%
  inner_join(capsicum_core, by = c("sample.id" = "id"))

# Find missing sample.id values
missing_ids <- anti_join(pepper_pc, capsicum_core, by = c("sample.id" = "id"))

# Print missing sample IDs
print(missing_ids$sample.id)

##############################################################################
library(ggplot2)
library(dplyr)

# Load data
#Pannel A
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_EDITS.csv")
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv")
#pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Matched_Pepper_Data_423.csv")
pca <- read.csv("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/Matched_Pepper_Data_397.csv")


# Count occurrences of each species
species_counts <- pca %>%
  group_by(species) %>%
  summarise(n = n())

# Create a named vector for new legend labels
species_labels <- setNames(
  paste0(species_counts$species, " (n=", species_counts$n, ")"),
  species_counts$species
)

# Basic PCA plot with updated legend
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2, color = species)) +
  geom_point(size = 3, alpha = 0.8) +  # Points for samples
  theme_minimal() +  # Clean theme
  labs(
    title = "",
    x = "PC1 (2.57%)",
    y = "PC2 (1.74%)",
    color = "Species"  # Legend title
  ) +
  scale_color_viridis_d(labels = species_labels)  # Updated legend labels

# Display the plot
print(pca_plot)


ggsave(
  "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/FIGURES/S2A.pdf",
  plot = last_plot(),
  width = 8,
  height = 6,
  units = "in"
)


#####################################################################################################################
#Pannel B
#####################################################################################################################
#colour by hcpc group in pca

library(dplyr)
library(ggplot2)

pca_plot <- read.csv(
  "Pepper_397_BEAGLE_metadata_NEW_HCPC_clusters.csv",
  check.names = FALSE
) %>%
  mutate(cluster = factor(cluster))

stopifnot(
  nrow(pca_plot) == 397,
  !anyDuplicated(pca_plot$sample.id),
  !anyNA(pca_plot$cluster)
)

cluster_counts <- pca_plot %>% count(cluster)
print(cluster_counts)

cluster_labels <- setNames(
  paste0("Cluster ", cluster_counts$cluster,
         " (n=", cluster_counts$n, ")"),
  cluster_counts$cluster
)

# Only draw ellipses for groups with at least four samples
ellipse_data <- pca_plot %>%
  group_by(cluster) %>%
  filter(n() >= 4) %>%
  ungroup()

ggplot(pca_plot, aes(EV1, EV2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(
    data = ellipse_data,
    type = "norm",
    linewidth = 1
  ) +
  scale_color_viridis_d(labels = cluster_labels) +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc.percent[1]),
    y = sprintf("PC2 (%.2f%%)", pc.percent[2]),
    color = "HCPC Cluster"
  ) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(
  "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/FIGURES/S2B_Pepper_397_BEAGLE_PCA_NEW_HCPC_clusters.pdf",
  plot = last_plot(),
  width = 8,
  height = 6,
  units = "in"
)

#####################################################################################################################
#CORE VISUALISATION
#Pannel C
#####################################################################################################################
library(ggplot2)
library(dplyr)


pca_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE"

training_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core/GS_all_phenotyped"

output_folder <- file.path(pca_folder, "FIGURES")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  file.path(
    pca_folder,
    "Pepper_397_annuum_BEAGLE_PCA_metadata.csv"
  ),
  check.names = FALSE
)

metadata$sample.id <- trimws(as.character(metadata$sample.id))

stopifnot(
  nrow(metadata) == 397,
  !anyDuplicated(metadata$sample.id),
  all(c("sample.id", "EV1", "EV2") %in% names(metadata))
)

environments <- c(
  "Control",
  "Heat_Stress_1",
  "Heat_Stress_2"
)

training_records <- bind_rows(
  lapply(environments, function(environment) {
    
    read.csv(
      file.path(
        training_folder,
        paste0("Training_IDs_by_trait_", environment, ".csv")
      ),
      check.names = FALSE
    )
  })
)

# Accessions used for at least one trait/environment
training_ids <- unique(
  trimws(as.character(training_records$Sample_ID))
)

metadata$Partition <- ifelse(
  metadata$sample.id %in% training_ids,
  "Training",
  "Prediction only"
)

metadata$Partition <- factor(
  metadata$Partition,
  levels = c("Prediction only", "Training")
)

counts <- table(metadata$Partition)
print(counts)

####################################################################
# 3. Plot
####################################################################

# Draw training accessions last
metadata <- metadata %>%
  arrange(Partition)

# Use variance percentages from the existing PCA, if available
x_label <- "PC1"
y_label <- "PC2"

if (
  exists("pc.percent") &&
  is.numeric(pc.percent) &&
  length(pc.percent) >= 2 &&
  all(is.finite(pc.percent[1:2]))
) {
  x_label <- sprintf("PC1 (%.2f%%)", pc.percent[1])
  y_label <- sprintf("PC2 (%.2f%%)", pc.percent[2])
}

p <- ggplot(
  metadata,
  aes(x = EV1, y = EV2, color = Partition)
) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(
    values = c(
      "Prediction only" = "grey80",
      "Training" = "red3"
    ),
    breaks = c("Prediction only", "Training"),
    labels = c(
      paste0(
        "Prediction only (n=", counts[["Prediction only"]], ")"
      ),
      paste0(
        "Training (n=", counts[["Training"]], ")"
      )
    ),
    drop = FALSE
  ) +
  labs(
    x = "PC1 (2.57%)",
    y = "PC2 (1.74%)",
    color = "GS use",
    caption = paste(
      "",
      ""
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(p)

ggsave(
  filename = file.path(output_folder, "S2C.pdf"),
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)

#####################################################################################################################
#pannel D
#####################################################################################################################
#  PCA highlighting yield rank-response groups
library(ggplot2)
library(dplyr)

####################################################################
# 1. Paths
####################################################################

base_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1"

pca_folder <- file.path(
  base_folder,
  "raw_vcf", "partitioned_to_annuum", "BEAGLE"
)

rank_folder <- file.path(
  base_folder,
  "GS_core", "GS_all_phenotyped", "S13"
)

figure_folder <- file.path(pca_folder, "FIGURES")
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

####################################################################
# 2. Load PCA metadata and selected accession lists
####################################################################

metadata <- read.csv(
  file.path(
    pca_folder,
    "Pepper_397_annuum_BEAGLE_PCA_metadata.csv"
  ),
  check.names = FALSE
)

heat_accessions <- read.csv(
  file.path(rank_folder, "Heat_favoured_accessions_core397.csv"),
  check.names = FALSE
)

control_accessions <- read.csv(
  file.path(rank_folder, "Control_favoured_accessions_core397.csv"),
  check.names = FALSE
)

metadata$sample.id <- trimws(as.character(metadata$sample.id))

heat_ids <- unique(trimws(as.character(heat_accessions$Sample_ID)))
control_ids <- unique(trimws(as.character(control_accessions$Sample_ID)))

stopifnot(
  nrow(metadata) == 397,
  !anyDuplicated(metadata$sample.id),
  length(intersect(heat_ids, control_ids)) == 0,
  all(c(heat_ids, control_ids) %in% metadata$sample.id)
)

####################################################################
# 3. Assign groups
####################################################################

metadata <- metadata %>%
  mutate(
    Response_group = case_when(
      sample.id %in% heat_ids ~ "Heat-favoured",
      sample.id %in% control_ids ~ "Control-favoured",
      TRUE ~ "Other accessions"
    ),
    Response_group = factor(
      Response_group,
      levels = c(
        "Other accessions",
        "Control-favoured",
        "Heat-favoured"
      )
    )
  ) %>%
  arrange(Response_group)

counts <- table(metadata$Response_group)
print(counts)

####################################################################
# 4. Plot
# Grey accessions first; selected groups drawn on top
####################################################################

p_D <- ggplot(metadata, aes(x = EV1, y = EV2)) +
  geom_point(
    data = filter(metadata, Response_group == "Other accessions"),
    aes(color = Response_group),
    size = 2,
    alpha = 0.5
  ) +
  geom_point(
    data = filter(metadata, Response_group != "Other accessions"),
    aes(color = Response_group),
    size = 2.8,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "Other accessions" = "grey85",
      "Control-favoured" = "royalblue4",
      "Heat-favoured" = "firebrick3"
    ),
    breaks = c(
      "Control-favoured",
      "Heat-favoured",
      "Other accessions"
    ),
    labels = c(
      paste0(
        "Control-favoured (n=", counts[["Control-favoured"]], ")"
      ),
      paste0(
        "Heat-favoured (n=", counts[["Heat-favoured"]], ")"
      ),
      paste0(
        "Other accessions (n=", counts[["Other accessions"]], ")"
      )
    ),
    drop = FALSE
  ) +
  labs(
    title = "",
    x = "PC1 (2.57%)",
    y = "PC2 (1.74%)",
    color = "Yield rank response"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(p_D)

####################################################################
# 5. Save
####################################################################

ggsave(
  filename = file.path(figure_folder, "S2D.pdf"),
  plot = p_D,
  width = 8,
  height = 6,
  units = "in"
)

write.csv(
  metadata,
  file.path(
    pca_folder,
    "Pepper_397_PCA_metadata_yield_response_groups.csv"
  ),
  row.names = FALSE
)

#####################################################################################################################
#pannel E
#####################################################################################################################
# Pepper Dendrogram with training/test Bar

####################################################################
# Panel E: annotate the existing HCPC tree with actual GS use
# No PCA or clustering refitting required
####################################################################

library(dplyr)
library(ggplot2)
library(factoextra)

####################################################################
# 1. Paths
####################################################################

base_folder <- paste0(
  "/Users/annamccormick/R/World_Veg_Project/",
  "Reviewer_comments_round1"
)

training_folder <- file.path(
  base_folder,
  "GS_core",
  "GS_all_phenotyped"
)

figure_folder <- file.path(
  base_folder,
  "raw_vcf",
  "partitioned_to_annuum",
  "BEAGLE",
  "FIGURES"
)

dir.create(
  figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

####################################################################
# 2. Read actual training IDs, as in panel C
####################################################################

environments <- c(
  "Control",
  "Heat_Stress_1",
  "Heat_Stress_2"
)

training_records <- dplyr::bind_rows(
  lapply(environments, function(environment) {
    
    filename <- file.path(
      training_folder,
      paste0(
        "Training_IDs_by_trait_",
        environment,
        ".csv"
      )
    )
    
    stopifnot(file.exists(filename))
    
    records <- read.csv(
      filename,
      check.names = FALSE
    )
    
    stopifnot("Sample_ID" %in% names(records))
    
    records
  })
)

# Union: used for at least one trait in at least one environment
training_ids <- unique(
  trimws(as.character(training_records$Sample_ID))
)

training_ids <- training_ids[
  !is.na(training_ids) & nzchar(training_ids)
]

####################################################################

library(FactoMineR)
library(dplyr)

pca_folder <- file.path(
  base_folder,
  "raw_vcf",
  "partitioned_to_annuum",
  "BEAGLE"
)

pepper_pca <- read.csv(
  file.path(
    pca_folder,
    "Pepper_397_BEAGLE_metadata_NEW_HCPC_clusters.csv"
  ),
  check.names = FALSE
)

pepper_pca$sample.id <- trimws(
  as.character(pepper_pca$sample.id)
)

pca_data <- pepper_pca %>%
  dplyr::select(EV1, EV2, EV3, EV4)

rownames(pca_data) <- pepper_pca$sample.id

stopifnot(
  nrow(pca_data) == 397L,
  !anyDuplicated(pepper_pca$sample.id),
  all(vapply(pca_data, is.numeric, logical(1))),
  all(is.finite(as.matrix(pca_data)))
)

res.pca <- FactoMineR::PCA(
  pca_data,
  scale.unit = FALSE,
  ncp = 4,
  graph = FALSE
)

set.seed(1234)

res.hcpc <- FactoMineR::HCPC(
  res.pca,
  nb.clust = 5,
  consol = TRUE,
  graph = FALSE
)
# 3. Match training status to the existing tree's leaf order
####################################################################

stopifnot(exists("res.hcpc"))

tree <- res.hcpc$call$t$tree

leaf_ids <- trimws(
  as.character(tree$labels[tree$order])
)

stopifnot(
  length(leaf_ids) == 397L,
  !anyDuplicated(leaf_ids),
  length(training_ids) == 285L,
  all(training_ids %in% leaf_ids)
)

color_bar_df <- data.frame(
  x = seq_along(leaf_ids),
  Sample_ID = leaf_ids,
  Type = factor(
    ifelse(
      leaf_ids %in% training_ids,
      "Training",
      "Prediction only"
    ),
    levels = c("Prediction only", "Training")
  )
)

counts <- table(color_bar_df$Type)

print(counts)

stopifnot(
  unname(counts["Training"]) == 285L,
  unname(counts["Prediction only"]) == 112L
)

####################################################################
# 4. Define the annotation bar position
####################################################################

tree_height <- max(tree$height)

stopifnot(
  is.finite(tree_height),
  tree_height > 0
)

bar_position_y <- -0.04 * tree_height
bar_height <- 0.02 * tree_height

####################################################################
# 5. Draw the dendrogram and corrected annotation bar
####################################################################
f <- factoextra::fviz_dend(
  res.hcpc,
  cex = 0.2,
  rect = FALSE,
  show_labels = FALSE
) +
  geom_tile(
    data = color_bar_df,
    aes(
      x = x,
      y = y,
      fill = Type
    ),
    inherit.aes = FALSE,
    width = 1,
    height = bar_height
  ) +
  scale_fill_manual(
    values = c(
      "Prediction only" = "grey80",
      "Training" = "red3"
    ),
    breaks = c("Prediction only", "Training"),
    labels = c(
      paste0(
        "Prediction only (n=",
        counts[["Prediction only"]],
        ")"
      ),
      paste0(
        "Training (n=",
        counts[["Training"]],
        ")"
      )
    ),
    name = "GS use",
    drop = FALSE
  ) +
  expand_limits(
    y = bar_position_y - bar_height
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(
      t = 1, r = 1, b = 1, l = 1,
      unit = "cm"
    )
  )

color_bar_df$y <- rep(
  bar_position_y,
  nrow(color_bar_df)
)

# Update the annotation layer with the corrected data
f$layers[[3]]$data <- color_bar_df

print(f)

ggsave(
  filename = file.path(
    figure_folder,
    "S2E_HCPC_actual_training285.pdf"
  ),
  plot = f,
  width = 12,
  height = 6,
  units = "in"
)

write.csv(
  pepper_pca,
  "Pepper_397_BEAGLE_metadata_NEW_HCPC_clusters.csv",
  row.names = FALSE
)


########################################################################################################################
# Supplemental Figure S3
########################################################################################################################
#Kinship matrix heatmap 

library(SNPRelate)
library(pheatmap)

setwd("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/")
genofile <- snpgdsOpen("pepper.gds")

# Genomic relationship matrix
grm_result <- snpgdsGRM(
  genofile,
  method = "EIGMIX",
  autosome.only = FALSE,
  remove.monosnp = TRUE,
  num.thread = 2,
  with.id = TRUE
)

snpgdsClose(genofile)

# Convert relationship matrix to genomic kinship
kinship_matrix <- grm_result$grm / 2
sample_ids <- grm_result$sample.id
rownames(kinship_matrix) <- colnames(kinship_matrix) <- sample_ids

# Match metadata
metadata <- read.csv(
  "Pepper_397_annuum_BEAGLE_PCA_metadata_80_20.csv",
  check.names = FALSE
)

stopifnot(
  length(sample_ids) == 397,
  !anyDuplicated(metadata$sample.id),
  all(sample_ids %in% metadata$sample.id)
)

metadata <- metadata[match(sample_ids, metadata$sample.id), ]

annotation <- data.frame(
  Population = metadata$species,
  row.names = sample_ids
)

stopifnot(!anyNA(annotation$Population))

groups <- unique(annotation$Population)
ann_colors <- list(
  Population = setNames(rainbow(length(groups)), groups)
)

# Symmetric colour scale centred on zero
limit <- max(abs(kinship_matrix))

heatmap_plot <- pheatmap(
  kinship_matrix,
  scale = "none",
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Genomic Kinship Matrix — C. annuum (n=397)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE
)

pdf("Pepper_397_BEAGLE_kinship_heatmap.pdf", width = 12, height = 12)
print(heatmap_plot)
dev.off()

########################################################################################################################
# Supplemental Figure S4
########################################################################################################################
# fastSTRUCTURE for the core collection n=397

library(ggplot2)
library(dplyr)
library(tidyr)
library(MetBrewer)
library(readxl)

folder <- "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/fastSTRUCTURE"

# Consistent ancestry colors
palette <- met.brewer("Hiroshige", 8)
k_colors <- setNames(
  palette[c(1, 6, 3, 8, 2, 5)],
  paste0("K", 1:6)
)

plots <- list()

for (K in 2:6) {
  
  pepper <- read_excel(
    file.path(folder, paste0("n397_pepper_k", K, ".xlsx")),
    sheet = 1
  )
  
  k_columns <- paste0("K", seq_len(K))
  
  # Sort individuals by ancestry
  pepper <- pepper %>%
    arrange(across(all_of(k_columns), ~ desc(.x)))
  
  pepper$sample.id <- factor(
    pepper$sample.id,
    levels = unique(pepper$sample.id)
  )
  
  long_data <- pepper %>%
    pivot_longer(
      cols = all_of(k_columns),
      names_to = "category",
      values_to = "value"
    ) %>%
    mutate(category = factor(category, levels = k_columns))
  
  final_plot <- ggplot(
    long_data,
    aes(x = sample.id, y = value, fill = category)
  ) +
    geom_col(width = 1) +
    scale_fill_manual(
      name = "Ancestry cluster",
      values = k_colors[k_columns]
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste0("K = ", K, " (n = ", nrow(pepper), ")"),
      x = "",
      y = "Ancestry Proportion"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1, size = 1
      ),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8)
    )
  
  print(final_plot)
  plots[[paste0("K", K)]] <- final_plot
  
  ggsave(
    filename = file.path(
      folder,
      paste0("n397_pepper_fastSTRUCTURE_K", K, ".pdf")
    ),
    plot = final_plot,
    width = 12,
    height = 4
  )
}

########################################################################################################################
# Supplemental Figure S5
########################################################################################################################
# Compare HCPC solutions: K = 2–10

library(FactoMineR)
library(cluster)
library(ggplot2)

# Same coordinates used for clustering
X <- as.matrix(res.pca$ind$coord[, 1:4, drop = FALSE])
D <- dist(X)

k_values <- 2:10

hcpc_solutions <- list()
diagnostic_rows <- list()

for (k in k_values) {
  
  set.seed(1234)
  
  fit <- FactoMineR::HCPC(
    res.pca,
    nb.clust = k,
    consol = TRUE,
    graph = FALSE
  )
  
  # Align assignments with coordinate row order
  idx <- match(rownames(X), rownames(fit$data.clust))
  stopifnot(!anyNA(idx))
  
  groups <- as.integer(
    factor(fit$data.clust$clust[idx])
  )
  
  # Within-cluster sum of squared Euclidean distances
  wss <- sum(vapply(
    split(seq_len(nrow(X)), groups),
    function(rows) {
      Xi <- X[rows, , drop = FALSE]
      sum(sweep(Xi, 2, colMeans(Xi), "-")^2)
    },
    numeric(1)
  ))
  
  sil <- cluster::silhouette(groups, D)
  
  diagnostic_rows[[as.character(k)]] <- data.frame(
    K = k,
    Within_cluster_SS = wss,
    Mean_silhouette = mean(sil[, "sil_width"]),
    Smallest_cluster = min(table(groups))
  )
  
  hcpc_solutions[[as.character(k)]] <- fit
}

hcpc_diagnostics <- do.call(rbind, diagnostic_rows)
rownames(hcpc_diagnostics) <- NULL

print(hcpc_diagnostics)

####################################################################
# Elbow plot
####################################################################

p_elbow <- ggplot(
  hcpc_diagnostics,
  aes(K, Within_cluster_SS)
) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_continuous(breaks = k_values) +
  theme_bw() +
  labs(
    x = "Number of clusters (K)",
    y = "Within-cluster sum of squares",
    title = "HCPC elbow diagnostic"
  )

print(p_elbow)

####################################################################
# Silhouette plot
####################################################################

p_silhouette <- ggplot(
  hcpc_diagnostics,
  aes(K, Mean_silhouette)
) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_continuous(breaks = k_values) +
  theme_bw() +
  labs(
    x = "Number of clusters (K)",
    y = "Average silhouette width",
    title = "HCPC silhouette diagnostic"
  )

print(p_silhouette)

cat(
  "Highest average silhouette among evaluated solutions: K =",
  hcpc_diagnostics$K[
    which.max(hcpc_diagnostics$Mean_silhouette)
  ],
  "\n"
)


########################################################################################################################
# Supplemental Figure S6
########################################################################################################################
#PCA for the global

library(SNPRelate)
library(ggplot2)

setwd("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/")

vcf.fn <- "pepper_10000_annuum_beagle.vcf.gz"
gds.fn <- "Pepper_global_7843_BEAGLE.gds"

# Convert once
if (!file.exists(gds.fn)) {
  snpgdsVCF2GDS(
    vcf.fn,
    gds.fn,
    method = "copy.num.of.ref"
  )
}

snpgdsSummary(gds.fn)
global_gds <- snpgdsOpen(gds.fn)

# LD pruning
set.seed(1234)
global_snpset <- snpgdsLDpruning(
  global_gds,
  ld.threshold = 0.2,
  autosome.only = FALSE
)

global_snp_ids <- unlist(global_snpset, use.names = FALSE)

# PCA
global_pca <- snpgdsPCA(
  global_gds,
  snp.id = global_snp_ids,
  num.thread = 2,
  autosome.only = FALSE
)

snpgdsClose(global_gds)

global_pc_percent <- global_pca$varprop * 100
head(round(global_pc_percent, 2))

global_tab <- data.frame(
  sample.id = global_pca$sample.id,
  EV1 = global_pca$eigenvect[, 1],
  EV2 = global_pca$eigenvect[, 2],
  EV3 = global_pca$eigenvect[, 3],
  EV4 = global_pca$eigenvect[, 4]
)

stopifnot(nrow(global_tab) == 7843)

write.csv(
  global_tab,
  "Pepper_global_7843_BEAGLE_PC.csv",
  row.names = FALSE
)

global_plot <- ggplot(global_tab, aes(EV1, EV2)) +
  geom_point(size = 1, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Global C. annuum collection (n=7,843)",
    x = sprintf("PC1 (%.2f%%)", global_pc_percent[1]),
    y = sprintf("PC2 (%.2f%%)", global_pc_percent[2])
  )

print(global_plot)

ggsave(
  "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/FIGURES/S5_Pepper_global_7843_BEAGLE_PCA.pdf",
  plot = global_plot,
  width = 8,
  height = 6
)



########################################################################################################################
# Supplemental Figure S7
########################################################################################################################
#kinship for global

library(SNPRelate)
library(pheatmap)

setwd("/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/raw_vcf/partitioned_to_annuum/BEAGLE/")

global_gds <- snpgdsOpen("Pepper_global_7843_BEAGLE.gds")

# Same genomic relationship estimator as the core
global_grm <- snpgdsGRM(
  global_gds,
  method = "EIGMIX",
  autosome.only = FALSE,
  remove.monosnp = TRUE,
  num.thread = 2,
  with.id = TRUE
)

snpgdsClose(global_gds)

global_kinship <- global_grm$grm / 2
global_ids <- global_grm$sample.id

rownames(global_kinship) <- colnames(global_kinship) <- global_ids
stopifnot(length(global_ids) == 7843)

# Save the matrix for reuse
saveRDS(global_kinship, "Pepper_global_7843_BEAGLE_kinship.rds")

# All samples are C. annuum
global_annotation <- data.frame(
  Population = rep("Capsicum annuum", length(global_ids)),
  row.names = global_ids
)

global_ann_colors <- list(
  Population = c("Capsicum annuum" = "red")
)

# Compute clustering once and use it for both axes
global_hc <- hclust(
  dist(global_kinship, method = "euclidean"),
  method = "complete"
)

global_heatmap <- pheatmap(
  global_kinship,
  scale = "none",
  annotation_col = global_annotation,
  annotation_row = global_annotation,
  annotation_colors = global_ann_colors,
  cluster_rows = global_hc,
  cluster_cols = global_hc,
  main = "Genomic Kinship Matrix — Global C. annuum (n=7,843)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  border_color = NA,
  silent = TRUE
)

# Rasterise only the matrix; retain other plot elements
matrix_index <- which(global_heatmap$gtable$layout$name == "matrix")

heatmap_colors <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(7, "RdYlBu"))
)(100)

heatmap_breaks <- seq(
  min(global_kinship),
  max(global_kinship),
  length.out = 101
)

heatmap_raster <- as.raster(
  pheatmap:::scale_colours(
    global_kinship[global_hc$order, global_hc$order],
    col = heatmap_colors,
    breaks = heatmap_breaks,
    na_col = "grey"
  )
)

global_heatmap$gtable$grobs[[matrix_index]] <-
  grid::rasterGrob(heatmap_raster, interpolate = FALSE)

# Display
grid::grid.newpage()
grid::grid.draw(global_heatmap$gtable)

# Save
pdf(
  "Pepper_global_7843_BEAGLE_kinship_heatmap.pdf",
  width = 12,
  height = 12
)
grid::grid.newpage()
grid::grid.draw(global_heatmap$gtable)
dev.off()

########################################################################################################################
# Supplemental Figure S8
########################################################################################################################
# fastSTRUCTURE for the global C. annuum collection n = 7,843

library(ggplot2)
library(dplyr)
library(tidyr)
library(MetBrewer)
library(readxl)

folder <- "/Users/annamccormick/r/World_Veg_Project/Reviewer_comments_round1/fastSTRUCTURE"

palette <- met.brewer("Hiroshige", 8)

k_colors <- setNames(
  palette[c(1, 6, 3, 8, 2, 5)],
  paste0("K", 1:6)
)

plots <- list()

for (K in 2:6) {
  
  pepper <- read_excel(
    file.path(folder, paste0("n7843_pepper_k", K, ".xlsx")),
    sheet = 1
  )
  
  k_columns <- paste0("K", seq_len(K))
  
  # Ensure ancestry proportions are numeric
  pepper <- pepper %>%
    mutate(across(all_of(k_columns), ~ as.numeric(.x)))
  
  # Check that proportions are valid and sum to approximately 1
  ancestry <- as.matrix(pepper[k_columns])
  
  if (any(!is.finite(ancestry)) ||
      any(ancestry < 0 | ancestry > 1) ||
      any(abs(rowSums(ancestry) - 1) > 0.01)) {
    stop("Check ancestry proportions in the K = ", K, " Excel file.")
  }
  
  # Sort by ancestry and assign one plotting position per row
  pepper <- pepper %>%
    arrange(across(all_of(k_columns), ~ desc(.x))) %>%
    mutate(individual = row_number())
  
  long_data <- pepper %>%
    pivot_longer(
      cols = all_of(k_columns),
      names_to = "category",
      values_to = "value"
    ) %>%
    mutate(category = factor(category, levels = k_columns))
  
  final_plot <- ggplot(
    long_data,
    aes(x = individual, y = value, fill = category)
  ) +
    geom_col(width = 1, colour = NA) +
    scale_fill_manual(
      name = "Ancestry cluster",
      values = k_colors[k_columns]
    ) +
    scale_x_continuous(
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("K = ", K, " (n = ", nrow(pepper), ")"),
      x = "",
      y = "Ancestry Proportion"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8)
    )
  
  print(final_plot)
  plots[[paste0("K", K)]] <- final_plot
  
  ggsave(
    filename = file.path(
      folder,
      paste0("n7843_pepper_fastSTRUCTURE_K", K, ".pdf")
    ),
    plot = final_plot,
    width = 12,
    height = 4
  )
}

#####################################################################################################################
# Supplemental Figure S9
#####################################################################################################################
# UpSet plots for 46 high-confidence SNPs
# Input: Table_S7_new.xlsx
####################################################################
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(UpSetR)

####################################################################
# 1. Select the new Table S7
####################################################################

message("Select Table_S7_new.xlsx")

input_file <- file.choose()

output_folder <- file.path(
  dirname(input_file),
  "UpSet_46_high_confidence_SNPs"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

####################################################################
# 2. Read and check the table
####################################################################

snp_table <- readxl::read_xlsx(
  input_file,
  sheet = 1
)

names(snp_table) <- trimws(names(snp_table))

required_columns <- c(
  "SNP",
  "Traits",
  "Models",
  "n_traits",
  "n_models_total"
)

missing_columns <- setdiff(
  required_columns,
  names(snp_table)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Remove empty spreadsheet rows
snp_table <- snp_table %>%
  dplyr::filter(
    !is.na(SNP),
    nzchar(trimws(as.character(SNP)))
  ) %>%
  dplyr::mutate(
    SNP = trimws(as.character(SNP)),
    Traits = trimws(as.character(Traits)),
    Models = trimws(as.character(Models)),
    n_traits = as.integer(n_traits),
    n_models_total = as.integer(n_models_total)
  )

if (nrow(snp_table) != 46L) {
  stop(
    "Expected 46 SNP rows, but found ",
    nrow(snp_table),
    ". Check the selected file."
  )
}

if (anyDuplicated(snp_table$SNP)) {
  stop("Duplicated SNP identifiers were found.")
}

if (
  anyNA(snp_table$Traits) ||
  anyNA(snp_table$Models) ||
  any(!nzchar(snp_table$Traits)) ||
  any(!nzchar(snp_table$Models))
) {
  stop("Some SNPs have missing or empty Traits or Models.")
}

cat("Unique SNPs loaded:", nrow(snp_table), "\n")

####################################################################
# 3. Build trait memberships
####################################################################

trait_membership <- snp_table %>%
  dplyr::transmute(
    SNP,
    Trait = Traits
  ) %>%
  tidyr::separate_rows(
    Trait,
    sep = ","
  ) %>%
  dplyr::mutate(
    Trait = trimws(Trait)
  ) %>%
  dplyr::filter(
    !is.na(Trait),
    nzchar(Trait)
  ) %>%
  dplyr::distinct(SNP, Trait)

####################################################################
# 4. Build model memberships
####################################################################

model_membership <- snp_table %>%
  dplyr::select(SNP, Models) %>%
  tidyr::separate_rows(
    Models,
    sep = ","
  ) %>%
  dplyr::transmute(
    SNP,
    model_key = tolower(trimws(Models))
  )

model_lookup <- c(
  "blink" = "BLINK",
  "farmcpu" = "FarmCPU",
  "mlm" = "MLM",
  "mlmm" = "MLMM"
)

unknown_models <- setdiff(
  unique(model_membership$model_key),
  names(model_lookup)
)

if (length(unknown_models) > 0) {
  stop(
    "Unrecognised model names: ",
    paste(unknown_models, collapse = ", ")
  )
}

model_membership <- model_membership %>%
  dplyr::mutate(
    Model = unname(model_lookup[model_key])
  ) %>%
  dplyr::select(SNP, Model) %>%
  dplyr::distinct(SNP, Model)

####################################################################
# 5. Verify membership counts against Table S7
####################################################################

membership_counts <- snp_table %>%
  dplyr::select(
    SNP,
    Expected_traits = n_traits,
    Expected_models = n_models_total
  ) %>%
  dplyr::left_join(
    trait_membership %>%
      dplyr::count(SNP, name = "Observed_traits"),
    by = "SNP"
  ) %>%
  dplyr::left_join(
    model_membership %>%
      dplyr::count(SNP, name = "Observed_models"),
    by = "SNP"
  )

readr::write_csv(
  membership_counts,
  file.path(output_folder, "Membership_counts_per_SNP.csv")
)

if (
  anyNA(membership_counts) ||
  any(
    membership_counts$Expected_traits !=
    membership_counts$Observed_traits
  ) ||
  any(
    membership_counts$Expected_models !=
    membership_counts$Observed_models
  )
) {
  stop(
    "Membership counts do not match Table S7. ",
    "Check Membership_counts_per_SNP.csv."
  )
}

stopifnot(
  dplyr::n_distinct(trait_membership$SNP) == 46L,
  dplyr::n_distinct(model_membership$SNP) == 46L
)

cat("All 46 SNPs have matching trait and model counts.\n")

####################################################################
# 6. Create UpSet inputs
####################################################################

trait_sets <- split(
  trait_membership$SNP,
  trait_membership$Trait
)

trait_sets <- lapply(trait_sets, unique)

# Order traits by their number of SNPs
trait_sets <- trait_sets[
  order(lengths(trait_sets), decreasing = TRUE)
]

model_sets <- split(
  model_membership$SNP,
  model_membership$Model
)

model_sets <- lapply(model_sets, unique)

model_order <- c(
  "BLINK",
  "FarmCPU",
  "MLM",
  "MLMM"
)

model_order <- model_order[
  model_order %in% names(model_sets)
]

model_sets <- model_sets[model_order]

if (length(trait_sets) < 2L || length(model_sets) < 2L) {
  stop("At least two sets are required for each UpSet plot.")
}

trait_input <- UpSetR::fromList(trait_sets)
model_input <- UpSetR::fromList(model_sets)

stopifnot(
  nrow(trait_input) == 46L,
  nrow(model_input) == 46L
)

####################################################################
# 7. Summarise exact combinations
####################################################################

summarise_combinations <- function(membership, member_column) {
  
  membership %>%
    dplyr::group_by(SNP) %>%
    dplyr::summarise(
      Combination = paste(
        sort(unique(.data[[member_column]])),
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    dplyr::group_by(Combination) %>%
    dplyr::summarise(
      N_SNPs = dplyr::n(),
      SNPs = paste(sort(SNP), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(N_SNPs),
      Combination
    )
}

trait_combinations <- summarise_combinations(
  trait_membership,
  "Trait"
)

model_combinations <- summarise_combinations(
  model_membership,
  "Model"
)

stopifnot(
  sum(trait_combinations$N_SNPs) == 46L,
  sum(model_combinations$N_SNPs) == 46L
)

cat(
  "Number of exact trait combinations:",
  nrow(trait_combinations),
  "\n"
)

cat(
  "Number of exact model combinations:",
  nrow(model_combinations),
  "\n"
)

####################################################################
# 8. Define plotting functions
####################################################################

draw_trait_plot <- function() {
  
  p <- UpSetR::upset(
    trait_input,
    sets = names(trait_sets),
    keep.order = TRUE,
    nsets = length(trait_sets),
    nintersects = NA,
    order.by = "freq",
    mainbar.y.label =
      "SNPs per exact trait combination (total n = 46)",
    sets.x.label = "SNPs per trait",
    main.bar.color = "#4E79A7",
    sets.bar.color = "grey45",
    matrix.color = "black",
    point.size = 3,
    line.size = 0.8,
    text.scale = 1.2,
    mb.ratio = c(0.55, 0.45)
  )
  
  print(p)
}

draw_model_plot <- function() {
  
  p <- UpSetR::upset(
    model_input,
    sets = names(model_sets),
    keep.order = TRUE,
    nsets = length(model_sets),
    nintersects = NA,
    order.by = "freq",
    mainbar.y.label =
      "SNPs per exact model combination (total n = 46)",
    sets.x.label = "SNPs per model",
    main.bar.color = "#B07AA1",
    sets.bar.color = "grey45",
    matrix.color = "black",
    point.size = 4,
    line.size = 1,
    text.scale = 1.3,
    mb.ratio = c(0.65, 0.35)
  )
  
  print(p)
}

####################################################################
# 9. Save vector PDFs
####################################################################

save_plot_pdf <- function(
    filename,
    width,
    height,
    draw_function
) {
  
  grDevices::pdf(
    filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
  draw_function()
}

trait_pdf <- file.path(
  output_folder,
  "Figure_S27A_Trait_combinations_46SNPs.pdf"
)

model_pdf <- file.path(
  output_folder,
  "Figure_S27B_Model_combinations_46SNPs.pdf"
)

save_plot_pdf(
  filename = trait_pdf,
  width = max(18, 5 + 0.4 * nrow(trait_combinations)),
  height = max(9, 4 + 0.30 * length(trait_sets)),
  draw_function = draw_trait_plot
)

save_plot_pdf(
  filename = model_pdf,
  width = 12,
  height = 8,
  draw_function = draw_model_plot
)

# Display plots in RStudio
if (interactive()) {
  draw_trait_plot()
  draw_model_plot()
}

####################################################################
# 10. Export supporting tables
####################################################################

readr::write_csv(
  snp_table,
  file.path(output_folder, "Selected_46_SNPs.csv")
)

readr::write_csv(
  trait_membership,
  file.path(output_folder, "Trait_memberships.csv")
)

readr::write_csv(
  model_membership,
  file.path(output_folder, "Model_memberships.csv")
)

readr::write_csv(
  trait_combinations,
  file.path(output_folder, "Trait_combinations.csv")
)

readr::write_csv(
  model_combinations,
  file.path(output_folder, "Model_combinations.csv")
)

writeLines(
  c(
    paste("Input file:", input_file),
    "Worksheet: 1",
    "Unique SNPs: 46",
    paste("Run date:", Sys.Date()),
    "",
    capture.output(sessionInfo())
  ),
  file.path(output_folder, "Input_and_sessionInfo.txt")
)

####################################################################
# 11. Print completion summary
####################################################################

cat(
  "\nFinished.\n",
  "\nSNPs in trait panel: 46",
  "\nSNPs in model panel: 46",
  "\n\nResults saved to:\n",
  output_folder,
  "\n"
)


#####################################################################################################################
# Supplemental Figure S10
#####################################################################################################################
# rrBLUP: BEAGLE-imputed C. annuum core
# Train using all available observed phenotypes for each trait
# Predict GEBVs for all genotyped accessions (core)
# Calculate Pas for control (core)
####################################################################

library(vcfR)
library(rrBLUP)

####################################################################
# 1. File paths
####################################################################

#folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core"
folder <- "/home/ahmccorm/scratch/world_veg/core_n397_PA"
pheno_folder <- folder

output_folder <- file.path(folder, "GS_all_phenotyped")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

####################################################################
# 2. Read BEAGLE-imputed genotypes
####################################################################

vcf <- read.vcfR(
  file.path(folder, "Pepper_397_annuum_LD_pruned_beagle.vcf.gz")
)

gt <- extract.gt(vcf, element = "GT")

# Handle phased and unphased diploid biallelic genotypes
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

# Individuals as rows; SNPs as columns
geno <- t(geno_numeric)

if (anyNA(geno)) {
  stop("Missing or unsupported genotype calls remain in the VCF.")
}

if (anyDuplicated(rownames(geno))) {
  stop("Duplicate sample IDs in the VCF.")
}

cat("Genotyped accessions:", nrow(geno), "\n")
cat("SNPs:", ncol(geno), "\n")

# Release large intermediate objects
rm(vcf, gt, geno_numeric)
gc()

####################################################################
# 3. Read phenotypes
####################################################################

phenotypes <- list(
  Control = read.csv(
    file.path(pheno_folder, "control_pheno_filtered_286.csv"),
    row.names = 1,
    check.names = FALSE
  ),
  Heat_Stress_1 = read.csv(
    file.path(pheno_folder, "heat_stress_1_filtered_286.csv"),
    row.names = 1,
    check.names = FALSE
  ),
  Heat_Stress_2 = read.csv(
    file.path(pheno_folder, "heat_stress_2_filtered_286.csv"),
    row.names = 1,
    check.names = FALSE
  )
)

####################################################################
# 4. Fit each trait using all available phenotype observations
####################################################################

run_rrBLUP_all_traits <- function(geno, pheno, environment,
                                  output_folder) {
  
  stopifnot(
    !anyDuplicated(rownames(pheno)),
    all(vapply(pheno, is.numeric, logical(1)))
  )
  
  # Match by accession ID and put phenotypes in genotype order
  matched_ids <- rownames(geno)[
    rownames(geno) %in% rownames(pheno)
  ]
  
  pheno <- pheno[matched_ids, , drop = FALSE]
  
  cat("\n", environment, "\n", sep = "")
  cat("Accessions with phenotype records:", length(matched_ids), "\n")
  
  if (length(matched_ids) < 3) {
    stop("Too few matching phenotype and genotype IDs.")
  }
  
  traits <- colnames(pheno)
  
  # One prediction column per trait, for every genotyped accession
  gebv <- matrix(
    NA_real_,
    nrow = nrow(geno),
    ncol = length(traits),
    dimnames = list(rownames(geno), traits)
  )
  
  summary_list <- list()
  training_list <- list()
  
  for (trait in traits) {
    
    # Use every observed phenotype for this particular trait
    observed <- is.finite(pheno[[trait]])
    train_ids <- matched_ids[observed]
    y <- pheno[[trait]][observed]
    
    training_list[[trait]] <- data.frame(
      Environment = environment,
      Trait = trait,
      Sample_ID = train_ids
    )
    
    status <- "Fitted"
    
    if (length(y) < 3) {
      
      status <- "Skipped: fewer than 3 observations"
      
    } else if (var(y) == 0) {
      
      status <- "Skipped: zero phenotype variance"
      
    } else {
      
      cat("Fitting", trait, "- n =", length(y), "\n")
      
      Z_train <- geno[train_ids, , drop = FALSE]
      
      # Remove SNPs that are invariant in this trait's training set
      keep <- colSums(
        sweep(Z_train, 2, Z_train[1, ], FUN = "!=")
      ) > 0
      
      if (!any(keep)) {
        
        status <- "Skipped: no variable SNPs"
        
      } else {
        
        # Centre markers using training-set means
        marker_means <- colMeans(Z_train[, keep, drop = FALSE])
        
        Z_train <- sweep(
          Z_train[, keep, drop = FALSE],
          2, marker_means, FUN = "-"
        )
        
        model <- mixed.solve(y = y, Z = Z_train)
        
        # Marker contributions for all accessions, centred identically
        gebv[, trait] <- as.numeric(
          geno[, keep, drop = FALSE] %*% model$u
        ) - sum(marker_means * as.numeric(model$u))
      }
    }
    
    summary_list[[trait]] <- data.frame(
      Environment = environment,
      Trait = trait,
      N_training = length(train_ids),
      N_without_observed_phenotype = nrow(geno) - length(train_ids),
      Status = status
    )
  }
  
  # Save all predictions together; no Excel merging needed
  write.csv(
    data.frame(
      Sample_ID = rownames(geno),
      gebv,
      row.names = NULL,
      check.names = FALSE
    ),
    file.path(
      output_folder,
      paste0("GEBVs_all_", nrow(geno), "_", environment, ".csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    do.call(rbind, summary_list),
    file.path(
      output_folder,
      paste0("Training_summary_", environment, ".csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    do.call(rbind, training_list),
    file.path(
      output_folder,
      paste0("Training_IDs_by_trait_", environment, ".csv")
    ),
    row.names = FALSE
  )
  
  return(invisible(gebv))
}

####################################################################
# 5. Run all three environments
####################################################################

for (environment in names(phenotypes)) {
  
  run_rrBLUP_all_traits(
    geno = geno,
    pheno = phenotypes[[environment]],
    environment = environment,
    output_folder = output_folder
  )
}

#########################################################################################################
# Step 3: Model Prediction Accuracies — rrBLUP
#########################################################################################################
library(rrBLUP)
#folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core"
folder <- "/home/ahmccorm/scratch/world_veg/core_n397_PA"
####################
# STEP ONE
####################

# Use the BEAGLE genotype matrix already created by your GS code
geno_data <- geno

# Load phenotype data
control_pheno_check <- read.csv(
  file.path(folder, "control_pheno_filtered_286.csv"),
  row.names = 1
)

heat_stress_1_check <- read.csv(
  file.path(folder, "heat_stress_1_filtered_286.csv"),
  row.names = 1
)

heat_stress_2_check <- read.csv(
  file.path(folder, "heat_stress_2_filtered_286.csv"),
  row.names = 1
)

# Identify common genotypes
common_genotypes <- intersect(
  rownames(geno_data),
  rownames(control_pheno_check)
)

geno_filtered <- geno_data[common_genotypes, , drop = FALSE]

geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
# STEP TWO
####################

# Use all matching phenotyped accessions instead of the old 150
train_genotypes <- common_genotypes

geno_train <- geno_matrix[train_genotypes, , drop = FALSE]

control_train <- control_pheno_check[train_genotypes, , drop = FALSE]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, , drop = FALSE]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, , drop = FALSE]

####################
# STEP THREE
####################

# Use your original function unchanged
source(file.path(folder, "xval_kfold_functions.R"))

# Run the same 10-fold cross-validation × 50 repetitions
xval_k10_rrblup <- k.xval(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_k10_rrblup,
  file.path(folder, "xval_rrblup_kfold_10.RData")
)


####################################################################
# Method Two: Gaussian Kernel
####################################################################

# Separate phenotype object so the rrBLUP input stays unchanged
control_train_kernel <- data.frame(
  Genotype = rownames(control_train),
  control_train,
  row.names = rownames(control_train),
  check.names = FALSE
)

stopifnot(
  identical(rownames(geno_train), rownames(control_train_kernel))
)

# A.mat expects -1/0/1 genotype coding
K <- A.mat(geno_train - 1)
k_dist <- dist(K)

# Run your original Gaussian cross-validation function
xval_k10_GAUSS <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = control_train_kernel,
  y.trainset = control_train_kernel,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_k10_GAUSS,
  file.path(folder, "xval_GAUSS_kfold_10.RData")
)

####################################################################
# Method Three: Exponential Kernel — original implementation
####################################################################

# Keep your original exponential kernel construction
K.Exp <- Kernel_computation(
  X = geno_train,
  name = "exponential",
  degree = NULL,
  nL = NULL
)

exp_dist <- dist(K.Exp)

exp_dist_mat <- as.matrix(exp_dist)

rownames(exp_dist_mat) <- colnames(exp_dist_mat) <-
  rownames(geno_train)

stopifnot(
  identical(
    rownames(exp_dist_mat),
    control_train_kernel$Genotype
  )
)

# Run your original exponential cross-validation function
xval_k10_EXP <- k.xval.EXP(
  g.in = geno_train,
  y.in = control_train_kernel,
  y.trainset = control_train_kernel,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_k10_EXP,
  file.path(folder, "xval_EXP_kfold_10.RData")
)


#####################################################################################################################
# Supplemental Figure S11/12
#####################################################################################################################
# PAs for Heat Stress 1 and Heat Stress 2 (core)
####################################################################

heat_phenotypes <- list(
  Heat_Stress_1 = heat_stress_1_check,
  Heat_Stress_2 = heat_stress_2_check
)

for (environment in names(heat_phenotypes)) {
  
  cat("\nStarting PAs for", environment, "\n")
  
  # Match this environment's phenotype records to the genotypes
  pheno_env <- heat_phenotypes[[environment]]
  
  matched_ids <- intersect(
    rownames(geno),
    rownames(pheno_env)
  )
  
  geno_env <- geno[matched_ids, , drop = FALSE]
  pheno_env <- pheno_env[matched_ids, , drop = FALSE]
  
  ##################################################################
  # rrBLUP
  ##################################################################
  
  result_rrblup <- k.xval(
    g.in = geno_env,
    y.in = pheno_env,
    y.trainset = pheno_env,
    k.fold = 10,
    reps = 50
  )
  
  saveRDS(
    result_rrblup,
    file.path(
      folder,
      paste0("xval_rrblup_", environment, "_kfold_10.RData")
    )
  )
  
  ##################################################################
  # Gaussian
  ##################################################################
  
  pheno_kernel <- data.frame(
    Genotype = rownames(pheno_env),
    pheno_env,
    row.names = rownames(pheno_env),
    check.names = FALSE
  )
  
  K_env <- A.mat(geno_env - 1)
  gaussian_dist <- dist(K_env)
  
  result_gauss <- k.xval.GAUSS(
    g.in = geno_env,
    y.in = pheno_kernel,
    y.trainset = pheno_kernel,
    k_dist = gaussian_dist,
    k.fold = 10,
    reps = 50
  )
  
  saveRDS(
    result_gauss,
    file.path(
      folder,
      paste0("xval_GAUSS_", environment, "_kfold_10.RData")
    )
  )
  
  ##################################################################
  # Exponential — your original implementation
  ##################################################################
  
  K_exp <- Kernel_computation(
    X = geno_env,
    name = "exponential",
    degree = NULL,
    nL = NULL
  )
  
  exponential_dist <- as.matrix(dist(K_exp))
  
  rownames(exponential_dist) <- colnames(exponential_dist) <-
    rownames(geno_env)
  
  result_exp <- k.xval.EXP(
    g.in = geno_env,
    y.in = pheno_kernel,
    y.trainset = pheno_kernel,
    k_dist = exponential_dist,
    k.fold = 10,
    reps = 50
  )
  
  saveRDS(
    result_exp,
    file.path(
      folder,
      paste0("xval_EXP_", environment, "_kfold_10.RData")
    )
  )
}


#####################################################################################################################
# Supplemental Figure S10/S11/S12 PLOTS
#####################################################################################################################

####################################################################
# Core collection: PA plots
# rrBLUP and Gaussian — Control and Heat Stress 1
####################################################################

library(ggplot2)
library(dplyr)

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core/PA"

environments <- c(
  Control = "",
  Heat_Stress_1 = "_Heat_Stress_1",
  Heat_Stress_2 = "_Heat_Stress_2"
)

environment_labels <- c(
  Control = "Control",
  Heat_Stress_1 = "Heat Stress 1",
  Heat_Stress_2 = "Heat Stress 2"
)

model_labels <- c(
  rrblup = "rrBLUP",
  GAUSS = "Gaussian Kernel",
  EXP = "Exponential Kernel"
)
plots_core <- list()
pa_tables_core <- list()

for (environment in names(environments)) {
  
  suffix <- environments[[environment]]
  
  ##################################################################
  # Load models
  ##################################################################
  
  model_results <- lapply(names(model_labels), function(method) {
    
    filename <- paste0(
      "xval_", method, suffix, "_kfold_10.RData"
    )
    
    results <- readRDS(
      file.path(folder, filename)
    )$xval.result
    
    results$model <- model_labels[[method]]
    
    return(results)
  })
  
  all_gs_traits <- bind_rows(model_results) %>%
    mutate(
      r.mean = as.numeric(as.character(r.mean)),
      r.sd = as.numeric(as.character(r.sd)),
      model = factor(
        model,
        levels = unname(model_labels)
      )
    ) %>%
    filter(
      !is.na(trait),
      is.finite(r.mean),
      is.finite(r.sd)
    )
  
  pa_tables_core[[environment]] <- all_gs_traits
  
  ##################################################################
  # Plot
  ##################################################################
  
  p <- ggplot(
    all_gs_traits,
    aes(y = r.mean, x = model, color = model)
  ) +
    theme_bw() +
    geom_errorbar(
      aes(
        ymin = r.mean - r.sd,
        ymax = r.mean + r.sd
      ),
      width = 0.3,
      position = position_dodge(0.2)
    ) +
    geom_point(size = 3) +
    facet_wrap(
      vars(trait),
      scales = "free_x",
      nrow = 5
    ) +
    geom_hline(
      yintercept = 0.5,
      color = "red",
      linewidth = 1.5,
      linetype = "longdash"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 15)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(
        "Core Collection: Prediction Accuracy (r) for ",
        n_distinct(all_gs_traits$trait),
        " Traits — ",
        environment_labels[[environment]]
      ),
      x = "Model",
      y = "Prediction Accuracy (r)"
    )
  
  print(p)
  plots_core[[environment]] <- p
  
  ggsave(
    filename = file.path(
      folder,
      paste0("Core_PA_", environment, "_rrBLUP_Gaussian.pdf")
    ),
    plot = p,
    width = 24,
    height = 13,
    units = "in"
  )
}

#####################################################################################################################
# Supplemental Figure S13
#####################################################################################################################
# top X above 0.5
# Traits with rrBLUP PA > 0.5 in all three environments

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

####################################################################
# 1. Paths
####################################################################

base_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1"

pa_folder <- file.path(base_folder, "GS_core", "PA")
gs_folder <- file.path(base_folder, "GS_core", "GS_all_phenotyped")

output_folder <- file.path(base_folder, "S11")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

environments <- c("Control", "Heat_Stress_1", "Heat_Stress_2")

####################################################################
# 2. Select traits using the new core rrBLUP PAs
####################################################################

pa_results <- dplyr::bind_rows(
  lapply(environments, function(environment) {
    
    suffix <- if (environment == "Control") {
      ""
    } else {
      paste0("_", environment)
    }
    
    filename <- paste0(
      "xval_rrblup", suffix, "_kfold_10.RData"
    )
    
    readRDS(file.path(pa_folder, filename))$xval.result %>%
      dplyr::transmute(
        Trait = as.character(trait),
        PA = as.numeric(as.character(r.mean)),
        Environment = environment
      ) %>%
      dplyr::filter(!is.na(Trait), nzchar(Trait))
  })
)

stopifnot(
  !anyDuplicated(pa_results[c("Trait", "Environment")])
)

top_traits <- pa_results %>%
  dplyr::group_by(Trait) %>%
  dplyr::summarise(
    Pass = dplyr::n_distinct(Environment) == 3 &&
      all(is.finite(PA) & PA > 0.5),
    .groups = "drop"
  ) %>%
  dplyr::filter(Pass) %>%
  dplyr::arrange(Trait) %>%
  dplyr::pull(Trait)

cat("Traits passing all three environments:", length(top_traits), "\n")
print(top_traits)

if (length(top_traits) == 0) {
  stop("No traits have rrBLUP PA > 0.5 in all three environments.")
}

write.csv(
  pa_results %>%
    dplyr::filter(Trait %in% top_traits),
  file.path(output_folder, "S11_Selected_traits_PAs.csv"),
  row.names = FALSE
)

####################################################################
# 3. Load new GEBVs
####################################################################

gebv_all <- dplyr::bind_rows(
  lapply(environments, function(environment) {
    
    filename <- paste0(
      "GEBVs_all_397_", environment, ".csv"
    )
    
    df <- read.csv(
      file.path(gs_folder, filename),
      check.names = FALSE
    )
    
    stopifnot(
      all(c("Sample_ID", top_traits) %in% names(df))
    )
    
    df$Sample_ID <- trimws(as.character(df$Sample_ID))
    
    stopifnot(!anyDuplicated(df$Sample_ID))
    
    df %>%
      dplyr::select(
        Sample_ID,
        tidyselect::all_of(top_traits)
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(top_traits),
        names_to = "Trait",
        values_to = "GEBV"
      ) %>%
      dplyr::mutate(Condition = environment)
  })
)

# Retain accession/trait combinations with finite GEBVs
# in all three environments before calculating ranks
gebv_all <- gebv_all %>%
  dplyr::group_by(Sample_ID, Trait) %>%
  dplyr::filter(
    dplyr::n_distinct(Condition) == 3,
    all(is.finite(GEBV))
  ) %>%
  dplyr::ungroup()

####################################################################
# 4. Rank accessions within each trait/environment
# Rank 1 = highest GEBV
####################################################################

rank_long <- gebv_all %>%
  dplyr::group_by(Trait, Condition) %>%
  dplyr::mutate(
    Rank = rank(-GEBV, ties.method = "average")
  ) %>%
  dplyr::ungroup()

####################################################################
# 5. Classify patterns in order: Heat 1 → Control → Heat 2
####################################################################

pattern_levels <- c(
  "Decreasing",
  "Increasing",
  "U-Shaped",
  "Inverted-U",
  "Irregular"
)

get_pattern <- function(r1, r2, r3) {
  
  if (any(!is.finite(c(r1, r2, r3)))) {
    return(NA_character_)
  }
  
  if (r1 > r2 && r2 > r3) return("Decreasing")
  if (r1 < r2 && r2 < r3) return("Increasing")
  if (r1 < r2 && r3 < r2) return("U-Shaped")
  if (r1 > r2 && r3 > r2) return("Inverted-U")
  
  "Irregular"
}

rank_wide <- rank_long %>%
  dplyr::select(Sample_ID, Trait, Condition, Rank) %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = Rank
  ) %>%
  dplyr::mutate(
    Pattern = mapply(
      get_pattern,
      Heat_Stress_1,
      Control,
      Heat_Stress_2
    ),
    Pattern = factor(Pattern, levels = pattern_levels)
  )

rank_plot_data <- rank_long %>%
  dplyr::left_join(
    rank_wide %>% dplyr::select(Sample_ID, Trait, Pattern),
    by = c("Sample_ID", "Trait")
  ) %>%
  dplyr::mutate(
    Condition = factor(
      Condition,
      levels = c("Heat_Stress_1", "Control", "Heat_Stress_2"),
      labels = c("Heat Stress 1", "Control", "Heat Stress 2")
    )
  )

####################################################################
# 6. Generate panels
####################################################################

pattern_colors <- c(
  "Decreasing" = "#1B9E77",
  "Increasing" = "#D95F02",
  "U-Shaped" = "#7570B3",
  "Inverted-U" = "#E7298A",
  "Irregular" = "#66A61E"
)

plot_list <- lapply(top_traits, function(trait_name) {
  
  plot_data <- rank_plot_data %>%
    dplyr::filter(Trait == trait_name)
  
  ggplot(
    plot_data,
    aes(
      x = Condition,
      y = Rank,
      group = Sample_ID,
      color = Pattern
    )
  ) +
    geom_line(alpha = 0.5) +
    geom_point(size = 0.8) +
    scale_y_reverse() +
    scale_color_manual(
      values = pattern_colors,
      limits = pattern_levels,
      drop = FALSE
    ) +
    labs(
      title = trait_name,
      x = NULL,
      y = "Rank (1 = highest GEBV)",
      color = "Pattern"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

names(plot_list) <- top_traits

####################################################################
# 7. Combine with one shared legend
####################################################################

final_plot <- patchwork::wrap_plots(
  plot_list,
  ncol = 5,
  guides = "collect"
) &
  theme(legend.position = "right")

print(final_plot)

####################################################################
# 8. Save PDF and rank-pattern tables
####################################################################

ggsave(
  filename = file.path(
    output_folder,
    "S11_Core397_rank_change_patterns.pdf"
  ),
  plot = final_plot,
  width = 22,
  height = max(8, ceiling(length(top_traits) / 5) * 3.5),
  units = "in",
  limitsize = FALSE
)

write.csv(
  rank_wide,
  file.path(output_folder, "S11_Core397_ranks_and_patterns.csv"),
  row.names = FALSE
)

pattern_counts <- rank_wide %>%
  dplyr::count(Trait, Pattern, .drop = FALSE)

write.csv(
  pattern_counts,
  file.path(output_folder, "S11_Core397_pattern_counts.csv"),
  row.names = FALSE
)


#####################################################################################################################
# Supplemental Figure S14
#####################################################################################################################

####################################################################
# Figure 14A
####################################################################

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

####################################################################
# Figure 14B: Core agronomic GEBV distributions by environment
# Dashed lines indicate medians
####################################################################

library(tidyverse)

####################################################################
# 1. Paths
####################################################################

base_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1"

core_pa_folder <- file.path(base_folder, "GS_core", "PA")
global_pa_folder <- file.path(base_folder, "GS_global", "PA")

core_gs_folder <- file.path(
  base_folder, "GS_core", "GS_all_phenotyped"
)

output_folder <- file.path(base_folder, "Figure2B")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

environments <- c(
  "Control",
  "Heat_Stress_1",
  "Heat_Stress_2"
)

condition_labels <- c(
  Control = "Control",
  Heat_Stress_1 = "Heat Stress 1",
  Heat_Stress_2 = "Heat Stress 2"
)

####################################################################
# 2. Load core and global rrBLUP PAs
####################################################################

load_pa <- function(folder, collection) {
  
  purrr::map_dfr(environments, function(environment) {
    
    suffix <- if (environment == "Control") {
      ""
    } else {
      paste0("_", environment)
    }
    
    filename <- paste0(
      "xval_rrblup", suffix, "_kfold_10.RData"
    )
    
    readRDS(file.path(folder, filename))$xval.result %>%
      dplyr::transmute(
        Trait = as.character(trait),
        PA = as.numeric(as.character(r.mean)),
        Collection = collection,
        Environment = environment
      ) %>%
      dplyr::filter(!is.na(Trait), nzchar(Trait))
  })
}

all_pa <- dplyr::bind_rows(
  load_pa(core_pa_folder, "Core"),
  load_pa(global_pa_folder, "Global")
)

stopifnot(
  !anyDuplicated(all_pa[c("Trait", "Collection", "Environment")])
)

####################################################################
# 3. Select traits passing all three conditions in both collections
####################################################################

qualifying_traits <- all_pa %>%
  dplyr::group_by(Collection, Trait) %>%
  dplyr::summarise(
    Pass = dplyr::n_distinct(Environment) == 3 &&
      all(is.finite(PA) & PA > 0.5),
    .groups = "drop"
  ) %>%
  dplyr::filter(Pass)

core_traits <- qualifying_traits %>%
  dplyr::filter(Collection == "Core") %>%
  dplyr::pull(Trait)

global_traits <- qualifying_traits %>%
  dplyr::filter(Collection == "Global") %>%
  dplyr::pull(Trait)

common_traits <- sort(intersect(core_traits, global_traits))

cat("Common qualifying traits:", length(common_traits), "\n")
print(common_traits)

if (length(common_traits) == 0) {
  stop("No traits pass the PA threshold in all six combinations.")
}

####################################################################
# 4. Load core GEBVs for selected traits
####################################################################

gebv_core_only <- purrr::map_dfr(
  environments,
  function(environment) {
    
    filename <- paste0(
      "GEBVs_all_397_", environment, ".csv"
    )
    
    df <- readr::read_csv(
      file.path(core_gs_folder, filename),
      show_col_types = FALSE
    )
    
    stopifnot(
      all(c("Sample_ID", common_traits) %in% names(df)),
      !anyDuplicated(df$Sample_ID)
    )
    
    df %>%
      dplyr::select(
        Sample_ID,
        tidyselect::all_of(common_traits)
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(common_traits),
        names_to = "Trait",
        values_to = "GEBV"
      ) %>%
      dplyr::mutate(
        Condition = condition_labels[[environment]]
      )
  }
) %>%
  dplyr::filter(is.finite(GEBV)) %>%
  dplyr::mutate(
    Trait = factor(Trait, levels = common_traits),
    Condition = factor(
      Condition,
      levels = c("Control", "Heat Stress 1", "Heat Stress 2")
    )
  )

####################################################################
# 5. Calculate medians
####################################################################
####################################################################
# 5. Calculate means and quartiles by trait and environment
####################################################################

distribution_summary <- gebv_core_only %>%
  dplyr::group_by(Trait, Condition) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Mean = mean(GEBV),
    Q25 = unname(quantile(GEBV, 0.25)),
    Median = median(GEBV),
    Q75 = unname(quantile(GEBV, 0.75)),
    .groups = "drop"
  )

reference_lines <- distribution_summary %>%
  tidyr::pivot_longer(
    cols = c(Mean, Q25, Median, Q75),
    names_to = "Statistic",
    values_to = "Position"
  ) %>%
  dplyr::mutate(
    Line_type = dplyr::case_when(
      Statistic == "Mean" ~ "Mean",
      Statistic == "Median" ~ "Median",
      TRUE ~ "25th / 75th percentile"
    ),
    Line_type = factor(
      Line_type,
      levels = c(
        "Mean",
        "Median",
        "25th / 75th percentile"
      )
    )
  )

####################################################################
# 6. Readable facet labels
####################################################################

pretty_labels <- c(
  "DASanthesis" = "Days after sowing to anthesis",
  "DASmaturity" = "Days after sowing to maturity",
  "DATanthesis" = "Days after transplanting to anthesis",
  "DATmaturity" = "Days after transplanting to maturity",
  "biomassfinalplant" = "Final plant biomass",
  "biomassfinalplot" = "Final plot biomass",
  "fruitlength" = "Fruit length",
  "fruitno" = "Fruit number",
  "fruitshapeindex" = "Fruit shape index",
  "fruitweight" = "Fruit weight",
  "fruitwidth" = "Fruit width",
  "height" = "Height",
  "heightfinal" = "Final height",
  "leafangledaily" = "Daily leaf angle",
  "leafareadaily" = "Daily leaf area",
  "leafareadailyplant" = "Daily plant leaf area",
  "leafareaprojecteddaily" = "Daily projected leaf area",
  "leafareaprojectedplant" = "Projected plant leaf area",
  "yield" = "Yield"
)

# Retain original name if a label has not been defined
facet_labels <- setNames(common_traits, common_traits)

labelled_traits <- intersect(
  common_traits,
  names(pretty_labels)
)

facet_labels[labelled_traits] <- pretty_labels[labelled_traits]

####################################################################
####################################################################
# 7. Plot with means and quartiles
####################################################################

condition_colors <- c(
  "Control" = "#1f77b4",
  "Heat Stress 1" = "#ff7f0e",
  "Heat Stress 2" = "red"
)

p <- ggplot(
  gebv_core_only,
  aes(x = GEBV, fill = Condition)
) +
  geom_density(alpha = 0.5) +
  
  geom_vline(
    data = reference_lines,
    aes(
      xintercept = Position,
      color = Condition,
      linetype = Line_type
    ),
    linewidth = 0.5,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  
  facet_wrap(
    ~ Trait,
    scales = "free",
    ncol = 4,
    labeller = labeller(Trait = facet_labels)
  ) +
  
  scale_fill_manual(
    values = condition_colors,
    name = "Condition"
  ) +
  
  scale_color_manual(
    values = condition_colors,
    name = "Condition"
  ) +
  
  scale_linetype_manual(
    values = c(
      "Mean" = "solid",
      "Median" = "dashed",
      "25th / 75th percentile" = "dotted"
    ),
    name = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = "none",
    linetype = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(color = "black")
    )
  ) +
  
  labs(
    title = "",
    x = "GEBV",
    y = "Density"
  )

print(p)

####################################################################
# 8. Save using your original dimensions and output folder
####################################################################

ggsave(
  filename = file.path(
    output_folder,
    "F2B_GEBV_core397_density.pdf"
  ),
  plot = p,
  width = 15,
  height = max(8, ceiling(length(common_traits) / 4) * 2.5),
  units = "in"
)

readr::write_csv(
  distribution_summary,
  file.path(
    output_folder,
    "S12B_Core397_GEBV_mean_quantiles.csv"
  )
)

readr::write_csv(
  all_pa %>%
    dplyr::filter(Trait %in% common_traits) %>%
    dplyr::arrange(Trait, Collection, Environment),
  file.path(output_folder, "F2B_Qualifying_trait_PAs.csv")
)

#####################################################################################################################
# Supplemental Figure S15
#####################################################################################################################
# Yield GEBV rank changes across conditions — core n = 397

library(ggplot2)
library(dplyr)
library(tidyr)

####################################################################
# 1. Paths and import
####################################################################

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core/GS_all_phenotyped"

output_folder <- file.path(folder, "S13")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

load_yield <- function(filename, value_name) {
  
  df <- read.csv(
    file.path(folder, filename),
    check.names = FALSE
  )
  
  stopifnot(
    all(c("Sample_ID", "yield") %in% names(df)),
    !anyDuplicated(df$Sample_ID)
  )
  
  df <- df %>%
    select(Sample_ID, yield) %>%
    mutate(Sample_ID = trimws(as.character(Sample_ID)))
  
  names(df)[2] <- value_name
  df
}

GEBVs_control <- load_yield(
  "GEBVs_all_397_Control.csv",
  "GEBV_Control"
)

GEBVs_heat1 <- load_yield(
  "GEBVs_all_397_Heat_Stress_1.csv",
  "GEBV_Heat1"
)

GEBVs_heat2 <- load_yield(
  "GEBVs_all_397_Heat_Stress_2.csv",
  "GEBV_Heat2"
)

# Match accessions by ID across all three files
GEBVs_all <- GEBVs_control %>%
  inner_join(GEBVs_heat1, by = "Sample_ID") %>%
  inner_join(GEBVs_heat2, by = "Sample_ID") %>%
  filter(
    if_all(
      c(GEBV_Control, GEBV_Heat1, GEBV_Heat2),
      ~ is.finite(.x)
    )
  )

cat("Accessions with yield GEBVs in all three conditions:",
    nrow(GEBVs_all), "\n")

####################################################################
# 2. Rank accessions
# Higher yield GEBV = better rank; rank 1 is best
####################################################################

GEBVs_all <- GEBVs_all %>%
  mutate(
    Rank_Heat1 = rank(-GEBV_Heat1, ties.method = "average"),
    Rank_Control = rank(-GEBV_Control, ties.method = "average"),
    Rank_Heat2 = rank(-GEBV_Heat2, ties.method = "average"),
    
    Rank_change_intensity =
      abs(Rank_Control - Rank_Heat1) +
      abs(Rank_Control - Rank_Heat2)
  )

# Same threshold as before: upper 20% across all accessions
threshold <- quantile(
  GEBVs_all$Rank_change_intensity,
  probs = 0.80,
  names = FALSE
)

####################################################################
# 3. Select the two rank-change patterns
####################################################################

# Better relative ranking in BOTH heat-stress conditions
extreme_V_genotypes <- GEBVs_all %>%
  filter(
    Rank_Heat1 < Rank_Control,
    Rank_Heat2 < Rank_Control,
    Rank_change_intensity >= threshold
  ) %>%
  arrange(desc(Rank_change_intensity))

# Worse relative ranking in BOTH heat-stress conditions
extreme_invertedV_genotypes <- GEBVs_all %>%
  filter(
    Rank_Heat1 > Rank_Control,
    Rank_Heat2 > Rank_Control,
    Rank_change_intensity >= threshold
  ) %>%
  arrange(desc(Rank_change_intensity))

cat("Heat-favoured rank pattern:", nrow(extreme_V_genotypes), "\n")
cat("Control-favoured rank pattern:",
    nrow(extreme_invertedV_genotypes), "\n")

####################################################################
# 4. Reshape ranks for plotting
####################################################################

GEBVs_long <- GEBVs_all %>%
  select(Sample_ID, Rank_Heat1, Rank_Control, Rank_Heat2) %>%
  pivot_longer(
    cols = c(Rank_Heat1, Rank_Control, Rank_Heat2),
    names_to = "Condition",
    values_to = "Rank"
  ) %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c("Rank_Heat1", "Rank_Control", "Rank_Heat2"),
      labels = c("Heat Stress 1", "Control", "Heat Stress 2")
    )
  )

v_genotype_plot_data <- GEBVs_long %>%
  filter(Sample_ID %in% extreme_V_genotypes$Sample_ID)

invertedV_plot_data <- GEBVs_long %>%
  filter(Sample_ID %in% extreme_invertedV_genotypes$Sample_ID)

####################################################################
# 5. Plot helper: same appearance as before
####################################################################

plot_rank_pattern <- function(data, title, line_color) {
  
  ggplot(
    data,
    aes(x = Condition, y = Rank, group = Sample_ID)
  ) +
    geom_line(
      color = line_color,
      alpha = 0.8,
      linewidth = 1.1
    ) +
    geom_point(color = line_color, size = 2) +
    scale_y_reverse() +
    labs(
      title = title,
      x = "Condition",
      y = "Yield GEBV rank (1 = best)"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      )
    )
}

####################################################################
# 6. Heat-favoured plot
####################################################################

p_heat <- plot_rank_pattern(
  v_genotype_plot_data,
  paste0(
    "Heat-favoured yield rankings (n = ",
    nrow(extreme_V_genotypes), ")"
  ),
  "firebrick3"
)

print(p_heat)

ggsave(
  filename = file.path(output_folder, "S13_Heat_favoured_core397.pdf"),
  plot = p_heat,
  width = 10,
  height = 7,
  units = "in"
)

####################################################################
# 7. Control-favoured plot
####################################################################

p_control <- plot_rank_pattern(
  invertedV_plot_data,
  paste0(
    "Control-favoured yield rankings (n = ",
    nrow(extreme_invertedV_genotypes), ")"
  ),
  "royalblue4"
)

print(p_control)

ggsave(
  filename = file.path(output_folder, "S13_Control_favoured_core397.pdf"),
  plot = p_control,
  width = 10,
  height = 7,
  units = "in"
)

####################################################################
# 8. Save selected accessions and all rankings
####################################################################

write.csv(
  extreme_V_genotypes,
  file.path(output_folder, "Heat_favoured_accessions_core397.csv"),
  row.names = FALSE
)

write.csv(
  extreme_invertedV_genotypes,
  file.path(output_folder, "Control_favoured_accessions_core397.csv"),
  row.names = FALSE
)

write.csv(
  GEBVs_all,
  file.path(output_folder, "Yield_GEBVs_and_ranks_core397.csv"),
  row.names = FALSE
)


#PA values for yield
library(dplyr)

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core/PA"

get_yield_pa <- function(filename, condition) {
  readRDS(file.path(folder, filename))$xval.result %>%
    filter(as.character(trait) == "yield") %>%
    transmute(
      Condition = condition,
      PA = as.numeric(as.character(r.mean)),
      SD = as.numeric(as.character(r.sd))
    )
}

yield_pa <- bind_rows(
  get_yield_pa(
    "xval_rrblup_kfold_10.RData",
    "Control"
  ),
  get_yield_pa(
    "xval_rrblup_Heat_Stress_1_kfold_10.RData",
    "Heat Stress 1"
  )
)

print(yield_pa)



#####################################################################################################################
# Supplemental Figure S16
#####################################################################################################################
# rrBLUP: BEAGLE-imputed C. annuum core
# Train using all available observed phenotypes for each trait
# Predict GEBVs for all genotyped accessions (global)
# Calculate Pas for control (global)

library(vcfR)
library(rrBLUP)

####################################################################
# 1. File paths
####################################################################

#folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core"
folder <- "/home/ahmccorm/scratch/world_veg/global_n7843_PA"
pheno_folder <- folder

output_folder <- file.path(folder, "GS_all_phenotyped")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

####################################################################
# 2. Read BEAGLE-imputed genotypes
####################################################################

vcf <- read.vcfR(
  file.path(folder, "pepper_10000_annuum_beagle.vcf.gz")
)

gt <- extract.gt(vcf, element = "GT")

# Handle phased and unphased diploid biallelic genotypes
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

# Individuals as rows; SNPs as columns
geno <- t(geno_numeric)

if (anyNA(geno)) {
  stop("Missing or unsupported genotype calls remain in the VCF.")
}

if (anyDuplicated(rownames(geno))) {
  stop("Duplicate sample IDs in the VCF.")
}

cat("Genotyped accessions:", nrow(geno), "\n")
cat("SNPs:", ncol(geno), "\n")

# Release large intermediate objects
rm(vcf, gt, geno_numeric)
gc()

####################################################################
# 3. Read phenotypes
####################################################################

phenotypes <- list(
  Control = read.csv(
    file.path(pheno_folder, "control_pheno_filtered_286.csv"),
    row.names = 1,
    check.names = FALSE
  ),
  Heat_Stress_1 = read.csv(
    file.path(pheno_folder, "heat_stress_1_filtered_286.csv"),
    row.names = 1,
    check.names = FALSE
  ),
  Heat_Stress_2 = read.csv(
    file.path(pheno_folder, "heat_stress_2_filtered_286.csv"),
    row.names = 1,
    check.names = FALSE
  )
)

####################################################################
# 4. Fit each trait using all available phenotype observations
####################################################################

run_rrBLUP_all_traits <- function(geno, pheno, environment,
                                  output_folder) {
  
  stopifnot(
    !anyDuplicated(rownames(pheno)),
    all(vapply(pheno, is.numeric, logical(1)))
  )
  
  # Match by accession ID and put phenotypes in genotype order
  matched_ids <- rownames(geno)[
    rownames(geno) %in% rownames(pheno)
  ]
  
  pheno <- pheno[matched_ids, , drop = FALSE]
  
  cat("\n", environment, "\n", sep = "")
  cat("Accessions with phenotype records:", length(matched_ids), "\n")
  
  if (length(matched_ids) < 3) {
    stop("Too few matching phenotype and genotype IDs.")
  }
  
  traits <- colnames(pheno)
  
  # One prediction column per trait, for every genotyped accession
  gebv <- matrix(
    NA_real_,
    nrow = nrow(geno),
    ncol = length(traits),
    dimnames = list(rownames(geno), traits)
  )
  
  summary_list <- list()
  training_list <- list()
  
  for (trait in traits) {
    
    # Use every observed phenotype for this particular trait
    observed <- is.finite(pheno[[trait]])
    train_ids <- matched_ids[observed]
    y <- pheno[[trait]][observed]
    
    training_list[[trait]] <- data.frame(
      Environment = environment,
      Trait = trait,
      Sample_ID = train_ids
    )
    
    status <- "Fitted"
    
    if (length(y) < 3) {
      
      status <- "Skipped: fewer than 3 observations"
      
    } else if (var(y) == 0) {
      
      status <- "Skipped: zero phenotype variance"
      
    } else {
      
      cat("Fitting", trait, "- n =", length(y), "\n")
      
      Z_train <- geno[train_ids, , drop = FALSE]
      
      # Remove SNPs that are invariant in this trait's training set
      keep <- colSums(
        sweep(Z_train, 2, Z_train[1, ], FUN = "!=")
      ) > 0
      
      if (!any(keep)) {
        
        status <- "Skipped: no variable SNPs"
        
      } else {
        
        # Centre markers using training-set means
        marker_means <- colMeans(Z_train[, keep, drop = FALSE])
        
        Z_train <- sweep(
          Z_train[, keep, drop = FALSE],
          2, marker_means, FUN = "-"
        )
        
        model <- mixed.solve(y = y, Z = Z_train)
        
        # Marker contributions for all accessions, centred identically
        gebv[, trait] <- as.numeric(
          geno[, keep, drop = FALSE] %*% model$u
        ) - sum(marker_means * as.numeric(model$u))
      }
    }
    
    summary_list[[trait]] <- data.frame(
      Environment = environment,
      Trait = trait,
      N_training = length(train_ids),
      N_without_observed_phenotype = nrow(geno) - length(train_ids),
      Status = status
    )
  }
  
  # Save all predictions together; no Excel merging needed
  write.csv(
    data.frame(
      Sample_ID = rownames(geno),
      gebv,
      row.names = NULL,
      check.names = FALSE
    ),
    file.path(
      output_folder,
      paste0("GEBVs_all_", nrow(geno), "_", environment, ".csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    do.call(rbind, summary_list),
    file.path(
      output_folder,
      paste0("Training_summary_", environment, ".csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    do.call(rbind, training_list),
    file.path(
      output_folder,
      paste0("Training_IDs_by_trait_", environment, ".csv")
    ),
    row.names = FALSE
  )
  
  return(invisible(gebv))
}

####################################################################
# 5. Run all three environments
####################################################################

for (environment in names(phenotypes)) {
  
  run_rrBLUP_all_traits(
    geno = geno,
    pheno = phenotypes[[environment]],
    environment = environment,
    output_folder = output_folder
  )
}

#########################################################################################################
# Step 3: Model Prediction Accuracies — rrBLUP
#########################################################################################################
library(rrBLUP)
#folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_core"
folder <- "/home/ahmccorm/scratch/world_veg/global_n7843_PA"
####################
# STEP ONE
####################

# Use the BEAGLE genotype matrix already created by your GS code
geno_data <- geno

# Load phenotype data
control_pheno_check <- read.csv(
  file.path(folder, "control_pheno_filtered_286.csv"),
  row.names = 1
)

heat_stress_1_check <- read.csv(
  file.path(folder, "heat_stress_1_filtered_286.csv"),
  row.names = 1
)

heat_stress_2_check <- read.csv(
  file.path(folder, "heat_stress_2_filtered_286.csv"),
  row.names = 1
)

# Identify common genotypes
common_genotypes <- intersect(
  rownames(geno_data),
  rownames(control_pheno_check)
)

geno_filtered <- geno_data[common_genotypes, , drop = FALSE]

geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
# STEP TWO
####################

# Use all matching phenotyped accessions 
train_genotypes <- common_genotypes

geno_train <- geno_matrix[train_genotypes, , drop = FALSE]

control_train <- control_pheno_check[train_genotypes, , drop = FALSE]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, , drop = FALSE]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, , drop = FALSE]

####################
# STEP THREE
####################

# Use your original function unchanged
source(file.path(folder, "xval_kfold_functions.R"))

# Run the same 10-fold cross-validation × 50 repetitions
xval_k10_rrblup <- k.xval(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_k10_rrblup,
  file.path(folder, "xval_rrblup_kfold_10.RData")
)


####################################################################
# Method Two: Gaussian Kernel
####################################################################

# Separate phenotype object so the rrBLUP input stays unchanged
control_train_kernel <- data.frame(
  Genotype = rownames(control_train),
  control_train,
  row.names = rownames(control_train),
  check.names = FALSE
)

stopifnot(
  identical(rownames(geno_train), rownames(control_train_kernel))
)

# A.mat expects -1/0/1 genotype coding
K <- A.mat(geno_train - 1)
k_dist <- dist(K)

# Run your original Gaussian cross-validation function
xval_k10_GAUSS <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = control_train_kernel,
  y.trainset = control_train_kernel,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_k10_GAUSS,
  file.path(folder, "xval_GAUSS_kfold_10.RData")
)

####################################################################
# Method Three: Exponential Kernel — original implementation
####################################################################

# Keep your original exponential kernel construction
K.Exp <- Kernel_computation(
  X = geno_train,
  name = "exponential",
  degree = NULL,
  nL = NULL
)

exp_dist <- dist(K.Exp)

exp_dist_mat <- as.matrix(exp_dist)

rownames(exp_dist_mat) <- colnames(exp_dist_mat) <-
  rownames(geno_train)

stopifnot(
  identical(
    rownames(exp_dist_mat),
    control_train_kernel$Genotype
  )
)

# Run your original exponential cross-validation function
xval_k10_EXP <- k.xval.EXP(
  g.in = geno_train,
  y.in = control_train_kernel,
  y.trainset = control_train_kernel,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)

saveRDS(
  xval_k10_EXP,
  file.path(folder, "xval_EXP_kfold_10.RData")
)


#####################################################################################################################
# Supplemental Figures S16/S17/18
#####################################################################################################################
# PLOTs
# PAs for control (global)
# PAs for heat-stress #1 (global)
# PAs for heat-stress #2 (global)

####################################################################
# Global collection: PA plots for all three environments
####################################################################

library(ggplot2)
library(dplyr)

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/GS_global/PA"

# Control filenames have no environment suffix
environments <- c(
  Control = "",
  Heat_Stress_1 = "_Heat_Stress_1",
  Heat_Stress_2 = "_Heat_Stress_2"
)

environment_labels <- c(
  Control = "Control",
  Heat_Stress_1 = "Heat Stress 1",
  Heat_Stress_2 = "Heat Stress 2"
)

model_labels <- c(
  rrblup = "rrBLUP",
  GAUSS = "Gaussian Kernel",
  EXP = "Exponential Kernel"
)

plots <- list()
pa_tables <- list()

for (environment in names(environments)) {
  
  suffix <- environments[[environment]]
  
  ##################################################################
  # Load all three models
  ##################################################################
  
  model_results <- lapply(names(model_labels), function(method) {
    
    filename <- paste0(
      "xval_", method, suffix, "_kfold_10.RData"
    )
    
    results <- readRDS(
      file.path(folder, filename)
    )$xval.result
    
    results$model <- model_labels[[method]]
    
    return(results)
  })
  
  all_gs_traits <- bind_rows(model_results) %>%
    mutate(
      r.mean = as.numeric(as.character(r.mean)),
      r.sd = as.numeric(as.character(r.sd)),
      model = factor(
        model,
        levels = unname(model_labels)
      )
    ) %>%
    filter(
      !is.na(trait),
      is.finite(r.mean),
      is.finite(r.sd)
    )
  
  pa_tables[[environment]] <- all_gs_traits
  
  ##################################################################
  # Plot using your previous layout
  ##################################################################
  
  p <- ggplot(
    all_gs_traits,
    aes(y = r.mean, x = model, color = model)
  ) +
    theme_bw() +
    geom_errorbar(
      aes(
        ymin = r.mean - r.sd,
        ymax = r.mean + r.sd
      ),
      width = 0.3,
      position = position_dodge(0.2)
    ) +
    geom_point(size = 3) +
    facet_wrap(
      vars(trait),
      scales = "free_x",
      nrow = 5
    ) +
    geom_hline(
      yintercept = 0.5,
      color = "red",
      linewidth = 1.5,
      linetype = "longdash"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 15)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0(
        "Prediction Accuracy (r) for ",
        n_distinct(all_gs_traits$trait),
        " Traits — ",
        environment_labels[[environment]]
      ),
      x = "Model",
      y = "Prediction Accuracy (r)"
    )
  
  print(p)
  plots[[environment]] <- p
  
  ggsave(
    filename = file.path(
      folder,
      paste0("Global_PA_", environment, ".pdf")
    ),
    plot = p,
    width = 24,
    height = 13,
    units = "in"
  )
}






#####################################################################################################################
# Supplemental Figure S19
#####################################################################################################################
####################################################################
# rrBLUP traits with PA > 0.5 in all three conditions
####################################################################

# Combine the three environment tables created above
rr_only <- bind_rows(pa_tables, .id = "condition") %>%
  filter(model == "rrBLUP") %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Control", "Heat_Stress_1", "Heat_Stress_2"),
      labels = c("Control", "Heat Stress 1", "Heat Stress 2")
    )
  )

# Require a valid result from every condition
rr_traits_all3 <- rr_only %>%
  group_by(trait) %>%
  summarise(
    n_conditions = n_distinct(condition),
    all_above = all(is.finite(r.mean) & r.mean > 0.5),
    .groups = "drop"
  ) %>%
  filter(n_conditions == 3, all_above) %>%
  pull(trait)

print(rr_traits_all3)
cat("Number of qualifying traits:", length(rr_traits_all3), "\n")

rr_filtered <- rr_only %>%
  filter(trait %in% rr_traits_all3)

if (nrow(rr_filtered) == 0) {
  stop("No rrBLUP traits have PA > 0.5 in all three conditions.")
}

p_rr <- ggplot(
  rr_filtered,
  aes(x = condition, y = r.mean)
) +
  theme_bw() +
  geom_hline(
    yintercept = 0.5,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_errorbar(
    aes(
      ymin = r.mean - r.sd,
      ymax = r.mean + r.sd
    ),
    width = 0.25,
    linewidth = 0.8
  ) +
  geom_point(size = 3, color = "#1F77B4") +
  facet_wrap(~ trait, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste0(
      "rrBLUP: ",
      length(rr_traits_all3),
      " Traits with PA > 0.5 in All Three Conditions"
    ),
    x = "Condition",
    y = "Prediction Accuracy (r)"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_rr)

ggsave(
  filename = file.path(
    folder,
    "Global_rrBLUP_traits_all3_above_0.5.pdf"
  ),
  plot = p_rr,
  width = 18,
  height = 10,
  units = "in"
)




#####################################################################################################################
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
# Supplemental Figure S21
#####################################################################################################################
# Quality traits global collection
library(vcfR)
library(rrBLUP)

####################################################################
# 1. Paths
####################################################################

folder <- "/home/ahmccorm/scratch/world_veg/quality/global_n7843_quality"

vcf_file <- "pepper_10000_annuum_beagle.vcf.gz"
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


##############
#PLOT
##############
library(ggplot2)
library(dplyr)

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/quality/global"

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
all_models_global <- bind_rows(
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
p_global <- ggplot(
  all_models_global,
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
    title = "Prediction Accuracy (r) for Quality Traits — Global Collection",
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  coord_cartesian(ylim = c(0, 1))

print(p_global)

ggsave(
  filename = file.path(folder, "prediction_accuracy_23traits_global.pdf"),
  plot = p_global,
  width = 20,
  height = 12,
  units = "in"
)

#above 0.5
bind_rows(
  Core = all_models_core,
  Global = all_models_global,
  .id = "Collection"
) %>%
  group_by(Collection, model) %>%
  summarise(
    Traits_above_0.5 = n_distinct(trait[r.mean > 0.5]),
    Total_traits = n_distinct(trait),
    .groups = "drop"
  ) %>%
  print()





#####################################################################################################################
# Supplemental Figure S22
#####################################################################################################################
# Relationship between climate and quality trait GEBV distributions in core and global Capsicum collections 

###################################################
# Figure A – Climate vs Quality GEBV correlations
###################################################
library(tidyverse)
library(pheatmap)
library(ggfortify)

folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/S21"

# Read new files, rename accession IDs, and remove Set if present
load_gebv <- function(filename) {
  read_csv(
    file.path(folder, filename),
    show_col_types = FALSE
  ) %>%
    rename(Line = Sample_ID) %>%
    mutate(Line = trimws(as.character(Line))) %>%
    select(-any_of("Set"))
}

gebv_control <- load_gebv("GEBVs_all_397_Control.csv")
gebv_heat1   <- load_gebv("GEBVs_all_397_Heat_Stress_1.csv")
gebv_heat2   <- load_gebv("GEBVs_all_397_Heat_Stress_2.csv")

gebv_quality <- load_gebv("GEBVs_Quality_all_397.csv")

gebv_quality_long <- gebv_quality %>%
  pivot_longer(
    cols = -Line,
    names_to = "QualityTrait",
    values_to = "QualityGEBV"
  ) %>%
  mutate(QualityGEBV = as.numeric(QualityGEBV))




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



####################################################################
# Figure B: Quality GEBV distributions
# Traits with rrBLUP PA > 0.5 in BOTH core and global collections
####################################################################

library(tidyverse)

####################################################################
# 1. Paths
####################################################################

base_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1"

folder <- file.path(base_folder, "S21")

####################################################################
# 2. Load new quality GEBVs
####################################################################

load_quality_csv <- function(filename, group_name) {
  readr::read_csv(
    file.path(folder, filename),
    show_col_types = FALSE
  ) %>%
    dplyr::rename(Line = Sample_ID) %>%
    dplyr::select(-tidyselect::any_of(c("Set", "Group"))) %>%
    dplyr::mutate(
      Line = trimws(as.character(Line)),
      Group = group_name
    )
}
gebv_quality_core <- load_quality_csv(
  "GEBVs_Quality_all_397.csv",
  "Core"
)

gebv_quality_global <- load_quality_csv(
  "GEBVs_Quality_all_7843.csv",
  "Global"
)

gebv_quality_all_long <- bind_rows(
  gebv_quality_core,
  gebv_quality_global
) %>%
  pivot_longer(
    cols = -c(Line, Group),
    names_to = "Trait",
    values_to = "GEBV"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Core", "Global"))
  ) %>%
  filter(is.finite(GEBV))

####################################################################
# 3. Load new rrBLUP prediction accuracies
####################################################################

rrblup_core <- readRDS(
  file.path(
    base_folder,
    "quality", "core",
    "xval_rrblup_23traits_kfold_10.RDS"
  )
)$xval.result %>%
  mutate(Source = "Core")

rrblup_global <- readRDS(
  file.path(
    base_folder,
    "quality", "global",
    "xval_rrblup_23traits_kfold_10.RDS"
  )
)$xval.result %>%
  mutate(Source = "Global")

rrblup_all <- bind_rows(
  rrblup_core,
  rrblup_global
) %>%
  mutate(
    r.mean = as.numeric(as.character(r.mean)),
    trait = as.character(trait)
  )

####################################################################
# 4. Select traits with PA > 0.5 in BOTH collections
####################################################################

core_high_traits <- rrblup_all %>%
  filter(
    Source == "Core",
    !is.na(trait),
    is.finite(r.mean),
    r.mean > 0.5
  ) %>%
  pull(trait)

global_high_traits <- rrblup_all %>%
  filter(
    Source == "Global",
    !is.na(trait),
    is.finite(r.mean),
    r.mean > 0.5
  ) %>%
  pull(trait)

common_high_accuracy_traits <- sort(
  intersect(core_high_traits, global_high_traits)
)

cat(
  "Traits with rrBLUP PA > 0.5 in both collections:",
  length(common_high_accuracy_traits),
  "\n"
)

print(common_high_accuracy_traits)

if (length(common_high_accuracy_traits) == 0) {
  stop("No traits have rrBLUP PA > 0.5 in both collections.")
}

# Check that selected traits exist in both GEBV tables
missing_traits <- setdiff(
  common_high_accuracy_traits,
  intersect(
    names(gebv_quality_core),
    names(gebv_quality_global)
  )
)

if (length(missing_traits) > 0) {
  stop(
    "Selected traits missing from GEBV files: ",
    paste(missing_traits, collapse = ", ")
  )
}

gebv_quality_filtered <- gebv_quality_all_long %>%
  filter(Trait %in% common_high_accuracy_traits) %>%
  mutate(
    Trait = factor(
      Trait,
      levels = common_high_accuracy_traits
    )
  )

####################################################################
# 5. Calculate means and quartiles
####################################################################

quality_summary <- gebv_quality_filtered %>%
  group_by(Trait, Group) %>%
  summarise(
    N = n(),
    Mean = mean(GEBV),
    Q25 = unname(quantile(GEBV, 0.25)),
    Median = median(GEBV),
    Q75 = unname(quantile(GEBV, 0.75)),
    .groups = "drop"
  )

reference_lines <- quality_summary %>%
  pivot_longer(
    cols = c(Mean, Q25, Median, Q75),
    names_to = "Statistic",
    values_to = "Position"
  ) %>%
  mutate(
    Line_type = case_when(
      Statistic == "Mean" ~ "Mean",
      Statistic == "Median" ~ "Median",
      TRUE ~ "25th / 75th percentile"
    ),
    Line_type = factor(
      Line_type,
      levels = c(
        "Mean",
        "Median",
        "25th / 75th percentile"
      )
    )
  )

####################################################################
# 6. Clean facet labels and plot
####################################################################

label_vec <- setNames(
  gsub("_", " ", common_high_accuracy_traits),
  common_high_accuracy_traits
)

# Match ggplot's original default colours for Core and Global
collection_colors <- c(
  Core = "#F8766D",
  Global = "#00BFC4"
)

p <- ggplot(
  gebv_quality_filtered,
  aes(x = GEBV, fill = Group)
) +
  geom_density(alpha = 0.4) +
  
  geom_vline(
    data = reference_lines,
    aes(
      xintercept = Position,
      color = Group,
      linetype = Line_type
    ),
    linewidth = 0.5,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  
  facet_wrap(
    ~ Trait,
    scales = "free",
    ncol = 4,
    labeller = labeller(Trait = label_vec)
  ) +
  
  scale_fill_manual(
    values = collection_colors,
    name = "Collection"
  ) +
  
  scale_color_manual(
    values = collection_colors,
    name = "Collection"
  ) +
  
  scale_linetype_manual(
    values = c(
      "Mean" = "solid",
      "Median" = "dashed",
      "25th / 75th percentile" = "dotted"
    ),
    name = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = "none",
    linetype = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(color = "black")
    )
  ) +
  
  labs(
    title = "",
    x = "GEBV",
    y = "Density"
  )

print(p)

####################################################################
# 7. Save at original dimensions
####################################################################

ggsave(
  filename = file.path(
    folder,
    "S20B_Quality_GEBV_distributions.pdf"
  ),
  plot = p,
  width = 20,
  height = 10,
  units = "in"
)

readr::write_csv(
  quality_summary,
  file.path(
    folder,
    "S20B_Quality_GEBV_mean_quantiles.csv"
  )
)

rrblup_all %>%
  filter(trait %in% common_high_accuracy_traits) %>%
  select(Source, trait, r.mean) %>%
  arrange(trait, Source) %>%
  write_csv(
    file.path(folder, "S20B_Common_high_accuracy_traits.csv")
  )

#####################################################################################################################
# Supplemental Figure S23
#####################################################################################################################
# Agronomic GEBV distributions: core vs global
# Traits with rrBLUP PA > 0.5 in all three conditions in BOTH sets
library(tidyverse)

####################################################################
# 1. Paths
####################################################################

base_folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1"

core_pa_folder <- file.path(base_folder, "GS_core", "PA")
global_pa_folder <- file.path(base_folder, "GS_global", "PA")

core_gs_folder <- file.path(
  base_folder, "GS_core", "GS_all_phenotyped"
)

global_gs_folder <- file.path(
  base_folder, "GS_global", "GS"
)

output_folder <- file.path(base_folder, "S21")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

environments <- c(
  "Control",
  "Heat_Stress_1",
  "Heat_Stress_2"
)

condition_labels <- c(
  Control = "Control",
  Heat_Stress_1 = "Heat Stress 1",
  Heat_Stress_2 = "Heat Stress 2"
)

####################################################################
# 2. Load rrBLUP prediction accuracies
####################################################################

load_pa <- function(folder, collection) {
  
  purrr::map_dfr(environments, function(environment) {
    
    suffix <- if (environment == "Control") {
      ""
    } else {
      paste0("_", environment)
    }
    
    filename <- paste0(
      "xval_rrblup", suffix, "_kfold_10.RData"
    )
    
    path <- file.path(folder, filename)
    
    if (!file.exists(path)) {
      stop("PA file not found: ", path)
    }
    
    readRDS(path)$xval.result %>%
      dplyr::transmute(
        trait = as.character(trait),
        r.mean = as.numeric(as.character(r.mean)),
        set = collection,
        condition = environment
      ) %>%
      dplyr::filter(!is.na(trait), nzchar(trait))
  })
}

all_pa <- dplyr::bind_rows(
  load_pa(core_pa_folder, "Core"),
  load_pa(global_pa_folder, "Global")
)

# Require one PA result per trait, collection, and condition
stopifnot(
  !anyDuplicated(all_pa[c("set", "condition", "trait")])
)

####################################################################
# 3. Select traits passing all three conditions in BOTH collections
####################################################################

qualifying_traits <- all_pa %>%
  dplyr::group_by(set, trait) %>%
  dplyr::summarise(
    passes_all = dplyr::n_distinct(condition) == 3 &&
      all(is.finite(r.mean) & r.mean > 0.5),
    .groups = "drop"
  ) %>%
  dplyr::filter(passes_all)

core_traits <- qualifying_traits %>%
  dplyr::filter(set == "Core") %>%
  dplyr::pull(trait)

global_traits <- qualifying_traits %>%
  dplyr::filter(set == "Global") %>%
  dplyr::pull(trait)

common_traits <- sort(
  intersect(core_traits, global_traits)
)

cat(
  "Core traits passing all three conditions:",
  length(core_traits), "\n"
)

cat(
  "Global traits passing all three conditions:",
  length(global_traits), "\n"
)

cat("Common qualifying traits:", length(common_traits), "\n")

print(common_traits)

if (length(common_traits) == 0) {
  stop(
    "No traits meet the PA threshold in all six ",
    "collection/condition combinations."
  )
}

####################################################################
# 4. Save qualifying trait PAs
####################################################################

readr::write_csv(
  all_pa %>%
    dplyr::filter(trait %in% common_traits) %>%
    dplyr::select(set, condition, trait, r.mean) %>%
    dplyr::arrange(trait, set, condition),
  file.path(
    output_folder,
    "S21_Common_traits_prediction_accuracies.csv"
  )
)

####################################################################
# 5. Load GEBVs and reshape selected traits
####################################################################

load_gebvs <- function(folder, collection, n_accessions) {
  
  purrr::map_dfr(environments, function(environment) {
    
    filename <- paste0(
      "GEBVs_all_", n_accessions, "_", environment, ".csv"
    )
    
    path <- file.path(folder, filename)
    
    if (!file.exists(path)) {
      stop("GEBV file not found: ", path)
    }
    
    data <- readr::read_csv(
      path,
      show_col_types = FALSE
    )
    
    required_columns <- c("Sample_ID", common_traits)
    missing_columns <- setdiff(required_columns, names(data))
    
    if (length(missing_columns) > 0) {
      stop(
        "Columns missing from ", filename, ": ",
        paste(missing_columns, collapse = ", ")
      )
    }
    
    stopifnot(!anyDuplicated(data$Sample_ID))
    
    data %>%
      dplyr::select(
        Sample_ID,
        tidyselect::all_of(common_traits)
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(common_traits),
        names_to = "Trait",
        values_to = "GEBV"
      ) %>%
      dplyr::mutate(
        Condition = condition_labels[[environment]],
        Group = collection
      )
  })
}

gebv_filtered <- dplyr::bind_rows(
  load_gebvs(core_gs_folder, "Core", 397),
  load_gebvs(global_gs_folder, "Global", 7843)
) %>%
  dplyr::filter(is.finite(GEBV)) %>%
  dplyr::mutate(
    Trait = factor(Trait, levels = common_traits),
    Condition = factor(
      Condition,
      levels = c("Control", "Heat Stress 1", "Heat Stress 2")
    ),
    Group = factor(
      Group,
      levels = c("Core", "Global")
    )
  )

####################################################################
####################################################################
# 6. Calculate means and quartiles
####################################################################

distribution_summary <- gebv_filtered %>%
  dplyr::group_by(Trait, Condition, Group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Mean = mean(GEBV),
    Q25 = unname(quantile(GEBV, 0.25)),
    Median = median(GEBV),
    Q75 = unname(quantile(GEBV, 0.75)),
    .groups = "drop"
  )

reference_lines <- distribution_summary %>%
  tidyr::pivot_longer(
    cols = c(Mean, Q25, Median, Q75),
    names_to = "Statistic",
    values_to = "Position"
  ) %>%
  dplyr::mutate(
    Line_type = dplyr::case_when(
      Statistic == "Mean" ~ "Mean",
      Statistic == "Median" ~ "Median",
      TRUE ~ "25th / 75th percentile"
    ),
    Line_type = factor(
      Line_type,
      levels = c(
        "Mean",
        "Median",
        "25th / 75th percentile"
      )
    )
  )

####################################################################
# 7. Plot: one trait per row, three conditions across columns
####################################################################

collection_colors <- c(
  Core = "#F8766D",
  Global = "#00BFC4"
)

p <- ggplot(
  gebv_filtered,
  aes(x = GEBV, fill = Group)
) +
  geom_density(alpha = 0.4) +
  
  geom_vline(
    data = reference_lines,
    aes(
      xintercept = Position,
      color = Group,
      linetype = Line_type
    ),
    linewidth = 0.5,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  
  facet_wrap(
    vars(Trait, Condition),
    scales = "free",
    ncol = 3
  ) +
  
  scale_fill_manual(
    values = collection_colors,
    name = "Collection"
  ) +
  
  scale_color_manual(
    values = collection_colors,
    name = "Collection"
  ) +
  
  scale_linetype_manual(
    values = c(
      "Mean" = "solid",
      "Median" = "dashed",
      "25th / 75th percentile" = "dotted"
    ),
    name = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  
  labs(
    title = "Common High-Accuracy Traits: GEBV Distributions",
    subtitle = paste0(
      length(common_traits),
      " traits with rrBLUP PA > 0.5 in all three conditions ",
      "in both collections"
    ),
    x = "GEBV",
    y = "Density"
  ) +
  
  theme(
    legend.position = "top",
    legend.box = "vertical",
    strip.text = element_text(size = 10)
  ) +
  
  guides(
    fill = guide_legend(order = 1),
    color = "none",
    linetype = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(color = "black")
    )
  )

print(p)

####################################################################
# 8. Save using original PDF dimensions
####################################################################

ggsave(
  filename = file.path(
    output_folder,
    "S21_Agronomic_GEBV_density_high_accuracy_traits.pdf"
  ),
  plot = p,
  width = 18,
  height = max(8, 1.5 * length(common_traits)),
  units = "in",
  limitsize = FALSE
)

readr::write_csv(
  distribution_summary,
  file.path(
    output_folder,
    "S21_Agronomic_GEBV_mean_quantiles.csv"
  )
)

#####################################################################################################################
# Supplemental Figure S24
#####################################################################################################################
# genomic heritability
####################################################################

setwd(
  "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/heritability/"
)

library(readxl)
library(dplyr)
library(vcfR)
library(sommer)

####################################################################
# 1. Settings
####################################################################

vcf_file <- "Pepper_397_annuum_LD_pruned_beagle.vcf.gz"

control_file <- "Control.xlsx"
heat1_file <- "Heat_stress_1.xlsx"
heat2_file <- "Heat_stress_2.xlsx"
sheet_name <- "data"

min_maf <- 0.05
max_missing <- 0.10

####################################################################
# 2. Read and standardize phenotypes
####################################################################

read_pheno <- function(path, treatment) {
  
  df <- read_excel(path, sheet = sheet_name)
  names(df) <- tolower(names(df))
  
  if (!("g2pname" %in% names(df))) {
    stop("No g2pname column in ", path)
  }
  
  df$g2pname <- trimws(as.character(df$g2pname))
  
  df <- df[
    !is.na(df$g2pname) & nzchar(df$g2pname),
    ,
    drop = FALSE
  ]
  
  exclude <- c(
    "g2pname", "rep", "genotype", "plantno.", "treatment"
  )
  
  for (cn in names(df)) {
    if (!(cn %in% exclude)) {
      df[[cn]] <- suppressWarnings(
        as.numeric(as.character(df[[cn]]))
      )
    }
  }
  
  df$treatment <- treatment
  df
}

control <- read_pheno(control_file, "Control")
hs1 <- read_pheno(heat1_file, "HS1")
hs2 <- read_pheno(heat2_file, "HS2")

all_pheno <- bind_rows(control, hs1, hs2)

####################################################################
# 3. Read VCF and convert genotypes to 0/1/2
####################################################################

vcf <- read.vcfR(vcf_file, verbose = FALSE)

gt <- extract.gt(vcf, element = "GT")
gt[] <- gsub("|", "/", gt, fixed = TRUE)

gt_numeric <- matrix(
  NA_real_,
  nrow = nrow(gt),
  ncol = ncol(gt),
  dimnames = dimnames(gt)
)

gt_numeric[which(gt == "0/0")] <- 0
gt_numeric[which(gt == "0/1" | gt == "1/0")] <- 1
gt_numeric[which(gt == "1/1")] <- 2

colnames(gt_numeric) <- trimws(colnames(gt_numeric))

stopifnot(!anyDuplicated(colnames(gt_numeric)))

rm(vcf, gt)
gc()

####################################################################
# 4. Filter SNPs by missingness and MAF
####################################################################

missing_rate <- rowMeans(is.na(gt_numeric))

gt_numeric <- gt_numeric[
  missing_rate <= max_missing,
  ,
  drop = FALSE
]

# Allele frequency calculated from observed calls
allele_frequency <- rowMeans(gt_numeric, na.rm = TRUE) / 2
maf <- pmin(allele_frequency, 1 - allele_frequency)

keep_maf <- is.finite(maf) & maf >= min_maf

gt_numeric <- gt_numeric[keep_maf, , drop = FALSE]

if (nrow(gt_numeric) == 0) {
  stop("No SNPs remain after filtering.")
}

cat("SNPs retained after missingness/MAF filters:",
    nrow(gt_numeric), "\n")

# Mean-impute any remaining missing calls
row_means <- rowMeans(gt_numeric, na.rm = TRUE)
missing_indices <- which(is.na(gt_numeric), arr.ind = TRUE)

if (nrow(missing_indices) > 0) {
  gt_numeric[missing_indices] <- row_means[
    missing_indices[, 1]
  ]
}

####################################################################
# 5. Remove invalid/invariant markers and build a finite GRM
####################################################################

# Individuals as rows; SNPs as columns
M <- t(gt_numeric)

marker_sd <- apply(M, 2, sd)

keep <- colSums(!is.finite(M)) == 0 &
  is.finite(marker_sd) &
  marker_sd > 0

cat("SNPs removed before scaling:", sum(!keep), "\n")
cat("SNPs retained for GRM:", sum(keep), "\n")

if (!any(keep)) {
  stop("No valid variable SNPs remain.")
}

M <- M[, keep, drop = FALSE]

M_scaled <- scale(M, center = TRUE, scale = TRUE)

stopifnot(all(is.finite(M_scaled)))

grm <- tcrossprod(M_scaled) / ncol(M_scaled)

# Enforce symmetry and add a small diagonal adjustment
grm <- (grm + t(grm)) / 2
diag(grm) <- diag(grm) + 1e-6

stopifnot(
  all(is.finite(grm)),
  identical(rownames(grm), colnames(grm))
)

cat("GRM accessions:", nrow(grm), "\n")
cat("Non-finite GRM entries:", sum(!is.finite(grm)), "\n")

rm(gt_numeric, M, M_scaled)
gc()

####################################################################
# 6. Heritability function
####################################################################

calc_h2 <- function(pheno_df, trait, grm) {
  
  dat <- data.frame(
    id = trimws(as.character(pheno_df$g2pname)),
    y = as.numeric(pheno_df[[trait]])
  )
  
  dat <- dat[
    is.finite(dat$y) & dat$id %in% rownames(grm),
    ,
    drop = FALSE
  ]
  
  common <- sort(unique(dat$id))
  
  if (length(common) < 10 || nrow(dat) < 2) {
    return(NA_real_)
  }
  
  if (var(dat$y) == 0) {
    return(NA_real_)
  }
  
  dat$id <- factor(dat$id, levels = common)
  grm_sub <- grm[common, common, drop = FALSE]
  
  fit <- tryCatch(
    mmer(
      y ~ 1,
      random = ~ vs(id, Gu = grm_sub),
      rcov = ~ units,
      data = dat,
      verbose = FALSE
    ),
    error = function(e) {
      message(
        "Failed: ", trait, " — ", conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(fit)) {
    return(NA_real_)
  }
  
  # Reject fits explicitly reported as non-converged
  if (isFALSE(fit$convergence)) {
    message("Model did not converge for: ", trait)
    return(NA_real_)
  }
  
  vc <- summary(fit)$varcomp
  
  g_row <- grep("id", rownames(vc))
  e_row <- grep("units", rownames(vc))
  
  if (
    length(g_row) != 1 ||
    length(e_row) != 1 ||
    !("VarComp" %in% colnames(vc))
  ) {
    message("Cannot identify variance components for: ", trait)
    return(NA_real_)
  }
  
  Vg <- as.numeric(vc[g_row, "VarComp"])
  Ve <- as.numeric(vc[e_row, "VarComp"])
  
  if (
    !all(is.finite(c(Vg, Ve))) ||
    Vg < 0 ||
    Ve < 0 ||
    Vg + Ve <= 0
  ) {
    return(NA_real_)
  }
  
  Vg / (Vg + Ve)
}

####################################################################
# 7. Identify traits and run each treatment
####################################################################

exclude <- c(
  "g2pname", "rep", "genotype", "plantno.", "treatment"
)

trait_cols <- names(all_pheno)[
  vapply(all_pheno, is.numeric, logical(1))
]

trait_cols <- setdiff(trait_cols, exclude)

cat("Traits to analyse:", length(trait_cols), "\n")

run_treatment <- function(df, trt) {
  
  out <- vapply(
    trait_cols,
    function(trait) {
      cat("Fitting:", trt, "—", trait, "\n")
      calc_h2(df, trait, grm)
    },
    numeric(1)
  )
  
  data.frame(
    trait = trait_cols,
    treatment = trt,
    h2 = as.numeric(out)
  )
}

res_control <- run_treatment(
  filter(all_pheno, treatment == "Control"),
  "Control"
)

res_hs1 <- run_treatment(
  filter(all_pheno, treatment == "HS1"),
  "HS1"
)

res_hs2 <- run_treatment(
  filter(all_pheno, treatment == "HS2"),
  "HS2"
)

res_all <- bind_rows(res_control, res_hs1, res_hs2)

####################################################################
# 8. Check and save estimates
####################################################################

res_all %>%
  group_by(treatment) %>%
  summarise(
    Estimated = sum(is.finite(h2)),
    Missing = sum(!is.finite(h2)),
    .groups = "drop"
  ) %>%
  print()

write.csv(
  res_all,
  "Figure_S22_heritability_results.csv",
  row.names = FALSE
)

plot_data <- res_all %>%
  filter(is.finite(h2)) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c("Control", "HS1", "HS2")
    )
  )

if (nrow(plot_data) == 0) {
  stop("No finite heritability estimates available for plotting.")
}

####################################################################
# 9. Plot and save two-panel PDF
####################################################################

draw_heritability <- function() {
  
  par(mfrow = c(1, 2))
  
  hist(
    plot_data$h2,
    breaks = seq(0, 1, by = 0.05),
    main = "A. Heritability across treatments",
    xlab = expression(h^2),
    xlim = c(0, 1)
  )
  
  abline(v = 0.25, col = "red", lwd = 2, lty = 2)
  
  boxplot(
    h2 ~ treatment,
    data = plot_data,
    main = "B. Heritability by treatment",
    ylab = expression(h^2),
    xlab = "",
    ylim = c(0, 1)
  )
  
  abline(h = 0.25, col = "red", lwd = 2, lty = 2)
}

# Display plots
old_par <- par(no.readonly = TRUE)
draw_heritability()
par(old_par)

# Save PDF
pdf("Figure_S22_heritability.pdf", width = 10, height = 5)
draw_heritability()
dev.off()



#####################################################################################################################
# Supplemental Figure S25
#####################################################################################################################
# Per-trait genomic heritability
####################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

folder <- paste0(
  "/Users/annamccormick/R/World_Veg_Project/",
  "Reviewer_comments_round1/heritability"
)

####################################################################
# 1. Read estimates
####################################################################

heritability <- read.csv(
  file.path(folder, "Figure_S22_heritability_results.csv"),
  check.names = FALSE
) %>%
  dplyr::mutate(
    trait = as.character(trait),
    h2 = as.numeric(as.character(h2)),
    treatment = factor(
      treatment,
      levels = c("Control", "HS1", "HS2"),
      labels = c("Control", "Heat Stress 1", "Heat Stress 2")
    )
  )

stopifnot(
  !anyDuplicated(heritability[c("trait", "treatment")]),
  !anyNA(heritability$treatment)
)

# Retain missing estimates in the exported table
# but exclude them from plotting
plot_data <- heritability %>%
  dplyr::filter(is.finite(h2))

if (nrow(plot_data) == 0) {
  stop("No finite heritability estimates.")
}

if (any(plot_data$h2 < 0 | plot_data$h2 > 1)) {
  stop("Some heritability estimates fall outside 0–1.")
}

####################################################################
# 2. Order traits by mean heritability
####################################################################

trait_order <- plot_data %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(
    mean_h2 = mean(h2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_h2, trait) %>%
  dplyr::pull(trait)

plot_data <- plot_data %>%
  dplyr::mutate(
    trait = factor(trait, levels = trait_order)
  )

####################################################################
# 3. Plot each trait in each environment
####################################################################

p_h2 <- ggplot(
  plot_data,
  aes(x = h2, y = trait)
) +
  geom_point(
    size = 2.4,
    color = "#4E79A7"
  ) +
  facet_wrap(
    ~ treatment,
    nrow = 1,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.3
    ),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(
      size = 12,
      face = "bold"
    )
  ) +
  labs(
    x = expression(
      "Observation-level genomic heritability (" * h^2 * ")"
    ),
    y = NULL,
    caption = paste(
      "Traits ordered by mean heritability across environments.",
      "Missing estimates are omitted."
    )
  )

print(p_h2)

ggsave(
  filename = file.path(
    folder,
    "Figure_S22_Per_trait_heritability.pdf"
  ),
  plot = p_h2,
  width = 13,
  height = max(8, 2 + 0.22 * length(trait_order)),
  units = "in"
)

####################################################################
# 4. Supplementary table: one row per trait
####################################################################

h2_table <- heritability %>%
  dplyr::select(trait, treatment, h2) %>%
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = h2,
    names_prefix = "h2_"
  ) %>%
  dplyr::arrange(trait)

write.csv(
  h2_table,
  file.path(folder, "Table_Per_trait_genomic_heritability.csv"),
  row.names = FALSE
)

####################################################################
# 5. Summaries for the main Results text
####################################################################

h2_summary <- plot_data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(
    N_traits = dplyr::n(),
    Minimum = min(h2),
    Median = median(h2),
    Maximum = max(h2),
    .groups = "drop"
  )

print(h2_summary)

write.csv(
  h2_summary,
  file.path(folder, "Heritability_summary_by_environment.csv"),
  row.names = FALSE
)


#####################################################################################################################
# Supplemental Figure S26
#####################################################################################################################
#LD decay on HPC
####################################################################
# LD decay: unpruned C. annuum core, n = 397
# Input: quality-filtered SNPs sampled across 12 chromosomes
# Calculate genotype-correlation r² for pairs within 1 Mb
####################################################################

library(ggplot2)
library(dplyr)

folder <- paste0(
  "/home/ahmccorm/kantar_koastore/anna/",
  "World_Veg_Collab_Pepper/Filtered_vcf/LD_core397"
)

setwd(folder)

output_prefix <- "LD_core397_unpruned"
max_distance <- 1000000
max_pairs <- 100000
min_individuals <- 10

set.seed(12345)

out_file <- function(suffix) {
  paste0(output_prefix, suffix)
}

####################################################################
# 1. Read sampled genotypes
####################################################################

sample_ids <- trimws(readLines("sample_ids.txt"))

stopifnot(
  length(sample_ids) == 397,
  !anyDuplicated(sample_ids)
)

cat("Reading sampled genotype table...\n")

sampled <- read.delim(
  "sampled_genotypes.tsv",
  header = FALSE,
  colClasses = "character",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(ncol(sampled) == length(sample_ids) + 2)

chromosome <- sampled[[1]]
position <- as.numeric(sampled[[2]])

stopifnot(all(is.finite(position)))

gt <- as.matrix(sampled[, -(1:2), drop = FALSE])
gt[] <- gsub("|", "/", gt, fixed = TRUE)

# Convert diploid calls to allele counts
geno <- matrix(
  NA_integer_,
  nrow = nrow(gt),
  ncol = ncol(gt),
  dimnames = list(NULL, sample_ids)
)

geno[which(gt == "0/0")] <- 0L
geno[which(gt == "0/1" | gt == "1/0")] <- 1L
geno[which(gt == "1/1")] <- 2L

# Missing calls remain NA; no genotype imputation
cat("Samples:", ncol(geno), "\n")
cat("Sampled SNPs:", nrow(geno), "\n")
cat("Missing/unsupported calls:", sum(is.na(geno)), "\n")

rm(sampled, gt)
gc()

####################################################################
# 2. Pairwise r² using jointly observed accessions
####################################################################

calculate_r2 <- function(x, y) {
  
  valid <- is.finite(x) & is.finite(y)
  
  if (sum(valid) < min_individuals) {
    return(NA_real_)
  }
  
  x <- x[valid]
  y <- y[valid]
  
  if (sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  
  r2 <- cor(x, y)^2
  
  if (!is.finite(r2)) {
    return(NA_real_)
  }
  
  pmin(1, pmax(0, r2))
}

####################################################################
# 3. Uniformly sample eligible SNP pairs per chromosome
####################################################################

calculate_chromosome <- function(chrom) {
  
  cat("\nProcessing:", chrom, "\n")
  
  idx <- which(chromosome == chrom)
  idx <- idx[order(position[idx])]
  
  pos <- position[idx]
  g <- geno[idx, , drop = FALSE]
  
  if (length(pos) < 2) {
    return(NULL)
  }
  
  # Partners must have a greater position and be within 1 Mb
  first_partner <- findInterval(pos, pos) + 1L
  last_partner <- findInterval(pos + max_distance, pos)
  
  n_partners <- pmax(
    0L,
    last_partner - first_partner + 1L
  )
  
  cumulative <- cumsum(as.double(n_partners))
  total_pairs <- sum(as.double(n_partners))
  
  cat("SNPs:", length(pos), "\n")
  cat("Eligible pairs:", total_pairs, "\n")
  
  if (total_pairs == 0) {
    return(NULL)
  }
  
  pair_ranks <- sample.int(
    total_pairs,
    min(max_pairs, total_pairs)
  )
  
  snp1 <- findInterval(pair_ranks - 1, cumulative) + 1L
  previous_total <- c(0, cumulative)[snp1]
  
  snp2 <- first_partner[snp1] +
    pair_ranks - previous_total - 1L
  
  r2_values <- vapply(
    seq_along(pair_ranks),
    function(i) {
      calculate_r2(g[snp1[i], ], g[snp2[i], ])
    },
    numeric(1)
  )
  
  result <- data.frame(
    chromosome = chrom,
    distance = pos[snp2] - pos[snp1],
    r2 = r2_values,
    snp1_pos = pos[snp1],
    snp2_pos = pos[snp2]
  )
  
  result <- result[is.finite(result$r2), , drop = FALSE]
  
  cat("Valid sampled pairs:", nrow(result), "\n")
  
  # Save each chromosome as it finishes
  write.csv(
    result,
    out_file(paste0("_chr", make.names(chrom), "_pairs.csv")),
    row.names = FALSE
  )
  
  result
}

all_ld_results <- dplyr::bind_rows(
  lapply(unique(chromosome), calculate_chromosome)
)

if (nrow(all_ld_results) == 0) {
  stop("No valid SNP-pair correlations were calculated.")
}

write.csv(
  all_ld_results,
  out_file("_raw_results.csv"),
  row.names = FALSE
)

rm(geno)
gc()

####################################################################
# 4. Original distance bins
####################################################################

distance_breaks <- c(
  seq(0, 5000, by = 250),
  seq(5500, 10000, by = 500),
  seq(11000, 20000, by = 1000),
  seq(22000, 30000, by = 2000),
  seq(32500, 50000, by = 2500),
  seq(55000, 100000, by = 5000),
  seq(110000, 200000, by = 10000),
  seq(220000, 400000, by = 20000),
  seq(425000, 600000, by = 25000),
  seq(650000, 1000000, by = 50000)
)

all_ld_results$distance_bin <- cut(
  all_ld_results$distance,
  breaks = distance_breaks,
  include.lowest = TRUE
)

ld_summary <- all_ld_results %>%
  dplyr::filter(!is.na(distance_bin)) %>%
  dplyr::group_by(chromosome, distance_bin) %>%
  dplyr::summarise(
    mean_r2 = mean(r2),
    median_r2 = median(r2),
    sd_r2 = sd(r2),
    n_pairs = dplyr::n(),
    mean_distance = mean(distance),
    .groups = "drop"
  )

genome_wide_summary <- all_ld_results %>%
  dplyr::filter(!is.na(distance_bin)) %>%
  dplyr::group_by(distance_bin) %>%
  dplyr::summarise(
    mean_r2 = mean(r2),
    median_r2 = median(r2),
    sd_r2 = sd(r2),
    n_pairs = dplyr::n(),
    mean_distance = mean(distance),
    .groups = "drop"
  )

write.csv(
  ld_summary,
  out_file("_summary.csv"),
  row.names = FALSE
)

write.csv(
  genome_wide_summary,
  out_file("_genome_wide_summary.csv"),
  row.names = FALSE
)

####################################################################
# 5. Fit exponential decay models
# Distance is fitted in kb for numerical stability
####################################################################

fit_decay <- function(df) {
  
  df <- df %>%
    dplyr::filter(
      is.finite(mean_r2),
      is.finite(mean_distance),
      mean_r2 > 0
    ) %>%
    dplyr::mutate(distance_kb = mean_distance / 1000)
  
  if (nrow(df) < 4) {
    return(NULL)
  }
  
  fit <- tryCatch(
    nls(
      mean_r2 ~ a * exp(-b * distance_kb) + c,
      data = df,
      start = list(
        a = max(df$mean_r2) - min(df$mean_r2),
        b = 0.01,
        c = min(df$mean_r2)
      ),
      algorithm = "port",
      lower = c(a = 0, b = 1e-10, c = 0),
      control = nls.control(maxiter = 500)
    ),
    error = function(e) {
      message("Decay fit failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  p <- coef(fit)
  
  if (any(!is.finite(p)) || p["a"] <= 0) {
    return(NULL)
  }
  
  half_kb <- unname(log(2) / p["b"])
  
  threshold_kb <- NA_real_
  
  if (p["c"] < 0.2 && p["a"] + p["c"] > 0.2) {
    threshold_kb <- unname(
      -log((0.2 - p["c"]) / p["a"]) / p["b"]
    )
  }
  
  # Avoid reporting crossings outside the observed bin-mean range
  observed_range <- range(df$distance_kb)
  
  within_range <- function(x) {
    is.finite(x) &&
      x >= observed_range[1] &&
      x <= observed_range[2]
  }
  
  data.frame(
    a = unname(p["a"]),
    b_per_kb = unname(p["b"]),
    c = unname(p["c"]),
    half_decay_distance_bp = if (within_range(half_kb)) {
      half_kb * 1000
    } else {
      NA_real_
    },
    r2_0.2_distance_bp = if (within_range(threshold_kb)) {
      threshold_kb * 1000
    } else {
      NA_real_
    }
  )
}

chromosome_parameters <- dplyr::bind_rows(
  lapply(unique(ld_summary$chromosome), function(chrom) {
    
    fit <- fit_decay(
      ld_summary[ld_summary$chromosome == chrom, , drop = FALSE]
    )
    
    if (is.null(fit)) {
      return(NULL)
    }
    
    data.frame(chromosome = chrom, fit)
  })
)

write.csv(
  chromosome_parameters,
  out_file("_parameters.csv"),
  row.names = FALSE
)

gw_fit <- fit_decay(genome_wide_summary)

####################################################################
# 6. Prepare fitted curve and annotations
####################################################################

curve_data <- data.frame(
  distance_kb = numeric(),
  predicted_r2 = numeric()
)

subtitle_text <- "Exponential decay fit unavailable"

if (!is.null(gw_fit)) {
  
  write.csv(
    gw_fit,
    out_file("_genome_wide_parameters.csv"),
    row.names = FALSE
  )
  
  print(gw_fit)
  
  curve_data <- data.frame(
    distance_kb = exp(seq(
      log(min(genome_wide_summary$mean_distance) / 1000),
      log(max(genome_wide_summary$mean_distance) / 1000),
      length.out = 500
    ))
  )
  
  curve_data$predicted_r2 <- gw_fit$a *
    exp(-gw_fit$b_per_kb * curve_data$distance_kb) +
    gw_fit$c
  
  format_distance <- function(x) {
    if (is.finite(x)) {
      paste0(round(x), " bp")
    } else {
      "not estimated within observed range"
    }
  }
  
  subtitle_text <- paste0(
    "Half-decay above asymptote: ",
    format_distance(gw_fit$half_decay_distance_bp),
    "\nr² = 0.2: ",
    format_distance(gw_fit$r2_0.2_distance_bp)
  )
}

####################################################################
# 7. Plot helper
####################################################################

make_plot <- function(points, raw = FALSE) {
  
  p <- ggplot(
    points,
    aes(x = distance_kb, y = plot_r2)
  ) +
    geom_point(
      color = "darkblue",
      size = if (raw) 1 else 3,
      alpha = if (raw) 0.2 else 0.6
    ) +
    geom_line(
      data = curve_data,
      aes(x = distance_kb, y = predicted_r2),
      inherit.aes = FALSE,
      color = if (raw) "red" else "steelblue",
      linewidth = 1.2
    ) +
    geom_hline(
      yintercept = 0.2,
      color = "orange",
      linetype = "dotted"
    ) +
    scale_x_log10(
      breaks = c(0.1, 1, 10, 100, 1000)
    ) +
    labs(
      title = if (raw) {
        "Genome-wide LD decay: sampled SNP pairs"
      } else {
        "Genome-wide average LD decay"
      },
      subtitle = subtitle_text,
      x = "Distance (kb, log scale)",
      y = if (raw) "r²" else "Mean r²",
      caption = "C. annuum core: 397 accessions; unpruned SNPs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(gw_fit)) {
    
    if (is.finite(gw_fit$half_decay_distance_bp)) {
      p <- p + geom_vline(
        xintercept = gw_fit$half_decay_distance_bp / 1000,
        color = "red",
        linetype = "dashed"
      )
    }
    
    if (is.finite(gw_fit$r2_0.2_distance_bp)) {
      p <- p + geom_vline(
        xintercept = gw_fit$r2_0.2_distance_bp / 1000,
        color = "orange",
        linetype = "dashed"
      )
    }
  }
  
  p
}

####################################################################
# 8. Save binned and raw-pair plots
####################################################################

binned_points <- genome_wide_summary %>%
  dplyr::mutate(
    distance_kb = mean_distance / 1000,
    plot_r2 = mean_r2
  )

p_binned <- make_plot(binned_points)

raw_points <- all_ld_results[
  sample.int(
    nrow(all_ld_results),
    min(5000, nrow(all_ld_results))
  ),
  ,
  drop = FALSE
] %>%
  dplyr::mutate(
    distance_kb = distance / 1000,
    plot_r2 = r2
  )

p_raw <- make_plot(raw_points, raw = TRUE)

# Base PDF plotmath can avoid Unicode font warnings
p_binned <- p_binned + labs(y = expression(Mean~r^2))
p_raw <- p_raw + labs(y = expression(r^2))

ggsave(
  out_file("_genome_wide_binned.pdf"),
  plot = p_binned,
  width = 10,
  height = 8
)

ggsave(
  out_file("_genome_wide_raw.pdf"),
  plot = p_raw,
  width = 10,
  height = 8
)

cat("\nLD analysis completed.\n")

print(
  head(
    genome_wide_summary[
      order(genome_wide_summary$mean_distance),
    ],
    15
  )
)

writeLines(
  capture.output(sessionInfo()),
  out_file("_sessionInfo.txt")
)




##### LD stats
folder <- "/Users/annamccormick/R/World_Veg_Project/Reviewer_comments_round1/LD"

chr <- read.csv(file.path(folder, "LD_core397_unpruned_parameters.csv"))
gw <- read.csv(file.path(folder, "LD_core397_unpruned_genome_wide_parameters.csv"))
gw$chromosome <- "Genome-wide"

fits <- rbind(gw[, names(chr)], chr)

crossing <- ifelse(
  is.finite(fits$r2_0.2_distance_bp),
  sprintf("%.1f", fits$r2_0.2_distance_bp / 1000),
  ifelse(fits$c >= 0.2, "No fitted crossing", "Not estimated within range")
)

ld_table <- data.frame(
  "Region" = fits$chromosome,
  "a" = sprintf("%.4f", fits$a),
  "b (kb^-1)" = sprintf("%.6f", fits$b_per_kb),
  "c" = sprintf("%.4f", fits$c),
  "Half-decay distance (kb)" = sprintf("%.1f", fits$half_decay_distance_bp / 1000),
  "Distance to r² = 0.2 (kb)" = crossing,
  check.names = FALSE
)

write.csv(
  ld_table,
  file.path(folder, "Table_SX_LD_decay.csv"),
  row.names = FALSE
)


#####################################################################################################################
# Supplemental Figure S27
#####################################################################################################################
library(dplyr)
library(pheatmap)
library(grid)

setwd("~/R/World_Veg_Project/anna_pheno_correlations/")
out_png <- "correlation_heatmap_phenotype_treatment_clustered_FIN.png"


# Read files
control <- read.csv("control_pheno_filtered_286.csv", stringsAsFactors = FALSE) %>%
  mutate(treatment = "Control")

hs1 <- read.csv("heat_stress_1_filtered_286.csv", stringsAsFactors = FALSE) %>%
  mutate(treatment = "Heat_Stress_1")

hs2 <- read.csv("heat_stress_2_filtered_286.csv", stringsAsFactors = FALSE) %>%
  mutate(treatment = "Heat_Stress_2")


# Define ID + trait columns

id_col <- "X"  # change if your ID col differs

# columns you do NOT want treated as phenotypes
non_trait <- c(id_col, "treatment", "g2pname", "genotype")

# traits shared across all three datasets
trait_cols <- Reduce(intersect, list(
  setdiff(names(control), non_trait),
  setdiff(names(hs1), non_trait),
  setdiff(names(hs2), non_trait)
))



# Make WIDE tables (trait_treatment columns)

make_wide <- function(df) {
  tr <- unique(df$treatment)
  
  df %>%
    dplyr::select(dplyr::all_of(id_col), dplyr::all_of(trait_cols)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(trait_cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) %>%
    dplyr::rename_with(
      ~ paste0(.x, "_", tr),
      dplyr::all_of(trait_cols)
    )
}

wide_ctrl <- make_wide(control)
wide_hs1  <- make_wide(hs1)
wide_hs2  <- make_wide(hs2)

wide_all <- wide_ctrl %>%
  inner_join(wide_hs1, by = id_col) %>%
  inner_join(wide_hs2, by = id_col)

Xmat <- wide_all %>%
  dplyr::select(-dplyr::all_of(id_col)) %>%
  as.matrix()



# Filter bad columns (missingness + zero variance)

missing_threshold <- 0.5   # 
min_non_na <- 10           # prevents unstable correlations

keep <- apply(Xmat, 2, function(v) {
  finite <- is.finite(v)
  if (mean(!finite) > missing_threshold) return(FALSE)
  if (sum(finite) < min_non_na) return(FALSE)
  if (sd(v[finite]) == 0) return(FALSE)
  TRUE
})

Xmat2 <- Xmat[, keep, drop = FALSE]



# Correlation among phenotype×treatment columns

cor_mat <- cor(Xmat2, use = "pairwise.complete.obs")
cor_mat[!is.finite(cor_mat)] <- 0


#  Symmetric clustering 

d  <- as.dist(1 - cor_mat)
hc <- hclust(d, method = "complete")


# Draw with pheatmap, then save 

ph <- pheatmap(
  cor_mat,
  cluster_rows = hc,
  cluster_cols = hc,
  color  = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = NA,
  fontsize_row = 4,
  fontsize_col = 4,
  angle_col = 90,
  main = "Clustered Correlation Matrix - Phenotype × Treatment Combinations",
  silent = TRUE
)

png(out_png, width = 4048, height = 4048, res = 300)
grid.newpage()
grid.draw(ph$gtable)
dev.off()

message("Saved: ", normalizePath(out_png))

#save as PDF
out_pdf <- "correlation_heatmap_phenotype_treatment_clustered_FIN.pdf"

pdf(out_pdf, width = 12, height = 12)  # dimensions in inches

grid.newpage()
grid.draw(ph$gtable)

dev.off()

message("Saved PDF: ", normalizePath(out_pdf))

#checking trait counts
# Print matrix dimensions
n_traits <- ncol(cor_mat)

cat("Number of phenotype × treatment variables retained:\n")
cat("X-axis:", n_traits, "\n")
cat("Y-axis:", n_traits, "\n\n")

cat("Matrix dimensions:", paste(dim(cor_mat), collapse = " x "), "\n")
#219 x 219



#PDF
library(corrplot)
setwd("~/R/World_Veg_Project/anna_pheno_correlations/")

pdf(file.path("full_correlation_plot_219x219_all_by_all.pdf"),
    width = 15,   # inches
    height = 15)

corrplot(cor_mat,
         method = "color",
         type = "full",
         order = "original",
         tl.cex = 0.2,
         tl.col = "black",
         title = paste0("Complete Correlation Matrix - ",
                        nrow(cor_mat),
                        " Variables\n(73 phenotypic traits × 3 treatments)"),
         mar = c(0, 0, 3, 0))

dev.off()



#####################################################################################################################
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
# Supplemental Figure S29
#####################################################################################################################
# Global STI: BEAGLE genotypes → GS → prediction accuracies
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

folder <- "/home/ahmccorm/scratch/world_veg/STI/global"

output_folder <- file.path(folder, "outputs")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(
  folder, "pepper_10000_annuum_beagle.vcf.gz"
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
    
    cat(
      "Fitting final STI GS model:", trait,
      "- training n =", length(y), "\n"
    )
    
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
# zero variance, or an unsuccessful final GS fit
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
  file.path(
    output_folder,
    "xval_rrblup_STI_global7843_kfold10.RDS"
  )
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
  file.path(
    output_folder,
    "xval_GAUSS_STI_global7843_kfold10.RDS"
  )
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
  file.path(
    output_folder,
    "xval_EXP_STI_global7843_kfold10.RDS"
  )
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
  file.path(output_folder, "STI_PA_summary_global7843.csv"),
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
      "Global STI Prediction Accuracies — ",
      n_distinct(all_models$trait), " Traits"
    ),
    x = "Model",
    y = "Prediction Accuracy (r)"
  ) +
  coord_cartesian(ylim = c(0, 1))

ggsave(
  filename = file.path(
    output_folder,
    "STI_PA_all_traits_global7843.pdf"
  ),
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

cat(
  "\nTraits with PA > 0.5 for all three models:",
  length(high_traits), "\n"
)

write.csv(
  data.frame(Trait = high_traits),
  file.path(output_folder, "STI_traits_all_models_above_0.5.csv"),
  row.names = FALSE
)

cat("\nGlobal STI GS and PA analyses completed.\n")

#####################################################################################################################
#####################################################################################################################
                                                      # FIN #
#####################################################################################################################
#####################################################################################################################




