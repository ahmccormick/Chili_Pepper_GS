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
#Figure S2. Population structure for the core (n=423) capsicum collection

#core population structure/training plots
#SNPRelate
#########
library(SNPRelate)
setwd("~/R/cannabis_GEAV/Outputs/")

vcf.fn <- "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_filtered_and_LD_pruned.vcf.gz"

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
write.csv(tab, "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC.csv")


library(ggplot2)
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC.csv")
pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Pepper_PC.csv")

# Basic PCA plot: EV1 vs. EV2
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2)) +
  geom_point(size = 3, alpha = 0.8) + # Points for samples
  theme_minimal() +                   # Clean theme
  labs(
    title = "Pepper PCA 1",
    x = "EV1 (2.76 %)",
    y = "EV2 (2.05%)"
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
pepper_pc <- read_csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Pepper_PC.csv")

# Load the Global Capsicum Core Collection Excel file (first sheet)
capsicum_core <- read_excel("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Global capsicum core collection.xlsx")

# Print column names to check structure
colnames(pepper_pc)
colnames(capsicum_core)

# Perform an inner join to match based on the sample.id column
matched_data <- pepper_pc %>%
  inner_join(capsicum_core, by = c("sample.id" = "id")) 

# Save the matched data to a new CSV file
#write_csv(matched_data, "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv")
write_csv(matched_data, "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv") 


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
pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Matched_Pepper_Data_423.csv")
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
    x = "EV1 (2.76%)",
    y = "EV2 (2.05%)",
    color = "Species"  # Legend title
  ) +
  scale_color_viridis_d(labels = species_labels)  # Updated legend labels

# Display the plot
print(pca_plot)

#############################################################################
# Count occurrences of each category in 'wild'
wild_counts <- pca %>%
  group_by(wild) %>%
  summarise(n = n())

# Create a named vector for new legend labels
wild_labels <- setNames(
  paste0(wild_counts$wild, " (n=", wild_counts$n, ")"),
  wild_counts$wild
)

# Basic PCA plot with updated legend for 'wild'
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2, color = as.factor(wild))) +
  geom_point(size = 3, alpha = 0.8) +  # Points for samples
  theme_minimal() +  # Clean theme
  labs(
    title = "Column: wild",
    x = "EV1 (2.76%)",
    y = "EV2 (2.05%)",
    color = "Wild Status"  # Legend title
  ) +
  scale_color_viridis_d(labels = wild_labels)  # Updated legend labels

# Display the plot
print(pca_plot)

##############
#Pannel B
#colour by hcpc group in pca
#####################################################################################################################
library(ggplot2)
library(dplyr)

# Load matched metadata with cluster column
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv")
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/pepper_metadata_with_clusters.csv")
pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Matched_Pepper_Data_423.csv")
pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/pepper_metadata_with_clusters.csv")


# Make sure cluster is treated as a factor
pca$cluster <- as.factor(pca$cluster)

# Count samples per cluster for legend
cluster_counts <- pca %>%
  group_by(cluster) %>%
  summarise(n = n())

# Create new legend labels
cluster_labels <- setNames(
  paste0("Cluster ", cluster_counts$cluster, " (n=", cluster_counts$n, ")"),
  cluster_counts$cluster
)

# Plot PCA colored by cluster
ggplot(pca, aes(x = EV1, y = EV2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA of Pepper Accessions Colored by HCPC Cluster",
    x = "EV1 (2.76%)",
    y = "EV2 (2.05%)",
    color = "Cluster"
  ) +
  scale_color_viridis_d(labels = cluster_labels) +  # Uses discrete viridis colors
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


ggplot(pca, aes(x = EV1, y = EV2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", linetype = "solid", size = 1, alpha = 0.6) +  # Ellipses
  theme_minimal() +
  labs(
    title = "",
    x = "EV1 (2.76%)",
    y = "EV2 (2.05%)",
    color = "Cluster"
  ) +
  scale_color_viridis_d(labels = cluster_labels) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


##############
#CORE VISUALISATION
#Pannel C
##############
# Read in the training core set
#core_ids <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
core_ids <- read.csv("~/R/World_Veg_Project/filtered_data_2/Training_Genotypes_n150_run1.csv")$Training_Genotypes
pepper_data <- read.csv("~/R/World_Veg_Project/filtered_data_2/pepper_metadata_with_clusters.csv")

# Create a logical column for Core vs. Non-core
pepper_data$Core <- ifelse(pepper_data$sample.id %in% core_ids, "Core", "Non-core")

# Convert to factor for consistent plotting
pepper_data$Core <- factor(pepper_data$Core, levels = c("Non-core", "Core"))
# Extract PCA coordinates
pca_coords <- res.hcpc$data.clust
pca_coords$sample.id <- rownames(pca_coords)

# Merge with core training label from pepper_data
pca_annotated <- left_join(pca_coords, pepper_data[, c("sample.id", "Core")], by = "sample.id")
ggplot(pca_annotated, aes(x = EV1, y = EV2, color = Core)) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(values = c("grey80", "red3")) +
  labs(
    title = "",
    x = "EV1 (2.76%)",
    y = "EV2 (2.05%)",
    color = "Training Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



#############
#pannel D
#pca_df <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/pepper_metadata_with_clusters.csv")
#hot_28 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/Extreme_V_SampleIDs.csv") %>%
  rename(sample.id = x) %>%
  mutate(Group = "Heat-loving")
#cold_57 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/Extreme_InvertedV_SampleIDs.csv") %>%
  rename(sample.id = x) %>%
  mutate(Group = "Cold-loving")

#############
  library(dplyr)

pca_df <- read.csv("~/R/World_Veg_Project/filtered_data_2/pepper_metadata_with_clusters.csv")

hot_28 <- read_csv("~/R/World_Veg_Project/filtered_data_2/Extreme_V_SampleIDs.csv",
                   col_names = "sample.id") %>%
  mutate(Group = "Heat-loving")

cold_57 <- read_csv("~/R/World_Veg_Project/filtered_data_2/Extreme_InvertedV_SampleIDs.csv",
                    col_names = "sample.id") %>%
  mutate(Group = "Cold-loving")


extreme_lines <- bind_rows(hot_28, cold_57)

pca_df <- left_join(pca_df, extreme_lines, by = "sample.id")


pca_df$Group <- recode(pca_df$Group,
                       "Cold-loving" = "Cool responsive",
                       "Heat-loving" = "Heat responsive")

ggplot(pca_df, aes(x = EV1, y = EV2)) +
  geom_point(aes(color = is.na(Group)), alpha = 0.3, size = 2) +
  geom_point(data = filter(pca_df, !is.na(Group)),
             aes(color = Group), size = 3, alpha = 0.9) +
  scale_color_manual(values = c(
    "Cool responsive" = "blue4",
    "Heat responsive" = "red3",
    `TRUE` = "grey80"
  )) +
  labs(
    title = "",
    x = "EV1 (2.76%)",
    y = "EV2 (2.05%)",
    color = "Group"
  ) +
  theme_minimal(base_size = 14)


##############
# Pepper Dendrogram with Species Bar
##############
#pannel E
# Load necessary libraries
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)

# Load the dataset
pepper_pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Matched_Pepper_Data_423.csv")
#pepper_pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_423.csv")
#pepper_pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/Matched_Pepper_Data_284.csv")

# Select only numeric columns for PCA
pca_data <- pepper_pca %>% select_if(is.numeric)

# Set row names as sample ID if available
row.names(pca_data) <- pepper_pca$sample.id  

# Perform PCA
res.pca <- PCA(pca_data, graph = FALSE)

# Perform HCPC with k = 5
res.hcpc <- HCPC(res.pca, nb.clust = 5, graph = FALSE)

# Plot the dendrogram
fviz_dend(res.hcpc, 
          cex = 0.6,         # Adjust label size
          rect = TRUE,       # Add rectangles around clusters
          rect_fill = TRUE,  # Color the clusters
          rect_border = "black")  # Black border for clarity



# Extract clustering results
pepper_pca$cluster <- res.hcpc$data.clust$clust

# Define species colors
species_colors <- c(
  "Capsicum annuum" = "blue",
  "Capsicum baccatum" = "yellow",
  "Capsicum chacoense" = "brown",
  "Capsicum chinense" = "red",
  "Capsicum frutescens" = "green"
)

# Convert species to factor
pepper_pca$species <- factor(pepper_pca$species, levels = names(species_colors))

# Create a color bar dataframe
color_bar_df <- data.frame(Type = pepper_pca$species, x = 1:nrow(pepper_pca))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed

# Add species color bar WITHOUT cluster rectangles
f <- fviz_dend(res.hcpc, 
               cex = 0.2, 
               rect = FALSE,    # Remove black dashed cluster rectangles
               rect_fill = FALSE,   
               rect_border = "transparent") +  # Ensure no unwanted borders
  geom_tile(data = color_bar_df, aes(x = x, y = bar_position_y, fill = Type), 
            width = 1, height = 0.2) +
  scale_fill_manual(values = species_colors) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"))

# Plot the dendrogram
print(f)


# View how many samples are in each cluster
table(pepper_pca$cluster)

# Cluster 1
cluster1_samples <- pepper_pca %>% filter(cluster == 1) %>% pull(sample.id)

# Cluster 2
cluster2_samples <- pepper_pca %>% filter(cluster == 2) %>% pull(sample.id)

# Cluster 3
cluster3_samples <- pepper_pca %>% filter(cluster == 3) %>% pull(sample.id)

# Cluster 4
cluster4_samples <- pepper_pca %>% filter(cluster == 4) %>% pull(sample.id)

# Cluster 5
cluster5_samples <- pepper_pca %>% filter(cluster == 5) %>% pull(sample.id)


head(pepper_pca[, c("sample.id", "cluster")])


write.csv(pepper_pca, "~/R/World_Veg_Collab_Pepper/outputs/pepper_metadata_with_clusters.csv", row.names = FALSE)

##############
# K optimal?

library("FactoMineR")
library("factoextra")
pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC.csv")
pca2 <- pca[3:6]

head(pca2)
# Standardizing the PCA data before clustering
my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")
a

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")
b



########################################################################################################################
# Supplemental Figure S3
########################################################################################################################
#Figure S3. Kinship matrix heatmap for the core collection (n=423) coloured by species.

setwd("~/R/World_Veg_Project/anna_kinship_matrix/")

library(SNPRelate)
library(pheatmap)

# 2. Define File Names 
vcf.fn <- "~/R/World_Veg_Project/anna_kinship_matrix/Pepper_filtered_and_LD_pruned.vcf.gz"
gds.fn <- "~/R/World_Veg_Project/anna_kinship_matrix/Pepper_filtered_and_LD_pruned.gds"

# 3. Convert VCF to GDS 
snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")

# 4. Open GDS File 
genofile <- snpgdsOpen(gds.fn)

# 5. Calculate Kinship Matrix 
ibs <- snpgdsIBS(genofile, num.thread = 2)
kinship_matrix <- 1 - ibs$ibs
rownames(kinship_matrix) <- colnames(kinship_matrix) <- ibs$sample.id

# 6. Read Metadata and Prepare Annotation 
# Replace with your file name and format (e.g., CSV, tab-delimited)
metadata <- read.csv("Matched_Pepper_Data_423.csv")
metadata <- metadata[match(ibs$sample.id, metadata$sample.id), ]  # ensure same order
rownames(metadata) <- metadata$sample.id
annotation <- data.frame(Population = metadata$species)
rownames(annotation) <- metadata$sample.id

# Optional: set custom colors for groups
group_colors <- unique(metadata$species)
colors <- setNames(rainbow(length(group_colors)), group_colors)
ann_colors <- list(Population = colors)

# Plot Heatmap with Colored Annotations 


pheatmap(kinship_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Kinship Matrix Heatmap by Population",
         labels_row = NA,
         labels_col = NA,
         show_rownames = F, show_colnames = F,
         annotation_names_row = F, annotation_names_col = F)

out_pdf <- "kinship_heatmap_by_population.pdf"

pheatmap(kinship_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Kinship Matrix Heatmap by Population",
         labels_row = NA,
         labels_col = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         filename = out_pdf,
         width = 12,
         height = 12)


########################################################################################################################
# Supplemental Figure S4
########################################################################################################################
#Figure S4. fastSTRUCTURE for the core collection n=423.
# Pepper K2 Plot with Wild Status
##############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MetBrewer)  # For color palettes

# Load the MetBrewer palette
hiroshige_palette <- met.brewer("Hiroshige", 8)  
k_colors <- hiroshige_palette[c(1, 6)]  # Colors for K1 and K2

# Load the Pepper data
pepper_k2 <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_k2.csv")

# Sort by K1 and K2 percentages
pepper_k2 <- pepper_k2 %>%
  arrange(desc(K1), desc(K2))

# Convert to long format after sorting
long_data <- pivot_longer(pepper_k2, cols = c("K1", "K2"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(pepper_k2$Name))

# Assign colors to each Species (Fixed primary colors)
species_colors <- c(
  "Capsicum annuum" = "blue",
  "Capsicum baccatum" = "yellow",
  "Capsicum chacoense" = "brown",
  "Capsicum chinense" = "red",
  "Capsicum frutescens" = "green"
)

# Assign colors for Wild Status
wild_colors <- c(
  "wild" = "orange",
  "_" = "black"
)

# Generate the plot with K1/K2 stacked bars, species bar, and wild status bar
final_plot <- ggplot() +
  # Plot K1 and K2 stacked bars
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add a colored bar underneath for species
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k2, aes(x = Name, y = -0.05, fill = species), height = 0.1) +
  scale_fill_manual(name = "Species", values = species_colors) +  # Use manually matched colors
  
  # Add a third colored bar for Wild Status
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k2, aes(x = Name, y = -0.10, fill = wild_status), height = 0.1) +
  scale_fill_manual(name = "Wild Status", values = wild_colors) +  # Wild = orange, Not wild = black
  
  labs(title = "K=2", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, 'cm'),
    legend.spacing.y = unit(0.2, 'cm')  # Fixed unit error
  )

# Print the plot
print(final_plot)

ggsave(filename = "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_fastSTRUCTURE_K2_with_wild.pdf",
       plot = final_plot, 
       width = 8,    
       height = 3)
##############
# Pepper K3 Plot with Wild Status
##############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MetBrewer)  # For color palettes

# Load the MetBrewer palette
hiroshige_palette <- met.brewer("Hiroshige", 8)  
k_colors <- hiroshige_palette[c(1, 4, 6)]  # Colors for K1, K2, and K3

# Load the Pepper data (K=3 file)
pepper_k3 <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_k3.csv")

# Sort by K1, K2, and K3 percentages
pepper_k3 <- pepper_k3 %>%
  arrange(desc(K1), desc(K2), desc(K3))

# Convert to long format after sorting
long_data <- pivot_longer(pepper_k3, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(pepper_k3$Name))

# Assign colors to each Species (Fixed primary colors)
species_colors <- c(
  "Capsicum annuum" = "blue",
  "Capsicum baccatum" = "yellow",
  "Capsicum chacoense" = "brown",
  "Capsicum chinense" = "red",
  "Capsicum frutescens" = "green"
)

# Assign colors for Wild Status
wild_colors <- c(
  "wild" = "orange",
  "_" = "black"
)

# Generate the plot with K1/K2/K3 stacked bars, species bar, and wild status bar
final_plot <- ggplot() +
  # Plot K1, K2, and K3 stacked bars
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1, K2, and K3
  
  # Add a colored bar underneath for species
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k3, aes(x = Name, y = -0.05, fill = species), height = 0.1) +
  scale_fill_manual(name = "Species", values = species_colors) +  # Use manually matched colors
  
  # Add a third colored bar for Wild Status
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k3, aes(x = Name, y = -0.10, fill = wild_status), height = 0.1) +
  scale_fill_manual(name = "Wild Status", values = wild_colors) +  # Wild = orange, Not wild = black
  
  labs(title = "K=3", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, 'cm'),
    legend.spacing.y = unit(0.2, 'cm')  # Fixed unit error
  )

# Print the plot
print(final_plot)

# Save the updated plot
ggsave(filename = "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_fastSTRUCTURE_K3_with_wild.pdf",
       plot = final_plot, 
       width = 8,    
       height = 3)

##############
# Pepper K4 Plot with Wild Status
##############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MetBrewer)  # For color palettes

# Load the MetBrewer palette
hiroshige_palette <- met.brewer("Hiroshige", 8)  
k_colors <- hiroshige_palette[c(1, 3, 5, 7)]  # Colors for K1, K2, K3, and K4

# Load the Pepper data (K=4 file)
pepper_k4 <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_k4.csv")

# Sort by K1, K2, K3, and K4 percentages
pepper_k4 <- pepper_k4 %>%
  arrange(desc(K1), desc(K2), desc(K3), desc(K4))

# Convert to long format after sorting
long_data <- pivot_longer(pepper_k4, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(pepper_k4$Name))

# Assign colors to each Species (Fixed primary colors)
species_colors <- c(
  "Capsicum annuum" = "blue",
  "Capsicum baccatum" = "yellow",
  "Capsicum chacoense" = "brown",
  "Capsicum chinense" = "red",
  "Capsicum frutescens" = "green"
)

# Assign colors for Wild Status
wild_colors <- c(
  "wild" = "orange",
  "_" = "black"
)

# Generate the plot with K1/K2/K3/K4 stacked bars, species bar, and wild status bar
final_plot <- ggplot() +
  # Plot K1, K2, K3, and K4 stacked bars
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1, K2, K3, and K4
  
  # Add a colored bar underneath for species
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k4, aes(x = Name, y = -0.05, fill = species), height = 0.1) +
  scale_fill_manual(name = "Species", values = species_colors) +  # Use manually matched colors
  
  # Add a third colored bar for Wild Status
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k4, aes(x = Name, y = -0.10, fill = wild_status), height = 0.1) +
  scale_fill_manual(name = "Wild Status", values = wild_colors) +  # Wild = orange, Not wild = black
  
  labs(title = "K=4", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, 'cm'),
    legend.spacing.y = unit(0.2, 'cm')  # Fixed unit error
  )

# Print the plot
print(final_plot)

# Save the updated plot
ggsave(filename = "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_fastSTRUCTURE_K4_with_wild.pdf",
       plot = final_plot, 
       width = 8,    
       height = 3)


##############
# Pepper K5 Plot with Wild Status
##############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MetBrewer)  # For color palettes

# Load the MetBrewer palette
hiroshige_palette <- met.brewer("Hiroshige", 8)  
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 8)]  # Colors for K1, K2, K3, K4, and K5

# Load the Pepper data (K=5 file)
pepper_k5 <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_k5.csv")

# Sort by K1, K2, K3, K4, and K5 percentages
pepper_k5 <- pepper_k5 %>%
  arrange(desc(K1), desc(K2), desc(K3), desc(K4), desc(K5))

# Convert to long format after sorting
long_data <- pivot_longer(pepper_k5, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(pepper_k5$Name))

# Assign colors to each Species (Fixed primary colors)
species_colors <- c(
  "Capsicum annuum" = "blue",
  "Capsicum baccatum" = "yellow",
  "Capsicum chacoense" = "brown",
  "Capsicum chinense" = "red",
  "Capsicum frutescens" = "green"
)

# Assign colors for Wild Status
wild_colors <- c(
  "wild" = "orange",
  "_" = "black"
)

# Generate the plot with K1/K2/K3/K4/K5 stacked bars, species bar, and wild status bar
final_plot <- ggplot() +
  # Plot K1, K2, K3, K4, and K5 stacked bars
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1, K2, K3, K4, and K5
  
  # Add a colored bar underneath for species
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k5, aes(x = Name, y = -0.05, fill = species), height = 0.1) +
  scale_fill_manual(name = "Species", values = species_colors) +  # Use manually matched colors
  
  # Add a third colored bar for Wild Status
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k5, aes(x = Name, y = -0.10, fill = wild_status), height = 0.1) +
  scale_fill_manual(name = "Wild Status", values = wild_colors) +  # Wild = orange, Not wild = black
  
  labs(title = "K=5", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, 'cm'),
    legend.spacing.y = unit(0.2, 'cm')  # Fixed unit error
  )

# Print the plot
print(final_plot)

# Save the updated plot
ggsave(filename = "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_fastSTRUCTURE_K5_with_wild.pdf",
       plot = final_plot, 
       width = 8,    
       height = 3)
##############
# Pepper K6 Plot with Wild Status
##############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MetBrewer)  # For color palettes

# Load the MetBrewer palette
hiroshige_palette <- met.brewer("Hiroshige", 8)  
k_colors <- hiroshige_palette[c(1, 3, 4, 5, 6, 8)]  # Colors for K1, K2, K3, K4, K5, and K6

# Load the Pepper data (K=6 file)
pepper_k6 <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_k6.csv")

# Sort by K1, K2, K3, K4, K5, and K6 percentages
pepper_k6 <- pepper_k6 %>%
  arrange(desc(K1), desc(K2), desc(K3), desc(K4), desc(K5), desc(K6))

# Convert to long format after sorting
long_data <- pivot_longer(pepper_k6, cols = c("K1", "K2", "K3", "K4", "K5", "K6"), 
                          names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(pepper_k6$Name))

# Assign colors to each Species (Fixed primary colors)
species_colors <- c(
  "Capsicum annuum" = "blue",
  "Capsicum baccatum" = "yellow",
  "Capsicum chacoense" = "brown",
  "Capsicum chinense" = "red",
  "Capsicum frutescens" = "green"
)

# Assign colors for Wild Status
wild_colors <- c(
  "wild" = "orange",
  "_" = "black"
)

# Generate the plot with K1/K2/K3/K4/K5/K6 stacked bars, species bar, and wild status bar
final_plot <- ggplot() +
  # Plot K1, K2, K3, K4, K5, and K6 stacked bars
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 to K6
  
  # Add a colored bar underneath for species
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k6, aes(x = Name, y = -0.05, fill = species), height = 0.1) +
  scale_fill_manual(name = "Species", values = species_colors) +  # Use manually matched colors
  
  # Add a third colored bar for Wild Status
  new_scale_fill() +  # Separate fill scale
  geom_tile(data = pepper_k6, aes(x = Name, y = -0.10, fill = wild_status), height = 0.1) +
  scale_fill_manual(name = "Wild Status", values = wild_colors) +  # Wild = orange, Not wild = black
  
  labs(title = "K=6", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, 'cm'),
    legend.spacing.y = unit(0.2, 'cm')  # Fixed unit error
  )

# Print the plot
print(final_plot)

# Save the updated plot
ggsave(filename = "~/R/World_Veg_Collab_Pepper/filtered_data/metadata/pepper_fastSTRUCTURE_K6_with_wild.pdf",
       plot = final_plot, 
       width = 8,    
       height = 3)


########################################################################################################################
# Supplemental Figure S5
########################################################################################################################
#PCA for the 10k 

#SNPRelate
#########
library(SNPRelate)
#vcf.fn <- "~/R/World_Veg_Collab_Pepper/filtered_data/pepper_10000.vcf.gz"
vcf.fn <- "~/R/World_Veg_Collab_Pepper/filtered_data/pepper_10000.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "pepper_10000_a.gds", method="copy.num.of.ref")
snpgdsSummary("pepper_10000_a.gds")
genofile <- snpgdsOpen("pepper_10000_a.gds")
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
#write.csv(tab, "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")

library(ggplot2)
pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")

# Basic PCA plot: EV1 vs. EV2
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2)) +
  geom_point(size = 3, alpha = 0.8) + # Points for samples
  theme_minimal() +                   # Clean theme
  labs(
    title = "Pepper PCA 10,250",
    x = "EV1 (18.02%)",
    y = "EV2 (13.77%)"
  )

# Display the plot
print(pca_plot)

####
#pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")
pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Pepper_PC_10250.csv")
# Remove trailing zero from PCA sample IDs
pca$sample.id <- sub("0$", "", pca$sample.id)

# Metadata 
#meta <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/10038_Tripodi_2021.csv")
meta <- read.csv("~/R/World_Veg_Project/filtered_data_2/10038_Tripodi_2021.csv")
colnames(meta)
colnames(meta)[colnames(meta) == "Sample_ID_G2P.Sol_database"] <- "sample.id"
merged <- merge(pca, meta, by = "sample.id")
#write.csv(merged, "~/R/World_Veg_Project/filtered_data_2/Pepper_PCA_with_Metadata.csv", row.names = FALSE)

library(ggplot2)

ggplot(merged, aes(x = EV1, y = EV2, color = Organism)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "",
    x = "EV1 (18.02%)",
    y = "EV2 (13.77%)"
  ) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))


########################################################################################################################
# Supplemental Figure S6
########################################################################################################################

#kinship for core
#SNPRelate
#########
library(SNPRelate)
#vcf.fn <- "~/R/World_Veg_Collab_Pepper/filtered_data/pepper_10000.vcf.gz"
vcf.fn <- "~/R/World_Veg_Project/filtered_data_2/pepper_10000.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "pepper_10000_a.gds", method="copy.num.of.ref")
snpgdsSummary("pepper_10000_a.gds")
genofile <- snpgdsOpen("pepper_10000_a.gds")
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
#write.csv(tab, "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")

library(ggplot2)
pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")

# Basic PCA plot: EV1 vs. EV2
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2)) +
  geom_point(size = 3, alpha = 0.8) + # Points for samples
  theme_minimal() +                   # Clean theme
  labs(
    title = "Pepper PCA 10,250",
    x = "EV1 (18.02 %)",
    y = "EV2 (13.77%)"
  )

# Display the plot
print(pca_plot)

####
pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")
pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")
# Remove trailing zero from PCA sample IDs
pca$sample.id <- sub("0$", "", pca$sample.id)

# Metadata 
#meta <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/10038_Tripodi_2021.csv")
meta <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/10038_Tripodi_2021.csv")
colnames(meta)
colnames(meta)[colnames(meta) == "Sample_ID_G2P.Sol_database"] <- "sample.id"
merged <- merge(pca, meta, by = "sample.id")

library(ggplot2)

ggplot(merged, aes(x = EV1, y = EV2, color = Organism)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA of Pepper Core Collection (n=10,124)",
    x = "PC1",
    y = "PC2"
  ) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

#########################
#second 
#########################
# 10K kinship matrix
setwd("~/R/World_Veg_Project/anna_kinship_matrix/")

library(SNPRelate)
library(pheatmap)

# 2. Define File Names 
vcf.fn <- "~/R/World_Veg_Project/anna_kinship_matrix/pepper_10000.vcf.gz"
gds.fn <- "~/R/World_Veg_Project/anna_kinship_matrix/pepper_10000.gds"

# 3. Convert VCF to GDS 
snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")

# 4. Open GDS File 
genofile <- snpgdsOpen(gds.fn)

# 5. Calculate Kinship Matrix 
ibs <- snpgdsIBS(genofile, num.thread = 2)
kinship_matrix <- 1 - ibs$ibs
rownames(kinship_matrix) <- colnames(kinship_matrix) <- ibs$sample.id

# 6. Read Metadata and Prepare Annotation 
# metadata already read in
metadata <- read.csv("Pepper_PCA_species_merged_clean.csv",
                     stringsAsFactors = FALSE,
                     check.names = FALSE)

# Clean IDs
metadata$sample.id <- trimws(as.character(metadata$sample.id))
vcf_ids <- trimws(as.character(ibs$sample.id))

# Keep one row per sample.id (choose the first occurrence)
metadata_uniq <- metadata[!is.na(metadata$sample.id) & metadata$sample.id != "", ]
metadata_uniq <- metadata_uniq[!duplicated(metadata_uniq$sample.id), ]

# Match to VCF order
idx <- match(vcf_ids, metadata_uniq$sample.id)

cat("VCF samples:", length(vcf_ids), "\n")
cat("Matched:", sum(!is.na(idx)), "\n")
cat("Unmatched:", sum(is.na(idx)), "\n")


# Build annotation using Organism
meta2 <- metadata_uniq[idx, ]  # aligned to vcf_ids; rows may be NA where unmatched

annotation <- data.frame(Population = meta2$Organism)
rownames(annotation) <- vcf_ids
annotation$Population[is.na(annotation$Population)] <- "UNKNOWN"

# Colors
pops <- unique(annotation$Population)
ann_colors <- list(Population = setNames(rainbow(length(pops)), pops))



# 7. Plot Heatmap with Colored Annotations 
out_pdf <- "kinship_heatmap_by_population_10k.pdf"

pheatmap(kinship_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         filename = out_pdf,
         width = 12,
         height = 12)


########################################################################################################################
# Supplemental Figure S7
########################################################################################################################
# Figure S7.  fastSTRUCTURE for the global collection (n=10,250) using X SNPs and coloured by species. 


#10k faststructure
library(readxl)   # to read .xlsx files
dplyr::glimpse    # if dplyr not loaded yet, run: library(dplyr)
library(dplyr)

# PCA file
pca <- read.csv("~/R/World_Veg_Project/filtered_data_2/Pepper_PC_10250.csv")

# Species metadata
meta <- read_excel("~/R/World_Veg_Project/filtered_data_2/10417_metadata_species.xlsx")

head(pca)
head(meta)
colnames(pca)
colnames(meta)

#need to add 0 to end of meta IDs to match correctly
meta <- meta %>%
  mutate(Sample_ID_fixed = paste0(Sample_ID, "0"))


meta_unique <- meta %>%
  distinct(Sample_ID_fixed, .keep_all = TRUE)

merged <- pca %>%
  left_join(meta_unique %>% select(Sample_ID_fixed, Organism),
            by = c("sample.id" = "Sample_ID_fixed"))
nrow(merged)


write.csv(
  merged,
  "~/R/World_Veg_Project/filtered_data_2/Pepper_PCA_species_merged_clean.csv",
  row.names = FALSE
)


###############
###############################################
# Load libraries
###############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MetBrewer)
library(readxl)

###############################################
# K2
###############################################
pepper_k2 <- read_excel(
  "~/R/World_Veg_Project/fastSTRUCTURE_10k/10k_pepper_k2.xlsx",
  sheet = 1
)

###############################################
# Clean species names
###############################################
# Replace underscores with spaces
pepper_k2$Organism <- gsub("_", " ", pepper_k2$Organism)

# Replace NA with "Unknown"
pepper_k2$Organism[is.na(pepper_k2$Organism)] <- "Unknown"

###############################################
# Sort samples by ancestry (optional but better visually)
###############################################
pepper_k2 <- pepper_k2 %>%
  arrange(desc(K1), desc(K2))

###############################################
# Reshape to long format for ggplot
###############################################
long_data <- pivot_longer(
  pepper_k2,
  cols = c("K1", "K2"),
  names_to = "category",
  values_to = "value"
)

###############################################
# Reorder individuals for plotting
###############################################
long_data$sample.id <- factor(long_data$sample.id,
                              levels = unique(pepper_k2$sample.id))

###############################################
# K1/K2 colors (MetBrewer)
###############################################
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 6)]

###############################################
# Species colors:
# 1) Keep explicit colors for key species
# 2) Auto-generate colors for new species
# 3) Add Unknown as grey
###############################################

# Explicit colors (your fixed set)
explicit_species_colors <- c(
  "Capsicum annuum"     = "blue",
  "Capsicum baccatum"   = "yellow",
  "Capsicum chacoense"  = "brown",
  "Capsicum chinense"   = "red",
  "Capsicum frutescens" = "green"
)

# All species present in your data
all_species <- sort(unique(pepper_k2$Organism))

# Species not in your explicit list
extra_species <- setdiff(all_species, names(explicit_species_colors))

# Auto colors for extra species
extra_colors <- colorRampPalette(met.brewer("Hokusai3"))(length(extra_species))
names(extra_colors) <- extra_species

# Combine explicit + generated
species_colors <- c(explicit_species_colors, extra_colors)

# Add Unknown explicitly
if ("Unknown" %in% all_species) {
  species_colors["Unknown"] <- "grey30"
}

###############################################
# Build the STRUCTURE-style plot
###############################################
final_plot <- ggplot() +
  
  # ----- K1/K2 ancestry bars -----
geom_bar(
  data = long_data,
  aes(x = sample.id, y = value, fill = category),
  stat = "identity"
) +
  scale_fill_manual(name = "K Categories", values = k_colors) +
  
  # ----- Species tile bar below -----
new_scale_fill() +
  geom_tile(
    data = pepper_k2,
    aes(x = sample.id, y = -0.05, fill = Organism),
    height = 0.10
  ) +
  scale_fill_manual(name = "Species", values = species_colors) +
  
  # ----- Labels and theme -----
labs(
  title = "K = 2",
  x = "",
  y = "Ancestry Proportion"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 1
    ),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )

###############################################
# Print the plot
###############################################
print(final_plot)

###############################################
# Save PDF
###############################################
ggsave(
  filename = "~/R/World_Veg_Project/fastSTRUCTURE_10k/pepper_fastSTRUCTURE_K2_speciesOnly.pdf",
  plot = final_plot,
  width = 12,
  height = 4
)


###############################################
# K3
###############################################
pepper_k3 <- read_excel(
  "~/R/World_Veg_Project/fastSTRUCTURE_10k/10k_pepper_k3.xlsx",
  sheet = 1
)

###############################################
# Clean species names
###############################################
pepper_k3$Organism <- gsub("_", " ", pepper_k3$Organism)
pepper_k3$Organism[is.na(pepper_k3$Organism)] <- "Unknown"

###############################################
# Sort samples by K1 → K2 → K3
###############################################
pepper_k3 <- pepper_k3 %>%
  arrange(desc(K1), desc(K2), desc(K3))

###############################################
# Convert to long format
###############################################
long_data <- pivot_longer(
  pepper_k3,
  cols = c("K1", "K2", "K3"),
  names_to = "category",
  values_to = "value"
)

###############################################
# Reorder samples
###############################################
long_data$sample.id <- factor(long_data$sample.id,
                              levels = unique(pepper_k3$sample.id))

###############################################
# Colors for K1/K2/K3 (3 colors)
###############################################
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 3, 6)]   # choose 3 nice spaced colors

###############################################
# Species colors (explicit + auto)
###############################################
# Explicit fixed colors
explicit_species_colors <- c(
  "Capsicum annuum"     = "blue",
  "Capsicum baccatum"   = "yellow",
  "Capsicum chacoense"  = "brown",
  "Capsicum chinense"   = "red",
  "Capsicum frutescens" = "green"
)

# All species present
all_species <- sort(unique(pepper_k3$Organism))

# Species not in explicit list
extra_species <- setdiff(all_species, names(explicit_species_colors))

# Auto colors for extra species
extra_colors <- colorRampPalette(met.brewer("Hokusai3"))(length(extra_species))
names(extra_colors) <- extra_species

# Combine all species colors
species_colors <- c(explicit_species_colors, extra_colors)

# Add Unknown if present
if ("Unknown" %in% all_species) {
  species_colors["Unknown"] <- "grey30"
}

###############################################
# Build the STRUCTURE-style plot
###############################################
final_plot <- ggplot() +
  
  # ---- K1/K2/K3 stacked bars ----
geom_bar(
  data = long_data,
  aes(x = sample.id, y = value, fill = category),
  stat = "identity"
) +
  scale_fill_manual(name = "K Categories", values = k_colors) +
  
  # ---- Species tile bar ----
new_scale_fill() +
  geom_tile(
    data = pepper_k3,
    aes(x = sample.id, y = -0.05, fill = Organism),
    height = 0.10
  ) +
  scale_fill_manual(name = "Species", values = species_colors) +
  
  # ---- Labels and themes ----
labs(
  title = "K = 3",
  x = "",
  y = "Ancestry Proportion"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 1
    ),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )

###############################################
# Print the plot
###############################################
print(final_plot)

###############################################
# Save PDF
###############################################
ggsave(
  filename = "~/R/World_Veg_Project/fastSTRUCTURE_10k/pepper_fastSTRUCTURE_K3_speciesOnly.pdf",
  plot = final_plot,
  width = 12,
  height = 4
)


###############################################
# K4
###############################################
pepper_k4 <- read_excel(
  "~/R/World_Veg_Project/fastSTRUCTURE_10k/10k_pepper_k4.xlsx",
  sheet = 1
)

###############################################
# Clean species names
###############################################
# Replace underscores with spaces
pepper_k4$Organism <- gsub("_", " ", pepper_k4$Organism)

# Replace NA with "Unknown"
pepper_k4$Organism[is.na(pepper_k4$Organism)] <- "Unknown"

###############################################
# Sort samples by K1 → K2 → K3 → K4
###############################################
pepper_k4 <- pepper_k4 %>%
  arrange(desc(K1), desc(K2), desc(K3), desc(K4))

###############################################
# Convert to long format
###############################################
long_data <- pivot_longer(
  pepper_k4,
  cols = c("K1", "K2", "K3", "K4"),
  names_to = "category",
  values_to = "value"
)

###############################################
# Reorder samples
###############################################
long_data$sample.id <- factor(
  long_data$sample.id,
  levels = unique(pepper_k4$sample.id)
)

###############################################
# Colors for K1–K4 clusters (4 colors)
###############################################
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 3, 5, 7)]   # 4 well-spaced colors

###############################################
# Species colors (explicit + auto)
###############################################

# Explicit fixed colors for main species
explicit_species_colors <- c(
  "Capsicum annuum"     = "blue",
  "Capsicum baccatum"   = "yellow",
  "Capsicum chacoense"  = "brown",
  "Capsicum chinense"   = "red",
  "Capsicum frutescens" = "green"
)

# All species in your dataset
all_species <- sort(unique(pepper_k4$Organism))

# Species not in your explicit list
extra_species <- setdiff(all_species, names(explicit_species_colors))

# Automatically generate colors for these extra species
extra_colors <- colorRampPalette(met.brewer("Hokusai3"))(length(extra_species))
names(extra_colors) <- extra_species

# Combine explicit + extra species colors
species_colors <- c(explicit_species_colors, extra_colors)

# Add Unknown if present
if ("Unknown" %in% all_species) {
  species_colors["Unknown"] <- "grey30"
}

###############################################
# Build the STRUCTURE-style plot
###############################################
final_plot <- ggplot() +
  
  # ---- K1–K4 stacked ancestry bars ----
geom_bar(
  data = long_data,
  aes(x = sample.id, y = value, fill = category),
  stat = "identity"
) +
  scale_fill_manual(name = "K Categories", values = k_colors) +
  
  # ---- Species tile bar (below STRUCTURE plot) ----
new_scale_fill() +
  geom_tile(
    data = pepper_k4,
    aes(x = sample.id, y = -0.05, fill = Organism),
    height = 0.10
  ) +
  scale_fill_manual(name = "Species", values = species_colors) +
  
  # ---- Labels and theme ----
labs(
  title = "K = 4",
  x = "",
  y = "Ancestry Proportion"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 1
    ),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )

###############################################
# Print the plot
###############################################
print(final_plot)

###############################################
# Save PDF
###############################################
ggsave(
  filename = "~/R/World_Veg_Project/fastSTRUCTURE_10k/pepper_fastSTRUCTURE_K4_speciesOnly.pdf",
  plot = final_plot,
  width = 12,
  height = 4
)

###############################################
# K5
###############################################
pepper_k5 <- read_excel(
  "~/R/World_Veg_Project/fastSTRUCTURE_10k/10k_pepper_k5.xlsx",
  sheet = 1
)

###############################################
# Clean species names
###############################################
pepper_k5$Organism <- gsub("_", " ", pepper_k5$Organism)
pepper_k5$Organism[is.na(pepper_k5$Organism)] <- "Unknown"

###############################################
# Sort samples by K1 → K2 → K3 → K4 → K5
###############################################
pepper_k5 <- pepper_k5 %>%
  arrange(desc(K1), desc(K2), desc(K3), desc(K4), desc(K5))

###############################################
# Convert to long format
###############################################
long_data <- pivot_longer(
  pepper_k5,
  cols = c("K1", "K2", "K3", "K4", "K5"),
  names_to = "category",
  values_to = "value"
)

###############################################
# Reorder samples
###############################################
long_data$sample.id <- factor(
  long_data$sample.id,
  levels = unique(pepper_k5$sample.id)
)

###############################################
# Colors for K1–K5 clusters (5 colors)
###############################################
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 3, 4, 5, 6)]  # 5 distinct colors

###############################################
# Species colors (explicit + auto)
###############################################

# Explicit fixed colors for main species
explicit_species_colors <- c(
  "Capsicum annuum"     = "blue",
  "Capsicum baccatum"   = "yellow",
  "Capsicum chacoense"  = "brown",
  "Capsicum chinense"   = "red",
  "Capsicum frutescens" = "green"
)

# All species in the dataset
all_species <- sort(unique(pepper_k5$Organism))

# Species not in explicit list
extra_species <- setdiff(all_species, names(explicit_species_colors))

# Auto-colors for extra species
extra_colors <- colorRampPalette(met.brewer("Hokusai3"))(length(extra_species))
names(extra_colors) <- extra_species

# Combine explicit + auto-generated colors
species_colors <- c(explicit_species_colors, extra_colors)

# Add Unknown (if present)
if ("Unknown" %in% all_species) {
  species_colors["Unknown"] <- "grey30"
}

###############################################
# Build the STRUCTURE-style plot
###############################################
final_plot <- ggplot() +
  
  # ---- K1–K5 stacked bars ----
geom_bar(
  data = long_data,
  aes(x = sample.id, y = value, fill = category),
  stat = "identity"
) +
  scale_fill_manual(name = "K Categories", values = k_colors) +
  
  # ---- Species tile bar ----
new_scale_fill() +
  geom_tile(
    data = pepper_k5,
    aes(x = sample.id, y = -0.05, fill = Organism),
    height = 0.10
  ) +
  scale_fill_manual(name = "Species", values = species_colors) +
  
  # ---- Labels and theme ----
labs(
  title = "K = 5",
  x = "",
  y = "Ancestry Proportion"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 1
    ),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )

###############################################
# Print the plot
###############################################
print(final_plot)

###############################################
# Save PDF
###############################################
ggsave(
  filename = "~/R/World_Veg_Project/fastSTRUCTURE_10k/pepper_fastSTRUCTURE_K5_speciesOnly.pdf",
  plot = final_plot,
  width = 12,
  height = 4
)


###############################################
# K6
###############################################
pepper_k6 <- read_excel(
  "~/R/World_Veg_Project/fastSTRUCTURE_10k/10k_pepper_k6.xlsx",
  sheet = 1
)

###############################################
# Clean species names
###############################################
pepper_k6$Organism <- gsub("_", " ", pepper_k6$Organism)
pepper_k6$Organism[is.na(pepper_k6$Organism)] <- "Unknown"

###############################################
# Sort samples by K1 → K2 → K3 → K4 → K5 → K6
###############################################
pepper_k6 <- pepper_k6 %>%
  arrange(desc(K1), desc(K2), desc(K3), desc(K4), desc(K5), desc(K6))

###############################################
# Convert to long format
###############################################
long_data <- pivot_longer(
  pepper_k6,
  cols = c("K1", "K2", "K3", "K4", "K5", "K6"),
  names_to = "category",
  values_to = "value"
)

###############################################
# Reorder samples
###############################################
long_data$sample.id <- factor(
  long_data$sample.id,
  levels = unique(pepper_k6$sample.id)
)

###############################################
# Colors for K1–K6 clusters (6 colors)
###############################################
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 2, 3, 5, 6, 7)]  # 6 distinct colors

###############################################
# Species colors (explicit + auto)
###############################################

# Explicit fixed colors for main species
explicit_species_colors <- c(
  "Capsicum annuum"     = "blue",
  "Capsicum baccatum"   = "yellow",
  "Capsicum chacoense"  = "brown",
  "Capsicum chinense"   = "red",
  "Capsicum frutescens" = "green"
)

# All species in the dataset
all_species <- sort(unique(pepper_k6$Organism))

# Species not in your explicit list
extra_species <- setdiff(all_species, names(explicit_species_colors))

# Auto-colors for extra species
extra_colors <- colorRampPalette(met.brewer("Hokusai3"))(length(extra_species))
names(extra_colors) <- extra_species

# Combine explicit + auto-generated colors
species_colors <- c(explicit_species_colors, extra_colors)

# Add Unknown (if present)
if ("Unknown" %in% all_species) {
  species_colors["Unknown"] <- "grey30"
}

###############################################
# Build the STRUCTURE-style plot
###############################################
final_plot <- ggplot() +
  
  # ---- K1–K6 stacked bars ----
geom_bar(
  data = long_data,
  aes(x = sample.id, y = value, fill = category),
  stat = "identity"
) +
  scale_fill_manual(name = "K Categories", values = k_colors) +
  
  # ---- Species tile bar ----
new_scale_fill() +
  geom_tile(
    data = pepper_k6,
    aes(x = sample.id, y = -0.05, fill = Organism),
    height = 0.10
  ) +
  scale_fill_manual(name = "Species", values = species_colors) +
  
  # ---- Labels and theme ----
labs(
  title = "K = 6",
  x = "",
  y = "Ancestry Proportion"
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 1
    ),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )

###############################################
# Print the plot
###############################################
print(final_plot)

###############################################
# Save PDF
###############################################
ggsave(
  filename = "~/R/World_Veg_Project/fastSTRUCTURE_10k/pepper_fastSTRUCTURE_K6_speciesOnly.pdf",
  plot = final_plot,
  width = 12,
  height = 4
)


########################################################################################################################
# Supplemental Figure S8
########################################################################################################################
#Figure S8. 73 trait prediction accuracies for control timepoint for the core collection.


#########################################################################################################
#GS for Heat 1/2/ctrl - on n=423 dataset
#using n=150/286 as training from the core 423
#########################################################################################################
######################Outline for document#########################
# Step 1: 
#1A - Pre-processing - conversion of .vcf file to genotype matrix
#1B - Convert genotype matrix to rrBLUP format
#1C - Transpose genotype matrix for downstream use
#1D - Genotype import & cleaning/imputation 
####################################################################
# Step 2:
# Genotype import
# Phenotype import
# Set training and test populations
# Run GS
####################################################################
# Step 3:
# Prediction accuracy using
# rrBLUP
# Gaussian Kernel
# Exponential Kernel
####################################################################

####################################################################
#1A FROM VCF TO GENOTYPE MATRIX - rows=individuals  columns =SNPs, in 0/1/2 format
####################################################################
library(vcfR)
library(adegenet)
vcf <- read.vcfR("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_filtered_and_LD_pruned.vcf.gz")
genotype_matrix <- extract.gt(vcf, element = "GT")

dim(genotype_matrix)  # Should return (SNPs × Samples)
write.csv(genotype_matrix, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_raw.csv", row.names = TRUE)

genotype_matrix_raw <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_raw.csv", 
                                row.names = 1, check.names = FALSE)

# Check dimensions
dim(genotype_matrix_raw)  # Should match (SNPs × Samples)

# Preview 
head(genotype_matrix_raw[, 1:5])  #first 5 samples

####################################################################
#1B - Convert genotype matrix to rrBLUP format
####################################################################
# Function to convert genotypes to numeric format for rrBLUP
convert_geno <- function(x) {
  ifelse(x == "0/0", 0, 
         ifelse(x == "0/1" | x == "1/0", 1, 
                ifelse(x == "1/1", 2, NA)))
}

# Apply conversion to the entire matrix
geno_numeric <- apply(genotype_matrix_raw, 2, convert_geno)

# Convert to dataframe
geno_df <- as.data.frame(geno_numeric)

# Check dimensions after conversion
dim(geno_df)  #(340734 × 423)

# Save the converted genotype matrix
write.csv(geno_df, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format.csv", row.names = TRUE)

####################################################################
#1C - Transpose genotype matrix for downstream use
####################################################################
library(rrBLUP)

# genotype matrix in rrBLUP format
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format.csv", 
                    row.names = 1, check.names = FALSE)

# Check dimensions to confirm SNPs x Samples
dim(geno_df)

# Preview first few values
head(geno_df[, 1:5])

geno_df <- t(geno_df)  # Transpose so genotypes (samples) become rows
geno_df <- as.data.frame(geno_df)  # Convert back to a dataframe
head(geno_df[, 1:5])
cat("Corrected genotype matrix dimensions:", dim(geno_df), "\n")


write.csv(geno_df, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_transposed.csv", row.names = TRUE)


##################
#1D - GENOTYPE IMPORT & CLEANING/IMPUTATION
##################

geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_transposed.csv", 
                    row.names = 1, check.names = FALSE)

# Define imputation function (replace missing values with the most common allele)
impute_geno <- function(x) {
  x[is.na(x)] <- as.numeric(names(which.max(table(x, useNA = "no"))))  # Replace NA with mode
  return(x)
}

# Apply the imputation function to each SNP column
cat("Imputing missing values in genotype data...\n")
geno_df_imputed <- apply(geno_df, 2, impute_geno)

# Convert back to dataframe
geno_df_imputed <- as.data.frame(geno_df_imputed)

# Confirm missing values after imputation
num_missing_after <- sum(is.na(geno_df_imputed))
cat("Missing values after imputation:", num_missing_after, "\n") # Should be 0

# Save the cleaned genotype matrix
write.csv(geno_df_imputed, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = TRUE)


############################################################################################################
# Step 2
############################################################################################################
library(rrBLUP)
library(dplyr)

##################
# GENOTYPE IMPORT (IMPUTED)
##################
geno_df_imputed <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv",
                            row.names = 1, check.names = FALSE)

##################
# PHENOTYPE IMPORT
##################

control_pheno_check <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)

##################
#TRAINING
##################
set.seed(123)  # For reproducibility

# Identify lines with complete phenotype data across ALL conditions
complete_lines <- rownames(control_pheno_check)[
  rowSums(is.na(control_pheno_check)) == 0 &
    rowSums(is.na(heat_stress_1_check)) == 0 &
    rowSums(is.na(heat_stress_2_check)) == 0
]

# Randomly sample 150 genotypes only from the fully complete lines
train_genotypes <- sample(complete_lines, 150)

##################
#TEST
##################
# Define test set as remaining lines (all genotyped lines except those in training)
test_genotypes <- setdiff(rownames(geno_df_imputed), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_df_imputed[train_genotypes, ]
geno_test <- geno_df_imputed[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

# Check missing values again
cat("Missing values in control_train:", sum(is.na(control_train)), "\n")
cat("Missing values in heat_stress_1_train:", sum(is.na(heat_stress_1_train)), "\n")
cat("Missing values in heat_stress_2_train:", sum(is.na(heat_stress_2_train)), "\n")



#convert genotype matrix to numeric for rrBLUP compatibility
geno_train <- as.matrix(geno_train)
mode(geno_train) <- "numeric"

geno_test <- as.matrix(geno_test)
mode(geno_test) <- "numeric"

############################################################################################################
library(rrBLUP)

# Function to train rrBLUP and predict GEBVs for multiple traits
run_rrBLUP_all_traits <- function(geno_train, geno_test, pheno_train, output_prefix) {
  
  traits <- colnames(pheno_train)  # Get trait names
  GEBVs_train_list <- list()
  GEBVs_test_list <- list()
  
  for (trait in traits) {
    cat("Training rrBLUP model for", trait, "...\n")
    
    y <- pheno_train[[trait]]  # Extract trait values
    y[is.na(y)] <- mean(y, na.rm = TRUE)  # Impute missing values (just in case)
    
    # Skip if variance is zero (to prevent errors)
    if (var(y, na.rm = TRUE) == 0) {
      cat("Skipping", trait, "due to zero variance.\n")
      next
    }
    
    # Train rrBLUP model
    rrblup_model <- mixed.solve(y = y, Z = geno_train)
    
    # Predict GEBVs for training and test sets
    GEBVs_train <- geno_train %*% rrblup_model$u
    GEBVs_test <- geno_test %*% rrblup_model$u
    
    # Convert to data frames and store
    GEBVs_train_list[[trait]] <- as.data.frame(GEBVs_train)
    GEBVs_test_list[[trait]] <- as.data.frame(GEBVs_test)
    
    # Rename columns
    colnames(GEBVs_train_list[[trait]]) <- paste0("GEBV_", trait, "_Train")
    colnames(GEBVs_test_list[[trait]]) <- paste0("GEBV_", trait, "_Test")
  }
  
  # Combine GEBVs for all traits
  GEBVs_train_final <- do.call(cbind, GEBVs_train_list)
  GEBVs_test_final <- do.call(cbind, GEBVs_test_list)
  
  # Add row names
  rownames(GEBVs_train_final) <- rownames(geno_train)
  rownames(GEBVs_test_final) <- rownames(geno_test)
  
  # Save results
  write.csv(GEBVs_train_final, paste0(output_prefix, "_Train.csv"), row.names = TRUE)
  write.csv(GEBVs_test_final, paste0(output_prefix, "_Test.csv"), row.names = TRUE)
  
  cat("GEBV predictions saved for all traits.\n")
}

#####
#Save data - rrBLUP for all 73 traits - GS results
#####

# Run rrBLUP for all traits in all conditions
run_rrBLUP_all_traits(geno_train, geno_test, control_train, 
                      "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control")
run_rrBLUP_all_traits(geno_train, geno_test, heat_stress_1_train, 
                      "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1")
run_rrBLUP_all_traits(geno_train, geno_test, heat_stress_2_train, 
                      "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2")

#the above ran test and training gebvs separate so then fused for each to complete out 423 manually in excel and named the below
#in output folder
#GEBVs_all_n423_73trait_control.csv
#GEBVs_all_n423_73trait_Heatstress_1.csv
#GEBVs_all_n423_73trait_Heatstress_2.csv


# Save the 150 training genotypes to a CSV file
write.csv(data.frame(Training_Genotypes = train_genotypes), 
          "~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv", 
          row.names = FALSE)

#########################################################################################################
#Step 3
#########################################################################################################
#Model Prediction Accuracies
#########################################################################################################

####################
#Method one: rrBLUP
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
#STEP ONE
####################

# Load Genotype Data (in rrBLUP format)
#geno_data <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Environmental and Phenotypic Data
control_pheno_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/heat_stress_2_filtered_286.csv", row.names = 1)

# Identify common genotypes
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))

# Filter the genotype matrix to match phenotype data
geno_filtered <- geno_data[common_genotypes, ]

# Ensure genotype data is in numeric format
geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
#STEP TWO
####################
# Load the previously saved training genotypes
#train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
train_genotypes <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# Define the test set as the remaining genotypes
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes, ]
geno_test <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]


####################
#STEP THREE
####################

# Load custom k-fold cross-validation function
source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
xval_k10_rrblup <- k.xval(g.in = geno_train, y.in = control_train, y.trainset = control_train, k.fold = 10, reps = 50)
#saveRDS(xval_k10_rrblup, "~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")
saveRDS(xval_k10_rrblup, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")


####################
#Method Two: Gaussian Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
#STEP ONE
####################

# Load Genotype Data (in rrBLUP format)
#geno_data <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Environmental and Phenotypic Data
control_pheno_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/heat_stress_2_filtered_286.csv", row.names = 1)

# Identify common genotypes
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))

# Filter the genotype matrix to match phenotype data
geno_filtered <- geno_data[common_genotypes, ]

# Ensure genotype data is in numeric format
geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
#STEP TWO
####################
# Load the previously saved training genotypes
#train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
train_genotypes <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# Define the test set as the remaining genotypes
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes, ]
geno_test <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

##############################################################################################################################################################
#add genotype column names to control_train as kin.blup needs it
control_train$Genotype <- rownames(control_train)

control_train <- control_train %>%
  relocate(Genotype)

print(colnames(control_train)[1])  # This should be "Genotype"

####################
#STEP THREE
####################

# Load custom k-fold cross-validation function
source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
#xval_k10_rrblup <- k.xval(g.in = geno_train, y.in = control_train, y.trainset = control_train, k.fold = 10, reps = 50)
#saveRDS(xval_k10_rrblup, "~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")
#saveRDS(xval_k10_rrblup, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")

# Compute relationship matrix using Gaussian Kernel
K <- A.mat(geno_train)
k_dist <- dist(K)

# Run k-fold cross-validation for Gaussian Kernel
xval_k10_GAUSS <- k.xval.GAUSS(g.in = geno_train, y.in = control_train, y.trainset = control_train, k_dist = k_dist, k.fold = 10, reps = 50)
#saveRDS(xval_k10_GAUSS, "~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_kfold_10.RData")
saveRDS(xval_k10_GAUSS, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies2/xval_GAUSS_kfold_10.RData")


####################
#Method Three: Exponential Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
#STEP ONE
####################

# Load Genotype Data (in rrBLUP format)
#geno_data <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
geno_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load Environmental and Phenotypic Data
control_pheno_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/heat_stress_2_filtered_286.csv", row.names = 1)

# Identify common genotypes
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))

# Filter the genotype matrix to match phenotype data
geno_filtered <- geno_data[common_genotypes, ]

# Ensure genotype data is in numeric format
geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

####################
#STEP TWO
####################
# Load the previously saved training genotypes
#train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/Training_Genotypes_n150_run1.csv")$Training_Genotypes
train_genotypes <- read.csv("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/Training_Genotypes_n150_run1.csv")$Training_Genotypes

# Define the test set as the remaining genotypes
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes, ]
geno_test <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets for training only (Now guaranteed to be complete)
control_train <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

##############################################################################################################################################################
#add genotype column names to control_train as kin.blup needs it
control_train$Genotype <- rownames(control_train)

control_train <- control_train %>%
  relocate(Genotype)

print(colnames(control_train)[1])  # This should be "Genotype"

####################
#STEP THREE
####################

# Load custom k-fold cross-validation function
source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/xval_kfold_functions.R")

# Run k-fold cross-validation for rrBLUP
#xval_k10_rrblup <- k.xval(g.in = geno_train, y.in = control_train, y.trainset = control_train, k.fold = 10, reps = 50)
#saveRDS(xval_k10_rrblup, "~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")
#saveRDS(xval_k10_rrblup, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")

# Compute relationship matrix using Gaussian Kernel
K <- A.mat(geno_train)
k_dist <- dist(K)

# Run k-fold cross-validation for Gaussian Kernel
#xval_k10_GAUSS <- k.xval.GAUSS(g.in = geno_train, y.in = control_train, y.trainset = control_train, k_dist = k_dist, k.fold = 10, reps = 50)
#saveRDS(xval_k10_GAUSS, "~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_kfold_10.RData")
#saveRDS(xval_k10_GAUSS, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_GAUSS_kfold_10.RData")

# Compute relationship matrix using Exponential Kernel
K.Exp <- Kernel_computation(X = geno_train, name = "exponential", degree = NULL, nL = NULL)
exp_dist <- dist(K.Exp)

##############################################################################################################################################################

rownames(geno_train) <- as.character(trimws(rownames(geno_train)))
control_train$Genotype <- as.character(trimws(control_train$Genotype))

print(setdiff(rownames(geno_train), control_train$Genotype))  # Should return character(0)
print(setdiff(control_train$Genotype, rownames(geno_train)))  # Should return character(0)

exp_dist_mat <- as.matrix(exp_dist)  # Convert from dist object
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)  # Assign correct row/col names

print(setdiff(rownames(exp_dist_mat), rownames(geno_train)))  # Should return character(0)
print(setdiff(rownames(exp_dist_mat), control_train$Genotype))  # Should return character(0)

##############################################################################################################################################################
# Run k-fold cross-validation for Exponential Kernel
xval_k10_EXP <- k.xval.EXP(g.in = geno_train, y.in = control_train, y.trainset = control_train, k_dist = exp_dist_mat, k.fold = 10, reps = 50)
#saveRDS(xval_k10_EXP, "~/R/World_Veg_Collab_Pepper/outputs/xval_EXP_kfold_10.RData")
saveRDS(xval_k10_EXP, "/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies3/xval_EXP_kfold_10.RData")

#########################################################################################################

####################
#Plotting Prediction Accuracies
####################
# Load cross-validation results
rrblup_kfold10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_rrblup_kfold_10.RData")$xval.result
gauss_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_GAUSS_kfold_10.RData")$xval.result
EXP_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_EXP_kfold_10.RData")$xval.result

#rrblup_kfold10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")$xval.result
#gauss_kfold_10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_GAUSS_kfold_10.RData")$xval.result
#EXP_kfold_10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_EXP_kfold_10.RData")$xval.result
#bayescpi_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_BayesCpi_kfold_10.RData")$xval.result

# Assign model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
#bayescpi_kfold_10$model <- "BayesCpi"

# Combine all results
#all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10, bayescpi_kfold_10)
all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

# Convert standard deviation to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

# Select traits to analyze
all_gs_traits <- all_models %>% filter(trait %in% unique(all_models$trait))  # Automatically selects all available traits

####################
#STEP FIVE
####################
#Plot prediction accuracy
ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +  # Adjust the number of rows as needed
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +  # Baseline for interpretation
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits in Chili Dataset",
       x = "Model",
       y = "Prediction Accuracy (r)")

#################################
# Supplemental Figure 8 PLOT
####################
#Plotting Prediction Accuracies
####################
# Load cross-validation results
rrblup_kfold10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_rrblup_kfold_10.RData")$xval.result
gauss_kfold_10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_GAUSS_kfold_10.RData")$xval.result
EXP_kfold_10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_EXP_kfold_10.RData")$xval.result

#rrblup_kfold10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_rrblup_kfold_10.RData")$xval.result
#gauss_kfold_10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_GAUSS_kfold_10.RData")$xval.result
#EXP_kfold_10 <- readRDS("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_EXP_kfold_10.RData")$xval.result
#bayescpi_kfold_10 <- readRDS("~/R/World_Veg_Collab_Pepper/outputs/xval_BayesCpi_kfold_10.RData")$xval.result

# Assign model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
#bayescpi_kfold_10$model <- "BayesCpi"

# Combine all results
#all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10, bayescpi_kfold_10)
all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

# Convert standard deviation to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

# Select traits to analyze
all_gs_traits <- all_models %>% filter(trait %in% unique(all_models$trait))  # Automatically selects all available traits

############################
# Convert columns to numeric (fixing character/factor issues)
all_gs_traits$r.mean <- as.numeric(all_gs_traits$r.mean)
all_gs_traits$r.sd <- as.numeric(all_gs_traits$r.sd)

# Check if there are any NA values after conversion
summary(all_gs_traits$r.mean)
summary(all_gs_traits$r.sd)

# Remove rows where r.mean or r.sd is NA (caused by conversion errors)
all_gs_traits <- all_gs_traits %>% filter(!is.na(r.mean), !is.na(r.sd))

p<- ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits",
       x = "Model",
       y = "Prediction Accuracy (r)")
p

ggsave("~/R/World_Veg_Project/PANNELS/prediction_accuracy_plot_1.pdf", plot = p, width = 24, height = 13, units = "in")


#########################################
# Figure X
#########################################
#follows on from the above 
#Filter traits with mean prediction accuracy above 0.5
library(dplyr)

all_gs_traits_filtered <- all_gs_traits %>% filter(r.mean > 0.5)

rrblup_traits_above_0.5 <- rrblup_traits_above_0.5 %>%
  mutate(facet_label = paste0(trait, " (PA = ", round(r.mean, 2), ")"))

p<- ggplot(rrblup_traits_above_0.5, aes(x = model, y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
                width = 0.3, position = position_dodge(0.2), color = "red") +
  geom_point(size = 3, color = "red") +
  facet_wrap(vars(facet_label), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  ylim(0, 1) +
  labs(
    title = "",
    x = "",
    y = "Prediction Accuracy (r)"
  )

p

ggsave("~/R/World_Veg_Collab_Pepper/FIGS/prediction_accuracy_plot_1a.pdf", plot = p, width = 18, height = 13, units = "in")


########################################################################################################################
# Supplemental Figure S9
########################################################################################################################
#Figure S9. 73 trait prediction accuracies for heat stress #1 timepoint for the core collection.


####################
#Plotting Prediction Accuracies - Heat stress #1
####################
# Load cross-validation results
rrblup_kfold10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_rrblup_kfold_10_heat1_model_1.RData")$xval.result
gauss_kfold_10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_GAUSS_kfold_10_heat1_model2.RData")$xval.result
EXP_kfold_10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_EXP_kfold_10_Heat1_model3.RData")$xval.result

# Assign model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"


# Combine all results
all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

# Convert standard deviation to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

# Select traits to analyze
all_gs_traits <- all_models %>% filter(trait %in% unique(all_models$trait))  # Automatically selects all available traits

############################
# Convert columns to numeric (fixing character/factor issues)
all_gs_traits$r.mean <- as.numeric(all_gs_traits$r.mean)
all_gs_traits$r.sd <- as.numeric(all_gs_traits$r.sd)

# Check if there are any NA values after conversion
summary(all_gs_traits$r.mean)
summary(all_gs_traits$r.sd)

# Remove rows where r.mean or r.sd is NA (caused by conversion errors)
all_gs_traits <- all_gs_traits %>% filter(!is.na(r.mean), !is.na(r.sd))

p <- ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits",
       x = "Model",
       y = "Prediction Accuracy (r)")

p

ggsave("~/R/World_Veg_Project/PANNELS/prediction_accuracy_plot_2.pdf", plot = p, width = 24, height = 13, units = "in")



########################################################################################################################
# Supplemental Figure S10
########################################################################################################################
#Figure S10. 73 trait prediction accuracies for heat stress #2 timepoint for the core collection.

####################
#Plotting Prediction Accuracies - Heat stress #2
####################
# Load cross-validation results
rrblup_kfold10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_rrblup_kfold_10_heat2_model_1.RData")$xval.result
gauss_kfold_10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_GAUSS_kfold_10_heat2_model2.RData")$xval.result
EXP_kfold_10 <- readRDS("~/R/World_Veg_Project/filtered_data_2/xval_EXP_kfold_10_Heat2_model3.RData")$xval.result

# Assign model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"


# Combine all results
all_models <- bind_rows(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

# Convert standard deviation to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

# Select traits to analyze
all_gs_traits <- all_models %>% filter(trait %in% unique(all_models$trait))  # Automatically selects all available traits

############################
# Convert columns to numeric (fixing character/factor issues)
all_gs_traits$r.mean <- as.numeric(all_gs_traits$r.mean)
all_gs_traits$r.sd <- as.numeric(all_gs_traits$r.sd)

# Check if there are any NA values after conversion
summary(all_gs_traits$r.mean)
summary(all_gs_traits$r.sd)

# Remove rows where r.mean or r.sd is NA (caused by conversion errors)
all_gs_traits <- all_gs_traits %>% filter(!is.na(r.mean), !is.na(r.sd))

p <- ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits",
       x = "Model",
       y = "Prediction Accuracy (r)")

p

ggsave("~/R/World_Veg_Project/PANNELS/prediction_accuracy_plot_3.pdf", plot = p, width = 24, height = 13, units = "in")






########################################################################################################################
# Supplemental Figure S11
########################################################################################################################
#Figure S11. Line rank changes for traits (n=13) with prediction accuracies above 0.5 for rrBLUP across all three environmental conditions for the core collection. 

############################################################
# Rank-change pattern plots for 13 traits across 3 environments
# FINAL VERSION — one legend, correct colours, publication-ready
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)

############################################################
# Pattern classification function
############################################################
get_pattern <- function(r1, r2, r3) {
  if (any(is.na(c(r1, r2, r3)))) return(NA)
  if (r1 > r2 & r2 > r3) return("Decreasing")
  if (r1 < r2 & r2 < r3) return("Increasing")
  if (r1 < r2 & r3 < r2) return("U-Shaped")
  if (r1 > r2 & r3 > r2) return("Inverted-U")
  return("Irregular")
}

############################################################
# Load GEBVs
############################################################
control <- read.csv("~/R/World_Veg_Project/data/outputs/GEBVs_all_n423_73trait_control.csv",
                    row.names = 1, check.names = FALSE)
heat1   <- read.csv("~/R/World_Veg_Project/data/outputs/GEBVs_all_n423_73trait_Heatstress_1.csv",
                    row.names = 1, check.names = FALSE)
heat2   <- read.csv("~/R/World_Veg_Project/data/outputs/GEBVs_all_n423_73trait_Heatstress_2.csv",
                    row.names = 1, check.names = FALSE)

############################################################
# Traits with PA > 0.5
############################################################
top_traits <- c(
  "GEBV_biomassfinalplant", "GEBV_biomassfinalplot", "GEBV_fruitlength", "GEBV_fruitno",
  "GEBV_fruitshapeindex", "GEBV_fruitweight", "GEBV_fruitwidth", "GEBV_leafangledaily",
  "GEBV_leafareadaily", "GEBV_leafareadailyplant",
  "GEBV_leafareaprojecteddaily", "GEBV_leafareaprojectedplant", "GEBV_yield"
)

############################################################
# Generate plots (NO legends here)
############################################################
plot_list <- list()

for (trait in top_traits) {
  
  # Long format GEBVs
  gebvs_long <- data.frame(
    Line = rownames(control),
    control = control[[trait]],
    heat1   = heat1[[trait]],
    heat2   = heat2[[trait]]
  ) %>%
    pivot_longer(cols = c("heat1", "control", "heat2"),
                 names_to = "Condition", values_to = "GEBV") %>%
    mutate(
      Condition = factor(Condition, levels = c("heat1", "control", "heat2"))
    ) %>%
    group_by(Condition) %>%
    mutate(Rank = rank(-GEBV)) %>%  # higher GEBV = higher rank
    ungroup()
  
  # Classify pattern per line
  rank_wide <- gebvs_long %>%
    select(Line, Condition, Rank) %>%
    pivot_wider(names_from = Condition, values_from = Rank) %>%
    mutate(Pattern = mapply(get_pattern, heat1, control, heat2))
  
  # Fixed factor levels
  rank_wide$Pattern <- factor(rank_wide$Pattern,
                              levels = c("Decreasing", "Increasing", "U-Shaped", "Inverted-U", "Irregular")
  )
  
  gebvs_long <- left_join(gebvs_long, rank_wide[, c("Line", "Pattern")], by = "Line")
  
  gebvs_long$Pattern <- factor(gebvs_long$Pattern,
                               levels = c("Decreasing", "Increasing", "U-Shaped", "Inverted-U", "Irregular")
  )
  
  # Plot WITHOUT legend
  p <- ggplot(gebvs_long,
              aes(x = Condition, y = Rank, group = Line, color = Pattern)) +
    geom_line(alpha = 0.5) +
    geom_point(size = 0.8) +
    scale_y_reverse() +
    scale_color_brewer(
      palette = "Dark2",
      drop = FALSE,
      limits = c("Decreasing", "Increasing", "U-Shaped", "Inverted-U", "Irregular")
    ) +
    labs(title = trait, x = NULL, y = "Rank", color = "Pattern") +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  plot_list[[trait]] <- p
}

############################################################
# Extract ONE legend from ONE plot
############################################################
get_legend <- function(myplot) {
  g <- ggplot_gtable(ggplot_build(myplot))
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[idx]]
}

# Create a dummy plot with legend ON
legend_plot <- plot_list[[1]] +
  theme(legend.position = "right")

shared_legend <- get_legend(legend_plot)

############################################################
# Combine panels + legend
############################################################
# Combine subplots into multipanel patchwork
panel_grid <- wrap_plots(plot_list, ncol = 5)

# Convert patchwork → grob (CRITICAL FIX)
panel_grid_grob <- patchworkGrob(panel_grid)

# Combine with legend
final_plot <- grid.arrange(
  panel_grid_grob,
  shared_legend,
  ncol = 2,
  widths = c(18, 2)
)


# Save output
ggsave("~/R/World_Veg_Project/PANNELS/rank_change_top13_patterns_FINAL.pdf",
       plot = final_plot, width = 22, height = 18)

############################################################
# Export PDF (vector)
############################################################
ggsave("~/R/World_Veg_Project/PANNELS/rank_change_top13_pattern_colored.pdf",
       plot = final_plot, width = 20, height = 16)



########################################################################################################################
# Supplemental Figure S12
########################################################################################################################
#Figure S12. 73 trait prediction accuracies for control timepoint in the global collection.


#########################################################################################################
#GS for Heat 1/2/ctrl - on n= 10,000  dataset
#using n=259 as training from the core 423
#########################################################################################################
#SNPRelate
#########
library(SNPRelate)
vcf.fn <- "~/R/World_Veg_Collab_Pepper/filtered_data/pepper_10000.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "pepper_10000_a.gds", method="copy.num.of.ref")
snpgdsSummary("pepper_10000_a.gds")
genofile <- snpgdsOpen("pepper_10000_a.gds")
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
write.csv(tab, "~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")


library(ggplot2)
pca <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/Pepper_PC_10250.csv")

# Basic PCA plot: EV1 vs. EV2
pca_plot <- ggplot(pca, aes(x = EV1, y = EV2)) +
  geom_point(size = 3, alpha = 0.8) + # Points for samples
  theme_minimal() +                   # Clean theme
  labs(
    title = "Pepper PCA 10,250",
    x = "EV1 (18.02 %)",
    y = "EV2 (13.77%)"
  )

# Display the plot
print(pca_plot)

################################################################################################
#OVERLAPS BETWEEN 10K+ vcf (10,250 samples) AND 289 PHENOTYPES

# Load the files
df1 <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/Pepper_PC_10250.csv")
#df1 <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/10038_Tripodi_2021.csv")
df2 <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/control_pheno_filtered_286.csv")

# Check column names to identify the sample ID column
colnames(df1)
colnames(df2)

# Assuming the sample ID columns are:
df1$Sample_ID_G2P_Sol_database
df2$Line
# Adjust if the column names are different!

# Find overlapping sample IDs
#overlap <- intersect(df1$Sample_ID_G2P_Sol_database, df2$Line)
overlap <- intersect(df1$sample.id, df2$Line)
# View overlapping samples
print(overlap)
length(overlap)  # how many overlap

# Save the overlapping sample IDs to a CSV file
write.csv(overlap, "~/R/World_Veg_Collab_Pepper/10250_GS/overlapping_sample_ids_282.csv", row.names = FALSE)


#df3 <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/Matched_Pepper_Data_423.csv")
#sample.id

# Assuming both use "sample.id" as the ID column
#overlap_1_3 <- intersect(df1$Sample_ID_G2P_Sol_database, df3$sample.id)
#overlap_1_3 <- intersect(df1$sample.id, df3$sample.id)

# View results
#print(overlap_1_3)
#length(overlap_1_3)  # Number of overlaps


################################################################################################
######################Outline for document#########################
# Step 1: 
#1A - Pre-processing - conversion of .vcf file to genotype matrix
#1B - Convert genotype matrix to rrBLUP format
#1C - Transpose genotype matrix for downstream use
#1D - Genotype import & cleaning/imputation 
####################################################################
# Step 2:
# Genotype import
# Phenotype import
# Set training and test populations
# Run GS
####################################################################
# Step 3:
# Prediction accuracy using
# rrBLUP
# Gaussian Kernel
# Exponential Kernel
####################################################################

####################################################################
#1A FROM VCF TO GENOTYPE MATRIX - rows=individuals  columns =SNPs, in 0/1/2 format
####################################################################
library(vcfR)
library(adegenet)
vcf <- read.vcfR("~/R/World_Veg_Collab_Pepper/filtered_data/pepper_10000.vcf.gz")
genotype_matrix <- extract.gt(vcf, element = "GT")

dim(genotype_matrix)  # Should return (SNPs × Samples)
write.csv(genotype_matrix, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_10k.csv", row.names = TRUE)

genotype_matrix_raw <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/genotype_matrix_10k.csv", 
                                row.names = 1, check.names = FALSE)

# Check dimensions
dim(genotype_matrix_raw)  # Should match (SNPs × Samples)

# Preview 
head(genotype_matrix_raw[, 1:5])  #first 5 samples

####################################################################
#1B - Convert genotype matrix to rrBLUP format
####################################################################
# Function to convert genotypes to numeric format for rrBLUP
convert_geno <- function(x) {
  ifelse(x == "0/0", 0, 
         ifelse(x == "0/1" | x == "1/0", 1, 
                ifelse(x == "1/1", 2, NA)))
}

# Apply conversion to the entire matrix
geno_numeric <- apply(genotype_matrix_raw, 2, convert_geno)

# Convert to dataframe
geno_df <- as.data.frame(geno_numeric)

# Check dimensions after conversion
#dim(geno_df)  #(340734 × 423)
dim(geno_df)  #(23462 × 10250)
# Save the converted genotype matrix
write.csv(geno_df, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format_10k_samples.csv", row.names = TRUE)


################################################################################################
####################################################################
#1C - Transpose genotype matrix for downstream use
####################################################################
library(rrBLUP)

# genotype matrix in rrBLUP format
geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format_10k_samples.csv", 
                    row.names = 1, check.names = FALSE)

# Check dimensions to confirm SNPs x Samples
dim(geno_df)

# Preview first few values
head(geno_df[, 1:5])

geno_df <- t(geno_df)  # Transpose so genotypes (samples) become rows
geno_df <- as.data.frame(geno_df)  # Convert back to a dataframe
head(geno_df[, 1:5])
cat("Corrected genotype matrix dimensions:", dim(geno_df), "\n")


write.csv(geno_df, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format_10k_samples_transposed.csv", row.names = TRUE)


##################
#1D - GENOTYPE IMPORT & CLEANING/IMPUTATION
##################

geno_df <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/genotype_matrix_rrBLUP_format_10k_samples_transposed.csv", 
                    row.names = 1, check.names = FALSE)

# Count missing values before imputation
num_missing_before <- sum(is.na(geno_df))


# Define imputation function (replace missing values with the most common allele)
impute_geno <- function(x) {
  x[is.na(x)] <- as.numeric(names(which.max(table(x, useNA = "no"))))  # Replace NA with mode
  return(x)
}

# Apply the imputation function to each SNP column
cat("Imputing missing values in genotype data...\n")
geno_df_imputed <- apply(geno_df, 2, impute_geno)

# Convert back to dataframe
geno_df_imputed <- as.data.frame(geno_df_imputed)

# Confirm missing values after imputation
num_missing_after <- sum(is.na(geno_df_imputed))
cat("Missing values after imputation:", num_missing_after, "\n") # Should be 0

# Save the cleaned genotype matrix
write.csv(geno_df_imputed, "~/R/World_Veg_Collab_Pepper/filtered_data/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv", row.names = TRUE)

# Total missing values
num_missing_before <- sum(is.na(geno_df))
total_elements <- prod(dim(geno_df))
missing_pct <- (num_missing_before / total_elements) * 100

cat("Missing values before imputation:", num_missing_before, "\n")
cat("Percent missing:", round(missing_pct, 2), "%\n")


################################################################################################
####################################################################
# Step 2:
# Genotype import
# Phenotype import
# Set training and test populations
# Run GS
####################################################################

############################################################################################################
# Step 2
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
control_pheno_check       <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_check       <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_check       <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)

##################
# TRAINING SET FILTERING
##################
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/overlapping_sample_ids_282.csv", header = TRUE)$x

# Filter to only those with complete phenotype data across all three timepoints
complete_lines <- train_genotypes[train_genotypes %in% rownames(control_pheno_check) &
                                    train_genotypes %in% rownames(heat_stress_1_check) &
                                    train_genotypes %in% rownames(heat_stress_2_check)]

complete_lines <- complete_lines[
  rowSums(is.na(control_pheno_check[complete_lines, ])) == 0 &
    rowSums(is.na(heat_stress_1_check[complete_lines, ])) == 0 &
    rowSums(is.na(heat_stress_2_check[complete_lines, ])) == 0
]

train_genotypes <- complete_lines
write.csv(train_genotypes, "~/R/World_Veg_Collab_Pepper/10250_GS/training_set_259_complete.csv", row.names = FALSE)

# Define test set
test_genotypes <- setdiff(rownames(geno_df_imputed), train_genotypes)

##################
# SUBSETTING DATA
##################
geno_train <- geno_df_imputed[train_genotypes, ]
geno_test  <- geno_df_imputed[test_genotypes, ]

rownames(geno_train) <- train_genotypes
rownames(geno_test)  <- test_genotypes

# Convert genotype to numeric matrices
geno_train <- as.matrix(sapply(geno_train, as.numeric))
rownames(geno_train) <- train_genotypes

geno_test  <- as.matrix(sapply(geno_test, as.numeric))
rownames(geno_test)  <- test_genotypes

# Subset phenotype data
control_train       <- control_pheno_check[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_check[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_check[train_genotypes, ]

##################
# rrBLUP FUNCTION
##################
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

##################
# RUN MODEL FOR EACH CONDITION
##################
run_rrBLUP_all_traits(geno_train, geno_test, control_train,       "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k")
run_rrBLUP_all_traits(geno_train, geno_test, heat_stress_1_train, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k")
run_rrBLUP_all_traits(geno_train, geno_test, heat_stress_2_train, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k")


################################################################################################
# Read files
control_train <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_Train.csv", row.names = 1)
control_test  <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_Test.csv",  row.names = 1)

# Add 'Line' column
control_train$Line <- rownames(control_train)
control_test$Line  <- rownames(control_test)

# Move Line to the front
control_train <- control_train[, c("Line", setdiff(names(control_train), "Line"))]
control_test  <- control_test[, c("Line", setdiff(names(control_test), "Line"))]

# Strip _Train and _Test suffixes from colnames (except 'Line')
colnames(control_train)[-1] <- gsub("_Train$", "", colnames(control_train)[-1])
colnames(control_test)[-1]  <- gsub("_Test$", "",  colnames(control_test)[-1])

# Double-check that names now match
stopifnot(all(colnames(control_train) == colnames(control_test)))

# Now safely merge
control_all <- rbind(control_train, control_test)

# Save
write.csv(control_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL.csv", row.names = FALSE)


# Read files
hs1_train <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_Train.csv", row.names = 1)
hs1_test  <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_Test.csv",  row.names = 1)

# Add 'Line' column
hs1_train$Line <- rownames(hs1_train)
hs1_test$Line  <- rownames(hs1_test)

# Move Line to the front
hs1_train <- hs1_train[, c("Line", setdiff(names(hs1_train), "Line"))]
hs1_test  <- hs1_test[,  c("Line", setdiff(names(hs1_test),  "Line"))]

# Strip _Train and _Test
colnames(hs1_train)[-1] <- gsub("_Train$", "", colnames(hs1_train)[-1])
colnames(hs1_test)[-1]  <- gsub("_Test$",  "", colnames(hs1_test)[-1])

# Confirm columns match
stopifnot(all(colnames(hs1_train) == colnames(hs1_test)))

# Merge
hs1_all <- rbind(hs1_train, hs1_test)

# Save
write.csv(hs1_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_ALL.csv", row.names = FALSE)




# Read files
hs2_train <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_Train.csv", row.names = 1)
hs2_test  <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_Test.csv",  row.names = 1)

# Add 'Line' column
hs2_train$Line <- rownames(hs2_train)
hs2_test$Line  <- rownames(hs2_test)

# Move Line to the front
hs2_train <- hs2_train[, c("Line", setdiff(names(hs2_train), "Line"))]
hs2_test  <- hs2_test[,  c("Line", setdiff(names(hs2_test),  "Line"))]

# Strip _Train and _Test
colnames(hs2_train)[-1] <- gsub("_Train$", "", colnames(hs2_train)[-1])
colnames(hs2_test)[-1]  <- gsub("_Test$",  "", colnames(hs2_test)[-1])

# Confirm columns match
stopifnot(all(colnames(hs2_train) == colnames(hs2_test)))

# Merge
hs2_all <- rbind(hs2_train, hs2_test)

# Save
write.csv(hs2_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_ALL.csv", row.names = FALSE)

#######################################################################
# Load your GEBV merged file
control_all <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL.csv")

# Load Tripodi reference file
tripodi <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/10038_Tripodi_2021.csv")

# Check column names to find the identifier column
colnames(tripodi)

# Extract the relevant columns as character vectors
control_lines <- as.character(control_all$Line)
tripodi_lines <- as.character(tripodi$Sample_ID_G2P_Sol_database)

# Fix the Line column: remove a trailing '0' *if* it exists
control_all$Line_fixed <- sub("0$", "", control_all$Line)

# Check overlap
tripodi_ids <- as.character(tripodi$Sample_ID_G2P_Sol_database)
n_matches <- sum(control_all$Line_fixed %in% tripodi_ids)
cat("Number of matches after removing trailing zero:", n_matches, "\n")

# View matching lines (optional)
matched_lines <- intersect(control_all$Line_fixed, tripodi_ids)
head(matched_lines)


# Replace Line with the corrected version
control_all$Line <- control_all$Line_fixed
control_all$Line_fixed <- NULL

# Save the corrected file
write.csv(control_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL_fixed.csv", row.names = FALSE)

#######################################################################
# Read in files
gebv <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL_fixed.csv")
tripodi <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/10038_Tripodi_2021.csv")

# Ensure character type
gebv$Line_fixed <- as.character(gebv$Line_fixed)
tripodi$Sample_ID_G2P_Sol_database <- as.character(tripodi$Sample_ID_G2P_Sol_database)

# Merge on matching IDs
merged <- merge(gebv, tripodi[, c("Sample_ID_G2P_Sol_database", "Organism")],
                by.x = "Line_fixed", by.y = "Sample_ID_G2P_Sol_database",
                all.x = TRUE)

# Save the result
write.csv(merged, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL_fixed_with_organism.csv", row.names = FALSE)



library(ggplot2)
library(dplyr)

# Load data
gebv_data <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL_fixed_with_organism.csv")

# Filter NAs in Organism
gebv_data <- gebv_data %>% filter(!is.na(Organism))

# Create labels with sample sizes
organism_counts <- gebv_data %>%
  group_by(Organism) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Organism, "\n(n=", n, ")"))

# Join back for custom x-axis labels
gebv_data <- gebv_data %>%
  left_join(organism_counts, by = "Organism")

# Plot
ggplot(gebv_data, aes(x = label, y = GEBV_yield, fill = Organism)) +
  geom_boxplot() +
  labs(title = "GEBV for Yield by Species",
       x = "Organism",
       y = "GEBV (Yield)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Hide legend since we now label on axis


#######################################################################

# Load GEBV Heat 1 file
hs1_all <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_ALL.csv")

# Load Tripodi metadata
tripodi <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/10038_Tripodi_2021.csv")

# Fix Line column: remove trailing 0 if present
hs1_all$Line_fixed <- sub("0$", "", hs1_all$Line)

# Replace Line column and remove temporary column
hs1_all$Line <- hs1_all$Line_fixed
hs1_all$Line_fixed <- NULL

# Save fixed file
write.csv(hs1_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_ALL_fixed.csv", row.names = FALSE)

# Merge in Organism info
hs1_all <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_ALL_fixed.csv")
tripodi$Sample_ID_G2P_Sol_database <- as.character(tripodi$Sample_ID_G2P_Sol_database)
hs1_all$Line <- as.character(hs1_all$Line)

merged_hs1 <- merge(hs1_all, tripodi[, c("Sample_ID_G2P_Sol_database", "Organism")],
                    by.x = "Line", by.y = "Sample_ID_G2P_Sol_database", all.x = TRUE)

# Save
write.csv(merged_hs1, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_ALL_fixed_with_organism.csv", row.names = FALSE)


#######################################################################
# Load GEBV Heat 2 file
hs2_all <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_ALL.csv")

# Fix Line column
hs2_all$Line_fixed <- sub("0$", "", hs2_all$Line)
hs2_all$Line <- hs2_all$Line_fixed
hs2_all$Line_fixed <- NULL

# Save fixed file
write.csv(hs2_all, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_ALL_fixed.csv", row.names = FALSE)

# Merge with Organism info
hs2_all <- read.csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_ALL_fixed.csv")
hs2_all$Line <- as.character(hs2_all$Line)

merged_hs2 <- merge(hs2_all, tripodi[, c("Sample_ID_G2P_Sol_database", "Organism")],
                    by.x = "Line", by.y = "Sample_ID_G2P_Sol_database", all.x = TRUE)

# Save
write.csv(merged_hs2, "~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_ALL_fixed_with_organism.csv", row.names = FALSE)
#######################################################################
#PLOTS
#######################################################################
library(dplyr)
library(ggplot2)
library(readr)

# Load and label each condition
control <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Control_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Control")

heat1 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_1_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Heat Stress 1")

heat2 <- read_csv("~/R/World_Veg_Collab_Pepper/outputs/GEBVs_all_73_Heat_Stress_2_10k_ALL_fixed_with_organism.csv") %>%
  mutate(Condition = "Heat Stress 2")

# Combine them
all_gebv <- bind_rows(control, heat1, heat2)

# Remove NAs in Organism
all_gebv <- all_gebv %>% filter(!is.na(Organism))

# Add (n=) labels
organism_labels <- all_gebv %>%
  group_by(Organism) %>%
  summarise(n = n()) %>%
  mutate(Organism_labeled = paste0(Organism, "\n(n=", n, ")"))

# Join back to label organisms
all_gebv <- left_join(all_gebv, organism_labels, by = "Organism")


ggplot(all_gebv, aes(x = Organism_labeled, y = GEBV_yield, fill = Organism)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Condition) +
  labs(title = "GEBV for Yield by Species and Condition",
       x = "Organism",
       y = "GEBV (Yield)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

############################################################################################################
# Step 3: Model Prediction Accuracies using rrBLUP + k-fold CV
############################################################################################################
#Model #1 - rrBLUP
############################################################################################################
library(rrBLUP)
library(dplyr)

###########################
# Step 1: Load Genotype Data
###########################

geno_df <- read.csv(
  "~/R/World_Veg_Collab_Pepper/10250_GS/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv",
  row.names = 1, check.names = FALSE
)

geno_matrix <- as.matrix(sapply(geno_df, as.numeric))
rownames(geno_matrix) <- rownames(geno_df)

###########################
# Step 2: Load Phenotype Data
###########################

control_pheno       <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)

###########################
# Step 3: Use Final 259 Training Set from Earlier
###########################

train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/training_set_259_complete.csv", header = TRUE)$x
test_genotypes  <- setdiff(rownames(geno_matrix), train_genotypes)

geno_train <- geno_matrix[train_genotypes, ]
geno_test  <- geno_matrix[test_genotypes, ]

control_train       <- control_pheno[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_pheno[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_pheno[train_genotypes, ]

###########################
# Step 4: Load Custom Cross-validation Function
###########################

#source("~/R/World_Veg_Collab_Pepper/scripts/xval_kfold_functions.R")  # <-- update if stored elsewhere
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R") 
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R")


###########################
# Step 5: Run k-fold Cross Validation
###########################

# You can change y.trainset to be any of the 3 conditions (run separately)
xval_k10_rrblup_control <- k.xval(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k.fold = 10,
  reps = 50
)

# Save the result
saveRDS(xval_k10_rrblup_control, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_rrblup_kfold_control_10k.RDS")


# Cross-validation for Heat Stress 1
xval_k10_rrblup_heat1 <- k.xval(
  g.in = geno_train,
  y.in = heat_stress_1_train,    # <- use training phenotype subset
  y.trainset = heat_stress_1_train,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_rrblup_heat1, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_rrblup_kfold_heat1_10k.RDS")

# Cross-validation for Heat Stress 2
xval_k10_rrblup_heat2 <- k.xval(
  g.in = geno_train,
  y.in = heat_stress_2_train,    # <- use training phenotype subset
  y.trainset = heat_stress_2_train,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_rrblup_heat2, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_rrblup_kfold_heat2_10k.RDS")


####################
#Plotting Prediction Accuracies
####################
# Load cross-validation results
rrblup_kfold10 <- readRDS("~/R/World_Veg_Collab_Pepper/10250_GS/xval_rrblup_kfold_control_10k.RDS")$xval.result
rrblup_kfold10 <- readRDS("~/R/World_Veg_Collab_Pepper/10250_GS/xval_rrblup_kfold_heat1_10k.RDS")$xval.result
rrblup_kfold10 <- readRDS("~/R/World_Veg_Collab_Pepper/10250_GS/xval_rrblup_kfold_heat2_10k.RDS")$xval.result

# Assign model name
rrblup_kfold10$model <- "rrBLUP"

# Combine (if adding others later)
all_models <- bind_rows(rrblup_kfold10)

# Ensure numeric types
all_models$r.mean <- as.numeric(all_models$r.mean)
all_models$r.sd    <- as.numeric(all_models$r.sd)

# Subset to all traits
all_gs_traits <- all_models %>%
  filter(trait %in% unique(all_models$trait))

# Plot
ggplot(all_gs_traits, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
                width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 5) +
  geom_hline(yintercept = 0.5, color = "red", linewidth = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Prediction Accuracy (r) for 73 Traits",
       x = "Model",
       y = "Prediction Accuracy (r)")


######################################################################################################
#Model #2
######################################################################################################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

# STEP 1: Load Genotype Data
geno_df <- read.csv(
  "~/R/World_Veg_Collab_Pepper/10250_GS/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv",
  row.names = 1, check.names = FALSE
)
geno_matrix <- as.matrix(sapply(geno_df, as.numeric))
rownames(geno_matrix) <- rownames(geno_df)

# STEP 2: Load Phenotype Data
control_pheno       <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)

# STEP 3: Load the 259 Final Training Set
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/training_set_259_complete.csv")$x
test_genotypes <- setdiff(rownames(geno_matrix), train_genotypes)

geno_train <- geno_matrix[train_genotypes, ]
geno_test  <- geno_matrix[test_genotypes, ]

# Subset phenotype datasets
control_train       <- control_pheno[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_pheno[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_pheno[train_genotypes, ]

# Add Genotype column required for kin.blup-based functions
control_train$Genotype <- rownames(control_train)
heat_stress_1_train$Genotype <- rownames(heat_stress_1_train)
heat_stress_2_train$Genotype <- rownames(heat_stress_2_train)

control_train <- control_train %>% relocate(Genotype)
heat_stress_1_train <- heat_stress_1_train %>% relocate(Genotype)
heat_stress_2_train <- heat_stress_2_train %>% relocate(Genotype)

# STEP 4: Load your custom k-fold function
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R")

# STEP 5: Calculate Gaussian Kernel
K <- A.mat(geno_train)
k_dist <- dist(K)

# STEP 6: Cross-validation

# CONTROL
xval_k10_gauss_control <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_gauss_control, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_gauss_kfold_control_10k.RDS")

# HEAT STRESS 1
xval_k10_gauss_heat1 <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = heat_stress_1_train,
  y.trainset = heat_stress_1_train,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_gauss_heat1, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_gauss_kfold_heat1_10k.RDS")

# HEAT STRESS 2
xval_k10_gauss_heat2 <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = heat_stress_2_train,
  y.trainset = heat_stress_2_train,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_gauss_heat2, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_gauss_kfold_heat2_10k.RDS")


######################################################################################################
#Model #3
######################################################################################################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

###########################
# Step 1: Load Genotype Data
###########################
geno_df <- read.csv(
  "~/R/World_Veg_Collab_Pepper/10250_GS/genotype_matrix_rrBLUP_format_10k_samples_transposed_imputed.csv",
  row.names = 1, check.names = FALSE
)
geno_matrix <- as.matrix(sapply(geno_df, as.numeric))
rownames(geno_matrix) <- rownames(geno_df)

###########################
# Step 2: Load Phenotype Data
###########################
control_pheno       <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/control_pheno_filtered_286.csv", row.names = 1)
heat_stress_1_pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_1_filtered_286.csv", row.names = 1)
heat_stress_2_pheno <- read.csv("~/R/World_Veg_Collab_Pepper/filtered_data/heat_stress_2_filtered_286.csv", row.names = 1)

###########################
# Step 3: Use Final 259 Training Set from Earlier
###########################
train_genotypes <- read.csv("~/R/World_Veg_Collab_Pepper/10250_GS/training_set_259_complete.csv", header = TRUE)$x
test_genotypes  <- setdiff(rownames(geno_matrix), train_genotypes)

geno_train <- geno_matrix[train_genotypes, ]
geno_test  <- geno_matrix[test_genotypes, ]

control_train       <- control_pheno[train_genotypes, ]
heat_stress_1_train <- heat_stress_1_pheno[train_genotypes, ]
heat_stress_2_train <- heat_stress_2_pheno[train_genotypes, ]

############################
# Step 4: Add Genotype Column (Filtered + Required for kinship-based models)
###########################

control_train$Genotype       <- rownames(control_train)
heat_stress_1_train$Genotype <- rownames(heat_stress_1_train)
heat_stress_2_train$Genotype <- rownames(heat_stress_2_train)

control_train       <- control_train %>% relocate(Genotype)
heat_stress_1_train <- heat_stress_1_train %>% relocate(Genotype)
heat_stress_2_train <- heat_stress_2_train %>% relocate(Genotype)


###########################
# Step 5: Load Custom Cross-validation Function
###########################
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R")

###########################
# Step 6: Compute Exponential Kernel
###########################
K.Exp <- Kernel_computation(X = geno_train, name = "exponential", degree = NULL, nL = NULL)
exp_dist_mat <- as.matrix(dist(K.Exp))
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)

###########################
# Step 7: Run k-fold Cross Validation for Each Timepoint
###########################
# Control
xval_k10_EXP_control <- k.xval.EXP(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_EXP_control, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_EXP_kfold_control_10k.RDS")

# Heat Stress 1
xval_k10_EXP_heat1 <- k.xval.EXP(
  g.in = geno_train,
  y.in = heat_stress_1_train,
  y.trainset = heat_stress_1_train,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_EXP_heat1, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_EXP_kfold_heat1_10k.RDS")

# Heat Stress 2
xval_k10_EXP_heat2 <- k.xval.EXP(
  g.in = geno_train,
  y.in = heat_stress_2_train,
  y.trainset = heat_stress_2_train,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)
saveRDS(xval_k10_EXP_heat2, "~/R/World_Veg_Collab_Pepper/10250_GS/xval_EXP_kfold_heat2_10k.RDS")

################################################################################################################################
library(dplyr)
library(ggplot2)

################
# Load rrBLUP results
################
rrblup_control <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_rrblup_kfold_control_10k.RDS")$xval.result %>%
  mutate(model = "rrBLUP", condition = "Control")

rrblup_heat1 <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_rrblup_kfold_heat1_10k.RDS")$xval.result %>%
  mutate(model = "rrBLUP", condition = "Heat1")

rrblup_heat2 <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_rrblup_kfold_heat2_10k.RDS")$xval.result %>%
  mutate(model = "rrBLUP", condition = "Heat2")

################
# Load Gauss results
################
gauss_control <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_gauss_kfold_control_10k.RDS")$xval.result %>%
  mutate(model = "Gauss", condition = "Control")

gauss_heat1 <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_gauss_kfold_heat1_10k.RDS")$xval.result %>%
  mutate(model = "Gauss", condition = "Heat1")

gauss_heat2 <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_gauss_kfold_heat2_10k.RDS")$xval.result %>%
  mutate(model = "Gauss", condition = "Heat2")

################
# Load EXP (Exponential Kernel) results (your new model)
################
exp_control <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_EXP_kfold_control_10k.RDS")$xval.result %>%
  mutate(model = "EXP", condition = "Control")

exp_heat1 <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_EXP_kfold_heat1_10k.RDS")$xval.result %>%
  mutate(model = "EXP", condition = "Heat1")

exp_heat2 <- readRDS("~/R/World_Veg_Project/filtered_data_2/10250_GS/xval_EXP_kfold_heat2_10k.RDS")$xval.result %>%
  mutate(model = "EXP", condition = "Heat2")

################
# Combine all models together
################
all_models <- bind_rows(
  rrblup_control, rrblup_heat1, rrblup_heat2,
  gauss_control, gauss_heat1, gauss_heat2,
  exp_control, exp_heat1, exp_heat2
)

# Convert to numeric
all_models$r.mean <- as.numeric(all_models$r.mean)
all_models$r.sd    <- as.numeric(all_models$r.sd)

################
# Plot function
################
plot_by_condition <- function(df, condition_name) {
  df %>%
    filter(condition == condition_name) %>%
    ggplot(aes(y = r.mean, x = model, color = model)) +
    theme_bw() +
    geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd),
                  width = 0.3, position = position_dodge(0.3)) +
    geom_point(size = 3, position = position_dodge(0.3)) +
    facet_wrap(vars(trait), scales = "free_y", ncol = 15) +
    geom_hline(yintercept = 0.5, color = "red", linewidth = 1.2, linetype = "dashed") +
    ylim(0, 1.0) +
    theme(axis.text.x = element_text(angle = 90, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 8)) +
    labs(title = paste("Prediction Accuracy for", condition_name),
         x = "Model", y = "Prediction Accuracy (r)")
}

################
# Plot each separately
################
########################################################################################################################
# Supplemental Figure S12
########################################################################################################################
#Figure S12. 73 trait prediction accuracies for control timepoint in the global collection.

#plot_by_condition(all_models, "Control")
p <- plot_by_condition(all_models, "Control")

ggsave(
  "~/R/World_Veg_Project/PANNELS/prediction_accuracy_plot_Control_10k.pdf",
  plot = p,
  width = 24,
  height = 13,
  units = "in"
)

########################################################################################################################
# Supplemental Figure S13
########################################################################################################################
#Figure S13. 73 trait prediction accuracies for heat stress #1 timepoint in the global collection.

plot_by_condition(all_models, "Heat1")

p <- plot_by_condition(all_models, "Heat1")

ggsave(
  "~/R/World_Veg_Project/PANNELS/prediction_accuracy_plot_Heat1_10k.pdf",
  plot = p,
  width = 24,
  height = 13,
  units = "in"
)

########################################################################################################################
# Supplemental Figure S14
########################################################################################################################
#Figure S14. 73 trait prediction accuracies for heat stress #2 timepoint in the global collection.

#plot_by_condition(all_models, "Heat2")

p <- plot_by_condition(all_models, "Heat2")

ggsave(
  "~/R/World_Veg_Project/PANNELS/prediction_accuracy_plot_Heat2_10k.pdf",
  plot = p,
  width = 24,
  height = 13,
  units = "in"
)

########################################################################################################################
# Supplemental Figure S15
########################################################################################################################
#Figure S15. Global results for the 17 traits with prediction accuracies above 0.5 for rrBLUP.

#load above data
rr_only <- all_models %>%
  filter(model == "rrBLUP")

rr_traits_all3 <- rr_only %>%
  group_by(trait) %>%
  summarise(min_r = min(r.mean, na.rm = TRUE)) %>%
  filter(min_r >= 0.5) %>%
  pull(trait)

rr_traits_all3

rr_filtered <- rr_only %>% filter(trait %in% rr_traits_all3)

p_rr <- rr_filtered %>%
  ggplot(aes(x = condition, y = r.mean)) +
  theme_bw() +
  
  # 0.5 threshold line
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 1) +
  
  # ERROR BARS
  geom_errorbar(aes(ymin = r.mean - r.sd,
                    ymax = r.mean + r.sd),
                width = 0.25, linewidth = 0.8) +
  
  # POINTS ONLY (rrBLUP)
  geom_point(size = 3, color = "#1F77B4") +
  
  facet_wrap(~ trait, ncol = 5, scales = "free_y") +
  ylim(0, 1.0) +
  labs(
    title = "rrBLUP Prediction Accuracy ≥ 0.5 in All Three Conditions",
    x = "Condition",
    y = "Prediction Accuracy (r)"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "none"
  )

p_rr

ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/S15_rrBLUP_traits_all3_above_0.5_10k.pdf",
  plot = p_rr,
  width = 18,
  height = 10,
  units = "in"
)


########################################################################################################################
# Supplemental Figure S16
########################################################################################################################
#Figure S16. 23 quality trait prediction accuracies for chili core collection (18 traits over 0.5).

##########################################
#McLeod Phenotype data - 23 traits (for n=329 lines)
##########################################
#Downloaded pheno data from 
#https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/S7JVEM
pheno_data <- read.csv("~/R/World_Veg_Collab_Pepper/McLeod_pheno_GS/G2P-SOL_Pepper_CC_Pheno_data_2023.csv", sep = ";")
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
install.packages("rrBLUP")
library(rrBLUP)
library(dplyr)

##################
# GENOTYPE IMPORT (IMPUTED)
##################
geno_df_imputed <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv",
                            row.names = 1, check.names = FALSE)

##################
# PHENOTYPE IMPORT
##################
pheno_McLeod <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)

training_set_lines <-read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", header = TRUE)

# Get vector of training line names from genotype file
training_line_names <- training_set_lines[[1]]

# Subset phenotype data to match those training lines
pheno_train <- pheno_McLeod[rownames(pheno_McLeod) %in% training_line_names, ] #gives 132 overlap or 31%

##################
# TRAINING SET FILTERING
##################
# Train = lines that have phenotype data
train_genotypes <- rownames(pheno_train)

# Test = all remaining lines without phenotype data
test_genotypes <- setdiff(rownames(geno_df_imputed), train_genotypes)


geno_train <- geno_df_imputed[train_genotypes, ]
geno_test  <- geno_df_imputed[test_genotypes, ]

# Make sure matrices are numeric
geno_train <- as.matrix(sapply(geno_train, as.numeric))
rownames(geno_train) <- train_genotypes

geno_test  <- as.matrix(sapply(geno_test, as.numeric))
rownames(geno_test)  <- test_genotypes



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



run_rrBLUP_all_traits(geno_train, geno_test, pheno_train, "/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_23traits_n423")
################

# Load both files
gebv_train <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_23traits_n423_Train.csv", row.names = 1)
gebv_test  <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_23traits_n423_Test.csv", row.names = 1)

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
write.csv(gebv_all, "/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv", row.names = TRUE)


################
#PLOTTING 23X GEBV HEATMAP
################
# Load your combined 23-trait GEBVs
gebv_23traits <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/outputs/GEBVs_McLeod_23traits_n423_ALL.csv", row.names = 1)


library(pheatmap)

# Drop non-numeric columns if needed
gebv_mat <- gebv_23traits %>% select(where(is.numeric))

# Optionally scale the data (row-wise)
pheatmap(gebv_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 5,
         main = "Heatmap of GEBVs across Traits and Lines")



####SCALE BY GEBV RANGES TO COMPARE
# Standardize the matrix (z-score by column)
scaled_mat <- scale(gebv_mat)

pheatmap(scaled_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Z-scored GEBVs across Traits and Lines")


#########################################################################################################
#Step 3
#########################################################################################################
#Model Prediction Accuracies - ran on HPC
#########################################################################################################

####################
#Method one: rrBLUP
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
#install.packages("hibayes")
library(ggplot2)


####################
#STEP ONE
####################
# Load genotype matrix (imputed, same as before)
geno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)

# Load phenotype matrix for 23 traits
pheno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)

# Load training set line names
train_genotypes <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", header = TRUE)[[1]]

# Filter genotype matrix to match phenotype
common_genotypes <- intersect(rownames(geno_data), rownames(pheno_data))
geno_filtered <- geno_data[common_genotypes, ]
pheno_filtered <- pheno_data[common_genotypes, ]

geno_matrix <- as.matrix(geno_filtered)
mode(geno_matrix) <- "numeric"

# Filter training genotypes to only those present in the data
train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))

# Subset genotype matrices for training and testing
geno_train <- geno_matrix[train_genotypes_filtered, ]
geno_test  <- geno_matrix[setdiff(rownames(geno_matrix), train_genotypes_filtered), ]

# Subset phenotype matrix to match training set
pheno_train <- pheno_filtered[train_genotypes_filtered, ]
pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

####################
#STEP TWO
####################

pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

####################
#STEP THREE
####################
# Load custom cross-validation function
#source("/home/ahmccorm/kantar_koastore/anna/World_Veg_Collab_Pepper/Genomic_Selection/PredictionAccuracies/xval_kfold_functions.R")
source("/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

# Drop Genotype column before cross-validation
pheno_train_numeric <- pheno_train %>% select(-Genotype)

# Run cross-validation
xval_k10_rrblup <- k.xval(
  g.in = geno_train,
  y.in = pheno_train_numeric,
  y.trainset = pheno_train_numeric,
  k.fold = 10,
  reps = 50
)
# Save output
saveRDS(xval_k10_rrblup, "/Users/annamccormick/R/World_Veg_Project/data/outputs/xval_rrblup_23traits_kfold_10.RDS")


####################
#Method Two: Gaussian Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

# Load genotype and phenotype data
geno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv", row.names = 1, check.names = FALSE)
pheno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", row.names = 1)
train_genotypes <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", header = TRUE)[[1]]

# Filter and format
common_genotypes <- intersect(rownames(geno_data), rownames(pheno_data))
geno_matrix <- as.matrix(geno_data[common_genotypes, ])
mode(geno_matrix) <- "numeric"
pheno_filtered <- pheno_data[common_genotypes, ]

train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))
geno_train <- geno_matrix[train_genotypes_filtered, ]
pheno_train <- pheno_filtered[train_genotypes_filtered, ]

# Add 'Genotype' column for kinship models
pheno_train$Genotype <- rownames(pheno_train)
pheno_train <- pheno_train %>% relocate(Genotype)

# Also drop 'Genotype' for matrix-only methods
pheno_train_numeric <- pheno_train %>% select(-Genotype)

####################
source("/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

# Compute Gaussian distance matrix (Euclidean distance on A.mat)
K <- A.mat(geno_train)
k_dist <- dist(K)

# Run k-fold CV
xval_k10_GAUSS <- k.xval.GAUSS(
  g.in = geno_train,
  y.in = pheno_train,
  y.trainset = pheno_train,
  k_dist = k_dist,
  k.fold = 10,
  reps = 50
)

saveRDS(xval_k10_GAUSS, "/Users/annamccormick/R/World_Veg_Project/data/outputs/xval_GAUSS_kfold_10.RDS")


####################
# Method Three: Exponential Kernel
####################
library(rrBLUP)
library(dplyr)
library(hibayes)
library(ggplot2)

####################
# STEP ONE: Load & Filter Data
####################

# Load genotype matrix
geno_data <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/genotype_matrix_imputed.csv", 
                      row.names = 1, check.names = FALSE)

# Load phenotype matrix (control condition)
control_pheno_check <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/McLeod_per_line_pheno_data_n329_23traits.csv", 
                                row.names = 1)

# Load training genotypes
train_genotypes <- read.csv("/Users/annamccormick/R/World_Veg_Project/data/inputs/Training_Genotypes_n150_run1.csv", 
                            header = TRUE)[[1]]

# Filter for genotypes in both genotype and phenotype datasets
common_genotypes <- intersect(rownames(geno_data), rownames(control_pheno_check))
geno_matrix <- as.matrix(geno_data[common_genotypes, ])
mode(geno_matrix) <- "numeric"
control_pheno_filtered <- control_pheno_check[common_genotypes, ]

# Filter training set
train_genotypes_filtered <- intersect(train_genotypes, rownames(geno_matrix))
geno_train <- geno_matrix[train_genotypes_filtered, ]
control_train <- control_pheno_filtered[train_genotypes_filtered, ]

# Add 'Genotype' column for kin.blup compatibility
control_train$Genotype <- rownames(control_train)
control_train <- control_train %>% relocate(Genotype)

####################
# STEP TWO: Compute Exponential Kernel
####################
source("/Users/annamccormick/R/World_Veg_Project/data/inputs/xval_kfold_functions.R")

# Compute exponential kernel matrix
K.Exp <- Kernel_computation(X = geno_train, name = "exponential", degree = NULL, nL = NULL)
exp_dist <- dist(K.Exp)
exp_dist_mat <- as.matrix(exp_dist)

# Set row/col names explicitly
rownames(exp_dist_mat) <- colnames(exp_dist_mat) <- rownames(geno_train)

# Sanity checks
stopifnot(all(rownames(exp_dist_mat) == rownames(geno_train)))
stopifnot(all(rownames(exp_dist_mat) == control_train$Genotype))

####################
# STEP THREE: Run Cross-Validation
####################
xval_k10_EXP <- k.xval.EXP(
  g.in = geno_train,
  y.in = control_train,
  y.trainset = control_train,
  k_dist = exp_dist_mat,
  k.fold = 10,
  reps = 50
)

# Save output
saveRDS(xval_k10_EXP, "/Users/annamccormick/R/World_Veg_Project/data/outputs/xval_EXP_kfold_10.RDS")




###################################################################################################
#Plot prediction accuracies
###################################################################################################

library(ggplot2)
library(dplyr)
library(readr)

# Load cross-validation results (these should contain the $xval.result list)
rrblup_results <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_rrblup_23traits_kfold_10.RDS")$xval.result
gauss_results  <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_GAUSS_kfold_10.RDS")$xval.result
exp_results    <- readRDS("~/R/World_Veg_Project/data/outputs/PA_23trait_core423/xval_EXP_kfold_10.RDS")$xval.result

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
p<-ggplot(all_models, aes(x = model, y = r.mean, color = model)) +
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

p

ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/prediction_accuracy_23traits_423.pdf",
  plot = p,
  width = 20,
  height = 12,
  units = "in"
)



########################################################################################################################
# Supplemental Figure S17
########################################################################################################################
#Figure S17. 23 quality trait prediction accuracies for chili global 10k collection (17 traits over 0.5).




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
p<-ggplot(all_models, aes(x = model, y = r.mean, color = model)) +
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

p

ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/prediction_accuracy_23traits_10k.pdf",
  plot = p,
  width = 20,
  height = 12,
  units = "in"
)





########################################################################################################################
# Supplemental Figure S18
########################################################################################################################
#Figure S18. Agronomic trait GEBV distributions in core and global collections.

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

p<- ggplot(gebv_filtered, aes(x = GEBV, fill = Group)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Trait + Condition, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Common High-Accuracy Traits: GEBV Distributions by Condition and Group",
    x = "GEBV",
    y = "Density",
    fill = "Group"
  )
ggsave(
  filename = "~/R/World_Veg_Project/PANNELS/GEBV_density_high_accuracy_traits.pdf",
  plot = p,
  width = 18,
  height = 20,
  units = "in"
)




########################################################################################################################
# Supplemental Figure S19
########################################################################################################################
setwd("~/R/World_Veg_Project/anna_heritability/")

library(readxl)
library(dplyr)
library(vcfR)
library(sommer)

#packageVersion("vcfR")
#packageVersion("sommer")

vcf_file <- "Pepper_filtered_and_LD_pruned.vcf.gz"
control_file <- "Control.xlsx"
heat1_file   <- "Heat_stress_1.xlsx"
heat2_file   <- "Heat_stress_2.xlsx"
sheet_name <- "data"

min_maf <- 0.05
max_missing <- 0.10   # match his params


# Read + standardize phenotypes
read_pheno <- function(path, treatment){
  df <- read_excel(path, sheet = sheet_name)
  names(df) <- tolower(names(df))
  
  # Ensure required columns exist
  # Expect: g2pname, rep (or similar)
  if (!("g2pname" %in% names(df))) stop("No g2pname column in ", path)
  
  # standardize rep name to "rep"
  rep_col <- grep("^rep$", names(df), ignore.case = TRUE, value = TRUE)
  if (length(rep_col) == 1) names(df)[names(df) == rep_col] <- "rep"
  
  df$g2pname <- trimws(as.character(df$g2pname))
  df <- df[!is.na(df$g2pname), ]
  
  # convert everything except these to numeric
  exclude <- c("g2pname","rep","genotype","plantno.","treatment")
  for (cn in names(df)) {
    if (!(cn %in% exclude)) {
      df[[cn]] <- suppressWarnings(as.numeric(as.character(df[[cn]])))
    }
  }
  
  df$treatment <- treatment
  df
}

control <- read_pheno(control_file, "Control")
hs1     <- read_pheno(heat1_file,   "HS1")
hs2     <- read_pheno(heat2_file,   "HS2")

all_pheno <- bind_rows(control, hs1, hs2)


# Build GRM 
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT")  # variants x samples

# Convert GT strings -> 0/1/2 with proper handling
gt_numeric <- matrix(NA_real_, nrow = nrow(gt), ncol = ncol(gt),
                     dimnames = dimnames(gt))

# vectorized mapping is faster than nested loops
map_gt <- function(x){
  out <- rep(NA_real_, length(x))
  out[x %in% c("0/0","0|0")] <- 0
  out[x %in% c("0/1","1/0","0|1","1|0")] <- 1
  out[x %in% c("1/1","1|1")] <- 2
  out
}

for (i in seq_len(nrow(gt))) {
  gt_numeric[i, ] <- map_gt(gt[i, ])
  if (i %% 20000 == 0) cat("Processed", i, "variants\n")
}

# Filter SNPs by missingness
missing_rate <- rowMeans(is.na(gt_numeric))
gt_numeric <- gt_numeric[missing_rate <= max_missing, , drop=FALSE]

# Filter SNPs by MAF
calc_maf <- function(x){
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0)
  p <- sum(x) / (2 * length(x))
  min(p, 1 - p)
}
maf <- apply(gt_numeric, 1, calc_maf)
gt_numeric <- gt_numeric[maf >= min_maf, , drop=FALSE]

cat("SNPs retained:", nrow(gt_numeric), "\n")

# Mean-impute missing per SNP
row_means <- rowMeans(gt_numeric, na.rm = TRUE)
for (i in seq_len(nrow(gt_numeric))) {
  miss <- is.na(gt_numeric[i, ])
  if (any(miss)) gt_numeric[i, miss] <- row_means[i]
}

# Samples as rows
M <- t(gt_numeric)

# Scale markers then GRM
M_scaled <- scale(M, center = TRUE, scale = TRUE)
grm <- tcrossprod(M_scaled) / ncol(M_scaled)
diag(grm) <- diag(grm) + 1e-6


# h² function (boss-style): no Rep, no Ve/rbar
calc_h2 <- function(pheno_df, trait, grm){
  dat <- pheno_df[, c("g2pname","rep", trait)]
  colnames(dat) <- c("id","rep","y")
  dat$id <- as.character(dat$id)
  dat$y  <- suppressWarnings(as.numeric(dat$y))
  dat <- dat[!is.na(dat$y), , drop=FALSE]
  
  common <- intersect(unique(dat$id), rownames(grm))
  if (length(common) < 10) return(NA_real_)
  
  dat <- dat[dat$id %in% common, , drop=FALSE]
  grm_sub <- grm[common, common, drop=FALSE]
  dat$id <- factor(dat$id)
  
  fit <- try(mmer(
    y ~ 1,
    random = ~ vs(id, Gu = grm_sub),
    rcov   = ~ units,
    data   = dat,
    verbose = FALSE
  ), silent=TRUE)
  
  if (inherits(fit,"try-error")) return(NA_real_)
  
  # sommer variance comps
  vc <- summary(fit)$varcomp
  g_row <- grep("id", rownames(vc), value=TRUE)
  g_row <- setdiff(g_row, "units")
  if (length(g_row) < 1) return(NA_real_)
  
  Vg <- vc[g_row[1], "VarComp"]
  Ve <- vc["units","VarComp"]
  as.numeric(Vg / (Vg + Ve))
}

# Trait list (numeric columns only; exclude metadata)
exclude <- c("g2pname","rep","genotype","plantno.","treatment")
trait_cols <- names(all_pheno)[sapply(all_pheno, is.numeric)]
trait_cols <- setdiff(trait_cols, exclude)


# Run by treatment, save results
run_treatment <- function(df, trt){
  out <- sapply(trait_cols, calc_h2, pheno_df = df, grm = grm)
  data.frame(trait = trait_cols, treatment = trt, h2 = as.numeric(out))
}

res_control <- run_treatment(filter(all_pheno, treatment=="Control"), "Control")
res_hs1     <- run_treatment(filter(all_pheno, treatment=="HS1"),     "HS1")
res_hs2     <- run_treatment(filter(all_pheno, treatment=="HS2"),     "HS2")

res_all <- bind_rows(res_control, res_hs1, res_hs2)

######################
# A Plot
######################
hist(res_all$h2[is.finite(res_all$h2)], breaks=20,
     main="Distribution of narrow-sense h² (all treatments)",
     xlab="h²")

abline(v = 0.25, col = "red", lwd = 2, lty = 2)

######################
# B plot
######################
boxplot(h2 ~ treatment, data = res_all,
        main="h² by treatment", ylab="h²", xlab="")

abline(h = 0.25, col = "red", lwd = 2, lty = 2)

#################
#pdf save
pdf("Figure_S19_heritability.pdf", width = 10, height = 5)

par(mfrow = c(1, 2))  # 1 row, 2 columns

######################
# A Plot
######################
hist(res_all$h2[is.finite(res_all$h2)], breaks = 20,
     main = "Distribution of narrow-sense h² (all treatments)",
     xlab = "h²")

abline(v = 0.25, col = "red", lwd = 2, lty = 2)

######################
# B Plot
######################
boxplot(h2 ~ treatment, data = res_all,
        main = "h² by treatment", ylab = "h²", xlab = "")

abline(h = 0.25, col = "red", lwd = 2, lty = 2)

dev.off()


########################################################################################################################
# Supplemental Figure S20
########################################################################################################################
# LD decay across the genome in the core collection
vcf_file <- "~/R/World_Veg_Project/anna_LD_heritability/Capsicum_n323_filtered.vcf.gz"  

args <- commandArgs(trailingOnly = TRUE)

show_help <- function() {
  cat("
Generic LD Decay Analysis Script

Usage: Rscript LD_decay_generic.R [options]

VCF file is hardcoded in the script. Edit the 'vcf_file' variable at the top of the script.

Optional:
  --max-distance         Maximum distance for LD calculation in bp (default: 1000000)
  --sample-size          Maximum number of SNP pairs per chromosome (default: 100000)
  --max-snps             Maximum SNPs to sample per chromosome (default: 20000)
  --min-individuals      Minimum individuals required for LD calculation (default: 10)
  --output-prefix        Prefix for output files (default: 'ld_decay')
  --random-seed          Random seed for reproducibility (default: 12345, use NULL for random)
  --help                 Show this help message
\n")
  quit(save = "no", status = 0)
}

if ("--help" %in% args || "-h" %in% args) show_help()

if (!file.exists(vcf_file)) {
  stop("Error: VCF file '", vcf_file, "' not found!\n",
       "Edit 'vcf_file' at top of script.")
}

# Defaults (match your current defaults)
max_distance    <- 1000000
sample_size     <- 100000
max_snps        <- 20000
min_individuals <- 10
output_prefix   <- "ld_decay"
random_seed     <- 12345  # set to NULL for non-reproducible runs

# Parse optional args
i <- 1
while (i <= length(args)) {
  key <- args[i]
  val <- if (i < length(args)) args[i + 1] else NA
  
  if (key == "--max-distance")      { max_distance    <- as.numeric(val); i <- i + 2; next }
  if (key == "--sample-size")       { sample_size     <- as.numeric(val); i <- i + 2; next }
  if (key == "--max-snps")          { max_snps        <- as.numeric(val); i <- i + 2; next }
  if (key == "--min-individuals")   { min_individuals <- as.numeric(val); i <- i + 2; next }
  if (key == "--output-prefix")     { output_prefix   <- val;            i <- i + 2; next }
  if (key == "--random-seed") {
    # allow explicit "NULL"
    if (!is.na(val) && toupper(val) == "NULL") random_seed <- NULL else random_seed <- as.numeric(val)
    i <- i + 2
    next
  }
  
  cat("Warning: Unknown argument '", key, "' ignored\n", sep = "")
  i <- i + 1
}

cat("\n==========================================\n")
cat("LD Decay Analysis Configuration\n")
cat("==========================================\n")
cat("VCF file:           ", vcf_file, "\n")
cat("Max distance:       ", max_distance, " bp\n")
cat("Sample size:        ", sample_size, " SNP pairs per chromosome\n")
cat("Max SNPs sampled:   ", max_snps, " per chromosome\n")
cat("Min individuals:    ", min_individuals, "\n")
cat("Output prefix:      ", output_prefix, "\n")
cat("Random seed:        ", ifelse(is.null(random_seed), "NULL (random)", random_seed), "\n")
cat("==========================================\n\n")

# Packages
required_packages <- c("vcfR", "ggplot2", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package '", pkg, "' is not installed. Install with install.packages('", pkg, "')")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

use_metbrewer <- FALSE
if (requireNamespace("MetBrewer", quietly = TRUE)) {
  suppressPackageStartupMessages(library(MetBrewer))
  use_metbrewer <- TRUE
} else {
  cat("Note: MetBrewer not found. Using base palette.\n")
}

# Read VCF
cat("Reading VCF...\n")
vcf <- tryCatch(read.vcfR(vcf_file, verbose = FALSE),
                error = function(e) stop("Error reading VCF: ", e$message))

n_variants <- nrow(vcf@fix)
n_samples  <- ncol(vcf@gt) - 1
cat("VCF loaded: ", n_variants, " variants, ", n_samples, " samples\n\n")

if (n_variants < 10) stop("Too few variants (", n_variants, "). Need >= 10.")
if (n_samples < min_individuals) stop("Too few samples (", n_samples, "). Need >= ", min_individuals, ".")

cat("Extracting genotypes...\n")
gt_matrix <- tryCatch(extract.gt(vcf, element = "GT", as.numeric = TRUE),
                      error = function(e) stop("Error extracting GT: ", e$message))

chrom_info <- data.frame(
  chromosome = as.character(vcf@fix[, "CHROM"]),
  position   = as.numeric(vcf@fix[, "POS"]),
  snp_id     = seq_len(nrow(vcf@fix)),
  stringsAsFactors = FALSE
)

# Set seed
if (!is.null(random_seed)) {
  set.seed(random_seed)
  cat("Random seed set:", random_seed, "\n\n")
} else {
  cat("No random seed set (non-reproducible)\n\n")
}

calculate_r2 <- function(snp1, snp2, min_ind = 10) {
  valid <- !is.na(snp1) & !is.na(snp2)
  if (sum(valid) < min_ind) return(NA_real_)
  
  x <- snp1[valid]; y <- snp2[valid]
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  
  r2 <- suppressWarnings(cor(x, y)^2)
  if (!is.finite(r2) || r2 < 0 || r2 > 1) return(NA_real_)
  r2
}

estimate_ld_decay_chromosome <- function(chrom_name,
                                         max_dist = max_distance,
                                         max_pairs = sample_size,
                                         max_snp = max_snps,
                                         min_ind = min_individuals) {
  
  cat("Processing chromosome:", chrom_name, "\n")
  chrom_snps <- chrom_info[chrom_info$chromosome == chrom_name, , drop = FALSE]
  if (nrow(chrom_snps) < 10) {
    cat("  Skipping:", chrom_name, "- insufficient SNPs (", nrow(chrom_snps), ")\n", sep="")
    return(NULL)
  }
  
  chrom_gt  <- gt_matrix[chrom_snps$snp_id, , drop = FALSE]
  positions <- chrom_snps$position
  
  n_snps <- nrow(chrom_gt)
  if (n_snps > max_snp) {
    idx <- sample.int(n_snps, max_snp)
    chrom_gt  <- chrom_gt[idx, , drop = FALSE]
    positions <- positions[idx]
  }
  
  n_sampled <- nrow(chrom_gt)
  if (n_sampled < 2) return(NULL)
  
  # maximum possible pairs among sampled SNPs
  max_possible <- (n_sampled * (n_sampled - 1)) / 2
  max_pairs <- min(max_pairs, max_possible)
  
  # accumulate rows in a list (fast)
  out <- vector("list", length = max_pairs)
  k <- 0L
  
  for (i in 1:(n_sampled - 1)) {
    if (k >= max_pairs) break
    for (j in (i + 1):n_sampled) {
      if (k >= max_pairs) break
      d <- abs(positions[j] - positions[i])
      if (d <= max_dist && d > 0) {
        r2 <- calculate_r2(chrom_gt[i, ], chrom_gt[j, ], min_ind = min_ind)
        if (is.finite(r2)) {
          k <- k + 1L
          out[[k]] <- data.frame(
            chromosome = chrom_name,
            distance   = d,
            r2         = r2,
            snp1_pos   = positions[i],
            snp2_pos   = positions[j],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (k == 0L) {
    cat("  Chromosome", chrom_name, ": No valid LD values\n")
    return(NULL)
  }
  
  ld <- do.call(rbind, out[seq_len(k)])
  cat("  Chromosome", chrom_name, ": calculated", nrow(ld), "LD values\n")
  ld
}

# Run per chromosome
unique_chromosomes <- unique(chrom_info$chromosome)
cat("Found", length(unique_chromosomes), "chromosome(s):", paste(unique_chromosomes, collapse = ", "), "\n\n")

all_ld_list <- lapply(unique_chromosomes, estimate_ld_decay_chromosome)
all_ld_list <- all_ld_list[!vapply(all_ld_list, is.null, logical(1))]

if (length(all_ld_list) == 0) stop("No LD values could be calculated. Check VCF/params.")
all_ld_results <- do.call(rbind, all_ld_list)

cat("\nTotal LD calculations:", nrow(all_ld_results), "\n")

# Save raw results

setwd("~/R/World_Veg_Project/anna_LD_heritability/")

raw_results_file <- paste0(output_prefix, "_raw_results.csv")
write.csv(all_ld_results, raw_results_file, row.names = FALSE)
cat("Saved:", raw_results_file, "\n")

# ---- distance bins ----
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

make_bin_labels <- function(brks) {
  lab <- character(length(brks) - 1)
  for (i in seq_len(length(brks) - 1)) {
    a <- brks[i]; b <- brks[i + 1]
    if (b < 1000) lab[i] <- sprintf("%d-%dbp", a, b)
    else if (a < 1000) lab[i] <- sprintf("%dbp-%.1fkb", a, b / 1000)
    else lab[i] <- sprintf("%.1f-%.1fkb", a / 1000, b / 1000)
  }
  lab
}

distance_labels <- make_bin_labels(distance_breaks)
all_ld_results$distance_bin <- cut(all_ld_results$distance,
                                   breaks = distance_breaks,
                                   labels = distance_labels,
                                   include.lowest = TRUE)

# Summary by chromosome/bin
ld_summary <- all_ld_results %>%
  dplyr::filter(!is.na(distance_bin)) %>%
  dplyr::group_by(chromosome, distance_bin) %>%
  dplyr::summarise(
    mean_r2       = mean(r2, na.rm = TRUE),
    median_r2     = median(r2, na.rm = TRUE),
    sd_r2         = sd(r2, na.rm = TRUE),
    n_pairs       = dplyr::n(),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  )

summary_file <- paste0(output_prefix, "_summary.csv")
write.csv(ld_summary, summary_file, row.names = FALSE)
cat("Saved:", summary_file, "\n")

# Genome-wide bin summary
genome_wide_summary <- all_ld_results %>%
  dplyr::filter(!is.na(distance_bin)) %>%
  dplyr::group_by(distance_bin) %>%
  dplyr::summarise(
    mean_r2       = mean(r2, na.rm = TRUE),
    median_r2     = median(r2, na.rm = TRUE),
    sd_r2         = sd(r2, na.rm = TRUE),
    n_pairs       = dplyr::n(),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(chromosome = "Genome-wide")

gw_summary_file <- paste0(output_prefix, "_genome_wide_summary.csv")
write.csv(genome_wide_summary, gw_summary_file, row.names = FALSE)
cat("Saved:", gw_summary_file, "\n")

# Fit decay models
fit_decay_model <- function(df) {
  df <- df[df$mean_r2 > 0 & is.finite(df$mean_r2) & is.finite(df$mean_distance), , drop = FALSE]
  if (nrow(df) < 4) return(NULL)
  
  tryCatch({
    m <- nls(mean_r2 ~ a * exp(-b * mean_distance) + c,
             data = df,
             start = list(a = max(df$mean_r2), b = 1e-6, c = 0.01),
             lower = list(a = 0, b = 1e-10, c = 0),
             algorithm = "port",
             control = nls.control(maxiter = 100, warnOnly = TRUE))
    
    p <- coef(m)
    if (any(!is.finite(p)) || p["a"] <= 0 || p["b"] <= 0) return(NULL)
    
    max_reasonable <- max(df$mean_distance, na.rm = TRUE) * 10
    
    half_decay <- log(2) / p["b"]
    if (!is.finite(half_decay) || half_decay < 0 || half_decay > max_reasonable) half_decay <- NA_real_
    
    r2_0.2 <- NA_real_
    if ((0.2 - p["c"]) > 0 && (0.2 - p["c"]) < p["a"]) {
      r2_0.2 <- -log((0.2 - p["c"]) / p["a"]) / p["b"]
      if (!is.finite(r2_0.2) || r2_0.2 < 0 || r2_0.2 > max_reasonable) r2_0.2 <- NA_real_
    }
    
    list(a = unname(p["a"]), b = unname(p["b"]), c = unname(p["c"]),
         half_decay_distance = half_decay, r2_0.2_distance = r2_0.2)
  }, error = function(e) NULL)
}

# Genome-wide parameters
genome_wide_params <- NULL
gw_fit <- fit_decay_model(genome_wide_summary)
if (!is.null(gw_fit)) {
  genome_wide_params <- data.frame(
    chromosome = "Genome-wide",
    a = gw_fit$a, b = gw_fit$b, c = gw_fit$c,
    half_decay_distance = gw_fit$half_decay_distance,
    r2_0.2_distance     = gw_fit$r2_0.2_distance
  )
  gw_params_file <- paste0(output_prefix, "_genome_wide_parameters.csv")
  write.csv(genome_wide_params, gw_params_file, row.names = FALSE)
  cat("Saved:", gw_params_file, "\n")
}

# Per-chromosome parameters (use list -> bind once)
cat("\nFitting per-chromosome decay models...\n")
chroms <- unique(ld_summary$chromosome)
param_list <- list()

for (chrom in chroms) {
  dat <- ld_summary[ld_summary$chromosome == chrom, , drop = FALSE]
  fit <- fit_decay_model(dat)
  if (is.null(fit)) next
  param_list[[chrom]] <- data.frame(
    chromosome = chrom,
    a = fit$a, b = fit$b, c = fit$c,
    half_decay_distance = fit$half_decay_distance,
    r2_0.2_distance     = fit$r2_0.2_distance,
    stringsAsFactors = FALSE
  )
}

decay_parameters <- if (length(param_list) > 0) do.call(rbind, param_list) else data.frame()

if (nrow(decay_parameters) > 0) {
  params_file <- paste0(output_prefix, "_parameters.csv")
  write.csv(decay_parameters, params_file, row.names = FALSE)
  cat("Saved:", params_file, "\n")
}

# Smooth curves for plotting
smooth_decay_data <- data.frame()
if (nrow(decay_parameters) > 0) {
  for (i in seq_len(nrow(decay_parameters))) {
    chrom <- decay_parameters$chromosome[i]
    p <- decay_parameters[i, ]
    
    max_dist <- max(ld_summary$mean_distance[ld_summary$chromosome == chrom], na.rm = TRUE)
    if (!is.finite(max_dist) || max_dist <= 1000) next
    
    distances <- seq(1000, max_dist, length.out = 100)
    pred <- p$a * exp(-p$b * distances) + p$c
    pred <- pmax(pred, 0)
    
    smooth_decay_data <- rbind(smooth_decay_data, data.frame(
      chromosome = chrom,
      distance = distances,
      predicted_r2 = pred,
      half_decay_distance = p$half_decay_distance
    ))
  }
}

# Palette
n_chroms <- length(unique(ld_summary$chromosome))
color_palette <- if (use_metbrewer && n_chroms > 0) MetBrewer::met.brewer("VanGogh1", n = n_chroms) else grDevices::rainbow(n_chroms)

#Plot helpers 
plot_genomewide <- function(points_df, x_kb, y, smooth_df, params, title, subtitle, outfile, raw = FALSE) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = points_df,
                        ggplot2::aes(x = !!x_kb, y = !!y),
                        color = "darkblue",
                        size = if (raw) 1 else 3,
                        alpha = if (raw) 0.2 else 0.6) +
    ggplot2::geom_line(data = smooth_df,
                       ggplot2::aes(x = distance/1000, y = predicted_r2),
                       color = if (raw) "red" else "steelblue",
                       linewidth = 1.4) +
    ggplot2::geom_hline(yintercept = 0.2, color = "orange", linetype = "dotted", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = params$half_decay_distance/1000, color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Distance (kb, log scale)", y = if (raw) "r²" else "Mean r²") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray50"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.na(params$r2_0.2_distance)) {
    p <- p + ggplot2::geom_vline(xintercept = params$r2_0.2_distance/1000, color = "orange", linetype = "dashed", linewidth = 1)
  }
  
  ggplot2::ggsave(outfile, p, width = 10, height = 8, dpi = 300)
  cat("Saved:", outfile, "\n")
}


# Genome-wide plots
if (!is.null(genome_wide_params)) {
  max_gw_dist <- max(genome_wide_summary$mean_distance, na.rm = TRUE)
  if (is.finite(max_gw_dist) && max_gw_dist > 1000) {
    gw_distances <- seq(1000, max_gw_dist, length.out = 100)
    gw_predicted <- genome_wide_params$a * exp(-genome_wide_params$b * gw_distances) + genome_wide_params$c
    gw_predicted <- pmax(gw_predicted, 0)
    gw_smooth_data <- data.frame(distance = gw_distances, predicted_r2 = gw_predicted)
    
    subtitle_binned <- paste0(
      "Half-decay: ", round(genome_wide_params$half_decay_distance), " bp",
      if (!is.na(genome_wide_params$r2_0.2_distance)) paste0(" | r²=0.2: ", round(genome_wide_params$r2_0.2_distance), " bp") else "",
      " | ", nrow(genome_wide_summary), " bins"
    )
    
    # ADD THIS LINE
    genome_wide_summary$mean_distance_kb <- genome_wide_summary$mean_distance / 1000
    
    plot_genomewide(
      points_df = genome_wide_summary,
      x_kb = rlang::sym("mean_distance_kb"),   # CHANGED
      y = rlang::sym("mean_r2"),
      smooth_df = gw_smooth_data,
      params = genome_wide_params,
      title = "Genome-wide Average LD Decay (Binned)",
      subtitle = subtitle_binned,
      outfile = paste0(output_prefix, "_genome_wide_1.pdf"),
      raw = FALSE
    )
    
    # raw sampled points plot
    max_raw_points <- 5000
    sampled_raw <- if (nrow(all_ld_results) > max_raw_points) all_ld_results[sample(nrow(all_ld_results), max_raw_points), ] else all_ld_results
    
    subtitle_raw <- paste0(
      "Half-decay: ", round(genome_wide_params$half_decay_distance), " bp",
      if (!is.na(genome_wide_params$r2_0.2_distance)) paste0(" | r²=0.2: ", round(genome_wide_params$r2_0.2_distance), " bp") else "",
      " | ",
      if (nrow(all_ld_results) > max_raw_points) paste0(max_raw_points, " sampled points of ", nrow(all_ld_results)) else paste0(nrow(all_ld_results), " SNP pairs")
    )
    
    points_raw <- sampled_raw
    points_raw$distance_kb <- points_raw$distance / 1000
    
    plot_genomewide(
      points_df = points_raw,
      x_kb = rlang::sym("distance_kb"),        # CHANGED
      y = rlang::sym("r2"),
      smooth_df = gw_smooth_data,
      params = genome_wide_params,
      title = "Genome-wide LD Decay (Raw Data)",
      subtitle = subtitle_raw,
      outfile = paste0(output_prefix, "_genome_wide_2.pdf"),
      raw = TRUE
    )
  }
}





########################################################################################################################
# Supplemental Figure S21
########################################################################################################################
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

###################################################################################################
# SUPPLEMENTAL TABLES S18 and S19
###################################################################################################
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(dplyr)

gff_path <- "/Users/annamccormick/R/World_Veg_Project/Reference_sequence_cornell/Zhangshugang.gff3"
gff <- import(gff_path)

genes <- gff[gff$type == "gene"]
mrna  <- gff[gff$type == "mRNA"]
exons <- gff[gff$type == "exon"]
cds   <- gff[gff$type == "CDS"]

collapse_cl <- function(x) {
  vapply(x, function(z) paste(as.character(z), collapse = "; "), character(1))
}

annotate_snps <- function(snps_df, genes, mrna, exons, cds) {
  
  snps_df <- snps_df %>%
    mutate(Chromosome_gff = sprintf("Chr%02d", as.integer(Chromosome)))
  
  snp_gr <- GRanges(
    seqnames = snps_df$Chromosome_gff,
    ranges   = IRanges(snps_df$Position, snps_df$Position)
  )
  
  # overlaps
  gene_hits <- findOverlaps(snp_gr, genes, ignore.strand = TRUE)
  mrna_hits <- findOverlaps(snp_gr, mrna,  ignore.strand = TRUE)
  exon_hits <- findOverlaps(snp_gr, exons, ignore.strand = TRUE)
  cds_hits  <- findOverlaps(snp_gr, cds,   ignore.strand = TRUE)
  
  # init
  snps_df$gene_id   <- NA_character_
  snps_df$gene_name <- NA_character_
  
  snps_df$mrna_id   <- NA_character_
  snps_df$mrna_note <- NA_character_
  
  snps_df$in_exon <- FALSE
  snps_df$in_cds  <- FALSE
  
  # gene fields (ID/Name exist even if Note doesn't)
  qg <- queryHits(gene_hits); sg <- subjectHits(gene_hits)
  snps_df$gene_id[qg]   <- mcols(genes)$ID[sg]
  snps_df$gene_name[qg] <- mcols(genes)$Name[sg]
  
  # mRNA functional annotation (THIS is where Note lives in your file)
  qm <- queryHits(mrna_hits); sm <- subjectHits(mrna_hits)
  snps_df$mrna_id[qm]   <- mcols(mrna)$ID[sm]
  snps_df$mrna_note[qm] <- collapse_cl(mcols(mrna)$Note[sm])
  
  # exon/CDS flags
  snps_df$in_exon[unique(queryHits(exon_hits))] <- TRUE
  snps_df$in_cds[unique(queryHits(cds_hits))]   <- TRUE
  
  # classification
  snps_df$annotation <- "Intergenic"
  snps_df$annotation[!is.na(snps_df$gene_id)] <- "Genic"
  snps_df$annotation[snps_df$in_exon]         <- "Exonic"
  snps_df$annotation[snps_df$in_cds]          <- "CDS"
  
  # "best available" functional label
  snps_df$functional_note <- snps_df$mrna_note  # primary source in your GFF
  
  # nearest gene + distance (for intergenic SNPs)
  nearest_gene_idx <- nearest(snp_gr, genes, ignore.strand = TRUE)
  snps_df$nearest_gene_id <- mcols(genes)$ID[nearest_gene_idx]
  snps_df$dist_to_nearest_gene_bp <- distance(snp_gr, genes[nearest_gene_idx])
  
  # nearest mRNA note (more interpretable than nearest_gene_id)
  nearest_mrna_idx <- nearest(snp_gr, mrna, ignore.strand = TRUE)
  snps_df$nearest_mrna_note <- collapse_cl(mcols(mrna)$Note[nearest_mrna_idx])
  
  snps_df
}

cold <- read.csv("~/R/World_Veg_Project/top_SNP_annotations/significant_snps_cold57.csv",
                 stringsAsFactors = FALSE)
hot  <- read.csv("~/R/World_Veg_Project/top_SNP_annotations/significant_snps_hot28.csv",
                 stringsAsFactors = FALSE)

cold_annot <- annotate_snps(cold, genes, mrna, exons, cds)
hot_annot  <- annotate_snps(hot,  genes, mrna, exons, cds)

write.csv(cold_annot,
          "~/R/World_Veg_Project/top_SNP_annotations/significant_snps_cold57_annotated.csv",
          row.names = FALSE)
write.csv(hot_annot,
          "~/R/World_Veg_Project/top_SNP_annotations/significant_snps_hot28_annotated.csv",
          row.names = FALSE)


####################################################################################################################################################################
################################################################################# FIN #############################################################################
####################################################################################################################################################################













