##############
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
