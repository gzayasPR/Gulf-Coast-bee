library(tidyverse)
library(ggpubr)


rm(list=ls())

args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
pca_Scores <- args[2]


setwd(out.dir)
pca_results <- read.csv(pca_Scores)

pca_results <- pca_results %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

# Define a color palette
color_palette <- c("AL" = "#4CAF50", "FL" = "#2196F3")

# Base plot theme for consistency
base_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pc12_ANGSD <- ggplot(pca_results, aes(x = PC1_ANGSD, y = PC2_ANGSD, color = Location, shape = Sex)) +
  geom_point(size = 4) +
  labs(title = "PCA ANGSD", x = "Principal Component 1", y = "Principal Component 2") +
  scale_color_manual(values = color_palette, name = "Location") +
  scale_shape_manual(values = c(17, 16), name = "Sex") +
  base_theme

pc12_DeepVariant <- ggplot(pca_results, aes(x = PC1_DeepVariant, y = PC2_DeepVariant, color = Location, shape = Sex)) +
  geom_point(size = 4) +
  labs(title = "PCA DeepVariant", x = "Principal Component 1", y = "Principal Component 2") +
  scale_color_manual(values = color_palette, name = "Location") +
  scale_shape_manual(values = c(17, 16), name = "Sex") +
  base_theme

pc12_GATK <- ggplot(pca_results, aes(x = PC1_GATK, y = PC2_GATK, color = Location, shape = Sex)) +
  geom_point(size = 4) +
  labs(title = "PCA GATK", x = "Principal Component 1", y = "Principal Component 2") +
  scale_color_manual(values = color_palette, name = "Location") +
  scale_shape_manual(values = c(17, 16), name = "Sex") +
  base_theme

combined_plot <- ggarrange(pc12_DeepVariant, pc12_ANGSD, pc12_GATK,
                           ncol = 3, labels = c("A", "B", "C"), legend = "top",
                           common.legend = TRUE, align = "hv")

# Save the plot with improved resolution
ggsave("PCA_Comparison.png", combined_plot, width = 18, height = 8, dpi = 300)