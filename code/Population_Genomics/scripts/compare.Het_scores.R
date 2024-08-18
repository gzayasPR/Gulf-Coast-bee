library(tidyverse)
library(ggpubr)


rm(list=ls())

args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
Het_Scores <- args[2]


setwd(out.dir)
Het_results <- read.csv(Het_Scores)

Het_results <- Het_results %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

# Define a color palette
color_palette <- c("AL" = "#4CAF50", "FL" = "#2196F3")
Het_results$OHET_GATK <- Het_results$OHET_GATK*100 
Het_results$OHET_DeepVariant <- Het_results$OHET_DeepVariant*100 
Het_results$OHET_ANGSD <- Het_results$OHET_ANGSD*100 
# Order Het_results based on the BioSample column
Het_results <- Het_results %>%
  arrange(Sample)
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
# Determine the y-axis limits
y_limits <- range(c(Het_results$OHET_ANGSD, Het_results$OHET_DeepVariant, Het_results$OHET_GATK), na.rm = TRUE)
y_limits[2] <- y_limits[2] +2
y_limits[1] <- y_limits[1] -2
Het_ANGSD <- ggplot(Het_results, aes(x = Location, y = OHET_ANGSD, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "ANGSD",
       x = "Location",
       y = "Observed Heterozygosity (%)")+
  scale_fill_manual(values = color_palette, name = "Location") +
  ylim(y_limits) +  # Set consistent y-axis limits
  base_theme

Het_DeepVariant <- ggplot(Het_results, aes(x = Location, y = OHET_DeepVariant, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "DeepVariant",
       x = "Location",
       y = "Observed Heterozygosity (%)")+
  scale_fill_manual(values = color_palette, name = "Location") +
  ylim(y_limits) +  # Set consistent y-axis limits
  base_theme

Het_GATK <- ggplot(Het_results, aes(x = Location, y = OHET_GATK, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "GATK",
       x = "Location",
       y = "Observed Heterozygosity (%)")+
  scale_fill_manual(values = color_palette, name = "Location") + 
  ylim(y_limits) +  # Set consistent y-axis limits
  base_theme

combined_plot <- ggarrange(Het_DeepVariant, Het_ANGSD, Het_GATK,
                           ncol = 3, labels = c("A", "B", "C"), legend = "top",
                           common.legend = TRUE, align = "hv")

# Save the plot with improved resolution
ggsave("Het_Comparison.png", combined_plot, width = 18, height = 8, dpi = 300)
