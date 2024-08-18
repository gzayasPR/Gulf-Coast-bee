library(tidyverse)
library(ggpubr)


rm(list=ls())

args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
Het_Scores <- args[2]


setwd(out.dir)
INB_results <- read.csv(Het_Scores)

INB_results <- INB_results %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

# Define a color palette
color_palette <- c("AL" = "#4CAF50", "FL" = "#2196F3")
INB_results$F__GATK <- as.numeric(INB_results$F_GATK)
INB_results$F__DeepVariant <-  as.numeric(INB_results$F_DeepVariant)
INB_results$F__ANGSD <-  as.numeric(INB_results$F_ANGSD)
print(INB_results)
# Order INB_results based on the BioSample column
INB_results <- INB_results %>%
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
y_limits <- range(c(INB_results$F_ANGSD, INB_results$F_DeepVariant, INB_results$F_GATK), na.rm = TRUE)
print(y_limits)
y_limits[2] <- y_limits[2] + 0.025
y_limits[1] <- 0
print(y_limits)
Het_ANGSD <- ggplot(INB_results, aes(x = Location, y = F_ANGSD, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "ANGSD",
       x = "Location",
       y = "Inbreeding Coefficient (F)")+
  scale_fill_manual(values = color_palette, name = "Location") +
coord_cartesian(ylim =y_limits,expand=TRUE) +  # Set consistent y-axis limits
  base_theme

Het_DeepVariant <- ggplot(INB_results, aes(x = Location, y = F_DeepVariant, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "DeepVariant",
       x = "Location",
       y = "Inbreeding Coefficient (F)")+
  scale_fill_manual(values = color_palette, name = "Location") +
coord_cartesian(ylim =y_limits,expand=TRUE) +  # Set consistent y-axis limits
  base_theme

Het_GATK <- ggplot(INB_results, aes(x = Location, y = F_GATK, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "GATK",
       x = "Location",
       y = "Inbreeding Coefficient (F)")+
  scale_fill_manual(values = color_palette, name = "Location") + 
coord_cartesian(ylim =y_limits,expand=TRUE) +  # Set consistent y-axis limits
  base_theme

combined_plot <- ggarrange(Het_DeepVariant, Het_ANGSD, Het_GATK,
                           ncol = 3, labels = c("A", "B", "C"), legend = "top",
                           common.legend = TRUE, align = "hv")

# Save the plot with improved resolution
ggsave("Inbreeding_Comparison.png", combined_plot, width = 18, height = 8, dpi = 300)
