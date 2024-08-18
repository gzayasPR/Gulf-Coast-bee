library(tidyverse)
rm(list=ls())

args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
het.path<- args[2]
meta.path <- args[3]


setwd(out.dir)
# Read heterozygosity data
het_data <- read.table(het.path, header=TRUE)
names(het_data )[1] <- "BioSample"
head(het_data)
het_data$Observed_Heterozygosity <-  ((het_data$N_SITES - het_data$O.HOM.) / het_data$N_SITES)*100
het_data$Expected_Heterozygosity <- ((het_data$N_SITES - het_data$E.HOM.)  / het_data$N_SITES)*100
# Read metadata
metadata <- read.csv(meta.path)
names(metadata )[1] <- "BioSample"


metadata <- metadata %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))
# Merge heterozygosity data with metadata
merged_data <- merge(het_data, metadata, by='BioSample')
merged_data <- merged_data %>% filter(Sex == "female")
# Plot observed vs expected heterozygosity

# Plot the ratios with boxplot separated by location
plot <- ggplot(merged_data, aes(x = Location, y = Observed_Heterozygosity, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "Observed Heterozygosity by Location",
       x = "Location",
       y = "Observed Heterozygosity (%)") +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("AL" = "#4CAF50", "FL" = "#2196F3")) +
  theme(panel.background = element_rect(fill = "white", color = "white"))

# Save the plot
ggsave("Observed Heterozygosity_boxplot.png", plot, width = 11, height = 6)



# Plot the ratios with boxplot separated by location
plot <- ggplot(merged_data, aes(x = Location, y = Expected_Heterozygosity, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "Expected Heterozygosity by Location",
       x = "Location",
       y = "Expected Heterozygosity (%)") +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("AL" = "#4CAF50", "FL" = "#2196F3")) +
  theme(panel.background = element_rect(fill = "white", color = "white"))

# Save the plot
ggsave("Expected Heterozygosity_boxplot.png", plot, width = 11, height = 6)

write.csv(merged_data,"Heterozygosity.csv")