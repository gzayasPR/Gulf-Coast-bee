# Load necessary libraries
library(tidyverse)
library(stringr)
library(ggplot2)

# Clear workspace
rm(list=ls())

# Get command line arguments
args <- commandArgs(TRUE)

# Set project directory and paths
out.dir <- args[1]
meta.path <- args[2]
setwd(out.dir)

# Load metadata
meta_data <- read.csv(meta.path)
names(meta_data)[1] <- "BioSample"

# Recode locations
meta_data <- meta_data %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

# Function to process a single file and calculate the Heterozygosity
process_file <- function(file) {
  data <- scan(file, quiet = TRUE)
  Heterozygosity <- data[2] / sum(data) * 100  # Convert to percentage
  return(Heterozygosity)
}

# Get a list of all .est.ml files in the directory
file_list <- list.files(pattern = "\\.est\\.ml$")

# Apply the process_file function to each file and store the results
Heterozygosity <- sapply(file_list, process_file)

# Combine the results into a data frame
results <- data.frame(
  File = file_list,
  Heterozygosity = Heterozygosity
)
results$BioSample <- sub(".est.ml", "", results$File)

# Print the results
print(results)

# Save the results to a CSV file
write.csv(results, "Heterozygosity.csv", row.names = FALSE)

# Merge the results with metadata
merged_data <- merge(results, meta_data, by = "BioSample")

# Filter to include only females
female_data <- merged_data %>% filter(Sex == "female")

# Plot the ratios with boxplot separated by location
plot <- ggplot(female_data, aes(x = Location, y = Heterozygosity, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "Observed Heterozygosity by Location",
       x = "Location",
       y = "Heterozygosity (%)") +
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
ggsave("Heterozygosity_boxplot.png", plot, width = 11, height = 6)

# Optionally, save the filtered data to a CSV file
write.csv(female_data, "female_data.csv", row.names = FALSE)
