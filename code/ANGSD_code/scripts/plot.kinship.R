library(ggplot2)
library(dplyr)
library(tidyverse)
library(pheatmap)

# Clear environment
rm(list=ls())

# Set project directory and paths from command line arguments
args <- commandArgs(TRUE)
out.dir <- args[1]
meta.path <- args[2]
ids.path <- args[3]
setwd(out.dir)
print(args)

# Load the kinship data
kinship_data <- read.table("newres", header=TRUE)

# Load the metadata
meta_data <- read.csv(meta.path)
names(meta_data)[1] <- "BioSample"
meta_data <- meta_data %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

# Load the IDs
ids <- read.table(ids.path)[1]
names(ids) <- "BioSample"
ids$index <- 0:(nrow(ids) - 1)

# Merge the kinship data with the IDs and metadata
kinship_data <- kinship_data %>%
  left_join(ids, by = c("a" = "index")) %>%
  rename(BioSampleA = BioSample) %>%
  left_join(ids, by = c("b" = "index")) %>%
  rename(BioSampleB = BioSample)

kinship_data <- kinship_data %>%
  left_join(meta_data, by = c("BioSampleA" = "BioSample")) %>%
  rename(LocationA = Location, SexA = Sex) %>%
  left_join(meta_data, by = c("BioSampleB" = "BioSample")) %>%
  rename(LocationB = Location, SexB = Sex)

# Create a custom factor for location comparisons
kinship_data$Location_Comparison <- paste(kinship_data$LocationA, kinship_data$LocationB, sep="x")
kinship_data$Location_Comparison[kinship_data$Location_Comparison == "ALxFL"] <- "FLxAL"
kinship_data$Location_Comparison <- factor(kinship_data$Location_Comparison, levels= c("FLxFL","ALxAL", "FLxAL"))

# Order BioSampleA and BioSampleB based on Location_Comparison
ordered_biosamples <- kinship_data %>%
  arrange(LocationA,LocationB) %>%
  select(BioSampleA, BioSampleB) %>%
  unique() %>%
  unlist() %>%
  unique()

kinship_data$BioSampleA <- factor(kinship_data$BioSampleA, levels = ordered_biosamples)
kinship_data$BioSampleB <- factor(kinship_data$BioSampleB, levels = ordered_biosamples)

# Create a wide format matrix for the kinship data
kinship_matrix <- kinship_data %>%
  select(BioSampleA, BioSampleB, KING) %>%
  spread(key = BioSampleB, value = KING, fill = 0) %>%
  column_to_rownames(var = "BioSampleA")

breaks <- c(seq(min(kinship_matrix), 0, length.out=51), seq(0, max(kinship_matrix), length.out=50)[-1])

# Create a heatmap with the midpoint at 0
plot <- pheatmap(as.matrix(kinship_matrix),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         main = "Kinship Coefficients Heatmap",
         cluster_rows = FALSE,
         cluster_cols = FALSE)
# Save the plot
ggsave("Kinship_heatmap_ordered.png",plot , width = 11, height = 6)

# Set the correct column order
col_order <- c("BioSampleA", "BioSampleB", "a", "b", "LocationA", "LocationB", "SexA", "SexB", "nSites", "rab", "Fa", "Fb", "theta", "inbred_relatedness_1_2", "inbred_relatedness_2_1", "fraternity", "identity", "zygosity", "2of3_IDB", "FDiff", "loglh", "nIter", "bestoptimll", "coverage", "2dsfs", "R0", "R1", "KING", "2dsfs_loglike", "2dsfsf_niter")
str(kinship_data)
# Reorder columns
kinship_data.save <- kinship_data[col_order]

# Save the reordered data to a CSV file
write.csv(kinship_data.save, "kinship_ordered.csv")
