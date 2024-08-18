# Load necessary libraries
library(tidyverse)
library(stringr)
library(ape)

rm(list=ls())
args <- commandArgs(TRUE)

# Set project directory and paths
out.dir <- args[1]
name <- args[2]
meta.path <- args[3]
ids.path <- args[4]
setwd(out.dir)

# Load metadata
meta_data <- read.csv(meta.path)
names(meta_data)[1] <- "BioSample"
meta_data$color <- meta_data$Location
meta_data <- meta_data %>%
  mutate(Location = recode(Location,
                           "USA-AL_Baldwin_co" = "AL",
                           "USA-FL_Escambia_co" = "FL"))

ids <- read.table(ids.path)[1]
names(ids) <- "BioSample"
ids$index <- 1:nrow(ids)

# Merge metadata with ids to get location information
ids <- merge(ids, meta_data, by="BioSample")

# Set the path to the tree file
tree_file <- "Tree.tree"

# Read the tree file
tree <- read.tree(tree_file)

# Replace tip labels with BioSample
tree$tip.label <- ids$BioSample[match(tree$tip.label, ids$index)]

# Create a vector of colors based on the location
location_colors <- ids %>% 
  distinct(Location) %>%
  mutate(color = rainbow(n()))

# Merge colors with ids
ids <- merge(ids, location_colors, by="Location")
# Optionally, you can save the plot to a file
pdf("tree_plot.pdf")
plot.phylo(tree, main="Phylogenetic Tree from ANGSD", tip.color=ids$color.y[match(tree$tip.label, ids$BioSample)])
dev.off()
