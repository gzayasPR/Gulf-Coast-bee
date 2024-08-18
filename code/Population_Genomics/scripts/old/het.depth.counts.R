# Install and load necessary packages if not already installed
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR")
if (!requireNamespace("genetics", quietly = TRUE)) install.packages("genetics")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("bedr", quietly = TRUE)) install.packages("bedr")
if (!requireNamespace("gaston", quietly = TRUE)) install.packages("gaston")
library(vcfR)
library(genetics)
library(parallel)
library(bedr)
library(gaston)
library(tidyverse)
rm(list=ls())
# Set project directory and paths
proj.dir <- "C:/Users/gzayas97/Documents/USDA Internship/"
setwd(proj.dir)
data.dir <- paste(proj.dir, "/data/", sep="")
results.dir <- paste(proj.dir, "/results/", sep="")
vcf.path <- paste(data.dir, "/Hesperapis_oraria.vcf", sep="")
vcf.path <- paste(data.dir, "/Hesperapis_oraria.vcf", sep="")
vcf <- read.vcfR(vcf.path)

# Extract the genotype (GT) and depth (DP) information
meta.path <- paste(data.dir,"/Samples.Metadata.csv",sep="")
meta.data <- read.csv(meta.path)
females.id <- meta.data$BioSample[meta.data$Sex == "female"]

# Get the column indices of the individuals
individual_indices <- which(colnames(vcf@gt) %in% females.id)

# Subset the VCF object
females.vcf <- vcf[, c(1, individual_indices)]

# Extract the genotype matrix
genotype_matrix <- extract.gt(females.vcf)
depth_matrix <- extract.gt(females.vcf, element = "DP", as.numeric = TRUE)
str(genotype_matrix)
# Filter genotypes based on minimum depth
genotype_matrix[depth_matrix > 50 ] <- NA
str(genotype_matrix)

dp.vector <- c(depth_matrix)
gt.vector <- c(genotype_matrix)


# Combine read depth and heterozygosity into a data frame
data <- data.frame(ReadDepth = dp.vector, Genotype = gt.vector)
# Create a logical vector to identify heterozygous calls
data$heterozygous <- data$Genotype %in% c("0/1","0|1","1/0","1|0")
data$heterozygous <- ifelse(is.na(data$Genotype), NA,data$Genotype %in% c("0/1","0|1","1/0","1|0"))
summary_data <- data %>%
  group_by(ReadDepth) %>%
  summarise(ProportionHet = mean(heterozygous, na.rm = TRUE),
            Count = n(),
            Count.het = sum(heterozygous, na.rm = TRUE))

na.omit(summary_data)
combined_plot <- ggplot(na.omit(summary_data)) +
  geom_bar(aes(x = ReadDepth, y = Count/1000000), stat = "identity", fill = "lightblue", alpha = 0.7) +
  geom_point(aes(x = ReadDepth, y = ProportionHet), size = 3, color = "red") +
  geom_line(aes(x = ReadDepth, y = ProportionHet), color = "red") +
  scale_x_continuous(
    breaks = seq(5, 50, by = 5)
  ) +
  scale_y_continuous(
    name = "Number of genotypes (millions)",
    sec.axis = sec_axis(~., name = "Proportion of Heterozygous genotypes")
  ) +
  labs(x = "Read depth") +
  theme_minimal()+
  theme(plot.background = element_rect(fill='white'))


# Print the combined plot
print(combined_plot)

setwd(results.dir)
dir.create("DepthvsHet",showWarnings = F)
setwd("DepthvsHet")

ggsave("Read Depth vs Genotype Counts and Het. Proportions.png",combined_plot)



combined_plot <- ggplot(na.omit(summary_data)) +
  geom_bar(aes(x = ReadDepth, y = Count.het /100000), stat = "identity", fill = "lightblue", alpha = 0.7) +
  geom_point(aes(x = ReadDepth, y = ProportionHet), size = 3, color = "red") +
  geom_line(aes(x = ReadDepth, y = ProportionHet), color = "red") +
  scale_x_continuous(
    breaks = seq(5, 50, by = 5)
  ) +
  scale_y_continuous(
    name = "Number of het genotypes (Hundred Thousands)",
    sec.axis = sec_axis(~., name = "Proportion Heterozygous genotypes")
  ) +
  labs(x = "Read depth") +
  theme_minimal()+
  theme(plot.background = element_rect(fill='white'))


# Print the combined plot
print(combined_plot)

