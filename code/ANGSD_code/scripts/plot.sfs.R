# Load necessary libraries
library(tidyverse)
library(stringr)
library(ape)

rm(list=ls())
args <- commandArgs(TRUE)

# Set project directory and paths
out.dir <- args[1]
setwd(out.dir)

sfs<-scan("smallFolded.sfs")
barplot(sfs[-1])
# Optionally, you can save the plot to a file
pdf("sfs.pdf")
barplot(sfs[-1])
dev.off()
