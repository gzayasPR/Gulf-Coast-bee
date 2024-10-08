if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR",repos="http://cran.r-project.org")
if (!requireNamespace("genetics", quietly = TRUE)) install.packages("genetics",repos="http://cran.r-project.org")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel",repos="http://cran.r-project.org")
if (!requireNamespace("bedr", quietly = TRUE)) install.packages("bedr",repos="http://cran.r-project.org")
if (!requireNamespace("gaston", quietly = TRUE)) install.packages("gaston",repos="http://cran.r-project.org")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table",repos="http://cran.r-project.org")
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet",repos="http://cran.r-project.org")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2",repos="http://cran.r-project.org")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape",repos="http://cran.r-project.org")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr",repos="http://cran.r-project.org")
if (!requireNamespace("detectRUNS", quietly = TRUE)) install.packages("detectRUNS",repos="http://cran.r-project.org")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse",repos="http://cran.r-project.org")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap",repos="http://cran.r-project.org")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos="http://cran.r-project.org")
library(BiocManager)
BiocManager::install("SNPRelate")