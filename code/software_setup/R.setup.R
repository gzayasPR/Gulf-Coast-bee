# Set the library path to use the user-specific R library directory
user_lib <- Sys.getenv("R_LIBS_USER")
.libPaths(user_lib)

# Check if the directory is writable
if (!dir.exists(user_lib) || !file.access(user_lib, 2) == 0) {
  stop("Library path is not writable: ", user_lib)
}

# Install necessary packages if not already installed
packages <- c("vcfR","genetics", "parallel", "bedr", "gaston", "data.table", 
              "adegenet", "ggplot2", "ape", "ggpubr", "detectRUNS", 
              "tidyverse", "pheatmap")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org", lib = user_lib)
  }
}

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.r-project.org", lib = user_lib)
}