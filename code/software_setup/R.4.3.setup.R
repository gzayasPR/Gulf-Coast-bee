# Create the local library if it doesn't exist
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)

# Install necessary packages if not already installed
packages <- c("tidyverse", "ggplot2")


install.packages(packages, lib = Sys.getenv("R_LIBS_USER"), repos = "http://cran.us.r-project.org", 
                 configure.args="--with-udunits2-lib=$UDUNITS_ROOT/lib --with-udunits2-include=$UDUNITS_ROOT/include")

install.packages(packages,  repos = "http://cran.us.r-project.org")