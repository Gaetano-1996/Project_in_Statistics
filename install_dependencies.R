#==============================================================================#
############################ INSTALL DEPENDENCIES ############################
#==============================================================================#

# R script to install all required packages for the causal discovery analysis
# Run this script once before running the main analysis

cat("Installing required packages for causal discovery analysis...\n")

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(missing_packages) > 0) {
    cat("Installing:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
  } else {
    cat("All packages already installed.\n")
  }
}

# CRAN packages
cran_packages <- c(
  # Core data manipulation and visualization
  "tidyverse",      # Data manipulation and visualization
  "data.table",     # Fast data manipulation
  "haven",          # Read Stata files
  
  # Missing data analysis and visualization
  "mice",           # Multiple imputation
  "ggmice",         # MICE visualization
  "visdat",         # Missing data visualization
  "skimr",          # Summary statistics
  
  # Causal discovery
  "pcalg",          # PC algorithm for causal discovery
  "micd",           # Mixed-type data conditional independence
  
  # Parallel processing
  "future",         # Parallel computing framework
  "parallel"        # Base parallel processing
)

# Specialized packages (may require Bioconductor)
bioc_packages <- c(
  "tpc"             # Temporal PC algorithm
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install_if_missing(cran_packages)

# Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install specialized packages
cat("Installing specialized packages...\n")
for (pkg in bioc_packages) {
  if (!(pkg %in% installed.packages()[,"Package"])) {
    cat("Installing", pkg, "...\n")
    # Try CRAN first, then other sources
    tryCatch({
      install.packages(pkg)
    }, error = function(e) {
      cat("Package", pkg, "not found on CRAN. Please install manually.\n")
      cat("For tpc package, you may need to install from source or GitHub.\n")
    })
  }
}

cat("Package installation completed!\n")
cat("You can now run the main analysis with: source('main.R')\n")