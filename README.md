# Why do people have children? 
## Insights from a cross-national study
### Project in Statistics 24/25, Gaetano Tedesco 

This repository contains a comprehensive causal discovery analysis investigating fertility decisions across developed countries using the Generations and Gender Survey (GGS) data. The project employs multiple imputation combined with temporal PC algorithms on mixed data to identify causal relationships influencing the decision to have children.

## 🎯 Research Question

**What are the causal factors that influence people's decisions to have children across different developed countries?**

This study uses causal discovery methods to move beyond correlational analyses and identify plausible causal relationships in fertility decision-making, accounting for temporal ordering and cultural differences across countries.

## 📊 Data Availability

**Important Note**: The GGS microdata used in this analysis are **not publicly available** and are not included in this repository due to data protection regulations. 

To replicate this analysis, you would need to:

1. Apply for data access through the [GGS Data Archive](https://www.ggp-i.org/)
2. Obtain the following files:
   - Wave 1 survey data
   - Wave 2 survey data  
   - (Stata file with variable labels)
3. Place these files in your designated data directory and update the `data_dir` path in `main.R`

## 📁 Repository Structure

```
Project_in_Statistics/
├── README.md                   # This file
├── install_dependencies.R      # Package installation script
├── main.R                      # Main analysis pipeline
├── src.R                       # Function library with full documentation
├── no_context                  # Experiment without age and sex
├── figures/                    # Analysis output visualizations
│   ├── eda/                    # Exploratory data analysis plots
│   ├── chains/                 # MICE convergence diagnostics
│   └── graphs/                 # Causal discovery results
└── summaries/                  # Statistical summary files
    ├── overall_summary.csv     # Dataset-wide statistics
    └── summary_{country}.csv   # Country-specific summaries
```

## 🔧 Installation and Setup

### Prerequisites
- R version 4.5.0 or higher
- RStudio (recommended)
- Sufficient computational resources (analysis is computationally intensive)

### Quick Setup

1. **Clone this repository:**
   ```bash
   git clone <repository-url>
   cd Project_in_Statistics
   ```

2. **Install dependencies:**
   ```r
   source("install_dependencies.R")
   ```

3. **Verify installation:**
   ```r
   source("src.R")  # Should load without errors
   ```

## 🚀 Running the Analysis

### Full Analysis Pipeline

1. **Set up your data directory** in `main.R`:
   ```r
   data_dir = "path/to/your/ggs/data"
   save_data = "path/to/processed/data"
   ```

2. **Run the complete analysis:**
   ```r
   source("main.R")
   ```

### Analysis Steps

The `main.R` script orchestrates the following pipeline:

1. **Data Preparation** (`prepare_ggs_data`)
   - Imports and harmonizes GGS Wave 1 and Wave 2 data
   - Creates outcome variable (new child between waves)
   - Processes questionnaire items across multiple domains and generate variables
   - Filters observations with >25% missing values

2. **Exploratory Data Analysis** (`perform_eda`, `generate_country_summaries`)
   - Generates comprehensive summary statistics
   - Creates missing data visualizations
   - Produces country-specific comparisons

3. **Temporal Tier Definition** (`create_tiers_vector`)
   - Assigns variables to temporal tiers (T0-T5)
   - Ensures proper causal ordering for discovery

4. **Causal Discovery** (`ggs_discovery`)
   - Applies multiple imputation (MICE) with 25 imputed datasets
   - Runs temporal PC algorithm by country (mixed data CI test from `micd`)
   - Generates causal graphs and convergence diagnostics

### Customizing the Analysis

Key parameters in `main.R`:
```r
# Parallel processing
n.core = parallel::detectCores() 

# Multiple imputation settings
m = 25               # Number of imputed datasets
maxit = 25           # MICE iterations
parallelseed = 123   # random seed for parallel computing

# Causal discovery settings  
alpha = 0.1                   # Significance level
context.var = NULL            # allow defintion of context variables
tiers = tiers                 # tiers specification
indepTest = micd::mixMItest   # Conditional independence test (mix data + MI) 
 
```

## 📈 Output Files

### Visualizations
- **EDA plots**: Missing data patterns, outcome distributions
- **Convergence plots**: MICE diagnostic chains by country
- **Causal graphs**: Full networks and neighborhood subgraphs

### Summary Files
- **CSV summaries**: Comprehensive statistics by country
- **RDS datasets**: Processed data for further analysis

## 🔬 Methodology

### Temporal Structure
Variables are organized into temporal tiers:

- **T0**: Background (sex, age)
- **T1**: Demographics (education, migration)  
- **T2**: Past events (previous births, deceased children)
- **T3**: Current status (n. of kids, health, relationships, housing)
- **T4**: Intentions (fertility intentions)
- **T5**: Outcome (new child)

### Statistical Approach
- **Multiple Imputation**: Random forest imputation for missing data
- **Causal Discovery**: Temporal PC algorithm with conditional Gaussian independence testing
- **Country-Specific Analysis**: Separate models to capture cultural differences

## 💻 System Requirements

### Computational Considerations
- **Memory**: Analysis requires substantial RAM (≥8GB recommended)
- **Processing**: Multi-core CPU beneficial for parallel processing
- **Time**: Full analysis may take 1-4 hours depending on hardware

### Dependencies
See `install_dependencies.R` for complete package list. 

Key packages:

- `tidyverse` - Data manipulation and visualization
- `mice` / `ggmice` - Multiple imputation
- `pcalg` / `tpc` / `micd` - Causal discovery algorithms
- `future` / `parallel` - Parallel processing

---

## 📋 R Session Information

The analysis was performed on a Mac Book Pro 14', M3 MAX (14c CPU, 40c GPU), 36 Gb of RAM with the following R environment:

```
R version 4.5.0 (2025-04-11)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.5

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C/UTF-8/C/C/C/C

time zone: Europe/Copenhagen
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] skimr_2.1.5       haven_2.5.5       data.table_1.17.4 future_1.49.0    
 [5] visdat_0.6.0      tpc_1.0           micd_1.1.1        ggmice_0.1.0     
 [9] mice_3.18.0       pcalg_2.7-12      lubridate_1.9.4   forcats_1.0.0    
[13] stringr_1.5.1     dplyr_1.1.4       purrr_1.0.4       readr_2.1.5      
[17] tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.2     tidyverse_2.0.0  

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1    farver_2.1.2        fastmap_1.2.0      
 [4] fastICA_1.2-7       digest_0.6.37       rpart_4.1.24       
 [7] timechange_0.3.0    lifecycle_1.0.4     cluster_2.1.8.1    
[10] survival_3.8-3      magrittr_2.0.3      compiler_4.5.0     
[13] rlang_1.1.6         tools_4.5.0         igraph_2.1.4       
[16] knitr_1.50          repr_1.1.7          RColorBrewer_1.1-3 
[19] abind_1.4-8         withr_3.0.2         BiocGenerics_0.54.0
[22] nnet_7.3-20         grid_4.5.0          stats4_4.5.0       
[25] jomo_2.7-6          globals_0.18.0      scales_1.4.0       
[28] iterators_1.0.14    MASS_7.3-65         cli_3.6.5          
[31] reformulas_0.4.1    generics_0.1.4      RcppParallel_5.1.10
[34] robustbase_0.99-4-1 tzdb_0.5.0          bdsmatrix_1.3-7    
[37] minqa_1.2.8         splines_4.5.0       BiocManager_1.30.25
[40] base64enc_0.1-3     vctrs_0.6.5         boot_1.3-31        
[43] glmnet_4.1-9        Matrix_1.7-3        jsonlite_2.0.0     
[46] hms_1.1.3           RBGL_1.84.0         mitml_0.4-5        
[49] clue_0.3-66         listenv_0.9.1       foreach_1.5.2      
[52] parallelly_1.45.0   glue_1.8.0          ggm_2.5.1          
[55] nloptr_2.2.1        DEoptimR_1.1-3-1    pan_1.9            
[58] codetools_0.2-20    stringi_1.8.7       shape_1.4.6.1      
[61] gtable_0.3.6        sfsmisc_1.1-20      lme4_1.1-37        
[64] pillar_1.10.2       htmltools_0.5.8.1   graph_1.86.0       
[67] R6_2.6.1            zigg_0.0.2          Rdpack_2.6.4       
[70] evaluate_1.0.3      lattice_0.22-7      rbibutils_2.3      
[73] backports_1.5.0     Rfast_2.1.5.1       broom_1.0.8        
[76] corpcor_1.6.10      Rcpp_1.0.14         nlme_3.1-168       
[79] xfun_0.52           pkgconfig_2.0.3
```