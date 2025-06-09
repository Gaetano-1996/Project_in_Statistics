#==============================================================================#
###################### WHY DO PEOPLE HAVE CHILDREN ? ###########################
#==============================================================================#
# 
# This script consolidates the entire causal discovery analysis pipeline:
# 1. Data preparation from raw GGS files
# 2. Exploratory data analysis 
# 3. Causal discovery using temporal PC algorithm with multiple imputation
#
# All functions are sourced from src.R
# Results are made available in global environment and saved on disk.


# SETUP ====================================================================####

# Clear environment
rm(list = ls())

# Load functions
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("src.R")

# Define directory structure
plots_dir = "plots2"
eda_dir = "plots2/eda"
chains_dir = "plots2/chains/"
graphs_dir = "plots2/graphs/"
summaries_dir = "summaries2"


# Create directory structure
create_directories()

# DATA PREPARATION =========================================================####

# Prepare GGS data 
# Returns df_final and saves df_mix.rds, df_final.rds
data_dir = "~/Desktop/Research/Causal Discovery/Data csv"
save_data = "~/Desktop/Research/Causal Discovery/data_clean"
df_final = prepare_ggs_data(data_dir, 
                              save_data)

# EXPLORATORY DATA ANALYSIS ================================================####

# Perform overall EDA
eda_results = perform_eda(df_final, eda_dir = eda_dir, summaries_dir = summaries_dir)

# Generate country-specific summaries
country_stats = generate_country_summaries(df_final, eda_dir = eda_dir, summaries_dir = summaries_dir)

# TIER DEFINITION ==========================================================####


# Get variable names (excluding country for tier assignment)
df_names = names(df_final %>% select(-country))

# Create tiers vector
tiers_numeric = create_tiers_vector(df_names)

cat("Tier structure created:\n")
cat("- T0 (Background):", sum(tiers_numeric == 0), "variables\n")
cat("- T1 (Demographics):", sum(tiers_numeric == 1), "variables\n") 
cat("- T2 (Past events):", sum(tiers_numeric == 2), "variables\n")
cat("- T3 (Current status):", sum(tiers_numeric == 3), "variables\n")
cat("- T4 (Intentions):", sum(tiers_numeric == 4), "variables\n")
cat("- T5 (Outcome):", sum(tiers_numeric == 5), "variables\n")


# CAUSAL DISCOVERY ==========================================================####


# Run causal discovery with multiple imputation and temporal PC
discovery_results = ggs_discovery(
  data = df_final,
  tiers = tiers_numeric,
  alpha = 0.1,
  m = 2,  # Number of imputations (reduced for faster execution)
  maxit = 2, # Maximum number of iteration for MICE (reduced for faster execution)
  n.core = parallel::detectCores(),
  plot.graph = TRUE,
  save.graph = F,
  plot.chains = TRUE,  # Set to TRUE if you want to see convergence plots
  save.chains = F,
  neighborhood = c(1, 2),  # Save both 1-hop and 2-hop neighborhood graphs
  mincor = 0.1,
  chains_path = chains_dir,
  graphs_path = graphs_dir,
  verbose = 1
)
