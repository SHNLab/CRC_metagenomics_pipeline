###############################################################################
# Script Name:  aggregate_model_abundance.R
# Author     :  Aman Sharma
# Date       :  2025-10-07
# Institute  :  National Agri-Food Biomanufacturing Institute (NABI), Mohali, India
# Description:  This script processes microbial abundance CSV files from a
#               "Filtered" folder. For each sample, it aggregates abundances 
#               by metabolic model, matches the sample to its group
#               ("healthy" or "diseased") based on a metadata CSV, and exports 
#               the results as Excel files with standardized naming.
#               Intended for downstream analysis of metagenomic datasets.
###############################################################################

# Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(writexl)

# ------------------------ USER PARAMETERS ------------------------
# Path to sample metadata file
sample_metadata_file <- "sample_ids.csv"   # Columns: sample_id, group

# Folder containing filtered CSV abundance files
input_folder <- "Filtered"

# Output folder (optional, can be same as working directory)
output_folder <- "Processed"

# Create output folder if it does not exist
if (!dir.exists(output_folder)) dir.create(output_folder)

# ------------------------ STEP 1: Read Sample Metadata ------------------------
sample_info <- read_csv(sample_metadata_file) %>%
  mutate(number = row_number())   # Generate numeric IDs for naming

# ------------------------ STEP 2: List All Input Files ------------------------
files <- list.files(path = input_folder, pattern = "*.csv", full.names = TRUE)

# ------------------------ STEP 3: Process Each File ------------------------
for (file in files) {
  
  # Extract base sample name (e.g., SRR5898924)
  base_name <- str_split(basename(file), "_")[[1]][1]
  
  # Read CSV and select relevant columns (2nd and 3rd)
  df <- read_csv(file, col_types = cols()) %>%
    select(2, 3)
  
  colnames(df) <- c("Abundance", "ModelFile")
  
  # Aggregate abundances by ModelFile
  df_agg <- df %>%
    group_by(ModelFile) %>%
    summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
    filter(!is.na(ModelFile) & ModelFile != "")   # Remove empty ModelFile rows
  
  # Match with sample metadata
  sample_row <- sample_info %>% filter(sample_id == base_name)
  
  if (nrow(sample_row) == 0) {
    warning(paste("No metadata match found for sample:", base_name))
    next
  }
  
  # Extract group (healthy/diseased) and numeric ID
  status <- sample_row$group[1]
  num <- sample_row$number[1]
  
  # Define output filename based on group
  out_file <- switch(
    tolower(status),
    "healthy" = paste0("hi_p", num, ".xlsx"),
    "crc"     = paste0("di_p", num, ".xlsx"),
    {
      warning(paste("Unknown group for sample:", base_name))
      next
    }
  )
  
  # Prepend output folder path
  out_file_path <- file.path(output_folder, out_file)
  
  # Save aggregated data to XLSX
  write_xlsx(df_agg, out_file_path)
  
  cat("Saved:", out_file_path, "\n")
}

cat("All files processed successfully.\n")
