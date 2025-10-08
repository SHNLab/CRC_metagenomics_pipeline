###############################################################################
# Script Name:  filter_and_sort_abundance.R
# Author     :  Aman Sharma
# Date       :  2025-10-07
# Institute  :  National Agri-Food Biomanufacturing Institute (NABI), Mohali, India
# Description:  This script reads previously aggregated microbial abundance XLSX 
#               files (healthy and diseased samples), filters out low-abundance
#               models (TotalAbundance ≤ 0.001), sorts them in descending order 
#               of abundance, and saves the processed files into a separate folder 
#               ("filtered_files"). Intended for downstream analysis and visualization.
###############################################################################

# Load required libraries
library(readxl)
library(writexl)

# ------------------------ USER PARAMETERS ------------------------
# Threshold for filtering low-abundance entries
abundance_threshold <- 0.001

# Output folder for filtered files
out_dir <- "filtered_files"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Pattern to select input files (healthy/diabetic CRC samples)
file_pattern <- "^[dh]i_.*\\.xlsx$"

# ------------------------ STEP 1: List All Input Files ------------------------
files <- list.files(pattern = file_pattern)

# ------------------------ STEP 2: Process Each File ------------------------
for (file in files) {
  
  # Read the Excel file
  df <- read_excel(file)
  
  # Filter for TotalAbundance above threshold and sort descending
  df_filtered <- df %>%
    dplyr::filter(TotalAbundance > abundance_threshold) %>%
    dplyr::arrange(dplyr::desc(TotalAbundance))
  
  # Construct output file path
  out_file <- file.path(out_dir, file)
  
  # Save filtered and sorted data
  write_xlsx(df_filtered, out_file)
  
  cat("Processed:", file, "→", out_file, "\n")
}

cat("All files processed successfully.\n")
