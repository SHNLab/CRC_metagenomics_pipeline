# Title: Species-Level Differential Abundance Analysis
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This script identifies microbial species that are differentially
#              abundant between CRC and Healthy groups across different diets
#              and visualizes the results as a summary heatmap.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(pheatmap)
library(RColorBrewer)
library(rstatix)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "/Volumes/Backups of /bacarena/publication_style_analysis_FINAL_16spet")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. METADATA AND DATA PREPARATION
# ----------------------------------------------------
# --- Generate Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:132)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)
print("✅ Metadata generated.")

# --- Load and Consolidate Microbial Abundance Data ---
all_abundances_files <- fs::dir_ls(path = diet_path, glob = "**/updated_abundances*.csv", recurse = TRUE)
if (length(all_abundances_files) == 0) stop("CRITICAL ERROR: No abundance files found.")

all_abundances <- all_abundances_files %>%
  purrr::map_dfr(~ {
    file_path <- .x
    file_name <- fs::path_file(file_path)
    diet_name <- fs::path_split(file_path)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE) %>%
      mutate(
        diet = diet_name,
        sample_id = sample_id,
        Species = stringr::str_replace(Model, "\\.mat$", "") %>%
          stringr::str_replace_all("_", " ")
      )
  })

# --- Final Merged Data ---
merged_abundance_data <- all_abundances %>%
  inner_join(metadata, by = "sample_id") %>%
  filter(RelAbundance > 0)
print("✅ Microbial abundance data loaded and merged.")


# 4. DIFFERENTIAL ABUNDANCE ANALYSIS
# ----------------------------------------------------
# For each diet and species, compare abundance between CRC and Healthy
species_stat_results <- merged_abundance_data %>%
  group_by(diet, Species) %>%
  filter(n_distinct(group) == 2) %>%
  rstatix::wilcox_test(RelAbundance ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

readr::write_csv(species_stat_results, fs::path(output_dir, "species_differential_abundance_stats.csv"))
print("✅ Differential abundance statistics calculated.")


# 5. HEATMAP VISUALIZATION OF SIGNIFICANT SPECIES
# ----------------------------------------------------
# First, identify the most significant species to display on the heatmap.
# We will select species that are significant (p.adj < 0.05) in at least ONE diet comparison.
significant_species_list <- species_stat_results %>%
  filter(p.adj < 0.05) %>%
  pull(Species) %>%
  unique()

if(length(significant_species_list) > 0) {
  # Prepare the data matrix for the heatmap
  heatmap_data <- merged_abundance_data %>%
    filter(Species %in% significant_species_list) %>%
    group_by(diet, group, Species) %>%
    summarise(mean_abundance = mean(RelAbundance), .groups = "drop") %>%
    unite("diet_group", diet, group, sep = "_") %>%
    pivot_wider(names_from = diet_group, values_from = mean_abundance, values_fill = 0) %>%
    column_to_rownames("Species")
  
  # Annotation for the heatmap columns
  annotation_col <- data.frame(
    Condition = str_extract(colnames(heatmap_data), "(CRC|Healthy)"),
    Diet = str_remove(colnames(heatmap_data), "_(CRC|Healthy)")
  )
  rownames(annotation_col) <- colnames(heatmap_data)
  
  # Generate and save the heatmap
  png(fs::path(output_dir, "Figure5_Species_Heatmap.png"), width = 12, height = 14, units = "in", res = 1200)
  pheatmap(
    mat = heatmap_data,
    main = "Differentially Abundant Species (Z-score)",
    annotation_col = annotation_col,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    scale = "row", # Scale by species (Z-score)
    fontsize = 14,
    fontsize_row = 10,
    cutree_cols = 2
  )
  dev.off()
  print("✅ Figure 5 (Species Heatmap) saved.")
  
} else {
  print("No significantly different species found at the p.adj < 0.05 threshold.")
}

# --- Script End ---
print("--- 🎉 Species-level analysis complete! ---")