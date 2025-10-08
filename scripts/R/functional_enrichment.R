# Title: Functional Enrichment Analysis of Metabolite Sets
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This script performs a GSEA-like functional enrichment analysis
#              on sets of metabolites to identify enriched or depleted functions
#              in CRC vs. Healthy across different diets.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have this new library:
# install.packages("fgsea")

library(tidyverse)
library(fs)
library(fgsea) # For the enrichment analysis
library(rstatix)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_FUNCTIONAL")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. METADATA AND DATA PREPARATION
# ----------------------------------------------------
# --- Generate Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:133)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)

# --- Load Metabolite Data ---
all_scfa_files <- fs::dir_ls(path = diet_path, glob = "**/SCFA_Harmful_Metabolite_Profile*.csv", recurse = TRUE)
if (length(all_scfa_files) == 0) stop("CRITICAL ERROR: No SCFA/Metabolite files found.")
all_metabolite_data <- all_scfa_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    diet_name <- fs::path_split(.x)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(.x, show_col_types = FALSE, progress = FALSE) %>%
      select(starts_with("Metabolite_"), starts_with("WeightedFlux_")) %>%
      pivot_longer(everything(), names_to = c(".value", "id"), names_pattern = "(.+)_(\\d+)") %>%
      filter(!is.na(Metabolite) & !is.na(WeightedFlux)) %>%
      mutate(
        sample_id = sample_id,
        diet = diet_name,
        Metabolite = str_to_title(str_replace_all(Metabolite, "_", " "))
      ) %>%
      select(sample_id, diet, Metabolite, WeightedFlux)
  })
merged_flux_data <- all_metabolite_data %>%
  inner_join(metadata, by = "sample_id")
print("✅ Metabolite flux data loaded.")


# 4. DEFINE FUNCTIONAL METABOLITE SETS
# ----------------------------------------------------
functional_sets <- list(
  `Beneficial SCFAs` = c("Butyrate", "Propionate", "Acetate", "Valerate"),
  `Harmful Metabolites` = c("Hydrogen Sulfide", "Ammonia", "Indole", "Putrescine", "Trimethylamine"),
  `Other Organic Acids` = c("Lactate", "Succinate", "Formate")
)
print("✅ Functional metabolite sets defined.")


# 5. PREPARE RANKED LIST FOR GSEA (ROBUST METHOD)
# ----------------------------------------------------
# We will rank metabolites based on their log2 fold change. This is a robust way to avoid the previous errors.
ranked_metabolites <- merged_flux_data %>%
  group_by(diet, Metabolite) %>%
  # Calculate mean flux for each group
  summarise(
    mean_healthy = mean(WeightedFlux[group == "Healthy"], na.rm = TRUE),
    mean_crc = mean(WeightedFlux[group == "CRC"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Calculate the ranking metric (log2 fold change)
  mutate(rank_metric = log2((mean_crc + 1e-9) / (mean_healthy + 1e-9))) %>%
  # **CRITICAL FIX**: Remove any non-finite numbers that would crash the next step
  filter(is.finite(rank_metric)) %>%
  select(diet, Metabolite, rank_metric)
print("✅ Metabolite data ranked for enrichment analysis.")


# 6. RUN THE ENRICHMENT ANALYSIS
# ----------------------------------------------------
all_diets <- unique(ranked_metabolites$diet)
gsea_results <- all_diets %>%
  purrr::map_dfr(~{
    current_diet <- .x
    
    # Create the named vector required by fgsea
    ranks_for_diet <- ranked_metabolites %>%
      filter(diet == current_diet) %>%
      select(Metabolite, rank_metric) %>%
      deframe()
    
    # Run the GSEA
    fgsea(pathways = functional_sets, stats = ranks_for_diet, minSize=2) %>%
      as_tibble() %>%
      mutate(diet = current_diet)
  })

write_csv(gsea_results, fs::path(output_dir, "functional_enrichment_results.csv"))
print("✅ Functional enrichment analysis complete.")


# 7. VISUALIZE THE ENRICHMENT RESULTS
# ----------------------------------------------------
# A bar chart showing the Normalized Enrichment Score (NES).
# Positive NES: Enriched in CRC | Negative NES: Enriched in Healthy.
enrichment_plot <- ggplot(gsea_results, aes(x = diet, y = NES, fill = pathway)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Functional Enrichment of Metabolite Sets by Diet",
    subtitle = "Positive NES indicates enrichment in CRC; Negative NES indicates enrichment in Healthy",
    x = "Diet",
    y = "Normalized Enrichment Score (NES)",
    fill = "Metabolite Set"
  ) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = fs::path(output_dir, "Figure_Functional_Enrichment_Summary.png"),
       plot = enrichment_plot, width = 12, height = 8, dpi = 1200)
print("✅ Functional enrichment plot saved.")