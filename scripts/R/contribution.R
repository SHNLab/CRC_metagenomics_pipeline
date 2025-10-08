# Title: Definitive In-Depth Functional Contribution Analysis
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-17
# Description: This definitive script calculates the functional contribution of
#              species to the production of ALL metabolites, visualizes the
#              results as publication-quality bubble plots with a dark pastel
#              color scheme, and generates a comprehensive summary text file.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have these packages installed:
# install.packages(c("tidyverse", "fs", "rstatix", "RColorBrewer"))

library(tidyverse)
library(fs)
library(rstatix)
library(RColorBrewer) # For advanced color scales


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
# Create a new, dedicated folder for these enhanced results
output_dir <- fs::path(base_path, "publication_style_analysis_FUNCTIONAL_V3")
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
print("✅ Metadata generated.")

# --- Load Abundance Data ---
all_abundances_files <- fs::dir_ls(path = diet_path, glob = "**/updated_abundances*.csv", recurse = TRUE)
if (length(all_abundances_files) == 0) stop("CRITICAL ERROR: No abundance files found.")
all_abundances <- all_abundances_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    diet_name <- fs::path_split(.x)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(.x, show_col_types = FALSE, progress = FALSE) %>%
      mutate(
        diet = diet_name,
        sample_id = sample_id,
        Species = stringr::str_replace(Model, "\\.mat$", "") %>% stringr::str_replace_all("_", " ")
      )
  })
merged_abundance_data <- all_abundances %>%
  inner_join(metadata, by = "sample_id") %>%
  mutate(sample_diet_id = paste(sample_id, diet, sep = "___"))
print("✅ Microbial abundance data loaded.")

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
  inner_join(metadata, by = "sample_id") %>%
  mutate(sample_diet_id = paste(sample_id, diet, sep = "___"))
print("✅ Metabolite flux data loaded.")


# 4. CALCULATE AND VISUALIZE FUNCTIONAL CONTRIBUTIONS FOR ALL METABOLITES
# ----------------------------------------------------
# --- Initialize the summary report ---
report_lines <- list()
report_lines[[1]] <- "--- Functional Contribution Analysis Summary ---"

# --- Get the list of all metabolites with positive flux to analyze ---
metabolites_to_analyze <- merged_flux_data %>%
  filter(WeightedFlux > 0) %>%
  pull(Metabolite) %>%
  unique()

# --- Loop through each metabolite ---
for (current_metabolite in metabolites_to_analyze) {
  
  cat(paste("\n--- Processing:", current_metabolite, "---\n"))
  
  # Filter for the current metabolite
  metabolite_flux <- merged_flux_data %>%
    filter(Metabolite == current_metabolite) %>%
    select(sample_diet_id, total_metabolite_flux = WeightedFlux)
  
  # Calculate Contribution Score = Species Abundance * Total Metabolite Flux
  contribution_data <- merged_abundance_data %>%
    inner_join(metabolite_flux, by = "sample_diet_id") %>%
    mutate(contribution_score = RelAbundance * total_metabolite_flux)
  
  # --- Find Top 10 Overall Contributors for plotting ---
  top_contributors <- contribution_data %>%
    group_by(Species) %>%
    summarise(mean_contribution = mean(contribution_score, na.rm = TRUE)) %>%
    slice_max(order_by = mean_contribution, n = 10) %>%
    pull(Species)
  
  # --- Summarize Data for plotting and reporting ---
  summary_data <- contribution_data %>%
    filter(Species %in% top_contributors) %>%
    group_by(diet, group, Species) %>%
    summarise(
      mean_contribution = mean(contribution_score, na.rm = TRUE),
      mean_abundance = mean(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- Add Results to the Summary Report ---
  report_lines[[length(report_lines) + 1]] <- paste0("\n## Metabolite: ", current_metabolite)
  top_healthy <- summary_data %>% filter(group == "Healthy") %>% arrange(desc(mean_contribution)) %>% head(3)
  top_crc <- summary_data %>% filter(group == "CRC") %>% arrange(desc(mean_contribution)) %>% head(3)
  
  report_lines[[length(report_lines) + 1]] <- "  **Top 3 Contributors in HEALTHY group:**"
  for(i in 1:nrow(top_healthy)) {
    report_lines[[length(report_lines) + 1]] <- paste0("    ", i, ". ", top_healthy$Species[i], " (Avg. Contribution: ", round(top_healthy$mean_contribution[i], 2), ")")
  }
  report_lines[[length(report_lines) + 1]] <- "  **Top 3 Contributors in CRC group:**"
  for(i in 1:nrow(top_crc)) {
    report_lines[[length(report_lines) + 1]] <- paste0("    ", i, ". ", top_crc$Species[i], " (Avg. Contribution: ", round(top_crc$mean_contribution[i], 2), ")")
  }
  
  
  # --- Generate the Bubble Plot ---
  bubble_plot <- ggplot(summary_data, aes(x = diet, y = Species, size = mean_abundance, color = mean_contribution)) +
    # **FIX**: Use shape 16 for solid circles without an outline
    geom_point(shape = 16, alpha = 1) +
    facet_wrap(~ group) +
    # **FIX**: Use a dark pastel color palette from RColorBrewer
    scale_color_distiller(palette = "Spectral", direction = 1, name = "Mean Contribution") +
    scale_size(range = c(2, 15), name = "Mean Relative\nAbundance") +
    labs(
      title = paste("Functional Contribution of Top Species to", current_metabolite, "Production"),
      subtitle = "Comparing CRC vs. Healthy across different diets",
      x = "Diet",
      y = "Bacterial Species"
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 16),
      panel.grid.major = element_line(linetype = "dashed")
    )
  
  # Create a subdirectory for the plots
  plot_dir <- fs::path(output_dir, "Metabolite_Contribution_Plots")
  if (!fs::dir_exists(plot_dir)) dir_create(plot_dir)
  
  ggsave(filename = fs::path(plot_dir, paste0("Contribution_Plot_", current_metabolite, ".png")),
         plot = bubble_plot, width = 16, height = 10, dpi = 1200)
  
  print(paste("✅ Plot saved for", current_metabolite))
}

# 5. WRITE THE FINAL SUMMARY REPORT
# ----------------------------------------------------
output_file_path <- fs::path(output_dir, "Functional_Contribution_Summary.txt")
writeLines(unlist(report_lines), con = output_file_path)

print(paste("\n--- 🎉 Success! Your comprehensive report has been saved to:", output_file_path, "---"))
print(paste("--- All individual metabolite plots are saved in:", fs::path(output_dir, "Metabolite_Contribution_Plots"), "---"))

# --- Script End ---