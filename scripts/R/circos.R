# Title: Definitive Script for Phylum-to-Metabolite Circos Diagram (FINAL)
# Author: abhilasha 
# Date: 2025-09-22
# Institute: NABI, Mohali, India
# Description: This script creates an enhanced Circos diagram. This final version
#              uses the correct, up-to-date arguments for circlize functions.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(circlize)
library(RColorBrewer)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_CIRCOS_enhanced")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. LOAD AND PREPARE DATA (No changes in this section)
# ----------------------------------------------------
# --- Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:133)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)

# --- Taxonomy Mapping ---
taxonomy_file <- fs::path(base_path, "taxonomy_mapping.csv")
if (!fs::file_exists(taxonomy_file)) {
  stop("--> 'taxonomy_mapping.csv' not found. Please create this file.")
}
taxonomy_map <- read_csv(taxonomy_file, show_col_types = FALSE)

# --- Abundance Data ---
all_abundances_files <- fs::dir_ls(path = diet_path, glob = "**/updated_abundances*.csv", recurse = TRUE)
all_abundances <- all_abundances_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(.x, show_col_types = FALSE, progress = FALSE) %>%
      mutate(
        diet = fs::path_split(.x)[[1]] %>% purrr::pluck(length(.) - 2),
        sample_id = sample_id,
        Species = stringr::str_replace(Model, "\\.mat$", "") %>% stringr::str_replace_all("_", " ")
      )
  })
merged_abundance_data <- all_abundances %>%
  inner_join(metadata, by = "sample_id") %>%
  inner_join(taxonomy_map, by = "Species")
print("✅ Microbial abundance data loaded and mapped to phyla.")

# --- Metabolite Data ---
all_scfa_files <- fs::dir_ls(path = diet_path, glob = "**/SCFA_Harmful_Metabolite_Profile*.csv", recurse = TRUE)
all_metabolite_data <- all_scfa_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(.x, show_col_types = FALSE, progress = FALSE) %>%
      select(starts_with("Metabolite_"), starts_with("WeightedFlux_")) %>%
      pivot_longer(everything(), names_to = c(".value", "id"), names_pattern = "(.+)_(\\d+)") %>%
      filter(!is.na(Metabolite)) %>%
      mutate(
        diet = fs::path_split(.x)[[1]] %>% purrr::pluck(length(.) - 2),
        sample_id = sample_id,
        Metabolite = str_to_title(str_replace_all(Metabolite, "_", " "))
      )
  })
merged_flux_data <- all_metabolite_data %>%
  inner_join(metadata, by = "sample_id")
print("✅ Metabolite flux data loaded.")


# 4. CALCULATE CONTRIBUTIONS (No changes in this section)
# ----------------------------------------------------
key_metabolites <- c("Ammonia", "Butyrate", "Succinate", "Hydrogen Sulfide", "Trimethylamine", "Propionate", "Acetate")

contribution_data <- merged_abundance_data %>%
  inner_join(
    merged_flux_data %>% filter(Metabolite %in% key_metabolites),
    by = c("sample_id", "diet", "group"),
    relationship = "many-to-many"
  ) %>%
  mutate(contribution_score = RelAbundance * WeightedFlux) %>%
  group_by(group, Phylum, Metabolite) %>%
  summarise(total_contribution = sum(contribution_score, na.rm = TRUE), .groups = "drop") %>%
  filter(total_contribution > 0)
print("✅ Phylum-to-Metabolite contributions calculated.")


# 5. GENERATE ENHANCED CIRCOS DIAGRAMS
# ----------------------------------------------------
for (current_group in c("CRC", "Healthy")) {
  
  cat(paste("\n--- Generating ENHANCED outputs for:", current_group, "---\n"))
  
  # Prepare data
  adjacency_matrix <- contribution_data %>%
    filter(group == current_group) %>%
    select(Phylum, Metabolite, total_contribution) %>%
    pivot_wider(names_from = Metabolite, values_from = total_contribution, values_fill = 0) %>%
    column_to_rownames("Phylum")
  
  # Filter out phyla with negligible contributions
  original_phyla_count <- nrow(adjacency_matrix)
  phylum_totals <- rowSums(adjacency_matrix)
  adjacency_matrix <- adjacency_matrix[phylum_totals > 1e-6, ]
  
  if(nrow(adjacency_matrix) < original_phyla_count) {
    removed_count <- original_phyla_count - nrow(adjacency_matrix)
    cat(paste0("INFO: Removed ", removed_count, " phyla with negligible contribution for the '", current_group, "' group.\n"))
  }
  
  # Save data to text file
  output_text_file <- fs::path(output_dir, paste0("Figure_Data_", current_group, ".txt"))
  write.csv(adjacency_matrix, file = output_text_file, quote = FALSE)
  print(paste("✅ Data for", current_group, "group saved to:", fs::path_file(output_text_file)))
  
  # Define Colors
  phylum_colors <- brewer.pal(max(3, nrow(adjacency_matrix)), "Paired")[1:nrow(adjacency_matrix)]
  names(phylum_colors) <- rownames(adjacency_matrix)
  metabolite_colors <- brewer.pal(max(3, ncol(adjacency_matrix)), "Pastel1")[1:ncol(adjacency_matrix)]
  names(metabolite_colors) <- colnames(adjacency_matrix)
  grid_colors <- c(phylum_colors, metabolite_colors)
  
  # --- Start plotting device ---
  png(fs::path(output_dir, paste0("Figure_Circos_Diagram_Enhanced_", current_group, ".png")),
      width = 12, height = 12, units = "in", res = 1200)
  
  circos.par(
    gap.after = c(rep(2, nrow(adjacency_matrix) - 1), 10, rep(2, ncol(adjacency_matrix) - 1), 10),
    start.degree = 90
  )
  
  chordDiagram(
    as.matrix(adjacency_matrix),
    grid.col = grid_colors,
    transparency = 0.25,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.2)
  )
  
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector_name = CELL_META$sector.index
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      circos.text(
        mean(xlim),
        ylim[1] + 0.1,
        sector_name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.2,
        font = 2
      )
      circos.axis(
        h = "top",
        labels.cex = 0.7,
        major.tick.length = 1, #<-- CORRECTED: Replaced deprecated argument
        sector.index = sector_name,
        track.index = 2
      )
      circos.barplot(
        value = xlim[2],
        pos = 0,
        col = grid_colors[sector_name],
        border = "white"
        #<-- CORRECTED: Removed unused 'width' argument
      )
    },
    bg.border = NA
  )
  
  title(paste("Phylum-to-Metabolite Contributions in", current_group, "Group (Enhanced)"))
  dev.off()
  circos.clear()
  
  print(paste("✅ Enhanced Circos diagram saved for", current_group, "group."))
}

# --- Script End ---