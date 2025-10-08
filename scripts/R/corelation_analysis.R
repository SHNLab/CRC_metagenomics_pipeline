# Title: Definitive Diet-Specific Integrative Correlation Network Analysis
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This definitive script builds and visualizes a separate correlation
#              network for each diet with enhanced, bold, black labels and saves
#              a corresponding CSV file for each network. It includes a robust
#              fix to handle cases with only one significant metabolite.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have these libraries installed:
# install.packages(c("tidyverse", "fs", "rstatix", "Hmisc", "igraph", "ggraph", "RColorBrewer"))

library(tidyverse)
library(fs)
library(rstatix)
library(Hmisc)
library(igraph)
library(ggraph)
library(RColorBrewer)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_DIET_NETWORKS")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. METADATA AND DATA PREPARATION (Self-contained)
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


# 4. LOOP THROUGH EACH DIET TO GENERATE A SEPARATE NETWORK
# ----------------------------------------------------
all_diets <- unique(merged_abundance_data$diet)

for (current_diet in all_diets) {
  
  cat(paste("\n--- Processing Network for:", current_diet, "Diet ---\n"))
  
  abundance_subset <- merged_abundance_data %>% filter(diet == current_diet)
  flux_subset <- merged_flux_data %>% filter(diet == current_diet)
  
  species_stat_results <- abundance_subset %>%
    group_by(Species) %>%
    filter(n_distinct(group) == 2 && (var(RelAbundance[group == "CRC"], na.rm = TRUE) > 0 || var(RelAbundance[group == "Healthy"], na.rm = TRUE) > 0)) %>%
    rstatix::wilcox_test(RelAbundance ~ group) %>%
    rstatix::adjust_pvalue(method = "fdr")
  
  metabolite_stat_results <- flux_subset %>%
    group_by(Metabolite) %>%
    filter(n_distinct(group) == 2 && (var(WeightedFlux[group == "CRC"], na.rm = TRUE) > 0 || var(WeightedFlux[group == "Healthy"], na.rm = TRUE) > 0)) %>%
    rstatix::wilcox_test(WeightedFlux ~ group) %>%
    rstatix::adjust_pvalue(method = "fdr")
  
  sig_species <- species_stat_results %>% filter(p.adj < 0.05) %>% pull(Species) %>% unique()
  sig_metabolites <- metabolite_stat_results %>% filter(p.adj < 0.05) %>% pull(Metabolite) %>% unique()
  
  print(paste("Found", length(sig_species), "significant species and", length(sig_metabolites), "significant metabolites for", current_diet, "diet."))
  
  if(length(sig_species) > 1 && length(sig_metabolites) > 0) {
    wide_species <- abundance_subset %>%
      filter(Species %in% sig_species) %>%
      pivot_wider(id_cols = sample_id, names_from = Species, values_from = RelAbundance, values_fill = 0)
    wide_metabolites <- flux_subset %>%
      filter(Metabolite %in% sig_metabolites) %>%
      pivot_wider(id_cols = sample_id, names_from = Metabolite, values_from = WeightedFlux, values_fill = 0)
    corr_data_matrix <- inner_join(wide_species, wide_metabolites, by = "sample_id")
    
    # **FIX**: Handle correlation calculation differently if there is only 1 significant metabolite
    if (length(sig_metabolites) > 1) {
      correlation_results <- Hmisc::rcorr(
        as.matrix(corr_data_matrix %>% select(all_of(sig_species))),
        as.matrix(corr_data_matrix %>% select(all_of(sig_metabolites))),
        type = "spearman"
      )
      rho_matrix <- correlation_results$r
      p_matrix <- correlation_results$P
    } else {
      # Manual correlation if only one metabolite
      rho_vec <- c()
      p_vec <- c()
      for(s in sig_species) {
        test <- cor.test(corr_data_matrix[[s]], corr_data_matrix[[sig_metabolites]], method = "spearman")
        rho_vec <- c(rho_vec, test$estimate)
        p_vec <- c(p_vec, test$p.value)
      }
      rho_matrix <- matrix(rho_vec, nrow = length(sig_species), ncol = 1, dimnames = list(sig_species, sig_metabolites))
      p_matrix <- matrix(p_vec, nrow = length(sig_species), ncol = 1, dimnames = list(sig_species, sig_metabolites))
    }
    
    corr_table <- rho_matrix %>%
      as.data.frame() %>%
      rownames_to_column("from") %>%
      pivot_longer(-from, names_to = "to", values_to = "rho") %>%
      left_join(p_matrix %>% as.data.frame() %>% rownames_to_column("from") %>% pivot_longer(-from, names_to = "to", values_to = "p_value"), by = c("from", "to"))
    
    edge_list <- corr_table %>% filter(abs(rho) > 0.5, p_value < 0.05)
    
    # **NEW**: Save the edge list CSV for the current diet
    write_csv(edge_list, fs::path(output_dir, paste0("Network_Edge_List_", current_diet, ".csv")))
    print(paste("✅ Network edge list saved for", current_diet, "diet."))
    
    if(nrow(edge_list) == 0) {
      print("No strong correlations found for this diet. Skipping network plot.")
      next
    }
    
    nodes_in_network <- unique(c(edge_list$from, edge_list$to))
    node_list <- tibble(name = nodes_in_network) %>%
      mutate(type = if_else(name %in% sig_species, "Microbe", "Metabolite"))
    
    network_graph <- igraph::graph_from_data_frame(d = edge_list, vertices = node_list, directed = FALSE)
    
    # --- Visualize the Network with Improved Labels ---
    network_plot <- ggraph(network_graph, layout = "fr") +
      geom_edge_link(aes(color = rho, width = abs(rho)), alpha = 0.7) +
      scale_edge_color_gradient2(low = "#6D94C5", mid = "#E8DFCA", high = "#F08787", name = "Spearman Rho", limits = c(-1, 1)) +
      scale_edge_width_continuous(range = c(0.5, 3), name = "Correlation Strength") +
      geom_node_point(aes(fill = type, size = degree(network_graph)), shape = 21, color = "black") +
      scale_fill_manual(values = c("Microbe" = "#0F828C", "Metabolite" = "#ff9f43"), name = "Node Type") +
      scale_size(range = c(5, 15), name = "Connections (Degree)") +
      # **FIX**: Enhanced node labels
      geom_node_text(aes(label = name), repel = TRUE, size = 14/.pt, fontface = "bold", color = "black", bg.color = "white", bg.r = 0.1) +
      theme_graph(base_family = "sans") +
      labs(title = paste("Integrative Network for", current_diet, "Diet"))
    
    ggsave(filename = fs::path(output_dir, paste0("Figure_Network_", current_diet, ".png")),
           plot = network_plot, width = 20, height = 15, dpi = 1200)
    
    print(paste("✅ Figure (Network Plot) saved for", current_diet, "diet."))
  } else {
    print("⚠️ Network analysis skipped: not enough significant features to build a network.")
  }
}

# --- Script End ---
print("\n--- 🎉 All diet-specific network analyses complete! ---")