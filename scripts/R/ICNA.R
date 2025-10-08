# Title: Integrative Correlation Network Analysis (Enhanced Visualization)
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This definitive script builds and visualizes a correlation network,
#              with improved aesthetics and readability for publication-quality figures.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(rstatix)
library(Hmisc)
library(igraph)
library(ggraph)
library(RColorBrewer) # For more color options


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_NETWORK_REVISED") # New output directory for revised plots
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


# 4. IDENTIFY SIGNIFICANT FEATURES (WITH FIX)
# ----------------------------------------------------
# --- Species ---
species_stat_results <- merged_abundance_data %>%
  group_by(diet, Species) %>%
  # **FIX**: Ensure at least one group has variance before testing
  filter(n_distinct(group) == 2 &&
           (var(RelAbundance[group == "CRC"], na.rm = TRUE) > 0 || var(RelAbundance[group == "Healthy"], na.rm = TRUE) > 0)) %>%
  rstatix::wilcox_test(RelAbundance ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

# --- Metabolites ---
metabolite_stat_results <- merged_flux_data %>%
  group_by(diet, Metabolite) %>%
  # **FIX**: Ensure at least one group has variance before testing
  filter(n_distinct(group) == 2 &&
           (var(WeightedFlux[group == "CRC"], na.rm = TRUE) > 0 || var(WeightedFlux[group == "Healthy"], na.rm = TRUE) > 0)) %>%
  rstatix::wilcox_test(WeightedFlux ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

sig_species <- species_stat_results %>% filter(p.adj < 0.05) %>% pull(Species) %>% unique()
sig_metabolites <- metabolite_stat_results %>% filter(p.adj < 0.05) %>% pull(Metabolite) %>% unique()
print(paste("Found", length(sig_species), "significant species and", length(sig_metabolites), "significant metabolites."))


# 5. CREATE THE NETWORK
# ----------------------------------------------------
# Define correlation threshold
correlation_threshold <- 0.5 # INCREASE THIS VALUE (e.g., to 0.6 or 0.7) to simplify the network

if(length(sig_species) > 1 && length(sig_metabolites) > 1) {
  # --- Prepare Correlation Data ---
  wide_species <- merged_abundance_data %>%
    filter(Species %in% sig_species) %>%
    pivot_wider(id_cols = sample_diet_id, names_from = Species, values_from = RelAbundance, values_fill = 0)
  wide_metabolites <- merged_flux_data %>%
    filter(Metabolite %in% sig_metabolites) %>%
    pivot_wider(id_cols = sample_diet_id, names_from = Metabolite, values_from = WeightedFlux, values_fill = 0)
  corr_data_matrix <- inner_join(wide_species, wide_metabolites, by = "sample_diet_id")
  
  # --- Calculate Spearman Correlation ---
  correlation_results <- Hmisc::rcorr(
    as.matrix(corr_data_matrix %>% select(all_of(sig_species))),
    as.matrix(corr_data_matrix %>% select(all_of(sig_metabolites))),
    type = "spearman"
  )
  rho_matrix <- correlation_results$r
  p_matrix <- correlation_results$P
  
  # --- Build Edge and Node Lists ---
  corr_table <- rho_matrix %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to", values_to = "rho") %>%
    left_join(p_matrix %>% as.data.frame() %>% rownames_to_column("from") %>% pivot_longer(-from, names_to = "to", values_to = "p_value"), by = c("from", "to"))
  
  # Filter for strong, significant correlations
  edge_list <- corr_table %>% filter(abs(rho) > correlation_threshold, p_value < 0.05)
  write_csv(edge_list, fs::path(output_dir, "network_edge_list.csv"))
  
  # Create a node list of only the nodes present in our filtered edge list
  nodes_in_network <- unique(c(edge_list$from, edge_list$to))
  node_list <- tibble(name = nodes_in_network) %>%
    mutate(type = if_else(name %in% sig_species, "Microbe", "Metabolite"))
  
  # --- Create Graph Object ---
  network_graph <- igraph::graph_from_data_frame(d = edge_list, vertices = node_list, directed = FALSE)
  print("✅ Network graph object created.")
  
  # 6. VISUALIZE THE NETWORK
  # ----------------------------------------------------
  # Define custom colors
  node_colors <-c("Microbe" = "#0F828C", "Metabolite" = "#ff9f43") # Green for Microbes, Yellow for Metabolites
  edge_gradient_colors <- c("Negative" = "#6D94C5", "Neutral" = "#E8DFCA", "Positive" = "#F08787") # Blue for negative, Red for positive
  
  network_plot <- ggraph(network_graph, layout = "stress") + # Experiment with "stress", "fr", "kk", "nicely"
    # Edges
    geom_edge_link(aes(color = rho, width = abs(rho)), alpha = 0.8) +
    scale_edge_color_gradient2(low = edge_gradient_colors["Negative"], mid = edge_gradient_colors["Neutral"], high = edge_gradient_colors["Positive"], midpoint = 0, name = "Spearman Rho") +
    scale_edge_width_continuous(range = c(0.5, 3), name = "Correlation Strength") + # Vary width by strength
    # Nodes
    geom_node_point(aes(fill = type, size = degree(network_graph)), shape = 21, color = "black", stroke = 0.8) +
    scale_fill_manual(values = node_colors, name = "Node Type") +
    scale_size(range = c(5, 15), name = "Connections (Degree)") + # Adjusted size range
    # Labels
    geom_node_text(aes(label = name), repel = TRUE, size = 7, family = "sans", bg.color = "white", bg.r = 0.1) + # Added repel and background for readability
    # Theme
    theme_graph(base_family = "sans", foreground = "grey70", plot_margin = margin(10, 10, 10, 10)) + # Light background for text, more margin
    labs(title = "Integrative Network of Significant Species and Metabolites")
  
  ggsave(filename = fs::path(output_dir, "Figure10_Correlation_Network_Revised.png"),
         plot = network_plot, width = 20, height = 15, dpi = 1200)
  
  print("✅ Figure 10 (Correlation Network) saved with revised aesthetics.")
} else {
  print("⚠️ Network analysis skipped: not enough significant features to build a network.")
}

# --- Script End ---
print("--- 🎉 Network analysis complete! Check your new output folder for revised plots! ---")