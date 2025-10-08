# Results

This folder contains the processed outputs of the metagenomic analysis pipeline for colorectal cancer (CRC) patients and healthy controls. It is structured to provide easy access to taxonomic, functional, and metabolic data for downstream analyses.

## Folder Structure

- **abundances/**  
  Contains species- and genus-level relative abundance tables derived from Kraken2 and Bracken outputs.

- **FBA/**  
  Contains outputs from flux balance analysis (FBA) simulations, including predicted short-chain fatty acid (SCFA) profiles and updated microbial abundances for each dietary condition.

- **kraken2_output/**  
  Contains Kraken2 classification outputs and reports for individual samples.

## File Description

### abundances/
- Species- and genus-level abundance tables for individual samples.
- Merged abundance tables for CRC and healthy cohorts.

### FBA/
- SCFA profiles and abundance updates for six different dietary simulations.
- Each diet has a separate subfolder containing its respective outputs.

### kraken2_output/
- Kraken2 output files for each sample.
- Kraken2 reports required for Bracken abundance estimation.

## Usage
These results can be used for:
1. Comparative microbial community analysis between CRC patients and healthy individuals.
2. Assessing the metabolic output of microbial communities under different dietary scenarios.
3. Integrative analysis of taxonomic and functional potential in gut microbiomes.

