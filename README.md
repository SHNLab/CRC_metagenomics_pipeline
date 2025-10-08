# Microbiome and Metabolic Modeling Project

This repository contains the data, results, and scripts for the comparative analysis of colorectal cancer (CRC) patients and healthy individuals, integrating metagenomic, genome-scale metabolic modeling, and FBA simulations.

---

## Folder Structure

data/ # Raw and processed data used for analysis
results/ # Outputs generated from analysis and simulations
scripts/ # R scripts and FBA code for processing and simulations


---

## Folder Descriptions

### data/
This folder contains all raw and processed data necessary for reproducible microbiome and metabolic modeling analysis.  
**Note:** Only a subset of 5 CRC and 5 healthy samples is included here for demonstration purposes, to reduce repository size and make it easier to download and navigate. Full datasets can be retrieved from the NCBI SRA as indicated in the methodology.  

Subfolders include:  
- `abundance/` – Species- and genus-level abundance tables from Kraken2/Bracken  
- `multi_qc/` - Contains aggregated quality control (QC) reports generated using **MultiQC**, summarizing outputs from various tools (e.g., FastQC, Trimmomatic, MetaSPAdes, Bowtie2). Provides an overview of read quality, assembly statistics, and alignment metrics across all samples.
- `diets/` – Diet-specific flux files and exchange reaction constraints for FBA  
- `model/` – Genome-scale metabolic models and patient-specific community models  
- `raw_fastq/` – Raw paired-end sequencing reads (.fastq)  
- `read me.txt` – Original notes describing folder contents  

### results/
This folder contains the outputs generated from processing the data, including:  
- `abundances/` – Processed species- and genus-level abundance tables for each sample  
- `FBA/` – Short-chain fatty acid (SCFA) profiles and diet-specific simulation outputs  
- `kraken2_output/` – Kraken2 read-level taxonomic assignments and summary reports  

### scripts/
This folder contains scripts used to process the data and run simulations, divided into:  
- `R/` – R scripts for abundance processing, statistical analysis, and visualizations  
- `FBA/` – MATLAB or R-based scripts for running flux balance analysis simulations  

---

## Usage
1. `data/` provides all necessary input files for reproducible analysis.  
2. `scripts/` contains the processing and modeling scripts that take input data to generate results.  
3. `results/` stores all outputs, which can be directly used for analysis, figures, and downstream interpretation.

**Note:** The subset of data in this repository is intended for demonstration, testing, and reproducibility purposes. Users analyzing the full cohort should download the complete datasets from the NCBI SRA as described in the methodology.
