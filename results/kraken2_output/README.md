# Kraken2 Output

This folder contains taxonomic classification results generated using **Kraken2** for both colorectal cancer (CRC) patients and healthy controls. The outputs provide read-level taxonomic assignments and summary reports necessary for downstream abundance estimation and microbiome analysis.

---

## Folder Structure
kraken2_output/
├── CRC/
│ ├── SRR8865572/
│ │ └── SRR8865572_report.txt # Summary report of taxonomic assignments
│ ├── SRR8865573/
│ └── ...
└── Healthy/
├── SRR5898908/
│ └── SRR5898908_report.txt
├── SRR5898909/
└── ...


File Description

*_report.txt: Summary of reads assigned to each taxonomic rank (kingdom, phylum, genus, species).

Usage

Input for Bracken to generate species- and genus-level relative abundance tables.

Used for downstream taxonomic profiling and diversity analyses.