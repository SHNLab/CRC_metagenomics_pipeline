# Raw FASTQ data

This folder contains the metadata and instructions for downloading the raw metagenomic sequencing data.

## Data source
- CRC patients: BioProject PRJNA531273 (30 samples)
- Healthy controls: BioProject PRJNA397112 (110 samples)
- Sequencing platform: Illumina NextSeq 500, paired-end WGS

## Download instructions
1. Install SRA Toolkit (v3.1.0).
2. Download SRR files using prefetch:
		prefetch --verify yes <SRR_ID>
3. Convert SRR to FASTQ using fasterq-dump:
		fasterq-dump <SRR_ID>.sra -O <output_directory>
4. Organize files in subfolders by sample ID:
		raw_fastq/CRC/<SRR_ID>/<SRR_ID>_1.fastq
		raw_fastq/CRC/<SRR_ID>/<SRR_ID>_2.fastq
		raw_fastq/Healthy/<SRR_ID>/<SRR_ID>_1.fastq
		raw_fastq/Healthy/<SRR_ID>/<SRR_ID>_2.fastq


## Metadata
- A CSV file containing all sample IDs and group (CRC or Healthy) is included.
