# Diet Data

This folder contains dietary input files used for Flux Balance Analysis (FBA) simulations of individual microbiome models. Each file represents a specific diet type and defines flux constraints for metabolites or nutrients relevant to the simulation.

## Folder Structure
diet/
├── FIBER_FLUX.tsv
├── HFD_FLUX.tsv
├── VEGAN_FLUX.tsv
├── PROTEIN_FLUX.tsv
├── MEDT_FLUX.tsv
├── UNHEAL_FLUX.tsv
└── GLUTENFREE_FLUX.tsv

## File Description

- Each file is a **tab-separated values (TSV)** file containing the following columns:
  - `exchange_reaction` – Identifier of the exchange reaction in the AGORA model (e.g., `EX_glc__D_e` for glucose uptake).
  - `fluxval` –  Allowed flux (minimum uptake).

## Usage

1. Each diet file is used to **set constraints on exchange reactions** in FBA or BacArena simulations.
2. By applying different diet files, you can simulate the impact of dietary interventions on:
   - Microbial community composition
   - Production of metabolites such as short-chain fatty acids (SCFAs)
   - Metabolic cross-feeding between species
3. Ensure the **exchange reactions in the diet file match the reactions present in the model** to avoid errors.

## Notes

- Naming convention should clearly reflect the diet type (e.g., `HFD_FLUX.tsv` for high-fat diet).
- Keep all diet files in this folder for **consistency and reproducibility**.