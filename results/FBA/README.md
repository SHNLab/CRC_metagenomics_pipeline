# Flux Balance Analysis (FBA) Results

This directory contains the FBA simulation outputs categorized by diet type.

Each diet folder includes two main subdirectories:
1. **Updated_abundances_<DIET>/** — Predicted microbial abundances after FBA simulation.
2. **Metabolite_profile_<DIET>/** — Corresponding metabolite and SCFA flux profiles.

---

## Directory Structure
FBA/
├── FIBER/
│ ├── Updated_abundances_FIBER/
│ └── Metabolite_profile_FIBER/
├── MEDT/
│ ├── Updated_abundances_MEDT/
│ └── Metabolite_profile_MEDT/
├── HFD/
│ ├── Updated_abundances_HFD/
│ └── Metabolite_profile_HFD/
├── GLUTENFREE/
│ ├── Updated_abundances_GLUTENFREE/
│ └── Metabolite_profile_GLUTENFREE/
├── PROTEIN/
│ ├── Updated_abundances_PROTEIN/
│ └── Metabolite_profile_PROTEIN/
└── VEGAN/
├── Updated_abundances_VEGAN/
└── Metabolite_profile_VEGAN/
---

## Notes
- Flux values are derived from patient-specific microbial community models constrained with diet-specific flux inputs.
- Simulations were executed using MICOM/COBRA-based optimization pipelines.

---

## Recommended Use
Use these data for:
- Comparing microbial community adaptation across diets.
- Correlating diet-induced changes in abundances with SCFA production.
- Integrating with metadata (clinical, phenotypic) for downstream statistical or pathway enrichment analysis.
