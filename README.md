# Drought_intensity-frequency_root_microbiome
Data and R scripts from a factorial drought experiment investigating the effects of drought intensity and frequency on root-associated microbiomes.

## Overview
This repository contains data and R code for analyzing how drought intensity, drought frequency, plant community composition and their interactions shape root-associated microbial communities and their links to plant traits and performance in a temperate dry grassland experiment.

The workflow integrates bacterial and fungal community analyses, functional predictions (PICRUSt2, FUNGuild), indicator species analysis, and structural equation modeling (SEM).

---

## Repository structure
--code/                     # R scripts for full analysis workflow  
--data/  
-----data_raw/              # Raw input tables (ASV, taxonomy)  
-----phyloseq/              # Processed phyloseq objects (.rds)  
-----picrust2/              # PICRUSt2 input and annotation files  
--results/  
------figures/              # Final figures (PNG and PDF)  
------tables/               # Statistical outputs and processed tables  
--.gitignore  
--LICENSE  
--README.md  

---

## Analysis workflow
The analysis is organized into sequential R scripts located in `code/`. Scripts are numbered to reflect the workflow:

- **00–03**: Data preprocessing and alpha diversity analysis  
- **04–06**: Beta diversity (PCoA, PERMANOVA) and visualization  
- **07–08**: Taxonomic composition and top taxa analyses  
- **09–10**: Functional prediction (Bacteria) using PICRUSt2  
- **11–12**: Fungal trophic guild analysis (FUNGuild)  
- **13–14**: Indicator ASV analyses  
- **15–17**: Plant–microbe integration and structural equation modeling (SEM)  

Run scripts in numerical order to reproduce all results.

---

## Data description
The repository includes:

- Amplicon data:
  - ASV tables (bacteria and fungi)
  - Taxonomic assignments

- Processed objects:
  - Phyloseq objects for downstream analyses

- Functional annotations:
  - PICRUSt2 predicted pathways (MetaCyc-based)
  - NSTI scores and pathway classification

- Fungal functional groups:
  - FUNGuild trophic mode assignments

- SEM-integrated dataset (including plant traits)

---

## Outputs
All outputs are stored in `results/`:

- Figures (`results/figures/`):
  - Alpha and beta diversity
  - Taxonomic composition
  - Functional pathways
  - Indicator ASVs
  - Plant–microbe relationships

- Tables (`results/tables/`):
  - Statistical model results (ANOVA, PERMANOVA, SEM)
  - Diversity metrics
  - Taxonomic summaries
  - Functional pathway analyses

---

## Reproducibility
To reproduce the analyses:

1. Clone or download this repository  
2. Open the project in R or RStudio  
3. Run scripts in the `code/` folder in order (00 → 17)  

All results and figures can be regenerated from the provided data and scripts.

---

## Requirements
- R (≥ 4.0 recommended)

Main packages used include:
- phyloseq, vegan  
- dplyr, tidyr, ggplot2  
- piecewiseSEM, car, broom  
- Hmisc, corrplot, ggpubr  

---

## Associated study
Yang, G. et al.  
**Drought intensity rather than frequency structures root-associated microbial communities and their links to plant performance** 

---

## Data availability
All data required to reproduce the analyses are included in this repository.

Sequencing data are available at:
NCBI BioProject: PRJNA1455234

---

## License
This project is licensed under the MIT License.

---

## Contact
Gang Yang, [yangg@natur.cuni.cz]
