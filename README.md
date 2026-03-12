# Douglas-fir Genomic Selection

Genomic prediction framework for coastal Douglas-fir (*Pseudotsuga menziesii*) in British Columbia, integrating single-step GBLUP, GxE reaction norm modeling, and SNP panel optimization.

## Scope

This repository tracks **R and batch scripts only** — no data, no figures, no model outputs.
All large files (GDS, genotype matrices, phenotype files, ClimateNA extracts) live outside the repo.

## Pipeline components

| Component | Description |
|-----------|-------------|
| **SNP QC** | Staged quality control pipeline (v4.x) for freeBayes GBS variant calls against a scaffold-level reference. Reads from SeqArray GDS; outputs filtered 0/1/2 and dosage matrices. |
| **Site analysis** | breedR spatial mixed models for EP.708 progeny trial sites. Stage 1 BLUP extraction with AR1×AR1 or P-spline spatial adjustment. |
| **Genomic prediction** | BGLR stage 2 models: ABLUP, ssGBLUP, reaction norm (Jarquín), enviromic kernel. |
| **SNP panel optimization** | AIM-selected, GWAS-selected, and random panels at multiple density levels with nested cross-validation. |

## Data (not tracked)

Primary data paths (local):
- `C:/Users/bratclif/dfir_gds/` — GDS files, QC outputs
- `D:/OneDrive - NRCan RNCan/gs/doug-fir/data/` — phenotype masterfiles, pedigree, ClimateNA

## Requirements

- R ≥ 4.3
- Key packages: SeqArray, SNPRelate, breedR, BGLR, AGHmatrix, EnvRtype, vcfR, ggplot2, dplyr
