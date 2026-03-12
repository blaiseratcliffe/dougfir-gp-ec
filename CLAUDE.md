# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Genomic prediction framework for coastal Douglas-fir (*Pseudotsuga menziesii*) in British Columbia. The repo tracks **R and batch scripts only** — no data, figures, or model outputs. All large files (GDS, genotype matrices, phenotypes, ClimateNA extracts) live outside the repo on OneDrive.

## Running the pipeline

The QC pipeline is the primary entry point:

```bat
run_qc.bat
```

This calls `Rscript snp_qc_pipeline_v2.R` and writes a timestamped log to `data/qc/qc_report_*.txt`. The batch file hardcodes `R-4.4.2` and the work-machine OneDrive path.

Individual R scripts are run via `Rscript <script>.R` or `source()` within an R session. There is no build system, package manager, or test suite.

## Pipeline architecture

The scripts form a linear pipeline. Each script's output is the next script's input:

1. **`pre-qc_diag.R`** — Pre-QC diagnostics on raw VCF via vcfR. Computes sample/site missingness, het rates, depth, allele balance, batch effects, PCA, and KGD-like GRM. Outputs CSVs + diagnostic plots to guide threshold selection. Does not filter.

2. **`vcf_to_gds.R`** — Converts raw VCF to SeqArray GDS (full-fidelity: all INFO/FORMAT fields). Optionally recompresses to LZMA. Creates a derived SNPRelate GDS (biallelic SNPs, QC-filtered) for HWE/GRM. Uses `storage.option = "ZIP_RA"` as a string (not `seqStorageOption()`) to avoid a known parallel-import bug.

3. **`validate_gds_1.R`** — Validates VCF header completeness and confirms SeqArray GDS retains all declared INFO/FORMAT keys.

4. **`validate_gds_2.R`** — Spot-checks VCF vs GDS fidelity by comparing random variant subsets (chrom/pos, INFO, FORMAT fields).

5. **`build_pedigree_metadata.R`** — Constructs the 3-generation pedigree (G0 founders -> G1 parents -> G2 progeny) from `All_sites.csv` and `founders_climateNA.csv`. Outputs `pedigree_full.csv` and `pedigree_genotyped.csv`.

6. **`snp_qc_pipeline_v2.R`** (v4.6.1) — Main QC pipeline. Reads from SeqArray GDS. Staged filtering: A (genotype DP/GQ masking) -> B (dead samples) -> C (site call rate + QUAL) -> D (excess het / paralog, G1 only) -> E (refined sample QC) -> F (MAF). Post-QC: pedigree panel, IBS0 verification, Mendelian error check, duplicate detection. Outputs filtered 0/1/2 matrix and optional dosage matrix.

7. **`imputation_pipeline.R`** — Takes post-QC GDS and runs four imputation methods: mean imputation, LinkImputeR (Java), gtImputation (Python/SOM), and KGD (R). Builds GRM for each. Requires external tools in `tools/` directory (JDK, Python, LinkImputeR.jar, gtImputation repo, KGD repo).

## Key data paths (not in repo)

- **Work machine**: `D:\OneDrive - NRCan RNCan\gs\doug-fir\` (GDS files, QC outputs, tools/)
- **Home machine**: `D:\OneDrive\projects\df\`
- **GDS on NVMe** (work): `C:\Users\bratclif\dfir_gds\`

Scripts hardcode the work-machine paths. When running at home, paths need manual adjustment.

## Key R packages

- **Bioconductor**: SeqArray, SNPRelate
- **CRAN**: vcfR, dplyr, readr, stringr, ggplot2, scales, gridExtra, AGHmatrix, tibble
- **Requires R >= 4.3**

## Biological context

- Diploid conifer with high genetic load and rapid LD decay — excess heterozygosity is expected and not always paralog signal.
- 3-generation pedigree: 45 G0 founders (ungenotyped) -> ~1,377 G1 (genotyped) -> 490 G2 (genotyped).
- Variant calls from freeBayes GBS against a scaffold-level reference. freeBayes uses RO/AO (not standard AD) for allele depths.
- QC thresholds (v4.6) were deliberately relaxed from v4.5 to recover more SNPs (~12K was too few). The het-ratio filter was tightened to 1.25 to compensate.

## Git workflow

The repo uses a deny-all `.gitignore` that whitelists only `*.R`, `*.bat`, `.gitignore`, `README.md`, and `LICENSE`. The `tools/` directory is explicitly blocked. Sync scripts (`git_sync_home.bat`, `git_sync_work.bat`) handle push/pull between work and home machines via GitHub.
