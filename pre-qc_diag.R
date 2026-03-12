#!/usr/bin/env Rscript
# ==============================================================================
# VCF Deep Summary & Statistics — Pre-QC Diagnostics (v2.5)
# ==============================================================================
# Project : Douglas-fir Genomic Selection
# File    : UBC_071001_snps_RAW.vcf
#
# ------------------------------------------------------------------------------
# OVERVIEW / WHAT THIS SCRIPT DOES
# ------------------------------------------------------------------------------
# This script is a **pre-QC diagnostic pass** over a raw VCF. It is intended to
# help you choose sensible QC thresholds (sample and SNP filters) and to flag
# obvious data issues *before* you run your main filtering pipeline.
#
# It reads a VCF with diploid GT calls and computes:
#
#   A) BASIC FILE + VARIANT SUMMARY
#      • Sample count, variant count, per-chromosome site counts, VCF size
#      • Variant type breakdown (SNP/Insertion/Deletion/MNP/Complex/multiallelic)
#      • Ts/Tv ratio for **canonical biallelic SNPs** (A/C/G/T only)
#
#   B) GENOTYPE-LEVEL + SAMPLE-LEVEL QC DIAGNOSTICS
#      • Overall call rate and genotype class proportions (0/0, 0/1, 1/1, other)
#      • Per-sample missingness (all variants)
#      • Per-sample heterozygosity (canonical biallelic SNPs only)
#      • Outlier detection: samples > 3 SD high missingness / high het / low het
#
#   C) SITE-LEVEL QC DIAGNOSTICS
#      • Per-site missingness / call rate (all variants)
#      • For canonical biallelic SNPs only:
#          - Allele frequency and MAF
#          - HWE chi-square screen (rough diagnostic; not exact test)
#          - KGD-inspired "Fin" metric + excess_het_flag
#
#   D) DEPTH / QUALITY METRICS (IF PRESENT IN FORMAT)
#      • FORMAT/DP: overall + per-sample mean DP + per-site mean DP
#      • FORMAT/GQ: overall + per-sample mean GQ (if present)
#      • FORMAT/AD: allele balance diagnostics (if present)
#
#   E) ALLELE BALANCE + CONTAMINATION DIAGNOSTICS (AD-based)
#      • Per-site median allele balance at het calls
#      • Per-sample mean allele balance at het calls
#      • Minor allele contamination at homozygous sites
#
#   F) SITE FREQUENCY SPECTRUM + SNP CLUSTERING
#      • Folded SFS with singleton/doubleton enrichment check
#      • SNP proximity / clustering diagnostic
#
#   G) DEPTH-DEPENDENT GENOTYPE BIAS
#      • Het rate stratified by DP bins
#
#   H) PLATE LAYOUT / BATCH EFFECT DIAGNOSTICS
#      • Empty well flagging
#      • Batch-effect boxplots by sequencing plate and S-group
#      • DNA input quality vs QC metrics
#      • Well-position effects (spatial heatmap)
#
#   I) PCA DIAGNOSTIC
#      • PCA on dosage matrix colored by plate, S-group, missingness, depth
#
#   J) KGD-LIKE GRM DIAGNOSTIC (RELATEDNESS SANITY CHECK)
#      • Missingness-adjusted weighted GRM on a SNP subset
#      • Off-diagonal histogram + diagonal vs depth
#
#   K) PLOTS + CSV OUTPUTS
#      • Histograms, scatterplots, heatmaps supporting threshold selection
#      • Two CSVs: per-sample and per-site summary tables
#
# ------------------------------------------------------------------------------
# ASSUMPTIONS / INTERPRETATION NOTES
# ------------------------------------------------------------------------------
# • Diploid nuclear genotypes are assumed.
# • "Canonical biallelic SNPs" means single-base A/C/G/T REF and ALT, not
#   multi-allelic.
# • HWE is a rough chi-square screen; interpret cautiously in structured data.
# • Fin / excess_het_flag is a heuristic for paralog detection.
# • GRM diagnostic is a sanity check, not a production GRM pipeline.
# • Allele balance diagnostics require the AD (allelic depth) FORMAT field.
#
# ------------------------------------------------------------------------------
# CHANGES IN v2.5:
#   • Plate layout integration (batch effects, empty well flagging, well-position)
#   • PCA diagnostic (colored by plate, S-group, missingness, depth)
#   • Allele balance at het sites (AD-based) + per-site/per-sample summaries
#   • Minor allele contamination at homozygous sites (AD-based)
#   • Folded site frequency spectrum with singleton enrichment check
#   • SNP clustering / proximity diagnostic
#   • Depth-dependent genotype bias (het rate vs DP bins)
#   • DNA input quality (concentration/volume) vs QC metrics
# ==============================================================================

cat("
╔══════════════════════════════════════════════════════════════════════╗
║  VCF Deep Summary & Pre-QC Diagnostics (v2.5)                       ║
╚══════════════════════════════════════════════════════════════════════╝
\n")

# ── 0. Setup ─────────────────────────────────────────────────────────────────
required_pkgs <- c("vcfR", "ggplot2", "scales", "gridExtra")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(vcfR)
library(ggplot2)
library(scales)
library(gridExtra)

# ── Helper functions ─────────────────────────────────────────────────────────
safe_mean   <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_sd     <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)
safe_min    <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max    <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

which_min_safe <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  which.min(replace(x, is.na(x), Inf))
}

which_max_safe <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  which.max(replace(x, is.na(x), -Inf))
}

is_missing_gt_matrix <- function(gt) {
  dot_mask <- grepl("\\.", gt)
  dim(dot_mask) <- dim(gt)
  dimnames(dot_mask) <- dimnames(gt)
  is.na(gt) | dot_mask
}

print_meta_lines <- function(meta, label) {
  if (is.null(meta) || length(meta) == 0) {
    cat(sprintf("  No %s fields found.\n", label))
    return(invisible(NULL))
  }
  if (!is.null(dim(meta))) {
    cat(sprintf("  %s fields present: %d\n", label, nrow(meta)))
    for (i in seq_len(nrow(meta))) {
      cat(sprintf("    • %s\n", paste(meta[i, ], collapse = " | ")))
    }
  } else {
    cat(sprintf("  %s fields present: %d\n", label, length(meta)))
    for (i in seq_along(meta)) {
      cat(sprintf("    • %s\n", paste(meta[[i]], collapse = " | ")))
    }
  }
  invisible(NULL)
}

# ── 1. Configuration ─────────────────────────────────────────────────────────
vcf_path <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/UBC_071001_snps_RAW.vcf"

out_dir <- file.path(dirname(vcf_path), "pre_qc_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── KGD-inspired diagnostics toggles ─────────────────────────────────────────
MAKE_FIN_PLOT  <- TRUE
FIN_RATIO_THR  <- 1.5
FIN_MAF_MIN    <- 0.05
FIN_HOBS_MIN   <- 0.20

FLAG_SITE_DP_Q <- 0.995

DP_BAR_TOP_N         <- 80
DP_BAR_SHOW_TOP_HIGH <- TRUE

# KGD-like GRM diagnostic
RUN_KGD_GRM_DIAG   <- TRUE
KGD_SAMPLE_MISS_MAX <- 0.20
KGD_SAMPLE_DP_MIN   <- 20
KGD_MAX_SAMPLES     <- Inf
KGD_MAX_SNPS        <- 50000
KGD_SEED            <- 1

# ── v2.5: Plate layout diagnostics ──────────────────────────────────────────
RUN_PLATE_DIAG    <- TRUE
plate_layout_path <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/__UBC_071001_PlateLayout.xls"
# readxl is needed for .xls files; install if plate diagnostics are on
if (RUN_PLATE_DIAG && grepl("\\.(xls|xlsx)$", plate_layout_path, ignore.case = TRUE)) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl", repos = "https://cloud.r-project.org")
  }
}

# ── v2.5: PCA diagnostic ────────────────────────────────────────────────────
RUN_PCA_DIAG     <- TRUE
PCA_MAX_SNPS     <- 50000
PCA_SEED         <- 1

# ── v2.5: Allele balance diagnostics (requires AD in FORMAT) ────────────────
RUN_AB_DIAG        <- TRUE
AB_HET_EXPECTED    <- 0.5
AB_HET_OUTLIER_THR <- 0.15   # flag sites where |median_AB - 0.5| > this
AB_HOM_CONTAM_THR  <- 0.05   # flag samples where mean minor-allele frac at hom > this

# ── v2.5: Site frequency spectrum ───────────────────────────────────────────
RUN_SFS_DIAG  <- TRUE
SFS_MAX_MAC   <- 20

# ── v2.5: SNP clustering / proximity diagnostic ────────────────────────────
RUN_SNP_CLUSTER_DIAG <- TRUE
SNP_CLUSTER_WINDOWS  <- c(10, 50)
SNP_CLUSTER_THR      <- 3

# ── v2.5: Depth-dependent genotype bias ─────────────────────────────────────
RUN_DEPTH_HET_DIAG <- TRUE
DP_HET_BINS        <- c(0, 2, 5, 10, 20, 50, Inf)

# ── 2. Read VCF ──────────────────────────────────────────────────────────────
cat("── Reading VCF ─────────────────────────────────────────────────────\n")
t0 <- Sys.time()
vcf <- read.vcfR(vcf_path, verbose = TRUE)
t1 <- Sys.time()
cat(sprintf("   Read time: %.1f seconds\n\n", as.numeric(difftime(t1, t0, units = "secs"))))

# ── 3. Basic Dimensions ──────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 1: BASIC FILE DIMENSIONS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

n_samples  <- ncol(vcf@gt) - 1
n_variants <- nrow(vcf@fix)
chroms     <- unique(vcf@fix[, "CHROM"])
n_chroms   <- length(chroms)

cat(sprintf("  Total samples (individuals) : %s\n", format(n_samples, big.mark = ",")))
cat(sprintf("  Total variant sites         : %s\n", format(n_variants, big.mark = ",")))
cat(sprintf("  Unique chromosomes/scaffolds: %s\n", format(n_chroms, big.mark = ",")))
cat(sprintf("  VCF file size               : %.1f MB\n", file.info(vcf_path)$size / 1e6))

chrom_counts <- as.data.frame(table(vcf@fix[, "CHROM"]), stringsAsFactors = FALSE)
colnames(chrom_counts) <- c("CHROM", "n_variants")
chrom_counts <- chrom_counts[order(-chrom_counts$n_variants), ]

cat("\n  Top 20 chromosomes/scaffolds by variant count:\n")
print(head(chrom_counts, 20), row.names = FALSE)

# ── 4. Variant Type Breakdown ────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 2: VARIANT TYPE BREAKDOWN\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

ref_alleles <- vcf@fix[, "REF"]
alt_alleles <- vcf@fix[, "ALT"]
ref_uc <- toupper(ref_alleles)
alt_uc <- toupper(alt_alleles)

classify_variant <- function(ref, alt) {
  alts <- strsplit(alt, ",", fixed = TRUE)[[1]]
  types <- vapply(alts, function(a) {
    if (nchar(ref) == 1 && nchar(a) == 1) "SNP"
    else if (nchar(ref) > nchar(a)) "Deletion"
    else if (nchar(ref) < nchar(a)) "Insertion"
    else "MNP/Complex"
  }, character(1))
  if (length(alts) > 1) return(paste0("Multi-allelic (", paste(types, collapse = "/"), ")"))
  types
}

variant_types <- vapply(seq_len(n_variants), function(i) {
  classify_variant(ref_alleles[i], alt_alleles[i])
}, character(1))

type_table <- sort(table(variant_types), decreasing = TRUE)

cat("  Variant type counts:\n")
for (nm in names(type_table)) {
  cat(sprintf("    %-30s : %s  (%.1f%%)\n",
              nm, format(type_table[nm], big.mark = ","),
              100 * type_table[nm] / n_variants))
}

n_multiallelic <- sum(grepl(",", alt_alleles, fixed = TRUE), na.rm = TRUE)
cat(sprintf("\n  Multi-allelic sites: %s (%.2f%%)\n",
            format(n_multiallelic, big.mark = ","),
            100 * n_multiallelic / n_variants))

biallelic_snp_idx <- which(
  !is.na(ref_uc) & !is.na(alt_uc) &
    grepl("^[ACGT]$", ref_uc) & grepl("^[ACGT]$", alt_uc) &
    !grepl(",", alt_uc, fixed = TRUE)
)

cat("\n── Transition / Transversion Ratio ──\n")
if (length(biallelic_snp_idx) > 0) {
  ref_snp <- ref_uc[biallelic_snp_idx]
  alt_snp <- alt_uc[biallelic_snp_idx]
  transitions <- (ref_snp == "A" & alt_snp == "G") | (ref_snp == "G" & alt_snp == "A") |
    (ref_snp == "C" & alt_snp == "T") | (ref_snp == "T" & alt_snp == "C")
  transversions <- !transitions
  n_ts <- sum(transitions); n_tv <- sum(transversions)
  ts_tv <- if (n_tv > 0) n_ts / n_tv else NA_real_
  cat(sprintf("  Canonical biallelic SNPs : %s\n", format(length(biallelic_snp_idx), big.mark = ",")))
  cat(sprintf("  Transitions             : %s\n", format(n_ts, big.mark = ",")))
  cat(sprintf("  Transversions           : %s\n", format(n_tv, big.mark = ",")))
  cat(sprintf("  Ts/Tv ratio             : %s\n", if (is.na(ts_tv)) "NA" else sprintf("%.3f", ts_tv)))
} else {
  cat("  No canonical biallelic SNPs found for Ts/Tv calculation.\n")
}

# ── 5. Genotype Matrix Extraction ────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 3: GENOTYPE STATISTICS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
sample_names <- colnames(gt_matrix)

if (anyDuplicated(sample_names)) {
  dupes <- unique(sample_names[duplicated(sample_names)])
  stop(sprintf("Duplicate sample IDs found in VCF genotype columns: %s",
               paste(head(dupes, 20), collapse = ", ")))
}

missing_mask <- is_missing_gt_matrix(gt_matrix)

total_calls <- length(gt_matrix)
n_missing   <- sum(missing_mask)
called      <- sum(!missing_mask)
overall_call_rate <- if (total_calls > 0) 1 - (n_missing / total_calls) else NA_real_

hom_ref <- sum((gt_matrix == "0/0") | (gt_matrix == "0|0"), na.rm = TRUE)
het_all <- sum((gt_matrix == "0/1") | (gt_matrix == "1/0") |
                 (gt_matrix == "0|1") | (gt_matrix == "1|0"), na.rm = TRUE)
hom_alt <- sum((gt_matrix == "1/1") | (gt_matrix == "1|1"), na.rm = TRUE)
other_gt <- called - hom_ref - het_all - hom_alt

cat(sprintf("  Total genotype calls    : %s\n", format(total_calls, big.mark = ",")))
cat(sprintf("  Missing genotypes       : %s (%.2f%%)\n",
            format(n_missing, big.mark = ","), 100 * n_missing / total_calls))
cat(sprintf("  Overall call rate       : %.4f (%.2f%%)\n",
            overall_call_rate, 100 * overall_call_rate))

cat(sprintf("\n  Genotype class breakdown (of called genotypes):\n"))
if (called > 0) {
  cat(sprintf("    Hom-ref  (0/0) : %s  (%.2f%%)\n", format(hom_ref, big.mark = ","), 100 * hom_ref / called))
  cat(sprintf("    Het 0/1 class  : %s  (%.2f%%)\n", format(het_all, big.mark = ","), 100 * het_all / called))
  cat(sprintf("    Hom-alt  (1/1) : %s  (%.2f%%)\n", format(hom_alt, big.mark = ","), 100 * hom_alt / called))
  if (other_gt > 0) cat(sprintf("    Other/multi    : %s  (%.2f%%)\n", format(other_gt, big.mark = ","), 100 * other_gt / called))
} else {
  cat("    No called genotypes found.\n")
}

gt_bial <- if (length(biallelic_snp_idx) > 0) {
  gt_matrix[biallelic_snp_idx, , drop = FALSE]
} else {
  matrix(character(0), nrow = 0, ncol = ncol(gt_matrix), dimnames = list(character(0), sample_names))
}

gt_bial_missing <- if (nrow(gt_bial) > 0) is_missing_gt_matrix(gt_bial) else {
  matrix(FALSE, nrow = 0, ncol = ncol(gt_matrix), dimnames = list(character(0), sample_names))
}

gt_bial_het_mask <- if (nrow(gt_bial) > 0) {
  (gt_bial == "0/1") | (gt_bial == "1/0") | (gt_bial == "0|1") | (gt_bial == "1|0")
} else {
  matrix(FALSE, nrow = 0, ncol = ncol(gt_bial), dimnames = list(character(0), sample_names))
}

# ── 6. Per-Sample Statistics ─────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 4: PER-SAMPLE STATISTICS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

sample_missing   <- colSums(missing_mask)
sample_miss_rate <- sample_missing / n_variants

sample_called_bial    <- if (nrow(gt_bial) > 0) colSums(!gt_bial_missing) else rep(0L, n_samples)
sample_het_count_bial <- if (nrow(gt_bial) > 0) colSums(gt_bial_het_mask, na.rm = TRUE) else rep(0L, n_samples)
sample_het_rate_bial  <- ifelse(sample_called_bial > 0, sample_het_count_bial / sample_called_bial, NA_real_)

sample_called_all <- colSums(!missing_mask)
sample_homalt     <- colSums((gt_matrix == "1/1") | (gt_matrix == "1|1"), na.rm = TRUE)
sample_homalt_rate <- ifelse(sample_called_all > 0, sample_homalt / sample_called_all, NA_real_)

idx_min_miss <- which_min_safe(sample_miss_rate)
idx_max_miss <- which_max_safe(sample_miss_rate)
idx_min_het  <- which_min_safe(sample_het_rate_bial)
idx_max_het  <- which_max_safe(sample_het_rate_bial)

cat("  Per-sample missing rate (all variants):\n")
cat(sprintf("    Mean   : %.4f (%.2f%%)\n", safe_mean(sample_miss_rate), 100 * safe_mean(sample_miss_rate)))
cat(sprintf("    Median : %.4f (%.2f%%)\n", safe_median(sample_miss_rate), 100 * safe_median(sample_miss_rate)))
cat(sprintf("    Min    : %.4f (%s)\n", safe_min(sample_miss_rate), if (is.na(idx_min_miss)) "NA" else sample_names[idx_min_miss]))
cat(sprintf("    Max    : %.4f (%s)\n", safe_max(sample_miss_rate), if (is.na(idx_max_miss)) "NA" else sample_names[idx_max_miss]))
cat(sprintf("    SD     : %.4f\n", safe_sd(sample_miss_rate)))

cat("\n  Per-sample heterozygosity (canonical biallelic SNPs only):\n")
cat(sprintf("    Mean   : %.4f\n", safe_mean(sample_het_rate_bial)))
cat(sprintf("    Median : %.4f\n", safe_median(sample_het_rate_bial)))
cat(sprintf("    Min    : %.4f (%s)\n", safe_min(sample_het_rate_bial), if (is.na(idx_min_het)) "NA" else sample_names[idx_min_het]))
cat(sprintf("    Max    : %.4f (%s)\n", safe_max(sample_het_rate_bial), if (is.na(idx_max_het)) "NA" else sample_names[idx_max_het]))
cat(sprintf("    SD     : %.4f\n", safe_sd(sample_het_rate_bial)))

miss_mean <- safe_mean(sample_miss_rate); miss_sd <- safe_sd(sample_miss_rate)
het_mean  <- safe_mean(sample_het_rate_bial); het_sd <- safe_sd(sample_het_rate_bial)

miss_outliers  <- if (is.na(miss_sd)) character(0) else sample_names[!is.na(sample_miss_rate) & sample_miss_rate > (miss_mean + 3 * miss_sd)]
het_outliers_hi <- if (is.na(het_sd)) character(0) else sample_names[!is.na(sample_het_rate_bial) & sample_het_rate_bial > (het_mean + 3 * het_sd)]
het_outliers_lo <- if (is.na(het_sd)) character(0) else sample_names[!is.na(sample_het_rate_bial) & sample_het_rate_bial < (het_mean - 3 * het_sd)]

cat("\n  ⚠ Outlier samples (>3 SD from mean):\n")
cat(sprintf("    High missingness : %d samples\n", length(miss_outliers)))
cat(sprintf("    High het         : %d samples\n", length(het_outliers_hi)))
cat(sprintf("    Low het          : %d samples\n", length(het_outliers_lo)))

sample_stats <- data.frame(
  Sample = sample_names, N_called_all = sample_called_all,
  N_missing_all = sample_missing, Missing_rate_all = round(sample_miss_rate, 6),
  N_called_biallelic_snp = sample_called_bial, N_het_biallelic_snp = sample_het_count_bial,
  Het_rate_biallelic_snp = round(sample_het_rate_bial, 6),
  N_hom_alt_all = sample_homalt, Hom_alt_rate_all = round(sample_homalt_rate, 6),
  stringsAsFactors = FALSE
)

# ── 7. Per-Site Statistics ───────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 5: PER-SITE (PER-VARIANT) STATISTICS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

site_missing   <- rowSums(missing_mask)
site_miss_rate <- site_missing / n_samples
site_called    <- n_samples - site_missing

site_het_rate_bial <- rep(NA_real_, n_variants)
site_called_bial <- integer(0); site_het_bial <- integer(0)

if (nrow(gt_bial) > 0) {
  site_called_bial <- rowSums(!gt_bial_missing)
  site_het_bial    <- rowSums(gt_bial_het_mask, na.rm = TRUE)
  site_het_rate_bial[biallelic_snp_idx] <- ifelse(site_called_bial > 0, site_het_bial / site_called_bial, NA_real_)
}

cat("── Computing MAF + HWE + Fin metrics (canonical biallelic SNPs only) ──\n")

site_maf      <- rep(NA_real_, n_variants)
site_alt_freq <- rep(NA_real_, n_variants)
valid_pvals   <- numeric(0)

H_obs <- rep(NA_real_, n_variants); H_exp <- rep(NA_real_, n_variants)
fin_ratio <- rep(NA_real_, n_variants); excess_het_flag <- rep(FALSE, n_variants)

n_biallelic_snps <- length(biallelic_snp_idx)

# Allele count vectors (reused for SFS later)
bial_n_00 <- bial_n_01 <- bial_n_11 <- bial_n_total <- integer(0)
bial_allele_count <- bial_alt_allele_sum <- integer(0)

if (n_biallelic_snps > 0) {
  bial_n_00 <- rowSums((gt_bial == "0/0") | (gt_bial == "0|0"), na.rm = TRUE)
  bial_n_01 <- rowSums((gt_bial == "0/1") | (gt_bial == "1/0") | (gt_bial == "0|1") | (gt_bial == "1|0"), na.rm = TRUE)
  bial_n_11 <- rowSums((gt_bial == "1/1") | (gt_bial == "1|1"), na.rm = TRUE)
  
  bial_n_total       <- bial_n_00 + bial_n_01 + bial_n_11
  bial_allele_count  <- 2L * bial_n_total
  bial_alt_allele_sum <- bial_n_01 + 2L * bial_n_11
  
  snp_alt_freq <- ifelse(bial_allele_count > 0, bial_alt_allele_sum / bial_allele_count, NA_real_)
  snp_maf      <- pmin(snp_alt_freq, 1 - snp_alt_freq)
  
  site_alt_freq[biallelic_snp_idx] <- snp_alt_freq
  site_maf[biallelic_snp_idx]      <- snp_maf
  
  h_obs_vec <- ifelse(bial_n_total > 0, bial_n_01 / bial_n_total, NA_real_)
  h_exp_vec <- 2 * snp_alt_freq * (1 - snp_alt_freq)
  fin_vec   <- h_obs_vec / pmax(h_exp_vec, 1e-8)
  
  H_obs[biallelic_snp_idx]     <- h_obs_vec
  H_exp[biallelic_snp_idx]     <- h_exp_vec
  fin_ratio[biallelic_snp_idx] <- fin_vec
  
  maf_ok  <- !is.na(site_maf[biallelic_snp_idx]) & site_maf[biallelic_snp_idx] >= FIN_MAF_MIN
  hexp_ok <- !is.na(h_exp_vec) & h_exp_vec >= 0.02
  hobs_ok <- !is.na(h_obs_vec) & h_obs_vec >= FIN_HOBS_MIN
  excess_het_flag[biallelic_snp_idx] <- maf_ok & hexp_ok & hobs_ok & (fin_vec > FIN_RATIO_THR)
  
  # HWE chi-square screen
  valid_hwe <- bial_n_total >= 10
  p_ref <- ifelse(valid_hwe, (2 * bial_n_00 + bial_n_01) / (2 * bial_n_total), NA_real_)
  q_alt <- ifelse(valid_hwe, 1 - p_ref, NA_real_)
  exp_00 <- bial_n_total * p_ref^2
  exp_01 <- bial_n_total * 2 * p_ref * q_alt
  exp_11 <- bial_n_total * q_alt^2
  term_00 <- ifelse(valid_hwe & exp_00 > 0, (bial_n_00 - exp_00)^2 / exp_00, 0)
  term_01 <- ifelse(valid_hwe & exp_01 > 0, (bial_n_01 - exp_01)^2 / exp_01, 0)
  term_11 <- ifelse(valid_hwe & exp_11 > 0, (bial_n_11 - exp_11)^2 / exp_11, 0)
  chi2 <- ifelse(valid_hwe, term_00 + term_01 + term_11, NA_real_)
  hwe_pvals <- ifelse(valid_hwe, pchisq(chi2, df = 1, lower.tail = FALSE), NA_real_)
  valid_pvals <- hwe_pvals[!is.na(hwe_pvals)]
}

cat(sprintf("\n  Site missing rate (all variants):\n"))
cat(sprintf("    Mean   : %.4f (%.2f%%)\n", safe_mean(site_miss_rate), 100 * safe_mean(site_miss_rate)))
cat(sprintf("    Median : %.4f\n", safe_median(site_miss_rate)))
cat(sprintf("    Max    : %.4f\n", safe_max(site_miss_rate)))

if (n_biallelic_snps > 0) {
  cat(sprintf("\n  Minor Allele Frequency (canonical biallelic SNPs only):\n"))
  cat(sprintf("    Canonical biallelic SNPs: %s\n", format(n_biallelic_snps, big.mark = ",")))
  cat(sprintf("    Mean MAF                : %.4f\n", safe_mean(site_maf)))
  cat(sprintf("    Median MAF              : %.4f\n", safe_median(site_maf)))
}

cat("\n  Variants passing common call-rate thresholds (all variants):\n")
for (thr in c(0.50, 0.70, 0.80, 0.90, 0.95)) {
  n_pass <- sum((1 - site_miss_rate) >= thr)
  cat(sprintf("    Call rate ≥ %.0f%% : %s variants (%.1f%%)\n",
              100 * thr, format(n_pass, big.mark = ","), 100 * n_pass / n_variants))
}

rm(missing_mask); gc(verbose = FALSE)

# ── 8. Depth Statistics (DP) + AD Extraction ─────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 6: READ DEPTH STATISTICS (DP)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

dp_matrix <- tryCatch(extract.gt(vcf, element = "DP", as.numeric = TRUE), error = function(e) NULL)

sample_mean_dp <- NULL; site_mean_dp <- NULL

if (!is.null(dp_matrix) && sum(!is.na(dp_matrix)) > 0) {
  cat(sprintf("  Overall mean depth  : %.2f\n", mean(dp_matrix, na.rm = TRUE)))
  cat(sprintf("  Overall median depth: %.1f\n", median(dp_matrix, na.rm = TRUE)))
  cat(sprintf("  Min depth           : %d\n", min(dp_matrix, na.rm = TRUE)))
  cat(sprintf("  Max depth           : %d\n", max(dp_matrix, na.rm = TRUE)))
  cat(sprintf("  SD                  : %.2f\n", sd(dp_matrix, na.rm = TRUE)))
  
  sample_mean_dp <- colMeans(dp_matrix, na.rm = TRUE)
  site_mean_dp   <- rowMeans(dp_matrix, na.rm = TRUE)
  sample_mean_dp[is.nan(sample_mean_dp)] <- NA_real_
  site_mean_dp[is.nan(site_mean_dp)]     <- NA_real_
  
  idx_min_dp <- which_min_safe(sample_mean_dp)
  idx_max_dp <- which_max_safe(sample_mean_dp)
  
  cat(sprintf("\n  Per-sample mean depth:\n"))
  cat(sprintf("    Mean of means  : %.2f\n", safe_mean(sample_mean_dp)))
  cat(sprintf("    Min mean depth : %.2f (%s)\n", safe_min(sample_mean_dp),
              if (is.na(idx_min_dp)) "NA" else sample_names[idx_min_dp]))
  cat(sprintf("    Max mean depth : %.2f (%s)\n", safe_max(sample_mean_dp),
              if (is.na(idx_max_dp)) "NA" else sample_names[idx_max_dp]))
  
  cat(sprintf("\n  Per-site mean depth:\n"))
  cat(sprintf("    Mean of means  : %.2f\n", safe_mean(site_mean_dp)))
  cat(sprintf("    Median         : %.2f\n", safe_median(site_mean_dp)))
  
  cat("\n  Sites by mean depth range:\n")
  dp_breaks <- c(0, 2, 5, 10, 20, 50, 100, Inf)
  dp_labels <- c("<2", "2-5", "5-10", "10-20", "20-50", "50-100", ">100")
  dp_cut <- cut(site_mean_dp, breaks = dp_breaks, labels = dp_labels, right = FALSE)
  dp_tab <- table(dp_cut)
  for (i in seq_along(dp_tab)) {
    cat(sprintf("    %-10s : %s variants (%.1f%%)\n",
                names(dp_tab)[i], format(dp_tab[i], big.mark = ","), 100 * dp_tab[i] / n_variants))
  }
  sample_stats$Mean_DP <- round(sample_mean_dp[match(sample_stats$Sample, sample_names)], 2)
} else {
  cat("  DP field not found in genotype fields. Skipping depth analysis.\n")
}

high_depth_site_flag <- rep(FALSE, n_variants)
depth_site_q995 <- NA_real_
if (!is.null(site_mean_dp) && sum(!is.na(site_mean_dp)) > 0) {
  dp_ref <- site_mean_dp
  if (length(biallelic_snp_idx) > 0) dp_ref <- site_mean_dp[biallelic_snp_idx]
  depth_site_q995 <- as.numeric(quantile(dp_ref, probs = FLAG_SITE_DP_Q, na.rm = TRUE, names = FALSE))
  high_depth_site_flag <- !is.na(site_mean_dp) & (site_mean_dp >= depth_site_q995)
}

# ── v2.5: Extract AD matrix for allele balance diagnostics ──────────────────
ad_matrix_raw <- NULL
if (RUN_AB_DIAG) {
  cat("\n── Extracting AD (allelic depth) field ──\n")
  ad_matrix_raw <- tryCatch(
    extract.gt(vcf, element = "AD", as.numeric = FALSE),
    error = function(e) NULL
  )
  if (is.null(ad_matrix_raw) || sum(!is.na(ad_matrix_raw)) == 0) {
    cat("  AD field not found or empty. Allele balance diagnostics will be skipped.\n")
    ad_matrix_raw <- NULL
    RUN_AB_DIAG <- FALSE
  } else {
    cat(sprintf("  AD field extracted: %s entries with data.\n",
                format(sum(!is.na(ad_matrix_raw)), big.mark = ",")))
  }
}

# ── 9. QUAL Score Statistics ─────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 7: QUALITY SCORE (QUAL) STATISTICS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

qual_vals <- suppressWarnings(as.numeric(vcf@fix[, "QUAL"]))
n_qual_avail <- sum(!is.na(qual_vals))

if (n_qual_avail > 0) {
  cat(sprintf("  QUAL values available : %s / %s variants\n",
              format(n_qual_avail, big.mark = ","), format(n_variants, big.mark = ",")))
  cat(sprintf("  Mean QUAL   : %.2f\n", mean(qual_vals, na.rm = TRUE)))
  cat(sprintf("  Median QUAL : %.2f\n", median(qual_vals, na.rm = TRUE)))
  cat(sprintf("  Min QUAL    : %.2f\n", min(qual_vals, na.rm = TRUE)))
  cat(sprintf("  Max QUAL    : %.2f\n", max(qual_vals, na.rm = TRUE)))
  cat("\n  QUAL thresholds:\n")
  for (q in c(20, 30, 50, 100, 200, 500)) {
    n_pass <- sum(qual_vals >= q, na.rm = TRUE)
    cat(sprintf("    QUAL >= %4d : %s variants (%.1f%%)\n",
                q, format(n_pass, big.mark = ","), 100 * n_pass / n_variants))
  }
} else { cat("  No QUAL values available.\n") }

# ── 10. INFO Field Parsing ───────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 8: INFO FIELD OVERVIEW\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
info_fields <- queryMETA(vcf, element = "INFO")
print_meta_lines(info_fields, "INFO")

# ── 11. FORMAT Field Overview ────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 9: FORMAT FIELD OVERVIEW\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
format_fields <- queryMETA(vcf, element = "FORMAT")
print_meta_lines(format_fields, "FORMAT")

# ── 12. Genotype Quality (GQ) ────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 10: GENOTYPE QUALITY (GQ) STATISTICS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

gq_matrix <- tryCatch(extract.gt(vcf, element = "GQ", as.numeric = TRUE), error = function(e) NULL)
if (!is.null(gq_matrix) && sum(!is.na(gq_matrix)) > 0) {
  sample_mean_gq <- colMeans(gq_matrix, na.rm = TRUE)
  sample_mean_gq[is.nan(sample_mean_gq)] <- NA_real_
  idx_min_gq <- which_min_safe(sample_mean_gq)
  idx_max_gq <- which_max_safe(sample_mean_gq)
  cat(sprintf("  Overall mean GQ   : %.2f\n", mean(gq_matrix, na.rm = TRUE)))
  cat(sprintf("  Overall median GQ : %.1f\n", median(gq_matrix, na.rm = TRUE)))
  cat(sprintf("\n  Per-sample mean GQ:\n"))
  cat(sprintf("    Mean  : %.2f\n", safe_mean(sample_mean_gq)))
  cat(sprintf("    Min   : %.2f (%s)\n", safe_min(sample_mean_gq), if (is.na(idx_min_gq)) "NA" else sample_names[idx_min_gq]))
  cat(sprintf("    Max   : %.2f (%s)\n", safe_max(sample_mean_gq), if (is.na(idx_max_gq)) "NA" else sample_names[idx_max_gq]))
  cat("\n  GQ thresholds:\n")
  for (gq in c(10, 20, 30, 50)) {
    pct <- 100 * sum(gq_matrix >= gq, na.rm = TRUE) / sum(!is.na(gq_matrix))
    cat(sprintf("    GQ >= %2d : %.1f%% of genotypes\n", gq, pct))
  }
} else { cat("  GQ field not found. Skipping.\n") }
rm(gq_matrix); gc(verbose = FALSE)

# ── 13. Hardy-Weinberg Equilibrium Summary ───────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 11: HARDY-WEINBERG EQUILIBRIUM (CANONICAL BIALLELIC SNPs)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
if (n_biallelic_snps > 0) {
  cat(sprintf("  Valid HWE tests: %s\n", format(length(valid_pvals), big.mark = ",")))
  cat("  Note: rough diagnostic chi-square screen; interpret cautiously in structured data.\n")
} else { cat("  No canonical biallelic SNPs for HWE testing.\n") }

# ── 14. FILTER Field Summary ────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 12: FILTER FIELD SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
filter_vals <- vcf@fix[, "FILTER"]
filter_table <- sort(table(filter_vals, useNA = "ifany"), decreasing = TRUE)
cat("  FILTER values:\n")
for (nm in names(filter_table)) {
  cat(sprintf("    %-20s : %s (%.1f%%)\n", nm, format(filter_table[nm], big.mark = ","), 100 * filter_table[nm] / n_variants))
}

# ==============================================================================
# v2.5 NEW COMPUTATION SECTIONS (before plots)
# ==============================================================================

# ── 14b. Allele Balance Diagnostics (Section G + H) ─────────────────────────
site_median_AB  <- rep(NA_real_, n_variants)
site_AB_flag    <- rep(FALSE, n_variants)
sample_mean_AB_het        <- rep(NA_real_, n_samples)
sample_hom_minor_frac     <- rep(NA_real_, n_samples)
contam_flagged            <- character(0)

if (RUN_AB_DIAG && !is.null(ad_matrix_raw) && n_biallelic_snps > 0) {
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("  SECTION 13: ALLELE BALANCE DIAGNOSTICS (AD-based)\n")
  cat("══════════════════════════════════════════════════════════════════\n\n")
  
  # Parse AD field: "ref_count,alt_count" -> ref and alt depth matrices
  ad_bial <- ad_matrix_raw[biallelic_snp_idx, , drop = FALSE]
  
  parse_ad_ref <- function(x) {
    v <- suppressWarnings(as.integer(sub(",.*", "", x)))
    v[is.na(v)] <- NA_integer_; v
  }
  parse_ad_alt <- function(x) {
    v <- suppressWarnings(as.integer(sub(".*,", "", x)))
    v[is.na(v)] <- NA_integer_; v
  }
  
  ad_ref <- matrix(parse_ad_ref(ad_bial), nrow = nrow(ad_bial), ncol = ncol(ad_bial),
                   dimnames = dimnames(ad_bial))
  ad_alt <- matrix(parse_ad_alt(ad_bial), nrow = nrow(ad_bial), ncol = ncol(ad_bial),
                   dimnames = dimnames(ad_bial))
  ad_total <- ad_ref + ad_alt
  
  # Allele balance = ref / (ref + alt); compute only where ad_total > 0
  ab_ratio <- ifelse(ad_total > 0, ad_ref / ad_total, NA_real_)
  
  # --- Section G: AB at het sites ---
  het_mask_bial <- gt_bial_het_mask
  ab_het <- ab_ratio
  ab_het[!het_mask_bial | is.na(het_mask_bial)] <- NA_real_
  
  # Per-site median AB at het calls
  site_median_AB_bial <- apply(ab_het, 1, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) == 0) NA_real_ else median(vals)
  })
  site_median_AB[biallelic_snp_idx] <- site_median_AB_bial
  site_AB_flag[biallelic_snp_idx] <- !is.na(site_median_AB_bial) &
    abs(site_median_AB_bial - AB_HET_EXPECTED) > AB_HET_OUTLIER_THR
  
  # Per-sample mean AB at het calls
  sample_mean_AB_het <- apply(ab_het, 2, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) == 0) NA_real_ else mean(vals)
  })
  
  n_ab_flagged <- sum(site_AB_flag, na.rm = TRUE)
  cat(sprintf("  Per-site median AB at het sites:\n"))
  cat(sprintf("    Mean of medians : %.4f\n", safe_mean(site_median_AB_bial)))
  cat(sprintf("    Sites flagged (|AB-0.5| > %.2f) : %s (%.1f%%)\n",
              AB_HET_OUTLIER_THR, format(n_ab_flagged, big.mark = ","),
              100 * n_ab_flagged / n_biallelic_snps))
  
  cat(sprintf("\n  Per-sample mean AB at het sites:\n"))
  cat(sprintf("    Mean   : %.4f\n", safe_mean(sample_mean_AB_het)))
  cat(sprintf("    Min    : %.4f\n", safe_min(sample_mean_AB_het)))
  cat(sprintf("    Max    : %.4f\n", safe_max(sample_mean_AB_het)))
  
  # --- Section H: Minor allele contamination at hom sites ---
  hom_ref_mask <- (gt_bial == "0/0") | (gt_bial == "0|0")
  hom_alt_mask <- (gt_bial == "1/1") | (gt_bial == "1|1")
  
  # At hom-ref: minor allele fraction = alt / total
  # At hom-alt: minor allele fraction = ref / total
  minor_frac <- matrix(NA_real_, nrow = nrow(ad_bial), ncol = ncol(ad_bial),
                       dimnames = dimnames(ad_bial))
  hom_ref_ok <- hom_ref_mask & ad_total > 0 & !is.na(ad_total)
  hom_alt_ok <- hom_alt_mask & ad_total > 0 & !is.na(ad_total)
  minor_frac[hom_ref_ok] <- ad_alt[hom_ref_ok] / ad_total[hom_ref_ok]
  minor_frac[hom_alt_ok] <- ad_ref[hom_alt_ok] / ad_total[hom_alt_ok]
  
  sample_hom_minor_frac <- apply(minor_frac, 2, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) == 0) NA_real_ else mean(vals)
  })
  
  contam_flagged <- sample_names[!is.na(sample_hom_minor_frac) & sample_hom_minor_frac > AB_HOM_CONTAM_THR]
  
  cat(sprintf("\n  Minor allele fraction at homozygous sites (contamination check):\n"))
  cat(sprintf("    Mean per-sample minor allele frac : %.4f\n", safe_mean(sample_hom_minor_frac)))
  cat(sprintf("    Max per-sample minor allele frac  : %.4f\n", safe_max(sample_hom_minor_frac)))
  cat(sprintf("    Samples flagged (> %.2f)           : %d\n", AB_HOM_CONTAM_THR, length(contam_flagged)))
  
  rm(ad_bial, ad_ref, ad_alt, ad_total, ab_ratio, ab_het, minor_frac,
     hom_ref_mask, hom_alt_mask, hom_ref_ok, hom_alt_ok, het_mask_bial)
  gc(verbose = FALSE)
}
rm(ad_matrix_raw); gc(verbose = FALSE)

# ── 14c. Folded Site Frequency Spectrum (Section I) ──────────────────────────
sfs_mac_vec <- NULL
if (RUN_SFS_DIAG && n_biallelic_snps > 0) {
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("  SECTION 14: SITE FREQUENCY SPECTRUM (FOLDED)\n")
  cat("══════════════════════════════════════════════════════════════════\n\n")
  
  # Minor allele count per site
  mac_vec <- pmin(bial_alt_allele_sum, bial_allele_count - bial_alt_allele_sum)
  mac_vec <- mac_vec[bial_allele_count > 0]
  sfs_mac_vec <- mac_vec
  
  n_singletons  <- sum(mac_vec == 1, na.rm = TRUE)
  n_doubletons  <- sum(mac_vec == 2, na.rm = TRUE)
  n_polymorphic <- sum(mac_vec > 0, na.rm = TRUE)
  
  cat(sprintf("  Polymorphic biallelic SNPs : %s\n", format(n_polymorphic, big.mark = ",")))
  cat(sprintf("  Singletons (MAC=1)        : %s (%.1f%%)\n",
              format(n_singletons, big.mark = ","), 100 * n_singletons / max(n_polymorphic, 1)))
  cat(sprintf("  Doubletons (MAC=2)        : %s (%.1f%%)\n",
              format(n_doubletons, big.mark = ","), 100 * n_doubletons / max(n_polymorphic, 1)))
}

# ── 14d. SNP Clustering / Proximity Diagnostic (Section J) ──────────────────
site_snp_cluster <- list()
if (RUN_SNP_CLUSTER_DIAG) {
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("  SECTION 15: SNP CLUSTERING / PROXIMITY DIAGNOSTIC\n")
  cat("══════════════════════════════════════════════════════════════════\n\n")
  
  site_chrom <- vcf@fix[, "CHROM"]
  site_pos   <- suppressWarnings(as.integer(vcf@fix[, "POS"]))
  
  for (win in SNP_CLUSTER_WINDOWS) {
    col_name <- paste0("snp_cluster_", win, "bp")
    cluster_count <- integer(n_variants)
    
    # Process per chromosome for efficiency
    uchrom <- unique(site_chrom)
    for (ch in uchrom) {
      idx <- which(site_chrom == ch & !is.na(site_pos))
      if (length(idx) < 2) next
      pos_ch <- site_pos[idx]
      ord <- order(pos_ch)
      idx <- idx[ord]
      pos_ch <- pos_ch[ord]
      
      for (i in seq_along(idx)) {
        # Count SNPs within +/- win bp (excluding self)
        lower <- pos_ch[i] - win
        upper <- pos_ch[i] + win
        n_near <- sum(pos_ch >= lower & pos_ch <= upper) - 1L
        cluster_count[idx[i]] <- n_near
      }
    }
    site_snp_cluster[[col_name]] <- cluster_count
    
    n_flagged <- sum(cluster_count >= SNP_CLUSTER_THR, na.rm = TRUE)
    cat(sprintf("  Window %d bp: %s sites with >= %d nearby SNPs (%.1f%%)\n",
                win, format(n_flagged, big.mark = ","), SNP_CLUSTER_THR,
                100 * n_flagged / n_variants))
  }
}

# ── 14e. Depth-Dependent Genotype Bias (Section K) ──────────────────────────
dp_het_table <- NULL
if (RUN_DEPTH_HET_DIAG && !is.null(dp_matrix) && n_biallelic_snps > 0) {
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("  SECTION 16: DEPTH-DEPENDENT GENOTYPE BIAS\n")
  cat("══════════════════════════════════════════════════════════════════\n\n")
  
  dp_bial <- dp_matrix[biallelic_snp_idx, , drop = FALSE]
  gt_bial_flat <- as.vector(gt_bial)
  dp_bial_flat <- as.vector(dp_bial)
  
  gt_class <- rep(NA_character_, length(gt_bial_flat))
  gt_class[gt_bial_flat == "0/0" | gt_bial_flat == "0|0"] <- "hom_ref"
  gt_class[gt_bial_flat == "0/1" | gt_bial_flat == "1/0" | gt_bial_flat == "0|1" | gt_bial_flat == "1|0"] <- "het"
  gt_class[gt_bial_flat == "1/1" | gt_bial_flat == "1|1"] <- "hom_alt"
  
  dp_labels_het <- paste0("[", head(DP_HET_BINS, -1), ",", tail(DP_HET_BINS, -1), ")")
  dp_bin <- cut(dp_bial_flat, breaks = DP_HET_BINS, labels = dp_labels_het, right = FALSE)
  
  valid <- !is.na(gt_class) & !is.na(dp_bin)
  df_dh <- data.frame(dp_bin = dp_bin[valid], gt_class = gt_class[valid], stringsAsFactors = FALSE)
  
  dp_het_table <- as.data.frame.matrix(table(df_dh$dp_bin, df_dh$gt_class))
  
  cat("  Genotype rates by DP bin (biallelic SNPs):\n")
  cat(sprintf("  %-15s %10s %10s %10s %10s %10s\n", "DP_bin", "hom_ref", "het", "hom_alt", "total", "het_rate"))
  # Ensure all expected columns exist (0 if absent)
  for (gc_col in c("hom_ref", "het", "hom_alt")) {
    if (!gc_col %in% colnames(dp_het_table)) dp_het_table[[gc_col]] <- 0L
  }
  dp_het_table$total <- dp_het_table$hom_ref + dp_het_table$het + dp_het_table$hom_alt
  dp_het_table$het_rate <- ifelse(dp_het_table$total > 0, dp_het_table$het / dp_het_table$total, NA_real_)
  for (i in seq_len(nrow(dp_het_table))) {
    hr <- if (!is.na(dp_het_table$het_rate[i])) sprintf("%.4f", dp_het_table$het_rate[i]) else "NA"
    cat(sprintf("  %-15s %10s %10s %10s %10s %10s\n",
                rownames(dp_het_table)[i],
                format(dp_het_table$hom_ref[i], big.mark = ","),
                format(dp_het_table$het[i], big.mark = ","),
                format(dp_het_table$hom_alt[i], big.mark = ","),
                format(dp_het_table$total[i], big.mark = ","),
                hr))
  }
  
  rm(dp_bial, gt_bial_flat, dp_bial_flat, gt_class, dp_bin, df_dh)
  gc(verbose = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS (existing + v2.5 new)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  GENERATING DIAGNOSTIC PLOTS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

theme_diag <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey40"))

# Plot 01: MAF Distribution
p1 <- ggplot(data.frame(MAF = site_maf[!is.na(site_maf)]), aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "#2171b5", color = "white", linewidth = 0.2) +
  labs(title = "Minor Allele Frequency Distribution", subtitle = "Canonical biallelic SNPs only",
       x = "MAF", y = "Count") +
  geom_vline(xintercept = c(0.01, 0.05), linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 0.015, y = Inf, label = "MAF 0.01", vjust = 2, color = "red", size = 3) +
  annotate("text", x = 0.055, y = Inf, label = "MAF 0.05", vjust = 2, color = "red", size = 3) +
  theme_diag
ggsave(file.path(out_dir, "01_maf_distribution.png"), p1, width = 8, height = 5, dpi = 150)
cat("  ✓ 01_maf_distribution.png\n")

# Plot 02: Site Missing Rate
p2 <- ggplot(data.frame(MissRate = site_miss_rate), aes(x = MissRate)) +
  geom_histogram(bins = 50, fill = "#e6550d", color = "white", linewidth = 0.2) +
  labs(title = "Per-Site Missing Rate Distribution", subtitle = "All variants", x = "Missing Rate", y = "Count") +
  geom_vline(xintercept = c(0.10, 0.20, 0.30), linetype = "dashed", color = "grey30") +
  theme_diag
ggsave(file.path(out_dir, "02_site_missing_rate.png"), p2, width = 8, height = 5, dpi = 150)
cat("  ✓ 02_site_missing_rate.png\n")

# Plot 03: Sample Missing Rate
p3 <- ggplot(data.frame(MissRate = sample_miss_rate), aes(x = MissRate)) +
  geom_histogram(bins = 30, fill = "#756bb1", color = "white", linewidth = 0.2) +
  labs(title = "Per-Sample Missing Rate Distribution", subtitle = "All variants", x = "Missing Rate", y = "Count") +
  theme_diag
ggsave(file.path(out_dir, "03_sample_missing_rate.png"), p3, width = 8, height = 5, dpi = 150)
cat("  ✓ 03_sample_missing_rate.png\n")

# Plot 04: Sample Heterozygosity
p4 <- ggplot(data.frame(HetRate = sample_het_rate_bial), aes(x = HetRate)) +
  geom_histogram(bins = 30, fill = "#31a354", color = "white", linewidth = 0.2) +
  labs(title = "Per-Sample Heterozygosity Distribution", subtitle = "Canonical biallelic SNPs only",
       x = "Heterozygosity Rate", y = "Count") + theme_diag
if (!is.na(het_mean)) p4 <- p4 + geom_vline(xintercept = het_mean, linetype = "dashed", color = "black")
ggsave(file.path(out_dir, "04_sample_heterozygosity.png"), p4, width = 8, height = 5, dpi = 150)
cat("  ✓ 04_sample_heterozygosity.png\n")

# Plot 05: Het vs Missing
p5 <- ggplot(data.frame(MissRate = sample_miss_rate, HetRate = sample_het_rate_bial),
             aes(x = MissRate, y = HetRate)) +
  geom_point(alpha = 0.6, color = "#2171b5") +
  labs(title = "Sample Heterozygosity vs Missing Rate", x = "Missing Rate", y = "Heterozygosity Rate") +
  theme_diag
if (!is.na(het_mean) && !is.na(het_sd)) {
  p5 <- p5 + geom_hline(yintercept = het_mean + c(-3, 3) * het_sd, linetype = "dashed", color = "red", alpha = 0.5)
}
ggsave(file.path(out_dir, "05_sample_het_vs_missing.png"), p5, width = 8, height = 6, dpi = 150)
cat("  ✓ 05_sample_het_vs_missing.png\n")

# Plot 06: Depth Distribution + top-N barplots
if (!is.null(dp_matrix) && sum(!is.na(dp_matrix)) > 0) {
  dp_flat <- as.vector(dp_matrix); dp_flat <- dp_flat[!is.na(dp_flat)]
  if (length(dp_flat) > 1e6) dp_flat <- sample(dp_flat, 1e6)
  p6 <- ggplot(data.frame(DP = dp_flat), aes(x = DP)) +
    geom_histogram(bins = 100, fill = "#636363", color = "white", linewidth = 0.2) +
    labs(title = "Read Depth (DP) Distribution", x = "Read Depth", y = "Count") +
    coord_cartesian(xlim = c(0, quantile(dp_flat, 0.99))) +
    geom_vline(xintercept = c(5, 10), linetype = "dashed", color = "red", alpha = 0.6) + theme_diag
  ggsave(file.path(out_dir, "06_depth_distribution.png"), p6, width = 8, height = 5, dpi = 150)
  cat("  ✓ 06_depth_distribution.png\n")
  
  dp_df <- data.frame(Sample = sample_names, MeanDP = sample_mean_dp)
  dp_df <- dp_df[!is.na(dp_df$MeanDP), ]
  dp_df_low <- head(dp_df[order(dp_df$MeanDP), ], min(DP_BAR_TOP_N, nrow(dp_df)))
  dp_df_low$Sample <- factor(dp_df_low$Sample, levels = rev(dp_df_low$Sample))
  p6b <- ggplot(dp_df_low, aes(x = Sample, y = MeanDP)) +
    geom_bar(stat = "identity", fill = "#636363") + coord_flip() +
    labs(title = sprintf("Per-Sample Mean Read Depth (Lowest %d)", nrow(dp_df_low)), x = "", y = "Mean DP") +
    theme_diag + theme(axis.text.y = element_text(size = 7))
  ggsave(file.path(out_dir, "06b_sample_mean_depth_lowN.png"), p6b,
         width = 8, height = 0.15 * nrow(dp_df_low) + 3, dpi = 150, limitsize = FALSE)
  cat("  ✓ 06b_sample_mean_depth_lowN.png\n")
  
  if (isTRUE(DP_BAR_SHOW_TOP_HIGH)) {
    dp_df_high <- head(dp_df[order(-dp_df$MeanDP), ], min(DP_BAR_TOP_N, nrow(dp_df)))
    dp_df_high$Sample <- factor(dp_df_high$Sample, levels = rev(dp_df_high$Sample))
    p6c <- ggplot(dp_df_high, aes(x = Sample, y = MeanDP)) +
      geom_bar(stat = "identity", fill = "#636363") + coord_flip() +
      labs(title = sprintf("Per-Sample Mean Read Depth (Highest %d)", nrow(dp_df_high)), x = "", y = "Mean DP") +
      theme_diag + theme(axis.text.y = element_text(size = 7))
    ggsave(file.path(out_dir, "06c_sample_mean_depth_highN.png"), p6c,
           width = 8, height = 0.15 * nrow(dp_df_high) + 3, dpi = 150, limitsize = FALSE)
    cat("  ✓ 06c_sample_mean_depth_highN.png\n")
  }
}

# Sample depth dependence plots
if (!is.null(sample_mean_dp) && sum(!is.na(sample_mean_dp)) > 0) {
  df_depth_sample <- data.frame(Sample = sample_names, Mean_DP = sample_mean_dp,
                                Missing_rate = sample_miss_rate, Het_rate_bial = sample_het_rate_bial)
  p_dp_het <- ggplot(df_depth_sample, aes(x = Mean_DP, y = Het_rate_bial)) +
    geom_point(alpha = 0.6) + labs(title = "Sample Mean Depth vs Heterozygosity", x = "Mean DP", y = "Het Rate") + theme_diag
  ggsave(file.path(out_dir, "11_sample_depth_vs_het.png"), p_dp_het, width = 8, height = 6, dpi = 150)
  cat("  ✓ 11_sample_depth_vs_het.png\n")
  p_dp_miss <- ggplot(df_depth_sample, aes(x = Mean_DP, y = Missing_rate)) +
    geom_point(alpha = 0.6) + labs(title = "Sample Mean Depth vs Missingness", x = "Mean DP", y = "Missing Rate") + theme_diag
  ggsave(file.path(out_dir, "12_sample_depth_vs_missing.png"), p_dp_miss, width = 8, height = 6, dpi = 150)
  cat("  ✓ 12_sample_depth_vs_missing.png\n")
}

# Plot 07: QUAL Distribution
qual_plot_vals <- qual_vals[!is.na(qual_vals) & qual_vals > 0]
if (length(qual_plot_vals) > 0) {
  p7 <- ggplot(data.frame(QUAL = qual_plot_vals), aes(x = QUAL)) +
    geom_histogram(bins = 100, fill = "#de2d26", color = "white", linewidth = 0.2) +
    labs(title = "Variant Quality Score (QUAL) Distribution", x = "QUAL", y = "Count") +
    scale_x_log10(labels = comma) + geom_vline(xintercept = c(20, 30), linetype = "dashed", color = "grey30") + theme_diag
  ggsave(file.path(out_dir, "07_qual_distribution.png"), p7, width = 8, height = 5, dpi = 150)
  cat("  ✓ 07_qual_distribution.png\n")
}

# Plot 08: HWE p-value distribution
if (length(valid_pvals) > 0) {
  plot_pvals <- pmax(valid_pvals, .Machine$double.xmin)
  p8 <- ggplot(data.frame(pval = plot_pvals), aes(x = -log10(pval))) +
    geom_histogram(bins = 50, fill = "#fd8d3c", color = "white", linewidth = 0.2) +
    labs(title = "HWE p-value Distribution (Canonical Biallelic SNPs)", x = expression(-log[10](p)), y = "Count") +
    geom_vline(xintercept = -log10(1e-5), linetype = "dashed", color = "red") +
    annotate("text", x = -log10(1e-5) + 0.3, y = Inf, label = "p=1e-5", vjust = 2, color = "red", size = 3) + theme_diag
  ggsave(file.path(out_dir, "08_hwe_pvalue_distribution.png"), p8, width = 8, height = 5, dpi = 150)
  cat("  ✓ 08_hwe_pvalue_distribution.png\n")
}

# Plot 09: Variants per chromosome
chrom_top <- head(chrom_counts, 30)
chrom_top$CHROM <- factor(chrom_top$CHROM, levels = rev(chrom_top$CHROM))
p9 <- ggplot(chrom_top, aes(x = CHROM, y = n_variants)) +
  geom_bar(stat = "identity", fill = "#2171b5") + coord_flip() +
  labs(title = "Variants per Chromosome/Scaffold (Top 30)", x = "", y = "Number of Variants") +
  scale_y_continuous(labels = comma) + theme_diag
ggsave(file.path(out_dir, "09_variants_per_chrom.png"), p9, width = 8, height = 7, dpi = 150)
cat("  ✓ 09_variants_per_chrom.png\n")

# Plot 10: Diagnostic panel
panel_plots <- list(p1, p2, p3, p4)
p_panel <- arrangeGrob(grobs = panel_plots, ncol = 2, top = "Pre-QC Diagnostic Panel - Douglas-fir VCF")
ggsave(file.path(out_dir, "10_diagnostic_panel.png"), p_panel, width = 14, height = 10, dpi = 150)
cat("  ✓ 10_diagnostic_panel.png\n")

# Plot 13: Fin plot
if (MAKE_FIN_PLOT && n_biallelic_snps > 0) {
  df_fin <- data.frame(p_alt = site_alt_freq[biallelic_snp_idx], H_obs = H_obs[biallelic_snp_idx],
                       fin_ratio = fin_ratio[biallelic_snp_idx], excess_het_flag = excess_het_flag[biallelic_snp_idx])
  df_fin <- df_fin[!is.na(df_fin$p_alt) & !is.na(df_fin$H_obs), ]
  if (nrow(df_fin) > 200000) df_fin <- df_fin[sample.int(nrow(df_fin), 200000), ]
  p_fin <- ggplot(df_fin, aes(x = p_alt, y = H_obs)) +
    geom_point(alpha = 0.15, size = 0.7) +
    stat_function(fun = function(p) 2 * p * (1 - p), linewidth = 0.9, alpha = 0.85) +
    labs(title = "KGD Fin Plot: Observed vs Expected Heterozygosity",
         x = "Alt allele frequency (p)", y = "Observed heterozygosity (Hobs)") + theme_diag
  if (any(df_fin$excess_het_flag, na.rm = TRUE)) {
    p_fin <- p_fin + geom_point(data = df_fin[df_fin$excess_het_flag, ], aes(x = p_alt, y = H_obs),
                                alpha = 0.5, size = 0.9, color = "red")
  }
  ggsave(file.path(out_dir, "13_fin_plot_excess_het.png"), p_fin, width = 8, height = 6, dpi = 150)
  cat("  ✓ 13_fin_plot_excess_het.png\n")
}

# ── v2.5 NEW PLOTS ──────────────────────────────────────────────────────────

# Plot 26: Allele balance site histogram
if (RUN_AB_DIAG && sum(!is.na(site_median_AB)) > 0) {
  df_ab_site <- data.frame(medAB = site_median_AB[!is.na(site_median_AB)])
  p_ab_site <- ggplot(df_ab_site, aes(x = medAB)) +
    geom_histogram(bins = 50, fill = "#1b9e77", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(0.5 - AB_HET_OUTLIER_THR, 0.5 + AB_HET_OUTLIER_THR),
               linetype = "dotted", color = "red") +
    labs(title = "Per-Site Median Allele Balance at Het Calls",
         subtitle = "Canonical biallelic SNPs; red = flag threshold",
         x = "Median ref/(ref+alt)", y = "Count") + theme_diag
  ggsave(file.path(out_dir, "26_allele_balance_site_hist.png"), p_ab_site, width = 8, height = 5, dpi = 150)
  cat("  ✓ 26_allele_balance_site_hist.png\n")
  
  # Plot 27: Per-sample mean AB
  df_ab_samp <- data.frame(meanAB = sample_mean_AB_het[!is.na(sample_mean_AB_het)])
  p_ab_samp <- ggplot(df_ab_samp, aes(x = meanAB)) +
    geom_histogram(bins = 30, fill = "#d95f02", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    labs(title = "Per-Sample Mean Allele Balance at Het Calls",
         subtitle = "Canonical biallelic SNPs",
         x = "Mean ref/(ref+alt)", y = "Count") + theme_diag
  ggsave(file.path(out_dir, "27_allele_balance_sample_hist.png"), p_ab_samp, width = 8, height = 5, dpi = 150)
  cat("  ✓ 27_allele_balance_sample_hist.png\n")
}

# Plot 28: Hom contamination histogram
if (RUN_AB_DIAG && sum(!is.na(sample_hom_minor_frac)) > 0) {
  df_contam <- data.frame(frac = sample_hom_minor_frac[!is.na(sample_hom_minor_frac)])
  p_contam <- ggplot(df_contam, aes(x = frac)) +
    geom_histogram(bins = 30, fill = "#7570b3", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = AB_HOM_CONTAM_THR, linetype = "dashed", color = "red") +
    labs(title = "Minor Allele Fraction at Homozygous Sites (Contamination Check)",
         subtitle = "Per-sample mean; red = flag threshold",
         x = "Mean minor allele fraction", y = "Count") + theme_diag
  ggsave(file.path(out_dir, "28_hom_contamination_hist.png"), p_contam, width = 8, height = 5, dpi = 150)
  cat("  ✓ 28_hom_contamination_hist.png\n")
}

# Plot 29: Folded SFS
if (RUN_SFS_DIAG && !is.null(sfs_mac_vec)) {
  mac_plot <- sfs_mac_vec[sfs_mac_vec > 0]
  mac_capped <- ifelse(mac_plot > SFS_MAX_MAC, SFS_MAX_MAC + 1L, mac_plot)
  mac_labels <- c(as.character(1:SFS_MAX_MAC), paste0(">", SFS_MAX_MAC))
  mac_factor <- factor(mac_capped, levels = 1:(SFS_MAX_MAC + 1), labels = mac_labels)
  df_sfs <- data.frame(MAC = mac_factor)
  
  # Expected neutral SFS reference: proportional to 1/k
  obs_counts <- table(mac_factor)
  k_vals <- 1:SFS_MAX_MAC
  neutral_prop <- (1 / k_vals) / sum(1 / k_vals)
  neutral_expected <- neutral_prop * sum(obs_counts[1:SFS_MAX_MAC])
  
  p_sfs <- ggplot(df_sfs, aes(x = MAC)) +
    geom_bar(fill = "#e7298a", color = "white", linewidth = 0.2) +
    labs(title = "Folded Site Frequency Spectrum",
         subtitle = "Canonical biallelic SNPs; line = expected neutral (1/k)",
         x = "Minor Allele Count (MAC)", y = "Count") +
    theme_diag + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Overlay neutral expectation for MAC 1..SFS_MAX_MAC
  df_neutral <- data.frame(MAC = factor(1:SFS_MAX_MAC, levels = 1:(SFS_MAX_MAC + 1), labels = mac_labels),
                           expected = neutral_expected)
  p_sfs <- p_sfs + geom_point(data = df_neutral, aes(x = MAC, y = expected), color = "black", size = 2) +
    geom_line(data = df_neutral, aes(x = as.numeric(MAC), y = expected), color = "black", linewidth = 0.7)
  ggsave(file.path(out_dir, "29_folded_sfs.png"), p_sfs, width = 10, height = 5, dpi = 150)
  cat("  ✓ 29_folded_sfs.png\n")
}

# Plot 30: SNP clustering histogram
if (RUN_SNP_CLUSTER_DIAG && length(site_snp_cluster) > 0) {
  for (col_name in names(site_snp_cluster)) {
    cc <- site_snp_cluster[[col_name]]
    cc_capped <- pmin(cc, 10)
    df_clust <- data.frame(n_nearby = factor(cc_capped, levels = 0:10,
                                             labels = c(as.character(0:9), "10+")))
    p_clust <- ggplot(df_clust, aes(x = n_nearby)) +
      geom_bar(fill = "#66a61e", color = "white", linewidth = 0.2) +
      geom_vline(xintercept = SNP_CLUSTER_THR + 0.5, linetype = "dashed", color = "red") +
      labs(title = paste("SNP Clustering:", col_name),
           subtitle = sprintf("Count of nearby SNPs within window; red = flag threshold (%d)", SNP_CLUSTER_THR),
           x = "Number of nearby SNPs", y = "Sites") +
      theme_diag
    fname <- paste0("30_snp_clustering_", col_name, ".png")
    ggsave(file.path(out_dir, fname), p_clust, width = 8, height = 5, dpi = 150)
    cat(sprintf("  ✓ %s\n", fname))
  }
}

# Plot 31: Depth-dependent het bias
if (RUN_DEPTH_HET_DIAG && !is.null(dp_het_table) && "het_rate" %in% colnames(dp_het_table)) {
  df_dh_plot <- data.frame(DP_bin = factor(rownames(dp_het_table), levels = rownames(dp_het_table)),
                           het_rate = dp_het_table$het_rate, total = dp_het_table$total)
  p_dh <- ggplot(df_dh_plot, aes(x = DP_bin, y = het_rate)) +
    geom_col(fill = "#a6761d", color = "white", linewidth = 0.2) +
    geom_text(aes(label = paste0("n=", format(total, big.mark = ","))), vjust = -0.3, size = 2.5) +
    labs(title = "Depth-Dependent Genotype Bias",
         subtitle = "Het rate at biallelic SNPs stratified by genotype DP",
         x = "DP bin", y = "Het rate") +
    theme_diag
  ggsave(file.path(out_dir, "31_depth_het_bias.png"), p_dh, width = 8, height = 5, dpi = 150)
  cat("  ✓ 31_depth_het_bias.png\n")
}

# ── 15. Summary Report ───────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  SECTION 17: QC RECOMMENDATIONS BASED ON DIAGNOSTICS\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

cat("  Based on the pre-QC diagnostics, consider these thresholds:\n\n")
cat("  1. SAMPLE FILTERING\n")
if (!is.na(miss_mean) && !is.na(miss_sd)) {
  cat(sprintf("     - Remove samples with missing rate > %.2f (3 SD above mean)\n",
              min(0.50, miss_mean + 3 * miss_sd)))
}
cat(sprintf("     - Investigate %d high-het and %d low-het outlier samples\n",
            length(het_outliers_hi), length(het_outliers_lo)))
cat("     - Check per-sample depth for low-coverage individuals\n")
if (RUN_AB_DIAG && length(contam_flagged) > 0) {
  cat(sprintf("     - Investigate %d samples flagged for hom-site contamination\n", length(contam_flagged)))
}
cat("\n  2. VARIANT FILTERING\n")
cat(sprintf("     - Site call rate  : >= 80%% (removes %d variants)\n",
            n_variants - sum((1 - site_miss_rate) >= 0.80)))
if (n_biallelic_snps > 0) cat("     - MAF             : >= 0.01 or 0.05 depending on analysis\n")
cat("     - Fin/excess-het   : consider excluding excess_het_flag loci (paralog candidates)\n")
if (RUN_AB_DIAG) cat(sprintf("     - Allele balance   : consider excluding AB-flagged sites (%d flagged)\n",
                             sum(site_AB_flag, na.rm = TRUE)))
if (!is.null(site_mean_dp) && sum(!is.na(site_mean_dp)) > 0) {
  cat(sprintf("     - High site DP     : quantile %.3f (cutoff=%.2f)\n", FLAG_SITE_DP_Q, depth_site_q995))
}
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# KGD-like GRM diagnostic (FULL cleaned sample set)
# ══════════════════════════════════════════════════════════════════════════════
grm_K <- NULL          # stored for PCA reuse
grm_keep_names <- NULL

if (RUN_KGD_GRM_DIAG && n_biallelic_snps > 1000) {
  cat("\n── Running KGD-like GRM diagnostic (full cleaned sample set) ──\n")
  set.seed(KGD_SEED)
  
  keep_samp <- (sample_miss_rate <= KGD_SAMPLE_MISS_MAX)
  if (!is.null(sample_mean_dp) && sum(!is.na(sample_mean_dp)) > 0) {
    keep_samp <- keep_samp & (!is.na(sample_mean_dp)) & (sample_mean_dp >= KGD_SAMPLE_DP_MIN)
  }
  keep_names <- sample_names[keep_samp]
  cat(sprintf("  Samples passing diag filter: %d / %d\n", length(keep_names), length(sample_names)))
  
  if (length(keep_names) >= 50) {
    if (is.finite(KGD_MAX_SAMPLES) && length(keep_names) > KGD_MAX_SAMPLES) {
      keep_names <- sample(keep_names, KGD_MAX_SAMPLES)
    }
    snp_idx_sub <- biallelic_snp_idx
    if (is.finite(KGD_MAX_SNPS) && length(snp_idx_sub) > KGD_MAX_SNPS) {
      snp_idx_sub <- sample(snp_idx_sub, KGD_MAX_SNPS)
    }
    cat(sprintf("  SNPs used for diag: %d\n", length(snp_idx_sub)))
    
    gt_sub  <- gt_matrix[snp_idx_sub, keep_names, drop = FALSE]
    miss_sub <- is_missing_gt_matrix(gt_sub)
    
    g <- matrix(NA_real_, nrow = nrow(gt_sub), ncol = ncol(gt_sub), dimnames = dimnames(gt_sub))
    g[gt_sub == "0/0" | gt_sub == "0|0"] <- 0
    g[gt_sub == "0/1" | gt_sub == "1/0" | gt_sub == "0|1" | gt_sub == "1|0"] <- 1
    g[gt_sub == "1/1" | gt_sub == "1|1"] <- 2
    
    n_called_grm <- rowSums(!is.na(g))
    p_alt_grm    <- rowSums(g, na.rm = TRUE) / pmax(2 * n_called_grm, 1)
    
    X <- sweep(g, 1, 2 * p_alt_grm, FUN = "-"); X[is.na(X)] <- 0
    w <- 2 * p_alt_grm * (1 - p_alt_grm); w[!is.finite(w)] <- 0
    
    I_called <- !miss_sub & !is.na(g)
    Num <- crossprod(X)
    Z <- sweep(I_called, 1, sqrt(w), `*`)
    Den <- crossprod(Z)
    K <- Num / pmax(Den, 1e-8)
    
    grm_K <- K; grm_keep_names <- keep_names
    
    saveRDS(K, file.path(out_dir, "kgd_grm_cleaned_samples.rds"))
    cat("  ✓ kgd_grm_cleaned_samples.rds\n")
    
    K_off <- K[upper.tri(K)]
    df_k <- data.frame(K = K_off[is.finite(K_off)])
    p_k_hist <- ggplot(df_k, aes(x = K)) +
      geom_histogram(bins = 60, color = "white", linewidth = 0.2) +
      labs(title = "KGD-like GRM Diagnostic (off-diagonals)",
           subtitle = sprintf("Samples=%d, SNPs=%d", ncol(gt_sub), nrow(gt_sub)),
           x = "Relatedness (scaled)", y = "Count") + theme_diag
    ggsave(file.path(out_dir, "14_kgd_grm_offdiag_hist.png"), p_k_hist, width = 8, height = 5, dpi = 150)
    cat("  ✓ 14_kgd_grm_offdiag_hist.png\n")
    
    if (!is.null(sample_mean_dp) && sum(!is.na(sample_mean_dp)) > 0) {
      dp_keep <- sample_mean_dp[match(keep_names, sample_names)]
      df_diag <- data.frame(Sample = keep_names, Mean_DP = dp_keep, K_diag = diag(K))
      p_k_diag <- ggplot(df_diag, aes(x = Mean_DP, y = K_diag)) +
        geom_point(alpha = 0.6) +
        labs(title = "KGD-like GRM Diagnostic: Diagonal vs Mean DP", x = "Mean DP", y = "GRM diagonal") + theme_diag
      ggsave(file.path(out_dir, "15_kgd_grm_diag_vs_depth.png"), p_k_diag, width = 8, height = 6, dpi = 150)
      cat("  ✓ 15_kgd_grm_diag_vs_depth.png\n")
    }
    
    rm(gt_sub, miss_sub, g, X, I_called, Z, Num, Den, K)
    gc(verbose = FALSE)
  } else {
    warning("KGD_GRM_DIAG: too few samples after filter; skipping.")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# v2.5 Section F: PCA Diagnostic
# ══════════════════════════════════════════════════════════════════════════════
pca_scores <- NULL

if (RUN_PCA_DIAG && n_biallelic_snps > 500) {
  cat("\n── Running PCA diagnostic ──────────────────────────────────────\n")
  
  # Reuse GRM if available, otherwise compute a dosage-based PCA
  if (!is.null(grm_K) && !is.null(grm_keep_names)) {
    cat("  Reusing GRM for PCA (eigen-decomposition)...\n")
    eig <- eigen(grm_K, symmetric = TRUE)
    n_pcs <- min(20, ncol(grm_K))
    pca_scores <- data.frame(
      Sample = grm_keep_names,
      eig$vectors[, 1:n_pcs, drop = FALSE]
    )
    colnames(pca_scores) <- c("Sample", paste0("PC", 1:n_pcs))
    var_explained <- eig$values[1:n_pcs] / sum(eig$values[eig$values > 0]) * 100
  } else {
    cat("  Computing dosage-based PCA (GRM not available)...\n")
    set.seed(PCA_SEED)
    pca_samp <- sample_names[sample_miss_rate <= KGD_SAMPLE_MISS_MAX]
    if (!is.null(sample_mean_dp) && sum(!is.na(sample_mean_dp)) > 0) {
      pca_samp <- pca_samp[!is.na(sample_mean_dp[match(pca_samp, sample_names)]) &
                             sample_mean_dp[match(pca_samp, sample_names)] >= KGD_SAMPLE_DP_MIN]
    }
    if (length(pca_samp) >= 50) {
      pca_snps <- biallelic_snp_idx
      if (is.finite(PCA_MAX_SNPS) && length(pca_snps) > PCA_MAX_SNPS) {
        pca_snps <- sample(pca_snps, PCA_MAX_SNPS)
      }
      gt_pca <- gt_matrix[pca_snps, pca_samp, drop = FALSE]
      g_pca <- matrix(NA_real_, nrow = nrow(gt_pca), ncol = ncol(gt_pca), dimnames = dimnames(gt_pca))
      g_pca[gt_pca == "0/0" | gt_pca == "0|0"] <- 0
      g_pca[gt_pca == "0/1" | gt_pca == "1/0" | gt_pca == "0|1" | gt_pca == "1|0"] <- 1
      g_pca[gt_pca == "1/1" | gt_pca == "1|1"] <- 2
      
      p_pca <- rowSums(g_pca, na.rm = TRUE) / pmax(2 * rowSums(!is.na(g_pca)), 1)
      X_pca <- sweep(g_pca, 1, 2 * p_pca, FUN = "-"); X_pca[is.na(X_pca)] <- 0
      w_pca <- 2 * p_pca * (1 - p_pca); w_pca[!is.finite(w_pca)] <- 0
      I_pca <- !is.na(g_pca)
      Num_pca <- crossprod(X_pca)
      Z_pca <- sweep(I_pca, 1, sqrt(w_pca), `*`)
      Den_pca <- crossprod(Z_pca)
      K_pca <- Num_pca / pmax(Den_pca, 1e-8)
      
      eig <- eigen(K_pca, symmetric = TRUE)
      n_pcs <- min(20, ncol(K_pca))
      pca_scores <- data.frame(Sample = pca_samp, eig$vectors[, 1:n_pcs, drop = FALSE])
      colnames(pca_scores) <- c("Sample", paste0("PC", 1:n_pcs))
      var_explained <- eig$values[1:n_pcs] / sum(eig$values[eig$values > 0]) * 100
      grm_keep_names <- pca_samp
      
      rm(gt_pca, g_pca, X_pca, I_pca, Z_pca, Num_pca, Den_pca, K_pca)
      gc(verbose = FALSE)
    } else {
      cat("  Too few samples after filtering; skipping PCA.\n")
    }
  }
  
  if (!is.null(pca_scores)) {
    cat(sprintf("  PCA computed: %d samples, %d PCs\n", nrow(pca_scores), ncol(pca_scores) - 1))
    cat(sprintf("  Variance explained: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%\n",
                var_explained[1], var_explained[2], var_explained[3]))
    
    # Scree plot
    df_scree <- data.frame(PC = 1:min(10, length(var_explained)),
                           VarExp = var_explained[1:min(10, length(var_explained))])
    p_scree <- ggplot(df_scree, aes(x = PC, y = VarExp)) +
      geom_col(fill = "#2171b5") + geom_line(linewidth = 0.7) + geom_point(size = 2) +
      scale_x_continuous(breaks = df_scree$PC) +
      labs(title = "PCA Scree Plot", subtitle = "Proportion of variance explained per PC",
           x = "Principal Component", y = "% Variance Explained") + theme_diag
    ggsave(file.path(out_dir, "21_pca_scree.png"), p_scree, width = 8, height = 5, dpi = 150)
    cat("  ✓ 21_pca_scree.png\n")
    
    # Merge metadata for coloring
    pca_meta <- pca_scores
    pca_meta$Seq_Plate <- sub(".*_(P\\d+)_.*", "\\1", pca_meta$Sample)
    pca_meta$Missing_rate <- sample_miss_rate[match(pca_meta$Sample, sample_names)]
    if (!is.null(sample_mean_dp)) {
      pca_meta$Mean_DP <- sample_mean_dp[match(pca_meta$Sample, sample_names)]
    }
    
    pc1_lab <- sprintf("PC1 (%.1f%%)", var_explained[1])
    pc2_lab <- sprintf("PC2 (%.1f%%)", var_explained[2])
    
    # PCA by seq plate
    p_pca_plate <- ggplot(pca_meta, aes(x = PC1, y = PC2, color = Seq_Plate)) +
      geom_point(alpha = 0.6, size = 1.5) +
      labs(title = "PCA: PC1 vs PC2 by Sequencing Plate", x = pc1_lab, y = pc2_lab, color = "Seq Plate") +
      theme_diag + theme(legend.position = "right")
    ggsave(file.path(out_dir, "22_pca_by_seq_plate.png"), p_pca_plate, width = 10, height = 7, dpi = 150)
    cat("  ✓ 22_pca_by_seq_plate.png\n")
    
    # PCA by missingness
    p_pca_miss <- ggplot(pca_meta, aes(x = PC1, y = PC2, color = Missing_rate)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_viridis_c() +
      labs(title = "PCA: PC1 vs PC2 by Missingness", x = pc1_lab, y = pc2_lab, color = "Miss Rate") +
      theme_diag
    ggsave(file.path(out_dir, "25_pca_by_missingness.png"), p_pca_miss, width = 10, height = 7, dpi = 150)
    cat("  ✓ 25_pca_by_missingness.png\n")
    
    # PCA by mean DP
    if (!is.null(sample_mean_dp)) {
      p_pca_dp <- ggplot(pca_meta, aes(x = PC1, y = PC2, color = Mean_DP)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_viridis_c(option = "magma") +
        labs(title = "PCA: PC1 vs PC2 by Mean Depth", x = pc1_lab, y = pc2_lab, color = "Mean DP") +
        theme_diag
      ggsave(file.path(out_dir, "25b_pca_by_depth.png"), p_pca_dp, width = 10, height = 7, dpi = 150)
      cat("  ✓ 25b_pca_by_depth.png\n")
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# v2.5 Sections A-E: Plate Layout Diagnostics
# ══════════════════════════════════════════════════════════════════════════════
plate_info <- NULL

if (RUN_PLATE_DIAG && file.exists(plate_layout_path)) {
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("  SECTION 18: PLATE LAYOUT DIAGNOSTICS\n")
  cat("══════════════════════════════════════════════════════════════════\n\n")
  
  # --- Section A: Read & merge plate layout ---
  # Try reading as CSV first; if XLS, convert via readxl or fallback
  if (grepl("\\.csv$", plate_layout_path, ignore.case = TRUE)) {
    plate_raw <- read.csv(plate_layout_path, stringsAsFactors = FALSE, skip = 1)
  } else {
    # Try readxl for XLS/XLSX
    if (requireNamespace("readxl", quietly = TRUE)) {
      plate_raw <- as.data.frame(readxl::read_excel(plate_layout_path, skip = 1), stringsAsFactors = FALSE)
    } else {
      cat("  readxl not available; attempting CSV conversion.\n")
      csv_tmp <- tempfile(fileext = ".csv")
      # Fallback: try system libreoffice conversion
      sys_result <- tryCatch(system2("libreoffice",
                                     args = c("--headless", "--convert-to", "csv", "--outdir", dirname(csv_tmp), plate_layout_path),
                                     stdout = TRUE, stderr = TRUE), error = function(e) NULL)
      csv_converted <- sub("\\.[^.]+$", ".csv", basename(plate_layout_path))
      csv_path_try <- file.path(dirname(csv_tmp), csv_converted)
      if (file.exists(csv_path_try)) {
        plate_raw <- read.csv(csv_path_try, stringsAsFactors = FALSE, skip = 1)
      } else {
        cat("  WARNING: Could not read plate layout file. Skipping plate diagnostics.\n")
        RUN_PLATE_DIAG <- FALSE
        plate_raw <- NULL
      }
    }
  }
  
  if (RUN_PLATE_DIAG && !is.null(plate_raw)) {
    # Identify required columns
    rg_col <- grep("RG_Sample_Code", colnames(plate_raw), value = TRUE)[1]
    uc_col <- grep("Unique_Client_Code", colnames(plate_raw), value = TRUE)[1]
    conc_col <- grep("Concentration", colnames(plate_raw), value = TRUE)[1]
    vol_col  <- grep("Volume", colnames(plate_raw), value = TRUE)[1]
    total_col <- grep("^Total", colnames(plate_raw), value = TRUE)[1]
    
    if (is.na(rg_col) || is.na(uc_col)) {
      cat("  WARNING: RG_Sample_Code or Unique_Client_Code not found. Skipping.\n")
      RUN_PLATE_DIAG <- FALSE
    }
  }
  
  if (RUN_PLATE_DIAG && !is.null(plate_raw)) {
    plate_info <- data.frame(
      RG_Sample_Code = plate_raw[[rg_col]],
      Unique_Client_Code = plate_raw[[uc_col]],
      stringsAsFactors = FALSE
    )
    if (!is.na(conc_col))  plate_info$DNA_Concentration <- suppressWarnings(as.numeric(plate_raw[[conc_col]]))
    if (!is.na(vol_col))   plate_info$DNA_Volume <- suppressWarnings(as.numeric(plate_raw[[vol_col]]))
    if (!is.na(total_col)) plate_info$Total_DNA <- suppressWarnings(as.numeric(plate_raw[[total_col]]))
    
    # Parse seq plate from RG_Sample_Code: UBC_071001_P01_WA01 -> P01
    plate_info$Seq_Plate <- sub(".*_(P\\d+)_.*", "\\1", plate_info$RG_Sample_Code)
    
    # Parse well from RG_Sample_Code: ...W[A-H][01-12] -> row, col
    well_str <- sub(".*_W([A-H])(\\d+)$", "\\1_\\2", plate_info$RG_Sample_Code)
    plate_info$Well_Row <- sub("_.*", "", well_str)
    plate_info$Well_Col <- suppressWarnings(as.integer(sub(".*_", "", well_str)))
    
    # Parse from Unique_Client_Code: P{n}S{s}_{well}_{individual}
    plate_info$Client_Plate <- sub("^(P\\d+)S.*", "\\1", plate_info$Unique_Client_Code)
    plate_info$S_Group <- sub("^P\\d+(S[\\d&]+)_.*", "\\1", plate_info$Unique_Client_Code, perl = TRUE)
    plate_info$Individual_ID <- suppressWarnings(as.integer(sub(".*_(\\d+)$", "\\1", plate_info$Unique_Client_Code)))
    plate_info$Is_Empty <- !is.na(plate_info$Individual_ID) & plate_info$Individual_ID == 0
    
    # Merge with VCF sample names
    in_vcf <- plate_info$RG_Sample_Code %in% sample_names
    in_plate <- sample_names %in% plate_info$RG_Sample_Code
    
    cat(sprintf("  Plate layout samples      : %d\n", nrow(plate_info)))
    cat(sprintf("  Matched to VCF            : %d\n", sum(in_vcf)))
    cat(sprintf("  VCF samples not in layout : %d\n", sum(!in_plate)))
    cat(sprintf("  Layout samples not in VCF : %d\n", sum(!in_vcf)))
    
    # --- Section B: Empty well flagging ---
    empty_wells <- plate_info[plate_info$Is_Empty & in_vcf, ]
    n_empties_in_vcf <- nrow(empty_wells)
    
    cat(sprintf("\n  Empty wells (Individual_ID = 0): %d total, %d in VCF\n",
                sum(plate_info$Is_Empty), n_empties_in_vcf))
    
    if (n_empties_in_vcf > 0) {
      empty_miss <- sample_miss_rate[match(empty_wells$RG_Sample_Code, sample_names)]
      empty_dp   <- if (!is.null(sample_mean_dp)) sample_mean_dp[match(empty_wells$RG_Sample_Code, sample_names)] else rep(NA, n_empties_in_vcf)
      empty_call_rate <- 1 - empty_miss
      
      empty_report <- data.frame(
        RG_Sample_Code = empty_wells$RG_Sample_Code,
        Seq_Plate = empty_wells$Seq_Plate,
        Missing_rate = round(empty_miss, 4),
        Call_rate = round(empty_call_rate, 4),
        Mean_DP = round(empty_dp, 2),
        stringsAsFactors = FALSE
      )
      write.csv(empty_report, file.path(out_dir, "empty_well_report.csv"), row.names = FALSE)
      cat(sprintf("  ✓ empty_well_report.csv saved\n"))
      
      # Warn if any empties have meaningful call rates
      suspicious <- sum(empty_call_rate > 0.05, na.rm = TRUE)
      if (suspicious > 0) {
        cat(sprintf("  ⚠ WARNING: %d empty wells have call rate > 5%% — potential contamination!\n", suspicious))
      }
      
      cat(sprintf("  Empty well missingness: mean=%.4f, min=%.4f, max=%.4f\n",
                  safe_mean(empty_miss), safe_min(empty_miss), safe_max(empty_miss)))
    }
    
    # --- Section C: Plate batch-effect diagnostics ---
    # Merge QC metrics for samples in VCF
    plate_qc <- plate_info[in_vcf, ]
    plate_qc$Missing_rate <- sample_miss_rate[match(plate_qc$RG_Sample_Code, sample_names)]
    plate_qc$Het_rate     <- sample_het_rate_bial[match(plate_qc$RG_Sample_Code, sample_names)]
    if (!is.null(sample_mean_dp)) {
      plate_qc$Mean_DP <- sample_mean_dp[match(plate_qc$RG_Sample_Code, sample_names)]
    }
    
    cat("\n  ── Batch effect diagnostics ──\n")
    
    # Boxplot: missingness by seq plate
    p_batch_miss <- ggplot(plate_qc, aes(x = Seq_Plate, y = Missing_rate)) +
      geom_boxplot(outlier.size = 0.5) +
      labs(title = "Per-Sample Missingness by Sequencing Plate", x = "Seq Plate", y = "Missing Rate") +
      theme_diag + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "16_batch_boxplot_missingness.png"), p_batch_miss, width = 10, height = 5, dpi = 150)
    cat("  ✓ 16_batch_boxplot_missingness.png\n")
    
    # Boxplot: depth by seq plate
    if (!is.null(sample_mean_dp)) {
      p_batch_dp <- ggplot(plate_qc, aes(x = Seq_Plate, y = Mean_DP)) +
        geom_boxplot(outlier.size = 0.5) +
        labs(title = "Per-Sample Mean Depth by Sequencing Plate", x = "Seq Plate", y = "Mean DP") +
        theme_diag + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(out_dir, "17_batch_boxplot_depth.png"), p_batch_dp, width = 10, height = 5, dpi = 150)
      cat("  ✓ 17_batch_boxplot_depth.png\n")
    }
    
    # Boxplot: het by seq plate
    p_batch_het <- ggplot(plate_qc, aes(x = Seq_Plate, y = Het_rate)) +
      geom_boxplot(outlier.size = 0.5) +
      labs(title = "Per-Sample Heterozygosity by Sequencing Plate", x = "Seq Plate", y = "Het Rate") +
      theme_diag + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "18_batch_boxplot_het.png"), p_batch_het, width = 10, height = 5, dpi = 150)
    cat("  ✓ 18_batch_boxplot_het.png\n")
    
    # Boxplot: by S group
    p_batch_miss_s <- ggplot(plate_qc, aes(x = S_Group, y = Missing_rate)) +
      geom_boxplot(outlier.size = 0.5) +
      labs(title = "Per-Sample Missingness by S Group", x = "S Group", y = "Missing Rate") +
      theme_diag + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "16b_batch_boxplot_missingness_sgroup.png"), p_batch_miss_s, width = 8, height = 5, dpi = 150)
    cat("  ✓ 16b_batch_boxplot_missingness_sgroup.png\n")
    
    # Kruskal-Wallis tests
    kw_miss_plate <- tryCatch(kruskal.test(Missing_rate ~ Seq_Plate, data = plate_qc)$p.value, error = function(e) NA)
    kw_het_plate  <- tryCatch(kruskal.test(Het_rate ~ Seq_Plate, data = plate_qc)$p.value, error = function(e) NA)
    kw_miss_s     <- tryCatch(kruskal.test(Missing_rate ~ S_Group, data = plate_qc)$p.value, error = function(e) NA)
    cat(sprintf("\n  Kruskal-Wallis p-values (batch effect tests):\n"))
    cat(sprintf("    Missingness ~ Seq Plate : %s\n", if (is.na(kw_miss_plate)) "NA" else format.pval(kw_miss_plate, digits = 3)))
    cat(sprintf("    Het rate ~ Seq Plate    : %s\n", if (is.na(kw_het_plate)) "NA" else format.pval(kw_het_plate, digits = 3)))
    cat(sprintf("    Missingness ~ S Group   : %s\n", if (is.na(kw_miss_s)) "NA" else format.pval(kw_miss_s, digits = 3)))
    
    if (!is.null(sample_mean_dp)) {
      kw_dp_plate <- tryCatch(kruskal.test(Mean_DP ~ Seq_Plate, data = plate_qc)$p.value, error = function(e) NA)
      cat(sprintf("    Mean DP ~ Seq Plate     : %s\n", if (is.na(kw_dp_plate)) "NA" else format.pval(kw_dp_plate, digits = 3)))
    }
    
    # --- Section D: DNA input quality vs QC metrics ---
    if (!is.na(conc_col) && "DNA_Concentration" %in% colnames(plate_qc)) {
      plate_qc_dna <- plate_qc[!is.na(plate_qc$DNA_Concentration) & plate_qc$DNA_Concentration > 0, ]
      if (nrow(plate_qc_dna) > 10) {
        p_dna1 <- ggplot(plate_qc_dna, aes(x = DNA_Concentration, y = Missing_rate)) +
          geom_point(alpha = 0.4, size = 1) + geom_smooth(method = "loess", se = FALSE, color = "red") +
          labs(title = "DNA Concentration vs Missingness", x = "Concentration (ng/uL)", y = "Missing Rate") + theme_diag
        
        plots_dna <- list(p_dna1)
        if (!is.null(sample_mean_dp) && "Mean_DP" %in% colnames(plate_qc_dna)) {
          p_dna2 <- ggplot(plate_qc_dna, aes(x = DNA_Concentration, y = Mean_DP)) +
            geom_point(alpha = 0.4, size = 1) + geom_smooth(method = "loess", se = FALSE, color = "red") +
            labs(title = "DNA Concentration vs Mean Depth", x = "Concentration (ng/uL)", y = "Mean DP") + theme_diag
          plots_dna <- list(p_dna1, p_dna2)
        }
        if ("Total_DNA" %in% colnames(plate_qc_dna)) {
          plate_qc_dna2 <- plate_qc_dna[!is.na(plate_qc_dna$Total_DNA) & plate_qc_dna$Total_DNA > 0, ]
          if (nrow(plate_qc_dna2) > 10) {
            p_dna3 <- ggplot(plate_qc_dna2, aes(x = Total_DNA, y = Missing_rate)) +
              geom_point(alpha = 0.4, size = 1) + geom_smooth(method = "loess", se = FALSE, color = "red") +
              labs(title = "Total DNA vs Missingness", x = "Total DNA (ng)", y = "Missing Rate") + theme_diag
            plots_dna <- c(plots_dna, list(p_dna3))
          }
        }
        p_dna_panel <- arrangeGrob(grobs = plots_dna, ncol = min(length(plots_dna), 2))
        ggsave(file.path(out_dir, "19_dna_input_vs_qc.png"), p_dna_panel,
               width = 7 * min(length(plots_dna), 2), height = 6, dpi = 150)
        cat("  ✓ 19_dna_input_vs_qc.png\n")
      }
    }
    
    # --- Section E: Well-position effects ---
    plate_qc_well <- plate_qc[!is.na(plate_qc$Well_Row) & !is.na(plate_qc$Well_Col), ]
    if (nrow(plate_qc_well) > 50) {
      well_agg <- aggregate(Missing_rate ~ Well_Row + Well_Col, data = plate_qc_well, FUN = mean, na.rm = TRUE)
      
      p_well_miss <- ggplot(well_agg, aes(x = factor(Well_Col), y = Well_Row, fill = Missing_rate)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c(option = "inferno", direction = -1) +
        labs(title = "Well-Position Effect: Mean Missingness",
             subtitle = "Averaged across all plates; check edges for spatial artifacts",
             x = "Column", y = "Row", fill = "Mean\nMiss Rate") +
        theme_diag
      
      if (!is.null(sample_mean_dp) && "Mean_DP" %in% colnames(plate_qc_well)) {
        well_agg_dp <- aggregate(Mean_DP ~ Well_Row + Well_Col, data = plate_qc_well, FUN = mean, na.rm = TRUE)
        p_well_dp <- ggplot(well_agg_dp, aes(x = factor(Well_Col), y = Well_Row, fill = Mean_DP)) +
          geom_tile(color = "white") +
          scale_fill_viridis_c(option = "viridis") +
          labs(title = "Well-Position Effect: Mean Depth",
               subtitle = "Averaged across all plates",
               x = "Column", y = "Row", fill = "Mean DP") +
          theme_diag
        p_well_panel <- arrangeGrob(p_well_miss, p_well_dp, ncol = 2)
      } else {
        p_well_panel <- p_well_miss
      }
      ggsave(file.path(out_dir, "20_well_position_heatmap.png"), p_well_panel, width = 14, height = 5, dpi = 150)
      cat("  ✓ 20_well_position_heatmap.png\n")
    }
    
    # PCA colored by S group and client plate (if PCA available)
    if (!is.null(pca_scores)) {
      pca_plate <- merge(pca_scores, plate_info[, c("RG_Sample_Code", "S_Group", "Client_Plate")],
                         by.x = "Sample", by.y = "RG_Sample_Code", all.x = TRUE)
      
      if (sum(!is.na(pca_plate$S_Group)) > 10) {
        p_pca_s <- ggplot(pca_plate, aes(x = PC1, y = PC2, color = S_Group)) +
          geom_point(alpha = 0.6, size = 1.5) +
          labs(title = "PCA: PC1 vs PC2 by S Group",
               x = sprintf("PC1 (%.1f%%)", var_explained[1]),
               y = sprintf("PC2 (%.1f%%)", var_explained[2]), color = "S Group") +
          theme_diag
        ggsave(file.path(out_dir, "23_pca_by_s_group.png"), p_pca_s, width = 10, height = 7, dpi = 150)
        cat("  ✓ 23_pca_by_s_group.png\n")
      }
      
      if (sum(!is.na(pca_plate$Client_Plate)) > 10) {
        p_pca_cp <- ggplot(pca_plate, aes(x = PC1, y = PC2, color = Client_Plate)) +
          geom_point(alpha = 0.6, size = 1.5) +
          labs(title = "PCA: PC1 vs PC2 by Client Plate",
               x = sprintf("PC1 (%.1f%%)", var_explained[1]),
               y = sprintf("PC2 (%.1f%%)", var_explained[2]), color = "Client Plate") +
          theme_diag
        ggsave(file.path(out_dir, "24_pca_by_client_plate.png"), p_pca_cp, width = 10, height = 7, dpi = 150)
        cat("  ✓ 24_pca_by_client_plate.png\n")
      }
    }
  }
} else if (RUN_PLATE_DIAG) {
  cat(sprintf("  WARNING: Plate layout file not found: %s\n", plate_layout_path))
  cat("  Skipping plate layout diagnostics.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# SAVE OUTPUTS (updated with all v2.5 columns)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Saving output tables ────────────────────────────────────────\n")

# --- Update per_sample_stats with v2.5 columns ---
if (RUN_AB_DIAG) {
  sample_stats$Mean_AB_het <- round(sample_mean_AB_het[match(sample_stats$Sample, sample_names)], 6)
  sample_stats$Hom_minor_allele_frac <- round(sample_hom_minor_frac[match(sample_stats$Sample, sample_names)], 6)
}

# Plate layout columns
if (!is.null(plate_info)) {
  pi_match <- match(sample_stats$Sample, plate_info$RG_Sample_Code)
  sample_stats$Seq_Plate    <- plate_info$Seq_Plate[pi_match]
  sample_stats$S_Group      <- plate_info$S_Group[pi_match]
  sample_stats$Client_Plate <- plate_info$Client_Plate[pi_match]
  sample_stats$Well_Row     <- plate_info$Well_Row[pi_match]
  sample_stats$Well_Col     <- plate_info$Well_Col[pi_match]
  sample_stats$Individual_ID <- plate_info$Individual_ID[pi_match]
  sample_stats$Is_Empty     <- plate_info$Is_Empty[pi_match]
  if ("DNA_Concentration" %in% colnames(plate_info)) {
    sample_stats$DNA_Concentration <- plate_info$DNA_Concentration[pi_match]
  }
  if ("DNA_Volume" %in% colnames(plate_info)) {
    sample_stats$DNA_Volume <- plate_info$DNA_Volume[pi_match]
  }
  if ("Total_DNA" %in% colnames(plate_info)) {
    sample_stats$Total_DNA <- plate_info$Total_DNA[pi_match]
  }
}

# PCA scores
if (!is.null(pca_scores)) {
  pc_match <- match(sample_stats$Sample, pca_scores$Sample)
  n_pc_save <- min(5, ncol(pca_scores) - 1)
  for (pc_i in 1:n_pc_save) {
    sample_stats[[paste0("PC", pc_i)]] <- pca_scores[[paste0("PC", pc_i)]][pc_match]
  }
}

sample_stats <- sample_stats[order(-sample_stats$Missing_rate_all), ]
write.csv(sample_stats, file.path(out_dir, "per_sample_stats.csv"), row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", file.path(out_dir, "per_sample_stats.csv")))

# --- Build per_site_stats ---
site_stats <- data.frame(
  CHROM                    = vcf@fix[, "CHROM"],
  POS                      = suppressWarnings(as.integer(vcf@fix[, "POS"])),
  REF                      = ref_alleles,
  ALT                      = alt_alleles,
  QUAL                     = qual_vals,
  Type                     = variant_types,
  Call_rate_all            = round(1 - site_miss_rate, 6),
  MAF_biallelic_snp        = round(site_maf, 6),
  Het_rate_biallelic_snp   = round(site_het_rate_bial, 6),
  H_obs_biallelic_snp      = round(H_obs, 6),
  H_exp_biallelic_snp      = round(H_exp, 6),
  fin_ratio                = round(fin_ratio, 6),
  excess_het_flag          = excess_het_flag,
  Mean_DP                  = if (!is.null(site_mean_dp)) round(site_mean_dp, 2) else NA_real_,
  high_depth_site_flag     = high_depth_site_flag,
  depth_site_q995          = depth_site_q995,
  stringsAsFactors = FALSE
)

# v2.5 site columns
if (RUN_AB_DIAG) {
  site_stats$median_AB <- round(site_median_AB, 6)
  site_stats$AB_flag   <- site_AB_flag
}

if (RUN_SNP_CLUSTER_DIAG && length(site_snp_cluster) > 0) {
  for (col_name in names(site_snp_cluster)) {
    site_stats[[col_name]] <- site_snp_cluster[[col_name]]
  }
}

write.csv(site_stats, file.path(out_dir, "per_site_stats.csv"), row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", file.path(out_dir, "per_site_stats.csv")))

# ── Done ─────────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  ALL OUTPUTS SAVED TO: %s\n", out_dir))
cat("══════════════════════════════════════════════════════════════════\n")
cat("\n  Output files:\n")
cat("    • per_sample_stats.csv            — per-sample missingness, het, depth, AB, plate info, PCs\n")
cat("    • per_site_stats.csv              — per-variant MAF, call rate, het, fin, depth, AB, clustering\n")
cat("    • empty_well_report.csv           — empty well flagging (if plate layout provided)\n")
cat("    • 01_maf_distribution.png         — MAF histogram\n")
cat("    • 02_site_missing_rate.png        — site missing rate histogram\n")
cat("    • 03_sample_missing_rate.png      — sample missing rate histogram\n")
cat("    • 04_sample_heterozygosity.png    — sample het rate histogram\n")
cat("    • 05_sample_het_vs_missing.png    — het vs missing scatterplot\n")
cat("    • 06_depth_distribution.png       — DP histogram\n")
cat("    • 06b/06c_sample_mean_depth_*.png — per-sample mean depth (low/high N)\n")
cat("    • 07_qual_distribution.png        — QUAL score histogram\n")
cat("    • 08_hwe_pvalue_distribution.png  — HWE p-value histogram\n")
cat("    • 09_variants_per_chrom.png       — variants per chromosome\n")
cat("    • 10_diagnostic_panel.png         — 2x2 panel of core diagnostics\n")
cat("    • 11_sample_depth_vs_het.png      — sample mean DP vs heterozygosity\n")
cat("    • 12_sample_depth_vs_missing.png  — sample mean DP vs missingness\n")
cat("    • 13_fin_plot_excess_het.png      — KGD fin plot + flagged excess-het loci\n")
cat("    • 14_kgd_grm_offdiag_hist.png     — KGD-like GRM off-diagonal histogram\n")
cat("    • 15_kgd_grm_diag_vs_depth.png    — KGD-like GRM diagonal vs depth\n")
cat("    • 16_batch_boxplot_missingness.png — missingness by seq plate\n")
cat("    • 17_batch_boxplot_depth.png       — depth by seq plate\n")
cat("    • 18_batch_boxplot_het.png         — het rate by seq plate\n")
cat("    • 19_dna_input_vs_qc.png          — DNA concentration/total vs QC metrics\n")
cat("    • 20_well_position_heatmap.png     — well position spatial effects\n")
cat("    • 21_pca_scree.png                — PCA scree plot\n")
cat("    • 22_pca_by_seq_plate.png         — PCA colored by seq plate\n")
cat("    • 23_pca_by_s_group.png           — PCA colored by S group\n")
cat("    • 24_pca_by_client_plate.png      — PCA colored by client plate\n")
cat("    • 25_pca_by_missingness.png       — PCA colored by missingness\n")
cat("    • 25b_pca_by_depth.png            — PCA colored by mean depth\n")
cat("    • 26_allele_balance_site_hist.png  — per-site median allele balance\n")
cat("    • 27_allele_balance_sample_hist.png — per-sample mean allele balance\n")
cat("    • 28_hom_contamination_hist.png    — contamination check at hom sites\n")
cat("    • 29_folded_sfs.png               — folded site frequency spectrum\n")
cat("    • 30_snp_clustering_*.png          — SNP proximity clustering\n")
cat("    • 31_depth_het_bias.png            — depth-dependent het rate bias\n")
cat("    • kgd_grm_cleaned_samples.rds      — KGD-like GRM object\n")
cat("\n  ✔ Pre-QC diagnostics complete (v2.5).\n\n")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

# Final cleanup
to_rm <- c("vcf","gt_matrix","gt_bial","gt_bial_missing","gt_bial_het_mask",
           "sample_stats","site_stats","chrom_counts",
           "ref_alleles","alt_alleles","ref_uc","alt_uc","variant_types","filter_vals",
           "qual_vals","site_maf","site_alt_freq","site_miss_rate","site_called",
           "site_het_rate_bial","site_called_bial","site_het_bial","sample_missing",
           "sample_miss_rate","sample_called_all","sample_called_bial",
           "sample_het_count_bial","sample_het_rate_bial","sample_homalt",
           "sample_homalt_rate","sample_names","valid_pvals","site_mean_dp",
           "sample_mean_dp",
           "H_obs","H_exp","fin_ratio","excess_het_flag",
           "high_depth_site_flag","depth_site_q995",
           "dp_matrix","dp_flat","dp_df","dp_df_low","dp_df_high",
           "df_depth_sample","df_fin",
           "bial_n_00","bial_n_01","bial_n_11","bial_n_total",
           "bial_allele_count","bial_alt_allele_sum",
           "site_median_AB","site_AB_flag","sample_mean_AB_het","sample_hom_minor_frac",
           "sfs_mac_vec","site_snp_cluster","dp_het_table",
           "grm_K","grm_keep_names","pca_scores","plate_info","plate_raw",
           "p1","p2","p3","p4","p5","p6","p6b","p6c","p7","p8","p9","p_fin","p_panel",
           "p_dp_het","p_dp_miss","p_k_hist","p_k_diag","df_k","df_diag",
           "p_ab_site","p_ab_samp","p_contam","p_sfs","p_clust","p_dh",
           "p_scree","p_pca_plate","p_pca_miss","p_pca_dp","p_pca_s","p_pca_cp",
           "p_batch_miss","p_batch_dp","p_batch_het","p_batch_miss_s",
           "p_well_miss","p_well_dp","p_well_panel","p_dna_panel",
           "empty_report","plate_qc","plate_qc_dna","plate_qc_well",
           "contam_flagged","var_explained","pca_meta","pca_plate","eig")
rm(list = intersect(to_rm, ls()))
gc(verbose = FALSE)