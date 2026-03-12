# =============================================================================
# snp_qc_pipeline_v4_6_1_patched.R
#
# SNP quality control pipeline for the Douglas-fir genomic selection dataset.
# Reads from SeqArray GDS; all downstream QC logic identical to v4_5 except
# where noted below.
#
# ------------------------------------------------------------------------------
# CHANGES FROM v4_6
# ------------------------------------------------------------------------------
# v4.6.1 patch (4 fixes):
#   1. Marginal SNP mask: removed !snp_fail_hwe_excess condition so that
#      the Mendelian-error comparison actually populates in datasets where
#      the HWE test dominates the ratio filter.
#   2. Het-sensitivity table: now shows ratio-only, HWE-only, and combined
#      fail counts separately at each threshold so the user can see which
#      filter is driving removals.
#   3. Dosage variance ratio warning: gated on MAF >= 0.10; low-MAF
#      compression is expected dosage mechanics, not paralog signal.
#      MAF-stratified DVR summary added.
#   4. Ts/Tv warning: softened for non-WGS / transcriptome-reference
#      contexts where lower ratios are expected.
#
# ------------------------------------------------------------------------------
# CHANGES FROM v4_5
# ------------------------------------------------------------------------------
# 1. Threshold relaxation (motivated by v4_5 recovering only ~12K SNPs):
#    • SNP_CR: 0.90 → 0.70  (primary lever for SNP recovery)
#    • HET_EXCESS_RATIO_MAX: 1.50 → 1.25  (tighter paralog filter to
#      compensate for the relaxed call-rate gate — cleaner survivors)
#    • SAMPLE_CR_FAIL: 0.80 → 0.70  (modest; keeps borderline samples)
#    • SAMPLE_CR_WARN: 0.90 → 0.80
#    • PED_HET_EXCESS_RATIO_MAX: 1.25 → 1.15  (stricter pedigree panel)
#
# 2. Optional dosage matrix (USE_DOSAGE toggle):
#    • When enabled, Stage A extracts freeBayes FORMAT/RO (ref obs count)
#      and FORMAT/AO (alt obs count) from the GDS, and computes
#      continuous allele dosage = 2 × AO / (RO + AO).
#    • The dosage matrix is carried in parallel through all sample and
#      SNP filtering stages (same row/column subsetting as gt_012).
#    • Saved as snps_qc_pass_dosage.rds alongside the 012 matrix.
#    • All QC decisions (call rate, MAF, het, HWE, pedigree checks,
#      Mendelian errors, IBS0, duplicates) continue to use the discrete
#      gt_012 matrix — dosage is purely an output for downstream GRM
#      construction and genomic prediction.
#    • Falls back gracefully if RO/AO fields are absent from the GDS.
#
# 3. Log flushing: explicit flush(log_con) after each stage header to
#    prevent truncated log files on Windows.
#
# 4. Conifer-specific diagnostics:
#    • Het-ratio sensitivity analysis (Stage D): logs SNP counts at
#      alternative thresholds (1.20–1.50) so the user can evaluate
#      whether 1.25 is too aggressive for an obligate outcrosser with
#      high genetic load.
#    • Mendelian error stratified by het-ratio bin (section 9b2):
#      computes per-SNP Mendelian error rates for "marginal" SNPs
#      (ratio 1.25–1.40, removed by Stage D) and compares to pass-QC
#      SNPs. If marginal SNPs don't have elevated errors, the excess
#      het is biological and the threshold should be relaxed.
#    • Per-sample inbreeding coefficient F_hat (Stage E):
#      F_hat = 1 - obs_het / E[het], where E[het] accounts for which
#      loci each sample has calls at. Flags potential selfing,
#      contamination, or systematic reference bias.
#    • Dosage variance ratio (section 6 + 9d): compares observed
#      per-SNP dosage variance to HWE expectation (2p(1-p)). Ratio
#      < 1 flags collapsed paralogs; reported in snp_qc_summary.csv
#      and as a binned distribution in the diagnostic summary.
#    • Transition / transversion ratio (section 10): Ts/Tv as a
#      global quality check on variant calls.
#
# 5. MAF_GS relaxed from 0.05 to 0.01 so the GS panel includes rare
#    variants — in conifers with rapid LD decay, these contribute to
#    realized-relationship estimation in the GRM.
#
# ------------------------------------------------------------------------------
# OVERVIEW  (unchanged from v4_5)
# ------------------------------------------------------------------------------
# Staged filtering: A (genotype masking) → B (dead samples) → C (site QC) →
# D (excess het / paralog, G1 only) → E (refined sample QC) → F (MAF).
# Post-QC: pedigree panel, IBS0 verification, Mendelian error check,
# duplicate detection, diagnostic summary.
#
# ------------------------------------------------------------------------------
# REQUIRED PACKAGES
# ------------------------------------------------------------------------------
#   install.packages(c("dplyr", "readr", "stringr", "ggplot2"))
#   BiocManager::install("SeqArray")
#
# ------------------------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------------------------
#   UBC_071001_snps_RAW.full.seq.gds  — GDS converted from freeBayes VCF
#                                        (via seqVCF2GDS with fmt.import=NULL
#                                        to preserve DP/GQ/RO/AO FORMAT fields)
#   pedigree_full.csv
#   pedigree_genotyped.csv
#   All.sites.csv
#
# ------------------------------------------------------------------------------
# VERSION HISTORY
# ------------------------------------------------------------------------------
# v4_6_1 — Four targeted patches; see changelog above.
# v4_6 — Relaxed thresholds + optional dosage matrix + conifer diagnostics.
# v4_5 — Data source changed from VCF to GDS (SeqArray). QC logic from v4_4.
# v4_4 — Five new filters, parent search, Mendelian check, duplicate detection.
# v4_3 — Staged filtering restructure (A–F).
# v4_2 — Contig sanity check, scaffold collapsing.
# =============================================================================

library(SeqArray)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

stopifnot("R >= 4.0.0 required" = getRversion() >= "4.0.0")

# -----------------------------------------------------------------------------
# 0. Paths and logging
# -----------------------------------------------------------------------------
DATA_DIR       <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data"
GDS_FILE       <- "C:/Users/bratclif/dfir_gds/UBC_071001_snps_RAW.full.seq.gds"
PED_FILE       <- file.path(DATA_DIR, "pedigree_full.csv")
PED_GENO_FILE  <- file.path(DATA_DIR, "pedigree_genotyped.csv")
ALL_SITES_FILE <- file.path(DATA_DIR, "All.sites.csv")
OUT_DIR        <- file.path(DATA_DIR, "qc")
POSTQC_GDS_FILE <- file.path(OUT_DIR, "doug_fir_postQC.gds")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# LOG_FILE <- file.path(OUT_DIR, "qc_report.txt")
# log_con  <- file(LOG_FILE, open = "wt")
# options(warn = 1)
# 
# on.exit({ try(close(log_con), silent = TRUE) }, add = TRUE)

# Override cat() to write to both console and log file
# cat <- function(..., file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE) {
#   base::cat(..., file = stderr(), sep = sep, fill = fill, labels = labels, append = append)
#   base::cat(..., file = log_con, sep = sep, fill = fill, labels = labels, append = FALSE)
#   invisible(try(flush(log_con), silent = TRUE))
# }
# 
# flush_log <- function() invisible(try(flush(log_con), silent = TRUE))

# -----------------------------------------------------------------------------
# QC thresholds — staged filtering
# -----------------------------------------------------------------------------
# ── v4_6 relaxations (commented values = v4_5 originals) ──

# Stage A: per-genotype masking  (unchanged)
DP_MIN                    <- 3
DP_MAX_MULT               <- 2
GQ_MIN                    <- 0
SNP_MEAN_DP_MAX           <- 0

# Stage B: coarse sample pre-filter  (unchanged)
SAMPLE_CR_DEAD            <- 0.50
SAMPLE_DP_MIN             <- 5

# Stage C: site filtering pass 1
SNP_CR                    <- 0.70     # was 0.90 — primary lever for SNP recovery
QUAL_MIN                  <- 20
SNP_CLUSTER_BP            <- 0
SNP_CLUSTER_MAX           <- 3

# Stage D: excess het / paralog filter (G1 only)
HET_EXCESS_RATIO_MAX      <- 1.25    # was 1.50 — tighter to compensate for relaxed CR
HEXP_MIN                  <- 0.05
HWE_EXCESS_PVAL           <- 1e-5

# Stage E: refined sample QC
SAMPLE_CR_WARN            <- 0.80    # was 0.90
SAMPLE_CR_FAIL            <- 0.70    # was 0.80
HET_SD                    <- 3

# Stage F: MAF filter  (unchanged)
MAF_MIN                   <- 0.01
MAF_GS                    <- 0.01   # was 0.05 — include rare variants for GRM

# Pedigree checks  (unchanged)
IBS0_PASS                 <- 0.005
IBS0_WARN                 <- 0.020
MIN_IBS_MARKERS           <- 100
PARENT_SEARCH_TOP_N       <- 5

# Mendelian error check  (unchanged)
MENDEL_SNP_ERROR_THR      <- 0.05

# Duplicate / identity detection  (unchanged)
DUPLICATE_CONC_THR        <- 0.99
DUPLICATE_MAX_SNPS        <- 5000

# Pedigree-panel thresholds
PED_SNP_CR                <- 0.98
PED_MAF_MIN               <- 0.05
PED_HET_EXCESS_RATIO_MAX  <- 1.15   # was 1.25 — stricter pedigree panel
PED_CR_DIFF_MAX           <- 0.05
PED_PHYSICAL_THIN_BP      <- 50000

# ── Dosage mode toggle ──
USE_DOSAGE                <- TRUE

# -----------------------------------------------------------------------------
# Control samples to exclude from QC entirely
# -----------------------------------------------------------------------------
CONTROL_SAMPLE_IDS <- c(
  "UBC_071001_P23_WE04",
  "UBC_071001_P23_WE02",
  "UBC_071001_P23_WE07",
  "UBC_071001_P23_WE10",
  "UBC_071001_P23_WE08",
  "UBC_071001_P23_WE06",
  "UBC_071001_P23_WE05",
  "UBC_071001_P23_WE09",
  "UBC_071001_P23_WE01",
  "UBC_071001_P23_WE03"
)

# -----------------------------------------------------------------------------
# Helper functions  (identical to v4_5)
# -----------------------------------------------------------------------------
snp_call_rate <- function(mat) {
  if (nrow(mat) == 0L) return(rep(NA_real_, ncol(mat)))
  colMeans(!is.na(mat))
}

sample_call_rate <- function(mat) {
  if (ncol(mat) == 0L) return(rep(NA_real_, nrow(mat)))
  rowMeans(!is.na(mat))
}

sample_het <- function(mat) {
  result <- rowMeans(mat == 1L, na.rm = TRUE)
  result[rowSums(!is.na(mat)) == 0L] <- NA_real_
  result
}

snp_maf <- function(mat) {
  freq <- colMeans(mat, na.rm = TRUE) / 2
  freq[is.nan(freq)] <- NA_real_
  pmin(freq, 1 - freq)
}

snp_obs_het <- function(mat) {
  result <- colMeans(mat == 1L, na.rm = TRUE)
  result[is.nan(result)] <- NA_real_
  result
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

ibs0_rate <- function(g1, g2, min_markers = MIN_IBS_MARKERS) {
  keep <- !is.na(g1) & !is.na(g2)
  if (sum(keep) < min_markers) return(NA_real_)
  mean((g1[keep] == 0L & g2[keep] == 2L) | (g1[keep] == 2L & g2[keep] == 0L))
}

n_compared_markers <- function(g1, g2) {
  sum(!is.na(g1) & !is.na(g2))
}

flag_snp_clusters <- function(chrom, pos, window_bp, max_in_window) {
  n <- length(chrom)
  flagged <- logical(n)
  if (n == 0L) return(flagged)
  
  ord <- order(chrom, pos, na.last = TRUE)
  chrom_s <- chrom[ord]
  pos_s   <- pos[ord]
  
  i <- 1L
  while (i <= n) {
    if (is.na(pos_s[i])) { i <- i + 1L; next }
    j <- i
    while (j < n && identical(chrom_s[j + 1L], chrom_s[i]) &&
           !is.na(pos_s[j + 1L]) && (pos_s[j + 1L] - pos_s[i]) <= window_bp) {
      j <- j + 1L
    }
    cluster_size <- j - i + 1L
    if (cluster_size > max_in_window) {
      flagged[ord[i:j]] <- TRUE
    }
    i <- i + 1L
  }
  flagged
}

physical_thin_map <- function(map_df, min_bp) {
  if (nrow(map_df) == 0L || min_bp <= 0L) {
    return(map_df$snp_id)
  }
  map_ord <- map_df %>% arrange(chromosome, position, snp_id)
  keep <- logical(nrow(map_ord))
  last_chr <- NULL
  last_pos <- NA_integer_
  
  for (i in seq_len(nrow(map_ord))) {
    chr_i <- map_ord$chromosome[i]
    pos_i <- map_ord$position[i]
    if (is.na(pos_i)) {
      keep[i] <- FALSE
    } else if (!identical(chr_i, last_chr) || is.na(last_pos) || (pos_i - last_pos) >= min_bp) {
      keep[i] <- TRUE
      last_chr <- chr_i
      last_pos <- pos_i
    }
  }
  
  map_ord$snp_id[keep]
}

mendelian_errors_trio <- function(g_sire, g_dam, g_off) {
  n <- length(g_sire)
  err <- rep(NA, n)
  called <- !is.na(g_sire) & !is.na(g_dam) & !is.na(g_off)
  if (sum(called) == 0L) return(err)
  
  s <- g_sire[called]; d <- g_dam[called]; o <- g_off[called]
  min_off <- (s == 2L) + (d == 2L)
  max_off <- (s > 0L)  + (d > 0L)
  
  err[called] <- (o < min_off) | (o > max_off)
  err
}

pairwise_concordance <- function(gmat) {
  n <- nrow(gmat)
  called <- !is.na(gmat)
  n_joint <- tcrossprod(called * 1.0)
  
  n_match <- matrix(0.0, n, n)
  for (v in 0L:2L) {
    Mv <- (gmat == v & called) * 1.0
    n_match <- n_match + tcrossprod(Mv)
  }
  
  conc <- n_match / pmax(n_joint, 1)
  diag(conc) <- NA_real_
  dimnames(conc) <- list(rownames(gmat), rownames(gmat))
  conc
}

# Helper: safely reshape SeqArray FORMAT fields into [sample x variant] matrix.
as_sample_variant_matrix <- function(x, n_sample, n_variant, field_name) {
  if (is.matrix(x)) {
    if (identical(dim(x), c(n_sample, n_variant))) {
      return(x)
    }
    if (identical(dim(x), c(n_variant, n_sample))) {
      return(t(x))
    }
    stop(sprintf(
      "%s matrix has unexpected dimensions: got %s, expected %d x %d (sample x variant).",
      field_name, paste(dim(x), collapse = " x "), n_sample, n_variant
    ))
  }
  
  if (is.list(x) && all(c("length", "data") %in% names(x))) {
    ulen <- unique(x$length)
    if (length(ulen) != 1L || is.na(ulen) || ulen != 1L) {
      stop(sprintf(
        "%s is stored as variable-length data and could not be simplified safely.",
        field_name
      ))
    }
    expected_len <- as.numeric(n_sample) * as.numeric(n_variant)
    if (length(x$data) != expected_len) {
      stop(sprintf(
        "%s data length mismatch: got %d values, expected %.0f.",
        field_name, length(x$data), expected_len
      ))
    }
    return(matrix(x$data, nrow = n_sample, ncol = n_variant))
  }
  
  stop(sprintf("Unsupported structure returned for %s.", field_name))
}

# -----------------------------------------------------------------------------
# 1. Load pedigree / metadata
# -----------------------------------------------------------------------------
cat("Loading pedigree...\n")

ped_full <- read_csv(PED_FILE, col_types = cols(.default = col_character()))
ped_n_rows <- nrow(ped_full)
cat(sprintf("  pedigree_full rows: %d\n", ped_n_rows))
rm(ped_full)
gc()

ped_geno <- read_csv(PED_GENO_FILE, col_types = cols(.default = col_character())) %>%
  mutate(across(everything(), str_trim))

all_sites <- read_csv(ALL_SITES_FILE, col_types = cols(.default = col_character())) %>%
  mutate(across(everything(), str_trim))

# Remove control samples from metadata so they do not trigger later mismatch warnings
ped_geno <- ped_geno %>%
  filter(!genotype_id %in% CONTROL_SAMPLE_IDS)

all_sites <- all_sites %>%
  filter(!genotype.id %in% CONTROL_SAMPLE_IDS)

if (anyDuplicated(ped_geno$genotype_id[!is.na(ped_geno$genotype_id) & ped_geno$genotype_id != ""])) {
  dup_ids <- ped_geno %>%
    filter(!is.na(genotype_id), genotype_id != "") %>%
    count(genotype_id, sort = TRUE) %>%
    filter(n > 1L)
  print(dup_ids)
  stop("Duplicate genotype_id values found in pedigree_genotyped.csv")
}

if (anyDuplicated(all_sites$genotype.id[!is.na(all_sites$genotype.id) & all_sites$genotype.id != ""])) {
  dup_ids_all_sites <- all_sites %>%
    filter(!is.na(genotype.id), genotype.id != "") %>%
    count(genotype.id, sort = TRUE) %>%
    filter(n > 1L)
  warning("Duplicate genotype.id values found in All.sites.csv; using unique() for membership checks.", call. = FALSE)
  print(dup_ids_all_sites)
}

g1_ids <- ped_geno %>% filter(gen == "G1") %>% pull(genotype_id)
g2_ids <- ped_geno %>% filter(gen == "G2") %>% pull(genotype_id)
cat(sprintf("  G1: %d  |  G2: %d\n", length(g1_ids), length(g2_ids)))

# =============================================================================
# 2. Read GDS and restrict to canonical biallelic SNPs
# =============================================================================
cat("\nReading GDS...\n")
stopifnot(file.exists(GDS_FILE))

gds <- seqOpen(GDS_FILE)

gds_open <- TRUE
# on.exit({
#   if (gds_open) try(seqClose(gds), silent = TRUE)
# }, add = TRUE)

seqResetFilter(gds, verbose = FALSE)

all_sample_ids_raw <- seqGetData(gds, "sample.id")
all_variant_ids <- seqGetData(gds, "variant.id")

control_ids_found <- intersect(all_sample_ids_raw, CONTROL_SAMPLE_IDS)
all_sample_ids <- setdiff(all_sample_ids_raw, CONTROL_SAMPLE_IDS)

cat(sprintf("  Raw variants : %d\n", length(all_variant_ids)))
cat(sprintf("  Raw samples  : %d\n", length(all_sample_ids_raw)))
cat(sprintf("  Control samples excluded before QC : %d\n", length(control_ids_found)))
if (length(control_ids_found) > 0L) {
  cat("   ", paste(control_ids_found, collapse = ", "), "\n")
}
cat(sprintf("  Samples entering QC : %d\n", length(all_sample_ids)))

# -- Identify biallelic SNPs --
ref_raw <- seqGetData(gds, "$ref")
alt_raw <- seqGetData(gds, "$alt")

is_biallelic <- !grepl(",", alt_raw, fixed = TRUE)
is_snp <- grepl("^[ACGT]$", ref_raw) & grepl("^[ACGT]$", alt_raw)
keep_variant <- is_biallelic & is_snp

cat(sprintf("  Removed multiallelic / non-SNP variants : %d\n", sum(!keep_variant)))
cat(sprintf("  Biallelic SNPs retained before GT extraction : %d\n", sum(keep_variant)))

keep_variant_ids <- all_variant_ids[keep_variant]
seqSetFilter(
  gds,
  sample.id  = all_sample_ids,
  variant.id = keep_variant_ids,
  verbose    = FALSE
)

rm(ref_raw, alt_raw, keep_variant)
gc()

# -- Build SNP map from filtered variants --
vcf_id <- seqGetData(gds, "annotation/id")
chrom  <- as.character(seqGetData(gds, "chromosome"))
pos    <- as.integer(seqGetData(gds, "position"))
ref    <- as.character(seqGetData(gds, "$ref"))
alt    <- as.character(seqGetData(gds, "$alt"))
qual   <- suppressWarnings(as.numeric(seqGetData(gds, "annotation/qual")))

variant_id_gds <- seqGetData(gds, "variant.id")

snp_id_vec <- ifelse(
  is.na(vcf_id) | vcf_id == "." | vcf_id == "",
  paste0(chrom, ":", pos),
  vcf_id
)

map_raw <- data.frame(
  snp_id         = snp_id_vec,
  variant_id_gds = variant_id_gds,
  chromosome     = chrom,
  position       = pos,
  ref            = ref,
  alt            = alt,
  qual           = qual,
  stringsAsFactors = FALSE
)

n_dup_ids <- sum(duplicated(map_raw$snp_id))
if (n_dup_ids > 0L) {
  warning(sprintf(
    "%d duplicate SNP IDs detected; suffixing with _dup1, _dup2, ...",
    n_dup_ids
  ), call. = FALSE)
  map_raw$snp_id <- make.unique(map_raw$snp_id, sep = "_dup")
}

rm(vcf_id, variant_id_gds, chrom, pos, ref, alt, qual, snp_id_vec)
gc()

# -----------------------------------------------------------------------------
# 2b. Chromosome / contig sanity check
# -----------------------------------------------------------------------------
cat("\n-- Chromosome / contig summary --\n")

contig_snp_counts <- table(map_raw$chromosome)
n_contigs         <- length(contig_snp_counts)
contig_counts_sorted <- sort(contig_snp_counts, decreasing = TRUE)

cat(sprintf("  Total unique contigs : %d\n", n_contigs))
cat(sprintf("  Total SNPs mapped    : %d\n", sum(contig_snp_counts)))

MAJOR_CONTIG_CUMFRAC <- 0.90
cum_frac <- cumsum(as.numeric(contig_counts_sorted)) / sum(contig_snp_counts)
n_major  <- which(cum_frac >= MAJOR_CONTIG_CUMFRAC)[1L]
if (is.na(n_major)) n_major <- n_contigs
major_contigs <- names(contig_counts_sorted)[seq_len(n_major)]
minor_contigs <- setdiff(names(contig_snp_counts), major_contigs)

cat(sprintf("  Major contigs (>= %.0f%% cumulative SNPs) : %d\n",
            MAJOR_CONTIG_CUMFRAC * 100, length(major_contigs)))
cat(sprintf("  Minor scaffolds                          : %d\n", length(minor_contigs)))
cat(sprintf("  SNPs on major contigs : %d (%.1f%%)\n",
            sum(contig_snp_counts[major_contigs]),
            100 * sum(contig_snp_counts[major_contigs]) / sum(contig_snp_counts)))
cat(sprintf("  SNPs on minor scaffolds : %d (%.1f%%)\n",
            sum(contig_snp_counts[minor_contigs]),
            100 * sum(contig_snp_counts[minor_contigs]) / sum(contig_snp_counts)))

cat("\n  Top 20 contigs by SNP count:\n")
top_n <- min(20L, length(contig_counts_sorted))
for (i in seq_len(top_n)) {
  cat(sprintf("    %-30s %7d SNPs  (%.1f%%)\n",
              names(contig_counts_sorted)[i],
              contig_counts_sorted[i],
              100 * contig_counts_sorted[i] / sum(contig_snp_counts)))
}
if (n_contigs > 20L) {
  cat(sprintf("    ... plus %d more contigs\n", n_contigs - 20L))
}

if (length(major_contigs) > 50L) {
  warning(sprintf(
    paste0("50+ contigs needed to reach %.0f%% of SNPs; assembly may be highly ",
           "fragmented. Consider whether contig-level missingness is meaningful."),
    MAJOR_CONTIG_CUMFRAC * 100
  ), call. = FALSE)
}

chr_display <- setNames(
  ifelse(names(contig_snp_counts) %in% major_contigs,
         names(contig_snp_counts),
         "other_scaffolds"),
  names(contig_snp_counts)
)

rm(contig_counts_sorted, cum_frac, n_major, top_n)
gc()

# =============================================================================
# 3. Extract genotypes (and optional dosage) + per-genotype masks  [STAGE A]
# =============================================================================
cat("\n======================================================================\n")
cat("  STAGE A: Genotype extraction + per-genotype DP/GQ masking\n")
cat("======================================================================\n\n")


n_samples_gds  <- length(all_sample_ids)
n_variants_filt <- nrow(map_raw)

cat("  Extracting dosage matrix ($dosage_alt)...\n")
gt_012 <- seqGetData(gds, "$dosage_alt")
gt_012 <- as_sample_variant_matrix(gt_012, n_samples_gds, n_variants_filt, "$dosage_alt")
storage.mode(gt_012) <- "integer"

rownames(gt_012) <- all_sample_ids
colnames(gt_012) <- map_raw$snp_id

cat(sprintf("  Before masking - missing: %.2f%%\n", 100 * mean(is.na(gt_012))))

# -- Dosage extraction from RO / AO (freeBayes FORMAT fields) --
dosage_mat <- NULL
dosage_available <- FALSE

if (USE_DOSAGE) {
  cat("\n  [Dosage mode ON] Attempting to extract RO/AO from GDS...\n")
  
  ro_ok <- tryCatch({
    ro_raw <- seqGetData(gds, "annotation/format/RO")
    TRUE
  }, error = function(e) FALSE)
  
  ao_ok <- tryCatch({
    ao_raw <- seqGetData(gds, "annotation/format/AO")
    TRUE
  }, error = function(e) FALSE)
  
  if (ro_ok && ao_ok) {
    cat("  RO and AO fields found -- computing read-count dosage.\n")
    
    ro_mat <- as_sample_variant_matrix(ro_raw, n_samples_gds, n_variants_filt, "RO")
    ao_mat <- as_sample_variant_matrix(ao_raw, n_samples_gds, n_variants_filt, "AO")
    storage.mode(ro_mat) <- "numeric"
    storage.mode(ao_mat) <- "numeric"
    rm(ro_raw, ao_raw)
    
    total_reads <- ro_mat + ao_mat
    
    dosage_mat <- ifelse(
      total_reads > 0 & !is.na(gt_012),
      2.0 * ao_mat / total_reads,
      NA_real_
    )
    
    rownames(dosage_mat) <- all_sample_ids
    colnames(dosage_mat) <- map_raw$snp_id
    
    cat(sprintf("  Dosage matrix: %d samples x %d SNPs\n", nrow(dosage_mat), ncol(dosage_mat)))
    cat(sprintf("  Dosage mean (non-NA): %.4f\n", mean(dosage_mat, na.rm = TRUE)))
    cat(sprintf("  Dosage missing: %.2f%%\n", 100 * mean(is.na(dosage_mat))))
    
    site_mean_dosage <- colMeans(dosage_mat, na.rm = TRUE)
    site_mean_dosage[is.nan(site_mean_dosage)] <- NA_real_
    
    rm(ro_mat, ao_mat, total_reads)
    gc()
    
    dosage_available <- TRUE
  } else {
    if (!ro_ok) cat("  WARNING: annotation/format/RO not found in GDS.\n")
    if (!ao_ok) cat("  WARNING: annotation/format/AO not found in GDS.\n")
    cat("  Dosage mode requested but RO/AO unavailable -- falling back to 0/1/2 only.\n")
    USE_DOSAGE <- FALSE
  }
  
}

# -- DP extraction and masking --
need_dp <- DP_MIN > 0 || DP_MAX_MULT > 0 || SNP_MEAN_DP_MAX > 0 || SAMPLE_DP_MIN > 0
if (need_dp) {
  cat("  Extracting DP matrix (annotation/format/DP)...\n")
  dp_raw <- seqGetData(gds, "annotation/format/DP")
  dp_mat <- as_sample_variant_matrix(dp_raw, n_samples_gds, n_variants_filt, "DP")
  storage.mode(dp_mat) <- "numeric"
  rm(dp_raw)
  
  global_mean_dp <- mean(dp_mat, na.rm = TRUE)
  cat(sprintf("  Global mean genotype DP : %.2f\n", global_mean_dp))
  
  if (DP_MIN > 0) {
    mask_dp_low <- !is.na(dp_mat) & dp_mat < DP_MIN
    n_masked_dp_low <- sum(mask_dp_low & !is.na(gt_012))
    cat(sprintf("  Applying DP < %d mask... (%d genotypes masked)\n", DP_MIN, n_masked_dp_low))
    gt_012[mask_dp_low] <- NA_integer_
    if (dosage_available) dosage_mat[mask_dp_low] <- NA_real_
    rm(mask_dp_low)
  }
  
  effective_dp_max <- 0
  if (DP_MAX_MULT > 0) {
    effective_dp_max <- round(DP_MAX_MULT * global_mean_dp)
    mask_dp_high <- !is.na(dp_mat) & dp_mat > effective_dp_max
    n_masked_dp_high <- sum(mask_dp_high & !is.na(gt_012))
    cat(sprintf("  Applying DP > %d mask (%.0fx mean)... (%d genotypes masked)\n",
                effective_dp_max, DP_MAX_MULT, n_masked_dp_high))
    gt_012[mask_dp_high] <- NA_integer_
    if (dosage_available) dosage_mat[mask_dp_high] <- NA_real_
    rm(mask_dp_high)
  }
  
  snp_mean_dp <- colMeans(dp_mat, na.rm = TRUE)
  snp_mean_dp[is.nan(snp_mean_dp)] <- NA_real_
  
  sample_mean_dp_raw <- rowMeans(dp_mat, na.rm = TRUE)
  sample_mean_dp_raw[is.nan(sample_mean_dp_raw)] <- NA_real_
  names(sample_mean_dp_raw) <- all_sample_ids
  
  rm(dp_mat)
  gc()
} else {
  snp_mean_dp <- rep(NA_real_, n_variants_filt)
  sample_mean_dp_raw <- rep(NA_real_, n_samples_gds)
  names(sample_mean_dp_raw) <- all_sample_ids
  effective_dp_max <- 0
  global_mean_dp <- NA_real_
}

# -- GQ extraction and masking (disabled by default) --
if (GQ_MIN > 0) {
  cat(sprintf("  Extracting GQ matrix (annotation/format/GQ)...\n"))
  gq_raw <- seqGetData(gds, "annotation/format/GQ")
  gq_mat <- as_sample_variant_matrix(gq_raw, n_samples_gds, n_variants_filt, "GQ")
  storage.mode(gq_mat) <- "numeric"
  rm(gq_raw)
  
  mask_gq <- !is.na(gq_mat) & gq_mat < GQ_MIN
  n_masked_gq <- sum(mask_gq & !is.na(gt_012))
  cat(sprintf("  Applying GQ < %d mask... (%d genotypes masked)\n", GQ_MIN, n_masked_gq))
  gt_012[mask_gq] <- NA_integer_
  if (dosage_available) dosage_mat[mask_gq] <- NA_real_
  rm(gq_mat, mask_gq)
  gc()
}

cat(sprintf("  After masking  - missing: %.2f%%\n", 100 * mean(is.na(gt_012))))
if (dosage_available) {
  cat(sprintf("  Dosage after masking - missing: %.2f%%\n", 100 * mean(is.na(dosage_mat))))
}

# -- Close the GDS handle --
seqClose(gds)
gds_open <- FALSE
cat("  GDS handle closed (all data in memory).\n")

cat(sprintf("  Genotype matrix: %d samples x %d SNPs\n", nrow(gt_012), ncol(gt_012)))


raw_n_samples <- nrow(gt_012)
raw_n_snps    <- ncol(gt_012)

if (anyDuplicated(rownames(gt_012))) {
  dup_vcf_ids <- rownames(gt_012)[duplicated(rownames(gt_012))]
  print(sort(table(dup_vcf_ids), decreasing = TRUE))
  stop("Duplicate sample IDs found in the genotype matrix")
}

rm(all_sample_ids_raw, all_sample_ids, all_variant_ids, keep_variant_ids,
   control_ids_found, n_samples_gds, n_variants_filt)
gc()

# -----------------------------------------------------------------------------
# 4. Sample ID verification
# -----------------------------------------------------------------------------
cat("\n-- Sample ID verification --\n")

vcf_sample_ids <- rownames(gt_012)
known_gids <- all_sites %>%
  filter(!is.na(genotype.id), genotype.id != "", genotype.id != "NA") %>%
  pull(genotype.id) %>%
  unique()

in_vcf_not_ped <- setdiff(vcf_sample_ids, known_gids)
in_ped_not_vcf <- setdiff(known_gids, vcf_sample_ids)

cat(sprintf("  VCF samples in pedigree : %d / %d\n",
            length(intersect(vcf_sample_ids, known_gids)),
            length(vcf_sample_ids)))

if (length(in_vcf_not_ped) > 0L) {
  cat("  WARNING - in VCF but not pedigree:\n")
  cat("   ", paste(in_vcf_not_ped, collapse = ", "), "\n")
}
if (length(in_ped_not_vcf) > 0L) {
  cat("  WARNING - in pedigree but not VCF:\n")
  cat("   ", paste(head(in_ped_not_vcf, 20), collapse = ", "), "\n")
  if (length(in_ped_not_vcf) > 20L) {
    cat("   ... and", length(in_ped_not_vcf) - 20L, "more\n")
  }
}

rm(known_gids, in_vcf_not_ped, in_ped_not_vcf)
gc()

# =============================================================================
# 5a. Coarse sample pre-filter: dead + low-depth samples  [STAGE B]
# =============================================================================
cat("\n======================================================================\n")
cat("  STAGE B: Coarse sample pre-filter (dead + low-depth samples)\n")
cat("======================================================================\n\n")


raw_sample_cr <- sample_call_rate(gt_012)
dead_mask_cr <- raw_sample_cr < SAMPLE_CR_DEAD
n_dead_cr <- sum(dead_mask_cr, na.rm = TRUE)

cat(sprintf("  Call rate threshold for dead samples : %.0f%%\n", SAMPLE_CR_DEAD * 100))
cat(sprintf("  Samples with call rate < %.0f%%       : %d\n", SAMPLE_CR_DEAD * 100, n_dead_cr))

sample_dp_for_filter <- sample_mean_dp_raw[rownames(gt_012)]
dead_mask_dp <- if (SAMPLE_DP_MIN > 0 && !all(is.na(sample_dp_for_filter))) {
  !is.na(sample_dp_for_filter) & sample_dp_for_filter < SAMPLE_DP_MIN
} else {
  rep(FALSE, nrow(gt_012))
}
n_dead_dp <- sum(dead_mask_dp & !dead_mask_cr, na.rm = TRUE)

cat(sprintf("  Mean DP threshold                   : >= %d\n", SAMPLE_DP_MIN))
cat(sprintf("  Samples with mean DP < %d (new)     : %d\n", SAMPLE_DP_MIN, n_dead_dp))

dead_mask <- dead_mask_cr | dead_mask_dp
n_dead <- sum(dead_mask, na.rm = TRUE)

if (n_dead > 0L) {
  dead_ids <- rownames(gt_012)[dead_mask]
  cat(sprintf("  Removing %d samples total (CR + DP)...\n", n_dead))
  gt_012 <- gt_012[!dead_mask, , drop = FALSE]
  if (dosage_available) dosage_mat <- dosage_mat[!dead_mask, , drop = FALSE]
  cat(sprintf("  Samples remaining after coarse filter: %d\n", nrow(gt_012)))
} else {
  dead_ids <- character(0)
  cat("  No dead/low-depth samples to remove.\n")
}

rm(raw_sample_cr, dead_mask_cr, dead_mask_dp, dead_mask, sample_dp_for_filter)
gc()


# =============================================================================
# 5b. Site filtering pass 1  [STAGE C]
# =============================================================================
cat("\n======================================================================\n")
cat("  STAGE C: Site filtering pass 1 (call rate, QUAL, clusters, monomorphic)\n")
cat("======================================================================\n\n")


n_snps_pre_c <- ncol(gt_012)

snp_cr_vec_c  <- snp_call_rate(gt_012)
snp_maf_vec_c <- snp_maf(gt_012)

snp_fail_cr_c <- snp_cr_vec_c < SNP_CR
snp_monomorphic_c <- !is.na(snp_maf_vec_c) & snp_maf_vec_c == 0

snp_fail_dp_mean_c <- if (SNP_MEAN_DP_MAX > 0) {
  !is.na(snp_mean_dp) & snp_mean_dp > SNP_MEAN_DP_MAX
} else {
  rep(FALSE, n_snps_pre_c)
}

snp_qual_c <- map_raw$qual[match(colnames(gt_012), map_raw$snp_id)]
snp_fail_qual_c <- if (QUAL_MIN > 0) {
  !is.na(snp_qual_c) & snp_qual_c < QUAL_MIN
} else {
  rep(FALSE, n_snps_pre_c)
}

snp_chrom_c <- map_raw$chromosome[match(colnames(gt_012), map_raw$snp_id)]
snp_pos_c   <- map_raw$position[match(colnames(gt_012), map_raw$snp_id)]
snp_fail_cluster_c <- if (SNP_CLUSTER_BP > 0 && SNP_CLUSTER_MAX > 0) {
  flag_snp_clusters(snp_chrom_c, snp_pos_c, SNP_CLUSTER_BP, SNP_CLUSTER_MAX)
} else {
  rep(FALSE, n_snps_pre_c)
}

snp_remove_c <- snp_fail_cr_c | snp_monomorphic_c | snp_fail_dp_mean_c |
  snp_fail_qual_c | snp_fail_cluster_c

cat(sprintf("  SNPs entering Stage C               : %d\n", n_snps_pre_c))
cat(sprintf("  Fail call rate < %.0f%%               : %d\n", SNP_CR * 100, sum(snp_fail_cr_c, na.rm = TRUE)))
if (QUAL_MIN > 0) {
  cat(sprintf("  Fail QUAL < %d                      : %d\n", QUAL_MIN, sum(snp_fail_qual_c, na.rm = TRUE)))
}
cat(sprintf("  Monomorphic (MAF = 0)               : %d\n", sum(snp_monomorphic_c, na.rm = TRUE)))
if (SNP_CLUSTER_BP > 0) {
  cat(sprintf("  Fail SNP cluster (> %d in %dbp)      : %d\n",
              SNP_CLUSTER_MAX, SNP_CLUSTER_BP, sum(snp_fail_cluster_c, na.rm = TRUE)))
}
if (SNP_MEAN_DP_MAX > 0) {
  cat(sprintf("  Fail raw mean DP > %.0f              : %d\n", SNP_MEAN_DP_MAX, sum(snp_fail_dp_mean_c, na.rm = TRUE)))
}
cat(sprintf("  Total removed in Stage C            : %d\n", sum(snp_remove_c, na.rm = TRUE)))

stage_c_audit <- data.frame(
  snp_id            = colnames(gt_012),
  fail_stage_c      = snp_remove_c,
  fail_call_rate    = snp_fail_cr_c,
  fail_monomorphic  = snp_monomorphic_c,
  fail_mean_dp      = snp_fail_dp_mean_c,
  fail_qual         = snp_fail_qual_c,
  fail_cluster      = snp_fail_cluster_c,
  stringsAsFactors = FALSE
)

write_csv(
  stage_c_audit %>% filter(fail_stage_c),
  file.path(OUT_DIR, "stage_c_snp_failures.csv")
)

snps_removed_stage_c <- colnames(gt_012)[snp_remove_c]
gt_012 <- gt_012[, !snp_remove_c, drop = FALSE]
if (dosage_available) dosage_mat <- dosage_mat[, !snp_remove_c, drop = FALSE]
snp_mean_dp <- snp_mean_dp[!snp_remove_c]

cat(sprintf("  SNPs remaining after Stage C        : %d\n", ncol(gt_012)))

rm(snp_cr_vec_c, snp_maf_vec_c, snp_fail_cr_c, snp_monomorphic_c, snp_fail_dp_mean_c,
   snp_fail_qual_c, snp_fail_cluster_c, snp_qual_c, snp_chrom_c, snp_pos_c, snp_remove_c)
gc()


# =============================================================================
# 5c. Excess het + HWE excess-het direction filter on G1  [STAGE D]
# =============================================================================
cat("\n======================================================================\n")
cat("  STAGE D: Excess het + HWE excess-het filter (G1 only)\n")
cat("======================================================================\n\n")


g1_keep_d <- intersect(rownames(gt_012), g1_ids)
cat(sprintf("  G1 individuals available : %d\n", length(g1_keep_d)))

if (length(g1_keep_d) < 10L) {
  warning("Very few G1 individuals remain; SNP excess-het filter is unstable", call. = FALSE)
}

gt_g1_d <- gt_012[g1_keep_d, , drop = FALSE]

maf_g1_d <- snp_maf(gt_g1_d)
het_obs_g1_d <- snp_obs_het(gt_g1_d)
het_exp_g1_d <- 2 * maf_g1_d * (1 - maf_g1_d)
het_excess_ratio_d <- het_obs_g1_d / het_exp_g1_d
het_excess_ratio_d[is.nan(het_excess_ratio_d) | is.infinite(het_excess_ratio_d)] <- NA_real_

snp_fail_het_ratio <- !is.na(het_excess_ratio_d) &
  !is.na(het_exp_g1_d) &
  het_exp_g1_d >= HEXP_MIN &
  het_excess_ratio_d > HET_EXCESS_RATIO_MAX

n_het_excess_skip_d <- sum(is.na(het_excess_ratio_d) | is.na(het_exp_g1_d) | het_exp_g1_d < HEXP_MIN)

cat(sprintf("  Fail excess het ratio (Hobs/Hexp > %.2f) : %d\n",
            HET_EXCESS_RATIO_MAX, sum(snp_fail_het_ratio, na.rm = TRUE)))
cat(sprintf("  Excess-het skipped / NA (Hexp < %.2f)    : %d\n",
            HEXP_MIN, n_het_excess_skip_d))

n_00_d <- colSums(gt_g1_d == 0L, na.rm = TRUE)
n_01_d <- colSums(gt_g1_d == 1L, na.rm = TRUE)
n_11_d <- colSums(gt_g1_d == 2L, na.rm = TRUE)
n_total_d <- n_00_d + n_01_d + n_11_d

p_ref_d <- ifelse(n_total_d >= 10, (2 * n_00_d + n_01_d) / (2 * n_total_d), NA_real_)
q_alt_d <- ifelse(!is.na(p_ref_d), 1 - p_ref_d, NA_real_)

exp_00_d <- n_total_d * p_ref_d^2
exp_01_d <- n_total_d * 2 * p_ref_d * q_alt_d
exp_11_d <- n_total_d * q_alt_d^2

term_00_d <- ifelse(!is.na(exp_00_d) & exp_00_d > 0, (n_00_d - exp_00_d)^2 / exp_00_d, 0)
term_01_d <- ifelse(!is.na(exp_01_d) & exp_01_d > 0, (n_01_d - exp_01_d)^2 / exp_01_d, 0)
term_11_d <- ifelse(!is.na(exp_11_d) & exp_11_d > 0, (n_11_d - exp_11_d)^2 / exp_11_d, 0)

chi2_d <- ifelse(n_total_d >= 10, term_00_d + term_01_d + term_11_d, NA_real_)
hwe_pval_d <- ifelse(!is.na(chi2_d), pchisq(chi2_d, df = 1, lower.tail = FALSE), NA_real_)

is_excess_het_direction <- !is.na(het_obs_g1_d) & !is.na(het_exp_g1_d) & het_obs_g1_d > het_exp_g1_d

snp_fail_hwe_excess <- if (HWE_EXCESS_PVAL > 0) {
  !is.na(hwe_pval_d) & hwe_pval_d < HWE_EXCESS_PVAL & is_excess_het_direction
} else {
  rep(FALSE, ncol(gt_012))
}

n_hwe_only <- sum(snp_fail_hwe_excess & !snp_fail_het_ratio, na.rm = TRUE)
cat(sprintf("  Fail HWE excess-het (p < %s, G1)     : %d\n",
            formatC(HWE_EXCESS_PVAL, format = "e", digits = 0),
            sum(snp_fail_hwe_excess, na.rm = TRUE)))
cat(sprintf("  HWE excess-het only (not caught by ratio): %d\n", n_hwe_only))

snp_fail_het_d <- snp_fail_het_ratio | snp_fail_hwe_excess

cat(sprintf("  Total fail Stage D (ratio OR HWE)    : %d\n", sum(snp_fail_het_d, na.rm = TRUE)))

# ── PATCH 2: Het-ratio sensitivity analysis (v4.6.1) ────────────────────────
# Now shows ratio-only, HWE-only, and combined fail counts separately at each
# threshold so the user can see which filter is driving removals.
cat("\n  Het-ratio sensitivity (alternative thresholds):\n")
het_sensitivity_thresholds <- c(1.20, 1.25, 1.30, 1.35, 1.40, 1.50)
evaluable_mask <- !is.na(het_excess_ratio_d) & !is.na(het_exp_g1_d) & het_exp_g1_d >= HEXP_MIN
n_evaluable <- sum(evaluable_mask)
cat(sprintf("    Evaluable SNPs (Hexp >= %.2f) : %d\n", HEXP_MIN, n_evaluable))

n_hwe_fail_total <- sum(snp_fail_hwe_excess, na.rm = TRUE)
cat(sprintf("    HWE excess-het fails (all thresholds): %d\n", n_hwe_fail_total))

cat(sprintf("\n    %-10s  %10s  %10s  %10s  %10s  %s\n",
            "Threshold", "Ratio-only", "HWE-only", "Combined", "Survive", ""))
for (thr_i in het_sensitivity_thresholds) {
  fail_ratio_i <- evaluable_mask & het_excess_ratio_d > thr_i
  n_fail_ratio_i <- sum(fail_ratio_i, na.rm = TRUE)
  # HWE-only = fail HWE but would NOT fail at this ratio threshold
  n_hwe_only_i <- sum(snp_fail_hwe_excess & !fail_ratio_i, na.rm = TRUE)
  n_fail_combined_i <- sum(fail_ratio_i | snp_fail_hwe_excess, na.rm = TRUE)
  n_survive_i <- ncol(gt_012) - n_fail_combined_i
  marker <- if (thr_i == HET_EXCESS_RATIO_MAX) " <-- CURRENT" else ""
  cat(sprintf("    > %-7.2f  %10d  %10d  %10d  %10d%s\n",
              thr_i, n_fail_ratio_i, n_hwe_only_i, n_fail_combined_i, n_survive_i, marker))
}

# ── PATCH 1: Marginal SNPs for Mendelian-error stratification (v4.6.1) ──────
# v4.6.1 fix: removed the !snp_fail_hwe_excess condition. In datasets where
# the HWE chi-squared test catches nearly every ratio-fail SNP, requiring
# HWE-pass yielded zero marginal SNPs, defeating the comparison. Now the
# marginal set includes all SNPs with ratio in (threshold, 1.40] regardless
# of HWE status, since the purpose is to compare Mendelian error rates across
# the het-ratio gradient.
marginal_het_mask <- evaluable_mask &
  het_excess_ratio_d > HET_EXCESS_RATIO_MAX &
  het_excess_ratio_d <= 1.40

n_marginal <- sum(marginal_het_mask)
cat(sprintf("\n  Marginal SNPs saved for Mendelian stratification (ratio %.2f-1.40): %d\n",
            HET_EXCESS_RATIO_MAX, n_marginal))

if (n_marginal > 0L) {
  gt_marginal_het <- gt_012[, marginal_het_mask, drop = FALSE]
  marginal_het_info <- data.frame(
    snp_id           = colnames(gt_012)[marginal_het_mask],
    het_excess_ratio = het_excess_ratio_d[marginal_het_mask],
    stringsAsFactors = FALSE
  )
} else {
  gt_marginal_het <- NULL
  marginal_het_info <- data.frame(snp_id = character(), het_excess_ratio = numeric())
}

# Also save het-ratio for all SNPs that PASS Stage D, for binned comparison
het_ratio_all_pre_d <- data.frame(
  snp_id           = colnames(gt_012),
  het_excess_ratio = het_excess_ratio_d,
  fail_stage_d     = snp_fail_het_ratio | snp_fail_hwe_excess,
  stringsAsFactors = FALSE
)

het_stats_stage_d <- data.frame(
  snp_id           = colnames(gt_012),
  maf_g1           = maf_g1_d,
  het_obs_g1       = het_obs_g1_d,
  het_exp_g1       = het_exp_g1_d,
  het_excess_ratio = het_excess_ratio_d,
  hwe_pval_g1      = hwe_pval_d,
  fail_het_excess  = snp_fail_het_d,
  stringsAsFactors = FALSE
)

stage_d_audit <- data.frame(
  snp_id            = colnames(gt_012),
  fail_stage_d      = snp_fail_het_d,
  fail_het_ratio    = snp_fail_het_ratio,
  fail_hwe_excess   = snp_fail_hwe_excess,
  maf_g1            = maf_g1_d,
  het_obs_g1        = het_obs_g1_d,
  het_exp_g1        = het_exp_g1_d,
  het_excess_ratio  = het_excess_ratio_d,
  hwe_pval_g1       = hwe_pval_d,
  stringsAsFactors = FALSE
)

write_csv(
  stage_d_audit %>% filter(fail_stage_d),
  file.path(OUT_DIR, "stage_d_snp_failures.csv")
)

snps_removed_stage_d <- colnames(gt_012)[snp_fail_het_d]
gt_012 <- gt_012[, !snp_fail_het_d, drop = FALSE]
if (dosage_available) dosage_mat <- dosage_mat[, !snp_fail_het_d, drop = FALSE]
snp_mean_dp <- snp_mean_dp[!snp_fail_het_d]
het_stats_stage_d <- het_stats_stage_d[!snp_fail_het_d, ]

cat(sprintf("  SNPs remaining after Stage D         : %d\n", ncol(gt_012)))

rm(gt_g1_d, maf_g1_d, het_obs_g1_d, het_exp_g1_d, het_excess_ratio_d,
   snp_fail_het_ratio, snp_fail_hwe_excess, is_excess_het_direction,
   n_00_d, n_01_d, n_11_d, n_total_d, p_ref_d, q_alt_d,
   exp_00_d, exp_01_d, exp_11_d, term_00_d, term_01_d, term_11_d,
   chi2_d, hwe_pval_d, snp_fail_het_d, evaluable_mask)
gc()


# =============================================================================
# 5d. Refined sample QC  [STAGE E]
# =============================================================================
cat("\n======================================================================\n")
cat("  STAGE E: Refined sample QC (recomputed on clean sites)\n")
cat("======================================================================\n\n")


cat(sprintf("  Computing sample stats on %d clean SNPs...\n", ncol(gt_012)))

sample_qc <- data.frame(
  sample_id = rownames(gt_012),
  call_rate = sample_call_rate(gt_012),
  obs_het   = sample_het(gt_012)
) %>%
  left_join(
    ped_geno %>%
      select(genotype_id, gen, site) %>%
      distinct(genotype_id, .keep_all = TRUE) %>%
      rename(sample_id = genotype_id),
    by = "sample_id"
  ) %>%
  mutate(
    gen = if_else(is.na(gen) | gen == "", "UNKNOWN", gen),
    site = if_else(is.na(site) | site == "", "UNKNOWN", site),
    cr_flag = case_when(
      call_rate < SAMPLE_CR_FAIL ~ "FAIL",
      call_rate < SAMPLE_CR_WARN ~ "WARN",
      TRUE                       ~ "PASS"
    )
  )

global_het_mean <- safe_mean(sample_qc$obs_het)
global_het_sd   <- safe_sd(sample_qc$obs_het)
if (is.na(global_het_sd) || global_het_sd == 0) {
  global_het_sd <- 1
}

sample_qc <- sample_qc %>%
  group_by(gen) %>%
  mutate(
    group_n_nonmissing = sum(!is.na(obs_het)),
    het_mean_group = if_else(group_n_nonmissing >= 20L, safe_mean(obs_het), global_het_mean),
    het_sd_group   = if_else(group_n_nonmissing >= 20L, safe_sd(obs_het), global_het_sd),
    het_mean_group = if_else(is.na(het_mean_group), global_het_mean, het_mean_group),
    het_sd_group   = if_else(is.na(het_sd_group) | het_sd_group == 0, global_het_sd, het_sd_group),
    het_lower      = het_mean_group - HET_SD * het_sd_group,
    het_upper      = het_mean_group + HET_SD * het_sd_group,
    het_z          = (obs_het - het_mean_group) / het_sd_group
  ) %>%
  ungroup() %>%
  mutate(
    het_flag = case_when(
      is.na(obs_het)      ~ "NA",
      is.na(het_z)        ~ "NA",
      abs(het_z) > HET_SD ~ "OUTLIER",
      TRUE                ~ "PASS"
    )
  )

# -- Per-sample inbreeding coefficient (F_hat) --
cat("  Computing per-sample F_hat (inbreeding coefficient)...\n")

snp_p_e <- colMeans(gt_012, na.rm = TRUE) / 2
snp_hexp_e <- 2 * snp_p_e * (1 - snp_p_e)
called_mat_e <- !is.na(gt_012)
n_called_per_sample <- rowSums(called_mat_e)

sample_exp_het <- as.vector(called_mat_e %*% snp_hexp_e) /
  pmax(n_called_per_sample, 1L)
sample_exp_het[n_called_per_sample == 0L] <- NA_real_

f_hat <- 1 - (sample_qc$obs_het / sample_exp_het)
f_hat[is.nan(f_hat) | is.infinite(f_hat)] <- NA_real_

sample_qc$exp_het <- sample_exp_het
sample_qc$f_hat   <- f_hat

cat(sprintf("  F_hat summary (all samples):\n"))
cat(sprintf("    Mean: %.4f  Median: %.4f  SD: %.4f\n",
            mean(f_hat, na.rm = TRUE), median(f_hat, na.rm = TRUE),
            sd(f_hat, na.rm = TRUE)))

for (gen_i in sort(unique(sample_qc$gen))) {
  f_gen <- sample_qc$f_hat[sample_qc$gen == gen_i]
  cat(sprintf("    %s (n=%d): mean=%.4f  range=[%.4f, %.4f]\n",
              gen_i, sum(!is.na(f_gen)),
              mean(f_gen, na.rm = TRUE),
              min(f_gen, na.rm = TRUE),
              max(f_gen, na.rm = TRUE)))
}

n_high_f <- sum(!is.na(f_hat) & f_hat > 0.10)
if (n_high_f > 0L) {
  cat(sprintf("  WARNING: %d samples with F_hat > 0.10\n", n_high_f))
}

rm(snp_p_e, snp_hexp_e, called_mat_e, n_called_per_sample, sample_exp_het, f_hat)

cat(sprintf("  Call rate < %.0f%% (WARN) : %d samples\n",
            SAMPLE_CR_WARN * 100, sum(sample_qc$cr_flag == "WARN", na.rm = TRUE)))
cat(sprintf("  Call rate < %.0f%% (FAIL) : %d samples\n",
            SAMPLE_CR_FAIL * 100, sum(sample_qc$cr_flag == "FAIL", na.rm = TRUE)))
cat(sprintf("  Het outliers (>%.0f SD within generation) : %d samples\n",
            HET_SD, sum(sample_qc$het_flag == "OUTLIER", na.rm = TRUE)))

samples_remove_e <- sample_qc %>%
  filter(cr_flag == "FAIL" | het_flag == "OUTLIER") %>%
  pull(sample_id)

n_removed_stage_e <- length(samples_remove_e)

cat(sprintf("  Total flagged for removal (Stage E)  : %d\n", length(samples_remove_e)))

all_samples_removed <- union(dead_ids, samples_remove_e)
cat(sprintf("  Total samples removed (Stage B + E)  : %d\n", length(all_samples_removed)))

if (length(dead_ids) > 0L) {
  dead_sample_rows <- data.frame(
    sample_id = dead_ids,
    call_rate = NA_real_,
    obs_het   = NA_real_,
    cr_flag   = "DEAD",
    het_flag  = "NA",
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      ped_geno %>%
        select(genotype_id, gen, site) %>%
        distinct(genotype_id, .keep_all = TRUE) %>%
        rename(sample_id = genotype_id),
      by = "sample_id"
    ) %>%
    mutate(
      gen = if_else(is.na(gen) | gen == "", "UNKNOWN", gen),
      site = if_else(is.na(site) | site == "", "UNKNOWN", site)
    )
  sample_qc_full <- bind_rows(sample_qc, dead_sample_rows)
} else {
  sample_qc_full <- sample_qc
}
write_csv(sample_qc_full, file.path(OUT_DIR, "sample_qc_summary.csv"))

het_thresholds <- sample_qc %>%
  distinct(gen, het_lower, het_upper)

p_het <- ggplot(sample_qc, aes(x = call_rate, y = obs_het, colour = gen)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(
    data = het_thresholds,
    aes(yintercept = het_lower),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_hline(
    data = het_thresholds,
    aes(yintercept = het_upper),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_vline(xintercept = SAMPLE_CR_WARN, linetype = "dotted", colour = "orange") +
  geom_vline(xintercept = SAMPLE_CR_FAIL, linetype = "dashed", colour = "red") +
  facet_wrap(~ gen, scales = "free_y") +
  scale_colour_manual(values = c(G1 = "#2166ac", G2 = "#d6604d", UNKNOWN = "grey50"), drop = FALSE) +
  labs(
    title = "Sample QC: call rate vs heterozygosity (Stage E -- clean sites)",
    subtitle = sprintf("n=%d (excl. %d dead) | flagged=%d",
                       nrow(sample_qc), length(dead_ids), length(samples_remove_e)),
    x = "Call rate (recomputed on clean sites)",
    y = "Observed heterozygosity",
    colour = "Generation"
  ) +
  theme_bw()

ggsave(file.path(OUT_DIR, "sample_qc_hetplot.png"), p_het, width = 9, height = 5, dpi = 150)
cat("  Plot saved: sample_qc_hetplot.png\n")

samples_keep <- setdiff(rownames(gt_012), samples_remove_e)
gt_filt <- gt_012[samples_keep, , drop = FALSE]
if (dosage_available) dosage_filt <- dosage_mat[samples_keep, , drop = FALSE]
if (!is.null(gt_marginal_het)) {
  gt_marginal_het <- gt_marginal_het[intersect(samples_keep, rownames(gt_marginal_het)), , drop = FALSE]
}
cat(sprintf("\n  Samples retained after Stage E: %d / %d (excl. %d dead from Stage B)\n",
            nrow(gt_filt), raw_n_samples, length(dead_ids)))

rm(gt_012, samples_remove_e, p_het, het_thresholds, dead_ids)
if (dosage_available) rm(dosage_mat)
if (exists("dead_sample_rows")) rm(dead_sample_rows)
gc()

if (nrow(gt_filt) == 0L) stop("No samples remain after sample-level QC")
if (ncol(gt_filt) == 0L) stop("No SNPs remain after site filtering")

sample_keep_meta <- sample_qc %>%
  filter(sample_id %in% rownames(gt_filt)) %>%
  select(sample_id, gen, site)

sample_missing_by_group <- sample_qc %>%
  filter(sample_id %in% rownames(gt_filt)) %>%
  group_by(gen, site) %>%
  summarise(
    n_samples        = n(),
    mean_call_rate   = mean(call_rate, na.rm = TRUE),
    median_call_rate = median(call_rate, na.rm = TRUE),
    min_call_rate    = min(call_rate, na.rm = TRUE),
    n_fail_callrate  = sum(cr_flag == "FAIL", na.rm = TRUE),
    n_warn_callrate  = sum(cr_flag == "WARN", na.rm = TRUE),
    n_het_outlier    = sum(het_flag == "OUTLIER", na.rm = TRUE),
    .groups = "drop"
  )
write_csv(sample_missing_by_group, file.path(OUT_DIR, "sample_missingness_by_group.csv"))


# =============================================================================
# 5e. Site filtering pass 2: MAF on final sample set  [STAGE F]
# =============================================================================
cat("\n======================================================================\n")
cat("  STAGE F: Site filtering pass 2 (MAF on final sample set)\n")
cat("======================================================================\n\n")


n_snps_pre_f <- ncol(gt_filt)
snp_maf_final <- snp_maf(gt_filt)

snp_fail_maf_f <- !is.na(snp_maf_final) & snp_maf_final < MAF_MIN
snp_newly_mono <- !is.na(snp_maf_final) & snp_maf_final == 0

snp_remove_f <- snp_fail_maf_f | snp_newly_mono | is.na(snp_maf_final)

cat(sprintf("  SNPs entering Stage F               : %d\n", n_snps_pre_f))
cat(sprintf("  Fail MAF < %.2f                     : %d\n", MAF_MIN, sum(snp_fail_maf_f, na.rm = TRUE)))
cat(sprintf("  Newly monomorphic after sample QC   : %d\n", sum(snp_newly_mono & !snp_fail_maf_f, na.rm = TRUE)))
cat(sprintf("  All-NA after sample QC              : %d\n", sum(is.na(snp_maf_final) & !snp_fail_maf_f & !snp_newly_mono)))
cat(sprintf("  Total removed in Stage F            : %d\n", sum(snp_remove_f, na.rm = TRUE)))

stage_f_audit <- data.frame(
  snp_id              = colnames(gt_filt),
  fail_stage_f        = snp_remove_f,
  fail_maf            = snp_fail_maf_f,
  newly_monomorphic   = snp_newly_mono,
  maf_final_prefilter = snp_maf_final,
  stringsAsFactors = FALSE
)

write_csv(
  stage_f_audit %>% filter(fail_stage_f),
  file.path(OUT_DIR, "stage_f_snp_failures.csv")
)

snps_removed_stage_f <- colnames(gt_filt)[snp_remove_f]
gt_filt <- gt_filt[, !snp_remove_f, drop = FALSE]
if (dosage_available) dosage_filt <- dosage_filt[, !snp_remove_f, drop = FALSE]
snp_mean_dp <- snp_mean_dp[!snp_remove_f]
het_stats_stage_d <- het_stats_stage_d[!snp_remove_f, ]

snp_maf_final_pass <- snp_maf(gt_filt)

cat(sprintf("  SNPs remaining after Stage F (pass QC): %d\n", ncol(gt_filt)))

# -----------------------------------------------------------------------------
# 5f. All-SNP audit ledger
# -----------------------------------------------------------------------------
cat("\n-- Writing all-SNP audit ledger --\n")

snp_qc_audit_all <- map_raw %>%
  select(snp_id, chromosome, position, ref, alt, qual) %>%
  left_join(stage_c_audit, by = "snp_id") %>%
  left_join(stage_d_audit, by = "snp_id") %>%
  left_join(stage_f_audit, by = "snp_id") %>%
  mutate(
    fail_stage_c = coalesce(fail_stage_c, FALSE),
    fail_stage_d = coalesce(fail_stage_d, FALSE),
    fail_stage_f = coalesce(fail_stage_f, FALSE),
    fail_any     = fail_stage_c | fail_stage_d | fail_stage_f,
    final_status = if_else(fail_any, "FAIL", "PASS")
  )

write_csv(snp_qc_audit_all, file.path(OUT_DIR, "snp_qc_audit_all_snps.csv"))
cat(sprintf("  Wrote all-SNP QC audit ledger: %d rows\n", nrow(snp_qc_audit_all)))


# =============================================================================
# 6. Build comprehensive SNP QC summary
# =============================================================================
cat("\n-- Building SNP QC summary --\n")

g1_keep <- intersect(rownames(gt_filt), g1_ids)
g2_keep <- intersect(rownames(gt_filt), g2_ids)

snp_cr_final_all <- snp_call_rate(gt_filt)
snp_cr_g1 <- snp_call_rate(gt_filt[g1_keep, , drop = FALSE])
snp_cr_g2 <- snp_call_rate(gt_filt[g2_keep, , drop = FALSE])

snp_qc <- data.frame(
  snp_id            = colnames(gt_filt),
  call_rate         = snp_cr_final_all,
  maf               = snp_maf_final_pass,
  call_rate_g1      = snp_cr_g1,
  call_rate_g2      = snp_cr_g2,
  cr_diff_g1_g2     = abs(snp_cr_g1 - snp_cr_g2),
  mean_dp           = snp_mean_dp,
  stringsAsFactors = FALSE
) %>%
  left_join(het_stats_stage_d, by = "snp_id") %>%
  left_join(map_raw %>% select(snp_id, chromosome, position, ref, alt, qual), by = "snp_id") %>%
  mutate(
    pass_qc = TRUE,
    pass_gs = !is.na(maf) & maf >= MAF_GS
  )

# -- Dosage-based allele frequency and variance ratio (if available) --
if (dosage_available) {
  dosage_af <- colMeans(dosage_filt, na.rm = TRUE) / 2
  dosage_af[is.nan(dosage_af)] <- NA_real_
  snp_qc$dosage_af <- dosage_af
  snp_qc$dosage_maf <- pmin(dosage_af, 1 - dosage_af)
  
  dosage_var_obs <- apply(dosage_filt, 2, var, na.rm = TRUE)
  dosage_var_obs[is.nan(dosage_var_obs)] <- NA_real_
  dosage_var_exp <- 2 * dosage_af * (1 - dosage_af)
  dosage_var_ratio <- dosage_var_obs / dosage_var_exp
  dosage_var_ratio[is.nan(dosage_var_ratio) | is.infinite(dosage_var_ratio)] <- NA_real_
  
  snp_qc$dosage_var_obs   <- dosage_var_obs
  snp_qc$dosage_var_exp   <- dosage_var_exp
  snp_qc$dosage_var_ratio <- dosage_var_ratio
  
  rm(dosage_af, dosage_var_obs, dosage_var_exp, dosage_var_ratio)
}

snps_pass    <- snp_qc$snp_id
snps_pass_gs <- snp_qc %>% filter(pass_gs) %>% pull(snp_id)

cat(sprintf("  SNPs passing all QC (Stages A-F)     : %d\n", length(snps_pass)))
cat(sprintf("  SNPs passing QC + MAF >= %.2f (GS)   : %d\n", MAF_GS, length(snps_pass_gs)))
cat(sprintf("  SNPs with |CR_G1 - CR_G2| >= 0.10   : %d\n",
            sum(!is.na(snp_qc$cr_diff_g1_g2) & snp_qc$cr_diff_g1_g2 >= 0.10)))

write_csv(snp_qc, file.path(OUT_DIR, "snp_qc_summary.csv"))

rm(snp_cr_final_all, snp_cr_g1, snp_cr_g2, snp_maf_final,
   snp_maf_final_pass, snp_fail_maf_f, snp_newly_mono, snp_remove_f,
   het_stats_stage_d, snp_mean_dp)
gc()


# -----------------------------------------------------------------------------
# 7. Missingness by site x chromosome (minor scaffolds collapsed)
# -----------------------------------------------------------------------------
cat("\n-- Missingness by site x chromosome --\n")

snp_chr_display <- chr_display[map_raw$chromosome[match(colnames(gt_filt), map_raw$snp_id)]]
names(snp_chr_display) <- colnames(gt_filt)

site_levels <- sample_keep_meta %>%
  distinct(site) %>%
  arrange(site) %>%
  pull(site)
chr_display_levels <- sort(unique(snp_chr_display))

cat(sprintf("  Contig groups for missingness: %d (%d major + %s)\n",
            length(chr_display_levels), length(major_contigs),
            if (length(minor_contigs) > 0L) "1 collapsed scaffold group" else "no minor scaffolds"))

site_chr_list <- vector("list", length(site_levels) * length(chr_display_levels))
out_idx <- 0L

snps_by_chr <- split(names(snp_chr_display), snp_chr_display)

for (site_i in site_levels) {
  site_samples <- sample_keep_meta %>% filter(site == site_i) %>% pull(sample_id)
  if (length(site_samples) == 0L) next
  for (chr_i in chr_display_levels) {
    chr_snps <- snps_by_chr[[chr_i]]
    if (is.null(chr_snps) || length(chr_snps) == 0L) next
    sub_mat <- gt_filt[site_samples, chr_snps, drop = FALSE]
    out_idx <- out_idx + 1L
    site_chr_list[[out_idx]] <- data.frame(
      site = site_i,
      chromosome = chr_i,
      n_samples = nrow(sub_mat),
      n_snps = ncol(sub_mat),
      mean_call_rate = mean(!is.na(sub_mat)),
      missing_rate = mean(is.na(sub_mat))
    )
  }
}

sample_missing_by_site_chr <- if (out_idx > 0L) {
  bind_rows(site_chr_list[seq_len(out_idx)])
} else {
  data.frame(
    site = character(), chromosome = character(),
    n_samples = integer(), n_snps = integer(),
    mean_call_rate = numeric(), missing_rate = numeric()
  )
}
write_csv(sample_missing_by_site_chr, file.path(OUT_DIR, "sample_missingness_by_site_chromosome.csv"))
cat(sprintf("  Wrote site x chromosome missingness summary: %d rows\n", nrow(sample_missing_by_site_chr)))

rm(site_chr_list, site_levels, chr_display_levels, snp_chr_display, snps_by_chr)
gc()


# -----------------------------------------------------------------------------
# 7c. Build pedigree panel (strict QC + optional physical thinning)
# -----------------------------------------------------------------------------
cat("\n-- Building pedigree panel --\n")
cat(sprintf("  Pedigree panel thresholds: CR >= %.2f, MAF >= %.2f, Hobs/Hexp <= %.2f, |CR_G1-CR_G2| <= %.2f\n",
            PED_SNP_CR, PED_MAF_MIN, PED_HET_EXCESS_RATIO_MAX, PED_CR_DIFF_MAX))

ped_panel_base <- snp_qc %>%
  filter(
    pass_qc,
    !is.na(call_rate), call_rate >= PED_SNP_CR,
    !is.na(maf), maf >= PED_MAF_MIN,
    is.na(cr_diff_g1_g2) | cr_diff_g1_g2 <= PED_CR_DIFF_MAX,
    is.na(het_exp_g1) | het_exp_g1 < HEXP_MIN | (!is.na(het_excess_ratio) & het_excess_ratio <= PED_HET_EXCESS_RATIO_MAX)
  ) %>%
  arrange(chromosome, position, snp_id)

cat(sprintf("  Pedigree panel before thinning : %d SNPs\n", nrow(ped_panel_base)))

ped_panel_ids <- physical_thin_map(ped_panel_base, PED_PHYSICAL_THIN_BP)
ped_panel_map <- ped_panel_base %>%
  filter(snp_id %in% ped_panel_ids) %>%
  mutate(.ord = match(snp_id, ped_panel_ids)) %>%
  arrange(.ord) %>%
  select(-.ord)

cat(sprintf("  Pedigree panel after thinning  : %d SNPs\n", length(ped_panel_ids)))
if (PED_PHYSICAL_THIN_BP > 0) {
  cat(sprintf("  Physical thinning distance     : %d bp\n", PED_PHYSICAL_THIN_BP))
}

writeLines(ped_panel_ids, file.path(OUT_DIR, "snps_qc_pedigree_ids.txt"))
write_csv(ped_panel_map, file.path(OUT_DIR, "snps_qc_pedigree_map.csv"))

if (length(ped_panel_ids) < 5000L) {
  warning(sprintf(
    "Only %d pedigree-panel SNPs available; pedigree checks may have limited power.",
    length(ped_panel_ids)
  ), call. = FALSE)
}


# -----------------------------------------------------------------------------
# 8. Write filtered outputs (SNP map + GS SNP list)
# -----------------------------------------------------------------------------
cat("\n-- Writing filtered outputs --\n")

map_pass <- snp_qc %>%
  filter(pass_qc) %>%
  select(
    snp_id, chromosome, position, ref, alt, qual,
    maf, maf_g1, call_rate, call_rate_g1, call_rate_g2, cr_diff_g1_g2,
    mean_dp, het_obs_g1, het_exp_g1, het_excess_ratio, hwe_pval_g1,
    any_of("dosage_af"), any_of("dosage_maf"),
    any_of("dosage_var_obs"), any_of("dosage_var_exp"), any_of("dosage_var_ratio")
  )
write_csv(map_pass, file.path(OUT_DIR, "snps_qc_pass_map.csv"))
writeLines(snps_pass_gs, file.path(OUT_DIR, "snps_qc_pass_gs_ids.txt"))
cat(sprintf("  Wrote GS SNP list: %d SNPs\n", length(snps_pass_gs)))
rm(map_pass)
gc()

# -----------------------------------------------------------------------------
# 9. Pedigree verification via IBS0 (G1-G2 parent-offspring)
# -----------------------------------------------------------------------------
cat("\n-- Pedigree check: G1-G2 IBS0 verification --\n")

if (length(ped_panel_ids) == 0L) {
  warning("Pedigree panel is empty; skipping IBS0 pedigree checks.", call. = FALSE)
  check_results <- data.frame()
  pedigree_errors <- data.frame()
} else {
  gt_ped <- gt_filt[, ped_panel_ids, drop = FALSE]
  cat(sprintf("  SNPs used for IBS0 pedigree panel : %d\n", ncol(gt_ped)))
  
  g2_check_pairs <- ped_geno %>%
    filter(gen == "G2", genotype_id %in% rownames(gt_filt)) %>%
    transmute(
      offspring   = genotype_id,
      sire        = sire,
      dam         = dam,
      sire_in_vcf = !is.na(sire) & sire %in% rownames(gt_filt),
      dam_in_vcf  = !is.na(dam)  & dam  %in% rownames(gt_filt)
    ) %>%
    filter(sire_in_vcf | dam_in_vcf)
  
  cat(sprintf("  G2 offspring with >= 1 genotyped parent : %d\n", nrow(g2_check_pairs)))
  
  if (nrow(g2_check_pairs) > 0L) {
    check_results <- g2_check_pairs %>%
      mutate(
        n_markers_sire = mapply(
          function(off, par, ok) {
            if (!ok) return(NA_integer_)
            n_compared_markers(gt_ped[off, ], gt_ped[par, ])
          },
          offspring, sire, sire_in_vcf
        ),
        n_markers_dam = mapply(
          function(off, par, ok) {
            if (!ok) return(NA_integer_)
            n_compared_markers(gt_ped[off, ], gt_ped[par, ])
          },
          offspring, dam, dam_in_vcf
        ),
        ibs0_sire = mapply(
          function(off, par, ok) {
            if (!ok) return(NA_real_)
            ibs0_rate(gt_ped[off, ], gt_ped[par, ])
          },
          offspring, sire, sire_in_vcf
        ),
        ibs0_dam = mapply(
          function(off, par, ok) {
            if (!ok) return(NA_real_)
            ibs0_rate(gt_ped[off, ], gt_ped[par, ])
          },
          offspring, dam, dam_in_vcf
        )
      ) %>%
      mutate(
        po_check_sire = case_when(
          is.na(ibs0_sire)       ~ "NA",
          ibs0_sire <= IBS0_PASS ~ "PASS",
          ibs0_sire <= IBS0_WARN ~ "LOW",
          TRUE                   ~ "FAIL"
        ),
        po_check_dam = case_when(
          is.na(ibs0_dam)       ~ "NA",
          ibs0_dam <= IBS0_PASS ~ "PASS",
          ibs0_dam <= IBS0_WARN ~ "LOW",
          TRUE                  ~ "FAIL"
        )
      )
    
    cat(sprintf("\n  IBS0 results (n=%d offspring):\n", nrow(check_results)))
    cat(sprintf("    Sire: PASS=%d  LOW=%d  FAIL=%d  NA=%d\n",
                sum(check_results$po_check_sire == "PASS", na.rm = TRUE),
                sum(check_results$po_check_sire == "LOW",  na.rm = TRUE),
                sum(check_results$po_check_sire == "FAIL", na.rm = TRUE),
                sum(check_results$po_check_sire == "NA",   na.rm = TRUE)))
    cat(sprintf("    Dam : PASS=%d  LOW=%d  FAIL=%d  NA=%d\n",
                sum(check_results$po_check_dam == "PASS", na.rm = TRUE),
                sum(check_results$po_check_dam == "LOW",  na.rm = TRUE),
                sum(check_results$po_check_dam == "FAIL", na.rm = TRUE),
                sum(check_results$po_check_dam == "NA",   na.rm = TRUE)))
    
    pedigree_errors <- check_results %>%
      filter(po_check_sire == "FAIL" | po_check_dam == "FAIL")
    
    if (nrow(pedigree_errors) > 0L) {
      cat("\n  WARNING - potential pedigree errors:\n")
      print(pedigree_errors)
    }
    
    write_csv(check_results, file.path(OUT_DIR, "pedigree_check_ibs0.csv"))
    
    # 9a. Parent search for FAIL pairs
    failed_offspring <- unique(c(
      pedigree_errors %>% filter(po_check_sire == "FAIL") %>% pull(offspring),
      pedigree_errors %>% filter(po_check_dam  == "FAIL") %>% pull(offspring)
    ))
    
    if (length(failed_offspring) > 0L && length(g1_keep) > 0L) {
      cat(sprintf("\n-- Parent search: scanning %d offspring against %d G1 candidates --\n",
                  length(failed_offspring), length(g1_keep)))
      
      search_results_list <- vector("list", length(failed_offspring))
      
      for (fi in seq_along(failed_offspring)) {
        off_id <- failed_offspring[fi]
        off_gt <- gt_ped[off_id, ]
        
        off_row <- pedigree_errors %>% filter(offspring == off_id)
        sire_failed <- any(off_row$po_check_sire == "FAIL")
        dam_failed  <- any(off_row$po_check_dam  == "FAIL")
        recorded_sire <- off_row$sire[1L]
        recorded_dam  <- off_row$dam[1L]
        
        g1_ibs0 <- vapply(g1_keep, function(g1_id) {
          ibs0_rate(off_gt, gt_ped[g1_id, ])
        }, numeric(1))
        
        g1_n_markers <- vapply(g1_keep, function(g1_id) {
          n_compared_markers(off_gt, gt_ped[g1_id, ])
        }, integer(1))
        
        candidates <- data.frame(
          offspring      = off_id,
          candidate_g1   = g1_keep,
          ibs0           = g1_ibs0,
          n_markers      = g1_n_markers,
          stringsAsFactors = FALSE
        ) %>%
          filter(!is.na(ibs0)) %>%
          arrange(ibs0) %>%
          mutate(
            rank = row_number(),
            is_recorded_sire = candidate_g1 == recorded_sire,
            is_recorded_dam  = candidate_g1 == recorded_dam,
            failed_as = case_when(
              is_recorded_sire & sire_failed ~ "recorded_sire_FAIL",
              is_recorded_dam  & dam_failed  ~ "recorded_dam_FAIL",
              TRUE                           ~ ""
            )
          )
        
        recorded_ids <- c(
          if (sire_failed) recorded_sire else NULL,
          if (dam_failed)  recorded_dam  else NULL
        )
        top_n <- head(candidates, PARENT_SEARCH_TOP_N)
        recorded_rows <- candidates %>% filter(candidate_g1 %in% recorded_ids)
        search_results_list[[fi]] <- bind_rows(top_n, recorded_rows) %>%
          distinct(offspring, candidate_g1, .keep_all = TRUE) %>%
          arrange(ibs0)
      }
      
      parent_search <- bind_rows(search_results_list)
      
      if (nrow(parent_search) > 0L) {
        write_csv(parent_search, file.path(OUT_DIR, "pedigree_parent_search.csv"))
        cat(sprintf("  Wrote parent search results: %d rows for %d offspring\n",
                    nrow(parent_search), length(failed_offspring)))
        
        cat("\n  Parent search summary (top candidate vs recorded parent):\n")
        for (off_id in failed_offspring) {
          off_search <- parent_search %>% filter(offspring == off_id)
          off_errors <- pedigree_errors %>% filter(offspring == off_id)
          
          if (any(off_errors$po_check_sire == "FAIL")) {
            rec_sire <- off_errors$sire[1L]
            rec_sire_ibs0 <- off_search %>% filter(candidate_g1 == rec_sire) %>% pull(ibs0)
            rec_sire_ibs0 <- if (length(rec_sire_ibs0) > 0) sprintf("%.4f", rec_sire_ibs0[1]) else "NA"
            top_cand <- off_search %>% slice(1L)
            cat(sprintf("    %s sire: recorded=%s (IBS0=%s) | best=%s (IBS0=%.4f, rank=%d)\n",
                        off_id, rec_sire, rec_sire_ibs0,
                        top_cand$candidate_g1, top_cand$ibs0, top_cand$rank))
          }
          if (any(off_errors$po_check_dam == "FAIL")) {
            rec_dam <- off_errors$dam[1L]
            rec_dam_ibs0 <- off_search %>% filter(candidate_g1 == rec_dam) %>% pull(ibs0)
            rec_dam_ibs0 <- if (length(rec_dam_ibs0) > 0) sprintf("%.4f", rec_dam_ibs0[1]) else "NA"
            top_cand <- off_search %>% slice(1L)
            cat(sprintf("    %s dam:  recorded=%s (IBS0=%s) | best=%s (IBS0=%.4f, rank=%d)\n",
                        off_id, rec_dam, rec_dam_ibs0,
                        top_cand$candidate_g1, top_cand$ibs0, top_cand$rank))
          }
        }
      }
      
      rm(search_results_list, parent_search)
      gc()
    } else if (nrow(pedigree_errors) > 0L) {
      cat("  Parent search skipped: no G1 individuals available for comparison\n")
    }
  } else {
    cat("  Skipped - no G2 offspring with a genotyped parent in the filtered dataset\n")
    check_results <- data.frame()
    pedigree_errors <- data.frame()
  }
  
  rm(gt_ped)
  gc()
}


# -----------------------------------------------------------------------------
# 9b. Mendelian error check for complete trios
# -----------------------------------------------------------------------------
cat("\n-- Mendelian error check for complete trios --\n")

trio_results <- data.frame()
mendel_per_snp <- data.frame()
complete_trios <- data.frame()

if (length(ped_panel_ids) > 0L) {
  complete_trios <- ped_geno %>%
    filter(gen == "G2", genotype_id %in% rownames(gt_filt)) %>%
    transmute(
      offspring = genotype_id,
      sire      = sire,
      dam       = dam
    ) %>%
    filter(
      !is.na(sire) & sire %in% rownames(gt_filt),
      !is.na(dam)  & dam  %in% rownames(gt_filt)
    )
  
  cat(sprintf("  Complete trios (both parents genotyped) : %d\n", nrow(complete_trios)))
  
  if (nrow(complete_trios) > 0L) {
    gt_trio <- gt_filt[, ped_panel_ids, drop = FALSE]
    
    trio_err_list <- vector("list", nrow(complete_trios))
    snp_err_matrix <- matrix(NA, nrow = nrow(complete_trios), ncol = ncol(gt_trio),
                             dimnames = list(complete_trios$offspring, colnames(gt_trio)))
    
    for (ti in seq_len(nrow(complete_trios))) {
      off_id  <- complete_trios$offspring[ti]
      sire_id <- complete_trios$sire[ti]
      dam_id  <- complete_trios$dam[ti]
      
      errs <- mendelian_errors_trio(gt_trio[sire_id, ], gt_trio[dam_id, ], gt_trio[off_id, ])
      snp_err_matrix[ti, ] <- errs
      
      n_tested <- sum(!is.na(errs))
      n_errors <- sum(errs, na.rm = TRUE)
      
      trio_err_list[[ti]] <- data.frame(
        offspring  = off_id,
        sire      = sire_id,
        dam       = dam_id,
        n_tested  = n_tested,
        n_errors  = n_errors,
        error_rate = if (n_tested > 0) n_errors / n_tested else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    
    trio_results <- bind_rows(trio_err_list)
    
    cat(sprintf("  Median per-trio error rate  : %.4f\n", median(trio_results$error_rate, na.rm = TRUE)))
    cat(sprintf("  Mean per-trio error rate    : %.4f\n", mean(trio_results$error_rate, na.rm = TRUE)))
    cat(sprintf("  Max per-trio error rate     : %.4f (%s)\n",
                max(trio_results$error_rate, na.rm = TRUE),
                trio_results$offspring[which.max(trio_results$error_rate)]))
    
    median_err <- median(trio_results$error_rate, na.rm = TRUE)
    trio_results$flag <- ifelse(
      !is.na(trio_results$error_rate) & trio_results$error_rate > max(3 * median_err, 0.02),
      "HIGH", "OK"
    )
    n_high_trios <- sum(trio_results$flag == "HIGH", na.rm = TRUE)
    if (n_high_trios > 0L) {
      cat(sprintf("  Trios with elevated error rate : %d\n", n_high_trios))
    }
    
    write_csv(trio_results, file.path(OUT_DIR, "pedigree_mendelian_trios.csv"))
    
    snp_n_tested <- colSums(!is.na(snp_err_matrix))
    snp_n_errors <- colSums(snp_err_matrix, na.rm = TRUE)
    snp_err_rate <- ifelse(snp_n_tested > 0, snp_n_errors / snp_n_tested, NA_real_)
    
    mendel_per_snp <- data.frame(
      snp_id     = colnames(gt_trio),
      n_trios_tested = snp_n_tested,
      n_errors   = snp_n_errors,
      error_rate = snp_err_rate,
      flag       = ifelse(!is.na(snp_err_rate) & snp_err_rate > MENDEL_SNP_ERROR_THR, "HIGH", "OK"),
      stringsAsFactors = FALSE
    )
    
    n_high_snps <- sum(mendel_per_snp$flag == "HIGH", na.rm = TRUE)
    cat(sprintf("  SNPs with error rate > %.2f  : %d / %d (%.1f%%)\n",
                MENDEL_SNP_ERROR_THR, n_high_snps, nrow(mendel_per_snp),
                100 * n_high_snps / max(nrow(mendel_per_snp), 1)))
    
    write_csv(mendel_per_snp, file.path(OUT_DIR, "pedigree_mendelian_per_snp.csv"))
    
    rm(gt_trio, snp_err_matrix, trio_err_list)
    gc()
  } else {
    cat("  Skipped -- no complete trios available.\n")
  }
} else {
  cat("  Skipped -- pedigree panel is empty.\n")
}


# -----------------------------------------------------------------------------
# 9b2. Mendelian error stratified by het-ratio bin (conifer diagnostic)
# -----------------------------------------------------------------------------
cat("\n-- Mendelian error stratified by het-ratio bin (conifer diagnostic) --\n")


if (nrow(trio_results) > 0L && nrow(complete_trios) > 0L) {
  
  ped_het_ratios <- het_ratio_all_pre_d %>%
    filter(snp_id %in% ped_panel_ids, !fail_stage_d)
  
  if (nrow(ped_het_ratios) > 0L && nrow(mendel_per_snp) > 0L) {
    ped_mendel_binned <- mendel_per_snp %>%
      left_join(ped_het_ratios %>% select(snp_id, het_excess_ratio), by = "snp_id") %>%
      filter(!is.na(het_excess_ratio), n_trios_tested > 0L) %>%
      mutate(
        het_ratio_bin = cut(het_excess_ratio,
                            breaks = c(0, 0.90, 1.00, 1.10, 1.15, 1.20, 1.25, Inf),
                            labels = c("<0.90", "0.90-1.00", "1.00-1.10",
                                       "1.10-1.15", "1.15-1.20", "1.20-1.25", ">1.25"),
                            right = TRUE, include.lowest = TRUE)
      )
    
    cat("\n  Pass-QC pedigree panel: Mendelian error by het-ratio bin:\n")
    cat(sprintf("    %-12s  %6s  %10s  %10s\n", "Het ratio", "nSNPs", "Mean err", "Median err"))
    for (b in levels(ped_mendel_binned$het_ratio_bin)) {
      sub <- ped_mendel_binned %>% filter(het_ratio_bin == b)
      if (nrow(sub) > 0L) {
        cat(sprintf("    %-12s  %6d  %10.4f  %10.4f\n",
                    b, nrow(sub),
                    mean(sub$error_rate, na.rm = TRUE),
                    median(sub$error_rate, na.rm = TRUE)))
      }
    }
  }
  
  # -- Compute errors for MARGINAL SNPs (ratio threshold-1.40) --
  if (!is.null(gt_marginal_het) && ncol(gt_marginal_het) > 0L) {
    cat(sprintf("\n  Marginal SNPs (removed at %.2f, would survive at 1.40): %d SNPs\n",
                HET_EXCESS_RATIO_MAX, ncol(gt_marginal_het)))
    
    marginal_snp_err_matrix <- matrix(
      NA, nrow = nrow(complete_trios), ncol = ncol(gt_marginal_het),
      dimnames = list(complete_trios$offspring, colnames(gt_marginal_het))
    )
    
    for (ti in seq_len(nrow(complete_trios))) {
      off_id  <- complete_trios$offspring[ti]
      sire_id <- complete_trios$sire[ti]
      dam_id  <- complete_trios$dam[ti]
      
      if (all(c(off_id, sire_id, dam_id) %in% rownames(gt_marginal_het))) {
        errs <- mendelian_errors_trio(
          gt_marginal_het[sire_id, ],
          gt_marginal_het[dam_id, ],
          gt_marginal_het[off_id, ]
        )
        marginal_snp_err_matrix[ti, ] <- errs
      }
    }
    
    marg_snp_n_tested <- colSums(!is.na(marginal_snp_err_matrix))
    marg_snp_n_errors <- colSums(marginal_snp_err_matrix, na.rm = TRUE)
    marg_snp_err_rate <- ifelse(marg_snp_n_tested > 0,
                                marg_snp_n_errors / marg_snp_n_tested, NA_real_)
    
    marginal_mendel <- marginal_het_info %>%
      mutate(
        n_trios_tested = marg_snp_n_tested,
        n_errors       = marg_snp_n_errors,
        error_rate     = marg_snp_err_rate,
        het_ratio_bin  = cut(het_excess_ratio,
                             breaks = c(1.25, 1.30, 1.35, 1.40),
                             labels = c("1.25-1.30", "1.30-1.35", "1.35-1.40"),
                             right = TRUE, include.lowest = TRUE)
      )
    
    cat("\n  Marginal SNPs: Mendelian error by het-ratio sub-bin:\n")
    cat(sprintf("    %-12s  %6s  %10s  %10s  %8s\n",
                "Het ratio", "nSNPs", "Mean err", "Median err", ">5% err"))
    for (b in levels(marginal_mendel$het_ratio_bin)) {
      sub <- marginal_mendel %>% filter(het_ratio_bin == b, n_trios_tested > 0)
      if (nrow(sub) > 0L) {
        cat(sprintf("    %-12s  %6d  %10.4f  %10.4f  %8d\n",
                    b, nrow(sub),
                    mean(sub$error_rate, na.rm = TRUE),
                    median(sub$error_rate, na.rm = TRUE),
                    sum(sub$error_rate > MENDEL_SNP_ERROR_THR, na.rm = TRUE)))
      }
    }
    
    if (exists("ped_mendel_binned") && nrow(ped_mendel_binned) > 0L) {
      pass_high <- ped_mendel_binned %>%
        filter(het_excess_ratio > 1.15 & het_excess_ratio <= 1.25,
               n_trios_tested > 0)
      if (nrow(pass_high) > 0L) {
        cat(sprintf("\n  Comparison -- pass-QC SNPs with ratio 1.15-1.25:\n"))
        cat(sprintf("    n=%d  mean_err=%.4f  median_err=%.4f\n",
                    nrow(pass_high),
                    mean(pass_high$error_rate, na.rm = TRUE),
                    median(pass_high$error_rate, na.rm = TRUE)))
      }
      marginal_tested <- marginal_mendel %>% filter(n_trios_tested > 0)
      if (nrow(marginal_tested) > 0L) {
        cat(sprintf("  Marginal SNPs (ratio %.2f-1.40) overall:\n", HET_EXCESS_RATIO_MAX))
        cat(sprintf("    n=%d  mean_err=%.4f  median_err=%.4f\n",
                    nrow(marginal_tested),
                    mean(marginal_tested$error_rate, na.rm = TRUE),
                    median(marginal_tested$error_rate, na.rm = TRUE)))
      }
      
      if (nrow(pass_high) > 0L && nrow(marginal_tested) > 0L) {
        pass_mean  <- mean(pass_high$error_rate, na.rm = TRUE)
        marg_mean  <- mean(marginal_tested$error_rate, na.rm = TRUE)
        ratio_diff <- marg_mean / max(pass_mean, 1e-6)
        cat(sprintf("\n  Error rate ratio (marginal / pass-high): %.2f\n", ratio_diff))
        if (ratio_diff < 1.5) {
          cat("  -> Marginal SNPs have similar error rates to pass-QC high-ratio SNPs.\n")
          cat("    Consider relaxing HET_EXCESS_RATIO_MAX toward 1.40 to recover markers.\n")
        } else {
          cat("  -> Marginal SNPs show elevated error rates -- current threshold may be appropriate.\n")
        }
      }
    }
    
    write_csv(marginal_mendel, file.path(OUT_DIR, "het_ratio_marginal_mendel.csv"))
    
    rm(marginal_snp_err_matrix,
       marg_snp_n_tested, marg_snp_n_errors, marg_snp_err_rate, marginal_mendel)
    gc()
  } else {
    cat("  No marginal SNPs saved -- skipping marginal Mendelian analysis.\n")
  }
} else {
  cat("  Skipped -- no complete trios available for stratified analysis.\n")
}

if (exists("gt_marginal_het"))    rm(gt_marginal_het)
if (exists("marginal_het_info"))  rm(marginal_het_info)
if (exists("het_ratio_all_pre_d")) rm(het_ratio_all_pre_d)
if (exists("ped_mendel_binned"))  rm(ped_mendel_binned)
gc()


# -----------------------------------------------------------------------------
# 9c. Duplicate / identity check (all samples + G1-specific)
# -----------------------------------------------------------------------------
cat("\n-- Duplicate / identity check (all samples) --\n")

duplicate_pairs <- data.frame()

if (length(ped_panel_ids) > 0L && nrow(gt_filt) >= 2L) {
  dup_snps <- ped_panel_ids
  if (length(dup_snps) > DUPLICATE_MAX_SNPS) {
    set.seed(42)
    dup_snps <- sample(dup_snps, DUPLICATE_MAX_SNPS)
  }
  
  gt_dup <- gt_filt[, dup_snps, drop = FALSE]
  cat(sprintf("  Samples: %d  |  SNPs: %d\n", nrow(gt_dup), ncol(gt_dup)))
  
  cat("  Computing pairwise concordance matrix...\n")
  conc_mat <- pairwise_concordance(gt_dup)
  
  ut <- which(upper.tri(conc_mat), arr.ind = TRUE)
  conc_vals <- conc_mat[ut]
  
  high_conc <- !is.na(conc_vals) & conc_vals >= DUPLICATE_CONC_THR
  n_dup_pairs <- sum(high_conc)
  
  if (n_dup_pairs > 0L) {
    dup_idx <- ut[high_conc, , drop = FALSE]
    duplicate_pairs <- data.frame(
      sample_1    = rownames(conc_mat)[dup_idx[, 1]],
      sample_2    = rownames(conc_mat)[dup_idx[, 2]],
      concordance = conc_vals[high_conc],
      stringsAsFactors = FALSE
    ) %>%
      arrange(desc(concordance))
    
    duplicate_pairs <- duplicate_pairs %>%
      left_join(
        sample_keep_meta %>% select(sample_id, gen) %>% rename(sample_1 = sample_id, gen_1 = gen),
        by = "sample_1"
      ) %>%
      left_join(
        sample_keep_meta %>% select(sample_id, gen) %>% rename(sample_2 = sample_id, gen_2 = gen),
        by = "sample_2"
      )
    
    cat(sprintf("\n  WARNING -- %d sample pairs with concordance >= %.2f:\n", n_dup_pairs, DUPLICATE_CONC_THR))
    for (di in seq_len(min(n_dup_pairs, 20L))) {
      cat(sprintf("    %s x %s : %.4f  (%s/%s)\n",
                  duplicate_pairs$sample_1[di],
                  duplicate_pairs$sample_2[di],
                  duplicate_pairs$concordance[di],
                  duplicate_pairs$gen_1[di],
                  duplicate_pairs$gen_2[di]))
    }
    if (n_dup_pairs > 20L) {
      cat(sprintf("    ... plus %d more pairs (see CSV)\n", n_dup_pairs - 20L))
    }
    
    write_csv(duplicate_pairs, file.path(OUT_DIR, "duplicate_sample_pairs.csv"))
  } else {
    cat(sprintf("  No duplicate pairs detected (concordance threshold: %.2f)\n", DUPLICATE_CONC_THR))
  }
  
  g1_in_dup <- intersect(rownames(gt_dup), g1_ids)
  if (length(g1_in_dup) >= 2L) {
    g1_conc <- conc_mat[g1_in_dup, g1_in_dup]
    g1_ut_vals <- g1_conc[upper.tri(g1_conc)]
    g1_ut_vals <- g1_ut_vals[!is.na(g1_ut_vals)]
    cat(sprintf("\n  G1 pairwise concordance (n=%d founders, %d pairs):\n",
                length(g1_in_dup), length(g1_ut_vals)))
    cat(sprintf("    Mean: %.4f  Max: %.4f  Min: %.4f\n",
                mean(g1_ut_vals), max(g1_ut_vals), min(g1_ut_vals)))
    n_g1_dup <- sum(g1_ut_vals >= DUPLICATE_CONC_THR)
    if (n_g1_dup > 0L) {
      cat(sprintf("    WARNING -- %d G1 pairs above %.2f concordance (possible duplicates)\n",
                  n_g1_dup, DUPLICATE_CONC_THR))
    } else {
      cat("    No G1 duplicates detected.\n")
    }
  }
  
  rm(gt_dup, conc_mat, ut, conc_vals, high_conc)
  gc()
} else {
  cat("  Skipped -- insufficient samples or pedigree panel empty.\n")
}


# -----------------------------------------------------------------------------
# 9d. Save pass-QC genotype matrix (and optional dosage matrix)
# -----------------------------------------------------------------------------
cat("\n-- Saving pass-QC genotype matrix (RDS) --\n")

gt_pass <- gt_filt[, snps_pass, drop = FALSE]

saveRDS(gt_pass, file.path(OUT_DIR, "snps_qc_pass_012.rds"))
cat(sprintf("  Wrote 0/1/2 genotype matrix: %d samples x %d SNPs\n",
            nrow(gt_pass), ncol(gt_pass)))

if (dosage_available) {
  dosage_pass <- dosage_filt[, snps_pass, drop = FALSE]
  
  saveRDS(dosage_pass, file.path(OUT_DIR, "snps_qc_pass_dosage.rds"))
  cat(sprintf("  Wrote dosage matrix: %d samples x %d SNPs\n",
              nrow(dosage_pass), ncol(dosage_pass)))
  
  # -- Dosage diagnostic summary --
  cat("\n  Dosage diagnostic summary (pass-QC SNPs):\n")
  cat(sprintf("    Missing rate: %.2f%%\n", 100 * mean(is.na(dosage_pass))))
  cat(sprintf("    Mean dosage (non-NA): %.4f\n", mean(dosage_pass, na.rm = TRUE)))
  cat(sprintf("    SD dosage: %.4f\n", sd(as.vector(dosage_pass), na.rm = TRUE)))
  
  dosage_af_pass <- colMeans(dosage_pass, na.rm = TRUE) / 2
  hardcall_af_pass <- colMeans(gt_pass, na.rm = TRUE) / 2
  af_cor <- cor(dosage_af_pass, hardcall_af_pass, use = "complete.obs")
  cat(sprintf("    Dosage vs hard-call AF correlation: %.4f\n", af_cor))
  
  af_diff <- abs(dosage_af_pass - hardcall_af_pass)
  cat(sprintf("    Mean |dosage_AF - hardcall_AF|: %.4f\n", mean(af_diff, na.rm = TRUE)))
  cat(sprintf("    Max  |dosage_AF - hardcall_AF|: %.4f\n", max(af_diff, na.rm = TRUE)))
  cat(sprintf("    SNPs with |AF diff| > 0.05: %d\n", sum(af_diff > 0.05, na.rm = TRUE)))
  
  # ── PATCH 3: Dosage variance ratio diagnostic (v4.6.1) ────────────────────
  # v4.6.1 fix: DVR at low MAF is compressed by dosage mechanics (values
  # cluster near 0 for rare-variant sites), not by paralog signal. The
  # warning is now gated on MAF >= 0.10 where the ratio is interpretable.
  # A MAF-stratified summary replaces the undifferentiated global warning.
  cat("\n  Dosage variance ratio (observed / HWE-expected):\n")
  
  dvr_df <- data.frame(
    snp_id = snp_qc$snp_id[snp_qc$snp_id %in% snps_pass],
    dvr    = snp_qc$dosage_var_ratio[snp_qc$snp_id %in% snps_pass],
    maf    = snp_qc$maf[snp_qc$snp_id %in% snps_pass],
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(dvr))
  
  if (nrow(dvr_df) > 0L) {
    cat(sprintf("    Evaluable SNPs: %d\n", nrow(dvr_df)))
    cat(sprintf("    Mean: %.4f  Median: %.4f  SD: %.4f\n",
                mean(dvr_df$dvr), median(dvr_df$dvr), sd(dvr_df$dvr)))
    cat(sprintf("    Min: %.4f  Max: %.4f\n", min(dvr_df$dvr), max(dvr_df$dvr)))
    
    dvr_bins <- c(0, 0.50, 0.75, 0.90, 1.00, 1.10, 1.25, 1.50, Inf)
    dvr_labels <- c("<0.50", "0.50-0.75", "0.75-0.90", "0.90-1.00",
                    "1.00-1.10", "1.10-1.25", "1.25-1.50", ">1.50")
    dvr_cut <- cut(dvr_df$dvr, breaks = dvr_bins, labels = dvr_labels,
                   right = TRUE, include.lowest = TRUE)
    cat("\n    Distribution by bin (all SNPs):\n")
    for (b in dvr_labels) {
      nb <- sum(dvr_cut == b, na.rm = TRUE)
      cat(sprintf("      %-12s : %6d SNPs (%5.1f%%)\n",
                  b, nb, 100 * nb / nrow(dvr_df)))
    }
    
    # MAF-stratified DVR summary (v4.6.1)
    cat("\n    DVR by MAF bin (v4.6.1 stratified diagnostic):\n")
    cat(sprintf("    %-15s  %6s  %8s  %8s  %10s\n",
                "MAF range", "nSNPs", "Mean DVR", "Med DVR", "DVR<0.75"))
    maf_dvr_bins <- list(
      c(0, 0.05, "MAF 0-0.05"),
      c(0.05, 0.10, "MAF 0.05-0.10"),
      c(0.10, 0.20, "MAF 0.10-0.20"),
      c(0.20, 0.50, "MAF 0.20-0.50")
    )
    for (mb in maf_dvr_bins) {
      lo <- as.numeric(mb[1]); hi <- as.numeric(mb[2]); lab <- mb[3]
      sub_dvr <- dvr_df$dvr[dvr_df$maf > lo & dvr_df$maf <= hi]
      if (length(sub_dvr) > 0L) {
        n_low <- sum(sub_dvr < 0.75, na.rm = TRUE)
        cat(sprintf("    %-15s  %6d  %8.3f  %8.3f  %5d (%4.1f%%)\n",
                    lab, length(sub_dvr),
                    mean(sub_dvr, na.rm = TRUE),
                    median(sub_dvr, na.rm = TRUE),
                    n_low, 100 * n_low / length(sub_dvr)))
      }
    }
    
    # Paralog warning: only flag for common SNPs (MAF >= 0.10) where
    # DVR compression genuinely indicates collapsed paralogs
    dvr_common <- dvr_df %>% filter(maf >= 0.10)
    n_compressed_common <- sum(dvr_common$dvr < 0.75, na.rm = TRUE)
    n_common <- nrow(dvr_common)
    
    if (n_common > 0L && n_compressed_common > 0L) {
      cat(sprintf("\n    NOTE: %d / %d common SNPs (MAF >= 0.10) with DVR < 0.75 -- possible residual paralogs\n",
                  n_compressed_common, n_common))
    }
    
    # Contextualise the low-MAF compression as expected
    dvr_rare <- dvr_df %>% filter(maf < 0.10)
    n_compressed_rare <- sum(dvr_rare$dvr < 0.75, na.rm = TRUE)
    if (nrow(dvr_rare) > 0L && n_compressed_rare > 0L) {
      cat(sprintf("    INFO: %d / %d rare SNPs (MAF < 0.10) with DVR < 0.75 -- expected dosage compression at low MAF\n",
                  n_compressed_rare, nrow(dvr_rare)))
    }
  } else {
    cat("    No evaluable SNPs for dosage variance ratio.\n")
  }
  rm(dvr_df)
  if (exists("dvr_common")) rm(dvr_common)
  if (exists("dvr_rare")) rm(dvr_rare)
  
  rm(dosage_pass, dosage_af_pass, hardcall_af_pass, af_diff)
  gc()
}

rm(gt_pass)
gc()

# -----------------------------------------------------------------------------
# 9e. Export post-QC subset GDS (original content, retained samples/SNPs only)
# -----------------------------------------------------------------------------
cat("\n-- Exporting post-QC subset GDS --\n")

postqc_sample_ids <- rownames(gt_filt)

postqc_variant_ids <- map_raw$variant_id_gds[
  match(colnames(gt_filt), map_raw$snp_id)
]

if (any(is.na(postqc_variant_ids))) {
  missing_snps <- colnames(gt_filt)[is.na(postqc_variant_ids)]
  stop(sprintf(
    "Could not map %d retained SNP(s) back to original GDS variant.id. Example: %s",
    length(missing_snps),
    paste(head(missing_snps, 10), collapse = ", ")
  ))
}

if (file.exists(POSTQC_GDS_FILE)) {
  file.remove(POSTQC_GDS_FILE)
}

gds_export <- seqOpen(GDS_FILE)

seqResetFilter(gds_export, verbose = FALSE)
seqSetFilter(
  gds_export,
  sample.id  = postqc_sample_ids,
  variant.id = postqc_variant_ids,
  verbose    = FALSE
)

seqExport(
  gds_export,
  POSTQC_GDS_FILE,
  optimize = TRUE,
  digest   = TRUE,
  verbose  = TRUE
)

seqClose(gds_export)

cat(sprintf("  Wrote post-QC GDS: %s\n", POSTQC_GDS_FILE))
cat(sprintf("  Samples exported : %d\n", length(postqc_sample_ids)))
cat(sprintf("  SNPs exported    : %d\n", length(postqc_variant_ids)))

# Optional sanity check
gds_check <- seqOpen(POSTQC_GDS_FILE)

exp_sample_ids  <- seqGetData(gds_check, "sample.id")
exp_variant_ids <- seqGetData(gds_check, "variant.id")

cat(sprintf("  Verified exported GDS dimensions: %d samples x %d variants\n",
            length(exp_sample_ids), length(exp_variant_ids)))

cat(sprintf("  Sample ID match   : %s\n",
            if (setequal(exp_sample_ids, postqc_sample_ids)) "YES" else "NO"))
cat(sprintf("  Variant ID match  : %s\n",
            if (setequal(exp_variant_ids, postqc_variant_ids)) "YES" else "NO"))

seqClose(gds_check)

rm(postqc_sample_ids, postqc_variant_ids, exp_sample_ids, exp_variant_ids, gds_check)
gc()

# -----------------------------------------------------------------------------
# 10. Post-QC diagnostic summary
# -----------------------------------------------------------------------------
cat("\n======================================================================\n")
cat("  POST-QC DIAGNOSTIC SUMMARY\n")
cat("======================================================================\n\n")


postqc_n_samples <- nrow(gt_filt)
postqc_n_snps    <- length(snps_pass)

postqc_sample_cr  <- sample_call_rate(gt_filt[, snps_pass, drop = FALSE])
postqc_sample_het <- sample_het(gt_filt[, snps_pass, drop = FALSE])

cat(sprintf("  Final matrix: %d samples x %d SNPs\n", postqc_n_samples, postqc_n_snps))
cat(sprintf("  Overall missing rate: %.2f%%\n", 100 * mean(is.na(gt_filt[, snps_pass]))))

cat("\n  Per-sample call rate (post-QC):\n")
cat(sprintf("    Mean: %.4f  Median: %.4f  Min: %.4f  Max: %.4f\n",
            mean(postqc_sample_cr, na.rm = TRUE),
            median(postqc_sample_cr, na.rm = TRUE),
            min(postqc_sample_cr, na.rm = TRUE),
            max(postqc_sample_cr, na.rm = TRUE)))

cat("\n  Per-sample heterozygosity (post-QC):\n")
cat(sprintf("    Mean: %.4f  Median: %.4f  SD: %.4f\n",
            mean(postqc_sample_het, na.rm = TRUE),
            median(postqc_sample_het, na.rm = TRUE),
            sd(postqc_sample_het, na.rm = TRUE)))

postqc_snp_maf <- snp_maf(gt_filt[, snps_pass, drop = FALSE])
postqc_snp_cr  <- snp_call_rate(gt_filt[, snps_pass, drop = FALSE])
postqc_snp_het <- snp_obs_het(gt_filt[, snps_pass, drop = FALSE])

cat("\n  Per-SNP MAF (post-QC):\n")
cat(sprintf("    Mean: %.4f  Median: %.4f  Min: %.4f  Max: %.4f\n",
            mean(postqc_snp_maf, na.rm = TRUE),
            median(postqc_snp_maf, na.rm = TRUE),
            min(postqc_snp_maf, na.rm = TRUE),
            max(postqc_snp_maf, na.rm = TRUE)))

cat("\n  Per-SNP call rate (post-QC):\n")
cat(sprintf("    Mean: %.4f  Median: %.4f  Min: %.4f\n",
            mean(postqc_snp_cr, na.rm = TRUE),
            median(postqc_snp_cr, na.rm = TRUE),
            min(postqc_snp_cr, na.rm = TRUE)))

cat("\n  Per-SNP observed heterozygosity (post-QC):\n")
cat(sprintf("    Mean: %.4f  Median: %.4f  Max: %.4f\n",
            mean(postqc_snp_het, na.rm = TRUE),
            median(postqc_snp_het, na.rm = TRUE),
            max(postqc_snp_het, na.rm = TRUE)))

cat("\n  MAF spectrum (post-QC):\n")
maf_bins <- c(0.01, 0.05, 0.10, 0.20, 0.50)
for (bi in seq_along(maf_bins)) {
  lower <- if (bi == 1L) 0 else maf_bins[bi - 1L]
  upper <- maf_bins[bi]
  n_in_bin <- sum(!is.na(postqc_snp_maf) & postqc_snp_maf > lower & postqc_snp_maf <= upper)
  cat(sprintf("    MAF (%.2f, %.2f] : %d SNPs (%.1f%%)\n",
              lower, upper, n_in_bin, 100 * n_in_bin / postqc_n_snps))
}

# ── PATCH 4: Transition / transversion ratio (v4.6.1) ───────────────────────
# v4.6.1 fix: softened the warning text. Ts/Tv < 1.5 is expected for
# non-WGS data aligned to fragmented transcriptome-derived references,
# and does not necessarily indicate false-positive enrichment.
cat("\n  Transition / transversion ratio (post-QC):\n")

pass_map <- snp_qc %>% filter(pass_qc)

if (nrow(pass_map) > 0L) {
  ref_pass <- toupper(pass_map$ref)
  alt_pass <- toupper(pass_map$alt)
  
  is_transition <- (ref_pass == "A" & alt_pass == "G") |
    (ref_pass == "G" & alt_pass == "A") |
    (ref_pass == "C" & alt_pass == "T") |
    (ref_pass == "T" & alt_pass == "C")
  is_transversion <- !is_transition & nchar(ref_pass) == 1L & nchar(alt_pass) == 1L
  
  n_ts <- sum(is_transition, na.rm = TRUE)
  n_tv <- sum(is_transversion, na.rm = TRUE)
  ts_tv_ratio <- if (n_tv > 0L) n_ts / n_tv else NA_real_
  
  cat(sprintf("    Transitions    : %d\n", n_ts))
  cat(sprintf("    Transversions  : %d\n", n_tv))
  cat(sprintf("    Ts/Tv ratio    : %.3f\n", ts_tv_ratio))
  
  if (!is.na(ts_tv_ratio) && ts_tv_ratio < 1.0) {
    cat("    WARNING: Ts/Tv < 1.0 -- substantially below expectation; investigate variant call quality.\n")
  } else if (!is.na(ts_tv_ratio) && ts_tv_ratio < 1.5) {
    cat("    NOTE: Ts/Tv < 1.5 -- lower than WGS expectation (~2.0), but typical for\n")
    cat("    reduced-representation or capture data on fragmented/transcriptome references.\n")
    cat("    Not necessarily indicative of false-positive enrichment in this context.\n")
  } else if (!is.na(ts_tv_ratio) && ts_tv_ratio < 2.0) {
    cat("    Ts/Tv is moderate -- typical for whole-genome low-coverage data.\n")
  } else if (!is.na(ts_tv_ratio)) {
    cat("    Ts/Tv is consistent with high-quality SNP calls.\n")
  }
  
  rm(ref_pass, alt_pass, is_transition, is_transversion, pass_map)
}

cat("\n  Samples by generation:\n")
gen_counts <- sample_keep_meta %>%
  count(gen) %>%
  arrange(gen)
for (gi in seq_len(nrow(gen_counts))) {
  cat(sprintf("    %s : %d\n", gen_counts$gen[gi], gen_counts$n[gi]))
}

rm(postqc_sample_cr, postqc_sample_het, postqc_snp_maf, postqc_snp_cr, postqc_snp_het)
gc()


# -----------------------------------------------------------------------------
# 11. Final summary
# -----------------------------------------------------------------------------
cat("\n== FINAL QC SUMMARY ============================================\n")
cat(sprintf("  Raw samples                               : %d\n", raw_n_samples))
cat(sprintf("  Raw biallelic SNPs after hard filtering   : %d\n", raw_n_snps))
cat("  -- Staged filtering --\n")
cat(sprintf("  Stage A: Genotype masking (DP %d-%d)      : enabled\n", DP_MIN, effective_dp_max))
cat(sprintf("  Stage B: Dead/low-depth samples            : %d removed\n",
            sum(sample_qc_full$cr_flag == "DEAD", na.rm = TRUE)))
cat(sprintf("  Stage C: Site CR + QUAL + clusters + mono  : %d removed\n",
            length(snps_removed_stage_c)))
cat(sprintf("  Stage D: Excess het + HWE (G1)             : %d removed\n",
            length(snps_removed_stage_d)))
cat(sprintf("  Stage E: Sample CR < %.0f%% + het outliers  : %d removed\n",
            SAMPLE_CR_FAIL * 100, n_removed_stage_e))
cat(sprintf("  Stage F: MAF < %.2f                       : %d removed\n",
            MAF_MIN, length(snps_removed_stage_f)))
cat("  -- Final counts --\n")
cat(sprintf("  Samples retained                          : %d\n", nrow(gt_filt)))
cat(sprintf("  SNPs passing all QC (Stages A-F)          : %d\n", length(snps_pass)))
cat(sprintf("  SNPs passing QC + MAF >= %.2f (GS)        : %d\n", MAF_GS, length(snps_pass_gs)))
cat(sprintf("  Pedigree panel SNPs                       : %d\n", length(ped_panel_ids)))
cat(sprintf("  IBS0 offspring checks                     : %d\n", nrow(check_results)))
if (nrow(check_results) > 0L) {
  cat(sprintf("  Potential pedigree errors (IBS0 FAIL)     : %d\n", nrow(pedigree_errors)))
  if (nrow(pedigree_errors) > 0L) {
    cat(sprintf("  Parent search performed                   : %d offspring screened against all G1\n",
                length(unique(c(
                  pedigree_errors %>% filter(po_check_sire == "FAIL") %>% pull(offspring),
                  pedigree_errors %>% filter(po_check_dam  == "FAIL") %>% pull(offspring)
                )))))
  }
}
if (nrow(trio_results) > 0L) {
  cat(sprintf("  Mendelian error trios checked             : %d\n", nrow(trio_results)))
  cat(sprintf("  Trios with elevated error rate            : %d\n",
              sum(trio_results$flag == "HIGH", na.rm = TRUE)))
  cat(sprintf("  SNPs with Mendelian error > %.2f          : %d\n",
              MENDEL_SNP_ERROR_THR,
              sum(mendel_per_snp$flag == "HIGH", na.rm = TRUE)))
}
cat(sprintf("  Duplicate sample pairs detected           : %d\n", nrow(duplicate_pairs)))
if (dosage_available) {
  cat("  -- Dosage output --\n")
  cat("  Dosage matrix saved: snps_qc_pass_dosage.rds\n")
  cat(sprintf("  Dosage AF vs hard-call AF: included in snp_qc_summary.csv\n"))
}

cat("\nStaged filtering order:\n")
cat(sprintf("  A) Genotype-level DP masking (DP %d-%d) + GQ masking\n", DP_MIN, effective_dp_max))
cat(sprintf("  B) Remove dead (CR < %.0f%%) + low-depth (mean DP < %d) samples\n", SAMPLE_CR_DEAD * 100, SAMPLE_DP_MIN))
cat(sprintf("  C) Site filter pass 1: CR >= %.0f%% + QUAL >= %d + SNP cluster (%d in %dbp) + monomorphic\n",
            SNP_CR * 100, QUAL_MIN, SNP_CLUSTER_MAX, SNP_CLUSTER_BP))
cat(sprintf("  D) Excess het ratio (> %.2f) + HWE excess-het (p < %s) on G1\n",
            HET_EXCESS_RATIO_MAX, formatC(HWE_EXCESS_PVAL, format = "e", digits = 0)))
cat("  E) Refined sample QC: call rate + het outliers (recomputed on clean sites)\n")
cat("  F) MAF filter (recomputed on final sample set)\n")

cat("\nOutput files:\n")
cat("  qc/snps_qc_pass_012.rds\n")
if (dosage_available) cat("  qc/snps_qc_pass_dosage.rds\n")
cat("  qc/snps_qc_pass_map.csv\n")
cat("  qc/snps_qc_pass_gs_ids.txt\n")
cat("  qc/snps_qc_pedigree_ids.txt\n")
cat("  qc/snps_qc_pedigree_map.csv\n")
cat("  qc/sample_qc_summary.csv\n")
cat("  qc/sample_missingness_by_group.csv\n")
cat("  qc/sample_missingness_by_site_chromosome.csv\n")
cat("  qc/snp_qc_summary.csv\n")
cat("  qc/snp_qc_audit_all_snps.csv\n")
cat("  qc/stage_c_snp_failures.csv\n")
cat("  qc/stage_d_snp_failures.csv\n")
cat("  qc/stage_f_snp_failures.csv\n")
cat("  qc/pedigree_check_ibs0.csv\n")
cat("  qc/pedigree_parent_search.csv         (if FAIL pairs found)\n")
cat("  qc/pedigree_mendelian_trios.csv        (if complete trios exist)\n")
cat("  qc/pedigree_mendelian_per_snp.csv      (if complete trios exist)\n")
cat("  qc/het_ratio_marginal_mendel.csv       (if marginal SNPs + trios exist)\n")
cat("  qc/duplicate_sample_pairs.csv          (if duplicates detected)\n")
cat("  qc/sample_qc_hetplot.png\n")
cat("  qc/qc_report.txt\n")
cat("  qc/doug_fir_postQC.gds\n")

#cat("\nQC report written to:", LOG_FILE, "\n")
cat("Pipeline complete.\n")


# Final cleanup
rm(all_sites, ped_geno, sample_keep_meta, sample_missing_by_group, sample_missing_by_site_chr)
rm(gt_filt, g1_ids, g2_ids, g1_keep, g2_keep, snps_pass, snps_pass_gs, ped_panel_ids, ped_panel_map)
rm(sample_qc, sample_qc_full, snp_qc, snp_qc_audit_all, map_raw, check_results, pedigree_errors, ped_panel_base)
rm(contig_snp_counts, major_contigs, minor_contigs, chr_display, n_contigs)
rm(all_samples_removed, snps_removed_stage_c, snps_removed_stage_d, snps_removed_stage_f, n_removed_stage_e)
rm(stage_c_audit, stage_d_audit, stage_f_audit)
rm(sample_mean_dp_raw, effective_dp_max, global_mean_dp)
rm(trio_results, mendel_per_snp, duplicate_pairs)
if (dosage_available) { if (exists("dosage_filt")) rm(dosage_filt) }
if (exists("site_mean_dosage")) rm(site_mean_dosage)
if (exists("complete_trios")) rm(complete_trios)
if (exists("gen_counts")) rm(gen_counts)
gc()