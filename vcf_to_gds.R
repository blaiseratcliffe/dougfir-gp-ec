library(SeqArray)
library(SNPRelate)

# ==============================================================================
# Douglas-fir GS: Full-fidelity VCF -> GDS conversion + QC-derived SNPRelate GDS
# ==============================================================================
# Purpose
# -------
# This script converts a (potentially very large) VCF into an on-disk binary GDS
# representation using SeqArray, while retaining the full VCF data model:
#   - ALL INFO fields (info.import = NULL)
#   - ALL FORMAT fields (fmt.import  = NULL)
#
# It then optionally recompresses the SeqArray GDS to LZMA for smaller size and
# improved long-run disk efficiency, and finally creates a QC-filtered, SNP-only
# derived GDS for SNPRelate analyses (HWE and GRM).
#
# Why GDS?
# --------
# VCF is text-based and expensive to parse and hold in memory (especially GT
# matrices). SeqArray’s GDS stores genotypes and annotations in a compact,
# random-access format, enabling fast streaming computations without materializing
# giant in-memory objects.
#
# Key design choices (Option 1: R-native but efficient)
# -----------------------------------------------------
# 1) Parallel import with storage.option as a STRING
#    - We pass storage.option = "ZIP_RA" (a character string), NOT a
#      seqStorageOption(...) object. This avoids a known issue in some SeqArray
#      builds where parallel worker calls can error with:
#        "'arg' should be one of 'general', 'imputation'"
#
# 2) Full-fidelity import (retain INFO/FORMAT keys)
#    - info.import = NULL imports all declared INFO variables
#    - fmt.import  = NULL imports all declared FORMAT variables
#
# 3) Optional recompression to LZMA (after import)
#    - Importing with ZIP_RA is robust and fast under parallelism; then we
#      optionally recompress to LZMA to reduce file size.
#    - NOTE (SeqArray 1.46.3): seqRecompress() does NOT accept `digest=`.
#
# 4) Derived SNPRelate GDS for downstream genetics QC
#    - We apply QC thresholds (MAF/MAC/missingness) on the SeqArray object
#      without modifying the "master" GDS.
#    - Optionally restrict the derived dataset to biallelic SNPs only, which is
#      recommended for HWE/GRM assumptions.
#
# Outputs
# -------
# 1) seq_gds_full : master SeqArray GDS containing all imported INFO/FORMAT fields
# 2) snp_gds_qc   : derived SNPRelate GDS (QC-filtered; optionally biallelic SNPs)
# 3) hwe          : HWE results computed by SNPRelate on the derived SNP GDS
# 4) grm          : GRM matrix computed by SNPRelate on the derived SNP GDS
#
# Notes on data fidelity
# ----------------------
# - This script preserves INFO/FORMAT fields (keys) by importing them all.
# - Numeric precision: In Option 1, we do NOT force float64 storage during import
#   because that requires seqStorageOption(...) objects (which can trigger the
#   parallel scenario bug in some versions). In practice, this is usually fine
#   for standard VCF annotations, but extremely high-precision decimals may be
#   stored as float32 internally.
#
# Windows / OneDrive note
# -----------------------
# Your VCF is located under OneDrive. If you see slow reads or intermittent
# file-access issues (Files On-Demand), consider copying the VCF to NVMe first.
# ==============================================================================

# ---- paths (write outputs to NVMe C:) ----
vcf_path <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/UBC_071001_snps_RAW2.vcf"
out_dir  <- "C:/Users/bratclif/dfir_gds/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seq_gds_full <- file.path(out_dir, "UBC_071001_snps_RAW.full.seq.gds")
snp_gds_qc   <- file.path(out_dir, "UBC_071001_snps_QCfiltered.snp.gds")

# ---- cores ----
# Use N physical cores available for conversion and downstream computations.
n_cores <- 20L

# ---- import settings ----
# IMPORTANT: must be STRING for parallel import in this SeqArray build
storage_primary <- "ZIP_RA"

# Optional: recompress the master GDS after import (smaller + often faster to read long-term)
recompress_to_lzma <- FALSE

# Recommended: restrict derived SNPRelate dataset to biallelic SNPs for HWE/GRM assumptions
restrict_snp_gds_to_biallelic_snps <- TRUE

# ---- QC filter thresholds for derived SNP GDS ----
qc_maf          <- 0.01
qc_mac          <- 1L
qc_missing_rate <- 0.05

# ---- helper: identify biallelic SNPs from SeqArray allele strings ----
is_biallelic_snp <- function(allele_vec) {
  # allele values typically like "A,G" or "A,AT" (REF,ALT1[,ALT2...])
  # Keep only biallelic SNPs: exactly 2 alleles and both length 1
  spl <- strsplit(allele_vec, ",", fixed = TRUE)
  vapply(
    spl,
    FUN = function(a) length(a) == 2L && all(nchar(a) == 1L),
    FUN.VALUE = logical(1)
  )
}

# ---- helper: safely describe object dimensions ----
safe_dim_string <- function(x) {
  d <- dim(x)
  if (!is.null(d)) {
    return(paste(d, collapse = " x "))
  }
  paste(length(x), "(vector)")
}

# ---- helper: sanity-check critical master GDS nodes ----
validate_master_gds <- function(gds_path) {
  message("\nRunning master GDS sanity checks...")
  stopifnot(file.exists(gds_path))
  
  gds <- seqOpen(gds_path)
  on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)
  
  seqResetFilter(gds, verbose = FALSE)
  
  sample_ids  <- seqGetData(gds, "sample.id")
  variant_ids <- seqGetData(gds, "variant.id")
  
  message("  Samples : ", length(sample_ids))
  message("  Variants: ", length(variant_ids))
  
  fmt_info <- seqSummary(gds, "annotation/format", check = "none", verbose = FALSE)
  fmt_ids <- fmt_info$ID %||% character()
  
  message("  FORMAT fields: ",
          if (length(fmt_ids)) paste(fmt_ids, collapse = ", ") else "<none>")
  
  required_fmt <- c("DP", "GQ")
  missing_fmt <- setdiff(required_fmt, fmt_ids)
  if (length(missing_fmt) > 0L) {
    stop(
      "Master GDS is missing required FORMAT field(s): ",
      paste(missing_fmt, collapse = ", "),
      ". Rebuild with fmt.import = NULL."
    )
  }
  
  message("  Testing critical reads from master GDS...")
  
  gt <- tryCatch(
    seqGetData(gds, "$dosage_alt"),
    error = function(e) stop("Failed to read $dosage_alt from master GDS: ", conditionMessage(e))
  )
  message("    $dosage_alt dim: ", safe_dim_string(gt))
  
  dp <- tryCatch(
    seqGetData(gds, "annotation/format/DP"),
    error = function(e) stop("Failed to read annotation/format/DP from master GDS: ", conditionMessage(e))
  )
  message("    DP dim         : ", safe_dim_string(dp))
  
  gq <- tryCatch(
    seqGetData(gds, "annotation/format/GQ"),
    error = function(e) stop("Failed to read annotation/format/GQ from master GDS: ", conditionMessage(e))
  )
  message("    GQ dim         : ", safe_dim_string(gq))
  
  rm(gt, dp, gq)
  invisible(TRUE)
}

# simple infix helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- preflight checks ----
# Normalize the path and fail early if missing.
vcf_path <- normalizePath(vcf_path, winslash = "/", mustWork = TRUE)
stopifnot(file.exists(vcf_path))

# ==============================================================================
# Step 1 — Convert VCF -> SeqArray GDS (full-fidelity)
# ==============================================================================
# - Imports ALL INFO + ALL FORMAT fields declared in the VCF header
# - Uses parallel import on 24 cores
# - Uses ZIP_RA storage option for robustness/speed
if (file.exists(seq_gds_full)) file.remove(seq_gds_full)

t_import <- system.time({
  seqVCF2GDS(
    vcf.fn         = vcf_path,
    out.fn         = seq_gds_full,
    storage.option = storage_primary,  # STRING avoids parallel scenario bug
    info.import    = NULL,             # ALL INFO variables
    fmt.import     = NULL,             # ALL FORMAT variables
    scenario       = "general",
    parallel       = n_cores,
    raise.error    = TRUE,
    optimize       = TRUE,
    digest         = TRUE,             # supported in seqVCF2GDS; OK to keep
    verbose        = TRUE
  )
})

message("\nFull-fidelity SeqArray GDS written:\n  ", seq_gds_full)
print(t_import)

# ==============================================================================
# Step 2 — Optional recompression to LZMA (master GDS)
# ==============================================================================
# Recompressing does not change content; it only changes storage compression.
# SeqArray 1.46.3: seqRecompress() does NOT accept digest=.
if (isTRUE(recompress_to_lzma)) {
  message("\nRecompressing to LZMA (this may take a while, but is usually worth it on NVMe)...")
  seqRecompress(
    gds.fn   = seq_gds_full,
    compress = "LZMA",
    optimize = TRUE,
    verbose  = TRUE
  )
}

# ==============================================================================
# Step 3 — Sanity summary + validation of master GDS
# ==============================================================================
# Print a summary, then explicitly test the critical nodes needed downstream:
#   - $dosage_alt
#   - annotation/format/DP
#   - annotation/format/GQ
gds <- seqOpen(seq_gds_full)
on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)
print(seqSummary(gds))
seqClose(gds)
on.exit(NULL, add = FALSE)

validate_master_gds(seq_gds_full)

# ==============================================================================
# Step 4 — Create a QC-filtered, SNPRelate-compatible derived SNP GDS
# ==============================================================================
# We reopen the master, apply QC filters in-memory (filter mask), optionally
# restrict to biallelic SNPs, then export a derived SNP GDS for SNPRelate.
gds <- seqOpen(seq_gds_full)
on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)

# Apply QC filters (derived only; master file remains unchanged)
seqSetFilterCond(
  gds,
  maf          = qc_maf,
  mac          = qc_mac,
  missing.rate = qc_missing_rate,
  parallel     = n_cores,
  verbose      = TRUE
)

# Optionally restrict to biallelic SNPs for HWE/GRM assumptions
if (isTRUE(restrict_snp_gds_to_biallelic_snps)) {
  allele <- seqGetData(gds, "allele")
  vid    <- seqGetData(gds, "variant.id")
  keep   <- vid[is_biallelic_snp(allele)]
  seqSetFilter(gds, variant.id = keep, verbose = TRUE)
}

if (file.exists(snp_gds_qc)) file.remove(snp_gds_qc)
seqGDS2SNP(gds, snp_gds_qc, verbose = TRUE)

seqClose(gds)
on.exit(NULL, add = FALSE)

message("\nDerived SNP GDS written:\n  ", snp_gds_qc)

# ==============================================================================
# Step 5 — SNPRelate analyses on the derived SNP GDS
# ==============================================================================
# Compute:
# - HWE p-values (snpgdsHWE)
# - Genomic Relationship Matrix (snpgdsGRM; GCTA method)
geno <- snpgdsOpen(snp_gds_qc)
on.exit(try(snpgdsClose(geno), silent = TRUE), add = TRUE)

hwe <- snpgdsHWE(geno, with.id = TRUE)

grm <- snpgdsGRM(
  geno,
  method        = "GCTA",
  autosome.only = FALSE,   # <-- key fix for non-human / scaffolded genomes
  num.thread    = n_cores,
  useMatrix     = TRUE,
  verbose       = TRUE
)

snpgdsClose(geno)
on.exit(NULL, add = FALSE)

# ==============================================================================
# Return objects (useful when sourcing this script)
# ==============================================================================
invisible(list(
  seq_gds_full = seq_gds_full,
  snp_gds_qc   = snp_gds_qc,
  hwe          = hwe,
  grm          = grm,
  import_time  = t_import
))