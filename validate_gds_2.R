library(SeqArray)

# ==============================================================================
# Spot-check fidelity: VCF vs SeqArray GDS (INFO/FORMAT + genotype fields)
# ==============================================================================
# Purpose
# -------
# This script performs a *spot-check validation* that a SeqArray GDS retains the
# same information as the original VCF for a random subset of variants and a
# small subset of samples.
#
# It is intended as a practical integrity check after VCF -> GDS conversion:
#   - Confirms variant order alignment (by chrom/pos) between VCF and GDS
#   - Confirms selected INFO keys (e.g., DP) match between VCF and GDS
#   - Confirms selected FORMAT keys (e.g., GT/DP/GQ) match for chosen samples
#
# What it does (high level)
# -------------------------
# 1) Randomly selects `n` variant indices from the GDS's `variant.id` array.
#    - In SeqArray, `variant.id` is a sequential integer ID assigned in input
#      order during conversion.
# 2) Filters the GDS to those variants and a subset of samples.
# 3) Extracts the same variant rows from the VCF by *streaming the VCF* until the
#    requested record indices are reached (no random-access index required).
# 4) Compares:
#    - chrom/pos alignment
#    - INFO fields (as strings)
#    - FORMAT fields for selected samples (as strings)
# 5) Returns a list of mismatches ("issues") if any are detected.
#
# Why this approach works
# -----------------------
# - It avoids loading the entire VCF or full genotype matrices into RAM.
# - It does not require tabix indexing.
# - It provides strong evidence of conversion fidelity without needing a full
#   re-export and file diff (which is expensive for large datasets).
#
# Important limitations / caveats
# -------------------------------
# - This is a *spot check*, not a proof of full file equality.
# - INFO numeric formatting may differ between VCF text and GDS numeric storage
#   (e.g., "3" vs "3.0", or floating rounding). This script compares values as
#   character strings; if you see benign formatting differences, you can upgrade
#   the comparison to numeric tolerances for specific keys.
# - FORMAT fields in VCF are strings; in GDS they may be stored as integers,
#   floats, matrices, or lists (variable length). This script assumes the common
#   case where selected FORMAT fields are scalar per-sample per-variant (matrix).
#
# Outputs
# -------
# The main result is a list:
#   - picked_variants: integer variant indices checked
#   - compared_samples: sample IDs used for comparison
#   - n_issues: number of mismatches found
#   - issues: list of mismatch records with context
#
# Usage
# -----
# Edit the paths in the "run spotcheck" block at the bottom, and adjust:
#   - n (number of variants to spot-check)
#   - info_keys / fmt_keys (fields to compare)
#   - sample_sel (sample indices to compare)
# ==============================================================================


# ------------------------------------------------------------------------------
# extract_vcf_by_variant_index()
# ------------------------------------------------------------------------------
# Stream-read a VCF and extract a set of variant rows by their 1-based record
# index (i.e., the order of variants in the file after the header).
#
# This is used to fetch the exact same variant records that correspond to the
# GDS `variant.id` indices (which are assigned sequentially in import order).
#
# Parameters
# ----------
# vcf_path : character
#   Path to VCF (plain .vcf or gz .vcf.gz). Uses gzfile() which can also read
#   plain text VCF.
# variant_idx : integer vector
#   1-based indices of variant records in the VCF body (post-header).
# info_keys, fmt_keys : character vectors (currently unused as filters)
#   Included for future extensions (e.g., only parse specific keys).
# sample_idx : integer vector
#   Indices of samples (columns) to extract FORMAT values for. Default: first 5.
#
# Returns
# -------
# Named list keyed by variant index (character):
#   each element contains:
#     chrom, pos, id, ref, alt
#     info: named list of INFO key->value or flag->TRUE
#     fmt_names: FORMAT field names (from column 9)
#     fmt_mat: character matrix [selected_samples x fmt_names]
# ------------------------------------------------------------------------------
extract_vcf_by_variant_index <- function(vcf_path, variant_idx, info_keys = character(), fmt_keys = character(),
                                         sample_idx = 1:5) {
  
  variant_idx <- sort(unique(as.integer(variant_idx)))
  max_i <- max(variant_idx)
  
  con <- gzfile(vcf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  # Read header + sample names
  samp_ids <- character()
  repeat {
    ln <- readLines(con, n = 1L)
    if (length(ln) == 0L) stop("Unexpected EOF while reading header.")
    if (startsWith(ln, "#CHROM\t")) {
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(parts) > 9) samp_ids <- parts[10:length(parts)]
      break
    }
  }
  
  # Storage
  out <- vector("list", length(variant_idx))
  names(out) <- as.character(variant_idx)
  want <- variant_idx
  want_ptr <- 1L
  
  i <- 0L
  repeat {
    lines <- readLines(con, n = 10000L)
    if (length(lines) == 0L) break
    
    for (ln in lines) {
      if (ln == "" || ln[1] == "#") next
      i <- i + 1L
      if (i < want[want_ptr]) next
      if (i > max_i) break
      
      if (i == want[want_ptr]) {
        cols <- strsplit(ln, "\t", fixed = TRUE)[[1]]
        info <- cols[8]
        fmt  <- cols[9]
        samp <- cols[10:length(cols)]
        
        # INFO parse
        info_map <- list()
        if (!identical(info, ".")) {
          toks <- strsplit(info, ";", fixed = TRUE)[[1]]
          for (t in toks) {
            if (t == "") next
            if (grepl("=", t, fixed = TRUE)) {
              kv <- strsplit(t, "=", fixed = TRUE)[[1]]
              info_map[[kv[1]]] <- kv[2]
            } else {
              info_map[[t]] <- TRUE  # INFO flag
            }
          }
        }
        
        # FORMAT parse for selected samples
        fmt_names <- strsplit(fmt, ":", fixed = TRUE)[[1]]
        samp_sel <- samp[sample_idx]
        samp_names_sel <- samp_ids[sample_idx]
        
        fmt_mat <- matrix(
          NA_character_,
          nrow = length(samp_sel),
          ncol = length(fmt_names),
          dimnames = list(samp_names_sel, fmt_names)
        )
        
        for (r in seq_along(samp_sel)) {
          vals <- strsplit(samp_sel[r], ":", fixed = TRUE)[[1]]
          fmt_mat[r, seq_len(min(length(vals), length(fmt_names)))] <-
            vals[seq_len(min(length(vals), length(fmt_names)))]
        }
        
        out[[as.character(i)]] <- list(
          chrom = cols[1],
          pos   = as.integer(cols[2]),
          id    = cols[3],
          ref   = cols[4],
          alt   = cols[5],
          info = info_map,
          fmt_names = fmt_names,
          fmt_mat = fmt_mat
        )
        
        want_ptr <- want_ptr + 1L
        if (want_ptr > length(want)) break
      }
    }
    
    if (i >= max_i || want_ptr > length(want)) break
  }
  
  out
}


# ------------------------------------------------------------------------------
# spotcheck_vcf_vs_gds()
# ------------------------------------------------------------------------------
# Randomly sample `n` variants from a SeqArray GDS and compare selected INFO/FORMAT
# fields between:
#   - the GDS (annotation/info/* and annotation/format/*)
#   - the original VCF (parsed as strings)
#
# Parameters
# ----------
# vcf_path : character
#   Path to original VCF.
# seq_gds_path : character
#   Path to SeqArray GDS created from that VCF.
# n : integer
#   Number of variants to spot-check.
# info_keys : character vector
#   INFO keys to compare (e.g., c("DP")).
# fmt_keys : character vector
#   FORMAT keys to compare (e.g., c("GT","DP","GQ")).
# sample_sel : integer vector
#   Sample indices to compare (e.g., 1:5). This selects the same sample indices
#   from VCF and GDS.
# seed : integer
#   Seed for reproducible variant sampling.
#
# Returns
# -------
# list:
#   picked_variants: integer vector of sampled variant indices
#   compared_samples: sample IDs in the GDS subset used for comparison
#   n_issues: count of mismatches
#   issues: list of mismatch records (variant_id, problem, context)
# ------------------------------------------------------------------------------
spotcheck_vcf_vs_gds <- function(vcf_path, seq_gds_path, n = 50,
                                 info_keys = c("DP"),
                                 fmt_keys  = c("GT","DP","GQ"),
                                 sample_sel = 1:5,
                                 seed = 1) {
  
  gds <- seqOpen(seq_gds_path)
  on.exit(seqClose(gds), add = TRUE)
  
  variant.id <- seqGetData(gds, "variant.id")
  if (n > length(variant.id)) stop("n is larger than the number of variants in the GDS.")
  
  set.seed(seed)
  pick <- sort(sample(variant.id, n))
  
  # Filter GDS to the same subset of variants and samples
  seqSetFilter(gds, variant.id = pick, sample.sel = sample_sel, verbose = FALSE)
  
  # Positional sanity: confirm VCF record index aligns with GDS variant order
  gds_chr <- as.character(seqGetData(gds, "chromosome"))
  gds_pos <- seqGetData(gds, "position")
  
  # Pull selected INFO and FORMAT from GDS.
  # Note: Some FORMAT variables can be list-like; we assume scalar matrix fields here.
  gds_info <- lapply(info_keys, function(k) {
    tryCatch(seqGetData(gds, paste0("annotation/info/", k)), error = function(e) NULL)
  })
  names(gds_info) <- info_keys
  
  gds_fmt <- lapply(fmt_keys, function(k) {
    tryCatch(seqGetData(gds, paste0("annotation/format/", k)), error = function(e) NULL)
  })
  names(gds_fmt) <- fmt_keys
  
  # Extract matching rows from the VCF by variant index
  vcf_recs <- extract_vcf_by_variant_index(
    vcf_path, pick,
    info_keys = info_keys,
    fmt_keys  = fmt_keys,
    sample_idx = sample_sel
  )
  
  issues <- list()
  
  for (j in seq_along(pick)) {
    idx <- pick[j]
    v <- vcf_recs[[as.character(idx)]]
    
    if (is.null(v)) {
      issues[[length(issues) + 1]] <- list(
        variant_id = idx,
        problem = "VCF record not found during scan"
      )
      next
    }
    
    # Align check (chrom/pos)
    if (!identical(v$chrom, gds_chr[j]) || !identical(v$pos, gds_pos[j])) {
      issues[[length(issues) + 1]] <- list(
        variant_id = idx,
        problem = "chrom/pos mismatch between VCF and GDS",
        vcf_chrom = v$chrom, vcf_pos = v$pos,
        gds_chrom = gds_chr[j], gds_pos = gds_pos[j]
      )
      next
    }
    
    # INFO checks (string compare; numeric formatting may differ)
    for (k in info_keys) {
      gv <- gds_info[[k]]
      if (is.null(gv)) {
        issues[[length(issues) + 1]] <- list(
          variant_id = idx,
          problem = paste("GDS missing INFO", k)
        )
        next
      }
      
      vcf_val <- v$info[[k]]
      gds_val <- gv[j]
      
      if (!is.null(vcf_val)) {
        # INFO flag (present/absent)
        if (isTRUE(vcf_val) &&
            !(isTRUE(gds_val) || identical(as.character(gds_val), "TRUE"))) {
          issues[[length(issues) + 1]] <- list(
            variant_id = idx,
            problem = paste("INFO flag mismatch", k),
            vcf = vcf_val, gds = gds_val
          )
        }
        
        # INFO key=value
        if (!isTRUE(vcf_val)) {
          if (!identical(as.character(vcf_val), as.character(gds_val))) {
            issues[[length(issues) + 1]] <- list(
              variant_id = idx,
              problem = paste("INFO value mismatch", k),
              vcf = as.character(vcf_val),
              gds = as.character(gds_val)
            )
          }
        }
      }
    }
    
    # FORMAT checks (selected samples)
    for (k in fmt_keys) {
      gv <- gds_fmt[[k]]
      if (is.null(gv)) {
        issues[[length(issues) + 1]] <- list(
          variant_id = idx,
          problem = paste("GDS missing FORMAT", k)
        )
        next
      }
      if (!k %in% colnames(v$fmt_mat)) next
      
      vcf_vals <- v$fmt_mat[, k]
      
      # Common case: gv is a matrix [variant x sample]
      gds_vals <- gv[j, , drop = TRUE]
      
      if (!all(as.character(vcf_vals) == as.character(gds_vals))) {
        issues[[length(issues) + 1]] <- list(
          variant_id = idx,
          problem = paste("FORMAT mismatch", k),
          vcf = as.character(vcf_vals),
          gds = as.character(gds_vals)
        )
      }
    }
  }
  
  list(
    picked_variants   = pick,
    compared_samples  = seqGetData(gds, "sample.id"),
    n_issues          = length(issues),
    issues            = issues
  )
}

# ==============================================================================
# Run spotcheck (edit paths)
# ==============================================================================
res <- spotcheck_vcf_vs_gds(
  vcf_path     = "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/UBC_071001_snps_RAW.vcf",
  seq_gds_path = "C:/Users/bratclif/dfir_gds/UBC_071001_snps_RAW.full.seq.gds", 
  n = 50,
  info_keys = c("DP"),
  fmt_keys  = c("GT","DP","GQ"),
  sample_sel = 1:5,
  seed = 123
)

str(res$n_issues)
if (res$n_issues) str(res$issues[[1]])