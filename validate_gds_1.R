# ==============================================================================
# Validate VCF Header Completeness + Validate SeqArray GDS Retains INFO/FORMAT
# ==============================================================================
# What this script does
# ---------------------
# This script answers two questions:
#
#   1) Is the VCF "well-declared"?
#      - i.e., does every INFO key used in the VCF body appear in a header line
#        like: ##INFO=<ID=...>
#      - and does every FORMAT key used in the VCF body appear in a header line
#        like: ##FORMAT=<ID=...>
#
#   2) Did SeqArray retain the INFO/FORMAT keys during VCF -> GDS conversion?
#      - i.e., for every INFO key used in the VCF body, does the GDS contain a
#        corresponding node under: annotation/info/<KEY>
#      - and for every FORMAT key used in the VCF body, does the GDS contain a
#        corresponding node under: annotation/format/<KEY>
#
# Why this matters
# ----------------
# SeqArray imports INFO/FORMAT based on declarations in the VCF header. If the VCF
# body contains keys not declared in the header, tools may drop/ignore them.
# This script detects those cases and also confirms that the converted GDS
# includes the same keys that are actually used in the VCF.
#
# Performance notes
# -----------------
# - The VCF scan streams through the file and only parses columns 8 (INFO)
#   and 9 (FORMAT). It does NOT load genotypes into memory.
# - Progress reporting uses base::format() rather than sprintf("%,d") because
#   some Windows builds do not support the thousands-separator flag.
#
# Outputs
# -------
# - hdr:  INFO/FORMAT keys declared in the header
# - used: INFO/FORMAT keys actually used in the VCF body
# - gds:  INFO/FORMAT keys present in the SeqArray GDS
# - cmp:  set-differences showing problems (if any)
#
# Typical interpretation
# ----------------------
# - cmp$vcf_*_used_not_declared non-empty  => VCF header is incomplete / invalid-ish
# - cmp$gds_missing_*_used non-empty       => conversion likely dropped used keys
# - cmp$gds_missing_*_declared non-empty   => usually OK (declared but never used)
# ==============================================================================

library(SeqArray)

# ------------------------------------------------------------------------------
# parse_vcf_header_tags()
# ------------------------------------------------------------------------------
# Read the VCF header and extract:
# - all declared INFO IDs (##INFO=<ID=...>)
# - all declared FORMAT IDs (##FORMAT=<ID=...>)
# - sample IDs from the #CHROM line (optional convenience)
#
# Parameters
# ----------
# vcf_path : character
#   Path to VCF (plain .vcf or gz .vcf.gz). This uses gzfile(); gzfile() will
#   also read plain text files.
#
# Returns
# -------
# list with elements:
#   info_decl: character vector of INFO IDs declared in header
#   fmt_decl : character vector of FORMAT IDs declared in header
#   samp_ids : character vector of sample IDs (may be empty if none)
# ------------------------------------------------------------------------------
parse_vcf_header_tags <- function(vcf_path) {
  con <- gzfile(vcf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  info_decl <- character()
  fmt_decl  <- character()
  samp_ids  <- character()
  
  repeat {
    ln <- readLines(con, n = 1L)
    if (length(ln) == 0L) stop("Unexpected EOF while reading VCF header.")
    
    if (startsWith(ln, "##INFO=<ID=")) {
      id <- sub("^##INFO=<ID=([^,>]+).*$", "\\1", ln)
      info_decl <- c(info_decl, id)
      
    } else if (startsWith(ln, "##FORMAT=<ID=")) {
      id <- sub("^##FORMAT=<ID=([^,>]+).*$", "\\1", ln)
      fmt_decl <- c(fmt_decl, id)
      
    } else if (startsWith(ln, "#CHROM\t")) {
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(parts) > 9) samp_ids <- parts[10:length(parts)]
      break
    }
  }
  
  list(
    info_decl = unique(info_decl),
    fmt_decl  = unique(fmt_decl),
    samp_ids  = samp_ids
  )
}

# ------------------------------------------------------------------------------
# scan_vcf_used_tags()
# ------------------------------------------------------------------------------
# Stream-scan the VCF body and collect:
# - INFO keys used in column 8
# - FORMAT keys used in column 9
#
# This avoids parsing genotypes entirely. It only reads columns 8 and 9 from
# each variant row (still splits the row by tabs).
#
# Parameters
# ----------
# vcf_path : character
#   Path to VCF (plain .vcf or gz .vcf.gz).
# report_every : integer
#   Print a progress message every N variants.
#
# Returns
# -------
# list with elements:
#   n_variants: integer, number of variant lines scanned
#   info_used : character vector, INFO keys observed in body
#   fmt_used  : character vector, FORMAT keys observed in body
# ------------------------------------------------------------------------------
scan_vcf_used_tags <- function(vcf_path, report_every = 50000L) {
  con <- gzfile(vcf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  # Skip header
  repeat {
    ln <- readLines(con, n = 1L)
    if (length(ln) == 0L) stop("Unexpected EOF while skipping header.")
    if (startsWith(ln, "#CHROM\t")) break
  }
  
  used_info <- new.env(parent = emptyenv())
  used_fmt  <- new.env(parent = emptyenv())
  
  n <- 0L
  repeat {
    lines <- readLines(con, n = 10000L)
    if (length(lines) == 0L) break
    
    for (ln in lines) {
      if (ln == "" || ln[1] == "#") next
      n <- n + 1L
      
      cols <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(cols) < 9L) next
      
      # INFO column (8)
      info <- cols[8]
      if (!identical(info, ".")) {
        toks <- strsplit(info, ";", fixed = TRUE)[[1]]
        for (t in toks) {
          if (t == "") next
          key <- sub("=.*$", "", t)
          used_info[[key]] <- TRUE
        }
      }
      
      # FORMAT column (9)
      fmt <- cols[9]
      if (!identical(fmt, ".")) {
        keys <- strsplit(fmt, ":", fixed = TRUE)[[1]]
        for (k in keys) used_fmt[[k]] <- TRUE
      }
    }
    
    # Windows-safe progress formatting (avoid sprintf("%,d"))
    if (n %% report_every == 0L) {
      message(sprintf(
        "Scanned %s variants...",
        format(n, big.mark = ",", scientific = FALSE)
      ))
    }
  }
  
  list(
    n_variants = n,
    info_used  = sort(ls(used_info)),
    fmt_used   = sort(ls(used_fmt))
  )
}

# ------------------------------------------------------------------------------
# gds_tags()
# ------------------------------------------------------------------------------
# List INFO/FORMAT variables present in a SeqArray GDS file by inspecting:
#   annotation/info/*
#   annotation/format/*
#
# Parameters
# ----------
# seq_gds_path : character
#   Path to SeqArray GDS file.
#
# Returns
# -------
# list with:
#   info_gds: character vector of INFO keys present in GDS
#   fmt_gds : character vector of FORMAT keys present in GDS
# ------------------------------------------------------------------------------
gds_tags <- function(seq_gds_path) {
  gds <- seqOpen(seq_gds_path)
  on.exit(seqClose(gds), add = TRUE)
  
  info_df <- seqSummary(gds, "annotation/info")
  fmt_df  <- seqSummary(gds, "annotation/format")
  
  list(
    info_gds = sort(info_df$ID),
    fmt_gds  = sort(fmt_df$ID)
  )
}

# ------------------------------------------------------------------------------
# compare_tag_sets()
# ------------------------------------------------------------------------------
# Compute the important set differences:
# - used-in-body but not declared-in-header  => VCF header completeness failure
# - used-in-body but missing-in-GDS          => conversion retention failure
# - declared-in-header but missing-in-GDS    => usually OK if never used
#
# Returns a list of character vectors.
# ------------------------------------------------------------------------------
compare_tag_sets <- function(header, used, gds) {
  list(
    # VCF validity checks
    vcf_info_used_not_declared = setdiff(used$info_used, header$info_decl),
    vcf_fmt_used_not_declared  = setdiff(used$fmt_used,  header$fmt_decl),
    
    # Conversion fidelity checks (focus on USED keys)
    gds_missing_info_used = setdiff(used$info_used, gds$info_gds),
    gds_missing_fmt_used  = setdiff(used$fmt_used,  gds$fmt_gds),
    
    # Declared but not present in GDS (often OK if never used in body)
    gds_missing_info_declared = setdiff(header$info_decl, gds$info_gds),
    gds_missing_fmt_declared  = setdiff(header$fmt_decl,  gds$fmt_gds)
  )
}

# ------------------------------------------------------------------------------
# Example run (edit paths)
# ------------------------------------------------------------------------------
vcf_path <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/UBC_071001_snps_RAW.vcf"
seq_gds  <- "C:/Users/bratclif/dfir_gds/UBC_071001_snps_RAW.full.seq.gds"

hdr  <- parse_vcf_header_tags(vcf_path)
used <- scan_vcf_used_tags(vcf_path, report_every = 50000L)
gds  <- gds_tags(seq_gds)

cmp <- compare_tag_sets(hdr, used, gds)

cat("\n=== VCF body vs header ===\n")
cat("Variants scanned:", format(used$n_variants, big.mark = ",", scientific = FALSE), "\n")
cat("INFO used but NOT declared in header:", length(cmp$vcf_info_used_not_declared), "\n")
cat("FORMAT used but NOT declared in header:", length(cmp$vcf_fmt_used_not_declared), "\n")

if (length(cmp$vcf_info_used_not_declared)) print(cmp$vcf_info_used_not_declared)
if (length(cmp$vcf_fmt_used_not_declared))  print(cmp$vcf_fmt_used_not_declared)

cat("\n=== USED tags missing in GDS (this would indicate loss) ===\n")
cat("INFO used but missing in GDS:", length(cmp$gds_missing_info_used), "\n")
cat("FORMAT used but missing in GDS:", length(cmp$gds_missing_fmt_used), "\n")

if (length(cmp$gds_missing_info_used)) print(cmp$gds_missing_info_used)
if (length(cmp$gds_missing_fmt_used))  print(cmp$gds_missing_fmt_used)

cat("\n=== Declared-in-header tags missing in GDS (often OK if unused) ===\n")
cat("INFO declared but missing in GDS:", length(cmp$gds_missing_info_declared), "\n")
cat("FORMAT declared but missing in GDS:", length(cmp$gds_missing_fmt_declared), "\n")