###############################################################################
## Douglas-fir Genomic Selection Project (UBC_071001)
## SNP Imputation & GRM Construction Pipeline
## Methods: (1) Mean imputation, (2) LinkImputeR, (3) gtImputation, (4) KGD
##
## Input:  Post-QC GDS file from the v4.6 SeqArray/SNPRelate pipeline
## Output: Imputed genotype matrices + KGD GRM for downstream BGLR models
##         Timestamped log file (imputation_log_YYYYMMDD_HHMMSS.txt)
##
## IMPORTANT
## =========
## This script assumes doug_fir_postQC.gds is already the fully filtered,
## post-QC truth set. No additional biological/QC/paralog/HWE/HWD filtering
## is performed here. All post-QC SNPs are carried forward into imputation and
## GRM construction.
##
## A genotype-level DP mask is applied before imputation:
##   - low mask:  DP < 3
##   - high mask: DP > 66
##
## Masked genotypes are set to missing in memory and in the exported VCFs used
## by LinkImputeR and gtImputation.
##
## The only remaining locus handling is arithmetic: loci with zero or undefined
## G1 base-population variance are zero-weighted in GRM construction so the
## matrix math remains well-defined. This is not an additional QC screen.
###############################################################################
##
## INSTALLATION GUIDE
## ==================
##
## This script requires R packages, external software, and specific directory
## layout. Follow the steps below before first run.
##
##
## 1. R PACKAGES
## -------------
## Install from CRAN:
##
##   install.packages(c("vcfR", "ggplot2", "dplyr", "tidyr", "Matrix",
##                       "parallel"))
##
##   Note: vcfR is only used to read imputed VCF outputs from LinkImputeR
##   and gtImputation. If you only run mean imputation and/or KGD, vcfR
##   is not needed. The primary input is read via SeqArray (GDS format).
##
## Install from Bioconductor:
##
##   if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
##   BiocManager::install(c("SeqArray", "SNPRelate"))
##
## Install AGHmatrix (for GRM utilities, used in downstream models):
##
##   install.packages("AGHmatrix")
##
##
## 2. JAVA (required for LinkImputeR)
## -----------------------------------
## LinkImputeR requires Java 8 or later. No admin/system install needed —
## Java can be run from a portable extracted zip.
##
## PORTABLE (no admin — recommended):
##   1. Go to https://adoptium.net/temurin/releases/
##   2. Select: OS = Windows, Architecture = x64, Package = JDK, Type = .zip
##   3. Download and extract to your tools directory, e.g.:
##        D:/OneDrive - NRCan RNCan/gs/doug-fir/tools/jdk/
##   4. The java binary will be at something like:
##        tools/jdk/jdk-21.0.5+11/bin/java.exe
##   5. Set java_path in the config section below to point at it.
##      (javac and jar for building LinkImputeR are in the same bin/ folder.)
##
## WITH ADMIN (system install):
##   Windows: winget install EclipseAdoptium.Temurin.21.JDK
##   macOS:   brew install openjdk
##   Linux:   sudo apt install default-jdk
##   Then java_path = "java" (uses system PATH).
##
## Verify:  java -version  (should show 1.8+ / 8+)
##
##
## 3. PYTHON 3 (required for gtImputation)
## ----------------------------------------
## PORTABLE (no admin):
##   1. Go to https://www.python.org/downloads/windows/
##   2. Download "Windows embeddable package (64-bit)" (.zip, not installer)
##   3. Extract to your tools directory, e.g.:
##        D:/OneDrive - NRCan RNCan/gs/doug-fir/tools/python/
##   4. Enable pip in the embeddable package:
##      a. Edit python312._pth (or similar) in the extracted folder
##         and uncomment the "import site" line
##      b. Download get-pip.py:
##           curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
##      c. Run:  tools\python\python.exe get-pip.py
##   5. Install packages:
##        tools\python\python.exe -m pip install minisom numpy scipy
##   6. Set python_path in the config section below.
##
## WITH ADMIN (system install):
##   Windows: winget install Python.Python.3.12
##            Then: pip install minisom numpy scipy
##   macOS:   brew install python
##   Linux:   sudo apt install python3 python3-pip
##            Then: pip3 install minisom numpy scipy
##
## If pip is not available on a system install:
##   python -m ensurepip --upgrade
##   python -m pip install minisom numpy scipy
##
##
## 4. LinkImputeR (Java tool — genotype imputation for non-model organisms)
## ------------------------------------------------------------------------
## Reference: Money et al. (2017) BMC Genomics 18:523
##            Money et al. (2015) G3 5:2383-2390
##
## Source code: https://github.com/danielmoney/LinkImputeR
## Latest version: 1.2.4 (April 2020)
##
## The original download site (cultivatingdiversity.org) is no longer
## available. The GitHub repo contains source code only (no pre-built
## JAR) — it must be built from source.
##
## Build from source:
##
##   Prerequisites: JDK 8+ and Apache Ant (no admin needed for either)
##
##   PORTABLE (no admin):
##     JDK:  Already extracted per Section 2 above (javac/jar in bin/)
##     Ant:  Download the .zip from https://ant.apache.org/bindownload.cgi
##           Extract to e.g. tools/apache-ant-1.10.15/
##           No install needed — ant.bat is in the bin/ folder.
##
##   WITH ADMIN:
##     Windows: winget install EclipseAdoptium.Temurin.21.JDK
##              winget install Apache.Ant
##     macOS:   brew install openjdk ant
##     Linux:   sudo apt install default-jdk ant
##
##   Build steps:
##     git clone https://github.com/danielmoney/LinkImputeR.git
##     cd LinkImputeR
##
##     # If using portable JDK + Ant on Windows:
##     set JAVA_HOME=D:\...\tools\jdk\jdk-21.0.5+11
##     D:\...\tools\apache-ant-1.10.15\bin\ant.bat jar
##
##     # If using system-installed JDK + Ant:
##     ant jar
##
##   This produces dist/LinkImputeR.jar. If the project uses NetBeans:
##     ant -f nbproject/build-impl.xml jar
##
##   If no build.xml exists, compile manually (works with just the JDK):
##     mkdir build
##     tools\jdk\jdk-21.0.5+11\bin\javac -d build -sourcepath src src\Executable\LinkImputeR.java
##     tools\jdk\jdk-21.0.5+11\bin\jar cfe LinkImputeR.jar Executable.LinkImputeR -C build .
##
##   Place the built JAR at:
##     D:/OneDrive - NRCan RNCan/gs/doug-fir/tools/LinkImputeR.jar
##
##   NOTE: Since v1.1.3, no lib/ directory is needed — just the single
##   JAR file is sufficient.
##
## Test (using portable Java):
##   "D:/..../tools/jdk/jdk-21.0.5+11/bin/java" -cp "D:/.../tools/LinkImputeR/LinkImputeR.jar;D:/.../tools/LinkImputeR/lib/*" Executable.LinkImputeR
##   (Should print usage info, not an error)
##
## Memory: For ~100K SNPs x ~1,600 samples, allocate 16-32 GB:
##   java -Xmx24G -cp "LinkImputeR.jar;lib/*" Executable.LinkImputeR ...
##
## Note on freeBayes VCF format:
##   freeBayes uses RO (ref counts) and AO (alt counts) instead of the
##   standard AD field. This script automatically patches the exported
##   VCF to add an AD field computed as "RO,AO" per sample, so
##   LinkImputeR can parse allele depths using its standard AD reader.
##   No manual VCF conversion is needed.
##
##
## 5. gtImputation (Python tool — SOM-based imputation for non-model organisms)
## -----------------------------------------------------------------------------
## Reference: Mora-Márquez et al. (2024) Molecular Ecology Resources e13992
##
## Download:
##   git clone https://github.com/GGFHF/gtImputation.git
##
##   Place the cloned repo at:
##     D:/OneDrive - NRCan RNCan/gs/doug-fir/tools/gtImputation/
##
##   The script expects to find gtImputation.py at one of:
##     tools/gtImputation/Package/gtImputation.py
##     tools/gtImputation/gtImputation.py
##
## Python dependencies (install if not already):
##   pip install minisom numpy scipy
##
## Test:
##   python "D:/.../tools/gtImputation/Package/gtImputation.py" --help
##
## Performance warning:
##   gtImputation trains a separate SOM for each missing genotype.
##   With ~100K SNPs x ~1,600 samples and ~10% missingness, this means
##   hundreds of thousands of small neural networks. Expect hours to days
##   of runtime. Consider running overnight or on a compute node.
##   The script prints an estimate of missing genotype count before starting.
##
##
## 6. KGD (R code — depth-adjusted GRM from GBS read counts)
## ----------------------------------------------------------
## Reference: Dodds et al. (2015) BMC Genomics 16:1047
##            Dodds et al. (2019) G3 9:3239-3247
##
## Download:
##   git clone https://github.com/AgResearch/KGD.git
##
##   Place the cloned repo at:
##     D:/OneDrive - NRCan RNCan/gs/doug-fir/tools/KGD/
##
##   The script sources:
##     tools/KGD/GBS-Chip-Gmatrix.R   (main KGD functions)
##
## No compilation needed — KGD is pure R (with optional Rcpp acceleration).
##
## This pipeline builds the allele count (RA) matrix directly from the
## RO/AO fields extracted via SeqArray from the GDS file — no VCF export
## or vcf2ra.py conversion is needed for KGD.
##
##
## DIRECTORY LAYOUT (expected after setup)
## ----------------------------------------
##
##   D:/OneDrive - NRCan RNCan/gs/doug-fir/
##   ├── data/
##   │   └── qc/
##   │       ├── doug_fir_postQC.gds          <-- your post-QC GDS (primary input)
##   │       └── imputation/                  <-- created by this script
##   │           ├── imputation_log_*.txt             <-- timestamped run log
##   │           ├── doug_fir_postQC_export.vcf.gz
##   │           ├── doug_fir_postQC_export_DPmask3_66.vcf.gz
##   │           ├── doug_fir_postQC_export_DPmask3_66_AD.vcf.gz
##   │           ├── linkimputeR/
##   │           ├── gtImputation/
##   │           ├── kgd/
##   │           ├── generation_assignment.rds
##   │           ├── allele_freq_g1.rds
##   │           ├── allele_freq_g1_reads.rds
##   │           ├── geno_mean_imputed.rds
##   │           ├── GRM_mean_imputed.rds
##   │           ├── geno_linkimputeR.rds
##   │           ├── geno_linkimputeR_probs.rds
##   │           ├── GRM_linkimputeR.rds
##   │           ├── geno_gtImputation.rds
##   │           ├── GRM_gtImputation.rds
##   │           ├── GRM_kgd_depth_adjusted.rds
##   │           ├── KGD_full_result.rds
##   │           └── diagnostic_*.png
##   └── tools/
##       ├── jdk/
##       │   └── jdk-21.0.5+11/
##       │       └── bin/
##       │           ├── java.exe
##       │           ├── javac.exe
##       │           └── jar.exe
##       ├── python/
##       │   ├── python.exe
##       │   └── Lib/site-packages/
##       ├── apache-ant-1.10.15/
##       │   └── bin/ant.bat
##       ├── LinkImputeR.jar
##       ├── gtImputation/
##       │   └── Package/
##       │       └── gtImputation.py
##       └── KGD/
##           └── GBS-Chip-Gmatrix.R
##
##
## QUICK START
## -----------
##   1. Install R packages (see section 1 above)
##   2. Extract portable JDK and Python zips into tools/ (sections 2-3)
##      OR use system-installed Java/Python if you have admin
##   3. pip install minisom numpy scipy
##   4. Download/clone LinkImputeR, gtImputation, KGD into tools/
##   5. Set run_* flags below to choose which methods to execute
##   6. Adjust gds_file path to point to your post-QC GDS
##   7. Source this script:  source("doug_fir_imputation_pipeline.R")
##   8. Check the timestamped log in imputation/ for diagnostics,
##      warnings, section timings, and the end-of-run summary
##
###############################################################################

SCRIPT_VERSION <- "1.1.0"  # Added structured logging

# ===========================================================================
# 0. CONFIGURATION
# ===========================================================================

## ---- Paths ----
base_dir   <- "D:/OneDrive - NRCan RNCan/gs/doug-fir"
data_dir   <- file.path(base_dir, "data")
qc_dir     <- file.path(data_dir, "qc")
imp_dir    <- file.path(qc_dir, "imputation")

# Post-QC GDS file (authoritative already-filtered input)
gds_file   <- file.path(qc_dir, "doug_fir_postQC.gds")

# External tool paths — adjust to your machine
java_path       <- file.path(base_dir, "tools/jdk/jdk-25.0.2+10/bin/java")
# java_path     <- "java"

python_path     <- file.path("C:/Users/bratclif/AppData/Local/Programs/Python/Python313/python.exe")
# python_path   <- "python"

linkimputeR_jar <- file.path(base_dir, "tools/LinkImputeR/LinkImputeR.jar")
kgd_dir         <- file.path(base_dir, "tools/KGD")
gtimp_dir       <- file.path(base_dir, "tools/gtImputation")

## ---- Methods to run ----
run_mean_imputation <- FALSE
run_linkimputeR     <- FALSE
run_gtImputation    <- TRUE
run_kgd             <- FALSE

# Number of cores for parallel operations
ncores <- parallel::detectCores() - 2

## ---- Java classpath separator (OS-aware) ----
## Windows uses ";" and Unix/macOS use ":" in Java -cp arguments.
cp_sep <- if (.Platform$OS.type == "windows") ";" else ":"

## ---- Java heap size (GB) ----
## For ~106K post-QC SNPs x ~1,600 samples, 16 GB is sufficient.
## Increase if LinkImputeR runs out of memory on larger datasets.
java_heap_gb <- 16

## ---- Genotype-level DP mask (applied before imputation) ----
DP_MASK_LOW  <- 3
DP_MASK_HIGH <- 66

# VCF export paths — only created if LinkImputeR or gtImputation are enabled
vcf_export <- file.path(imp_dir, "doug_fir_postQC_export.vcf.gz")
vcf_export_masked <- file.path(
  imp_dir,
  sprintf("doug_fir_postQC_export_DPmask%d_%d.vcf.gz", DP_MASK_LOW, DP_MASK_HIGH)
)
vcf_export_ad <- file.path(
  imp_dir,
  sprintf("doug_fir_postQC_export_DPmask%d_%d_AD.vcf.gz", DP_MASK_LOW, DP_MASK_HIGH)
)

## ---- Generation assignment (G1 vs G2) ----
## G1 parents are at Lost, Adam, and Fleet trial sites.
## G2 progeny are at Jordan and BigTree trial sites.
## G0 founders are NOT genotyped.
##
## Allele frequencies for GRM centering are computed from G1 individuals only
## (base population proxy). GRMs are built on all individuals (G1 + G2) but
## centered on G1 frequencies. No additional SNP/sample QC filtering is
## performed in this script.

PED_GENO_FILE  <- file.path(data_dir, "pedigree_genotyped.csv")
SAMPLE_QC_FILE <- file.path(qc_dir, "sample_qc_summary.csv")

## ---- Create output directory ----
dir.create(imp_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Clean stale GRM objects if sourcing repeatedly in one session ----
for (obj in c("G_mean", "G_li", "G_gt", "G_kgd")) {
  if (exists(obj, inherits = FALSE)) rm(list = obj)
}


# ===========================================================================
# LOGGING UTILITIES
# ===========================================================================
#
# Lightweight structured logging: ISO-8601 timestamps, severity levels,
# section tags, dual output (console + persistent log file), and a warning
# accumulator that replays at end-of-run.
#
# Levels: INFO, WARN, ERROR, DEBUG, TIMING
# ===========================================================================

LOG_FILE <- file.path(imp_dir, sprintf("imputation_log_%s.txt",
                                       format(Sys.time(), "%Y%m%d_%H%M%S")))
LOG_WARNINGS  <- character(0)
LOG_ERRORS    <- character(0)
SECTION_TIMES <- list()

## Null-coalescing operator -----------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

## Core logging function --------------------------------------------------
log_msg <- function(msg, level = "INFO", section = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  prefix <- if (!is.null(section)) {
    sprintf("[%s] [%-6s] [%s]", ts, level, section)
  } else {
    sprintf("[%s] [%-6s]", ts, level)
  }
  line <- paste(prefix, msg)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
  if (level == "WARN") {
    LOG_WARNINGS <<- c(LOG_WARNINGS, paste0("[", section %||% "global", "] ", msg))
  }
  if (level == "ERROR") {
    LOG_ERRORS <<- c(LOG_ERRORS, paste0("[", section %||% "global", "] ", msg))
  }
  invisible(NULL)
}

## Section timer wrapper --------------------------------------------------
log_section <- function(name, expr) {
  log_msg(sprintf("========== START: %s ==========", name), "INFO", name)
  t0 <- proc.time()
  result <- tryCatch(force(expr), error = function(e) {
    elapsed <- (proc.time() - t0)["elapsed"]
    log_msg(sprintf("FAILED after %.1f sec: %s", elapsed, conditionMessage(e)),
            "ERROR", name)
    stop(e)
  })
  elapsed <- (proc.time() - t0)["elapsed"]
  SECTION_TIMES[[name]] <<- elapsed
  log_msg(sprintf("========== END: %s (%.1f sec) ==========", name, elapsed),
          "TIMING", name)
  invisible(result)
}

## Genotype matrix snapshot -----------------------------------------------
log_geno_snapshot <- function(mat, label, generation_vec = NULL, section = NULL) {
  n_snp <- nrow(mat); n_samp <- ncol(mat)
  miss_pct <- 100 * mean(is.na(mat))
  log_msg(sprintf("%s: %d SNPs x %d samples, overall missing = %.2f%%",
                  label, n_snp, n_samp, miss_pct), "INFO", section)
  if (!is.null(generation_vec) && length(generation_vec) == n_samp) {
    for (gen in c("G1", "G2")) {
      idx <- which(generation_vec == gen)
      if (length(idx) > 0) {
        gen_miss <- 100 * mean(is.na(mat[, idx, drop = FALSE]))
        log_msg(sprintf("  %s (%d samples): missing = %.2f%%",
                        gen, length(idx), gen_miss), "INFO", section)
      }
    }
    na_gen <- which(is.na(generation_vec))
    if (length(na_gen) > 0) {
      log_msg(sprintf("  Unassigned (%d samples): missing = %.2f%%",
                      length(na_gen),
                      100 * mean(is.na(mat[, na_gen, drop = FALSE]))),
              "INFO", section)
    }
  }
}

## GRM summary with eigenvalue diagnostics --------------------------------
log_grm_summary <- function(G, name, generation_vec = NULL, section = NULL) {
  d <- diag(G)
  o <- G[lower.tri(G)]
  log_msg(sprintf("GRM [%s]: dim = %d x %d", name, nrow(G), ncol(G)),
          "INFO", section)
  log_msg(sprintf("  Diagonal:  mean=%.4f  sd=%.4f  range=[%.4f, %.4f]",
                  mean(d), sd(d), min(d), max(d)), "INFO", section)
  log_msg(sprintf("  Off-diag:  mean=%.4f  sd=%.4f  range=[%.4f, %.4f]",
                  mean(o), sd(o), min(o), max(o)), "INFO", section)
  eig_vals <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  n_neg <- sum(eig_vals < 0)
  log_msg(sprintf("  Eigenvalues: min=%.6f  max=%.4f  n_negative=%d",
                  min(eig_vals), max(eig_vals), n_neg), "INFO", section)
  if (n_neg > 0) {
    log_msg(sprintf("  Most negative eigenvalue: %.6e", min(eig_vals)),
            "WARN", section)
  }
  if (!is.null(generation_vec) && length(generation_vec) == nrow(G)) {
    is_g1_idx <- which(generation_vec == "G1")
    is_g2_idx <- which(generation_vec == "G2")
    if (length(is_g1_idx) > 1) {
      d_g1 <- diag(G)[is_g1_idx]
      o_g1 <- G[is_g1_idx, is_g1_idx]; o_g1 <- o_g1[lower.tri(o_g1)]
      log_msg(sprintf("  G1 block: diag mean=%.4f sd=%.4f  off-diag mean=%.4f",
                      mean(d_g1), sd(d_g1), mean(o_g1)), "INFO", section)
    }
    if (length(is_g2_idx) > 1) {
      d_g2 <- diag(G)[is_g2_idx]
      o_g2 <- G[is_g2_idx, is_g2_idx]; o_g2 <- o_g2[lower.tri(o_g2)]
      log_msg(sprintf("  G2 block: diag mean=%.4f sd=%.4f  off-diag mean=%.4f",
                      mean(d_g2), sd(d_g2), mean(o_g2)), "INFO", section)
    }
    if (length(is_g1_idx) > 0 && length(is_g2_idx) > 0) {
      o_cross <- as.vector(G[is_g1_idx, is_g2_idx])
      log_msg(sprintf("  G1xG2 cross-block: mean=%.4f  range=[%.4f, %.4f]",
                      mean(o_cross), min(o_cross), max(o_cross)),
              "INFO", section)
    }
  }
}

## External tool runner with stderr capture -------------------------------
run_external <- function(cmd, args = character(0), label = "external",
                         section = NULL, capture_stderr = TRUE) {
  full_cmd <- paste(c(shQuote(cmd), args), collapse = " ")
  log_msg(sprintf("Executing [%s]: %s", label, full_cmd), "INFO", section)
  t0 <- proc.time()
  if (capture_stderr) {
    stderr_f <- tempfile(fileext = ".stderr.log")
    ret <- system2(cmd, args, stdout = "", stderr = stderr_f)
    err_lines <- readLines(stderr_f, warn = FALSE)
    unlink(stderr_f)
  } else {
    ret <- system(full_cmd, intern = FALSE)
    err_lines <- character(0)
  }
  elapsed <- (proc.time() - t0)["elapsed"]
  log_msg(sprintf("[%s] Return code: %d  Elapsed: %.1f sec",
                  label, ret, elapsed),
          ifelse(ret == 0, "INFO", "WARN"), section)
  if (length(err_lines) > 0) {
    n_show <- min(30, length(err_lines))
    for (el in tail(err_lines, n_show)) {
      log_msg(paste("  stderr:", el), "DEBUG", section)
    }
    if (length(err_lines) > n_show) {
      log_msg(sprintf("  ... (%d additional stderr lines suppressed)",
                      length(err_lines) - n_show), "DEBUG", section)
    }
  }
  invisible(ret)
}

## File size reporter -----------------------------------------------------
log_file_info <- function(filepath, label = NULL, section = NULL) {
  if (file.exists(filepath)) {
    sz <- file.size(filepath)
    sz_str <- if (sz > 1e9) sprintf("%.2f GB", sz / 1e9)
    else if (sz > 1e6) sprintf("%.2f MB", sz / 1e6)
    else if (sz > 1e3) sprintf("%.1f KB", sz / 1e3)
    else sprintf("%d bytes", sz)
    nm <- label %||% basename(filepath)
    log_msg(sprintf("File [%s]: %s (%s)", nm, filepath, sz_str),
            "INFO", section)
  } else {
    log_msg(sprintf("File NOT FOUND: %s", filepath), "WARN", section)
  }
}

## End-of-run summary -----------------------------------------------------
log_run_summary <- function() {
  log_msg("", "INFO")
  log_msg("############################################################", "INFO")
  log_msg("##               END-OF-RUN SUMMARY                       ##", "INFO")
  log_msg("############################################################", "INFO")
  if (length(SECTION_TIMES) > 0) {
    log_msg("", "INFO")
    log_msg("--- Section Timings ---", "TIMING")
    total_sec <- 0
    for (nm in names(SECTION_TIMES)) {
      s <- SECTION_TIMES[[nm]]
      total_sec <- total_sec + s
      if (s >= 3600) {
        log_msg(sprintf("  %-30s %7.1f sec  (%.2f hr)", nm, s, s / 3600), "TIMING")
      } else if (s >= 60) {
        log_msg(sprintf("  %-30s %7.1f sec  (%.1f min)", nm, s, s / 60), "TIMING")
      } else {
        log_msg(sprintf("  %-30s %7.1f sec", nm, s), "TIMING")
      }
    }
    log_msg(sprintf("  %-30s %7.1f sec  (%.1f min)", "TOTAL",
                    total_sec, total_sec / 60), "TIMING")
  }
  ## Freeze warning/error lists before replay to prevent feedback loop
  ## (log_msg at WARN level appends to LOG_WARNINGS, so replaying with
  ## WARN would contaminate the list during iteration)
  frozen_warnings <- LOG_WARNINGS
  frozen_errors   <- LOG_ERRORS
  
  if (length(frozen_warnings) > 0) {
    log_msg("", "INFO")
    log_msg(sprintf("--- %d Warnings Accumulated ---", length(frozen_warnings)), "INFO")
    for (w in frozen_warnings) log_msg(paste("  *", w), "INFO")
  } else {
    log_msg("No warnings accumulated.", "INFO")
  }
  if (length(frozen_errors) > 0) {
    log_msg("", "INFO")
    log_msg(sprintf("--- %d Errors Encountered ---", length(frozen_errors)), "INFO")
    for (e in frozen_errors) log_msg(paste("  *", e), "INFO")
  }
  log_msg("", "INFO")
  log_msg(sprintf("Log file: %s", LOG_FILE), "INFO")
  log_msg(sprintf("Completed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "INFO")
}


# ===========================================================================
# SESSION METADATA (logged once at startup)
# ===========================================================================

log_msg("=== Douglas-fir Imputation Pipeline ===")
log_msg(sprintf("Script version: %s", SCRIPT_VERSION))
log_msg(sprintf("Log file: %s", LOG_FILE))
log_msg(sprintf("Start time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
log_msg(sprintf("R version: %s", R.version.string))
log_msg(sprintf("Platform: %s", R.version$platform))
log_msg(sprintf("OS type: %s", .Platform$OS.type))
log_msg(sprintf("Hostname: %s", Sys.info()["nodename"]))
log_msg(sprintf("User: %s", Sys.info()["user"]))
log_msg(sprintf("Working directory: %s", getwd()))
log_msg(sprintf("Cores available: %d (using %d)", parallel::detectCores(), ncores))

## Log package versions
.log_pkg_version <- function(pkg) {
  v <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "not installed")
  log_msg(sprintf("  %-15s %s", pkg, v))
}
log_msg("R package versions:")
for (pkg in c("SeqArray", "SNPRelate", "vcfR", "Matrix", "ggplot2",
              "dplyr", "tidyr", "AGHmatrix", "BGLR")) {
  .log_pkg_version(pkg)
}
rm(.log_pkg_version)

## Log external tool versions (non-fatal if tools are unavailable)
if (run_linkimputeR || run_gtImputation) {
  java_ver <- tryCatch({
    v <- system2(java_path, "-version", stdout = TRUE, stderr = TRUE)
    paste(v, collapse = " | ")
  }, error = function(e) "unavailable")
  log_msg(sprintf("Java version: %s", java_ver))
}
if (run_gtImputation) {
  python_ver <- tryCatch({
    v <- system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
    paste(v, collapse = " | ")
  }, error = function(e) "unavailable")
  log_msg(sprintf("Python version: %s", python_ver))
}

## Log configuration
log_msg("Configuration:")
log_msg(sprintf("  GDS input: %s", gds_file))
log_msg(sprintf("  Output dir: %s", imp_dir))
log_msg(sprintf("  DP mask: low < %d, high > %d", DP_MASK_LOW, DP_MASK_HIGH))
log_msg(sprintf("  Java heap: %d GB", java_heap_gb))
log_msg("Methods selected:")
log_msg(sprintf("  [%s] Mean imputation (baseline)", ifelse(run_mean_imputation, "X", " ")))
log_msg(sprintf("  [%s] LinkImputeR", ifelse(run_linkimputeR, "X", " ")))
log_msg(sprintf("  [%s] gtImputation (SOM)", ifelse(run_gtImputation, "X", " ")))
log_msg(sprintf("  [%s] KGD depth-adjusted GRM", ifelse(run_kgd, "X", " ")))

## Log input file metadata
if (file.exists(gds_file)) {
  log_file_info(gds_file, "PostQC GDS")
  gds_md5 <- tryCatch(tools::md5sum(gds_file), error = function(e) NA_character_)
  if (!is.na(gds_md5)) log_msg(sprintf("  MD5: %s", gds_md5))
} else {
  log_msg(sprintf("FATAL: GDS file not found: %s", gds_file), "ERROR")
  stop("GDS file not found: ", gds_file)
}

log_msg("")


# ===========================================================================
# 1. LOAD POST-QC DATA
# ===========================================================================

log_section("Load Post-QC Data", {
  
  SEC <- "Load"
  
  library(SeqArray)
  library(Matrix)
  
  if (run_linkimputeR || run_gtImputation) {
    library(vcfR)
  }
  
  gds <- seqOpen(gds_file, readonly = TRUE)
  log_msg("GDS file opened successfully.", "INFO", SEC)
  
  ## Verify standard REF/ALT allele ordering (freeBayes VCF convention).
  .allele_check <- seqGetData(gds, "allele")
  stopifnot(
    "GDS allele field must contain REF,ALT pairs (comma-separated)" =
      all(grepl(",", head(.allele_check, 100)))
  )
  log_msg("Allele field verified: REF,ALT comma-separated pairs.", "INFO", SEC)
  rm(.allele_check)
  
  sample_ids <- seqGetData(gds, "sample.id")
  variant_ids <- seqGetData(gds, "variant.id")
  snp_chrom   <- seqGetData(gds, "chromosome")
  snp_pos     <- seqGetData(gds, "position")
  snp_ids     <- paste0(snp_chrom, "_", snp_pos)
  
  geno_t <- tryCatch(
    seqGetData(gds, "$dosage_alt"),
    error = function(e) {
      log_msg("$dosage_alt not available; computing as 2 - $dosage (ref allele).",
              "WARN", SEC)
      2L - seqGetData(gds, "$dosage")
    }
  )
  
  geno <- t(geno_t)
  rm(geno_t)
  rownames(geno) <- snp_ids
  colnames(geno) <- sample_ids
  
  n_samples <- ncol(geno)
  n_snps    <- nrow(geno)
  
  ## Input sanity checks
  ## Post-QC expectations: ~1,599 samples (after sample-level QC from ~1,872
  ## genotyped), ~106K SNPs (MAF >= 0.01 post-QC truth set). Bounds are
  ## deliberately wide to catch wrong-project GDS files without triggering
  ## on normal run-to-run variation from threshold changes.
  if (n_samples < 500 || n_samples > 2500) {
    log_msg(sprintf("Unexpected sample count: %d. Expected ~1,600 for UBC_071001.",
                    n_samples), "WARN", SEC)
  }
  if (n_snps < 5000 || n_snps > 200000) {
    log_msg(sprintf("Unexpected SNP count: %d. Expected ~100K-110K post-QC SNPs (MAF >= 0.01).",
                    n_snps), "WARN", SEC)
  }
  
  miss_rate <- sum(is.na(geno)) / length(geno)
  log_msg(sprintf("Loaded: %d samples, %d SNPs from GDS", n_samples, n_snps),
          "INFO", SEC)
  log_msg(sprintf("Genotype matrix: %d SNPs x %d samples", n_snps, n_samples),
          "INFO", SEC)
  log_msg(sprintf("Overall missing rate (pre-DP-mask): %.2f%%",
                  miss_rate * 100), "INFO", SEC)
  
  ## ---- Extract allele read depths ----
  log_msg("Extracting allele read depths (RO/AO)...", "INFO", SEC)
  
  ro_raw <- seqGetData(gds, "annotation/format/RO")
  ao_raw <- seqGetData(gds, "annotation/format/AO")
  
  extract_format_matrix <- function(raw, n_snps, n_samples, snp_ids, sample_ids) {
    if (is.matrix(raw)) {
      m <- t(raw)
    } else if (is.list(raw) && !is.null(raw$data)) {
      if (!is.null(raw$length) && !all(raw$length == 1L)) {
        stop("FORMAT field has variable-length entries (Number != 1). ",
             "Cannot safely reshape into a simple [SNP x sample] matrix.\n",
             "  Unique $length values: ", paste(unique(raw$length), collapse = ", "))
      }
      expected_len <- as.numeric(n_samples) * as.numeric(n_snps)
      if (length(raw$data) != expected_len) {
        stop("FORMAT $data length (", length(raw$data), ") does not match ",
             "expected n_samples * n_snps (", expected_len, ").\n",
             "  n_samples=", n_samples, ", n_snps=", n_snps)
      }
      m <- matrix(raw$data, nrow = n_samples, ncol = n_snps)
      m <- t(m)
    } else {
      stop("Unexpected return format from seqGetData for FORMAT field")
    }
    stopifnot(
      "Extracted matrix has wrong number of rows (SNPs)" = nrow(m) == n_snps,
      "Extracted matrix has wrong number of columns (samples)" = ncol(m) == n_samples
    )
    dimnames(m) <- list(snp_ids, sample_ids)
    m
  }
  
  ro <- extract_format_matrix(ro_raw, n_snps, n_samples, snp_ids, sample_ids)
  ao <- extract_format_matrix(ao_raw, n_snps, n_samples, snp_ids, sample_ids)
  storage.mode(ro) <- "numeric"
  storage.mode(ao) <- "numeric"
  rm(ro_raw, ao_raw)
  
  log_msg(sprintf("Allele depths extracted: RO and AO matrices (%d x %d)",
                  n_snps, n_samples), "INFO", SEC)
  
  ## Depth summary (informative for GBS data quality assessment)
  total_depth <- ro + ao
  td_nonzero  <- total_depth[total_depth > 0 & !is.na(total_depth)]
  if (length(td_nonzero) > 0) {
    log_msg(sprintf("Read depth (non-zero cells): median=%.1f  mean=%.1f  Q25=%.0f  Q75=%.0f  max=%d",
                    median(td_nonzero), mean(td_nonzero),
                    quantile(td_nonzero, 0.25), quantile(td_nonzero, 0.75),
                    max(td_nonzero)), "INFO", SEC)
  }
  rm(total_depth, td_nonzero)
  
  ## ---- Extract genotype DP and apply genotype-level DP mask ----
  log_msg("Extracting genotype DP matrix for masking...", "INFO", SEC)
  dp_raw <- seqGetData(gds, "annotation/format/DP")
  dp <- extract_format_matrix(dp_raw, n_snps, n_samples, snp_ids, sample_ids)
  storage.mode(dp) <- "numeric"
  rm(dp_raw)
  
  mask_dp_low  <- !is.na(dp) & dp < DP_MASK_LOW
  mask_dp_high <- !is.na(dp) & dp > DP_MASK_HIGH
  mask_dp_any  <- mask_dp_low | mask_dp_high
  
  n_mask_low  <- sum(mask_dp_low  & !is.na(geno))
  n_mask_high <- sum(mask_dp_high & !is.na(geno))
  n_mask_any  <- sum(mask_dp_any  & !is.na(geno))
  
  log_msg(sprintf("Genotype DP mask: low (DP < %d) -> %d calls masked",
                  DP_MASK_LOW, n_mask_low), "INFO", SEC)
  log_msg(sprintf("Genotype DP mask: high (DP > %d) -> %d calls masked",
                  DP_MASK_HIGH, n_mask_high), "INFO", SEC)
  log_msg(sprintf("Total genotype calls masked by DP: %d (%.2f%% of non-NA calls)",
                  n_mask_any,
                  100 * n_mask_any / sum(!is.na(geno))), "INFO", SEC)
  
  if (n_mask_any > 0) {
    geno[mask_dp_any] <- NA_integer_
    ro[mask_dp_any]   <- NA_real_
    ao[mask_dp_any]   <- NA_real_
  }
  
  log_msg(sprintf("Overall missing rate after DP masking: %.2f%%",
                  100 * mean(is.na(geno))), "INFO", SEC)
  
  ## Free DP matrix and mask intermediates
  rm(dp, mask_dp_low, mask_dp_high)
  gc(verbose = FALSE)
  
  seqClose(gds)
  log_msg("GDS file closed.", "INFO", SEC)
  
  ## ---- Assign generations (G1 vs G2) ----
  log_msg("Assigning generations...", "INFO", SEC)
  
  generation <- rep(NA_character_, n_samples)
  names(generation) <- sample_ids
  
  if (file.exists(PED_GENO_FILE)) {
    ped_geno <- read.csv(PED_GENO_FILE, stringsAsFactors = FALSE)
    ped_geno[] <- lapply(ped_geno, trimws)
    
    ped_match <- match(sample_ids, ped_geno$genotype_id)
    matched_gen <- ped_geno$gen[ped_match]
    
    generation[!is.na(matched_gen) & matched_gen == "G1"] <- "G1"
    generation[!is.na(matched_gen) & matched_gen == "G2"] <- "G2"
    
    n_matched <- sum(!is.na(ped_match))
    log_msg(sprintf("Generation source: pedigree_genotyped.csv (matched %d of %d samples)",
                    n_matched, n_samples), "INFO", SEC)
    
    rm(ped_geno, ped_match, matched_gen, n_matched)
    
  } else if (file.exists(SAMPLE_QC_FILE)) {
    sample_qc <- read.csv(SAMPLE_QC_FILE, stringsAsFactors = FALSE)
    
    qc_match <- match(sample_ids, sample_qc$sample_id)
    matched_gen <- sample_qc$gen[qc_match]
    
    generation[!is.na(matched_gen) & matched_gen == "G1"] <- "G1"
    generation[!is.na(matched_gen) & matched_gen == "G2"] <- "G2"
    
    n_matched <- sum(!is.na(qc_match))
    log_msg(sprintf("Generation source: sample_qc_summary.csv (matched %d of %d samples)",
                    n_matched, n_samples), "INFO", SEC)
    
    rm(sample_qc, qc_match, matched_gen, n_matched)
    
  } else {
    log_msg("FATAL: Neither pedigree_genotyped.csv nor sample_qc_summary.csv found.",
            "ERROR", SEC)
    stop("FATAL: Neither pedigree_genotyped.csv nor sample_qc_summary.csv found.\n",
         "  Checked:\n",
         "    ", PED_GENO_FILE, "\n",
         "    ", SAMPLE_QC_FILE, "\n",
         "  Run the QC pipeline first, or place pedigree_genotyped.csv in the data directory.")
  }
  
  is_g1 <- which(generation == "G1")
  is_g2 <- which(generation == "G2")
  is_na <- which(is.na(generation))
  
  log_msg(sprintf("Generation assignment: G1=%d  G2=%d  unassigned=%d",
                  length(is_g1), length(is_g2), length(is_na)), "INFO", SEC)
  
  if (length(is_na) > 0) {
    log_msg(sprintf("WARNING: %d samples could not be assigned to G1 or G2.",
                    length(is_na)), "WARN", SEC)
    log_msg("Unassigned samples EXCLUDED from AF estimation, INCLUDED in GRM.",
            "WARN", SEC)
    if (length(is_na) <= 10) {
      log_msg(sprintf("  Unassigned IDs: %s", paste(sample_ids[is_na], collapse = ", ")),
              "INFO", SEC)
    } else {
      log_msg(sprintf("  First 10 unassigned IDs: %s",
                      paste(head(sample_ids[is_na], 10), collapse = ", ")),
              "INFO", SEC)
    }
  }
  if (length(is_g1) == 0) {
    log_msg("FATAL: No G1 samples identified.", "ERROR", SEC)
    stop("FATAL: No G1 samples identified. Cannot compute base population allele frequencies.")
  }
  
  ## Log genotype snapshot AFTER DP masking AND generation assignment
  log_geno_snapshot(geno, "Post-DP-mask genotypes", generation, SEC)
  
  ## ---- Compute G1-only allele frequencies (base population) ----
  log_msg("Computing G1-only allele frequencies (base population)...", "INFO", SEC)
  
  p_g1 <- rowMeans(geno[, is_g1, drop = FALSE], na.rm = TRUE) / 2
  
  ro_g1_sum <- rowSums(ro[, is_g1, drop = FALSE], na.rm = TRUE)
  ao_g1_sum <- rowSums(ao[, is_g1, drop = FALSE], na.rm = TRUE)
  p_g1_reads <- ao_g1_sum / (ro_g1_sum + ao_g1_sum)
  
  n_nonfinite_p_g1 <- sum(!is.finite(p_g1))
  n_nonfinite_p_g1_reads <- sum(!is.finite(p_g1_reads))
  
  g1_poly_idx <- which(is.finite(p_g1) & p_g1 > 0 & p_g1 < 1)
  g1_reads_poly_idx <- which(is.finite(p_g1_reads) & p_g1_reads > 0 & p_g1_reads < 1)
  
  log_msg(sprintf("G1 AF computed from %d samples.", length(is_g1)), "INFO", SEC)
  log_msg(sprintf("  Finite G1 AF (genotypes): %d of %d SNPs",
                  sum(is.finite(p_g1)), n_snps), "INFO", SEC)
  log_msg(sprintf("  Finite G1 AF (reads):     %d of %d SNPs",
                  sum(is.finite(p_g1_reads)), n_snps), "INFO", SEC)
  log_msg(sprintf("  Polymorphic in G1 (genotypes): %d of %d SNPs",
                  length(g1_poly_idx), n_snps), "INFO", SEC)
  log_msg(sprintf("  Polymorphic in G1 (reads):     %d of %d SNPs",
                  length(g1_reads_poly_idx), n_snps), "INFO", SEC)
  
  if (n_nonfinite_p_g1 > 0) {
    log_msg(sprintf("%d SNPs have undefined G1 genotype-based AF (zero-weighted in GRM).",
                    n_nonfinite_p_g1), "WARN", SEC)
  }
  if (n_nonfinite_p_g1_reads > 0) {
    log_msg(sprintf("%d SNPs have undefined G1 read-based AF.", n_nonfinite_p_g1_reads),
            "WARN", SEC)
  }
  if (length(g1_poly_idx) > 0) {
    maf_vals <- pmin(p_g1[g1_poly_idx], 1 - p_g1[g1_poly_idx])
    log_msg(sprintf("  G1 MAF distribution: mean=%.4f  median=%.4f  Q25=%.4f  Q75=%.4f",
                    mean(maf_vals), median(maf_vals),
                    quantile(maf_vals, 0.25), quantile(maf_vals, 0.75)),
            "INFO", SEC)
    rm(maf_vals)
  }
  
  ## ---- Export VCF if external tools need it ----
  needs_vcf <- run_linkimputeR || run_gtImputation
  
  if (needs_vcf) {
    
    # Step 1: Export base VCF from GDS
    if (!file.exists(vcf_export)) {
      log_msg("Exporting GDS to VCF for external tools...", "INFO", SEC)
      gds <- seqOpen(gds_file, readonly = TRUE)
      
      available_fmt <- seqSummary(gds, "annotation/format", check = "none")
      required_fmt <- c("DP", "GQ", "RO", "AO")
      missing_fmt <- setdiff(required_fmt, available_fmt$id)
      if (length(missing_fmt) > 0) {
        seqClose(gds)
        log_msg(sprintf("Required FORMAT fields missing: %s",
                        paste(missing_fmt, collapse = ", ")), "ERROR", SEC)
        stop("Required FORMAT fields missing from GDS: ",
             paste(missing_fmt, collapse = ", "))
      }
      log_msg(sprintf("FORMAT fields verified: %s",
                      paste(required_fmt, collapse = ", ")), "INFO", SEC)
      
      seqGDS2VCF(gds, vcf_export, fmt.var = c("DP", "GQ", "RO", "AO"))
      seqClose(gds)
      log_file_info(vcf_export, "Base VCF export", SEC)
    } else {
      log_msg("VCF export already exists (skipping).", "INFO", SEC)
      log_file_info(vcf_export, "Base VCF export", SEC)
    }
    
    # Step 2: Apply DP-based GT masking to exported VCF
    if (!file.exists(vcf_export_masked)) {
      log_msg("Patching exported VCF to apply DP-based GT masking...", "INFO", SEC)
      t_mask <- Sys.time()
      
      all_lines <- readLines(gzfile(vcf_export))
      log_msg(sprintf("Read %d lines from exported VCF.", length(all_lines)),
              "INFO", SEC)
      
      is_meta   <- startsWith(all_lines, "##")
      is_colhdr <- startsWith(all_lines, "#CHROM")
      is_data   <- !is_meta & !is_colhdr
      
      meta_lines <- all_lines[is_meta]
      col_line   <- all_lines[is_colhdr]
      data_lines <- all_lines[is_data]
      n_data     <- length(data_lines)
      
      stopifnot(
        "Mismatch: VCF data lines vs DP mask rows" =
          n_data == nrow(mask_dp_any)
      )
      log_msg(sprintf("VCF data lines: %d (matches mask rows: %d)",
                      n_data, nrow(mask_dp_any)), "INFO", SEC)
      
      ## Split all lines in one vectorized call
      split_lines <- strsplit(data_lines, "\t", fixed = TRUE)
      n_fixed <- 9L
      samp_mat <- do.call(rbind, lapply(split_lines, function(x) x[-(1:n_fixed)]))
      fixed_strings <- vapply(split_lines, function(x) {
        paste(x[1:n_fixed], collapse = "\t")
      }, character(1))
      rm(split_lines)
      
      stopifnot(
        "Sample count mismatch between VCF and mask" =
          ncol(samp_mat) == ncol(mask_dp_any)
      )
      
      ## Replace GT field with "./." at masked positions
      mask_idx <- which(mask_dp_any, arr.ind = FALSE)
      n_gt_masked <- length(mask_idx)
      if (n_gt_masked > 0) {
        samp_mat[mask_idx] <- sub("^[^:]*", "./.", samp_mat[mask_idx])
      }
      log_msg(sprintf("GT fields masked in VCF: %d", n_gt_masked), "INFO", SEC)
      
      samp_strings <- do.call(paste, c(as.data.frame(samp_mat, stringsAsFactors = FALSE), sep = "\t"))
      new_data <- paste(fixed_strings, samp_strings, sep = "\t")
      rm(samp_mat, fixed_strings, samp_strings)
      
      out_lines <- c(meta_lines, col_line, new_data)
      con_out <- gzfile(vcf_export_masked, open = "wt")
      writeLines(out_lines, con_out)
      close(con_out)
      
      elapsed_mask <- round(difftime(Sys.time(), t_mask, units = "secs"), 1)
      log_file_info(vcf_export_masked, "DP-masked VCF", SEC)
      log_msg(sprintf("DP masking completed in %.1f seconds.", elapsed_mask),
              "TIMING", SEC)
      
      rm(all_lines, meta_lines, col_line, data_lines, new_data, out_lines)
      gc()
    } else {
      log_msg("DP-masked VCF already exists (skipping).", "INFO", SEC)
      log_file_info(vcf_export_masked, "DP-masked VCF", SEC)
    }
    
    # Step 3: Patch DP-masked VCF to add AD field for LinkImputeR
    if (run_linkimputeR && !file.exists(vcf_export_ad)) {
      log_msg("Patching DP-masked VCF to add AD field from RO+AO...", "INFO", SEC)
      t_patch <- Sys.time()
      
      all_lines <- readLines(gzfile(vcf_export_masked))
      log_msg(sprintf("Read %d lines from DP-masked VCF.", length(all_lines)),
              "INFO", SEC)
      
      is_meta   <- startsWith(all_lines, "##")
      is_colhdr <- startsWith(all_lines, "#CHROM")
      is_data   <- !is_meta & !is_colhdr
      
      meta_lines <- all_lines[is_meta]
      col_line   <- all_lines[is_colhdr]
      data_lines <- all_lines[is_data]
      n_data     <- length(data_lines)
      
      stopifnot(
        "Mismatch: VCF data lines vs ro/ao matrix rows" =
          n_data == nrow(ro)
      )
      
      log_msg(sprintf("Building AD strings for %d x %d cells...",
                      n_data, ncol(ro)), "INFO", SEC)
      
      ro_char <- matrix(as.character(as.integer(ro)), nrow = nrow(ro))
      ao_char <- matrix(as.character(as.integer(ao)), nrow = nrow(ao))
      ro_char[is.na(ro)] <- "."
      ao_char[is.na(ao)] <- "."
      
      ad_mat <- matrix(paste(ro_char, ao_char, sep = ","),
                       nrow = nrow(ro), ncol = ncol(ro))
      rm(ro_char, ao_char)
      
      split_lines <- strsplit(data_lines, "\t", fixed = TRUE)
      n_fixed <- 9L
      samp_mat <- do.call(rbind, lapply(split_lines, function(x) x[-(1:n_fixed)]))
      
      fixed_strings <- vapply(split_lines, function(x) {
        x[n_fixed] <- paste0(x[n_fixed], ":AD")
        paste(x[1:n_fixed], collapse = "\t")
      }, character(1))
      rm(split_lines)
      
      samp_ad_mat <- matrix(
        paste(samp_mat, ad_mat, sep = ":"),
        nrow = n_data, ncol = ncol(samp_mat)
      )
      rm(samp_mat, ad_mat)
      
      samp_strings <- do.call(paste, c(as.data.frame(samp_ad_mat, stringsAsFactors = FALSE), sep = "\t"))
      new_data <- paste(fixed_strings, samp_strings, sep = "\t")
      rm(samp_ad_mat, fixed_strings, samp_strings)
      
      ad_hdr <- '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the reference and alternate alleles in the order listed">'
      out_lines <- c(meta_lines, ad_hdr, col_line, new_data)
      
      con_out <- gzfile(vcf_export_ad, open = "wt")
      writeLines(out_lines, con_out)
      close(con_out)
      
      elapsed_patch <- round(difftime(Sys.time(), t_patch, units = "secs"), 1)
      log_file_info(vcf_export_ad, "AD-patched VCF", SEC)
      log_msg(sprintf("AD patching completed in %.1f seconds.", elapsed_patch),
              "TIMING", SEC)
      
      rm(all_lines, meta_lines, col_line, data_lines, new_data, out_lines)
      gc()
    } else if (run_linkimputeR) {
      log_msg("AD-patched VCF already exists (skipping).", "INFO", SEC)
      log_file_info(vcf_export_ad, "AD-patched VCF", SEC)
    }
    
  } else {
    log_msg("No external tools selected — skipping VCF export.", "INFO", SEC)
  }
})


# ===========================================================================
# UTILITY: PSD bending for GRMs (with logging)
# ===========================================================================

make_psd <- function(G, tol = 1e-6, name = "GRM", section = NULL) {
  eig <- eigen(G, symmetric = TRUE)
  n_neg <- sum(eig$values < tol)
  if (n_neg > 0) {
    most_neg <- min(eig$values)
    log_msg(sprintf("Bending %s: %d eigenvalues < %.1e (most negative: %.6e) -> set to %.1e",
                    name, n_neg, tol, most_neg, tol), "INFO", section)
    eig$values[eig$values < tol] <- tol
    G <- tcrossprod(eig$vectors * rep(sqrt(eig$values), each = nrow(eig$vectors)))
  } else {
    log_msg(sprintf("PSD check %s: already positive definite (min eigenvalue: %.6e).",
                    name, min(eig$values)), "INFO", section)
  }
  (G + t(G)) / 2
}


# ===========================================================================
# 2. METHOD 1: MEAN IMPUTATION (baseline)
# ===========================================================================

log_section("Mean Imputation", {
  
  SEC <- "MeanImp"
  
  if (run_mean_imputation) {
    
    geno_mean <- geno
    
    ## Fill missing calls from G1 base frequencies
    fill_vals <- ifelse(is.finite(p_g1), 2 * p_g1, 0)
    na_mask <- is.na(geno_mean)
    n_filled <- sum(na_mask)
    geno_mean[na_mask] <- fill_vals[row(geno_mean)][na_mask]
    rm(na_mask, fill_vals)
    
    log_msg(sprintf("Mean imputation: %d missing genotypes filled.", n_filled),
            "INFO", SEC)
    log_msg(sprintf("Remaining NAs: %d", sum(is.na(geno_mean))), "INFO", SEC)
    
    log_msg("Building GRM (VanRaden Method 1) centered on G1 frequencies...",
            "INFO", SEC)
    
    p_mean_use <- p_g1
    p_mean_use[!is.finite(p_mean_use)] <- 0
    
    M_mean <- sweep(t(geno_mean), 2, 2 * p_mean_use, "-")
    w_mean <- 2 * p_mean_use * (1 - p_mean_use)
    zero_var_mean <- !is.finite(w_mean) | w_mean <= 0
    
    if (any(zero_var_mean)) {
      log_msg(sprintf("Zero-weighting %d SNPs with zero/undefined G1 variance.",
                      sum(zero_var_mean)), "INFO", SEC)
      M_mean[, zero_var_mean] <- 0
    }
    
    denom_mean <- sum(w_mean[!zero_var_mean])
    log_msg(sprintf("GRM denominator: %.2f (from %d informative SNPs)",
                    denom_mean, sum(!zero_var_mean)), "INFO", SEC)
    
    if (denom_mean <= 0) {
      log_msg("FATAL: GRM denominator is zero.", "ERROR", SEC)
      stop("FATAL: GRM denominator is zero after zero-weighting non-informative loci.")
    }
    
    G_mean <- tcrossprod(M_mean) / denom_mean
    rownames(G_mean) <- colnames(G_mean) <- sample_ids
    
    log_grm_summary(G_mean, "Mean Imputation (pre-bending)", generation, SEC)
    G_mean <- make_psd(G_mean, name = "Mean_Imp", section = SEC)
    
    saveRDS(geno_mean, file.path(imp_dir, "geno_mean_imputed.rds"))
    saveRDS(G_mean,    file.path(imp_dir, "GRM_mean_imputed.rds"))
    log_file_info(file.path(imp_dir, "geno_mean_imputed.rds"), "Geno matrix", SEC)
    log_file_info(file.path(imp_dir, "GRM_mean_imputed.rds"), "GRM", SEC)
    
  } else {
    log_msg("Skipped (run_mean_imputation = FALSE)", "INFO", SEC)
  }
})


# ===========================================================================
# 3. METHOD 2: LinkImputeR
# ===========================================================================

log_section("LinkImputeR", {
  
  SEC <- "LinkImpR"
  
  if (run_linkimputeR) {
    
    if (!file.exists(linkimputeR_jar)) {
      log_msg(sprintf("LinkImputeR.jar not found at: %s", linkimputeR_jar),
              "WARN", SEC)
      log_msg("Build from source: https://github.com/danielmoney/LinkImputeR",
              "INFO", SEC)
      linkimputeR_ready <- FALSE
    } else {
      log_file_info(linkimputeR_jar, "LinkImputeR JAR", SEC)
      linkimputeR_ready <- TRUE
    }
    
    if (linkimputeR_ready) {
      
      linkimp_dir <- file.path(imp_dir, "linkimputeR")
      dir.create(linkimp_dir, recursive = TRUE, showWarnings = FALSE)
      
      ini_file <- file.path(linkimp_dir, "accuracy.ini")
      
      ini_content <- paste0(
        "[Input]\n",
        "filename = ", normalizePath(vcf_export_ad, winslash = "/", mustWork = FALSE), "\n",
        "save = ", file.path(linkimp_dir, "filtered_input.vcf.gz"), "\n",
        "maxdepth = 200\n",
        "\n",
        "[InputFilters]\n",
        "maf = 0.0\n",
        "positionmissing = 1.0\n",
        "\n",
        "[CaseFilters]\n",
        "maf = 0.0\n",
        "missing = 1.0\n",
        "\n",
        "[Global]\n",
        "depth = 2\n",
        "error = 0.01\n",
        "\n",
        "[Accuracy]\n",
        "maskmethod = bysnp\n",
        "accuracymethod = correct\n",
        "numbermasked = 10000\n",
        "mindepth = 5\n",
        "\n",
        "[Stats]\n",
        "root = ", file.path(linkimp_dir, "accuracy_stats_"), "\n",
        "level = pretty\n",
        "\n",
        "[Output]\n",
        "control = ", file.path(linkimp_dir, "impute_control.xml"), "\n",
        "\n",
        "[Log]\n",
        "file = ", file.path(linkimp_dir, "linkimputeR.log"), "\n",
        "level = brief\n"
      )
      
      writeLines(ini_content, ini_file)
      log_msg(sprintf("INI written to: %s", ini_file), "INFO", SEC)
      
      log_msg("Running LinkImputeR accuracy assessment...", "INFO", SEC)
      log_msg(sprintf("Estimated time: 30-60 minutes for ~%d SNPs x ~%d samples",
                      n_snps, n_samples), "INFO", SEC)
      
      ## Build command arguments for system2
      li_cp <- paste(linkimputeR_jar, file.path(dirname(linkimputeR_jar), "lib", "*"),
                     sep = cp_sep)
      li_acc_args <- c(
        sprintf("-Xmx%dG", java_heap_gb),
        "-cp", shQuote(li_cp),
        "Executable.LinkImputeR",
        "-s", shQuote(normalizePath(ini_file, winslash = "/"))
      )
      
      ret <- run_external(java_path, li_acc_args, "LinkImputeR-accuracy", SEC)
      
      sum_file <- file.path(linkimp_dir, "accuracy_stats_sum.dat")
      if (file.exists(sum_file)) {
        log_msg("--- LinkImputeR Accuracy Summary ---", "INFO", SEC)
        
        acc_table <- read.delim(sum_file, header = TRUE, stringsAsFactors = FALSE)
        ## Log each row of the accuracy table
        for (ri in seq_len(nrow(acc_table))) {
          log_msg(sprintf("  Case: %-20s Accuracy=%.4f  Samples=%d  Positions=%d",
                          acc_table$Name[ri], acc_table$Accuracy[ri],
                          acc_table$Samples[ri], acc_table$Positions[ri]),
                  "INFO", SEC)
        }
        
        valid_cases <- acc_table$Samples > 0 & acc_table$Positions > 0
        
        if (any(valid_cases)) {
          best_idx  <- which(valid_cases)[which.max(acc_table$Accuracy[valid_cases])]
          best_case <- acc_table$Name[best_idx]
          log_msg(sprintf("Best case: %s (accuracy=%.4f, samples=%d, positions=%d)",
                          best_case, acc_table$Accuracy[best_idx],
                          acc_table$Samples[best_idx], acc_table$Positions[best_idx]),
                  "INFO", SEC)
        } else {
          log_msg("All cases have 0 samples or positions. LinkImputeR could not process data.",
                  "WARN", SEC)
          best_case <- NULL
        }
      } else {
        log_msg(sprintf("Accuracy results not found at: %s", sum_file), "WARN", SEC)
        best_case <- NULL
      }
      
      if (!is.null(best_case)) {
        log_msg(sprintf("Running LinkImputeR final imputation with: %s", best_case),
                "INFO", SEC)
        
        linkimp_vcf <- file.path(linkimp_dir, "imputed_output.vcf.gz")
        control_xml <- file.path(linkimp_dir, "impute_control.xml")
        
        li_imp_args <- c(
          sprintf("-Xmx%dG", java_heap_gb),
          "-cp", shQuote(li_cp),
          "Executable.LinkImputeR",
          shQuote(normalizePath(control_xml, winslash = "/", mustWork = FALSE)),
          shQuote(best_case),
          shQuote(normalizePath(linkimp_vcf, winslash = "/", mustWork = FALSE))
        )
        
        ret <- run_external(java_path, li_imp_args, "LinkImputeR-impute", SEC)
        
        if (file.exists(linkimp_vcf)) {
          log_msg("Reading LinkImputeR imputed VCF...", "INFO", SEC)
          log_file_info(linkimp_vcf, "Imputed VCF", SEC)
          
          vcf_li <- read.vcfR(linkimp_vcf, verbose = FALSE)
          
          gt_li <- extract.gt(vcf_li, element = "GT")
          geno_li <- matrix(NA_integer_, nrow = nrow(gt_li), ncol = ncol(gt_li),
                            dimnames = dimnames(gt_li))
          geno_li[gt_li == "0/0" | gt_li == "0|0"] <- 0L
          geno_li[gt_li == "0/1" | gt_li == "1/0" | gt_li == "0|1" | gt_li == "1|0"] <- 1L
          geno_li[gt_li == "1/1" | gt_li == "1|1"] <- 2L
          
          log_msg(sprintf("Imputed matrix: %d SNPs x %d samples",
                          nrow(geno_li), ncol(geno_li)), "INFO", SEC)
          log_msg(sprintf("Remaining NAs after LinkImputeR: %d (%.2f%%)",
                          sum(is.na(geno_li)),
                          100 * mean(is.na(geno_li))), "INFO", SEC)
          
          li_snp_ids <- rownames(geno_li)
          p_li_g1 <- p_g1[match(li_snp_ids, snp_ids)]
          
          if (all(is.na(p_li_g1))) {
            log_msg("Direct SNP ID match failed. Trying positional matching...",
                    "WARN", SEC)
            li_fix <- vcf_li@fix[, c("CHROM", "POS")]
            li_snp_ids <- paste0(li_fix[, "CHROM"], "_", li_fix[, "POS"])
            p_li_g1 <- p_g1[match(li_snp_ids, snp_ids)]
            log_msg(sprintf("Positional match: %d of %d SNPs matched.",
                            sum(!is.na(p_li_g1)), length(p_li_g1)), "INFO", SEC)
          }
          
          ## Verify SNP order
          li_match_pos <- match(li_snp_ids, snp_ids)
          li_match_valid <- li_match_pos[!is.na(li_match_pos)]
          if (length(li_match_valid) > 1 && is.unsorted(li_match_valid)) {
            log_msg("LinkImputeR output SNPs NOT in same order as input GDS. GRM centering may be misaligned.",
                    "WARN", SEC)
          } else {
            log_msg("SNP order verified: matches input GDS.", "INFO", SEC)
          }
          rm(li_match_pos, li_match_valid)
          
          if (sum(is.na(geno_li)) > 0) {
            log_msg("Filling residual NAs with G1 mean imputation...", "INFO", SEC)
            fill_vals_li <- ifelse(is.finite(p_li_g1), 2 * p_li_g1, 0)
            na_mask_li <- is.na(geno_li)
            n_residual <- sum(na_mask_li)
            geno_li[na_mask_li] <- fill_vals_li[row(geno_li)][na_mask_li]
            log_msg(sprintf("Residual NAs filled: %d", n_residual), "INFO", SEC)
            rm(na_mask_li, fill_vals_li)
          }
          
          log_msg("Building GRM (VanRaden Method 1) from LinkImputeR SNPs...",
                  "INFO", SEC)
          
          p_li_use <- p_li_g1
          p_li_use[!is.finite(p_li_use)] <- 0
          
          M_li <- sweep(t(geno_li), 2, 2 * p_li_use, "-")
          w_li <- 2 * p_li_use * (1 - p_li_use)
          zero_var_li <- !is.finite(w_li) | w_li <= 0
          
          if (any(zero_var_li)) {
            log_msg(sprintf("Zero-weighting %d SNPs with zero/undefined G1 variance.",
                            sum(zero_var_li)), "INFO", SEC)
            M_li[, zero_var_li] <- 0
          }
          
          denom_li <- sum(w_li[!zero_var_li])
          log_msg(sprintf("GRM denominator: %.2f (from %d informative SNPs)",
                          denom_li, sum(!zero_var_li)), "INFO", SEC)
          
          if (denom_li <= 0) {
            log_msg("GRM denominator is zero for LinkImputeR output. Cannot build GRM.",
                    "WARN", SEC)
          } else {
            G_li <- tcrossprod(M_li) / denom_li
            rownames(G_li) <- colnames(G_li) <- colnames(geno_li)
            
            log_grm_summary(G_li, "LinkImputeR (pre-bending)", generation, SEC)
            G_li <- make_psd(G_li, name = "LinkImputeR", section = SEC)
            
            saveRDS(geno_li, file.path(imp_dir, "geno_linkimputeR.rds"))
            saveRDS(G_li,    file.path(imp_dir, "GRM_linkimputeR.rds"))
            log_file_info(file.path(imp_dir, "geno_linkimputeR.rds"), "Geno matrix", SEC)
            log_file_info(file.path(imp_dir, "GRM_linkimputeR.rds"), "GRM", SEC)
          }
          
          ip_li <- extract.gt(vcf_li, element = "IP")
          if (!is.null(ip_li)) {
            saveRDS(ip_li, file.path(imp_dir, "geno_linkimputeR_probs.rds"))
            log_msg("Imputation probabilities saved.", "INFO", SEC)
          }
          
          rm(vcf_li, gt_li)
        } else {
          log_msg("Imputed VCF not found. LinkImputeR may have failed.", "WARN", SEC)
        }
      }
    }
    
  } else {
    log_msg("Skipped (run_linkimputeR = FALSE)", "INFO", SEC)
  }
})


# ===========================================================================
# 4. METHOD 3: gtImputation (SOM-based neural network)
# ===========================================================================

log_section("gtImputation (SOM)", {
  
  SEC <- "gtImp"
  
  if (run_gtImputation) {
    
    gtimp_script <- file.path(gtimp_dir, "Package", "gtImputation.py")
    if (!file.exists(gtimp_script)) {
      gtimp_script <- file.path(gtimp_dir, "gtImputation.py")
    }
    
    if (!file.exists(gtimp_script)) {
      log_msg(sprintf("gtImputation.py not found at: %s", gtimp_script), "WARN", SEC)
      log_msg("Clone from https://github.com/GGFHF/gtImputation.git", "INFO", SEC)
      run_gtimp <- FALSE
    } else {
      log_file_info(gtimp_script, "gtImputation script", SEC)
      run_gtimp <- TRUE
    }
    
    if (run_gtimp) {
      
      gtimp_out_dir <- file.path(imp_dir, "gtImputation")
      dir.create(gtimp_out_dir, recursive = TRUE, showWarnings = FALSE)
      
      gtimp_input_vcf  <- normalizePath(vcf_export_masked, winslash = "/", mustWork = FALSE)
      gtimp_output_vcf <- file.path(gtimp_out_dir, "gtimp_imputed.vcf")
      
      yml_file <- file.path(gtimp_out_dir, "gtimp_config.yml")
      
      yml_content <- paste0(
        "# gtImputation configuration for Douglas-fir GS project\n",
        "# Generated: ", Sys.time(), "\n",
        "\n",
        "input:\n",
        "  vcf_file: '", gtimp_input_vcf, "'\n",
        "  format: VCF\n",
        "\n",
        "output:\n",
        "  vcf_file: '", normalizePath(gtimp_output_vcf, winslash = "/", mustWork = FALSE), "'\n",
        "\n",
        "som_parameters:\n",
        "  xdim: 3\n",
        "  ydim: 3\n",
        "  sigma: 1.0\n",
        "  learning_rate: 0.5\n",
        "  num_iteration: 1000\n",
        "\n",
        "imputation:\n",
        "  method: som\n",
        "  criterion: CK\n",
        "  min_r2: 0.1\n",
        "  max_snps: 15\n",
        "\n",
        "filtering:\n",
        "  maf: 0.0\n",
        "  max_missing: 1.0\n"
      )
      
      writeLines(yml_content, yml_file)
      log_msg(sprintf("Config written to: %s", yml_file), "INFO", SEC)
      
      n_missing_geno <- sum(is.na(geno))
      log_msg(sprintf("Dataset: %d SNPs x %d samples", n_snps, n_samples), "INFO", SEC)
      log_msg(sprintf("Missing genotypes to impute: %d (%.2f%%)",
                      n_missing_geno, 100 * n_missing_geno / length(geno)),
              "INFO", SEC)
      
      if (n_missing_geno > 100000) {
        est_hours_low  <- round(n_missing_geno * 0.5 / 3600, 1)
        est_hours_high <- round(n_missing_geno * 2.0 / 3600, 1)
        log_msg(sprintf("ESTIMATED RUNTIME: %.1f - %.1f hours (at 0.5-2.0 sec per SOM)",
                        est_hours_low, est_hours_high), "WARN", SEC)
        if (est_hours_high > 48) {
          log_msg("Runtime may exceed 48 hours. Consider compute node or alternative method.",
                  "WARN", SEC)
        }
      }
      
      ## Try config mode first, then direct CLI mode
      gtimp_args_config <- c(
        shQuote(normalizePath(gtimp_script, winslash = "/")),
        "--config", shQuote(normalizePath(yml_file, winslash = "/"))
      )
      
      gtimp_args_direct <- c(
        shQuote(normalizePath(gtimp_script, winslash = "/")),
        "--vcf", shQuote(gtimp_input_vcf),
        "--method", "som",
        "--xdim", "3", "--ydim", "3",
        "--sigma", "1.0",
        "--learning-rate", "0.5",
        "--num-iteration", "1000",
        "--criterion", "CK",
        "--min-r2", "0.1",
        "--max-snps", "15",
        "--output", shQuote(normalizePath(gtimp_output_vcf, winslash = "/", mustWork = FALSE))
      )
      
      log_msg("Attempting config mode...", "INFO", SEC)
      ret <- run_external(python_path, gtimp_args_config, "gtImputation-config", SEC)
      
      if (ret != 0) {
        log_msg("Config mode failed. Trying direct CLI mode...", "WARN", SEC)
        ret <- run_external(python_path, gtimp_args_direct, "gtImputation-direct", SEC)
      }
      
      if (file.exists(gtimp_output_vcf)) {
        log_msg("Reading gtImputation output VCF...", "INFO", SEC)
        log_file_info(gtimp_output_vcf, "Imputed VCF", SEC)
        
        vcf_gt <- read.vcfR(gtimp_output_vcf, verbose = FALSE)
        
        gt_gt <- extract.gt(vcf_gt, element = "GT")
        geno_gt <- matrix(NA_integer_, nrow = nrow(gt_gt), ncol = ncol(gt_gt),
                          dimnames = dimnames(gt_gt))
        geno_gt[gt_gt == "0/0" | gt_gt == "0|0"] <- 0L
        geno_gt[gt_gt == "0/1" | gt_gt == "1/0" | gt_gt == "0|1" | gt_gt == "1|0"] <- 1L
        geno_gt[gt_gt == "1/1" | gt_gt == "1|1"] <- 2L
        
        log_msg(sprintf("Imputed matrix: %d SNPs x %d samples",
                        nrow(geno_gt), ncol(geno_gt)), "INFO", SEC)
        log_msg(sprintf("Remaining NAs after gtImputation: %d (%.2f%%)",
                        sum(is.na(geno_gt)),
                        100 * mean(is.na(geno_gt))), "INFO", SEC)
        
        gt_snp_ids <- rownames(geno_gt)
        p_gt_g1 <- p_g1[match(gt_snp_ids, snp_ids)]
        
        if (all(is.na(p_gt_g1))) {
          log_msg("Direct SNP ID match failed. Trying positional matching...",
                  "WARN", SEC)
          gt_fix <- vcf_gt@fix[, c("CHROM", "POS")]
          gt_snp_ids <- paste0(gt_fix[, "CHROM"], "_", gt_fix[, "POS"])
          p_gt_g1 <- p_g1[match(gt_snp_ids, snp_ids)]
          log_msg(sprintf("Positional match: %d of %d SNPs matched.",
                          sum(!is.na(p_gt_g1)), length(p_gt_g1)), "INFO", SEC)
        }
        
        ## Verify SNP order
        gt_match_pos <- match(gt_snp_ids, snp_ids)
        gt_match_valid <- gt_match_pos[!is.na(gt_match_pos)]
        if (length(gt_match_valid) > 1 && is.unsorted(gt_match_valid)) {
          log_msg("gtImputation output SNPs NOT in same order as input GDS. GRM centering may be misaligned.",
                  "WARN", SEC)
        } else {
          log_msg("SNP order verified: matches input GDS.", "INFO", SEC)
        }
        rm(gt_match_pos, gt_match_valid)
        
        if (sum(is.na(geno_gt)) > 0) {
          log_msg("Filling residual NAs with G1 mean imputation...", "INFO", SEC)
          fill_vals_gt <- ifelse(is.finite(p_gt_g1), 2 * p_gt_g1, 0)
          na_mask_gt <- is.na(geno_gt)
          n_residual <- sum(na_mask_gt)
          geno_gt[na_mask_gt] <- fill_vals_gt[row(geno_gt)][na_mask_gt]
          log_msg(sprintf("Residual NAs filled: %d", n_residual), "INFO", SEC)
          rm(na_mask_gt, fill_vals_gt)
        }
        
        log_msg("Building GRM (VanRaden Method 1) from gtImputation SNPs...",
                "INFO", SEC)
        
        p_gt_use <- p_gt_g1
        p_gt_use[!is.finite(p_gt_use)] <- 0
        
        M_gt <- sweep(t(geno_gt), 2, 2 * p_gt_use, "-")
        w_gt <- 2 * p_gt_use * (1 - p_gt_use)
        zero_var_gt <- !is.finite(w_gt) | w_gt <= 0
        
        if (any(zero_var_gt)) {
          log_msg(sprintf("Zero-weighting %d SNPs with zero/undefined G1 variance.",
                          sum(zero_var_gt)), "INFO", SEC)
          M_gt[, zero_var_gt] <- 0
        }
        
        denom_gt <- sum(w_gt[!zero_var_gt])
        log_msg(sprintf("GRM denominator: %.2f (from %d informative SNPs)",
                        denom_gt, sum(!zero_var_gt)), "INFO", SEC)
        
        if (denom_gt <= 0) {
          log_msg("GRM denominator is zero for gtImputation output. Cannot build GRM.",
                  "WARN", SEC)
        } else {
          G_gt <- tcrossprod(M_gt) / denom_gt
          rownames(G_gt) <- colnames(G_gt) <- colnames(geno_gt)
          
          log_grm_summary(G_gt, "gtImputation (pre-bending)", generation, SEC)
          G_gt <- make_psd(G_gt, name = "gtImputation", section = SEC)
          
          saveRDS(geno_gt, file.path(imp_dir, "geno_gtImputation.rds"))
          saveRDS(G_gt,    file.path(imp_dir, "GRM_gtImputation.rds"))
          log_file_info(file.path(imp_dir, "geno_gtImputation.rds"), "Geno matrix", SEC)
          log_file_info(file.path(imp_dir, "GRM_gtImputation.rds"), "GRM", SEC)
        }
        
        rm(vcf_gt, gt_gt)
      } else {
        log_msg(sprintf("gtImputation output VCF not found at: %s", gtimp_output_vcf),
                "WARN", SEC)
        log_msg("Common issues: missing Python packages, VCF format, or memory exhaustion.",
                "INFO", SEC)
      }
    }
    
  } else {
    log_msg("Skipped (run_gtImputation = FALSE)", "INFO", SEC)
  }
})


# ===========================================================================
# 5. METHOD 4: KGD DOSAGE-BASED GRM (no imputation)
# ===========================================================================
#
# [PATCHED — complete rewrite of KGD integration]
#
# Allele coding note: our geno counts ALT alleles (0/1/2) and p_g1_reads
# is the ALT allele frequency. KGD internally counts REF alleles and uses
# REF frequency. However, VanRaden's GRM is invariant to allele coding:
# if X = genon - 2p then flipping allele coding gives -X, and
# tcrossprod(-X) = tcrossprod(X). The diagonal correction terms (P0*P1,
# Kdepth) are symmetric in p and 1-p.

log_section("KGD Depth-Adjusted GRM", {
  
  SEC <- "KGD"
  
  if (run_kgd) {
    
    kgd_out_dir <- file.path(imp_dir, "kgd")
    dir.create(kgd_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    ## ---- 5a. Export RA file ----
    log_msg("Writing vcf2ra.py-compatible RA file for reproducibility...", "INFO", SEC)
    ra_file <- file.path(kgd_out_dir, "allele_counts.tsv")
    
    ro_int <- matrix(as.integer(ro), nrow = nrow(ro), ncol = ncol(ro))
    ao_int <- matrix(as.integer(ao), nrow = nrow(ao), ncol = ncol(ao))
    ro_int[is.na(ro_int)] <- 0L
    ao_int[is.na(ao_int)] <- 0L
    
    ra_matrix <- matrix(paste(ro_int, ao_int, sep = ","),
                        nrow = n_snps, ncol = n_samples)
    rm(ro_int, ao_int)
    
    ra_row_strings <- apply(ra_matrix, 1, paste, collapse = "\t")
    rm(ra_matrix)
    data_lines_ra <- paste(snp_chrom, snp_pos, ra_row_strings, sep = "\t")
    rm(ra_row_strings)
    
    writeLines(c(paste(c("CHROM", "POS", sample_ids), collapse = "\t"),
                 data_lines_ra), ra_file)
    rm(data_lines_ra)
    log_file_info(ra_file, "RA file", SEC)
    log_msg(sprintf("RA file: %d SNPs x %d samples", n_snps, n_samples), "INFO", SEC)
    
    ## ---- 5b. Source KGD in an isolated environment ----
    kgd_source <- file.path(kgd_dir, "GBS-Chip-Gmatrix.R")
    
    if (!file.exists(kgd_source)) {
      log_msg(sprintf("FATAL: KGD source not found at: %s", kgd_source), "ERROR", SEC)
      stop("FATAL: KGD source not found at:\n  ", kgd_source)
    }
    log_file_info(kgd_source, "KGD source", SEC)
    
    kgd_env <- new.env(parent = globalenv())
    
    kgd_env$functions.only  <- TRUE
    kgd_env$gform           <- "vcf"
    kgd_env$outlevel        <- 9
    kgd_env$cex.pointsize   <- 1
    kgd_env$hirel.thresh    <- 0.9
    kgd_env$alleles.keep    <- TRUE
    kgd_env$use.Rcpp        <- TRUE
    kgd_env$nThreads        <- max(1, ncores)
    
    log_msg("Sourcing KGD (functions only, isolated env)...", "INFO", SEC)
    source(kgd_source, local = kgd_env)
    
    ## ---- 5c. Populate KGD globals ----
    log_msg("Populating KGD globals from in-memory matrices...", "INFO", SEC)
    
    kgd_env$nsnps     <- n_snps
    kgd_env$nind      <- n_samples
    kgd_env$seqID     <- sample_ids
    kgd_env$SNP_Names <- snp_ids
    
    kgd_env$genon <- t(geno)
    storage.mode(kgd_env$genon) <- "integer"
    
    ro_t <- t(ro)
    ao_t <- t(ao)
    ro_t[is.na(ro_t)] <- 0
    ao_t[is.na(ao_t)] <- 0
    
    kgd_env$depth <- ro_t + ao_t
    kgd_env$depth[is.na(kgd_env$genon)] <- 0
    storage.mode(kgd_env$depth) <- "integer"
    
    kgd_env$alleles <- matrix(0L, nrow = n_samples, ncol = 2 * n_snps)
    kgd_env$alleles[, seq(1, 2 * n_snps - 1, 2)] <- as.integer(ro_t)
    kgd_env$alleles[, seq(2, 2 * n_snps, 2)]     <- as.integer(ao_t)
    rm(ro_t, ao_t)
    
    kgd_env$p <- p_g1_reads
    kgd_env$p[!is.finite(kgd_env$p)] <- 0
    
    kgd_env$fcolo <- rep("black", n_samples)
    
    kgd_env$redosamples()
    
    log_msg(sprintf("KGD globals: nsnps=%d  nind=%d", kgd_env$nsnps, kgd_env$nind),
            "INFO", SEC)
    log_msg(sprintf("  genon: %dx%d  depth: %dx%d  alleles: %dx%d",
                    nrow(kgd_env$genon), ncol(kgd_env$genon),
                    nrow(kgd_env$depth), ncol(kgd_env$depth),
                    nrow(kgd_env$alleles), ncol(kgd_env$alleles)),
            "INFO", SEC)
    
    mean_depth_nz <- mean(kgd_env$depth[kgd_env$depth > 0])
    call_rate_kgd <- 1 - mean(kgd_env$depth == 0)
    log_msg(sprintf("  Mean depth (non-zero): %.2f  Call rate: %.4f",
                    mean_depth_nz, call_rate_kgd), "INFO", SEC)
    
    ## ---- 5d. Run calcG with full diagnostics ----
    old_wd <- getwd()
    setwd(kgd_out_dir)
    
    puse_kgd <- p_g1_reads
    if (any(!is.finite(puse_kgd))) {
      n_nonfinite <- sum(!is.finite(puse_kgd))
      log_msg(sprintf("%d SNPs have undefined G1 read-based AF; replacing with 0 in puse.",
                      n_nonfinite), "WARN", SEC)
      puse_kgd[!is.finite(puse_kgd)] <- 0
    }
    
    log_msg("Running KGD calcG() (calclevel=9, npc=4, G1-centered)...", "INFO", SEC)
    
    Gresult <- tryCatch({
      kgd_env$calcG(
        snpsubset  = seq_len(n_snps),
        sfx        = "postQC_all",
        puse       = puse_kgd,
        indsubset  = seq_len(n_samples),
        npc        = 4,
        calclevel  = 9
      )
    }, error = function(e) {
      setwd(old_wd)
      log_msg(sprintf("FATAL: KGD calcG() failed: %s", conditionMessage(e)),
              "ERROR", SEC)
      stop(e)
    })
    
    setwd(old_wd)
    
    ## ---- 5e. Extract G5, apply PSD bending, save ----
    G_kgd <- Gresult$G5
    rownames(G_kgd) <- colnames(G_kgd) <- sample_ids
    
    log_grm_summary(G_kgd, "KGD G5 (pre-bending)", generation, SEC)
    G_kgd <- make_psd(G_kgd, name = "KGD_G5", section = SEC)
    
    saveRDS(G_kgd,    file.path(imp_dir, "GRM_kgd_depth_adjusted.rds"))
    saveRDS(Gresult,  file.path(imp_dir, "KGD_full_result.rds"))
    log_file_info(file.path(imp_dir, "GRM_kgd_depth_adjusted.rds"), "GRM", SEC)
    log_file_info(file.path(imp_dir, "KGD_full_result.rds"), "Full KGD result", SEC)
    
    ## Log KGD diagnostic plot files
    kgd_plots <- list.files(kgd_out_dir, pattern = "\\.png$", full.names = TRUE)
    if (length(kgd_plots) > 0) {
      log_msg(sprintf("KGD diagnostic plots (%d files) written to: %s",
                      length(kgd_plots), kgd_out_dir), "INFO", SEC)
      for (pf in kgd_plots) log_msg(sprintf("  %s", basename(pf)), "DEBUG", SEC)
    }
    
    ## ---- 5f. Clean up KGD environment ----
    rm(kgd_env)
    log_msg("KGD environment removed (all KGD globals cleaned up).", "INFO", SEC)
    
  } else {
    log_msg("Skipped (run_kgd = FALSE)", "INFO", SEC)
  }
})


# ===========================================================================
# 6. DIAGNOSTIC COMPARISON
# ===========================================================================

log_section("Diagnostic Comparison", {
  
  SEC <- "Diag"
  
  library(ggplot2)
  
  ## ---- 6a. Compare GRM diagonals ----
  log_msg("Comparing GRM diagonal elements (self-relatedness)...", "INFO", SEC)
  
  diag_df <- data.frame(Sample = colnames(geno))
  
  add_diag_safely <- function(df, G, col_name, expected_n) {
    if (nrow(G) == expected_n) {
      df[[col_name]] <- diag(G)
    } else {
      log_msg(sprintf("%s GRM has %d samples vs expected %d — excluding.",
                      col_name, nrow(G), expected_n), "WARN", SEC)
    }
    df
  }
  
  if (exists("G_mean"))  diag_df <- add_diag_safely(diag_df, G_mean, "Mean_Imp", n_samples)
  if (exists("G_li"))    diag_df <- add_diag_safely(diag_df, G_li, "LinkImputeR", n_samples)
  if (exists("G_gt"))    diag_df <- add_diag_safely(diag_df, G_gt, "gtImputation", n_samples)
  if (exists("G_kgd"))   diag_df <- add_diag_safely(diag_df, G_kgd, "KGD", n_samples)
  
  n_methods <- ncol(diag_df) - 1
  if (n_methods >= 1) {
    diag_long <- tidyr::pivot_longer(diag_df, -Sample,
                                     names_to = "Method", values_to = "Gii")
    
    p1 <- ggplot(diag_long, aes(x = Gii, fill = Method)) +
      geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
      labs(title = "GRM Diagonal Elements by Imputation Method",
           x = expression(G[ii]~"(self-relatedness)"), y = "Count") +
      theme_minimal(base_size = 13)
    
    plot_file <- file.path(imp_dir, "diagnostic_Gii_comparison.png")
    ggsave(plot_file, p1, width = 9, height = 5, dpi = 200)
    log_file_info(plot_file, "Gii comparison plot", SEC)
  } else {
    log_msg("No GRMs computed — skipping diagonal comparison.", "INFO", SEC)
  }
  
  ## ---- 6b. Pairwise off-diagonal comparisons ----
  log_msg("Pairwise off-diagonal GRM comparisons...", "INFO", SEC)
  
  grm_list <- list()
  if (exists("G_mean"))  grm_list[["Mean_Imp"]]     <- G_mean
  if (exists("G_li"))    grm_list[["LinkImputeR"]]  <- G_li
  if (exists("G_gt"))    grm_list[["gtImputation"]] <- G_gt
  if (exists("G_kgd"))   grm_list[["KGD"]]          <- G_kgd
  
  if (length(grm_list) >= 2) {
    grm_dims <- sapply(grm_list, nrow)
    ref_dim  <- grm_dims[1]
    
    if (!all(grm_dims == ref_dim)) {
      log_msg("GRM dimensions differ across methods:", "WARN", SEC)
      for (nm in names(grm_dims)) {
        log_msg(sprintf("  %s: %dx%d", nm, grm_dims[nm], grm_dims[nm]),
                "INFO", SEC)
      }
      dim_counts <- table(grm_dims)
      target_dim <- as.integer(names(which.max(dim_counts)))
      grm_list <- grm_list[grm_dims == target_dim]
      log_msg(sprintf("Using %d GRMs with dimension %d.",
                      length(grm_list), target_dim), "INFO", SEC)
    }
  }
  
  if (length(grm_list) >= 2) {
    off_list <- lapply(grm_list, function(G) G[lower.tri(G)])
    set.seed(42)
    n_pairs  <- length(off_list[[1]])
    samp_idx <- sample(n_pairs, min(50000, n_pairs))
    
    method_names <- names(grm_list)
    plot_counter <- 0
    for (i in seq_along(method_names)) {
      for (j in seq_along(method_names)) {
        if (j <= i) next
        nm_x <- method_names[i]
        nm_y <- method_names[j]
        plot_counter <- plot_counter + 1
        
        cor_val <- cor(off_list[[nm_x]], off_list[[nm_y]], use = "complete.obs")
        log_msg(sprintf("Pairwise Gij correlation: %s vs %s  r=%.4f",
                        nm_x, nm_y, cor_val), "INFO", SEC)
        
        plot_df <- data.frame(
          x = off_list[[nm_x]][samp_idx],
          y = off_list[[nm_y]][samp_idx]
        )
        
        shared_range <- range(c(plot_df$x, plot_df$y), na.rm = TRUE)
        pp <- ggplot(plot_df, aes(x = x, y = y)) +
          geom_point(alpha = 0.05, size = 0.5) +
          geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
          labs(title = paste("Pairwise Relatedness:", nm_x, "vs", nm_y),
               x = bquote(G[ij] ~ "(" * .(nm_x) * ")"),
               y = bquote(G[ij] ~ "(" * .(nm_y) * ")")) +
          theme_minimal(base_size = 13) +
          coord_fixed(ratio = 1, xlim = shared_range, ylim = shared_range)
        
        fname <- paste0("diagnostic_Gij_", gsub(" ", "_", tolower(nm_x)),
                        "_vs_", gsub(" ", "_", tolower(nm_y)), ".png")
        ggsave(file.path(imp_dir, fname), pp, width = 7, height = 7, dpi = 200)
        log_file_info(file.path(imp_dir, fname), "Gij plot", SEC)
      }
    }
    log_msg(sprintf("Generated %d pairwise GRM comparison plots.", plot_counter),
            "INFO", SEC)
  } else {
    log_msg("Fewer than 2 GRMs computed — skipping pairwise comparisons.", "INFO", SEC)
  }
  
  ## ---- 6c. Summary statistics ----
  log_msg("", "INFO", SEC)
  log_msg("=== GRM Summary Statistics ===", "INFO", SEC)
  
  for (nm in names(grm_list)) {
    log_grm_summary(grm_list[[nm]], nm, generation, SEC)
  }
  
  ## ---- 6d. ssGBLUP compatibility note ----
  log_msg("", "INFO", SEC)
  log_msg("ssGBLUP note: For H-matrix compatibility, verify that", "INFO", SEC)
  log_msg("  mean(diag(G)[G1]) ~ mean(diag(A)[G1]) before blending.", "INFO", SEC)
})


# ===========================================================================
# 7. PREPARE FOR DOWNSTREAM BGLR
# ===========================================================================

log_section("Prepare BGLR Outputs", {
  
  SEC <- "BGLR"
  
  saveRDS(generation, file.path(imp_dir, "generation_assignment.rds"))
  saveRDS(p_g1,       file.path(imp_dir, "allele_freq_g1.rds"))
  saveRDS(p_g1_reads, file.path(imp_dir, "allele_freq_g1_reads.rds"))
  
  log_msg("Saved common output files:", "INFO", SEC)
  log_file_info(file.path(imp_dir, "generation_assignment.rds"),
                "Generation assignment", SEC)
  log_file_info(file.path(imp_dir, "allele_freq_g1.rds"),
                "G1 AF (genotypes)", SEC)
  log_file_info(file.path(imp_dir, "allele_freq_g1_reads.rds"),
                "G1 AF (reads)", SEC)
  
  if (exists("G_mean")) {
    log_file_info(file.path(imp_dir, "geno_mean_imputed.rds"), "Mean imp geno", SEC)
    log_file_info(file.path(imp_dir, "GRM_mean_imputed.rds"), "Mean imp GRM", SEC)
  }
  if (exists("G_li")) {
    log_file_info(file.path(imp_dir, "geno_linkimputeR.rds"), "LinkImputeR geno", SEC)
    log_file_info(file.path(imp_dir, "GRM_linkimputeR.rds"), "LinkImputeR GRM", SEC)
  }
  if (exists("G_gt")) {
    log_file_info(file.path(imp_dir, "geno_gtImputation.rds"), "gtImputation geno", SEC)
    log_file_info(file.path(imp_dir, "GRM_gtImputation.rds"), "gtImputation GRM", SEC)
  }
  if (exists("G_kgd")) {
    log_file_info(file.path(imp_dir, "GRM_kgd_depth_adjusted.rds"), "KGD GRM", SEC)
  }
  
  log_msg("All GRMs centered on G1 (base population) allele frequencies.", "INFO", SEC)
  log_msg("", "INFO", SEC)
  log_msg("Usage in BGLR:", "INFO", SEC)
  log_msg('  G <- readRDS("GRM_kgd_depth_adjusted.rds")', "INFO", SEC)
  log_msg('  fit <- BGLR(y = y, ETA = list(list(K = G, model = "RKHS")), ...)',
          "INFO", SEC)
  log_msg("", "INFO", SEC)
  log_msg("For ssGBLUP (H-matrix), blend G with A:", "INFO", SEC)
  log_msg("  Verify: mean(diag(G)[G1]) ~ mean(diag(A)[G1]) before blending.",
          "INFO", SEC)
})


# ===========================================================================
# 8. END-OF-RUN SUMMARY
# ===========================================================================

log_run_summary()