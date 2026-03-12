# =============================================================================
# build_pedigree_metadata.R
#
# Constructs the complete 3-generation Douglas-fir pedigree from All_sites.csv
# and founders_climateNA.csv. Outputs:
#
#   (1) pedigree_full.csv      — all individuals in pedigreemm/optiSel format
#                                columns: id, sire, dam, gen, site, genotyped
#   (2) pedigree_genotyped.csv — genotyped subset only (G1 + G2), same format
#   (3) pedigree_summary.txt   — counts, family structure, connectivity report
#
# ── Pedigree structure (verified 2026-03-04) ─────────────────────────────────
#
#  Generation  |  n    | Sites              | Genotyped | Notes
#  ------------|-------|--------------------|-----------|--------------------
#  G0 founders |  45   | Various BC provenances | No    | Ungenotyped seed
#               |       |                    |           | orchard clones
#  G1 parents  | 1,377 | lost / fleet / adam | Yes      | 56 selected as
#               |       |                    |           | G2 parents
#  G2 progeny  |   490 | bigtree / jordan   | Yes       | 3-generation tip
#
# ── External (ungenotyped) G1 parents ────────────────────────────────────────
#  Select ID | Site       | Cross (G0 fem × G0 male) | Connectivity
#  ----------|------------|--------------------------|------------------------
#  3216      | 32sechelt  | 300  × 578               | Connected via G0-300
#  3151      | 11tamahi   | 83   × 423               | Pedigree island
#  3132      | 03caycuse  | 241  × 453               | Pedigree island
#  3472      | 16muir     | 82   × 434               | Novel G0s (82, 434)
#  3522      | 28fleet    | 289  × 46                | Connected (common G0s)
#  3526      | 32sechelt  | 81   × 623               | Connected
#  3529      | 32sechelt  | 66   × 53                | Connected
#  3534      | 32sechelt  | 630  × 441               | Connected
#  3537      | 32sechelt  | 560  × 623               | Connected
#  3548      | 32sechelt  | 214  × 129               | Connected
#
# ── Reciprocal crosses (6 pairs, valid by design) ────────────────────────────
#  100×210, 344×441, 285×566, 53×231, 214×227, 129×33
#  These produce half-sib structure; treat as distinct families in models.
#
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(tibble)

# ── 0. Paths (edit as needed) ─────────────────────────────────────────────────
ALL_SITES_FILE  <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/All.sites.csv"
FOUNDERS_FILE   <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/founders_climateNA.csv"
OUT_FULL        <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/pedigree_full.csv"
OUT_GENOTYPED   <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/pedigree_genotyped.csv"
OUT_SUMMARY     <- "D:/OneDrive - NRCan RNCan/gs/doug-fir/data/pedigree_summary.txt"

# ── 1. Load inputs ─────────────────────────────────────────────────────────────
all_sites <- read_csv(ALL_SITES_FILE, col_types = cols(.default = col_character()))
founders  <- read_csv(FOUNDERS_FILE,  col_types = cols(.default = col_character()))

# Standardise: strip whitespace, normalise NA strings
all_sites <- all_sites %>%
  mutate(across(everything(), str_trim)) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("NA", ""), NA, .)))

cat("All_sites rows loaded:", nrow(all_sites), "\n")

# ── 2. Classify rows by site ──────────────────────────────────────────────────
G1_SITES  <- c("lost", "fleet", "adam")
G2_SITES  <- c("bigtree", "jordan")
EXT_SITES <- c("32sechelt", "11tamahi", "03caycuse", "16muir", "28fleet")

g1_rows  <- all_sites %>% filter(site %in% G1_SITES)
g2_rows  <- all_sites %>% filter(site %in% G2_SITES)
ext_rows <- all_sites %>% filter(site %in% EXT_SITES)

cat("G1 rows:", nrow(g1_rows),
    " | G2 rows:", nrow(g2_rows),
    " | External rows:", nrow(ext_rows), "\n")

# ── 3. Build G0 records ───────────────────────────────────────────────────────
# G0 individuals appear as fem/male parents in G1 rows and in external rows.
# They are identified by numeric IDs in the range 33–630.
# They are NOT genotyped — founders only.

g0_from_g1  <- g1_rows  %>% select(fem, male) %>% unlist() %>% unique()
g0_from_ext <- ext_rows %>% select(fem, male) %>% unlist() %>% unique()
g0_ids <- unique(c(g0_from_g1, g0_from_ext))
g0_ids <- g0_ids[!is.na(g0_ids)]

g0_pedigree <- tibble(
  id        = g0_ids,
  sire      = NA_character_,   # founders have no recorded parents
  dam       = NA_character_,
  gen       = "G0",
  site      = "founder",
  genotyped = FALSE,
  select_id = NA_character_,
  genotype_id = NA_character_
)

cat("G0 founders:", nrow(g0_pedigree), "\n")

# ── 4. Build G1 records ───────────────────────────────────────────────────────
# Each G1 row has a genotype.id (the lab sample ID) and optionally a select
# value (the numeric select ID used as parent in G2 crosses).
# sire/dam convention: use fem=dam, male=sire (consistent with pedigreemm).
# Note: 6 reciprocal cross pairs exist — kept as separate families.

g1_pedigree <- g1_rows %>%
  filter(!is.na(genotype.id)) %>%
  transmute(
    id          = genotype.id,
    sire        = male,           # G0 male parent
    dam         = fem,            # G0 female parent
    gen         = "G1",
    site        = site,
    genotyped   = TRUE,
    select_id   = select,         # NA if not selected as G2 parent
    genotype_id = genotype.id
  )

cat("G1 genotyped individuals:", nrow(g1_pedigree), "\n")
cat("  G1 selected as G2 parents:",
    sum(!is.na(g1_pedigree$select_id)), "\n")

# ── 5. Build external (ungenotyped) G1 parent records ────────────────────────
# These individuals were selected from external sites and used as parents in G2,
# but were NOT genotyped and have no genotype.id.
# Their pedigree parents (G0) are recorded in fem/male columns.
# ID = select value (3132, 3151, 3216, 3472, 3522, 3526, 3529, 3534, 3537, 3548)

ext_pedigree <- ext_rows %>%
  filter(!is.na(select)) %>%
  transmute(
    id          = select,
    sire        = male,
    dam         = fem,
    gen         = "G1_ext",      # G1 generation, external site, ungenotyped
    site        = site,
    genotyped   = FALSE,
    select_id   = select,
    genotype_id = NA_character_
  ) %>%
  distinct(id, .keep_all = TRUE)

cat("External ungenotyped G1 parents:", nrow(ext_pedigree), "\n")
cat("  IDs:", paste(sort(ext_pedigree$id), collapse = ", "), "\n")

# ── 6. Build G2 records ───────────────────────────────────────────────────────
# G2 parents are identified by their select IDs (3xxx) in fem/male columns.
# All G2 parents should be present in g1_pedigree (select_id column) or
# ext_pedigree (id column).

# Build a lookup: select_id → G1 genotype.id (for genotyped G1 parents)
select_to_gid <- g1_pedigree %>%
  filter(!is.na(select_id)) %>%
  select(select_id, genotype_id) %>%
  deframe()   # named vector: select_id → genotype_id

g2_pedigree <- g2_rows %>%
  filter(!is.na(genotype.id)) %>%
  transmute(
    id          = genotype.id,
    # Use genotype.id as parent ID where available; else keep select_id
    # (for external ungenotyped parents 3132, 3151, 3216, 3472, etc.)
    sire        = ifelse(male %in% names(select_to_gid),
                         select_to_gid[male], male),
    dam         = ifelse(fem  %in% names(select_to_gid),
                         select_to_gid[fem],  fem),
    # Also store the original select IDs for cross-referencing
    sire_select = male,
    dam_select  = fem,
    gen         = "G2",
    site        = site,
    genotyped   = TRUE,
    select_id   = NA_character_,
    genotype_id = genotype.id
  )

cat("G2 genotyped individuals:", nrow(g2_pedigree), "\n")

# Verify all G2 parent select IDs are resolved
g2_parent_selects <- unique(c(g2_rows$fem, g2_rows$male))
g2_parent_selects <- g2_parent_selects[!is.na(g2_parent_selects)]
all_g1_ids <- c(names(select_to_gid), ext_pedigree$id)
unresolved <- setdiff(g2_parent_selects, all_g1_ids)
if (length(unresolved) > 0) {
  warning("Unresolved G2 parent select IDs: ", paste(unresolved, collapse = ", "))
} else {
  cat("  All G2 parent select IDs resolved ✓\n")
}

# ── 7. Assemble full pedigree ─────────────────────────────────────────────────
# Standard columns for pedigreemm / optiSel / AGHmatrix:
#   id   = unique individual identifier
#   sire = father ID (NA for founders)
#   dam  = mother ID (NA for founders)
# Generation order must be: parents before offspring (founders first)

ped_full <- bind_rows(
  g0_pedigree  %>% select(id, sire, dam, gen, site, genotyped, select_id, genotype_id),
  g1_pedigree,
  ext_pedigree,
  g2_pedigree  %>% select(id, sire, dam, gen, site, genotyped, select_id, genotype_id)
)

# Verify all sire/dam values exist as id in the pedigree (except G0 founders)
all_ids <- ped_full$id
missing_sires <- ped_full %>%
  filter(!is.na(sire), !sire %in% all_ids) %>%
  select(id, sire, gen, site)
missing_dams  <- ped_full %>%
  filter(!is.na(dam),  !dam  %in% all_ids) %>%
  select(id, dam, gen, site)

if (nrow(missing_sires) > 0) {
  warning(nrow(missing_sires), " sire IDs not found in pedigree:")
  print(missing_sires)
}
if (nrow(missing_dams) > 0) {
  warning(nrow(missing_dams), " dam IDs not found in pedigree:")
  print(missing_dams)
}

cat("\nFull pedigree assembled:", nrow(ped_full), "individuals\n")

# ── 8. Genotyped-only subset ──────────────────────────────────────────────────
# For use directly in rrBLUP / sommer / ssGBLUP without G0 phantom rows.
# Ungenotyped parents (G0s, external G1s) are retained as NA parents here
# so the A-matrix builder can still recognise them as founders.

ped_genotyped <- ped_full %>%
  filter(genotyped == TRUE) %>%
  mutate(
    # Replace any sire/dam that is a G0 or external G1 (ungenotyped) with
    # a standardised founder ID so pedigreemm doesn't drop them.
    # Convention: keep as-is — G0 IDs like "100" are just founder labels.
    sire_final = sire,
    dam_final  = dam
  )

cat("Genotyped subset:", nrow(ped_genotyped),
    "(G1:", sum(ped_genotyped$gen == "G1"),
    ", G2:", sum(ped_genotyped$gen == "G2"), ")\n")

# ── 9. Write outputs ──────────────────────────────────────────────────────────
write_csv(ped_full,      OUT_FULL)
write_csv(ped_genotyped, OUT_GENOTYPED)

# ── 10. Summary report ────────────────────────────────────────────────────────
sink(OUT_SUMMARY)
cat("PEDIGREE METADATA SUMMARY\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("══ COUNTS ══════════════════════════════════════════════\n")
cat(sprintf("  G0 founders (ungenotyped)         : %4d\n", sum(ped_full$gen == "G0")))
cat(sprintf("  G1 genotyped (lost/fleet/adam)    : %4d\n", sum(ped_full$gen == "G1")))
cat(sprintf("  G1 external ungenotyped           : %4d\n", sum(ped_full$gen == "G1_ext")))
cat(sprintf("  G2 genotyped (bigtree/jordan)     : %4d\n", sum(ped_full$gen == "G2")))
cat(sprintf("  Total                             : %4d\n", nrow(ped_full)))

cat("\n══ G1 SITE BREAKDOWN ══════════════════════════════════\n")
g1_pedigree %>% count(site) %>%
  { for (i in seq_len(nrow(.))) cat(sprintf("  %-10s : %d\n", .[[i,1]], .[[i,2]])) }

cat("\n══ G2 SITE BREAKDOWN ══════════════════════════════════\n")
g2_pedigree %>% count(site) %>%
  { for (i in seq_len(nrow(.))) cat(sprintf("  %-10s : %d\n", .[[i,1]], .[[i,2]])) }

cat("\n══ G0 FAMILY STRUCTURE (G1 crosses) ═══════════════════\n")
g1_pedigree %>%
  mutate(cross = paste0(dam, "×", sire)) %>%
  count(cross, sort = TRUE) %>%
  head(20) %>%
  { for (i in seq_len(nrow(.))) cat(sprintf("  %-15s : %3d trees\n", .[[i,1]], .[[i,2]])) }

cat("\n══ EXTERNAL G1 PARENTS ═════════════════════════════════\n")
ext_pedigree %>%
  { for (i in seq_len(nrow(.)))
      cat(sprintf("  %s  site=%-12s  dam=%s  sire=%s\n",
                  .[i,"id"], .[i,"site"], .[i,"dam"], .[i,"sire"])) }

cat("\n══ G1 SELECTED AS G2 PARENTS ═══════════════════════════\n")
g1_pedigree %>%
  filter(!is.na(select_id)) %>%
  select(select_id, genotype_id, site) %>%
  arrange(select_id) %>%
  { for (i in seq_len(nrow(.)))
      cat(sprintf("  select=%s  gid=%-30s  site=%s\n",
                  .[i,"select_id"], .[i,"genotype_id"], .[i,"site"])) }

cat("\n══ G2 FAMILY STRUCTURE ═════════════════════════════════\n")
g2_rows %>%
  mutate(cross = paste0(fem, "×", male)) %>%
  count(cross, sort = TRUE) %>%
  { for (i in seq_len(nrow(.))) cat(sprintf("  %-15s : %3d trees\n", .[[i,1]], .[[i,2]])) }

cat("\n══ PEDIGREE CONNECTIVITY NOTES ═════════════════════════\n")
cat("  3216 (sechelt): G0 dam=300 bridges to main pedigree.\n")
cat("    300 produced G1 selects 3506, 3519, 3527 (all genotyped).\n")
cat("    3216's G2 progeny are genomic half-avuncular to those families.\n")
cat("  3132 (caycuse) and 3151 (tamahi): no shared G0 with main pedigree.\n")
cat("    8 G2 progeny at jordan form an isolated pedigree island.\n")
cat("  3472 (16muir): G0 parents 82 and 434 are unique to this branch.\n")
cat("    3472's G2 progeny connect to main pedigree only through\n")
cat("    non-3472 G1 co-parents (3518, 3552, 3495).\n")
cat("  6 reciprocal cross pairs in G1 (e.g. 100×210 and 210×100).\n")
cat("    Treat as separate full-sib families; note for A-matrix build.\n")

sink()
cat("Summary written to:", OUT_SUMMARY, "\n")
cat("\nDone. Output files:\n")
cat("  ", OUT_FULL, " (", nrow(ped_full), " rows)\n", sep = "")
cat("  ", OUT_GENOTYPED, " (", nrow(ped_genotyped), " rows)\n", sep = "")
cat("  ", OUT_SUMMARY, "\n", sep = "")
