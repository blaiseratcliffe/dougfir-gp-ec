# Literature Review Research Prompt — Douglas-fir Genomic Selection Manuscript

## Purpose

This document is a structured research prompt designed to compile a comprehensive literature review across five interconnected topics relevant to a Douglas-fir (*Pseudotsuga menziesii*) genomic selection study. The study uses exome capture-derived SNP data, single-step GBLUP (ssGBLUP), and environmental reaction norm models (Jarquín framework) for multi-site genomic prediction in coastal British Columbia.

The compiled output should serve as a reference base for the manuscript introduction, methods justification, and discussion sections, targeting journals such as *Heredity* or *Molecular Breeding*.

---

## Predecessor Studies (same dataset / breeding program)

The current manuscript builds directly on three prior publications from the same research group using the same EP.708 Douglas-fir breeding population and exome capture genotyping platform. The compiled literature review must position the new work relative to these findings:

### Thistlethwaite et al. (2017) — BMC Genomics 18:930
**"Genomic prediction accuracies in space and time for height and wood density of Douglas-fir using exome capture as the genotyping platform"**

- Training population: 1,372 trees, 37 full-sib families (54 parents), 3 sites (Adams, Fleet River, Lost Creek) in coastal BC, Ne ≈ 21
- Genotyping: Exome capture via RAPiD Genomics (40K probes designed from Douglas-fir transcriptome assembly; Howe et al. 2013), variant calling with freeBayes, yielding 69,551 SNPs after filtering (MAF >5%, HWD < −0.05, max depth <60, <40% missing)
- GS methods: RR-BLUP and GRR; cross-validated within-site, cross-site, multi-site, multi-site→single-site, and time-time (age 12→35 for height; height→wood density)
- Traits: HT12, HT35 (height at ages 12 and 35), WDres (wood density via Resistograph)
- Key findings:
  - High within-site prediction accuracies with EBVs (RR-BLUP: 0.79–0.91), comparable to ABLUP
  - Cross-site predictions surprisingly high, suggesting limited GxE
  - **When DEBVs (deregressed EBVs, parental averages removed) were used, prediction accuracies dropped to ~0 for all traits and cross-validation schemes** — demonstrating that GS models were tracking pedigree (family means), not marker-QTL LD
  - Height at age 12 identified as earliest acceptable age for forward prediction
  - Concluded that many more markers are needed to capture LD in the large Douglas-fir genome; exome-only markers could not resolve LD
  - The exome represents a very small fraction of the ~18.7 Gbp Douglas-fir genome; intergenic regions (not captured by exome) may contain important regulatory variants

### Thistlethwaite et al. (2019) — Heredity 122:848–863
**"Genomic selection of juvenile height across a single-generational gap in Douglas-fir"**

- Extended the 2017 study to cross-generational validation: F1 training population (1,321 trees, same 37 families, 3 sites) → F2 validation population (136 individuals at Jordan River, BC)
- EBVs estimated using full pedigree information (N = 36,311 from 11 F1 environments + P0 progenitors + 2 F2 environments)
- GS methods: RR-BLUP, GRR, and Bayes-B
- Four scenarios tested: F1 HTJ EBVs → F2 HTJ GEBVs; same with DEBVs; F1 HT35 EBVs → F2 HTJ GEBVs; same with DEBVs
- Key findings:
  - GS cross-generational accuracy for EBVs matched ABLUP (0.91–0.92 vs. 0.92)
  - **DEBVs again dropped to ~0**, confirming pedigree tracking not LD
  - Validation set showed two distinct familial clusters (EBV < 20 and EBV > 20), inflating overall correlations; within-cluster accuracies were much lower
  - Bayes-B outperformed RR-BLUP/GRR only in specific cluster subsets
  - Concluded: "Without capturing LD, GS cannot surpass the prediction of ABLUP. More markers or improved distribution of markers are required."

### Thistlethwaite, Gamal El-Dien, Ratcliffe et al. (2020) — PLoS ONE 15(6):e0232201
**"Linkage disequilibrium vs. pedigree: Genomic selection prediction accuracy in conifer species"**

- Compared Douglas-fir (full-sib, Ne ≈ 21, 1,321 trees, 37 families, exome capture, 56,454 SNPs) vs. Interior spruce (half-sib, Ne ≈ 93, 1,126 trees, 25 OP families, GBS, 62,190 SNPs)
- Random marker subsampling from 200–50,000 SNPs, each replicated 10 times, tested with RR-BLUP
- Key findings:
  - Prediction accuracy plateaued at ~10,000–15,000 markers for both species
  - Little variation in accuracy across random marker subsets — evidence that SNPs track relatedness, not specific marker-QTL LD
  - Douglas-fir full-sib accuracies consistently higher than Interior spruce half-sib (1.58× for height, 1.54× for wood density), reflecting Ne and family structure differences
  - **"Significant LD between markers and putative causal variants was not detected using 50,000 SNPs, and GS was enabled only through the tracking of relatedness"**
  - Genic regions (exome) may have faster LD decay than intergenic regions; intergenic LD patterns are virtually unknown in conifers
  - Recommended dramatically increasing marker density for cross-population GS utility

### What the current manuscript adds beyond these three papers

The new study extends this body of work by:
1. Implementing **ssGBLUP** (single-step GBLUP) rather than two-step RR-BLUP/GRR — integrating the ~1,870 genotyped individuals with the full ~15,000-tree pedigree network via the blended H matrix
2. Expanding from 3 genotyped sites to **5 genotyped sites** embedded in a 20-site phenotypic network (EP.708 series)
3. Adding **environmental reaction norm models** (Jarquín 2014 framework) via BGLR, using ClimateNA-derived covariates for all 20 sites
4. Including a **SNP panel optimization** component: AIM-selected, GWAS-selected, and random panels at multiple density levels
5. Using **three cross-validation schemes** (random 5-fold, leave-one-family-out, leave-one-site-out) to distinguish genetic vs. environmental prediction
6. Incorporating **spatial adjustment** (breedR AR1×AR1 or P-spline) in a formal two-stage pipeline, rather than using ASReml-derived EBVs directly

---

## Research Scope

For each topic below, identify:
- Foundational / seminal papers
- Key methodological advances (emphasis on 2018–2025)
- Active debates or unresolved questions
- Knowledge gaps relevant to Douglas-fir GS specifically

Cite specific authors and years throughout. Prioritize peer-reviewed sources.

---

## Topic 1: Genomic Prediction — Foundations and Current State

Compile literature on the development of genomic prediction / genomic selection (GS) from its origins through the current state of the art.

### Key areas to cover

- **Origins**: Meuwissen, Hayes & Goddard (2001) — the GS concept; early implementations in dairy cattle (e.g., VanRaden 2008; Hayes et al. 2009)
- **Statistical frameworks**: GBLUP, BayesA/B/Cπ/R, RKHS, Bayesian regression via BGLR (Pérez & de los Campos 2014); relative performance across genetic architectures (many small-effect loci vs. few large-effect QTL)
- **Genomic relationship matrix (G)**: VanRaden (2008) Methods 1 and 2; allele frequency centering; scaling; properties relative to the pedigree-based A matrix
- **Single-step GBLUP (ssGBLUP)**: Theory and implementation — Legarra et al. (2009), Aguilar et al. (2010), Christensen & Lund (2010); construction of the blended H matrix; weighting of G vs. A (τ and ω scaling); advantages for populations with mixed genotyped/ungenotyped individuals
- **Training population design**: Effects of training population size, relatedness to selection candidates, and number of generations on prediction accuracy; persistence of accuracy across generations
- **Marker density**: Minimum SNP requirements; diminishing returns at high density; interaction with effective population size (Ne) and LD extent
- **Additive vs. non-additive models**: Dominance and epistasis in GS; SCA estimation; practical importance for clonally deployed species vs. seedling deployment
- **Two-stage vs. one-stage modeling**: The trade-off between computational tractability (spatial adjustment → genomic prediction) and statistical efficiency (fully integrated models); bias implications of the two-stage approach

### Specific relevance to this project
The study uses ssGBLUP with a two-stage approach (breedR spatial adjustment, then BGLR genomic prediction), ~5–12K SNPs from exome capture sequencing (variant-called with freeBayes against a scaffold-level reference), a three-generation pedigree with ~1,870 genotyped individuals embedded in a ~15,000-tree phenotypic network across 20 trial sites.

---

## Topic 2: Genomic Prediction in Conifers

Summarize the trajectory of GS research and implementation in conifer tree breeding programs.

### Key species and studies to cover

- **Loblolly pine** (*Pinus taeda*): Resende et al. (2012a, 2012b) — among the first conifer GS studies; subsequent work by Zapata-Valenzuela, Isik, de los Campos, Munoz et al.
- **White spruce** (*Picea glauca*): Beaulieu et al. (2014); Lenz et al. (2017, 2020); Isabelle Gagnon's group
- **Interior spruce** (*Picea glauca × engelmannii*): Gamal El-Dien et al. (2015); Ratcliffe et al. (2015); directly compared with Douglas-fir in Thistlethwaite et al. (2020) — half-sib OP families, GBS genotyping, Ne ≈ 93, lower prediction accuracies than Douglas-fir full-sibs
- **Norway spruce** (*Picea abies*): Chen et al. (2018) — also used exome capture as genotyping platform, relevant comparison to the Douglas-fir exome capture approach; Lenz et al.
- **Scots pine** (*Pinus sylvestris*): Calleja-Rodriguez et al. (2019, 2020)
- **Maritime pine** (*Pinus pinaster*): Bartholomé et al. (2016); Isik et al.
- **Radiata pine** (*Pinus radiata*): Li et al.; Dungey et al.
- **Douglas-fir** (*Pseudotsuga menziesii*): Thistlethwaite et al. (2017, 2019, 2020) — the direct predecessors to the current study (see Predecessor Studies section above); Howe et al. (2013) transcriptome/SNP resource; any other published GS work in Douglas-fir. Key finding from the predecessor studies: GS prediction accuracy in Douglas-fir was driven entirely by pedigree relatedness tracking, not marker-QTL LD, even with ~70K exome capture SNPs. This motivates the current study's use of ssGBLUP and environmental covariates.

### Key themes

- Reported prediction accuracies for growth traits (height, diameter, volume) and wood quality traits (density, MFA, stiffness, Pilodyn)
- Conifer-specific challenges for GS:
  - Long generation intervals (15–30+ years to reproductive maturity)
  - Large effective population sizes (Ne often 50–500+ in breeding populations) causing rapid LD decay
  - High genetic load and inbreeding depression
  - Open-pollinated (OP) vs. controlled-cross (polymix, full-sib diallel) mating designs and their interaction with GS accuracy
  - Unbalanced, multi-site trial networks with incomplete genetic connectivity
  - Generally modest prediction accuracies compared to livestock (typically r = 0.2–0.5 for growth traits)
- Whether low-density marker panels (< 10K SNPs) are sufficient given rapid LD decay in conifers
- Operational deployment: has any conifer breeding program moved to routine GS-based selection?
- ssGBLUP vs. two-step GBLUP comparisons in conifers specifically
- **The pedigree vs. LD debate in conifer GS**: Multiple studies (Thistlethwaite et al. 2017, 2020; Lenz et al. 2017; Fuentes-Utrilla et al. 2017) have shown that GS prediction accuracy in conifers is largely driven by relatedness tracking rather than marker-QTL LD. Deregression of EBVs (removing parental averages) causes prediction accuracy to collapse. This has implications for: (a) the practical value of GS for within-family selection vs. among-family ranking; (b) the persistence of GS accuracy across generations; and (c) the marker density needed to capture LD in species with rapid LD decay and very large genomes

### Specific relevance to this project
The study is among the first to combine ssGBLUP with environmental reaction norm models in Douglas-fir. The EP.708 trial series uses controlled crosses (partial diallels), and the genotyped subset spans two generations across five trial sites. The three predecessor studies (Thistlethwaite et al. 2017, 2019, 2020) established that ~70K exome capture SNPs track pedigree but not LD in this population — a key motivation for the current study's ssGBLUP approach (which leverages both pedigree and genomic information simultaneously) and for exploring whether environmental covariates can extract additional predictive value beyond what relatedness tracking alone provides.

---

## Topic 3: Conifer Genomics — Genome Characteristics and Implications

Review the genomic features of conifers that create unique challenges for genotyping and genomic analysis.

### Genome architecture

- Exceptionally large genome sizes: 16–30 Gb (Douglas-fir ~16 Gb; loblolly pine ~22 Gb; sugar pine ~31 Gb)
- High repetitive element content (often >70% of the genome)
- Ancient whole-genome duplications (shared ancestral paleopolyploidy in Pinaceae) and residual paralogy
- Lack of chromosome-level reference assemblies for most species; scaffold-level or contig-level assemblies are the norm (Douglas-fir has a scaffold-level reference with thousands of contigs)
- Slow LD decay: typically within 1–10 kb in outbred conifer populations, compared to 10–100+ kb in inbred crop species
- High nucleotide diversity from obligate outcrossing, wind pollination, large census populations
- High genetic load, segregation distortion, and embryo lethality alleles

### Implications for genotyping and QC

- Paralog-collapsed reads: duplicated genomic regions mapping to the same reference location, inflating apparent heterozygosity
- Heterozygote miscalling at low sequencing depth (undercalling) and from paralog collapse (overcalling) — opposite biases that complicate QC
- SNP filtering strategies appropriate for conifers: MAF thresholds, call rate thresholds, excess heterozygosity ratios (and why aggressive filtering risks discarding genuine biological signal)
- Allele frequency estimation: importance of using the correct reference population (base/founder generation) to avoid the Wahlund effect when mixing generations
- HWE testing: expected departures from HWE in structured populations; generation-stratified testing
- GRM construction: sensitivity of VanRaden's G to allele frequency estimates; effects of paralog-contaminated markers on relationship estimation; implications of using exclusively genic SNPs (from exome capture) vs. genome-wide SNPs for relationship estimation and LD with causal variants
- Inbreeding coefficient (F) estimation in outcrossing species with high heterozygosity
- **Genic vs. intergenic LD patterns**: Thistlethwaite et al. (2020) noted that LD in conifer genic regions (the exome) may decay more rapidly than in intergenic regions (which are enriched for repetitive elements and have reduced recombination), but intergenic LD in conifers is virtually unstudied. This has direct implications for genotyping strategy: exome capture restricts SNPs to genic regions with potentially faster LD decay, possibly explaining the failure to detect marker-QTL LD with ~70K exome SNPs. Tan et al. (2017) found that intergenic SNPs provided slightly better prediction accuracy than genic SNPs in *Eucalyptus*.

### Key references
- Neale et al. (2014, 2017) — conifer genome reviews
- Birol et al. (2013) — white spruce genome assembly
- Stevens et al. (2016) — Douglas-fir draft genome
- Zimin et al. (2014, 2017) — loblolly pine and mega-genome assembly methods
- De La Torre et al. (2014) — conifer genome structure
- Pavy et al. — spruce gene catalogs

---

## Topic 4: Genotyping Methods — Exome Capture, GBS, and SNP Arrays

Compare genotyping approaches for non-model organisms with large, complex genomes, with emphasis on exome capture (the method used in this study).

### Exome capture / targeted sequence capture

This is the genotyping method used in the current study. Cover in depth:

- Probe design challenges in species without complete gene annotations; reliance on transcriptome assemblies and gene catalogs
- The exome capture workflow: probe hybridization, enrichment, sequencing; expected coverage profiles (more uniform than GBS but still variable across targets)
- Off-target capture and its utility (flanking intronic/intergenic sequence) vs. noise (repetitive elements, paralogs)
- Variant calling from capture data: freeBayes, GATK HaplotypeCaller; behavior with a scaffold-level reference assembly
- Depth-dependent genotype calling errors in capture data: heterozygote undercalling at low-coverage targets, paralog-collapsed reads at duplicated gene families inflating apparent heterozygosity — the same opposing biases seen in GBS but with different depth distributions
- QC pipelines for exome capture in conifers:
  - Staged QC architecture: removing failed samples before site-level QC, then refining sample QC on clean sites
  - The tension between stringent filtering (reducing genotyping error) and retaining sufficient markers
  - Excess heterozygosity filtering thresholds and the risk of discarding genuine biological signal in outcrossing conifers
  - Dodds et al. (2015) KGD concepts (fin plots, depth-adjusted relationship estimation) — applicability to capture data
- Cost-effectiveness for large breeding populations (hundreds to thousands of individuals)
- Examples in conifers: spruce exome capture panels (Pavy et al.; Suren et al. 2016), pine exome panels (Neves et al. 2013), Douglas-fir exome capture (the current study's platform — RAPiD Genomics, 40K probes from transcriptome assembly, ~70K SNPs post-filtering, as described in Thistlethwaite et al. 2017), Norway spruce exome capture (Chen et al. 2018)
- Comparison with GBS: more uniform per-locus coverage but restricted to genic regions; implications for capturing LD with causal variants in regulatory vs. coding regions
- Missing data handling: imputation approaches for capture data with moderate missingness

### Genotyping-by-sequencing (GBS)

Cover for comparative context (GBS is the most common alternative in conifer breeding programs):

- Elshire et al. (2011) — the original GBS protocol
- Reduced representation principle: restriction enzyme digestion to sample a reproducible fraction of the genome
- Expected coverage distributions: typically low and highly variable (median DP 5–15x)
- Depth-dependent genotype calling errors: heterozygote undercalling at low depth (allelic dropout); why per-genotype DP/GQ masking can be catastrophically aggressive at low coverage
- Variant callers for GBS data: freeBayes, GATK, TASSEL-GBS, Stacks
- QC challenges specific to GBS in complex genomes
- Missing data: higher missingness rates than capture or arrays; imputation requirements

### SNP arrays

- Axiom and Infinium array platforms used in conifers (e.g., the 50K spruce SNP array, loblolly pine arrays)
- Advantages: uniform call rates, reproducible genotyping, established QC pipelines
- Disadvantages: ascertainment bias, species-specificity, development cost

### Cross-cutting themes

- The trade-off between marker density, genotyping cost, and GS accuracy in species with rapid LD decay
- **Marker density plateau**: Thistlethwaite et al. (2020) found prediction accuracy plateaued at ~10,000–15,000 markers for both Douglas-fir and Interior spruce, with minimal variation across random marker subsets — evidence that markers capture relatedness rather than LD. The Grattapaglia & Resende (2011) deterministic formula (10×Ne×L markers) gives ~4,200 for Douglas-fir (Ne ≈ 21, L ≈ 20 Morgans) but this assumes LD-based prediction, which was not achieved
- How genotyping strategy interacts with genomic prediction model choice
- Douglas-fir-specific genotyping: what platforms or methods have been used? The current BC breeding program study uses exome capture with freeBayes variant calling against a scaffold-level reference; identify any prior Douglas-fir genotyping platform publications

---

## Topic 5: Environmental Covariates and Reaction Norm Models in Genomic Prediction

Review the integration of environmental information into genomic prediction frameworks.

### Foundational work

- **Jarquín et al. (2014)**: The enviromic kernel approach — constructing environmental relationship matrices (Ω or W matrices) and their Kronecker products with genomic relationship matrices for modeling GxE as a continuous reaction norm
- **Lopez-Cruz et al. (2015)**: Bayesian GxE reaction norms via BGLR
- **Pérez-Rodríguez et al.**: Extensions and comparisons of reaction norm models
- **Burgueño et al. (2012)**: Multi-environment GS models

### Environmental covariate sources and tools

- **ClimateNA** (Wang et al. 2016): Downscaled climate normals and variables for North American locations; widely used in forestry
- **EnvRtype** (Costa-Neto et al. 2021): R package for typifying environments, building environmental relationship kernels, and integrating with GS models
- **NASAPOWER, daymetr, and other weather data sources** for crop applications

### Reaction norm model frameworks

- Discrete-site GxE models (site as a factor): compound symmetry, factor analytic, and unstructured G×E covariance structures
- Continuous reaction norm models:
  - Linear Jarquín reaction norm: G × W interaction with environmental covariates entering linearly
  - Gaussian enviromic kernel: nonlinear environmental similarity matrices
- The distinction between modeling GxE via environmental covariates vs. estimating site-specific genetic effects
- Random regression / character-state approaches as alternatives

### Applications by domain

- **Crops**: Maize and wheat in the CIMMYT breeding program (Jarquín et al. 2014; Crossa et al. 2016, 2017; Costa-Neto et al. 2021); soybean, rice, other
- **Livestock**: Any enviromic or climate-covariate approaches in animal breeding
- **Forest trees**: Reaction norm or enviromic approaches in conifer GS:
  - Marchal & Cumbie (2024?) in loblolly pine (if published)
  - Any others in Douglas-fir, spruce, eucalyptus, or other tree species
  - Li et al. in radiata pine GxE modeling
  - Conventional GxE studies in conifers (non-genomic) that inform the expectation of GxE magnitude

### Key methodological questions

- **Environmental variable selection**: Which variables (temperature, precipitation, frost-free days, continentality, degree-days, etc.) are most predictive for tree growth traits? Dimensionality reduction approaches (PCA, partial correlations) for collinear climate variables
- **Narrow environmental gradients**: When sites span a limited climatic range (e.g., coastal BC with relatively homogeneous maritime climate), is there sufficient GxE signal to justify reaction norm modeling over simpler discrete-site models?
- **Cross-validation for GxE**: Leave-one-site-out vs. leave-one-family-out vs. random fold; what each tests (environmental extrapolation, genetic extrapolation, interpolation); which is most relevant for breeding program deployment
- **Model comparison**: Metrics for comparing ABLUP, GBLUP, ssGBLUP with discrete sites, ssGBLUP with linear reaction norm, ssGBLUP with Gaussian enviromic kernel
- **Enviromic kernel interpretation**: How to extract biological meaning from environmental similarity matrices; which environmental windows (growing season? dormancy period? establishment year?) matter most

---

## Output Format Requested

For each of the five topics, provide:

1. **Narrative synthesis** (not just a list of citations) organized by sub-theme, with key findings and conclusions stated explicitly
2. **Key citations** with author, year, journal, and a one-sentence summary of each paper's contribution
3. **Knowledge gaps and open questions** most relevant to a Douglas-fir GS study using exome capture, ssGBLUP, and environmental reaction norms
4. **Connections across topics** — how findings in one area inform or constrain another (e.g., how conifer genome complexity [Topic 3] affects genotyping method choice [Topic 4] which in turn constrains marker density for GS [Topics 1–2])
5. **Positioning relative to predecessor studies** — for each topic, note how the compiled findings relate to the key conclusions of Thistlethwaite et al. (2017, 2019, 2020), and identify where the current study's innovations (ssGBLUP, environmental reaction norms, expanded site network, SNP panel optimization) address gaps or limitations identified in those papers

Target length: comprehensive but focused. Approximately 3,000–5,000 words per topic, or the equivalent depth.
