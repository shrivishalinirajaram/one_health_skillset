# One Health Roadmap

# One Health Microbiome + Exposomics Skill Portfolio

This repository tracks my technical skill-building journey across 12 months (2 hours/day), focused on computational biology, exposomics, and microbiome science in the One Health framework.

---

## Categories

- [1. Biological Domain Knowledge](#1-biological-domain-knowledge)
- [2. Biostatistics & Quantitative Modeling](#2-biostatistics--quantitative-modeling)
- [3. Programming & Scripting](#3-programming--scripting)
- [4. Pipeline Development & Workflow Management](#4-pipeline-development--workflow-management)
- [5. Containerization & Reproducibility](#5-containerization--reproducibility)

---

## 1. Biological Domain Knowledge

Total Time: ~7 weeks  
Focus: Microbiome ecology, host–microbe–toxin interaction, One Health, exposure biology, CRISPR, TA systems

<details>
<summary>Know More</summary>

### 1.1. Microbiome Ecology

| Sub-Skill                           | Learn To...                                                           | Notes                                             |
| ----------------------------------- | --------------------------------------------------------------------- | ------------------------------------------------- |
| Microbial ecology principles        | Understand diversity, richness, evenness; define core/rare taxa       | Alpha/beta/gamma diversity; richness vs abundance |
| Functional guilds & metabolic roles | Interpret what microbes do (e.g., SCFA production, vitamin synthesis) | Use KEGG/MetaCyc modules to connect pathways      |
| Colonization resistance & dysbiosis | Identify when microbial communities shift toward pathogenic states    | Important for exposome impact assessment          |
| Cross-biome connections             | Compare oral–gut–skin–lung microbiomes                                | One Health relevance in zoonotic exchange         |


### 1.2. Host–Microbe Interactions

| Sub-Skill                      | Learn To...                                                           | Notes                                        |
| ------------------------------ | --------------------------------------------------------------------- | -------------------------------------------- |
| Mucosal immunity               | Understand immune education by commensals                             | Especially IL-10, Tregs, IgA                 |
| Microbe-induced signaling      | Analyze microbial metabolites (e.g., SCFA, bile acids) affecting host | Integration with metabolomics/exposomics     |
| Barrier function               | Understand epithelial integrity, tight junctions, and permeability    | Used in toxicology and inflammation contexts |
| Host gene expression responses | Link microbiome shifts to host transcriptomic changes                 | Required for multi-omics correlation work    |

### 1.3. One health systems thinking

| Sub-Skill                                | Learn To...                                                 | Notes                                                   |
| ---------------------------------------- | ----------------------------------------------------------- | ------------------------------------------------------- |
| Cross-species microbiome transmission    | Track microbial exchange between human–animal–environment   | Zoonoses, AMR transfer, ecological modeling             |
| Reservoirs and environmental persistence | Understand how soil, water, or air act as microbial vectors | Especially relevant for exposure sources                |
| Shared exposure consequences             | Compare how toxins affect different hosts                   | Use same metabolite pathways across hosts               |
| Surveillance and global monitoring       | Explore WHO/FAO/OIE approaches to microbial surveillance    | Helpful for integrative use cases or policy translation |

### 1.4. Exposomics and Environmental Toxicology

| Sub-Skill                        | Learn To...                                                                              | Notes                                                                      |
| -------------------------------- | ---------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| Exposure routes and kinetics     | Inhalation, dermal, oral; ADME principles                                                | Know how chemicals reach microbes and host tissues                         |
| Xenobiotic metabolism            | Identify microbial enzymes (e.g., azoreductase, nitroreductase) that transform chemicals | Impacts metabolomics interpretation                                        |
| Host–microbe–chemical crosstalk  | Recognize interaction effects (e.g., dysbiosis after exposure)                           | Synergistic or antagonistic outcomes matter                                |
| Microbial indicators of exposure | Use taxa as biomarkers of toxin presence                                                 | Differential abundance signatures or functional genes (e.g., efflux pumps) |

### 1.5. CRISPR, Toxin–Antitoxin (TA), and Other Mobile Elements

| Sub-Skill                           | Learn To...                                                     | Notes                                            |
| ----------------------------------- | --------------------------------------------------------------- | ------------------------------------------------ |
| CRISPR-Cas systems                  | Distinguish Class I vs II; identify arrays and cas operons      | Required for Spacerome, viral defense analysis   |
| Toxin–antitoxin systems             | Understand types (I–VI), addiction modules, and stress response | Metatranscriptomics of TA genes                  |
| Mobile genetic elements             | Detect plasmids, phages, ICEs                                   | Important in AMR spread, functional plasticity   |
| Horizontal gene transfer mechanisms | Conjugation, transformation, transduction                       | For exposure-driven evolution/adaptation studies |

### 1.6. Taxonomic, Phylogenetic & Evolutionary Foundations

| Sub-Skill                     | Learn To...                                                  | Notes                                                |
| ----------------------------- | ------------------------------------------------------------ | ---------------------------------------------------- |
| Prokaryotic taxonomy          | Understand NCBI/SILVA/GTDB hierarchies                       | Required for mapping and data standardization        |
| Strain-level diversity        | Differentiate between species, strains, and genotypes        | Applies in MAG analysis, shotgun metagenomics        |
| Phylogenetic inference        | Interpret phylogenies based on marker genes or whole genomes | Used for tree-based analysis, functional predictions |
| Evolution of microbial traits | Understand selective pressures, niche adaptation             | For interpreting functional enrichment results       |

</details>

---

## 2. Biostatistics & Quantitative Modeling

Total Time: ~8 weeks  
Focus: GLMs, ZINB, PERMANOVA, compositional stats, differential abundance, mixed models

<details>
<summary>Know More</summary>

### 2.1. Foundations of Statistical Thinking

| Sub-Skill                   | Learn To...                                                                        | Notes                                                  |
| --------------------------- | ---------------------------------------------------------------------------------- | ------------------------------------------------------ |
| Probability theory          | Interpret distributions, likelihoods, prior beliefs                                | Essential for understanding Bayesian models            |
| Key distributions           | Normal, Poisson, Binomial, Negative Binomial, Zero-Inflated, Dirichlet Multinomial | Needed for microbiome data and overdispersion modeling |
| Estimation & testing        | Perform point estimation, confidence intervals, hypothesis testing                 | Learn frequentist and Bayesian logic                   |
| Multiple testing correction | Apply FDR (Benjamini–Hochberg), Bonferroni, q-value                                | Required in multi-omic differential abundance tests    |

### 2.2. Generalized Linear and Mixed Models

| Sub-Skill                               | Learn To...                                                            | Notes                                                    |
| --------------------------------------- | ---------------------------------------------------------------------- | -------------------------------------------------------- |
| Linear models (LM)                      | Run `lm()` models with continuous outcomes                             | Basis for PCA, PERMANOVA, etc.                           |
| GLMs                                    | Use `glm()` for count/zero-inflated data (log, logit, identity link)   | Use family = `poisson`, `quasipoisson`, `binomial`, `nb` |
| Generalized linear mixed models (GLMMs) | Use `lme4::glmer`, `glmmTMB`, `nlme` to handle nested/repeated designs | For longitudinal exposome or repeated microbiome samples |
| Model diagnostics                       | Residuals, AIC, pseudo-R², overdispersion checks                       | Ensure you're not misinterpreting noise as signal        |

### 2.3. Zero-Inflation, Compositionality & Normalization

| Sub-Skill                 | Learn To...                                                | Notes                                                       |
| ------------------------- | ---------------------------------------------------------- | ----------------------------------------------------------- |
| Compositional data issues | Understand sum-constrained proportions, false correlations | Affects microbiome, metabolomics, exposome data             |
| Transformations           | Apply CLR, ILR, ALR, log+psuedocount, TSS                  | Use `microbiome::transform`, `compositions::clr`            |
| Normalization             | VST, RLE, CSS, GMPR, rarefaction                           | Learn when to use each based on your method                 |
| Zero-inflated modeling    | Choose ZINB, hurdle, or ZIBR based on dropout structure    | Use `glmmTMB`, `zinbwave`, `ZIBR`, `ANCOM-BC` appropriately |

### 2.4. Differential Abundance & Expression Modeling

| Sub-Skill                   | Learn To...                                                   | Notes                                                           |
| --------------------------- | ------------------------------------------------------------- | --------------------------------------------------------------- |
| Classical methods           | Apply t-test, ANOVA, Kruskal–Wallis, Wilcoxon                 | Only for simple comparisons — avoid overuse                     |
| RNA-seq-inspired models     | Use `DESeq2`, `edgeR`, `limma-voom`, `voomWithQualityWeights` | Base for protein, transcript, or taxonomic abundance shifts     |
| Compositional-aware methods | Use `ALDEx2`, `ANCOM-BC`, `MaAsLin2`, `metagenomeSeq`         | Choose based on effect size structure, covariates, and sparsity |
| Model comparison & post-hoc | Compare fit via AIC/BIC; apply post-hoc Tukey or emmeans      | Necessary when doing multi-group comparisons                    |

### 2.5. Multivariate & Ordination Methods

| Sub-Skill                 | Learn To...                                                   | Notes                                                   |
| ------------------------- | ------------------------------------------------------------- | ------------------------------------------------------- |
| PCA / PCoA / NMDS         | Reduce dimensionality and visualize beta-diversity            | Use `vegan::metaMDS`, `ape::pcoa`, `prcomp`             |
| Distance metrics          | Understand Bray-Curtis, Aitchison, Jaccard, Euclidean         | Choice affects ordination structure and interpretation  |
| PERMANOVA & Adonis        | Run `vegan::adonis` for group separation on distance matrices | Used in all microbiome diversity comparisons            |
| Procrustes & Mantel tests | Compare ordination structures across omics                    | Important for integrative exposomics/microbiome studies |

### 2.6. Longitudinal & Hierarchical Modeling

| Sub-Skill                | Learn To...                                                      | Notes                                                |
| ------------------------ | ---------------------------------------------------------------- | ---------------------------------------------------- |
| Repeated measures models | Mixed models (`lmer`, `glmmTMB`, `nlme`) for repeated timepoints | Use for exposome × time, diet interventions, etc.    |
| Time series smoothing    | Use `splines`, `smoothers`, `mgcv::gam` for trend detection      | Needed for exposome drift and adaptation modeling    |
| MaAsLin2, ZIBR, ZINBMM   | Use longitudinal microbiome/exposome-specific tools              | Part of MiNDSET development and interpretation logic |

### 2.7. Bayesian & Simulation-Based Inference

| Sub-Skill                     | Learn To...                                                      | Notes                                                                   |
| ----------------------------- | ---------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Bayesian modeling             | Use `brms`, `rstanarm`, `JAGS`, `Stan`                           | Full uncertainty modeling for complex biological systems                |
| Posterior estimation & priors | Understand credible intervals, shrinkage, regularization         | Required for ZIDM, eBay, and probabilistic exposome inference           |
| Simulated data testing        | Use `simstudy`, `synthpop`, or custom scripts for power analysis | Especially when working with underpowered animal/environmental datasets |

</details>

---

## 3. Programming & Scripting

Total Time: ~6 weeks  
Focus: R (tidyverse, Bioconductor), Python (pandas, Biopython), Bash (file parsing, job scripts)

<details>
<summary>Know More</summary>

### 3.1. R Programming for Statistical and Microbiome Analysis

| Sub-Skill                  | Learn To...                                                                   | Notes                                                |
| -------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------- |
| Tidyverse core             | Use `dplyr`, `tidyr`, `ggplot2`, `tibble`, `forcats` for clean pipelines      | Tidy data in, tidy data out                          |
| Bioconductor packages      | Use `phyloseq`, `DESeq2`, `edgeR`, `limma`, `vegan`, `microbiome`             | Essential for omics analysis                         |
| Plotting                   | Make publication-ready plots with `ggplot2`, `patchwork`, `cowplot`, `ggpubr` | Modular, themable, scalable plots                    |
| Data reshaping             | Use `pivot_longer`, `pivot_wider`, `separate`, `unite`                        | Clean metadata or taxonomic tables                   |
| Writing functions          | Wrap routines into reusable, modular functions                                | Critical for scaling code and pipelines              |
| Error handling & debugging | Use `tryCatch`, `message()`, `stop()`                                         | For large-scale batch processing and robust wrappers |
| Literate programming       | Write `RMarkdown`, `Quarto`, `.Rmd` reports                                   | For reproducible reports and notebooks               |

### 3.2. Python for Parsing, Preprocessing & Machine Learning

| Sub-Skill                | Learn To...                                                                    | Notes                                                       |
| ------------------------ | ------------------------------------------------------------------------------ | ----------------------------------------------------------- |
| Data manipulation        | Use `pandas`, `numpy`, `glob`, `os` to wrangle and process files               | Ideal for working with metadata, logs, sequence files       |
| Plotting                 | Use `matplotlib`, `seaborn`, `plotly` for dynamic plots                        | `sns.clustermap` for heatmaps, `plt.subplots` for panels    |
| Machine learning         | Use `scikit-learn` for RF, SVM, XGBoost, pipelines                             | For supervised learning & interpretable ML                  |
| Bioinformatics utilities | Use `BioPython`, `ete3`, `scikit-bio`, `PyPHLAWD`                              | FASTA/FASTQ parsing, phylogeny, alignment                   |
| API access & automation  | Use `requests`, `json`, `xml`, `BeautifulSoup` for scraping or data extraction | Programmatic data retrieval (e.g., SRA, MGnify, BioSamples) |
| Writing clean scripts    | Write `.py` modules, argparse CLI wrappers, logging                            | For building scalable tools and batch jobs                  |

### 3.3. Bash and Command-Line Proficiency

| Sub-Skill                 | Learn To...                                                        | Notes                                                        |
| ------------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------ |
| Navigation & file ops     | Use `cd`, `ls`, `mv`, `cp`, `rm`, `mkdir`, `find`, `xargs`, `tree` | For fast project traversal and data prep                     |
| File parsing              | Use `cut`, `awk`, `sed`, `sort`, `uniq`, `grep`, `head`, `tail`    | For log parsing, sample lists, metadata cleanup              |
| Job scheduling (HPC)      | Use `sbatch`, `qsub`, `squeue`, `sacct`, `#!/bin/bash` headers     | For running batch pipelines on clusters                      |
| Permissions & environment | Use `chmod`, `chown`, `PATH`, `export`, `.bashrc`                  | To avoid execution errors and dependency hell                |
| Writing shell scripts     | Automate routine jobs with `.sh` scripts and parameterized loops   | For reproducible project automation                          |
| Command-line utilities    | Install and use tools like `fastqc`, `multiqc`, `seqkit`, `jq`     | Widely used in pre/post processing for multi-omics workflows |

### 3.4. Inter-language Integration & Reusability

| Sub-Skill                         | Learn To...                                                               | Notes                                                   |
| --------------------------------- | ------------------------------------------------------------------------- | ------------------------------------------------------- |
| Call R from Python (`rpy2`)       | Seamlessly combine ggplot2 with scikit-learn models                       | Advanced, but powerful                                  |
| Call Python from R (`reticulate`) | Use machine learning tools within R pipelines                             | Use inside Quarto/Rmd                                   |
| Modular function writing          | Split analysis into reusable scripts or R/Py functions                    | Better for pipelines & reproducibility                  |
| Logging and messaging             | Use `logging` (Python), `message()` (R), `echo` (Bash) for robust outputs | Helps in debugging large runs                           |
| Config-driven scripts             | Load parameters via `.yaml`, `.json`, or `.toml`                          | Needed for Snakemake, Nextflow, or pipeline integration |

</details>

---

