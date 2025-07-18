# One Health Microbiome + Exposomics Skill Roadmap

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

## 4. Pipeline Development & Workflow Management

Total Time: ~5 weeks  
Focus: Snakemake, Nexflow, Portability and Sharing

<details>
<summary>Know More</summary>

### 4.1. Workflow Architecture Fundamentals

| Sub-Skill               | Learn To...                                                | Notes                                                     |
| ----------------------- | ---------------------------------------------------------- | --------------------------------------------------------- |
| DAG thinking            | Design workflows as Directed Acyclic Graphs (DAGs)         | Know how input/output files drive task dependencies       |
| Rule chaining           | Define jobs that depend on outputs from prior rules        | Crucial for reproducible logic                            |
| Inputs, outputs, params | Write pipeline rules with dynamic wildcards and parameters | Know the difference between rule-level and global configs |
| Rule modularity         | Separate reusable rules into modules/snippets              | Encourages reusability across pipelines                   |
| Configuration files     | Use YAML/JSON to parameterize your pipeline                | Keeps code clean and flexible across datasets             |

### 4.2. Snakemake

| Sub-Skill                 | Learn To...                                                        | Notes                                                 |
| ------------------------- | ------------------------------------------------------------------ | ----------------------------------------------------- |
| Rule definition           | Write rules with `input`, `output`, `params`, `shell`, `resources` | Core of Snakemake pipelines                           |
| Wildcards and checkpoints | Handle variable filenames, sample-specific outputs                 | Use `{sample}` or `{group}` patterns                  |
| Config files              | Load sample sheets and params from `.yaml`                         | Essential for dataset-specific runs                   |
| Logging & reports         | Write `log:` blocks and generate workflow reports                  | `snakemake --report` for HTML summaries               |
| Conda/env integration     | Use `conda:` block to manage per-rule environments                 | Promotes reproducibility and avoids version conflicts |
| Cluster execution         | Use `--cluster sbatch` and profile configs                         | Integrates seamlessly with SLURM and SGE              |

### 4.3. Nextflow

| Sub-Skill                        | Learn To...                                              | Notes                                               |
| -------------------------------- | -------------------------------------------------------- | --------------------------------------------------- |
| Processes and channels           | Write processes and connect them via channels            | Nextflow is channel-driven (dataflow paradigm)      |
| Input/output declaration         | Use `from`, `into`, `tuple`, `file()` for channel IO     | More explicit than Snakemake’s rule chaining        |
| Parameterization                 | Use `params.config` and `.nf` config files               | Supports environment switching, AWS profiles, etc.  |
| Docker & Singularity integration | Declare container for each process                       | Enables exact reproducibility and cloud portability |
| DSL2 modular structure           | Use modules, workflows, main.nf for large-scale projects | nf-core compliant structure                         |

### 4.4. Workflow Design Best Practices

| Sub-Skill                     | Learn To...                                                               | Notes                                          |
| ----------------------------- | ------------------------------------------------------------------------- | ---------------------------------------------- |
| Sample sheet integration      | Load `.csv` or `.tsv` sample metadata for looping rules                   | Important for automating per-sample operations |
| Reusability and encapsulation | Split logic into subworkflows or module files                             | Avoids massive monolithic scripts              |
| Dry-run and benchmarking      | Test rules without execution (`--dry-run`, `touch`)                       | Safe pre-run validation                        |
| Resource optimization         | Set `threads`, `memory`, and `runtime` per rule                           | Required for HPC scaling or SLURM efficiency   |
| Error handling and debugging  | Use `--rerun-incomplete`, `--printshellcmds`, and logging                 | For traceability and crash recovery            |
| Output structure              | Keep `results/`, `logs/`, `config/`, `scripts/`, and `workflow/` separate | Makes repos easier to navigate and share       |

### 4.5. Workflow Portability & Sharing

| Sub-Skill                           | Learn To...                                                                    | Notes                                                           |
| ----------------------------------- | ------------------------------------------------------------------------------ | --------------------------------------------------------------- |
| Containerized rules                 | Run pipelines with `--use-conda`, `--use-singularity`, `docker.enabled = true` | For portable pipelines                                          |
| Profiles for different systems      | Create cluster profiles for local, HPC, AWS                                    | Config switching via `--profile hpc` etc.                       |
| Publishing workflows                | Structure repos with `README.md`, `envs/`, `config/`, `workflow/`              | Enables GitHub distribution or nf-core compatibility            |
| Workflow documentation              | Write usage examples, environment setup instructions                           | Critical for collaboration, reproducibility, and reviewer trust |
| Integration with `make` or wrappers | Optional: use `Makefile` to wrap Snakemake or Nextflow commands                | For simplified command-line entry points                        |

</details>

---

## 5. Containerization & Environment Management

Total Time: ~5 weeks  
Focus: Conda, R, Python, Docker, Singularity

<details>
<summary>Know More</summary>

### 5.1. Conda & Environment Management (Cross-language)

| Sub-Skill                      | Learn To...                                                               | Notes                                    |
| ------------------------------ | ------------------------------------------------------------------------- | ---------------------------------------- |
| Creating environments          | Use `conda create -n env_name pkg1 pkg2`, `conda activate`                | Standard across Python/R pipelines       |
| Exporting & sharing            | `conda env export > env.yml`; recreate with `conda env create -f env.yml` | Always version-lock your environments    |
| Environment isolation          | Use different envs for different projects                                 | Prevent dependency conflicts             |
| Bioconda & channels            | Learn how to install tools from `bioconda`, `conda-forge`                 | e.g., `conda install -c bioconda fastqc` |
| Snakemake/Nextflow integration | Use `conda:` block per rule to autoinstall envs                           | Enables reproducibility per-step         |

### 5.2. R Environment Management with renv

| Sub-Skill               | Learn To...                                              | Notes                                        |
| ----------------------- | -------------------------------------------------------- | -------------------------------------------- |
| Project isolation       | Use `renv::init()` to create a local project environment | Tracks all installed packages in `renv.lock` |
| Dependency tracking     | Auto-record versions of every R package                  | Essential for reproducibility in notebooks   |
| Environment restoration | Use `renv::restore()` on a new system to rebuild the env | Like `conda env create`, but for R           |
| GitHub integration      | Commit `renv.lock`, ignore `renv/library`                | Makes repos portable but lightweight         |

### 5.3. Docker (Containers for Everything)

| Sub-Skill                   | Learn To...                                                               | Notes                                 |
| --------------------------- | ------------------------------------------------------------------------- | ------------------------------------- |
| Dockerfile creation         | Write `FROM`, `RUN`, `COPY`, `CMD` layers                                 | Builds a reproducible container image |
| Building and running images | `docker build -t mytool .`, `docker run -v $PWD:/data mytool`             | Mount volumes, pass arguments         |
| Versioned images            | Tag builds with version numbers (`:v1.0`, `:latest`)                      | Use in pipeline configs               |
| Base images                 | Choose wisely: `rocker/rstudio`, `continuumio/miniconda`, `biocontainers` | Ensures reproducibility across builds |
| Entrypoints & CMDs          | Create CLI wrappers using Python, R, or Bash                              | For command-line tool containers      |

### 5.4. Singularity (for HPC & Academic Clusters)

| Sub-Skill                        | Learn To...                                                       | Notes                                    |
| -------------------------------- | ----------------------------------------------------------------- | ---------------------------------------- |
| Running images                   | Use `singularity run mycontainer.sif`                             | Needed when Docker is unavailable on HPC |
| Building containers              | Convert Docker image: `singularity build out.sif docker://ubuntu` | Leverages DockerHub ecosystem            |
| Binding volumes                  | Use `--bind /path:/container/path` for data I/O                   | Required for filesystem access on HPC    |
| Snakemake & Nextflow integration | Use `--use-singularity` or process-level container declaration    | Fully portable pipelines                 |
| Managing versions                | Store `.sif` files with version info                              | Stable and audit-ready for publications  |

### 5.5. Container Registries & Distribution

| Sub-Skill                   | Learn To...                                                                | Notes                                     |
| --------------------------- | -------------------------------------------------------------------------- | ----------------------------------------- |
| DockerHub                   | Push/pull your Docker containers with `docker push`                        | Enables sharing tools publicly            |
| Biocontainers               | Use community-curated tools via `bioconda` or `docker://biocontainers/...` | Saves effort building from scratch        |
| GitHub Container Registry   | Optional: publish private/protected containers with GitHub Actions         | For organizations or controlled workflows |
| Registry tagging & security | Manage tokens, tags, and container metadata                                | Especially when deploying at scale        |

### 5.6. Integrated Environment Strategies

| Strategy                                              | Use Case                                       | Stack                                 |
| ----------------------------------------------------- | ---------------------------------------------- | ------------------------------------- |
| Conda-only                                            | Lightweight, flexible, simple pipelines        | `env.yml`                             |
| Docker + Conda                                        | Portability + dependency control               | Dockerfile + `env.yml`                |
| Singularity + Conda                                   | HPC-compatible reproducibility                 | `.sif` + `env.yml`                    |
| Full stack (Snakemake/Nextflow + Singularity + Conda) | Production-grade, fully reproducible workflows | Multi-config systems with YAML params |

</details>

---

## 6. Multi-Omics Data Processing

Total Time: ~12 weeks  
Focus: 16s, shotgun, WGS, metatranscriptomics, metaproteomics/metabolomics

<details>
<summary>Know More</summary>

### 6.1. Amplicon Sequencing (16S rRNA / ITS)

| Sub-Skill                      | Learn To...                                    | Tools                                            |
| ------------------------------ | ---------------------------------------------- | ------------------------------------------------ |
| Import and demultiplex         | Handle multiplexed read files, barcodes        | QIIME2, cutadapt, fastp                          |
| Denoising                      | Use ASV-level methods for error correction     | DADA2, Deblur                                    |
| Chimera filtering              | Identify and remove artifacts                  | DADA2 built-in, VSEARCH                          |
| Taxonomic classification       | Use reference databases to assign taxonomy     | SILVA, Greengenes, GTDB with Naive Bayes, IDTAXA |
| Phylogenetic tree construction | Align representative sequences and build trees | MAFFT + FastTree                                 |
| Diversity analyses             | Alpha/beta diversity, ordination, PERMANOVA    | QIIME2, phyloseq, vegan                          |

### 6.2. Shotgun Metagenomics

| Sub-Skill                      | Learn To...                               | Tools                              |
| ------------------------------ | ----------------------------------------- | ---------------------------------- |
| Preprocessing                  | QC, trimming, filtering low-quality reads | fastp, TrimGalore, FastQC, MultiQC |
| Host/decontamination filtering | Remove host or contaminant reads          | Bowtie2, BMTagger                  |
| Assembly (optional)            | Reconstruct contigs from reads            | MEGAHIT, SPAdes                    |
| Taxonomic profiling            | Generate species-level profiles           | MetaPhlAn3, Kraken2, Bracken       |
| Functional profiling           | Annotate gene families and pathways       | HUMAnN3 (UniRef90, MetaCyc)        |
| Binning (optional)             | Cluster contigs into MAGs                 | MetaBAT2, CONCOCT, MaxBin2         |

### 6.3. Metatranscriptomics

| Sub-Skill               | Learn To...                                            | Tools                        |
| ----------------------- | ------------------------------------------------------ | ---------------------------- |
| RNA-seq preprocessing   | Trimming, QC, rRNA depletion                           | TrimGalore, SortMeRNA, fastp |
| Mapping/quantification  | Map reads to reference or quantify via pseudoalignment | Salmon, Kallisto, BWA        |
| rRNA vs mRNA separation | Separate regulatory vs coding reads                    | SortMeRNA, custom filters    |
| Normalization           | TPM, RPKM, DESeq2 VST                                  | DESeq2, edgeR, limma-voom    |
| Differential expression | Compare across conditions                              | DESeq2, limma-voom           |
| Gene/pathway annotation | Map to GO, KEGG, eggNOG, MetaCyc                       | eggNOG-mapper, InterProScan  |

### 6.4. Metaproteomics

| Sub-Skill             | Learn To...                                              | Tools                           |
| --------------------- | -------------------------------------------------------- | ------------------------------- |
| Raw data processing   | MS file format conversion (.raw to .mzML)                | MSConvert                       |
| Database search       | Match spectra to UniProt sequences                       | FragPipe, MSFragger             |
| Quantification        | Protein/peptide intensities                              | IonQuant, MaxQuant              |
| Taxonomic assignment  | Peptide-based LCA profiling                              | UniPept, MetaProteomeAnalyzer   |
| Functional annotation | Map to GO, KEGG, EC, MetaCyc                             | eggNOG, InterProScan            |
| Output integration    | Generate peptide → protein → function abundance matrices | Custom R/Python parsing scripts |

### 6.5. Metabolomics / Exposomics

| Sub-Skill              | Learn To...                               | Tools                                                 |
| ---------------------- | ----------------------------------------- | ----------------------------------------------------- |
| Raw data preprocessing | Peak detection, alignment, deconvolution  | XCMS, MZmine                                          |
| Annotation & ID        | Match m/z features to chemical IDs        | GNPS, HMDB, MS-DIAL                                   |
| Normalization          | Log, Pareto, scaling by TIC               | `MetaboAnalystR`, `MSPrep`                            |
| Batch correction       | Combat, LOESS, internal standards         | `sva::ComBat`, `normStats`                            |
| Pathway enrichment     | Map m/z features to pathways              | Mummichog, MetaCyc, KEGG                              |
| Exposure integration   | Link to host response or microbiome shift | Statistical or ML-based integration with omics layers |

### 6.6. Cross-Omics Output Formatting

| Sub-Skill                        | Learn To...                                                     | Format                                            |
| -------------------------------- | --------------------------------------------------------------- | ------------------------------------------------- |
| Unified abundance tables         | Construct sample × feature matrix for taxa/proteins/metabolites | TSV, CSV, biom                                    |
| Metadata linkage                 | Join omics tables with sample-level metadata                    | `left_join`, `merge`, `sample_data()` in phyloseq |
| Feature mapping                  | Map gene/protein/peptide IDs to function, taxonomy              | UniProt ID mapping, eggNOG, InterProScan          |
| Output structuring for pipelines | Write reusable scripts for creating final tables                | R, Python, Snakemake rules                        |

</details>

---

## 7. Functional Annotation & Pathway Inference

Total Time: ~4 weeks  
Focus: Gene and protein functional annotation, pathway enrichment, metabolic reconstruction, ontology handling

<details>
<summary>Know More</summary>

### 7.1. Gene & Protein Functional Annotation

| Sub-Skill                          | Learn To...                                                                       | Tools & Resources                                           |
| ---------------------------------- | --------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| Map genes to orthologous groups    | Assign functional context across species                                          | **eggNOG-mapper**, **OrthoFinder**, **KEGG Orthology (KO)** |
| Assign GO terms (BP, MF, CC)       | Retrieve Gene Ontology Biological Process, Molecular Function, Cellular Component | **InterProScan**, **Blast2GO**, **UniProtKB**, **biomaRt**  |
| EC number annotation               | Assign Enzyme Commission codes to genes/proteins                                  | **eggNOG**, **KAAS**, **UniProt**                           |
| Filter high-confidence annotations | Remove generic or uninformative hits (e.g., “hypothetical protein”)               | Custom R/Python parsing logic                               |
| Functional redundancy analysis     | Quantify whether multiple features map to same function                           | Often important in microbial communities                    |

### 7.2. Pathway Annotation & Enrichment

| Sub-Skill                    | Learn To...                                                                                             | Tools & Resources                                                               |
| ---------------------------- | ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| Map functions to pathways    | Connect gene/protein/metabolite IDs to pathways                                                         | **KEGG**, **MetaCyc**, **Reactome**, **GOslim**                                 |
| Functional hierarchy mapping | Annotate multi-level (L1-L4) functions (e.g., metabolism → amino acid metabolism → lysine biosynthesis) | **HUMAnN3**, **eggNOG**, **UniRef90**, **KEGG modules**                         |
| Perform pathway enrichment   | Run Fisher’s exact test or GSEA-like approaches to find enriched functions                              | **clusterProfiler**, **gProfiler2**, **MetaboAnalyst**, **Pathview**, **DAVID** |
| Pathway coverage estimation  | Determine % of pathway covered by your observed genes                                                   | **HUMAnN3**, **MinPath**, **MetaPath**                                          |

### 7.3. Metabolic Reconstruction

| Sub-Skill                                      | Learn To...                                         | Tools & Resources                                        |
| ---------------------------------------------- | --------------------------------------------------- | -------------------------------------------------------- |
| Reconstruct metabolic maps                     | Build draft metabolic models from annotations       | **Pathway Tools**, **CarveMe**, **ModelSEED**, **AGORA** |
| Assign function to taxa (functional potential) | Assess whether taxon X encodes pathway Y            | **Tax4Fun2**, **PICRUSt2**, **HUMAnN3**                  |
| Link environment to function                   | Ask: “What does this community do in this context?” | Critical for exposomics and One Health research          |

### 7.4. Ontology Handling & Term Curation

| Sub-Skill                      | Learn To...                                                                                     | Tools & Resources                                          |
| ------------------------------ | ----------------------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| Standardize terms              | Use controlled vocabularies: GO, KEGG BRITE, MetaCyc                                            | **OntologyLookupService**, **OBO Foundry**, **Ontobee**    |
| Deduplicate IDs                | Merge multiple IDs mapping to the same function                                                 | Use UniProt cross-mapping tables                           |
| Visualize ontology graphs      | Show GO term hierarchies or relationships                                                       | **topGO**, **REVIGO**, **QuickGO**, **ggraph**, **igraph** |
| Custom annotation dictionaries | Create dictionaries of curated functional categories (e.g., “stress response”, “toxin-related”) | R/Python dictionaries, JSON term banks                     |

### 7.5. Output Structuring for Downstream Use

| Sub-Skill                           | Learn To...                                                 | Output Format                                        |
| ----------------------------------- | ----------------------------------------------------------- | ---------------------------------------------------- |
| Build functional abundance matrices | Sample × GO Term / EC / KO matrix                           | Used for DA, enrichment, ordination                  |
| Join taxonomy + function            | Create combined tables for taxon-function interpretation    | e.g., “Bacteroides → LPS biosynthesis”               |
| Store annotations as metadata       | Keep raw feature ID + annotation + abundance in tidy format | For use in MiNDSET, MoMo-MAP, Shiny dashboards, etc. |
| Create interpretable summaries      | Summarize functions by category, module, or subsystem       | Often shown in stacked barplots, pathway heatmaps    |

</details>

---

