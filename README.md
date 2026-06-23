# Papua New Guinea Rodent 18S Dietary Metabarcoding

**Dietary niche evolution in Papua New Guinea rodents assessed by 18S rRNA metabarcoding**  
819 individuals · 44 species · 2 mountain transects (Mt. Wilhelm & Saruwaged Range) · From sea level to 3700 m asl.

## Important note 
`PNG_18S_final_script_17_06_2026.Rmd` is the original, complete analysis pipeline. Note that this script reflects an iterative analytical workflow; parameters and intermediate objects were revised repeatedly during analysis, so some sections are not guaranteed to run sequentially from start to finish without modification.

To have a clean organisation, individual R scripts in `statistical_analyses/` were reorganised from the original Rmd using Claude Code, to provide a clearer, modular overview of each analysis. These scripts are intended to convey the analytical workflow, parameter choices, and overall logic of each analysis — they are not guaranteed to be directly executable as standalone, copy-paste-and-run scripts.

**If there is any contradiction between an individual script and the Rmd file, the Rmd file should be treated as the authoritative record of the analysis as performed.** 

Scripts are numbered in the order the analyses appear in the manuscript.

---

## Repository contents

```
.
├── PNG_18S_final_script_17_06_2026.Rmd       ← original analysis script 
├── README.md
├── bioinformatics_scripts/
│   ├── CRABS_get_ref.qsub                    ← CRABS reference database construction
│   ├── 01-demultiplexing_trimming/           ← Skewer demultiplexing + primer trimming
│   │   ├── primers_18S.fasta                 ← 18S primer sequences for Skewer
│   │   └── matrix_18S.txt                    ← barcode–primer match matrix for Skewer
│   ├── 02-merge_reads_create_ASV_table/      ← DADA2 quality filtering + ASV table
│   └── 03-assign-taxonomy/                   ← DADA2 taxonomy assignment (CRABS DB)
├── statistical_analyses/
│   ├── 01_phyloseq_build.R                   ← full phyloseq build chain (prerequisite)
│   ├── 02_sample_distribution.R              ← species altitude range + habitat figure (Fig 2)
│   ├── 03_abdomen.R                          ← phylosymbiosis via ABDOMEN (Fig 3)
│   ├── 04_niche_partitioning.R               ← indicspecies WIC/TNW/BIC/SI + PCoA (Fig 4)
│   ├── 05_asr_glmm.R                         ← ASR + Pagel's λ + GLMMs (Figs 8–9)
│   ├── 06_maaslin2.R                         ← MaAsLin2 differential abundance (Figs 5–6)
│   ├── 07_cooccurrence.R                     ← pairwise dietary co-occurrence LMMs (Fig 7), *script will be added*
│   └── 08_barplots.R                         ← dietary barplots per ensemble (Supp Figs 1–7)
├── data/
│   ├── phyloseqs/
│   │   ├── phyloseq_initial.R                ← phyloseq after taxonomy filtering (827 samples)
│   │   └── phyloseq_final.R                  ← final analysis phyloseq (819 samples, 44 species)
│   └── ddRAD_trees/
│       └── Joined_scaled_RH85.nwk            ← ddRAD host phylogeny (Rattini + Hydromyini)
└── supplementary_tables/
    ├── SuppTable1_species_diet_kingdom.csv   ← relative abundance per species at kingdom level
    ├── SuppTable2_species_diet_phylum.csv    ← relative abundance per species at phylum level
    ├── SuppTable3_species_diet_class.csv     ← relative abundance per species at class level
    ├── SuppTable4_species_diet_order.csv     ← relative abundance per species at order level
    ├── SuppTable5_species_diet_family.csv    ← relative abundance per species at family level
    ├── SuppTable6_ASV_individual_abundance.csv ← raw ASV counts per individual sample
    └── SuppTable7_sample_metadata.csv          ← full metadata for all 819 samples (162 fields)
```



## Bioinformatics pipeline

Scripts were run on an HPC cluster using PBS job scheduling (`.qsub` files).  

### Step 0 — Reference database construction (`CRABS_get_ref.qsub`)

A custom 18S reference database was built using **CRABS** (Creating Reference databases for Amplicon-Based Sequencing). The full pipeline is documented in `CRABS_get_ref.qsub`; only the download step ran automatically in the PBS job — the remaining steps were executed manually on the HPC.

| Sub-step | CRABS command | Key parameters |
|---|---|---|
| Download | `crabs db_download` | NCBI nucleotide; query: `18S OR SSU`; size 1–50,000 bp |
| Import | `crabs db_import` | `--seq_header species` or `accession` |
| Merge | `crabs db_merge` | `--uniq yes` |
| In silico PCR | `crabs insilico_pcr` | Primers: `GATYTGTCTGGTTVATTCCG` / `CATCACAGACCTGTTATYGC`; `--error 4.5` |
| Assign taxonomy | `crabs assign_tax` | NCBI `nucl_gb.accession2taxid`, `nodes.dmp`, `names.dmp` |
| Dereplicate | `crabs dereplicate` | `--method uniq_species` |
| Sequence cleanup | `crabs seq_cleanup` | `--minlen 80 --maxlen 500 --maxns 0 --nans 2` |
| Format for DADA2 | `crabs tax_format` | `--format dada2` |

The resulting reference database (`REF_dadaB_BlastFIlt.fasta`) was additionally filtered by BLAST similarity (≥85% threshold) before use in taxonomy assignment.

HPC resources: 1 node × 2 CPUs, 64 GB RAM, up to 24 h walltime (download step only).

### Step 1 — Demultiplexing (`01-demultiplexing_trimming/`)

Dual-index demultiplexing with **Skewer v0.2.2** using a primer FASTA file and primer match matrix (both provided in `01-demultiplexing_trimming/`):

```bash
skewer -x primers_18S.fasta -M matrix_18S.txt -b -m head -k 35 -d 0 -t 8 \
    sample_R1.fq.gz sample_R2.fq.gz -o output_prefix
```

Key parameters: `-m head` (trim from read head), `-k 35` (min length after trimming), `-d 0` (no mismatches in barcode).

### Step 2 — Primer trimming (`01-demultiplexing_trimming/`)

18S rRNA primer sequences trimmed from demultiplexed reads with Skewer. Inline barcodes represented as 4×N:

| Primer | Sequence |
|---|---|
| Forward (18SF) | `GATYTGTCTGGTTVATTCCG` |
| Reverse (18SR) | `CACAGACCTGTTATYGC` |

```bash
skewer -x NNNNGATYTGTCTGGTTVATTCCG -y NNNNCACAGACCTGTTATYGC \
    -m head -k 35 -d 0 -t 8 sample-pair1.fastq.gz sample-pair2.fastq.gz -o output
```

### Step 3 — Quality filtering, denoising & ASV table (`02-merge_reads_create_ASV_table/`)

Run with **DADA2** in R (`dada2_PNG_18S_novo[1-3]_new.R`):

| Step | Function | Key parameters |
|---|---|---|
| Quality filter | `filterAndTrim()` | `maxN=0`, `maxEE=c(2,2)`, `truncQ=2`, `minLen=50` |
| Dereplicate | `derepFastq()` | — |
| Denoise | `dada()` | `selfConsist=TRUE`, `MAX_CONSIST=25` |
| Merge pairs | `mergePairs()` | `minOverlap=10`, `maxMismatch=1`, `justConcatenate=FALSE` |
| ASV table | `makeSequenceTable()` | — |

HPC resources: 4 nodes × 12 CPUs, 64 GB RAM, up to 60 h walltime.

### Step 4 — Taxonomy assignment (`03-assign-taxonomy/`)

Taxonomic classification with `dada2::assignTaxonomy()` against a **CRABS**-filtered reference database (`REF_dadaB_BlastFIlt.fasta`). The CRABS pipeline retrieves sequences from GenBank/BOLD, filters by BLAST similarity (≥85% threshold), and formats them for DADA2.

```r
taxa <- assignTaxonomy(refseqs, "REF_dadaB_BlastFIlt.fasta",
                       multithread = FALSE, minBoot = 50)
```

HPC resources: 1 node × 24 CPUs, 256 GB RAM, up to 12 h walltime.

---

## Statistical analysis pipeline

### Data files

| File | Description |
|---|---|
| `data/phyloseqs/phyloseq_initial.R` | Phyloseq after duplicate merging + taxonomy filtering; 4726 taxa, 827 samples. Load with `load("data/phyloseqs/phyloseq_initial.R")`. |
| `data/phyloseqs/phyloseq_final.R` | Final analysis phyloseq: ensemble metadata, 18S ASV tree, cleaned habitats, *R. norvegicus* removed; 4723 taxa, 819 samples, 44 species. Load with `load("data/phyloseqs/phyloseq_final.R")` — object name is `PNG_18S_rodent_filtered_ensembles_tree_NEWMT_W_clean`. Alternatively, `01_phyloseq_build.R` saves it as `phyloseq_final.rds` (load with `readRDS()`). |
| `data/ddRAD_trees/Joined_scaled_RH85.nwk` | Scaled ddRAD host phylogeny (merged Rattini + Hydromyini); used by `02_abdomen.R` and `04_asr_glmm.R`. |

### Prerequisites

- R ≥ 4.x
- Key packages: `phyloseq`, `vegan`, `indicspecies`, `glmmTMB`, `MuMIn`, `Maaslin2`,  
  `phytools`, `ggtree`, `ape`, `dplyr`, `ggplot2`, `lme4`, `lmerTest`, `emmeans`,  
  `sjPlot`, `officer`, `webshot2`, `doSNOW`, `foreach`, `mvMORPH`, `RPANDA`, `rstan`,  
  `scatterpie`, `patchwork`, `ggrepel`
- `ABDOMEN.R` source file (required by `02_abdomen.R`): see [github.com/BPerezLamarque/ABDOMEN](https://github.com/BPerezLamarque/ABDOMEN)

### Running the scripts

Set the `data_dir`, `phyloseq_dir`, and `output_dir` path variables at the top of each script before running. Scripts are numbered in the order the analyses appear in the Methods section.

**Recommended order:**

| Script | Manuscript figures | Input | Description |
|---|---|---|---|
| `01_phyloseq_build.R` | — | DADA2 outputs | Builds `phyloseq_final.rds` from scratch (requires original sequencing data) |
| `02_sample_distribution.R` | Fig 2 | `phyloseq_final.R` | Species altitudinal range + habitat occurrence, faceted by mountain range |
| `03_abdomen.R` | Fig 3 | `phyloseq_final.rds`, ddRAD tree, `ABDOMEN.R` | Phylosymbiosis (Pagel's λ) via ABDOMEN Bayesian model; ancestral diet tree |
| `04_niche_partitioning.R` | Fig 4, Supp Fig 9 | `phyloseq_final.rds` | indicspecies WIC/TNW/BIC/SI; rarefied + phylogenetic distances; PCoA of niche overlap |
| `05_asr_glmm.R` | Figs 8–9, Supp Fig 10 | CSVs from `04`, ddRAD tree | ASR (fastAnc + chronos), Pagel's λ, species-level LM, individual-level Gamma-GLMM |
| `06_maaslin2.R` | Figs 5–6 | `phyloseq_final.rds` | MaAsLin2 differential abundance at dietary order level (two transects) |
| `07_cooccurrence.R` | Fig 7 | `phyloseq_final.rds` | Pairwise Bray-Curtis + Weighted UniFrac; 2×2 factorial LMM; emmeans + brackets |
| `08_barplots.R` | Supp Figs 1–7 | `phyloseq_final.rds` | Dietary composition barplots per ecological ensemble (top 25 orders) |

---

## Supplementary tables

| File | Description |
|---|---|
| `SuppTable1_species_diet_kingdom.csv` | Relative abundance of dietary items per species at kingdom level (Animal, Plant, Fungi) |
| `SuppTable2_species_diet_phylum.csv` | Relative abundance of dietary items per species at phylum level |
| `SuppTable3_species_diet_class.csv` | Relative abundance of dietary items per species at class level |
| `SuppTable4_species_diet_order.csv` | Relative abundance of dietary items per species at order level |
| `SuppTable5_species_diet_family.csv` | Relative abundance of dietary items per species at family level |
| `SuppTable6_ASV_individual_abundance.csv` | Raw ASV read counts with taxonomy per individual sample (not aggregated to species) |
| `SuppTable7_sample_metadata.csv` | Full sample metadata for all 819 individuals (162 fields): species ID, locality, altitude, habitat, spatial niche, tribe, taxonomic division, body mass, ecological ensemble, sequencing run, and all field/lab collection data |

---

## Citation

Solak et al. (2026). *Dietary niche evolution in Papua New Guinea rodents revealed by 18S rRNA metabarcoding.* Molecular Ecology. [in review]
