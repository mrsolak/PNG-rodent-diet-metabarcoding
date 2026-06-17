################################################################################
# 01_phyloseq_build.R
# Full phyloseq build chain: DADA2 outputs → phyloseq_final.rds
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# Run this script once to produce phyloseq_final.rds.
# All downstream analysis scripts (02–06) load that single RDS file.
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# PREREQUISITES
#   R 4.x with: phyloseq, ShortRead, Biostrings, vegan, dplyr, spgs, ape
################################################################################


# ==============================================================================
# 0. PATHS — set these to match your system before running
# ==============================================================================

# Root directory where per-run DADA2 outputs are stored
dada2_root  <- "/path/to/dada2_outputs"

# Host metadata files
meta_mis1   <- "/path/to/metadata/MiSeq_run1_metadata.txt"
meta_mis2   <- "/path/to/metadata/MiSeq_run2_metadata.txt"
meta_novo   <- "/path/to/metadata/NovaSeq_metadata.txt"
meta_host   <- "/path/to/metadata/PNG_host_metadata.txt"   # PNG_2025-01-22.txt

# ASV phylogenetic tree (built on HPC — see comments in Section 11)
asv_tree_path <- "/path/to/ASV_cleaned_fasttree.nwk"

# Output directory (same folder as analysis scripts or a data/ folder)
output_dir  <- "/path/to/output"
output_rds  <- file.path(output_dir, "phyloseq_final.rds")


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(ShortRead)
library(Biostrings)
library(phyloseq)
library(vegan)
library(dplyr)
library(ape)


# ==============================================================================
# 2. CUSTOM FUNCTIONS (Kreisinger lab, adapted)
# ==============================================================================

# Retain only ASVs present in BOTH PCR duplicates for a given sample.
# Removes sequencing artefacts that appear in only one of two replicates.
dupl.concensus <- function(PHYLOSEQ, NAMES) {
  IDS      <- as.character(data.frame(sample_data(PHYLOSEQ))[, NAMES])
  IDS.dupl <- IDS[duplicated(IDS)]
  PHYLOSEQ <- prune_samples(IDS %in% IDS.dupl, PHYLOSEQ)
  CATS     <- as.character(data.frame(sample_data(PHYLOSEQ))[, NAMES])
  CATS2    <- levels(factor(CATS))
  OTU_TAB  <- otu_table(PHYLOSEQ)
  rownames(OTU_TAB) <- CATS
  for (i in seq_along(CATS2)) {
    FILTER.act <- colSums(OTU_TAB[rownames(OTU_TAB) == CATS2[i], ] > 0) > 1
    OTU_TAB[rownames(OTU_TAB) == CATS2[i], ] <-
      t(apply(OTU_TAB[rownames(OTU_TAB) == CATS2[i], ], 1, function(x) x * FILTER.act))
  }
  rownames(OTU_TAB) <- sample_names(PHYLOSEQ)
  otu_table(PHYLOSEQ) <- OTU_TAB
  prune_taxa(taxa_sums(PHYLOSEQ) > 0, PHYLOSEQ)
}

# Merge PCR duplicate pairs into a single sample (summed counts).
merge.duplicates <- function(PHYLOSEQ, NAMES) {
  CATS      <- as.character(data.frame(sample_data(PHYLOSEQ))[, NAMES])
  sample_data(PHYLOSEQ)$duplic.id <- CATS
  SAMDAT    <- sample_data(PHYLOSEQ)
  SAMDAT.sub <- SAMDAT[duplicated(CATS) == FALSE, ]
  FASTA     <- refseq(PHYLOSEQ)
  rownames(SAMDAT.sub) <- SAMDAT.sub$duplic.id
  PHYLOSEQ.merge <- merge_samples(PHYLOSEQ, "duplic.id")
  sample_data(PHYLOSEQ.merge) <- SAMDAT.sub
  merge_phyloseq(PHYLOSEQ.merge, FASTA)
}


# ==============================================================================
# 3. PER-RUN PHYLOSEQ CONSTRUCTION
# ==============================================================================
# Loads DADA2 outputs (OTU table, taxonomy, reference sequences) for each of
# the 5 sequencing runs (MiSeq 1–2, NovaSeq 1–3) and attaches per-run metadata.

# ---- MiSeq run 1 ----
load(file.path(dada2_root, "seqtab_18S_Mis1/otutab_PNG_18S_mis1.R"))
seqtab1 <- otu_table(seqtab, taxa_are_rows = FALSE)
load(file.path(dada2_root, "seqtab_18S_Mis1/TAX_PNG_18S_mis1.R"))
HAPLO1 <- readDNAStringSet(file.path(dada2_root, "seqtab_18S_Mis1/REF_PNG_18S_mis1.fasta"))
PNG_18S_mis1 <- merge_phyloseq(seqtab1, tax_table(taxa), HAPLO1)

MET <- read.delim(meta_mis1, header = TRUE, stringsAsFactors = FALSE)
MET$Adaptor_F <- spgs::reverseComplement(MET$Adaptor_F)
MET$Adaptor_F <- toupper(MET$Adaptor_F)
SN  <- paste(MET$Adaptor_R, "_", MET$Adaptor_F, MET$primerF, MET$primer_R, sep = "")
MET <- sample_data(MET); sample_names(MET) <- SN
PNG_18S_mis1 <- merge_phyloseq(PNG_18S_mis1, MET)
save(PNG_18S_mis1, file = file.path(dada2_root, "PNG_18S_mis1_MET.R"))
rm(seqtab, seqtab1, taxa, HAPLO1, MET, SN)

# ---- MiSeq run 2 ----
load(file.path(dada2_root, "seqtab_18S_Mis2/otutab_PNG_18S_mis2.R"))
seqtab1 <- otu_table(seqtab, taxa_are_rows = FALSE)
load(file.path(dada2_root, "seqtab_18S_Mis2/TAX_PNG_18S_mis2.R"))
HAPLO1 <- readDNAStringSet(file.path(dada2_root, "seqtab_18S_Mis2/REF_PNG_18S_mis2.fasta"))
PNG_18S_mis2 <- merge_phyloseq(seqtab1, tax_table(taxa), HAPLO1)

MET <- read.delim(meta_mis2, header = TRUE, stringsAsFactors = FALSE)
SN  <- paste(MET$Adaptor_R, "_", MET$Adaptor_F, MET$primerF, MET$primer_R, sep = "")
MET <- sample_data(MET); sample_names(MET) <- SN
PNG_18S_mis2 <- merge_phyloseq(PNG_18S_mis2, MET)
save(PNG_18S_mis2, file = file.path(dada2_root, "PNG_18S_mis2.R"))
rm(seqtab, seqtab1, taxa, HAPLO1, MET, SN)

# ---- NovaSeq runs 1–3 (same metadata file for all three) ----
for (run in 1:3) {
  run_name <- paste0("PNG_18S_novo", run)
  load(file.path(dada2_root, paste0("seqtab_18S_novo", run, "/otutab_", run_name, ".R")))
  seqtab1 <- otu_table(seqtab, taxa_are_rows = FALSE)
  load(file.path(dada2_root, paste0("seqtab_18S_novo", run, "/TAX_", run_name, ".R")))
  HAPLO1 <- readDNAStringSet(file.path(dada2_root,
               paste0("seqtab_18S_novo", run, "/REF_", run_name, ".fasta")))
  ps_obj  <- merge_phyloseq(seqtab1, tax_table(taxa), HAPLO1)
  MET <- read.delim(meta_novo, header = TRUE, stringsAsFactors = FALSE)
  SN  <- paste(MET$Novogene_fastq_name_full)
  MET <- sample_data(MET); sample_names(MET) <- SN
  ps_obj <- merge_phyloseq(ps_obj, MET)
  assign(run_name, ps_obj)
  save(list = run_name, file = file.path(dada2_root, paste0(run_name, ".R")))
  rm(seqtab, seqtab1, taxa, HAPLO1, MET, SN, ps_obj)
}


# ==============================================================================
# 4. MERGE PLATFORMS + FULL HOST METADATA
# ==============================================================================

load(file.path(dada2_root, "PNG_18S_mis1_MET.R"))
load(file.path(dada2_root, "PNG_18S_mis2.R"))
load(file.path(dada2_root, "PNG_18S_novo1.R"))
load(file.path(dada2_root, "PNG_18S_novo2.R"))
load(file.path(dada2_root, "PNG_18S_novo3.R"))

PNG_18S_Miseq12    <- merge_phyloseq(PNG_18S_mis1, PNG_18S_mis2)
PNG_18S_Novoseq123 <- merge_phyloseq(PNG_18S_novo1, PNG_18S_novo2, PNG_18S_novo3)
MET2 <- read.delim(meta_host)

dir.create(file.path(dada2_root, "METs"), showWarnings = FALSE)

# MiSeq: merge with full host metadata table (PNG_2025-01-22.txt)
MET11 <- as.data.frame(sample_data(PNG_18S_Miseq12)); class(MET11) <- "data.frame"
write.csv(MET11, file = file.path(dada2_root, "METs/MET_miseq12.csv"))
MET11 <- read.csv(file.path(dada2_root, "METs/MET_miseq12.csv"))
MET3  <- merge(MET11, MET2, by = "sample_ID", all.x = TRUE)
MET3  <- sample_data(MET3); sample_names(MET3) <- sample_names(PNG_18S_Miseq12)
sample_data(PNG_18S_Miseq12) <- MET3
PNG_18S_Miseq12_HT <- PNG_18S_Miseq12
save(PNG_18S_Miseq12_HT, file = file.path(dada2_root, "PNG_18S_Miseq12_HT.R"))

# NovaSeq: merge with full host metadata table
MET1 <- as.data.frame(sample_data(PNG_18S_Novoseq123)); class(MET1) <- "data.frame"
write.csv(MET1, file = file.path(dada2_root, "METs/MET_novoseq123.csv"))
MET1 <- read.csv(file.path(dada2_root, "METs/MET_novoseq123.csv"))
MET3 <- merge(MET1, MET2, by = "sample_ID", all.x = TRUE)
MET3 <- sample_data(MET3); sample_names(MET3) <- sample_names(PNG_18S_Novoseq123)
sample_data(PNG_18S_Novoseq123) <- MET3
PNG_18S_Novoseq123_HT <- PNG_18S_Novoseq123
save(PNG_18S_Novoseq123_HT, file = file.path(dada2_root, "PNG_18S_Novoseq123_HT.R"))

rm(list = ls(pattern = "^PNG_18S_mis|^PNG_18S_novo|MET"))


# ==============================================================================
# 5. PROCRUSTES VALIDATION OF PCR DUPLICATES
# ==============================================================================
# Confirms PCR duplicate concordance before merging.
# Expected: Procrustes correlation ~0.75, p < 0.001.

load(file.path(dada2_root, "PNG_18S_Miseq12_HT.R"))
PS_check <- prune_samples(sample_sums(PNG_18S_Miseq12_HT) > 500, PNG_18S_Miseq12_HT)
A <- prune_samples(regexpr("_A", sample_data(PS_check)$ID) > 0, PS_check)
B <- prune_samples(regexpr("_B", sample_data(PS_check)$ID) > 0, PS_check)
INTER <- intersect(sample_data(A)$ID_D, sample_data(B)$ID_D)
A <- prune_samples(sample_data(A)$ID_D %in% INTER, A)
B <- prune_samples(sample_data(B)$ID_D %in% INTER, B)
if (sum(duplicated(sample_data(B)$ID_D)) > 0) {
  dup_idx <- which(duplicated(sample_data(B)$ID_D))[1]
  B <- prune_samples(sample_names(B) != sample_names(B)[dup_idx], B)
}
A_t <- transform_sample_counts(A, function(x) x / sum(x))
B_t <- transform_sample_counts(B, function(x) x / sum(x))
DA  <- as.matrix(vegdist(otu_table(A_t))^0.5); colnames(DA) <- rownames(DA) <- sample_data(A)$ID_D
DB  <- as.matrix(vegdist(otu_table(B_t))^0.5); colnames(DB) <- rownames(DB) <- sample_data(B)$ID_D
DB  <- DB[rownames(DA), colnames(DA)]
print(protest(ape::pcoa(DA)$vectors, ape::pcoa(DB)$vectors))
rm(PS_check, A, B, INTER, A_t, B_t, DA, DB, dup_idx)


# ==============================================================================
# 6. DUPLICATE MERGING (per platform)
# ==============================================================================

load(file.path(dada2_root, "PNG_18S_mis1_MET.R"))
load(file.path(dada2_root, "PNG_18S_mis2.R"))
load(file.path(dada2_root, "PNG_18S_novo1.R"))
load(file.path(dada2_root, "PNG_18S_novo2.R"))
load(file.path(dada2_root, "PNG_18S_novo3.R"))

# MiSeq 1
SUMMARY <- summary(as.factor(sample_data(PNG_18S_mis1)$sample_ID), maxsum = 10000)
PNG_18S_mis1.dupl <- prune_samples(
  sample_data(PNG_18S_mis1)$sample_ID %in% names(SUMMARY)[SUMMARY == 2], PNG_18S_mis1)
PNG_18S_mis1.dupl <- merge.duplicates(
  dupl.concensus(PNG_18S_mis1.dupl, "sample_ID"), "sample_ID")
save(PNG_18S_mis1.dupl, file = file.path(dada2_root, "PNG_18S_mis1.dupl.R"))

# MiSeq 2
SUMMARY <- summary(as.factor(sample_data(PNG_18S_mis2)$sample_ID), maxsum = 10000)
PNG_18S_mis2.dupl <- prune_samples(
  sample_data(PNG_18S_mis2)$sample_ID %in% names(SUMMARY)[SUMMARY == 2], PNG_18S_mis2)
PNG_18S_mis2.dupl <- merge.duplicates(
  dupl.concensus(PNG_18S_mis2.dupl, "sample_ID"), "sample_ID")
save(PNG_18S_mis2.dupl, file = file.path(dada2_root, "PNG_18S_mis2.dupl.R"))

# NovaSeq 1–3 combined
PNG_18S_novo123 <- merge_phyloseq(PNG_18S_novo1, PNG_18S_novo2, PNG_18S_novo3)
REMOVE <- c("FR1200.005_CC_EtOH_C", "FR1200.005_CC_EtOH_D",
            "FR1200.039_CC_EtOH_C", "FR1200.040_CC_EtOH_C", "PNGlost9_CC_EtOH_C")
PNG_18S_novo123 <- prune_samples(!sample_data(PNG_18S_novo123)$ID %in% REMOVE, PNG_18S_novo123)
SUMMARY <- summary(as.factor(sample_data(PNG_18S_novo123)$sample_ID), maxsum = 10000)
PNG_18S_novo123.dupl <- prune_samples(
  sample_data(PNG_18S_novo123)$sample_ID %in% names(SUMMARY)[SUMMARY == 2], PNG_18S_novo123)
PNG_18S_novo123.dupl <- merge.duplicates(
  dupl.concensus(PNG_18S_novo123.dupl, "sample_ID"), "sample_ID")
save(PNG_18S_novo123.dupl, file = file.path(dada2_root, "PNG_18S_novo123.dupl.R"))

rm(list = ls(pattern = "^PNG_18S_mis|^PNG_18S_novo|SUMMARY|REMOVE"))


# ==============================================================================
# 7. HOST METADATA RE-ATTACHMENT + MERGE ALL PLATFORMS
# ==============================================================================

load(file.path(dada2_root, "PNG_18S_mis1.dupl.R"))
load(file.path(dada2_root, "PNG_18S_mis2.dupl.R"))
load(file.path(dada2_root, "PNG_18S_novo123.dupl.R"))
MET2 <- read.delim(meta_host)

for (ps_name in c("PNG_18S_mis1.dupl", "PNG_18S_mis2.dupl", "PNG_18S_novo123.dupl")) {
  ps_obj <- get(ps_name)
  MET11  <- as.data.frame(sample_data(ps_obj)); class(MET11) <- "data.frame"
  MET3   <- merge(MET11, MET2, by = "sample_ID", all.x = TRUE)
  MET3   <- sample_data(MET3); sample_names(MET3) <- sample_names(ps_obj)
  sample_data(ps_obj) <- MET3
  new_name <- paste0(ps_name, "_HT")
  assign(new_name, ps_obj)
  save(list = new_name, file = file.path(dada2_root, paste0(new_name, ".R")))
}
rm(ps_obj, MET11, MET3, ps_name, new_name, MET2)

sample_data(PNG_18S_mis1.dupl_HT)$seqrun   <- "mis"
sample_data(PNG_18S_mis2.dupl_HT)$seqrun   <- "mis"
sample_data(PNG_18S_novo123.dupl_HT)$seqrun <- "novo"

PNG_18S_mis12_novo123.dupl.HT <- merge_phyloseq(
  PNG_18S_mis1.dupl_HT, PNG_18S_mis2.dupl_HT, PNG_18S_novo123.dupl_HT
)
PNG_18S_mis12_novo123.dupl.HT <- prune_samples(
  sample_sums(PNG_18S_mis12_novo123.dupl.HT) > 500, PNG_18S_mis12_novo123.dupl.HT)
PNG_18S_mis12_novo123.dupl.HT <- prune_taxa(
  taxa_sums(PNG_18S_mis12_novo123.dupl.HT) > 1, PNG_18S_mis12_novo123.dupl.HT)
save(PNG_18S_mis12_novo123.dupl.HT,
     file = file.path(dada2_root, "PNG_18S_mis12_novo123.dupl.HT.R"))


# ==============================================================================
# 8. REMOVE CONTROLS / NON-PNG SAMPLES / SEPARATE RODENTS
# ==============================================================================

rm(list = ls(pattern = "^PNG_18S_mis|^PNG_18S_novo"))
load(file.path(dada2_root, "PNG_18S_mis12_novo123.dupl.HT.R"))
PHYSEQ <- PNG_18S_mis12_novo123.dupl.HT

PHYSEQ <- prune_samples(sample_data(PHYSEQ)$extracted_tissue != "control_isolation", PHYSEQ)
PHYSEQ <- prune_samples(sample_data(PHYSEQ)$extracted_tissue != "Control_negative",  PHYSEQ)
PHYSEQ <- prune_samples(sample_data(PHYSEQ)$sample_ID != "MW3200.138", PHYSEQ)
PHYSEQ <- prune_samples(sample_data(PHYSEQ)$Dataset_2 == "PNG", PHYSEQ)
PHYSEQ <- prune_taxa(taxa_sums(PHYSEQ) > 0, PHYSEQ)

PNG_18S_rodent.dupl <- prune_samples(
  sample_data(PHYSEQ)$Order_Infraclass == "Rodentia", PHYSEQ)
PNG_18S_rodent.dupl <- prune_taxa(taxa_sums(PNG_18S_rodent.dupl) > 0, PNG_18S_rodent.dupl)
save(PNG_18S_rodent.dupl, file = file.path(dada2_root, "PNG_18S_rodent.dupl.R"))


# ==============================================================================
# 9. TAXONOMY FILTERING + DIET COLUMN
# ==============================================================================

rm(PHYSEQ, PNG_18S_mis12_novo123.dupl.HT)
load(file.path(dada2_root, "PNG_18S_rodent.dupl.R"))
PHYSEQ <- PNG_18S_rodent.dupl

PHYSEQ <- subset_taxa(PHYSEQ, !is.na(Phylum) & Phylum != "")

phyla_to_remove <- c(
  "Acanthocephala", "Apicomplexa", "Bacillariophyta", "Blastocladiomycota",
  "Cercozoa", "Ciliophora", "Chytridiomycota", "Endomyxa", "Evosea",
  "Mucoromycota", "Nematoda", "Oomycota", "Parabasalia", "Preaxostyla",
  "Tubulinea", "Zoopagomycota", "Olpidiomycota", "Sanchytriomycota",
  "Platyhelminthes", "Cnidaria", "Ctenophora", "Discosea", "Echinodermata",
  "Gastrotricha", "Nematomorpha", "Nemertea", "Porifera", "Prasinodermophyta",
  "Rhodophyta", "Rotifera", "Tardigrada", "Xenacoelomorpha", "Euglenozoa"
)
PHYSEQ <- subset_taxa(PHYSEQ, !(Phylum %in% phyla_to_remove))

PHYSEQ <- subset_taxa(PHYSEQ,
  !(Phylum == "Chordata" & (
    is.na(Class) | Class == "Ascidiacea" | is.na(Order) |
    Family == "Hominidae" | !(Class %in% c("Amphibia", "Lepidosauria"))
  ))
)

exclude_class <- c(
  "Atractiellomycetes", "Cystobasidiomycetes", "Exobasidiomycetes", "Microbotryomycetes",
  "Moniliellomycetes", "Pucciniomycetes", "Tritirachiomycetes", "Ustilaginomycetes",
  "Wallemiomycetes", "Saccharomycetes", "Eurotiomycetes", "Trematoda",
  "Agaricostilbomycetes", "Chlorophyceae", "Dipodascomycetes", "Dothideomycetes",
  "Hexanauplia", "Klebsormidiophyceae", "Malasseziomycetes", "Pichiomycetes",
  "Schizosaccharomycetes", "Sordariomycetes", "Spiculogloeomycetes",
  "Taphrinomycetes", "Thecostraca", "Trebouxiophyceae", "Tremellomycetes", "Ulvophyceae"
)
PHYSEQ <- prune_taxa(!(tax_table(PHYSEQ)[, "Class"] %in% exclude_class), PHYSEQ)

exclude_order <- c("Ixodida", "Mesostigmata", "Pseudoscorpiones", "Sarcoptiformes",
                   "Trombidiformes", "Enchytraeida", "Hirudinida", "Monhysterida",
                   "Plectida", "Rhabditida")
PHYSEQ <- prune_taxa(!(tax_table(PHYSEQ)[, "Order"] %in% exclude_order), PHYSEQ)
PHYSEQ <- prune_taxa(!(tax_table(PHYSEQ)[, "Family"] %in% c("Aeolosomatidae")), PHYSEQ)

tax_df <- as.data.frame(tax_table(PHYSEQ))
tax_df$Diet <- dplyr::case_when(
  tax_df$Phylum %in% c("Arthropoda", "Chordata", "Mollusca", "Annelida") ~ "Animal",
  tax_df$Phylum %in% c("Ascomycota", "Basidiomycota")                   ~ "Fungi",
  tax_df$Phylum %in% c("Streptophyta", "Chlorophyta")                   ~ "Plant",
  TRUE ~ "NA"
)
tax_df <- dplyr::select(tax_df, 1, Diet, dplyr::everything())
tax_table(PHYSEQ) <- as.matrix(tax_df)

PNG_18S_rodent_filtered <- PHYSEQ
save(PNG_18S_rodent_filtered, file = file.path(dada2_root, "PNG_18S_rodent_filtered.R"))


# ==============================================================================
# 10. ECOLOGICAL ENSEMBLE ANNOTATION
# ==============================================================================

PHYSEQ <- PNG_18S_rodent_filtered
sample_data(PHYSEQ)$Ecological_ensembles <- gsub(
  "-", "_", gsub("\\.", "_", sample_data(PHYSEQ)$Ecological_ensembles))

split_data <- strsplit(sample_data(PHYSEQ)$Ecological_ensembles, "_")
sample_data(PHYSEQ)$size      <- sapply(split_data, `[`, 1)
sample_data(PHYSEQ)$lifestyle <- sapply(split_data, `[`, 2)
sample_data(PHYSEQ)$diet      <- sapply(split_data, `[`, 3)

species_lifestyle <- c(
  "Anisomys imitatorFR" = "scansorial", "Anisomys imitatorMW" = "scansorial",
  "Hyomys goliath" = "terrestrial",     "Mallomys cf_aroaensis" = "terrestrial",
  "Mallomys hercules" = "terrestrial",  "Mallomys istapantap" = "terrestrial",
  "Mallomys rothschildi" = "scansorial","Uromys anak" = "scansorial"
)
sam_data <- sample_data(PHYSEQ)
sam_data$lifestyle <- ifelse(
  sam_data$Species_consensus_ddRAD_mt_others %in% names(species_lifestyle),
  species_lifestyle[sam_data$Species_consensus_ddRAD_mt_others], sam_data$lifestyle)
sample_data(PHYSEQ) <- sam_data

sample_data(PHYSEQ)$new_ensembles <- with(sample_data(PHYSEQ),
  paste(size, lifestyle, diet, sep = "_"))
sample_data(PHYSEQ)$new_ensembles <- ifelse(
  sample_data(PHYSEQ)$new_ensembles == "na_terrestrial_omnivore",
  paste(sample_data(PHYSEQ)$new_ensembles,
        tolower(sample_data(PHYSEQ)$Tribe), sep = "_"),
  sample_data(PHYSEQ)$new_ensembles)

PNG_18S_rodent_filtered_ensembles <- PHYSEQ
save(PNG_18S_rodent_filtered_ensembles,
     file = file.path(dada2_root, "PNG_18S_rodent_filtered_ensembles.R"))


# ==============================================================================
# 11. ASV PHYLOGENETIC TREE ATTACHMENT
# ==============================================================================
# The 18S ASV tree was constructed on an HPC cluster:
#   1. Export: writeFasta(refseq(PHYSEQ), "ASV.fasta")
#   2. Align with mothur against SILVA 132 seed alignment:
#      align.seqs(candidate=ASV.fasta, template=silva.seed_v132.align, processors=10)
#   3. Clean alignment:
#      - remove columns with 100% gaps
#      - remove ASVs with >60% gaps after column filtering
#   4. Build tree: FastTree -nt -gtr ASV_cleaned.align > ASV_cleaned_fasttree.nwk
#
# See bioinformatics_scripts/ for the full shell pipeline.
# Set asv_tree_path at the top of this script before running.

FASTTREE <- ape::read.tree(asv_tree_path)
PHYSEQ   <- PNG_18S_rodent_filtered_ensembles

shared_taxa <- intersect(taxa_names(PHYSEQ), FASTTREE$tip.label)
PHYSEQ   <- prune_taxa(shared_taxa, PHYSEQ)
FASTTREE <- ape::keep.tip(FASTTREE, shared_taxa)

PHYSEQ <- merge_phyloseq(
  sample_data(PHYSEQ), otu_table(PHYSEQ), tax_table(PHYSEQ), refseq(PHYSEQ), FASTTREE
)
PNG_18S_rodent_filtered_ensembles_tree <- PHYSEQ
save(PNG_18S_rodent_filtered_ensembles_tree,
     file = file.path(dada2_root, "PNG_18S_rodent_filtered_ensembles_tree.R"))


# ==============================================================================
# 12. FINAL CLEANUP
# ==============================================================================

PHYSEQ <- PNG_18S_rodent_filtered_ensembles_tree

# Standardise mountain range label
sample_data(PHYSEQ)$Metalokalita <- ifelse(
  sample_data(PHYSEQ)$Metalokalita %in% c("Baitabag", "Nagada", "Wanang"),
  "Mt_Wilhelm", sample_data(PHYSEQ)$Metalokalita)
sample_data(PHYSEQ)$Metalokalita <- factor(sample_data(PHYSEQ)$Metalokalita)

# Standardise locality codes
sample_data(PHYSEQ)$Loc[sample_data(PHYSEQ)$Loc == "WG200"] <- "MW150"
sample_data(PHYSEQ)$Loc[sample_data(PHYSEQ)$Loc == "BAI"]   <- "MW45"

# Remove Rattus norvegicus (single captive-origin individual)
PHYSEQ <- prune_samples(
  !sample_data(PHYSEQ)$Species_consensus_ddRAD_mt_others == "Rattus norvegicus", PHYSEQ)
PHYSEQ <- prune_taxa(taxa_sums(PHYSEQ) > 0, PHYSEQ)

# Normalise habitat categories
hab <- sample_data(PHYSEQ)$Habitat_polished_Daniel.Fanda
hab[hab %in% c("secondary", "primary", "shrub/primary")] <- "forest"
hab[hab %in% c("house", "garden")]                        <- "synanthropic"
sample_data(PHYSEQ)$Habitat_polished_Daniel.Fanda <- hab

# Rename lifestyle → spatial_niche
sdat <- data.frame(sample_data(PHYSEQ))
sdat <- dplyr::rename(sdat, spatial_niche = lifestyle)
sample_data(PHYSEQ) <- sample_data(sdat)

# Fix 4 samples with NA habitat (confirmed post-fieldwork)
samples_to_fix <- c("MW200.022", "BAI.032", "BAI.036", "MW2200.070")
new_values     <- c("synanthropic", "forest", "forest", "synanthropic")
for (i in seq_along(samples_to_fix))
  sample_data(PHYSEQ)[samples_to_fix[i], "Habitat_polished_Daniel.Fanda"] <- new_values[i]

cat("Final phyloseq: ", ntaxa(PHYSEQ), "taxa,", nsamples(PHYSEQ), "samples\n")
print(PHYSEQ)


# ==============================================================================
# 13. SAVE
# ==============================================================================

saveRDS(PHYSEQ, file = output_rds)
cat("Saved:", output_rds, "\n")


# ==============================================================================
# 14. SUPPLEMENTARY: EXCLUDED TAXA TABLE
# ==============================================================================

load(file.path(dada2_root, "PNG_18S_rodent.dupl.R"))
PNG_initial <- prune_samples(sample_sums(PNG_18S_rodent.dupl) >= 500, PNG_18S_rodent.dupl)

excluded_taxa_table <- as.data.frame(tax_table(PNG_initial)) %>%
  tibble::rownames_to_column("ASV_ID") %>%
  dplyr::filter(ASV_ID %in% setdiff(taxa_names(PNG_initial), taxa_names(PHYSEQ))) %>%
  dplyr::left_join(
    data.frame(Raw_Abundance = taxa_sums(PNG_initial)) %>% tibble::rownames_to_column("ASV_ID"),
    by = "ASV_ID") %>%
  dplyr::select(ASV_ID, Kingdom, Phylum, Class, Order, Raw_Abundance) %>%
  dplyr::arrange(dplyr::desc(Raw_Abundance))

write.csv(excluded_taxa_table,
          file.path(output_dir, "supplementary_excluded_taxa.csv"),
          row.names = FALSE)
