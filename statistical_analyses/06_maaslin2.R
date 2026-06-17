################################################################################
# 06_maaslin2.R
# MaAsLin2 differential abundance at dietary order level
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# Tests which dietary orders are differentially abundant across host predictors
# (spatial niche, altitude, tribe, body mass) using MaAsLin2 on Mt. Wilhelm
# and Finisterres samples separately. Excludes invasive Rattus (R. rattus group
# and incertae sedis group) to focus on native rodent community.
################################################################################


# ==============================================================================
# 0. PATHS
# ==============================================================================

data_dir     <- "/path/to/project/data"
phyloseq_dir <- file.path(data_dir, "phyloseqs")
output_dir   <- "/path/to/output/maaslin2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(phyloseq)
library(dplyr)
library(Maaslin2)


# ==============================================================================
# 2. CUSTOM FUNCTIONS
# ==============================================================================

# Propagate parent-rank name into unassigned child ranks
manage_unassigned <- function(PHYLOSEQ, UNASS_STRING = NA, ADD, AFTER = TRUE) {
  TAXall <- data.frame(tax_table(PHYLOSEQ), stringsAsFactors = FALSE)
  if (!is.na(UNASS_STRING)) TAXall[TAXall == UNASS_STRING] <- NA
  for (i in seq_len(nrow(TAXall))) {
    TAX <- as.character(TAXall[i, ])
    if (AFTER) {
      for (j in 2:length(TAX))
        if (is.na(TAX[j]))
          TAX[j] <- ifelse(grepl(ADD, TAX[j - 1], fixed = TRUE),
                           TAX[j - 1], paste0(TAX[j - 1], ADD))
    }
    TAXall[i, ] <- TAX
  }
  colnames_orig <- colnames(TAXall)
  TAXall <- tax_table(as.matrix(TAXall))
  colnames(TAXall) <- colnames_orig
  tax_table(PHYLOSEQ) <- TAXall
  PHYLOSEQ
}

# Rename duplicate Order names with unique suffixes
rename_taxa_uniquely <- function(physeq, taxrank) {
  taxa <- data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  taxa$BaseName <- taxa[[taxrank]]
  taxa$BaseName[is.na(taxa$BaseName) | taxa$BaseName == ""] <-
    paste0("Unassigned_", seq_len(sum(is.na(taxa$BaseName) | taxa$BaseName == "")))
  taxa <- taxa %>%
    dplyr::group_by(BaseName) %>%
    dplyr::mutate(
      Count      = dplyr::row_number(),
      UniqueName = dplyr::if_else(dplyr::n() > 1,
                                  paste(BaseName, Count, sep = "_"), BaseName)
    ) %>%
    dplyr::ungroup()
  taxa_names(physeq) <- taxa$UniqueName
  physeq
}


# ==============================================================================
# 3. LOAD + FILTER PHYLOSEQ
# ==============================================================================

ps <- readRDS(file.path(phyloseq_dir, "phyloseq_final.rds"))

# Exclude invasive Rattus (not native to PNG)
ps_native <- subset_samples(
  ps,
  !Hydromyini.Division.Rattini.Group %in% c("R. rattus group", "incertae sedis group")
)
ps_native <- filter_taxa(ps_native, function(x) sum(x) > 0, prune = TRUE)

# Split by mountain transect
ps_wilhelm    <- subset_samples(ps_native, Metalokalita == "Mt_Wilhelm")
ps_finisterres <- subset_samples(ps_native, Metalokalita == "Finisterres")
ps_wilhelm    <- filter_taxa(ps_wilhelm,    function(x) sum(x) > 0, prune = TRUE)
ps_finisterres <- filter_taxa(ps_finisterres, function(x) sum(x) > 0, prune = TRUE)


# ==============================================================================
# 4. PREPARE ORDER-LEVEL PHYLOSEQ (Mt. Wilhelm)
# ==============================================================================

prepare_order_ps <- function(ps_obj) {
  ps_order <- tax_glom(ps_obj, taxrank = "Order", NArm = FALSE)
  ps_order <- manage_unassigned(ps_order, UNASS_STRING = "NA", ADD = "", AFTER = TRUE)

  # Fix remaining empty Order entries
  ord <- tax_table(ps_order)[, "Order"]
  empty_idx <- which(ord == "")
  for (k in seq_along(empty_idx))
    tax_table(ps_order)[empty_idx[k], "Order"] <- paste0("Unknown_Order_", k)

  rename_taxa_uniquely(ps_order, taxrank = "Order")
}

ps_w_order <- prepare_order_ps(ps_wilhelm)
ps_f_order <- prepare_order_ps(ps_finisterres)


# ==============================================================================
# 5. FILTER: KEEP ORDERS PRESENT IN > 50 SAMPLES
# ==============================================================================

filter_by_prevalence <- function(ps_obj, min_samples = 50) {
  ps_bin  <- transform_sample_counts(ps_obj, function(x) ifelse(x == 0, 0, 1))
  keep    <- taxa_sums(ps_bin) > min_samples
  ps_filt <- prune_taxa(keep, ps_obj)
  prune_samples(sample_sums(ps_filt) > 0, ps_filt)
}

ps_w_filt <- filter_by_prevalence(ps_w_order, min_samples = 50)
ps_f_filt <- filter_by_prevalence(ps_f_order, min_samples = 50)

cat("Mt. Wilhelm — orders retained:", ntaxa(ps_w_filt),
    "| samples:", nsamples(ps_w_filt), "\n")
cat("Finisterres — orders retained:", ntaxa(ps_f_filt),
    "| samples:", nsamples(ps_f_filt), "\n")


# ==============================================================================
# 6. RUN MAASLIN2 — MT. WILHELM
# ==============================================================================

run_maaslin2 <- function(ps_obj, fixed_effects, out_dir) {
  OTU_tab <- as.data.frame(t(otu_table(ps_obj)))
  SD_df   <- as.data.frame(sample_data(ps_obj))
  class(SD_df) <- "data.frame"

  # Keep only samples with no NA in fixed effects
  complete_idx <- complete.cases(SD_df[, fixed_effects])
  OTU_tab <- OTU_tab[complete_idx, ]
  SD_df   <- SD_df[complete_idx, ]

  Maaslin2(
    input_data     = OTU_tab,
    input_metadata = SD_df,
    output         = out_dir,
    fixed_effects  = fixed_effects,
    random_effects = c("Species_consensus_ddRAD_mt_others"),
    normalization  = "NONE",
    transform      = "LOG",
    analysis_method = "LM",
    max_significance = 0.25,
    plot_heatmap   = TRUE,
    plot_scatter   = FALSE
  )
}

fixed_fx <- c("Alt", "spatial_niche", "Tribe", "mean_mass",
              "Habitat_polished_Daniel.Fanda")

cat("Running MaAsLin2 — Mt. Wilhelm...\n")
fit_wilhelm <- run_maaslin2(
  ps_w_filt, fixed_fx,
  file.path(output_dir, "MaAsLin2_Wilhelm")
)

cat("Running MaAsLin2 — Finisterres...\n")
fit_finisterres <- run_maaslin2(
  ps_f_filt, fixed_fx,
  file.path(output_dir, "MaAsLin2_Finisterres")
)


# ==============================================================================
# 7. COMBINE + FILTER SIGNIFICANT RESULTS
# ==============================================================================

read_maaslin_results <- function(out_dir, mountain) {
  res <- read.csv(file.path(out_dir, "all_results.tsv"), sep = "\t")
  res$mountain <- mountain
  res
}

results_combined <- dplyr::bind_rows(
  read_maaslin_results(file.path(output_dir, "MaAsLin2_Wilhelm"),    "Mt_Wilhelm"),
  read_maaslin_results(file.path(output_dir, "MaAsLin2_Finisterres"), "Finisterres")
)

sig_results <- results_combined %>%
  filter(qval < 0.25) %>%
  arrange(qval)

write.csv(results_combined, file.path(output_dir, "MaAsLin2_all_results.csv"),  row.names = FALSE)
write.csv(sig_results,      file.path(output_dir, "MaAsLin2_significant.csv"), row.names = FALSE)

cat("MaAsLin2 complete. Significant associations (q < 0.25):", nrow(sig_results), "\n")
cat("Results saved to:", output_dir, "\n")
