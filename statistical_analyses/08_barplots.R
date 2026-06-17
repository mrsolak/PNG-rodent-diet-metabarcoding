################################################################################
# 02_barplots.R
# Dietary composition barplots per ecological ensemble
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# Loads phyloseq_final.rds and produces stacked barplots showing the top 25
# dietary orders per ecological ensemble (Figures in manuscript).
################################################################################


# ==============================================================================
# 0. PATHS
# ==============================================================================

data_dir     <- "/path/to/project/data"
phyloseq_dir <- file.path(data_dir, "phyloseqs")
output_dir   <- "/path/to/output/barplots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(phyloseq)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(RColorBrewer)
library(viridis)


# ==============================================================================
# 2. LOAD PHYLOSEQ
# ==============================================================================

PHYSEQ <- readRDS(file.path(phyloseq_dir, "phyloseq_final.rds"))


# ==============================================================================
# 3. HELPER FUNCTION: fill missing taxonomic ranks propagated from parent
# ==============================================================================

manage_unassigned <- function(PHYLOSEQ, UNASS_STRING = NA, ADD, AFTER = TRUE) {
  TAXall <- data.frame(tax_table(PHYLOSEQ), stringsAsFactors = FALSE)
  if (!is.na(UNASS_STRING)) TAXall[TAXall == UNASS_STRING] <- NA
  for (i in 1:nrow(TAXall)) {
    TAX <- as.character(TAXall[i, ])
    if (AFTER) {
      for (j in 2:length(TAX))
        if (is.na(TAX[j]))
          TAX[j] <- ifelse(grepl(ADD, TAX[j - 1], fixed = TRUE), TAX[j - 1],
                           paste0(TAX[j - 1], ADD))
    }
    TAXall[i, ] <- TAX
  }
  colnames_orig <- colnames(TAXall)
  TAXall <- tax_table(as.matrix(TAXall))
  colnames(TAXall) <- colnames_orig
  tax_table(PHYLOSEQ) <- TAXall
  PHYLOSEQ
}


# ==============================================================================
# 4. PREPARE PHYLOSEQ
# ==============================================================================

# Fill unassigned orders with parent-rank suffix
PHYSEQ <- manage_unassigned(PHYSEQ, UNASS_STRING = "NA", ADD = "", AFTER = TRUE)

# Fix specific empty Order entries
ord <- tax_table(PHYSEQ)[, "Order"]
ord[ord == ""][1] <- "Class:Agaricomycetes1"
ord[ord == ""][2] <- "Class:Agaricomycetes2"
ord[ord == ""][3] <- "Class:Agaricomycetes3"
ord[ord == ""][4] <- "Class:Gastropoda1"
ord[ord == ""][5] <- "Class:Symphyla1"
tax_table(PHYSEQ)[, "Order"] <- ord

# Aggregate to Order level
PHYSEQ_order <- tax_glom(PHYSEQ, taxrank = "Order")

# Colour palette (27 colours for up to 25 orders + Others + NA)
base_colors <- c(
  "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "dodgerblue2",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gold1", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown", "purple", "darkred", "pink"
)

# Define ensembles
ensembles <- list(
  TER_OMN_RAT = subset_samples(PHYSEQ_order, new_ensembles == "na_terrestrial_omnivore_rattini"),
  TER_OMN_HYD = subset_samples(PHYSEQ_order, new_ensembles == "na_terrestrial_omnivore_hydromyini"),
  SMA_SCA_HER = subset_samples(PHYSEQ_order, new_ensembles == "small_scansorial_herbivore"),
  TER_CAR     = subset_samples(PHYSEQ_order, new_ensembles == "na_terrestrial_carnivore"),
  SCA_OMN     = subset_samples(PHYSEQ_order, new_ensembles == "na_scansorial_omnivore"),
  LAR_TER_HER = subset_samples(PHYSEQ_order, new_ensembles == "large_terrestrial_herbivore"),
  LAR_SCA_HER = subset_samples(PHYSEQ_order, new_ensembles == "large_scansorial_herbivore")
)


# ==============================================================================
# 5. PLOTTING FUNCTION
# ==============================================================================

plot_ensemble <- function(ps_obj, ensemble_name) {
  ps_obj <- prune_taxa(taxa_sums(ps_obj) > 0, ps_obj)
  ps_obj <- prune_samples(sample_sums(ps_obj) > 0, ps_obj)

  tax_tab <- as.data.frame(tax_table(ps_obj))
  otu_mat <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) otu_mat <- t(otu_mat)
  tax_tab$total_abundance <- rowSums(otu_mat)

  order_abundance <- tax_tab %>%
    group_by(Order) %>%
    summarize(order_total = sum(total_abundance), n_taxa = n()) %>%
    arrange(desc(order_total))

  n_top    <- min(25, nrow(order_abundance))
  top_orders  <- order_abundance$Order[1:n_top]
  to_merge <- setdiff(order_abundance$Order, top_orders)

  if (length(to_merge) > 0) {
    taxa_keep  <- rownames(tax_tab)[tax_tab$Order %in% top_orders]
    taxa_merge <- rownames(tax_tab)[tax_tab$Order %in% to_merge]
    others_row <- colSums(otu_mat[taxa_merge, , drop = FALSE])
    otu_mat_f  <- rbind(otu_mat[taxa_keep, , drop = FALSE], Others = others_row)
    tax_f      <- as.matrix(tax_table(ps_obj)[taxa_keep, , drop = FALSE])
    others_tax <- tax_f[1, , drop = FALSE]
    rownames(others_tax) <- "Others"
    others_tax[,] <- "Others"
    tax_f      <- rbind(tax_f, others_tax)
    ps_obj <- phyloseq(otu_table(otu_mat_f, taxa_are_rows = TRUE),
                       tax_table(tax_f), sample_data(ps_obj))
  }

  ps_obj <- transform_sample_counts(ps_obj, function(x) x / sum(x))

  sdf <- data.frame(
    sample_ID = sample_names(ps_obj),
    species   = as.character(sample_data(ps_obj)$Species_consensus_ddRAD_mt_others),
    mountain  = as.character(sample_data(ps_obj)$Metalokalita),
    altitude  = as.numeric(sample_data(ps_obj)$Alt),
    stringsAsFactors = FALSE
  )
  sdf <- dplyr::arrange(sdf, species, mountain, altitude)
  sample_data(ps_obj)$sample_ID <- factor(sample_names(ps_obj), levels = sdf$sample_ID)

  psmelt(ps_obj) %>%
    ggplot(aes(x = sample_ID, y = Abundance, fill = Order)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = setNames(base_colors[seq_along(unique(Order))], unique(Order))) +
    labs(title = ensemble_name, x = "Sample", y = "Relative abundance") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "right", legend.text = element_text(size = 7))
}


# ==============================================================================
# 6. GENERATE AND SAVE PLOTS
# ==============================================================================

for (ens_name in names(ensembles)) {
  p <- plot_ensemble(ensembles[[ens_name]], ens_name)
  ggsave(file.path(output_dir, paste0("barplot_", ens_name, ".pdf")),
         plot = p, width = 14, height = 5)
  cat("Saved:", ens_name, "\n")
}
