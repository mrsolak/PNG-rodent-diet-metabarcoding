################################################################################
# 03_niche_partitioning.R
# Individual dietary niche partitioning — indicspecies WIC/TNW/BIC/SI
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# Computes species-level niche breadth (TNW), within-individual consistency
# (WIC), and between-individual variation (BIC/Spec_index) using the
# indicspecies::nichevar() function on rarefied, phylogenetically-weighted data.
# Results are saved as CSV files for use by 04_asr_glmm.R.
# Also produces the PCoA of phylogenetic niche overlap (manuscript Figure 4).
#
# RUNTIME: parallelised over 8 cores; can take several hours for 44 species.
################################################################################


# ==============================================================================
# 0. PATHS
# ==============================================================================

data_dir       <- "/path/to/project/data"
phyloseq_dir   <- file.path(data_dir, "phyloseqs")
host_tree_path <- file.path(data_dir, "ddRAD_trees", "Joined_scaled_RH85.nwk")
output_dir     <- "/path/to/output/niche_partitioning"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(indicspecies)
library(phyloseq)
library(ape)
library(ShortRead)
library(foreach)
library(doSNOW)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)


# ==============================================================================
# 2. LOAD + PREPARE PHYLOSEQ
# ==============================================================================

ps <- readRDS(file.path(phyloseq_dir, "phyloseq_final.rds"))

# Rename ASVs to short IDs (required by indicspecies)
taxa_names(ps) <- paste0("ASC_", seq_len(ntaxa(ps)))

# Replace ASV tree with cleaned FastTree (built from aligned ASVs)
FASTTREE <- ape::read.tree(file.path(data_dir, "ASV_cleaned_fasttree.nwk"))  # ASV tree (HPC-generated FastTree)
FASTTREE$edge.length <- sqrt(FASTTREE$edge.length)   # sqrt transform suppresses long-branch effect
ps <- merge_phyloseq(sample_data(ps), otu_table(ps), tax_table(ps), refseq(ps), FASTTREE)

# Filter to species with >= 6 individuals, rarefy, transform to proportions
sample_data(ps)$Sp <- sample_data(ps)$Species_consensus_ddRAD_mt_others
SUMM   <- summary(as.factor(sample_data(ps)$Sp))
SUMM   <- SUMM[SUMM > 5]
ps_filt <- prune_samples(sample_data(ps)$Sp %in% names(SUMM), ps)
ps_rare <- rarefy_even_depth(ps_filt, rngseed = 42)
ps_prop <- transform_sample_counts(ps_rare, function(x) x / sum(x))

cat("Species retained (n >= 6):", length(SUMM), "\n")
cat("Samples after rarefaction:", nsamples(ps_rare), "\n")

# Cophenetic phylogenetic distance matrix (sqrt-transformed, rescaled 0–1)
PD <- cophenetic(phy_tree(ps_prop))
PD <- PD^0.5
PD <- PD / max(PD)
PD <- as.dist(PD)

Sp <- unique(sample_data(ps_prop)$Sp)


# ==============================================================================
# 3. SPECIES-LEVEL NICHE BREADTH — PHYLOGENETIC (main analysis)
# ==============================================================================

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoSNOW(cl)
pb   <- txtProgressBar(max = length(Sp), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

NV.phylo <- foreach(i = seq_along(Sp),
                    .packages = c("indicspecies", "ape", "phyloseq"),
                    .export   = c("ps_prop", "Sp", "PD"),
                    .options.snow = opts,
                    .combine  = "rbind") %dopar% {
  PHS.act <- prune_samples(sample_data(ps_prop)$Sp == Sp[i], ps_prop)
  DIET    <- data.frame(otu_table(PHS.act)); class(DIET) <- "data.frame"
  nichevar(P = DIET, mode = "single", D = PD, nboot = 100)
}
close(pb); stopCluster(cl)

NV.phylo.df <- data.frame(NV.phylo, Sp)
write.csv(NV.phylo.df,
          file.path(output_dir, "species_niche_breadth_phylo.csv"),
          row.names = FALSE)


# ==============================================================================
# 4. INDIVIDUAL-LEVEL NICHE METRICS — PHYLOGENETIC
# ==============================================================================

cl <- makeCluster(num_cores)
registerDoSNOW(cl)
pb   <- txtProgressBar(max = length(Sp), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

IN.phylo <- foreach(i = seq_along(Sp),
                    .packages = c("indicspecies", "ape", "phyloseq"),
                    .export   = c("ps_prop", "Sp", "PD"),
                    .options.snow = opts) %dopar% {
  PHS.act <- prune_samples(sample_data(ps_prop)$Sp == Sp[i], ps_prop)
  DIET    <- data.frame(otu_table(PHS.act)); class(DIET) <- "data.frame"
  nichevar(P = DIET, mode = "ind", D = PD, nboot = 100)$IS
}
close(pb); stopCluster(cl)

names(IN.phylo) <- Sp

# ==============================================================================
# 5. PAIRWISE OVERLAP — PHYLOGENETIC (Inverse BIC)
# ==============================================================================

cl <- makeCluster(num_cores)
registerDoSNOW(cl)
pb   <- txtProgressBar(max = length(Sp), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

OVind.phylo <- foreach(i = seq_along(Sp),
                       .packages = c("indicspecies", "ape", "phyloseq"),
                       .export   = c("ps_prop", "Sp", "PD"),
                       .options.snow = opts) %dopar% {
  PHS.act <- prune_samples(sample_data(ps_prop)$Sp == Sp[i], ps_prop)
  DIET    <- data.frame(otu_table(PHS.act)); class(DIET) <- "data.frame"
  nicheoverlap(P = DIET, mode = "single", D = PD)
}
close(pb); stopCluster(cl)

names(OVind.phylo) <- Sp


# ==============================================================================
# 6. BUILD SPECIES-LEVEL SUMMARY TABLE
# ==============================================================================

species_summary <- lapply(Sp, function(sp) {
  tnw_row <- NV.phylo.df[NV.phylo.df$Sp == sp, ]
  wic_vec <- IN.phylo[[sp]]
  ov_mat  <- OVind.phylo[[sp]]
  inv_bic <- if (!is.null(ov_mat) && nrow(ov_mat) > 1)
    mean(ov_mat[lower.tri(ov_mat)], na.rm = TRUE) else NA

  data.frame(
    species       = sp,
    TNW           = tnw_row$B,
    TNW_LC        = tnw_row$LC,
    TNW_UC        = tnw_row$UC,
    mean_WIC      = mean(wic_vec, na.rm = TRUE),
    Spec_index    = mean(wic_vec, na.rm = TRUE) / tnw_row$B,
    Inverse_BIC   = inv_bic,
    n_individuals = sum(sample_data(ps_prop)$Sp == sp),
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

write.csv(species_summary,
          file.path(output_dir, "species_niche_summary_phylo.csv"),
          row.names = FALSE)


# ==============================================================================
# 7. BUILD INDIVIDUAL-LEVEL TABLE
# ==============================================================================

ind_table <- do.call(rbind, lapply(Sp, function(sp) {
  sdat <- as.data.frame(sample_data(
    prune_samples(sample_data(ps_prop)$Sp == sp, ps_prop)))
  class(sdat) <- "data.frame"
  wic_vec <- IN.phylo[[sp]]
  sp_row  <- species_summary[species_summary$species == sp, ]
  if (length(wic_vec) != nrow(sdat)) return(NULL)
  data.frame(
    sdat,
    WIC         = wic_vec,
    TNW         = sp_row$TNW,
    Spec_index  = sp_row$Spec_index,
    Inverse_BIC = sp_row$Inverse_BIC,
    Division    = ifelse(sdat$Tribe == "Hydromyini", "Hydromyini", "Rattini"),
    stringsAsFactors = FALSE
  )
}))

write.csv(ind_table,
          file.path(output_dir, "individual_niche_table_phylo.csv"),
          row.names = FALSE)

cat("Niche partitioning results saved to:", output_dir, "\n")


# ==============================================================================
# 8. VISUALISATION — Species niche space plot (WIC/TNW vs TNW)
# ==============================================================================

host_tree <- ape::read.tree(host_tree_path)

sp_sum <- read.csv(file.path(output_dir, "species_niche_summary_phylo.csv"))

p_niche <- ggplot(sp_sum, aes(x = TNW, y = Spec_index, label = species, size = n_individuals)) +
  geom_point(alpha = 0.8, colour = "steelblue") +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_size_continuous(range = c(2, 8), name = "N individuals") +
  labs(x = "Total Niche Width (TNW)", y = "Specialisation index (WIC/TNW)",
       title = "Dietary niche partitioning") +
  theme_bw()

ggsave(file.path(output_dir, "niche_space_WIC_TNW.pdf"), p_niche, width = 8, height = 6)


# ==============================================================================
# 9. PCoA — PHYLOGENETIC NICHE OVERLAP (Figure 4)
# Reads the pairwise overlap table saved by the foreach loop above (OV.phylo.df).
# Distance = 1 - overlap (sqrt-transformed); pie charts show diet composition.
# ==============================================================================

library(tibble)
library(tidyr)
library(ggrepel)
library(phyloseq)

# ---- Load pairwise niche overlap table ----
OV <- read.delim(file.path(output_dir, "OV.phylo.df.txt"))

# Mirror to symmetric matrix
symmetric_df <- bind_rows(
  OV %>% select(X1, X2, O),
  OV %>% select(X1 = X2, X2 = X1, O = O)
) %>%
  distinct(X1, X2, .keep_all = TRUE)

square_matrix <- symmetric_df %>%
  pivot_wider(names_from = X1, values_from = O) %>%
  column_to_rownames("X2") %>%
  as.matrix()

all_species   <- sort(union(symmetric_df$X1, symmetric_df$X2))
square_matrix <- square_matrix[all_species, all_species]
diag(square_matrix) <- 1

DIST <- 1 - square_matrix
PCOA <- ape::pcoa(DIST^0.5)

sp_sum_pcoa <- read.csv(file.path(output_dir, "species_niche_summary_phylo.csv"))

DD <- data.frame(PCOA$vectors) %>%
  rownames_to_column("species") %>%
  left_join(sp_sum_pcoa %>% select(species, n_individuals, TNW),
            by = "species")

pct_var <- round(PCOA$values$Relative_eig[1:2] * 100, 1)

# ---- Diet composition per species for pie charts ----
PS_pie <- readRDS(file.path(phyloseq_dir, "phyloseq_final.rds"))
PS_pie <- tax_glom(PS_pie, taxrank = "Diet")

phyloseq_species     <- merge_samples(PS_pie, "Species_consensus_ddRAD_mt_others", fun = mean)
phyloseq_species_rel <- transform_sample_counts(phyloseq_species, function(x) x / sum(x))

diet_df <- psmelt(phyloseq_species_rel) %>%
  filter(!is.na(Diet)) %>%
  group_by(Sample, Diet) %>%
  summarise(rel_abundance = sum(Abundance), .groups = "drop") %>%
  rename(species = Sample) %>%
  pivot_wider(names_from = Diet, values_from = rel_abundance, values_fill = 0)

for (d in c("Animal", "Plant", "Fungi")) {
  if (!d %in% colnames(diet_df)) diet_df[[d]] <- 0
}
diet_df <- diet_df %>%
  mutate(total  = Animal + Plant + Fungi,
         Animal = Animal / total,
         Plant  = Plant  / total,
         Fungi  = Fungi  / total) %>%
  select(species, Animal, Plant, Fungi)

DD <- DD %>% left_join(diet_df, by = "species")

# ---- Custom pie polygon function ----
create_pie_slices <- function(data, x_col, y_col, diets, radius, n_points = 30) {
  pie_list <- list()
  for (i in seq_len(nrow(data))) {
    xc    <- data[[x_col]][i]
    yc    <- data[[y_col]][i]
    sp    <- data$species[i]
    props <- as.numeric(data[i, diets])
    props[is.na(props)] <- 0
    props[props < 0.01] <- 0
    if (sum(props) == 0) next
    props <- props / sum(props)
    cum   <- c(0, cumsum(props))
    for (j in seq_along(diets)) {
      if (props[j] == 0) next
      angles <- seq(2 * pi * cum[j], 2 * pi * cum[j + 1], length.out = n_points)
      pie_list[[length(pie_list) + 1]] <- data.frame(
        species = sp, Diet = diets[j],
        x = c(xc, xc + radius * cos(angles)),
        y = c(yc, yc + radius * sin(angles))
      )
    }
  }
  bind_rows(pie_list)
}

pcoa_range <- diff(range(DD$Axis.1))
pie_radius <- pcoa_range * 0.015

pie_data <- create_pie_slices(DD, "Axis.1", "Axis.2",
                               c("Animal", "Plant", "Fungi"),
                               radius = pie_radius, n_points = 200)

diet_colors <- c("Animal" = "#E31A1C", "Plant" = "#33A02C", "Fungi" = "#1F78B4")

p_pcoa <- ggplot() +
  geom_polygon(
    data = pie_data,
    aes(x = x, y = y, fill = Diet, group = interaction(species, Diet)),
    color = NA, alpha = 0.9
  ) +
  geom_point(
    data = DD, aes(x = Axis.1, y = Axis.2),
    shape = 1, size = pie_radius * 300, color = "black", stroke = 0.4
  ) +
  geom_text_repel(
    data = DD, aes(x = Axis.1, y = Axis.2, label = species),
    size = 2.6, max.overlaps = 20, box.padding = 1.2,
    point.padding = pie_radius * 150, min.segment.length = 0,
    segment.size = 0.3, segment.color = "grey40", seed = 42
  ) +
  scale_fill_manual(values = diet_colors, name = "Diet") +
  labs(
    title    = "PCoA \u2014 Phylogenetic niche overlap between species",
    subtitle = "Distance = 1 - overlap (sqrt-transformed) | Pie charts = diet composition",
    x        = sprintf("PCoA Axis 1 (%.1f%%)", pct_var[1]),
    y        = sprintf("PCoA Axis 2 (%.1f%%)", pct_var[2])
  ) +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8),
        legend.position = "bottom") +
  coord_fixed()

ggsave(file.path(output_dir, "Fig4_PCoA_niche_overlap.png"),
       p_pcoa, width = 12, height = 10, dpi = 300)
print(p_pcoa)

cat("\nAll niche partitioning and PCoA outputs saved to:", output_dir, "\n")
cat("Done.\n")
