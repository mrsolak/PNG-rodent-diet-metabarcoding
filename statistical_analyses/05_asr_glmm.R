################################################################################
# 04_asr_glmm.R
# Ancestral state reconstruction + GLMM for dietary niche metrics
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# PART 1: Ancestral State Reconstruction (ASR) of TNW, mean_WIC, Spec_index,
#         and Inverse_BIC on the ddRAD host phylogeny using phytools::fastAnc.
# PART 2: Statistical models:
#   (a) AIC model selection for TNW (species-level, Gaussian, MuMIn::dredge)
#   (b) Beta-GLMM for Spec_index and WIC/TNW (individual-level)
#   (c) Gamma-GLMM for WIC (individual-level, glmmTMB)
# PART 3: Summary tables (sjPlot::tab_model → HTML → Word)
#
# REQUIRES output from 03_niche_partitioning.R:
#   species_niche_summary_phylo.csv
#   individual_niche_table_phylo.csv
################################################################################


# ==============================================================================
# 0. PATHS
# ==============================================================================

niche_dir  <- "/path/to/output/niche_partitioning"  # from 03_niche_partitioning.R
output_dir <- "/path/to/output/ASR_GLMM"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

asr_out_dir  <- file.path(output_dir, "ASR")
glmm_out_dir <- file.path(output_dir, "GLMM")
dir.create(asr_out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(glmm_out_dir, recursive = TRUE, showWarnings = FALSE)

data_dir       <- "/path/to/project/data"
host_tree_path <- file.path(data_dir, "ddRAD_trees", "Joined_scaled_RH85.nwk")


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(ape)
library(phytools)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(viridis)
library(glmmTMB)
library(splines)
library(MuMIn)
library(broom.mixed)
library(sjPlot)
library(ggeffects)
library(performance)
library(officer)
library(webshot2)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

sp_data  <- read.csv(file.path(niche_dir, "species_niche_summary_phylo.csv"))
ind_data <- read.csv(file.path(niche_dir, "individual_niche_table_phylo.csv"))

host_tree <- ape::read.tree(host_tree_path)

# Filter to species with >= 6 individuals (matches 03_niche_partitioning.R threshold)
sp_data <- sp_data %>% filter(n_individuals >= 6)
cat("Species included in models:", nrow(sp_data), "\n")
cat("Species excluded (<6 ind):",
    length(setdiff(unique(ind_data$Species_consensus_ddRAD_mt_others), sp_data$species)), "\n")

# Filter individual table to matching species
ind_data <- ind_data %>%
  filter(Species_consensus_ddRAD_mt_others %in% sp_data$species)

# Prune phylogeny to retained species
host_tree_pruned <- ape::keep.tip(host_tree, sp_data$species)


# ==============================================================================
# 3. PART 1 — ANCESTRAL STATE RECONSTRUCTION
# ==============================================================================

traits <- c("TNW", "mean_WIC", "Spec_index", "Inverse_BIC")

asr_results <- list()

for (trait in traits) {
  vals <- setNames(sp_data[[trait]], sp_data$species)
  # Align values to tree tips
  vals_aln <- vals[host_tree_pruned$tip.label]

  asr_out <- phytools::fastAnc(host_tree_pruned, vals_aln, vars = TRUE, CI = TRUE)
  asr_results[[trait]] <- asr_out

  # Save node estimates
  node_df <- data.frame(
    node  = names(asr_out$ace),
    est   = asr_out$ace,
    lower = asr_out$CI95[, 1],
    upper = asr_out$CI95[, 2]
  )
  write.csv(node_df,
            file.path(asr_out_dir, paste0("ASR_nodes_", trait, ".csv")),
            row.names = FALSE)

  # Plot on tree
  pdf(file.path(asr_out_dir, paste0("ASR_tree_", trait, ".pdf")), width = 8, height = 10)
  phytools::contMap(host_tree_pruned, vals_aln, plot = FALSE) %>%
    plot(legend = 0.7 * max(nodeHeights(host_tree_pruned)),
         ftype = "i", fsize = 0.7, main = trait)
  dev.off()

  cat("ASR completed:", trait, "\n")
}


# ==============================================================================
# 4. PART 2a — SPECIES-LEVEL TNW MODEL SELECTION (Gaussian, dredge)
# ==============================================================================

# Merge trait data with species-level predictors
sp_model_data <- sp_data %>%
  left_join(
    ind_data %>%
      group_by(Species_consensus_ddRAD_mt_others) %>%
      summarise(
        spatial_niche = names(sort(table(spatial_niche), decreasing = TRUE))[1],
        Tribe         = names(sort(table(Tribe),         decreasing = TRUE))[1],
        mean_mass     = mean(mean_mass, na.rm = TRUE),
        Range_Type    = names(sort(table(Metalokalita),  decreasing = TRUE))[1],
        .groups = "drop"
      ),
    by = c("species" = "Species_consensus_ddRAD_mt_others")
  )

# Global model
options(na.action = "na.fail")
tnw_global <- lm(TNW ~ spatial_niche + Tribe + mean_mass + Range_Type,
                 data = sp_model_data)

tnw_dredge <- MuMIn::dredge(tnw_global, m.lim = c(0, 3))
write.csv(as.data.frame(tnw_dredge),
          file.path(glmm_out_dir, "TNW_model_selection.csv"),
          row.names = FALSE)
options(na.action = "na.omit")


# ==============================================================================
# 4. PART 2b — INDIVIDUAL-LEVEL BETA-GLMM FOR SPEC_INDEX
# ==============================================================================

# Prepare individual data
ind_clean <- ind_data %>%
  filter(!is.na(Alt), !is.na(mean_mass), !is.na(spatial_niche),
         !is.na(Tribe), !is.na(Habitat_polished_Daniel.Fanda),
         !is.na(Habitat_water_body_prox_Daniel.Fanda),
         Spec_index > 0, Spec_index < 1)

# Beta-GLMM: Spec_index ~ altitude (linear) + random effects
m_beta_linear <- glmmTMB::glmmTMB(
  Spec_index ~ scale(Alt) + spatial_niche + Tribe + scale(mean_mass) +
    (1 | Species_consensus_ddRAD_mt_others) + (1 | Locality),
  data   = ind_clean,
  family = glmmTMB::beta_family()
)

# Beta-GLMM: Spec_index ~ altitude (non-linear with natural spline)
m_beta_spline <- glmmTMB::glmmTMB(
  Spec_index ~ ns(Alt, df = 3) + spatial_niche + Tribe + scale(mean_mass) +
    (1 | Species_consensus_ddRAD_mt_others) + (1 | Locality),
  data   = ind_clean,
  family = glmmTMB::beta_family()
)

# Save model summaries
save(m_beta_linear, m_beta_spline,
     file = file.path(glmm_out_dir, "beta_glmm_Spec_index.RData"))


# ==============================================================================
# 4. PART 2c — INDIVIDUAL-LEVEL GAMMA-GLMM FOR WIC
# ==============================================================================

ind_wic <- ind_data %>%
  filter(!is.na(Alt), !is.na(mean_mass), !is.na(spatial_niche),
         !is.na(Tribe), WIC > 0)

m_gamma_wic <- glmmTMB::glmmTMB(
  WIC ~ scale(Alt) + spatial_niche + Tribe + scale(mean_mass) +
    (1 | Species_consensus_ddRAD_mt_others) + (1 | Locality),
  data   = ind_wic,
  family = Gamma(link = "log")
)

save(m_gamma_wic, file = file.path(glmm_out_dir, "gamma_glmm_WIC.RData"))


# ==============================================================================
# 5. PART 3 — PUBLICATION TABLES (HTML → Word)
# ==============================================================================

tab_beta <- sjPlot::tab_model(
  m_beta_linear, m_beta_spline,
  dv.labels = c("Spec_index (linear Alt)", "Spec_index (spline Alt)"),
  file = file.path(glmm_out_dir, "table_beta_Spec_index.html")
)

tab_gamma <- sjPlot::tab_model(
  m_gamma_wic,
  dv.labels = "WIC (Gamma GLMM)",
  file = file.path(glmm_out_dir, "table_gamma_WIC.html")
)

# Combine tables into a Word document
word_doc <- officer::read_docx()
for (html_file in list.files(glmm_out_dir, pattern = "\\.html$", full.names = TRUE)) {
  png_file <- sub("\\.html$", ".png", html_file)
  webshot2::webshot(html_file, png_file, zoom = 2)
  word_doc <- officer::body_add_img(word_doc, png_file, width = 6, height = 4)
  word_doc <- officer::body_add_break(word_doc)
}
print(word_doc, target = file.path(glmm_out_dir, "GLMM_tables_combined.docx"))

cat("ASR and GLMM results saved to:", output_dir, "\n")
