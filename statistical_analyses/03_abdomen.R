################################################################################
# 07_abdomen.R
# Phylosymbiosis analysis using the ABDOMEN model
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# ABDOMEN (Ancestral Diet By Ornstein–Uhlenbeck model of Multivariate
# Evolutionary Niche) requires the ABDOMEN.R source file to be available.
# See: https://github.com/BPerezLamarque/ABDOMEN
#
# Input:  phyloseq_final.rds  (or load phyloseq_final.R via load())
#         ddRAD host phylogeny (Newick format)
#         ABDOMEN.R source file
# Output: ABDOMEN fit + ancestral compositions + permutation test results
#
# Runtime: VERY LONG (hours–days depending on nb_cores and permutations)
# Runs on:  Stan via rstan; set nb_cores and chains appropriately
################################################################################

library(phyloseq)
library(ape)
library(dplyr)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(mvMORPH)
library(RPANDA)
library(rstan)
rstan_options(auto_write = TRUE)
library(RColorBrewer)
library(ggtree)
library(scatterpie)
options(scipen = 100)

# ==============================================================================
# PATHS — edit before running
# ==============================================================================

data_dir       <- "/path/to/project/data"
phyloseq_dir   <- file.path(data_dir, "phyloseqs")
host_tree_path <- file.path(data_dir, "ddRAD_trees", "Joined_scaled_RH85.nwk")
abdomen_script <- "/path/to/ABDOMEN.R"   # source file for ABDOMEN functions
output_dir     <- "/path/to/output/abdomen"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. SOURCE ABDOMEN
# ==============================================================================

source(abdomen_script)

# ==============================================================================
# 2. LOAD PHYLOSEQ AND PREPARE DATA
# ==============================================================================

ps <- readRDS(file.path(phyloseq_dir, "phyloseq_final.rds"))
# Alternative: load(file.path(phyloseq_dir, "phyloseq_final.R"))
# object name: PNG_18S_rodent_filtered_ensembles_tree_NEWMT_W_clean

# Add Division column (consolidate Rattus groups)
sample_data(ps)$Division <- sample_data(ps)$Hydromyini.Division.Rattini.Group
sample_data(ps)$Division <- ifelse(
  sample_data(ps)$Division %in% c("incertae sedis group", "R. rattus group "),
  "Rattus invasive",
  sample_data(ps)$Division
)
sample_data(ps)$Division <- ifelse(
  sample_data(ps)$Division == "R. leucopus group (New Guinean group) ",
  "Rattus endemic",
  sample_data(ps)$Division
)

# Aggregate to Diet level (Animal / Plant / Fungi)
PS_diet <- tax_glom(ps, taxrank = "Diet")

# ==============================================================================
# 3. LOAD HOST TREE AND CREATE SPECIES-LEVEL TREE
# ==============================================================================

host_tree_raw <- ape::read.tree(host_tree_path)

# Clean tip labels to sample IDs
host_tree_raw$tip.label <- sub("_.*$", "", host_tree_raw$tip.label)

ps_sample_data <- data.frame(sample_data(PS_diet))
ps_sample_data$SampleID <- rownames(ps_sample_data)

# Map tree tips to species
tree_tip_species <- sapply(host_tree_raw$tip.label, function(tip) {
  species <- ps_sample_data[ps_sample_data$SampleID == tip,
                             "Species_consensus_ddRAD_mt_others"]
  if (length(species) > 0) return(as.character(species[1]))
  return(NA)
})

tips_to_keep    <- host_tree_raw$tip.label[!is.na(tree_tip_species)]
host_tree_matched <- ape::keep.tip(host_tree_raw, tips_to_keep)
kept_species    <- tree_tip_species[!is.na(tree_tip_species)]

unique_species     <- unique(kept_species)
representative_tips <- sapply(unique_species, function(sp) {
  tips_to_keep[kept_species == sp][1]
})

tree <- ape::keep.tip(host_tree_matched, representative_tips)
species_for_tips <- sapply(tree$tip.label, function(tip) {
  sp <- ps_sample_data[ps_sample_data$SampleID == tip,
                        "Species_consensus_ddRAD_mt_others"]
  if (length(sp) > 0) return(as.character(sp[1]))
  return(NA)
})
tree$tip.label <- species_for_tips

# Ensure rooted, binary, ultrametric
if (!is.rooted(tree))    tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
if (!is.binary(tree))    tree <- multi2di(tree)
if (!is.ultrametric(tree)) tree <- chronos(tree)

# ==============================================================================
# 4. CALCULATE DIET PROPORTIONS PER SPECIES (Animal / Plant / Fungi)
# ==============================================================================

PS_species     <- merge_samples(PS_diet, "Species_consensus_ddRAD_mt_others", fun = mean)
PS_species_rel <- transform_sample_counts(PS_species, function(x) x / sum(x))

otu_table_species <- as(otu_table(PS_species_rel), "matrix")
if (ncol(otu_table_species) > nrow(otu_table_species))
  otu_table_species <- t(otu_table_species)

tax_table_species <- as(tax_table(PS_species_rel), "matrix")

diet_matrix <- matrix(0, nrow = nrow(otu_table_species), ncol = 3,
                      dimnames = list(rownames(otu_table_species),
                                      c("Animal", "Plant", "Fungi")))

for (sp in rownames(otu_table_species)) {
  for (asv in colnames(otu_table_species)) {
    if (asv %in% rownames(tax_table_species)) {
      diet_type <- tax_table_species[asv, "Diet"]
      if (!is.na(diet_type) && diet_type %in% c("Animal", "Plant", "Fungi"))
        diet_matrix[sp, diet_type] <- diet_matrix[sp, diet_type] + otu_table_species[sp, asv]
    }
  }
}

shared_species <- intersect(tree$tip.label, rownames(diet_matrix))
table <- diet_matrix[shared_species, ]
tree  <- keep.tip(tree, shared_species)
table <- table[tree$tip.label, ]

# ==============================================================================
# 5. ABDOMEN PARAMETERS
# ==============================================================================

name                <- "rodent_diet_proportions"
code_path           <- getwd()
detection_threshold <- 1e-05
seed                <- 3
prior_Z0            <- "empirical"
mean_prior_logY     <- 0
sd_prior_logY       <- 2
nb_cores            <- 10
chains              <- 4
warmup              <- 3000
iter                <- 8000

# ==============================================================================
# 6. RUN ABDOMEN
# ==============================================================================

fit_summary <- ABDOMEN(tree, table, name,
                       code_path           = code_path,
                       prior_Z0            = prior_Z0,
                       detection_threshold = detection_threshold,
                       seed                = seed,
                       mean_prior_logY     = mean_prior_logY,
                       sd_prior_logY       = sd_prior_logY,
                       nb_cores            = nb_cores,
                       chains              = chains,
                       warmup              = warmup,
                       iter                = iter)

# ==============================================================================
# 7. EXTRACT RESULTS
# ==============================================================================

original_lambda <- ABDOMEN_extract_lambda(tree, table, fit_summary,
                                          detection_threshold = detection_threshold)

Z0 <- ABDOMEN_extract_Z0(tree, table, fit_summary,
                          detection_threshold = detection_threshold)

Z0_nodes <- ABDOMEN_extract_Z0_nodes(tree, table, fit_summary,
                                      detection_threshold = detection_threshold)

R_matrices <- ABDOMEN_extract_R(tree, table, fit_summary,
                                 detection_threshold = detection_threshold)

# Generate default ABDOMEN plots
list_colors <- c("Animal" = "#f8766d", "Plant" = "#00ba38", "Fungi" = "#619cff")
ABDOMEN_process_output(tree, table, name, fit_summary,
                       code_path           = code_path,
                       detection_threshold = detection_threshold,
                       list_colors         = list_colors)

cat("Pagel's lambda:\n"); print(original_lambda)
cat("Root ancestral composition:\n"); print(Z0)

# Save fit_summary
save(fit_summary, file = file.path(output_dir, "fit_summary.RData"))
cat("fit_summary saved.\n")

# ==============================================================================
# 8. PERMUTATION TEST (optional — very slow)
# ==============================================================================

nb_permutations         <- 100
list_lambda_permutations <- c()

for (seed_perm in 1:nb_permutations) {
  set.seed(seed_perm)
  name_random  <- paste0("rodent_diet_permutation_", seed_perm)
  table_random <- table[sample(tree$tip.label), ]
  rownames(table_random) <- rownames(table)

  fit_summary_permut <- ABDOMEN(tree, table_random, name = name_random,
                                code_path           = code_path,
                                prior_Z0            = prior_Z0,
                                detection_threshold = detection_threshold,
                                seed                = seed_perm,
                                mean_prior_logY     = mean_prior_logY,
                                sd_prior_logY       = sd_prior_logY,
                                nb_cores            = nb_cores,
                                chains              = chains,
                                warmup              = warmup,
                                iter                = iter)

  list_lambda_permutations <- rbind(
    list_lambda_permutations,
    ABDOMEN_extract_lambda(tree, table_random, fit_summary_permut,
                           detection_threshold = detection_threshold)
  )
}

p_value <- length(which(list_lambda_permutations[, 1] >= original_lambda[1])) /
           nb_permutations

write.csv(list_lambda_permutations,
          file.path(output_dir, "abdomen_lambda_permutations.csv"),
          row.names = FALSE)

p_value_results <- data.frame(
  original_lambda = original_lambda[1],
  p_value         = p_value,
  nb_permutations = nb_permutations,
  significant     = ifelse(p_value < 0.05, "YES", "NO")
)
write.csv(p_value_results,
          file.path(output_dir, "abdomen_phylosymbiosis_pvalue.csv"),
          row.names = FALSE)

save(list_lambda_permutations, p_value, original_lambda, nb_permutations,
     Z0_nodes, R_matrices, Z0,
     file = file.path(output_dir, "abdomen_permutation_results.RData"))

save(Z0,       file = file.path(output_dir, "Z0.Rdata"))
save(Z0_nodes, file = file.path(output_dir, "Z0_nodes.Rdata"))

# ==============================================================================
# 9. POLISHED DIET TREE VISUALISATION
# ==============================================================================

diet_cols <- c("Animal", "Plant", "Fungi")

# Get species counts per species
sample_counts <- sample_data(ps) %>%
  as("data.frame") %>%
  dplyr::count(Species_consensus_ddRAD_mt_others) %>%
  dplyr::rename(Species = Species_consensus_ddRAD_mt_others, n_samples = n)

# Species divisions
species_divisions <- sample_data(ps) %>%
  as("data.frame") %>%
  group_by(Species_consensus_ddRAD_mt_others) %>%
  summarize(Division = first(Division)) %>%
  ungroup() %>%
  rename(Species = Species_consensus_ddRAD_mt_others)

division_colors <- c(
  "Coccymys division"  = "#E60000",
  "Hydromys division"  = "#004D99",
  "Mallomys division"  = "#33CC33",
  "Pogonomys division" = "#FFCC00",
  "Rattus endemic"     = "#660066",
  "Rattus invasive"    = "#FF6600",
  "Uromys division"    = "#00CCCC"
)

# Scale tree X/Y
p_tmp       <- ggtree(tree, layout = "rectangular")
tmp_df      <- p_tmp$data
x_range     <- max(tmp_df$x) - min(tmp_df$x)
y_range     <- diff(range(tmp_df$y))
scale_factor <- if (x_range > 0) y_range / x_range else 1

tree_scaled              <- tree
tree_scaled$edge.length  <- tree_scaled$edge.length * scale_factor

p_base      <- ggtree(tree_scaled, layout = "rectangular", size = 0.8)
tree_coords <- p_base$data

# Pie chart data
Z0_nodes_df <- as.data.frame(Z0_nodes)
if (!"node" %in% colnames(Z0_nodes_df))
  Z0_nodes_df$node <- as.numeric(rownames(Z0_nodes_df))

internal_data <- Z0_nodes_df %>%
  select(node, all_of(diet_cols)) %>%
  mutate(across(all_of(diet_cols), as.numeric))

tip_data <- as.data.frame(diet_matrix[tree$tip.label, ])
tip_data$node <- seq_along(tree$tip.label)
tip_data <- tip_data %>% select(node, all_of(diet_cols))

all_node_data <- bind_rows(tip_data, internal_data)

piechart_data <- tree_coords %>%
  left_join(all_node_data, by = "node") %>%
  filter(!is.na(Animal))

internal_pies <- piechart_data %>% filter(!isTip)
tip_pies      <- piechart_data %>% filter(isTip)

pie_radius   <- y_range * 0.01
label_offset <- pie_radius * 1.8

max_orig    <- max(tree$edge.length)
max_rounded <- ceiling(max_orig * 5) / 5

# Division backgrounds
tree_data_with_divisions <- p_base$data %>%
  left_join(species_divisions, by = c("label" = "Species"))

division_backgrounds <- tree_data_with_divisions %>%
  filter(isTip) %>%
  group_by(Division) %>%
  summarize(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    xmin = -0.05 * max(tree_coords$x),
    xmax = max(tree_coords$x) * 1.7
  )

division_labels <- division_backgrounds %>%
  mutate(
    y_center   = (ymin + ymax) / 2,
    x_position = max(tree_coords$x) * 1.55
  )

# Annotation text (lambda + p)
annotation_text <- paste0(
  "\u03bb = ", round(original_lambda[1], 2),
  ", p = ", ifelse(p_value < 0.001, "<0.001", round(p_value, 3))
)

p_tree <- ggtree(tree_scaled, layout = "rectangular", size = 0.8) +

  geom_rect(
    data = division_backgrounds,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill  = rep(division_colors[division_backgrounds$Division],
                length.out = nrow(division_backgrounds)),
    alpha = 0.2, inherit.aes = FALSE
  ) +

  geom_text(
    data = division_labels,
    aes(x = x_position, y = y_center, label = Division),
    size = 3.5, hjust = 0.3, fontface = "bold", inherit.aes = FALSE
  ) +

  geom_rootedge(rootedge = tree_scaled$edge.length[1] * 0.3,
                linetype = "dashed", size = 0.8) +

  geom_scatterpie(data = internal_pies,
                  aes(x = x, y = y, r = pie_radius),
                  cols = diet_cols, color = "black", size = 0.3) +

  geom_scatterpie(data = tip_pies,
                  aes(x = x, y = y, r = pie_radius),
                  cols = diet_cols, color = "black", size = 0.3) +

  geom_tiplab(
    aes(label = paste0(label, " (n=",
                       sample_counts$n_samples[match(label, sample_counts$Species)],
                       ")")),
    size = 3.0, offset = label_offset, hjust = 0, fontface = "italic"
  ) +

  annotate("text",
           x = 0.02 * max(tree_coords$x),
           y = max(tree_coords$y) * 0.98,
           label = annotation_text,
           size = 4, hjust = 0, vjust = 1, fontface = "bold") +

  scale_fill_manual(values = list_colors, name = "Diet Composition",
                    guide  = guide_legend(nrow = 1,
                                          override.aes = list(size = 4))) +

  xlim(0, max(tree_coords$x) * 1.8) +
  coord_cartesian(clip = "off") +
  ggtitle("Ancestral Reconstruction of Diet Composition") +
  theme_tree2() +
  theme(
    plot.title        = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 9),
    axis.title.x      = element_text(size = 12),
    plot.margin       = margin(10, 100, 30, 10),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white", color = "white")
  ) +
  scale_x_continuous(
    name   = "substitutions/site",
    breaks = seq(0, max_rounded, by = 0.2) * scale_factor,
    labels = seq(0, max_rounded, by = 0.2),
    expand = expansion(mult = c(0.02, 0.20))
  )

print(p_tree)

ggsave(file.path(output_dir, "rodent_diet_tree_ABDOMEN.png"),
       p_tree, width = 9, height = 18, dpi = 300, bg = "white",
       limitsize = FALSE)
ggsave(file.path(output_dir, "rodent_diet_tree_ABDOMEN.pdf"),
       p_tree, width = 18, height = 9, device = "pdf", limitsize = FALSE)

# ==============================================================================
# 10. EXTRACT ANCESTRAL DIET COMPOSITIONS FOR MANUSCRIPT
# ==============================================================================

# Root ancestor — from fit_summary$summary (has proper CIs)
z0_rows    <- grep("^Z0\\[", rownames(fit_summary$summary))
z0_summary <- fit_summary$summary[z0_rows, ]
rownames(z0_summary) <- diet_cols

root_df <- data.frame(
  ancestor   = "Root (all rodents)",
  diet       = diet_cols,
  mean_pct   = round(z0_summary[, "mean"]  * 100, 1),
  CI_low_95  = round(z0_summary[, "2.5%"]  * 100, 1),
  CI_high_95 = round(z0_summary[, "97.5%"] * 100, 1),
  CI_low_50  = round(z0_summary[, "25%"]   * 100, 1),
  CI_high_50 = round(z0_summary[, "75%"]   * 100, 1),
  has_CI     = TRUE,
  stringsAsFactors = FALSE
)

# Helper: find MRCA node in Z0_nodes
find_node <- function(Z0_nodes, required_species, label,
                      most_inclusive = FALSE) {
  matches <- sapply(Z0_nodes$MRCA, function(mrca_str) {
    mrca_species <- strsplit(mrca_str, "-")[[1]]
    all(required_species %in% mrca_species)
  })
  matched_rows  <- which(matches)
  if (length(matched_rows) == 0) {
    cat("WARNING: No node found for", label, "\n"); return(NULL)
  }
  mrca_lengths  <- sapply(Z0_nodes$MRCA[matched_rows],
                          function(x) length(strsplit(x, "-")[[1]]))
  best_row <- if (most_inclusive) matched_rows[which.max(mrca_lengths)]
              else                matched_rows[which.min(mrca_lengths)]
  cat("Found node for", label, ":", rownames(Z0_nodes)[best_row],
      "(node", Z0_nodes$node[best_row], ")\n")
  best_row
}

rattus_all_str      <- sort(c("Rattus rattus", "Rattus exulans",
                               "Rattus steini", "Rattus niobe", "Rattus verecundus"))
rattus_endemic_str  <- sort(c("Rattus steini", "Rattus niobe", "Rattus verecundus"))
rattus_invasive_str <- sort(c("Rattus rattus", "Rattus exulans"))

Z0_nodes_no_rattus <- Z0_nodes[!grepl("Rattus", Z0_nodes$MRCA), ]

hydromyini_row_nr   <- find_node(Z0_nodes_no_rattus,
                                  c("Abeomelomys spMW", "Hydromys sp",
                                    "Mallomys istapantap", "Pseudohydromys ellermani"),
                                  "Hydromyini ancestor", most_inclusive = TRUE)
hydromyini_row_full <- which(rownames(Z0_nodes) ==
                               rownames(Z0_nodes_no_rattus)[hydromyini_row_nr])

rattini_row  <- find_node(Z0_nodes, rattus_all_str,      "All Rattini ancestor")
endemic_row  <- find_node(Z0_nodes, rattus_endemic_str,  "Endemic Rattini ancestor")
invasive_row <- find_node(Z0_nodes, rattus_invasive_str, "Invasive Rattus ancestor")

extract_node_row <- function(Z0_nodes_df, row_idx, ancestor_label) {
  if (is.null(row_idx)) return(NULL)
  node_data <- Z0_nodes_df[row_idx, diet_cols]
  data.frame(
    ancestor   = ancestor_label,
    diet       = diet_cols,
    mean_pct   = round(as.numeric(node_data) * 100, 1),
    CI_low_95  = NA, CI_high_95 = NA,
    CI_low_50  = NA, CI_high_50 = NA,
    has_CI     = FALSE,
    stringsAsFactors = FALSE
  )
}

nodes_df <- bind_rows(
  extract_node_row(Z0_nodes, hydromyini_row_full, "Hydromyini ancestor"),
  extract_node_row(Z0_nodes, rattini_row,          "All Rattini ancestor"),
  extract_node_row(Z0_nodes, endemic_row,          "Endemic Rattini ancestor"),
  extract_node_row(Z0_nodes, invasive_row,         "Invasive Rattus ancestor")
)

all_results <- bind_rows(root_df, nodes_df)

cat("\n=== FORMATTED FOR MANUSCRIPT ===\n")
for (anc in unique(all_results$ancestor)) {
  sub <- all_results %>% filter(ancestor == anc)
  cat("\n", anc, ":\n")
  for (i in seq_len(nrow(sub))) {
    if (sub$has_CI[i]) {
      cat(sprintf("  %s: %.1f%% (95%% CI: %.1f%%\u2013%.1f%%)\n",
                  sub$diet[i], sub$mean_pct[i],
                  sub$CI_low_95[i], sub$CI_high_95[i]))
    } else {
      cat(sprintf("  %s: %.1f%%\n", sub$diet[i], sub$mean_pct[i]))
    }
  }
}

write.csv(all_results,
          file.path(output_dir, "ancestral_diet_compositions.csv"),
          row.names = FALSE)

if (!file.exists(file.path(output_dir, "fit_summary.RData")))
  save(fit_summary, file = file.path(output_dir, "fit_summary.RData"))

cat("\nAll outputs saved to:", output_dir, "\n")
cat("Done.\n")
