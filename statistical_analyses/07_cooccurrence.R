################################################################################
# 05_cooccurrence.R
# Individual-level dietary co-occurrence analysis
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# Input:  phyloseq_final.rds  (or load phyloseq_final.R via load())
# Output: CSVs + PNGs in output_dir
#
# Design: 2x2 factorial (same/diff species x same/diff locality)
# Metrics: Bray-Curtis + Weighted UniFrac (proportion-normalised)
# Model: distance ~ same_locality_f * same_species_f + (1|sp1) + (1|sp2) + (1|loc1)
################################################################################

library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(scales)

# ==============================================================================
# PATHS — edit before running
# ==============================================================================

data_dir     <- "/path/to/project/data"
phyloseq_dir <- file.path(data_dir, "phyloseqs")
output_dir   <- "/path/to/output/cooccurrence"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. LOAD PHYLOSEQ
# ==============================================================================

ps <- readRDS(file.path(phyloseq_dir, "phyloseq_final.rds"))
# Alternative: load(file.path(phyloseq_dir, "phyloseq_final.R"))
# object name: PNG_18S_rodent_filtered_ensembles_tree_NEWMT_W_clean

sample_df <- as.data.frame(sample_data(ps))
class(sample_df) <- "data.frame"
sample_df$SampleID <- rownames(sample_df)

sample_df <- sample_df %>%
  mutate(Division = case_when(
    Hydromyini.Division.Rattini.Group %in%
      c("incertae sedis group", "R. rattus group ") ~ "Rattus invasive",
    Hydromyini.Division.Rattini.Group ==
      "R. leucopus group (New Guinean group) "       ~ "Rattus endemic",
    TRUE ~ Hydromyini.Division.Rattini.Group
  ))

cat("Total samples before filtering:", nrow(sample_df), "\n")

# ==============================================================================
# 2. FILTER: EXCLUDE MW2400
# ==============================================================================

ps_filtered <- subset_samples(ps, Locality != "MW2400")
sample_df_filtered <- sample_df %>%
  filter(Locality != "MW2400",
         SampleID %in% sample_names(ps_filtered))

sample_df_filtered <- sample_df_filtered %>%
  mutate(Locality_Habitat = paste(Locality,
                                   Habitat_polished_Daniel.Fanda,
                                   sep = "_"))

cat("Samples after excluding MW2400:", nsamples(ps_filtered), "\n")
cat("Unique localities retained:", length(unique(sample_df_filtered$Locality)), "\n")
cat("Unique Locality+Habitat combinations:",
    length(unique(sample_df_filtered$Locality_Habitat)), "\n")

# ==============================================================================
# 3. PROPORTION NORMALISATION
# ==============================================================================

ps_prop <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
cat("Proportion normalisation applied.\n")

# ==============================================================================
# 4. COMPUTE DISTANCE MATRICES
# ==============================================================================

cat("Computing Bray-Curtis distances...\n")
bc_dist <- phyloseq::distance(ps_prop, method = "bray")
bc_mat  <- as.matrix(bc_dist)
cat("Bray-Curtis matrix:", nrow(bc_mat), "x", ncol(bc_mat), "\n")

cat("Computing Weighted UniFrac distances...\n")
wu_dist <- phyloseq::distance(ps_filtered, method = "wunifrac")
wu_mat  <- as.matrix(wu_dist)
cat("Weighted UniFrac matrix:", nrow(wu_mat), "x", ncol(wu_mat), "\n")

# ==============================================================================
# 5. BUILD PAIRS — 2x2 FACTORIAL (upper.tri indexing)
# ==============================================================================

cat("Building all pairs for 2x2 factorial design...\n")

sdf         <- sample_df_filtered %>% filter(SampleID %in% rownames(bc_mat))
sample_ids  <- sdf$SampleID
n_ind       <- length(sample_ids)

species_id  <- setNames(sdf$Species_consensus_ddRAD_mt_others, sdf$SampleID)
loc_hab_id  <- setNames(sdf$Locality_Habitat,                  sdf$SampleID)
locality_id <- setNames(sdf$Locality,                          sdf$SampleID)
division_id <- setNames(sdf$Division,                          sdf$SampleID)

idx      <- which(upper.tri(matrix(0, n_ind, n_ind)), arr.ind = TRUE)
ind1_vec <- sample_ids[idx[, 1]]
ind2_vec <- sample_ids[idx[, 2]]
sp1_vec  <- species_id[ind1_vec]
sp2_vec  <- species_id[ind2_vec]
lh1_vec  <- loc_hab_id[ind1_vec]
lh2_vec  <- loc_hab_id[ind2_vec]
loc1_vec <- locality_id[ind1_vec]
div1_vec <- division_id[ind1_vec]
div2_vec <- division_id[ind2_vec]

same_species_vec  <- sp1_vec == sp2_vec
same_locality_vec <- lh1_vec == lh2_vec

pairs_base <- data.frame(
  ind1          = ind1_vec,
  ind2          = ind2_vec,
  sp1           = sp1_vec,
  sp2           = sp2_vec,
  loc_hab1      = lh1_vec,
  loc_hab2      = lh2_vec,
  loc1          = loc1_vec,
  div1          = div1_vec,
  div2          = div2_vec,
  same_species  = same_species_vec,
  same_locality = same_locality_vec,
  stringsAsFactors = FALSE
) %>%
  mutate(
    group = case_when(
      same_species  & same_locality  ~ "Intraspecific\nSame locality",
      same_species  & !same_locality ~ "Intraspecific\nDiff. locality",
      !same_species & same_locality  ~ "Interspecific\nSame locality",
      !same_species & !same_locality ~ "Interspecific\nDiff. locality"
    ),
    group = factor(group, levels = c(
      "Intraspecific\nSame locality",
      "Intraspecific\nDiff. locality",
      "Interspecific\nSame locality",
      "Interspecific\nDiff. locality"
    )),
    same_species_f  = factor(same_species,
                             levels = c(TRUE, FALSE),
                             labels = c("Same species", "Diff. species")),
    same_locality_f = factor(same_locality,
                             levels = c(TRUE, FALSE),
                             labels = c("Same locality", "Diff. locality"))
  )

cat(sprintf("Total pairs:                    %d\n", nrow(pairs_base)))
cat(sprintf("  Intraspecific same locality:  %d\n",
            sum(pairs_base$same_species  & pairs_base$same_locality)))
cat(sprintf("  Intraspecific diff locality:  %d\n",
            sum(pairs_base$same_species  & !pairs_base$same_locality)))
cat(sprintf("  Interspecific same locality:  %d\n",
            sum(!pairs_base$same_species & pairs_base$same_locality)))
cat(sprintf("  Interspecific diff locality:  %d\n",
            sum(!pairs_base$same_species & !pairs_base$same_locality)))

# ==============================================================================
# 6. COLOUR SCHEME
# Locality = primary axis: same locality = dark, diff locality = light
# Species  = hue:          intraspecific  = blue, interspecific = red
# ==============================================================================

group_colours <- c(
  "Intraspecific\nSame locality"  = "#1B4F8A",   # dark blue
  "Intraspecific\nDiff. locality" = "#89B4DA",   # light blue
  "Interspecific\nSame locality"  = "#922B21",   # dark red
  "Interspecific\nDiff. locality" = "#F1948A"    # light red
)

# ==============================================================================
# 7. HELPER FUNCTIONS
# ==============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  return(sprintf("p = %.3f", p))
}

add_bracket <- function(plot, x1, x2, y_bracket, label,
                         tip_height = 0.002, text_size = 3.2) {
  plot +
    annotate("segment",
             x = x1, xend = x1,
             y = y_bracket - tip_height, yend = y_bracket,
             colour = "grey20", linewidth = 0.5) +
    annotate("segment",
             x = x2, xend = x2,
             y = y_bracket - tip_height, yend = y_bracket,
             colour = "grey20", linewidth = 0.5) +
    annotate("segment",
             x = x1, xend = x2,
             y = y_bracket, yend = y_bracket,
             colour = "grey20", linewidth = 0.5) +
    annotate("text",
             x     = (x1 + x2) / 2,
             y     = y_bracket,
             label = label,
             vjust = -0.4,
             size  = text_size,
             colour = "grey20")
}

make_barplot <- function(bar_data, dist_label, int_p_label,
                          p_same_loc, p_diff_loc,
                          error_label,
                          group_colours) {

  y_min       <- floor(min(bar_data$err_lo) * 20) / 20
  y_err_max   <- max(bar_data$err_hi)
  bracket_gap <- (y_err_max - y_min) * 0.06
  bracket1_y  <- y_err_max + bracket_gap
  bracket2_y  <- bracket1_y + bracket_gap * 1.8
  y_max       <- bracket2_y + bracket_gap * 2.5
  tip_h       <- (y_max - y_min) * 0.012

  p <- ggplot(bar_data, aes(x = group, y = mean_dist, fill = group)) +
    geom_col(width = 0.65, colour = "grey30", linewidth = 0.4) +
    geom_errorbar(aes(ymin = err_lo, ymax = err_hi),
                  width = 0.18, linewidth = 0.7, colour = "grey20") +
    geom_text(aes(y = err_hi, label = n_label),
              vjust = -0.4, size = 2.6, colour = "grey40") +
    scale_fill_manual(values = group_colours, guide = "none") +
    scale_y_continuous(limits = c(y_min, y_max),
                       oob    = scales::squish) +
    labs(x = NULL,
         y = paste0("Mean ", dist_label, "\ndissimilarity (", error_label, ")")) +
    annotate("text",
             x = 3.8, y = y_max * 0.995,
             label    = int_p_label,
             hjust    = 1, vjust = 1,
             size     = 3, fontface = "italic", colour = "grey30") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x  = element_text(size = 9,  colour = "grey20"),
      axis.title.y = element_text(size = 10),
      plot.margin  = margin(12, 8, 8, 8)
    )

  # Bracket 1: Intraspecific same (x=1) vs Interspecific same (x=3)
  p <- add_bracket(p, x1 = 1, x2 = 3, y_bracket = bracket1_y,
                   label = fmt_p(p_same_loc), tip_height = tip_h,
                   text_size = 3.2)

  # Bracket 2: Intraspecific diff (x=2) vs Interspecific diff (x=4)
  p <- add_bracket(p, x1 = 2, x2 = 4, y_bracket = bracket2_y,
                   label = fmt_p(p_diff_loc), tip_height = tip_h,
                   text_size = 3.2)
  p
}

# ==============================================================================
# 8. ANALYSIS FUNCTION
# ==============================================================================

run_metric_analysis <- function(dist_mat, dist_label, dist_label_short,
                                 pairs_base, output_dir,
                                 group_colours,
                                 max_pairs = 500000) {

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("METRIC:", dist_label, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  pairs <- pairs_base %>%
    mutate(distance = dist_mat[cbind(ind1, ind2)])

  # ---- LMM ----
  cat("Fitting LMM...\n")
  set.seed(42)
  if (nrow(pairs) > max_pairs) {
    cat(sprintf("  Subsampling to %d pairs\n", max_pairs))
    pairs_lmm <- pairs[sample(nrow(pairs), max_pairs), ]
  } else {
    pairs_lmm <- pairs
  }

  mod <- lmer(
    distance ~ same_locality_f * same_species_f +
      (1 | sp1) + (1 | sp2) + (1 | loc1),
    data    = pairs_lmm,
    REML    = TRUE,
    control = lmerControl(optimizer = "bobyqa",
                          optCtrl   = list(maxfun = 2e5))
  )

  cat("\n--- LMM Summary ---\n")
  print(summary(mod))

  anova_tab <- anova(mod)
  cat("\n--- ANOVA Table ---\n")
  print(anova_tab)

  # ---- Emmeans ----
  emm    <- emmeans(mod, ~ same_locality_f * same_species_f)
  emm_df <- as.data.frame(emm)
  names(emm_df)[names(emm_df) %in% c("asymp.LCL", "lower.CL")] <- "LCL"
  names(emm_df)[names(emm_df) %in% c("asymp.UCL", "upper.CL")] <- "UCL"
  emm_df <- emm_df %>%
    mutate(
      species_label  = factor(same_species_f,
                               levels = c("Same species", "Diff. species")),
      locality_label = factor(same_locality_f,
                               levels = c("Same locality", "Diff. locality"))
    )

  all_contrasts <- as.data.frame(pairs(emm, adjust = "tukey"))
  cat("\n--- All pairwise contrasts (Tukey) ---\n")
  print(all_contrasts)

  int_contrast <- contrast(emm, interaction = "pairwise")
  simple_loc   <- emmeans(mod, ~ same_locality_f | same_species_f)
  cat("\n--- Interaction Contrasts ---\n"); print(int_contrast)
  cat("\n--- Simple Effects of Locality ---\n"); print(pairs(simple_loc))

  # ---- Extract PI-requested pairwise p-values ----
  get_contrast_p <- function(contrasts_df, cell_a, cell_b) {
    match_rows <- (grepl(cell_a, contrasts_df$contrast, fixed = TRUE) &
                   grepl(cell_b, contrasts_df$contrast, fixed = TRUE))
    row <- contrasts_df[match_rows, ]
    if (nrow(row) == 0) {
      cat(sprintf("  WARNING: contrast '%s' vs '%s' not found\n", cell_a, cell_b))
      return(NA_real_)
    }
    row$p.value[1]
  }

  p_same_loc <- get_contrast_p(all_contrasts,
                                "Same locality Same species",
                                "Same locality Diff. species")
  p_diff_loc <- get_contrast_p(all_contrasts,
                                "Diff. locality Same species",
                                "Diff. locality Diff. species")

  cat(sprintf("\nBracket p-values:\n"))
  cat(sprintf("  Intra vs Inter (same locality): %s\n", fmt_p(p_same_loc)))
  cat(sprintf("  Intra vs Inter (diff locality): %s\n", fmt_p(p_diff_loc)))

  # ---- Interaction p ----
  int_row     <- rownames(anova_tab)[grep(":", rownames(anova_tab))]
  int_p       <- anova_tab[int_row, "Pr(>F)"]
  int_p_label <- ifelse(int_p < 0.001,
                        "Interaction p < 0.001",
                        sprintf("Interaction p = %.3f", int_p))

  # ---- Save model outputs ----
  sink(file.path(output_dir, paste0("lmm_summary_", dist_label_short, ".txt")))
  cat("=== MODEL FORMULA ===\n");            print(formula(mod))
  cat("\n=== SUMMARY ===\n");                print(summary(mod))
  cat("\n=== ANOVA ===\n");                  print(anova_tab)
  cat("\n=== EMMEANS ===\n");                print(emm)
  cat("\n=== ALL PAIRWISE CONTRASTS ===\n"); print(all_contrasts)
  cat("\n=== INTERACTION CONTRASTS ===\n");  print(int_contrast)
  cat("\n=== SIMPLE EFFECTS ===\n");         print(pairs(simple_loc))
  sink()

  write.csv(emm_df,
            file.path(output_dir,
                      paste0("emmeans_2x2_", dist_label_short, ".csv")),
            row.names = FALSE)
  write.csv(all_contrasts,
            file.path(output_dir,
                      paste0("pairwise_contrasts_", dist_label_short, ".csv")),
            row.names = FALSE)
  write.csv(pairs %>%
              select(ind1, ind2, sp1, sp2, loc_hab1, loc_hab2,
                     div1, div2, same_species, same_locality,
                     group, distance),
            file.path(output_dir,
                      paste0("all_pairs_", dist_label_short, ".csv")),
            row.names = FALSE)

  # ---- Summary stats per group: mean, CI, SD ----
  bar_stats <- pairs %>%
    group_by(group) %>%
    summarise(
      mean_dist = mean(distance),
      sd_dist   = sd(distance),
      se        = sd_dist / sqrt(n()),
      n         = n(),
      ci_lo     = mean_dist - 1.96 * se,
      ci_hi     = mean_dist + 1.96 * se,
      sd_lo     = mean_dist - sd_dist,
      sd_hi     = mean_dist + sd_dist,
      .groups   = "drop"
    ) %>%
    mutate(n_label = paste0("n=", formatC(n, format = "d", big.mark = ",")))

  # ---- Panel A1: mean +/- 95% CI ----
  bar_ci <- bar_stats %>%
    mutate(err_lo = ci_lo, err_hi = ci_hi)

  panel_A1 <- make_barplot(bar_ci, dist_label, int_p_label,
                            p_same_loc, p_diff_loc,
                            error_label   = "mean \u00b1 95% CI",
                            group_colours = group_colours)

  # ---- Panel A2: mean +/- SD ----
  bar_sd <- bar_stats %>%
    mutate(err_lo = pmax(sd_lo, 0),
           err_hi = pmin(sd_hi, 1))

  panel_A2 <- make_barplot(bar_sd, dist_label, int_p_label,
                            p_same_loc, p_diff_loc,
                            error_label   = "mean \u00b1 SD",
                            group_colours = group_colours)

  # ---- Panel B: Interaction plot (emmeans) ----
  line_colours <- c("Same species"  = "#1B4F8A",
                    "Diff. species" = "#922B21")
  point_shapes <- c("Same species"  = 16,
                    "Diff. species" = 17)

  panel_B <- ggplot(emm_df,
                    aes(x      = locality_label,
                        y      = emmean,
                        group  = species_label,
                        colour = species_label,
                        shape  = species_label)) +
    geom_line(linewidth = 1.1) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL),
                  width = 0.08, linewidth = 0.7) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = line_colours, name = NULL) +
    scale_shape_manual(values  = point_shapes,  name = NULL) +
    labs(x = NULL,
         y = paste0("Estimated marginal mean\n", dist_label, " dissimilarity")) +
    annotate("text",
             x = 1.5, y = max(emm_df$UCL),
             label    = int_p_label,
             hjust    = 0.5, vjust = -0.5,
             size     = 3.2, fontface = "bold", colour = "grey20") +
    theme_classic(base_size = 11) +
    theme(
      legend.position  = "top",
      legend.text      = element_text(size = 9),
      axis.text.x      = element_text(size = 10, colour = "grey20"),
      axis.title.y     = element_text(size = 10),
      plot.margin      = margin(8, 8, 8, 8)
    )

  # ---- Save individual panels ----
  ggsave(file.path(output_dir,
                   paste0("panel_A1_CI_",  dist_label_short, ".png")),
         panel_A1, width = 7, height = 5.5, dpi = 300)
  ggsave(file.path(output_dir,
                   paste0("panel_A2_SD_",  dist_label_short, ".png")),
         panel_A2, width = 7, height = 5.5, dpi = 300)
  ggsave(file.path(output_dir,
                   paste0("panel_B_interaction_", dist_label_short, ".png")),
         panel_B, width = 5, height = 5, dpi = 300)

  # ---- Species-pair summary ----
  sp_pair_summary <- pairs %>%
    filter(!same_species, same_locality) %>%
    group_by(sp1, sp2, div1, div2) %>%
    summarise(
      n_pairs      = n(),
      mean_dist    = round(mean(distance), 4),
      mean_overlap = round(1 - mean(distance), 4),
      .groups      = "drop"
    ) %>%
    arrange(mean_dist)

  write.csv(sp_pair_summary,
            file.path(output_dir,
                      paste0("species_pairs_", dist_label_short, ".csv")),
            row.names = FALSE)

  cat(sprintf("\nGroup means (%s):\n", dist_label))
  pairs %>%
    group_by(group) %>%
    summarise(mean_dist = round(mean(distance), 4),
              sd_dist   = round(sd(distance), 4),
              n = n(), .groups = "drop") %>%
    print()

  cat(sprintf("\nLMM interaction (%s): p = %.4f\n", int_row, int_p))

  return(list(
    mod         = mod,
    anova_tab   = anova_tab,
    emm_df      = emm_df,
    int_p       = int_p,
    int_p_label = int_p_label,
    p_same_loc  = p_same_loc,
    p_diff_loc  = p_diff_loc,
    panel_A1    = panel_A1,
    panel_A2    = panel_A2,
    panel_B     = panel_B,
    bar_stats   = bar_stats,
    pairs       = pairs
  ))
}

# ==============================================================================
# 9. RUN FOR BOTH METRICS
# ==============================================================================

res_bc <- run_metric_analysis(bc_mat,  "Bray-Curtis",      "BC",
                               pairs_base, output_dir, group_colours)
res_wu <- run_metric_analysis(wu_mat,  "Weighted UniFrac", "WUniFrac",
                               pairs_base, output_dir, group_colours)

# ==============================================================================
# 10. COMBINED FIGURES
# combined_CI: all panels using 95% CI error bars
# combined_SD: all panels using SD error bars
# ==============================================================================

make_combined <- function(res_bc, res_wu, panel_type, ci_label, output_dir) {

  panel_bc <- if (panel_type == "CI") res_bc$panel_A1 else res_bc$panel_A2
  panel_wu <- if (panel_type == "CI") res_wu$panel_A1 else res_wu$panel_A2

  combined <- (panel_bc | res_bc$panel_B) /
              (panel_wu | res_wu$panel_B) +
    plot_annotation(
      title    = "Individual-level dietary co-occurrence analysis",
      subtitle = paste0(
        "Top row: Bray-Curtis | Bottom row: Weighted UniFrac\n",
        "Left panels: mean +/- ", ci_label,
        " with brackets for key contrasts. ",
        "Right panels: LMM estimated marginal means +/- 95% CI.\n",
        "Colours: dark = same locality (co-occurring); ",
        "light = different locality. Blue = intraspecific; Red = interspecific."
      ),
      caption  = paste0(
        "Brackets show Tukey-adjusted p-values for: ",
        "(1) Intraspecific vs Interspecific within same locality, ",
        "(2) Intraspecific vs Interspecific within different localities. ",
        "Interaction p shown in corner of left panels and above right panels. ",
        "LMM includes random effects for sp1, sp2, and locality."
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9,  colour = "grey30"),
        plot.caption  = element_text(size = 8,  colour = "grey40",
                                     hjust = 0, margin = margin(t = 6))
      )
    )

  fname <- file.path(output_dir,
                     paste0("cooccurrence_combined_", panel_type, ".png"))
  ggsave(fname, combined, width = 12, height = 10, dpi = 300)
  cat(sprintf("Saved: %s\n", fname))
}

make_combined(res_bc, res_wu, "CI", "95% CI", output_dir)
make_combined(res_bc, res_wu, "SD", "SD",     output_dir)

# ==============================================================================
# 11. FINAL CONSOLE SUMMARY
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

for (res in list(list(r = res_bc, label = "Bray-Curtis"),
                 list(r = res_wu, label = "Weighted UniFrac"))) {
  cat(sprintf("%-20s | interaction p = %.4f\n", res$label, res$r$int_p))
  cat(sprintf("  Intra vs Inter (same loc): %s\n", fmt_p(res$r$p_same_loc)))
  cat(sprintf("  Intra vs Inter (diff loc): %s\n", fmt_p(res$r$p_diff_loc)))
}

cat("\nAll outputs saved to:", output_dir, "\n")
cat("Done.\n")
