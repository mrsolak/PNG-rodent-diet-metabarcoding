################################################################################
# 02_sample_distribution.R
# Species altitude range and habitat distribution figure (Figure 2)
#
# Papua New Guinea rodent 18S diet metabarcoding
# Solak et al. — Molecular Ecology
#
# NOTE: If any parameter or code differs from the main Rmd script
# (PNG_18S_final_script_17_06_2026.Rmd), treat the Rmd as ground truth.
#
# Produces a faceted plot showing each species' altitudinal range (black bar)
# and individual sample points coloured by habitat type, faceted by mountain
# range (Mt. Wilhelm / Finisterres). Background rectangles show taxonomic
# divisions. Corresponds to Figure 2 of the manuscript.
#
# Input:  data/phyloseqs/phyloseq_final.R  (load with load())
# Output: Fig2_sample_distribution.png in output_dir
################################################################################

library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)

# ==============================================================================
# PATHS — edit before running
# ==============================================================================

data_dir     <- "/path/to/project/data"
phyloseq_dir <- file.path(data_dir, "phyloseqs")
output_dir   <- "/path/to/output/sample_distribution"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. LOAD PHYLOSEQ
# ==============================================================================

load(file.path(phyloseq_dir, "phyloseq_final.R"))
PHYLOSEQ_NEWMT <- PNG_18S_rodent_filtered_ensembles_tree_NEWMT_W_clean

sample_md <- as(sample_data(PHYLOSEQ_NEWMT), "data.frame")

# ==============================================================================
# 2. SPECIES ORDER AND COLOUR DEFINITIONS
# ==============================================================================

species_order_lookup_reversed <- c(
  "Rattus rattus", "Rattus exulans", "Rattus niobe", "Rattus steini", "Rattus verecundus",
  "Anisomys imitatorFR", "Anisomys imitatorMW", "Lorentzimys lowFR", "Lorentzimys lowMW",
  "Lorentzimys midMW", "Lorentzimys highFR", "Lorentzimys highMW", "Macruromys major",
  "Hyomys goliath", "Pogonomys loriae", "Pogonomys cf_sylvestris", "Pogonomys spNorthEastFR",
  "Pogonomys spNorthEastMW", "Coccymys shawmayeri", "Abeomelomys spMW",
  "Mallomys istapantap", "Mallomys rothschildi", "Mallomys hercules", "Mallomys cf_aroaensis",
  "Hydromys sp", "Leptomys ernstmayri", "Microhydromys argenteus", "Microhydromys richardsoni",
  "Pseudohydromys ellermani", "Pseudohydromys fuscus", "Uromys anak", "Uromys sp1",
  "Uromys cf_caudimaculatus", "Protochromys fellowsii", "Melomys low", "Melomys midFR",
  "Paramelomys platyopsMW", "Paramelomys platyopsWG", "Paramelomys platyopsBAI", "Paramelomys platyopsFR",
  "Paramelomys rubex_MidFR", "Paramelomys rubex_MidMW", "Paramelomys rubex_HighFR", "Paramelomys rubex_HighMW"
)

habitat_colors <- c(
  "forest"       = "darkgreen",
  "grassland"    = "yellow",
  "ecotone"      = "darkblue",
  "synanthropic" = "darkred"
)

division_colors_palette <- c(
  "Coccymys division"  = "#E60000",
  "Hydromys division"  = "#004D99",
  "Mallomys division"  = "#33CC33",
  "Pogonomys division" = "#FFCC00",
  "Rattus endemic"     = "#660066",
  "Rattus invasive"    = "#FF6600",
  "Uromys division"    = "#00CCCC"
)
custom_division_order <- names(division_colors_palette)

# ==============================================================================
# 3. PREPARE DATA
# ==============================================================================

species_divisions <- sample_md %>%
  group_by(Species_consensus_ddRAD_mt_others) %>%
  summarize(Division = first(Hydromyini.Division.Rattini.Group)) %>%
  ungroup() %>%
  mutate(
    Division = case_when(
      Division %in% c("incertae sedis group", "R. rattus group ") ~ "Rattus invasive",
      Division == "R. leucopus group (New Guinean group) "         ~ "Rattus endemic",
      TRUE ~ as.character(Division)
    )
  ) %>%
  select(Species_consensus_ddRAD_mt_others, Division) %>%
  distinct()

plot_data_final <- sample_md %>%
  filter(!is.na(Species_consensus_ddRAD_mt_others),
         !is.na(Alt),
         !is.na(Metalokalita),
         !is.na(Habitat_polished_Daniel.Fanda)) %>%
  left_join(species_divisions, by = "Species_consensus_ddRAD_mt_others") %>%
  mutate(
    Species       = factor(Species_consensus_ddRAD_mt_others,
                           levels = species_order_lookup_reversed),
    Mountain_Range = factor(Metalokalita),
    Habitat        = factor(Habitat_polished_Daniel.Fanda),
    Division       = factor(Division, levels = custom_division_order)
  ) %>%
  filter(Species %in% species_order_lookup_reversed)

# Per-species overall altitudinal range (for range bar)
overall_range_data <- plot_data_final %>%
  group_by(Species, Mountain_Range) %>%
  summarise(Min_Alt = min(Alt), Max_Alt = max(Alt), .groups = "drop")

# Division background rectangles
background_rect_data <- plot_data_final %>%
  select(Species, Division) %>%
  distinct() %>%
  group_by(Division) %>%
  summarise(
    xmin = min(as.numeric(Species)) - 0.5,
    xmax = max(as.numeric(Species)) + 0.5
  ) %>%
  ungroup() %>%
  tidyr::crossing(Mountain_Range = unique(plot_data_final$Mountain_Range)) %>%
  mutate(fill_color = division_colors_palette[as.character(Division)])

custom_y_breaks <- c(0, 700, 1200, 1700, 2200, 2700, 3200, 3700, 4200, 4500)
vline_positions <- seq(from = 1.5,
                       to   = length(species_order_lookup_reversed) - 0.5,
                       by   = 1)

# ==============================================================================
# 4. PLOT
# ==============================================================================

p_habitat_range <- ggplot(plot_data_final, aes(x = Species, y = Alt)) +

  geom_rect(data = background_rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Division),
            inherit.aes = FALSE, alpha = 0.2) +

  scale_fill_manual(values = division_colors_palette, name = "Taxonomic Division") +

  geom_segment(data = overall_range_data,
               aes(x = Species, xend = Species, y = Min_Alt, yend = Max_Alt),
               color = "black", size = 4, inherit.aes = FALSE) +

  geom_point(aes(color = Habitat),
             position = position_jitter(width = 0.2, height = 0),
             alpha = 0.8, shape = 15, size = 4) +

  geom_vline(xintercept = vline_positions,
             linetype = "dashed", color = "gray60", size = 0.3) +

  scale_y_continuous(breaks = custom_y_breaks) +
  scale_color_manual(values = habitat_colors, name = "Habitat") +

  facet_grid(Mountain_Range ~ ., scales = "free_y", space = "free_y") +

  theme_minimal(base_size = 12) +
  labs(x = "Species", y = "Elevation [masl]",
       title = "Species Altitude Range and Habitat Occurrence") +
  theme(
    axis.text.x       = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      face = "italic", size = 8, color = "black"),
    axis.text.y       = element_text(size = 8),
    strip.text        = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
    panel.spacing.y   = unit(1, "lines"),
    axis.line.x       = element_blank(),
    legend.position   = "right"
  ) +
  guides(
    fill  = "none",
    color = guide_legend(override.aes = list(size = 10))
  )

print(p_habitat_range)

ggsave(file.path(output_dir, "Fig2_sample_distribution.png"),
       p_habitat_range, width = 20, height = 10, dpi = 300)

cat("Figure 2 saved to:", output_dir, "\n")
