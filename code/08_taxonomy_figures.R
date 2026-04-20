### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Taxonomy figures ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

## =========================================================
## 1. Load saved tables
## =========================================================

# Heatmap tables
bac_phylum_heat <- read.csv("./results/tables/taxonomy_heatmap_bacteria_phylum.csv", header = TRUE, check.names = FALSE)
fun_phylum_heat <- read.csv("./results/tables/taxonomy_heatmap_fungi_phylum.csv", header = TRUE, check.names = FALSE)
bac_genus_heat  <- read.csv("./results/tables/taxonomy_heatmap_bacteria_genus.csv", header = TRUE, check.names = FALSE)
fun_genus_heat  <- read.csv("./results/tables/taxonomy_heatmap_fungi_genus.csv", header = TRUE, check.names = FALSE)

# Barplot tables
bac_phylum_bar <- read.csv("./results/tables/taxonomy_barplot_bacteria_phylum.csv", header = TRUE, check.names = FALSE)
fun_phylum_bar <- read.csv("./results/tables/taxonomy_barplot_fungi_phylum.csv", header = TRUE, check.names = FALSE)
bac_genus_bar  <- read.csv("./results/tables/taxonomy_barplot_bacteria_genus.csv", header = TRUE, check.names = FALSE)
fun_genus_bar  <- read.csv("./results/tables/taxonomy_barplot_fungi_genus.csv", header = TRUE, check.names = FALSE)

# Model results
taxa_results <- read.csv("./results/tables/taxonomy_top5_taxa_beta_models.csv", header = TRUE, check.names = FALSE)

## =========================================================
## 2. Restore factor levels
## =========================================================

effect_levels <- c("Intensity", "Frequency", "Composition")
level_order <- c("Dry", "Wet", "Low", "High", "Conspecific", "Heterospecific")

restore_heatmap_factors <- function(df, tax_rank) {
  df$Effect <- factor(df$Effect, levels = effect_levels)
  df$Level  <- factor(df$Level, levels = level_order)
  
  tax_levels <- unique(df[[tax_rank]])
  df[[tax_rank]] <- factor(df[[tax_rank]], levels = tax_levels)
  
  df
}

restore_barplot_factors <- function(df) {
  df$Effect <- factor(df$Effect, levels = effect_levels)
  df$Level  <- factor(df$Level, levels = level_order)
  
  tax_levels <- unique(df$Taxon)
  df$Taxon <- factor(df$Taxon, levels = tax_levels)
  
  df
}

bac_phylum_heat <- restore_heatmap_factors(bac_phylum_heat, "Phylum")
fun_phylum_heat <- restore_heatmap_factors(fun_phylum_heat, "Phylum")
bac_genus_heat  <- restore_heatmap_factors(bac_genus_heat, "Genus")
fun_genus_heat  <- restore_heatmap_factors(fun_genus_heat, "Genus")

bac_phylum_bar <- restore_barplot_factors(bac_phylum_bar)
fun_phylum_bar <- restore_barplot_factors(fun_phylum_bar)
bac_genus_bar  <- restore_barplot_factors(bac_genus_bar)
fun_genus_bar  <- restore_barplot_factors(fun_genus_bar)

## =========================================================
## 3. Helper functions for heatmaps
## =========================================================

p_to_stars <- function(p) {
  dplyr::case_when(is.na(p)  ~ "", p < 0.001 ~ "***", p < 0.01  ~ "**",
                   p < 0.05  ~ "*", TRUE      ~ "")
}

make_sig_table <- function(results_df, taxa_vec, rank_name, tax_levels) {
  results_df %>%
    filter(Rank == rank_name, Taxon %in% taxa_vec,
           Term %in% c("Intensity", "Frequency", "Composition")) %>%
    mutate(Effect = factor(Term, levels = effect_levels),
           Taxon = factor(Taxon, levels = tax_levels),
           star = p_to_stars(p_adj), x_pos = 1.5) %>%
    dplyr::select(Taxon, Effect, x_pos, star)
}

make_highlight_table <- function(results_df, heat_df, taxa_vec, rank_name, tax_rank, tax_levels) {
  
  sig_main <- results_df %>%
    filter(Rank == rank_name, Taxon %in% taxa_vec,
           Term %in% c("Intensity", "Frequency", "Composition"),
           p_adj < 0.05) %>%
    transmute(Taxon, Effect = Term)
  
  highlight_df <- bind_rows(
    heat_df %>%
      filter(.data[[tax_rank]] %in% taxa_vec, Effect == "Intensity",
             Level %in% c("Dry", "Wet")) %>%
      group_by(Taxon = .data[[tax_rank]], Effect) %>%
      slice_max(order_by = mean_abundance, n = 1, with_ties = FALSE) %>%
      ungroup(),
    
    heat_df %>%
      filter(.data[[tax_rank]] %in% taxa_vec, Effect == "Frequency",
             Level %in% c("Low", "High")) %>%
      group_by(Taxon = .data[[tax_rank]], Effect) %>%
      slice_max(order_by = mean_abundance, n = 1, with_ties = FALSE) %>%
      ungroup(),
    
    heat_df %>%
      filter(.data[[tax_rank]] %in% taxa_vec, Effect == "Composition",
             Level %in% c("Conspecific", "Heterospecific")) %>%
      group_by(Taxon = .data[[tax_rank]], Effect) %>%
      slice_max(order_by = mean_abundance, n = 1, with_ties = FALSE) %>%
      ungroup()) %>%
    mutate(Effect = factor(Effect, levels = effect_levels),
           Taxon = factor(Taxon, levels = tax_levels)) %>%
    semi_join(sig_main, by = c("Taxon", "Effect"))
  
  highlight_df
}

plot_tax_heatmap <- function(heat_df, sig_df, highlight_df, tax_rank, title_text, colours) {
  ggplot(heat_df, aes(x = Level, y = .data[[tax_rank]], fill = log10(mean_abundance + 1e-6))) +
    geom_tile(color = "white", linewidth = 0.4, width = 0.75) +
    geom_tile(data = highlight_df, aes(x = Level, y = Taxon), inherit.aes = FALSE,
              fill = NA, color = "black", linewidth = 1.1, width = 0.75, height = 0.9) +
    geom_text(data = sig_df, aes(x = x_pos, y = Taxon, label = star), inherit.aes = FALSE,
              fontface = "bold", size = 5) +
    facet_wrap(~ Effect, scales = "free_x", nrow = 1) +
    scale_fill_gradientn(colours = colours,
                         values = scales::rescale(c(-6, -4, -3, -2, -1)),
                         name = "Mean relative abundance",
                         breaks = log10(c(0.1, 0.01, 0.0001)),
                         labels = c("0.1", "0.01", "0.0001")) +
    labs(x = NULL, y = tax_rank, title = title_text) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11, color = "black"),
          axis.text.y = element_text(face = ifelse(tax_rank == "Genus", "bold.italic", "bold"),
                                     size = 9, color = "black"),
          axis.title.y = element_text(face = "bold", size = 11),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 9))
}

## =========================================================
## 4. Build significance/highlight tables for heatmaps
## =========================================================

top5_bac_phyla  <- levels(bac_phylum_heat$Phylum)
top5_fun_phyla  <- levels(fun_phylum_heat$Phylum)
top5_bac_genera <- levels(bac_genus_heat$Genus)
top5_fun_genera <- levels(fun_genus_heat$Genus)

bac_phylum_sig <- make_sig_table(results_df = taxa_results,
                                 taxa_vec = top5_bac_phyla,
                                 rank_name = "Bacterial phylum",
                                 tax_levels = levels(bac_phylum_heat$Phylum))

fun_phylum_sig <- make_sig_table(results_df = taxa_results,
                                 taxa_vec = top5_fun_phyla,
                                 rank_name = "Fungal phylum",
                                 tax_levels = levels(fun_phylum_heat$Phylum))

bac_genus_sig <- make_sig_table(results_df = taxa_results,
                                taxa_vec = top5_bac_genera,
                                rank_name = "Bacterial genus",
                                tax_levels = levels(bac_genus_heat$Genus))

fun_genus_sig <- make_sig_table(results_df = taxa_results,
                                taxa_vec = top5_fun_genera,
                                rank_name = "Fungal genus",
                                tax_levels = levels(fun_genus_heat$Genus))

bac_phylum_highlight <- make_highlight_table(results_df = taxa_results,
                                             heat_df = bac_phylum_heat,
                                             taxa_vec = top5_bac_phyla,
                                             rank_name = "Bacterial phylum",
                                             tax_rank = "Phylum",
                                             tax_levels = levels(bac_phylum_heat$Phylum))

fun_phylum_highlight <- make_highlight_table(results_df = taxa_results,
                                             heat_df = fun_phylum_heat,
                                             taxa_vec = top5_fun_phyla,
                                             rank_name = "Fungal phylum",
                                             tax_rank = "Phylum",
                                             tax_levels = levels(fun_phylum_heat$Phylum))

bac_genus_highlight <- make_highlight_table(results_df = taxa_results,
                                            heat_df = bac_genus_heat,
                                            taxa_vec = top5_bac_genera,
                                            rank_name = "Bacterial genus",
                                            tax_rank = "Genus",
                                            tax_levels = levels(bac_genus_heat$Genus))

fun_genus_highlight <- make_highlight_table(results_df = taxa_results,
                                            heat_df = fun_genus_heat,
                                            taxa_vec = top5_fun_genera,
                                            rank_name = "Fungal genus",
                                            tax_rank = "Genus",
                                            tax_levels = levels(fun_genus_heat$Genus))

## =========================================================
## 5. Heatmap color palettes
## =========================================================

bac_cols_phy <- c("white", "#E3F2FD", "#90CAF9", "#42A5F5", "#1E88E5", "#0D47A1")
bac_cols_gen <- c("white", "#E0F7FA", "#80DEEA", "#26C6DA", "#0097A7", "#006064")
fun_cols_phy <- c("white", "#F3E5F5", "#CE93D8", "#AB47BC", "#8E24AA", "#4A148C")
fun_cols_gen <- c("white", "#F8BBD0", "#F06292", "#EC407A", "#C2185B", "#880E4F")

## =========================================================
## 6. Build heatmaps
## =========================================================

p_bac_phylum_heat <- plot_tax_heatmap(heat_df = bac_phylum_heat,
                                      sig_df = bac_phylum_sig,
                                      highlight_df = bac_phylum_highlight,
                                      tax_rank = "Phylum",
                                      title_text = "Bacteria - Phylum",
                                      colours = bac_cols_phy)

p_fun_phylum_heat <- plot_tax_heatmap(heat_df = fun_phylum_heat,
                                      sig_df = fun_phylum_sig,
                                      highlight_df = fun_phylum_highlight,
                                      tax_rank = "Phylum",
                                      title_text = "Fungi - Phylum",
                                      colours = fun_cols_phy)

p_bac_genus_heat <- plot_tax_heatmap(heat_df = bac_genus_heat,
                                     sig_df = bac_genus_sig,
                                     highlight_df = bac_genus_highlight,
                                     tax_rank = "Genus",
                                     title_text = "Bacteria - Genus",
                                     colours = bac_cols_gen)

p_fun_genus_heat <- plot_tax_heatmap(heat_df = fun_genus_heat,
                                     sig_df = fun_genus_sig,
                                     highlight_df = fun_genus_highlight,
                                     tax_rank = "Genus",
                                     title_text = "Fungi - Genus",
                                     colours = fun_cols_gen)

p_taxa_heatmaps <- ggarrange(p_bac_phylum_heat, p_fun_phylum_heat,
                             p_bac_genus_heat,  p_fun_genus_heat,
                             labels = c("(a)", "(b)", "(c)", "(d)"),
                             nrow = 2, ncol = 2)
p_taxa_heatmaps

## =========================================================
## 7. Helper functions for barplots
## =========================================================

bolditalic_taxon_labels <- function(labels) {
  lapply(labels, function(x) {
    if (x == "Other") {
      as.expression(bquote(bold(.(x))))
    } else {
      as.expression(bquote(bolditalic(.(x))))
    }
  })
}

plot_tax_bar <- function(df, title_text, palette, italic_legend = FALSE) {
  
  legend_labs <- if (italic_legend) {
    bolditalic_taxon_labels(levels(df$Taxon))
  } else {
    levels(df$Taxon)
  }
  
  ggplot(df, aes(x = Level, y = mean_abundance, fill = Taxon)) +
    geom_col(width = 0.7) +
    facet_wrap(~ Effect, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = palette, labels = legend_labs) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Mean relative abundance", title = title_text, fill = NULL) +
    guides(fill = guide_legend(nrow = 3)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11, color = "black"),
          axis.title.y = element_text(face = "bold", size = 10),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size = 9, color = "black"))
}

## =========================================================
## 8. Barplot color palettes
## =========================================================

bac_phy_cols <- brewer.pal(12, "Paired")[1:11]
bac_phy_cols[11] <- "grey80"

fun_phy_cols <- brewer.pal(5, "Set2")
fun_phy_cols <- c(fun_phy_cols, "grey80")

bac_gen_cols <- brewer.pal(11, "Set3")
bac_gen_cols[11] <- "grey80"

fun_gen_cols <- brewer.pal(11, "Spectral")
fun_gen_cols[11] <- "grey80"

# Ensure names match factor levels
names(bac_phy_cols) <- levels(bac_phylum_bar$Taxon)
names(fun_phy_cols) <- levels(fun_phylum_bar$Taxon)
names(bac_gen_cols) <- levels(bac_genus_bar$Taxon)
names(fun_gen_cols) <- levels(fun_genus_bar$Taxon)

## =========================================================
## 9. Build barplots
## =========================================================

p_bac_phylum_bar <- plot_tax_bar(df = bac_phylum_bar,
                                 title_text = "Bacteria - Phylum",
                                 palette = bac_phy_cols,
                                 italic_legend = FALSE)

p_fun_phylum_bar <- plot_tax_bar(df = fun_phylum_bar,
                                 title_text = "Fungi - Phylum",
                                 palette = fun_phy_cols,
                                 italic_legend = FALSE)

p_bac_genus_bar <- plot_tax_bar(df = bac_genus_bar,
                                title_text = "Bacteria - Genus",
                                palette = bac_gen_cols,
                                italic_legend = TRUE)

p_fun_genus_bar <- plot_tax_bar(df = fun_genus_bar,
                                title_text = "Fungi - Genus",
                                palette = fun_gen_cols,
                                italic_legend = TRUE)

p_taxa_bars <- ggarrange(p_bac_phylum_bar, p_fun_phylum_bar,
                         p_bac_genus_bar,  p_fun_genus_bar,
                         labels = c("(a)", "(b)", "(c)", "(d)"),
                         nrow = 2, ncol = 2)
p_taxa_bars
## =========================================================
## 10. Save figures
## =========================================================

# - Figure 3
ggsave(filename = "./results/figures/taxonomy_heatmaps_top5.pdf", plot = p_taxa_heatmaps,
       width = 16, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = "./results/figures/png/taxonomy_heatmaps_top5.png", plot = p_taxa_heatmaps,
       width = 16, height = 10, units = "in", dpi = 300)

# - Figure S4
ggsave(filename = "./results/figures/taxonomy_barplots_top10.pdf", plot = p_taxa_bars,
       width = 14, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = "./results/figures/png/taxonomy_barplots_top10.png", plot = p_taxa_bars,
       width = 14, height = 10, units = "in", dpi = 300)
