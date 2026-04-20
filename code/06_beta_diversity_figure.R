### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Beta-diversity figure ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(grid)
library(tibble)

## =========================================================
## 1. Helper functions
## =========================================================

p_to_stars <- function(p) {
  dplyr::case_when(is.na(p)  ~ "", p < 0.001 ~ "***", p < 0.01  ~ "**",
                   p < 0.05  ~ "*", TRUE      ~ "")
}

build_permanova_label <- function(perm_tbl, term_map = NULL) {
  perm_tbl <- as.data.frame(perm_tbl)
  
  if (!all(c("Term", "Pr(>F)") %in% names(perm_tbl))) {
    stop("PERMANOVA table must contain columns 'Term' and 'Pr(>F)'.")
  }
  
  perm_tbl <- perm_tbl %>%
    filter(!Term %in% c("Residual", "Total"), !is.na(`Pr(>F)`))
  
  if (!is.null(term_map)) {
    perm_tbl$Term_label <- dplyr::recode(perm_tbl$Term, !!!term_map)
  } else {
    perm_tbl$Term_label <- perm_tbl$Term
  }
  
  lines <- paste0(perm_tbl$Term_label, " ", p_to_stars(perm_tbl$`Pr(>F)`))
  paste(c("PERMANOVA", lines), collapse = "\n")
}

make_pcoa_plot <- function(ord_df, pcoa_var_tbl, dataset_name, permanova_label,
                           flip_x_annot = NULL, flip_y_annot = NULL) {
  ord_df <- ord_df %>%
    mutate(Intensity   = factor(Intensity, levels = c("Dry", "Wet")),
           Frequency   = factor(Frequency, levels = c("Low", "High")),
           Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))
  
  var_df <- pcoa_var_tbl %>%
    filter(Dataset == dataset_name)
  
  pc1 <- round(var_df$Variance_explained_percent[var_df$Axis == "PCoA1"], 1)
  pc2 <- round(var_df$Variance_explained_percent[var_df$Axis == "PCoA2"], 1)
  
  int_cols <- c(Dry = "#E7782F", Wet = "#1994E5")
  
  if (is.null(flip_x_annot)) {
    flip_x_annot <- min(ord_df$PCoA1, na.rm = TRUE)
  }
  if (is.null(flip_y_annot)) {
    flip_y_annot <- 0
  }
  
  ggplot(ord_df, aes(PCoA1, PCoA2)) +
    geom_point(aes(color = Intensity, fill = Intensity, shape = Frequency, alpha = Composition),
               size = 3, stroke = 0.8) +
    scale_color_manual(values = int_cols, name = "Intensity") +
    scale_fill_manual(values = int_cols, guide = "none") +
    scale_shape_manual(values = c(High = 21, Low = 24), name = "Frequency") +
    scale_alpha_manual(values = c(Conspecific = 1, Heterospecific = 0.3), name = "Composition") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    annotate("text", x = flip_x_annot, y = flip_y_annot, label = permanova_label,
             hjust = 0, vjust = 1.5, size = 3.5, fontface = "bold") +
    labs(x = paste0("PCoA1 (", pc1, "%)"), y = paste0("PCoA2 (", pc2, "%)"), title = dataset_name) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(face = "bold", size = 11),
          legend.title = element_text(face = "bold", size = 11),
          axis.line = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          axis.title = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5)) +
    guides(color = guide_legend(order = 1, override.aes = list(shape = 21, fill = int_cols, size = 3)),
           shape = guide_legend(order = 2, override.aes = list(fill = NA, color = "black")),
           alpha = guide_legend(order = 3, override.aes = list(shape = 21, color = "grey20", fill = "grey20")))
}

make_genera_panel <- function(header, taxa) {
  ggdraw() +
    draw_grob(grid::roundrectGrob(x = 0.5, y = 0.5, width = 0.85, height = 0.80,
                                  r = unit(0.03, "snpc"), gp = grid::gpar(col = "#888888", fill = NA, lwd = 1.2))) +
    draw_label(header, x = 0.5, y = 0.70, fontface = "bold", size = 10, hjust = 0.5) +
    draw_label(taxa, x = 0.5, y = 0.30, fontface = "italic", size = 9.5, hjust = 0.5, lineheight = 1.1)
}

make_donut_pair <- function(df, domain_label, phylum_colors, genera_low, genera_high,
                            left_label = "Low PCoA1", right_label = "High PCoA1") {
  plot_df <- df %>%
    mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum),
           PCoA1_end = factor(PCoA1_end, levels = c("Low_PCoA1", "High_PCoA1"),
                              labels = c(left_label, right_label))) %>%
    count(PCoA1_end, Phylum, name = "n") %>%
    group_by(PCoA1_end) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    filter(prop > 0)
  
  center_label_df <- plot_df %>%
    group_by(PCoA1_end) %>%
    summarise(asv_n = sum(n), .groups = "drop") %>%
    mutate(label = dplyr::case_when(PCoA1_end == left_label  ~ paste0(left_label, "\n(", domain_label, ")"),
                                    PCoA1_end == right_label ~ paste0(right_label, "\n(", domain_label, ")")),
           x = 1, y = 0.5)
  
  low_data <- plot_df %>%
    filter(PCoA1_end == left_label) %>%
    arrange(desc(Phylum)) %>%
    mutate(y_mid = cumsum(prop) - prop / 2)
  
  high_data <- plot_df %>%
    filter(PCoA1_end == right_label) %>%
    arrange(desc(Phylum)) %>%
    mutate(y_mid = cumsum(prop) - prop / 2)
  
  p_left <- ggplot(plot_df %>% filter(PCoA1_end == left_label),
                   aes(x = 2, y = prop, fill = Phylum)) +
    geom_col(width = 0.9, color = "white", linewidth = 0.5) +
    scale_fill_manual(values = phylum_colors, drop = FALSE) +
    geom_text(data = low_data %>% filter(prop > 0.05),
              aes(x = 2, y = y_mid, label = paste0(n, " ASVs")),
              inherit.aes = FALSE, color = "white", fontface = "bold", size = 3) +
    geom_text(data = center_label_df %>% filter(PCoA1_end == left_label),
              aes(x = x, y = y, label = label), inherit.aes = FALSE,
              fontface = "bold", size = 3.5, lineheight = 1.1) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "none")
  
  p_legend_src <- ggplot(plot_df, aes(x = 2, y = prop, fill = Phylum)) +
    geom_col(width = 0.9, color = "white", linewidth = 0.6) +
    scale_fill_manual(values = phylum_colors, name = "Phylum") +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 9),
          legend.key.size = unit(0.5, "cm"),
          legend.spacing.y = unit(0.15, "cm"))
  
  legend_only <- cowplot::get_legend(p_legend_src)
  
  p_right <- ggplot(plot_df %>% filter(PCoA1_end == right_label),
                    aes(x = 2, y = prop, fill = Phylum)) +
    geom_col(width = 0.9, color = "white", linewidth = 0.5) +
    scale_fill_manual(values = phylum_colors, drop = FALSE) +
    geom_text(data = high_data %>% filter(prop > 0.05),
              aes(x = 2, y = y_mid, label = paste0(n, " ASVs")),
              inherit.aes = FALSE, color = "white", fontface = "bold", size = 3) +
    geom_text(data = center_label_df %>% filter(PCoA1_end == right_label),
              aes(x = x, y = y, label = label), inherit.aes = FALSE,
              fontface = "bold", size = 3.5, lineheight = 1.1) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "none")
  
  p_genera_low <- make_genera_panel("Representative genera", genera_low)
  p_genera_high <- make_genera_panel(
    ifelse(grepl(",", genera_high), "Representative genera", "Representative genus"),
    genera_high)
  p_blank <- ggdraw()
  
  top_row <- plot_grid(p_left, legend_only, p_right, nrow = 1, rel_widths = c(1, 0.4, 1))
  bottom_row <- plot_grid(p_genera_low, p_blank, p_genera_high, nrow = 1, rel_widths = c(1, 0.4, 1))
  
  plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(3, 0.5))
}

## =========================================================
## 2. Load saved analysis outputs
## =========================================================

pcoa_bac <- read.csv("./results/tables/pcoa_coordinates_bacteria.csv", header = TRUE, check.names = FALSE)
pcoa_fun <- read.csv("./results/tables/pcoa_coordinates_fungi.csv", header = TRUE, check.names = FALSE)
pcoa_var <- read.csv("./results/tables/pcoa_variance_explained.csv", header = TRUE, check.names = FALSE)

perm_bac <- read.csv("./results/tables/permanova_bacteria_sample_level.csv", header = TRUE, check.names = FALSE)
perm_fun <- read.csv("./results/tables/permanova_fungi_sample_level.csv", header = TRUE, check.names = FALSE)

pcoa1_bac_asvs <- read.csv("./results/tables/pcoa1_associated_asvs_bacteria.csv", header = TRUE, check.names = FALSE)
pcoa1_fun_asvs <- read.csv("./results/tables/pcoa1_associated_asvs_fungi.csv", header = TRUE, check.names = FALSE)

## =========================================================
## 3. Build PERMANOVA labels
## =========================================================

term_map <- c("Intensity" = "I",
              "Frequency" = "F",
              "Composition" = "C",
              "Intensity:Frequency" = "I × F",
              "Intensity:Composition" = "I × C",
              "Frequency:Composition" = "F × C",
              "Intensity:Frequency:Composition" = "I × F × C")

perm_label_bac <- build_permanova_label(perm_bac, term_map = term_map)
perm_label_fun <- build_permanova_label(perm_fun, term_map = term_map)

## =========================================================
## 4. Main PCoA panels
## =========================================================

p_bac <- make_pcoa_plot(ord_df = pcoa_bac, pcoa_var_tbl = pcoa_var,
                        dataset_name = "Bacteria", permanova_label = perm_label_bac,
                        flip_x_annot = min(pcoa_bac$PCoA1, na.rm = TRUE), flip_y_annot = 0)

p_fun <- make_pcoa_plot(ord_df = pcoa_fun, pcoa_var_tbl = pcoa_var,
                        dataset_name = "Fungi", permanova_label = perm_label_fun,
                        flip_x_annot = min(pcoa_fun$PCoA1, na.rm = TRUE), flip_y_annot = 0)

p_pcoa <- ggarrange(p_bac, p_fun, nrow = 1, labels = c("(a)", "(b)"), common.legend = TRUE,
                    legend = "bottom")
p_pcoa

## =========================================================
## 5. Donut panels for PCoA1 ends
## =========================================================

# Bacterial phylum colors
all_phyla_bac <- sort(unique(ifelse(is.na(pcoa1_bac_asvs$Phylum) | pcoa1_bac_asvs$Phylum == "", "Unclassified", pcoa1_bac_asvs$Phylum)))
phylum_colors_bac <- c("#7BAF9E", "#C4977A", "#8FA8C8", "#C4A882", "#A98DB5",
                       "#B5C46A", "#D4847A", "#6AAFB5", "#C4B87A", "#9AB58F",
                       "#B57A8F", "#7A9BC4", "#C4C47A", "#8F7AB5", "#7AB5A0")
phylum_colors_bac <- rep(phylum_colors_bac, length.out = length(all_phyla_bac))
names(phylum_colors_bac) <- all_phyla_bac

bac_donut <- make_donut_pair(df = pcoa1_bac_asvs, domain_label = "Bacteria",
                             phylum_colors = phylum_colors_bac,
                             genera_low = "Streptomyces, Nocardioides",
                             genera_high = "BIrii41, Turneriella")

bac_donut <- ggarrange(NULL, bac_donut, nrow = 1, widths = c(0.1, 1))

# Fungal phylum colors
all_phyla_fun <- sort(unique(ifelse(is.na(pcoa1_fun_asvs$Phylum) | pcoa1_fun_asvs$Phylum == "", "Unclassified", pcoa1_fun_asvs$Phylum)))
phylum_colors_fun <- c("Ascomycota" = "#7A9BC4")
if (!all(all_phyla_fun %in% names(phylum_colors_fun))) {
  extra_phyla <- setdiff(all_phyla_fun, names(phylum_colors_fun))
  extra_cols <- rep("#BDBDBD", length(extra_phyla))
  names(extra_cols) <- extra_phyla
  phylum_colors_fun <- c(phylum_colors_fun, extra_cols)
}
phylum_colors_fun <- phylum_colors_fun[all_phyla_fun]

fun_donut <- make_donut_pair(df = pcoa1_fun_asvs, domain_label = "Fungi",
                             phylum_colors = phylum_colors_fun,
                             genera_low = "Darksidea, Dactylella",
                             genera_high = "Pseudomassariosphaeria")

fun_donut <- ggarrange(NULL, fun_donut, nrow = 1, widths = c(0.1, 1))

p_donut <- ggarrange(bac_donut, fun_donut, nrow = 1, labels = c("(c)", "(d)"))

## =========================================================
## 6. Combine and save final figure
## =========================================================

p_all <- ggarrange(p_pcoa, p_donut, nrow = 2, heights = c(1, 0.6))

# Save - Figure 2
ggsave(filename = "./results/figures/beta_diversity_figure.pdf",
       plot = p_all, width = 14, height = 10, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/beta_diversity_figure.png",
       plot = p_all, width = 14, height = 10, units = "in", dpi = 300, bg = "white")
