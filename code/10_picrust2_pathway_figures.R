### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: PICRUSt2 pathway figures ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)

## =========================================================
## 1. Load saved tables
## =========================================================

plot_cat_eff_I <- read.csv("./results/tables/picrust2_effect_plot_intensity.csv",
                           header = TRUE, check.names = FALSE)

plot_cat_eff_F <- read.csv("./results/tables/picrust2_effect_plot_frequency.csv",
                           header = TRUE, check.names = FALSE)

plot_cat_IF <- read.csv("./results/tables/picrust2_effect_plot_intensity_frequency.csv",
                        header = TRUE, check.names = FALSE)

significance_summary <- read.csv("./results/tables/picrust2_pathway_significance_summary.csv",
                                 header = TRUE, check.names = FALSE)

## =========================================================
## 2. Restore factor levels
## =========================================================

class_l1_levels <- c("Biosynthesis / cellular components",
                     "Substrate utilization / nutrient cycling",
                     "Energy metabolism",
                     "Cellular processes",
                     "Stress response",
                     "Other")

plot_cat_eff_I <- plot_cat_eff_I %>%
  mutate(Class_L1 = factor(Class_L1, levels = class_l1_levels),
         Direction = factor(Direction, levels = c("Dry", "Wet", "Neutral")),
         Label = factor(Label, levels = unique(Label)))

plot_cat_eff_F <- plot_cat_eff_F %>%
  mutate(Class_L1 = factor(Class_L1, levels = class_l1_levels),
         Direction = factor(Direction, levels = c("Low", "High", "Neutral")),
         Label = factor(Label, levels = unique(Label)))

plot_cat_IF <- plot_cat_IF %>%
  mutate(Class_L1 = factor(Class_L1, levels = class_l1_levels),
         Direction = factor(Direction, levels = c("Dry", "Wet", "Neutral")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Label = factor(Label, levels = unique(Label)))

## =========================================================
## 3. Titles from significance summary
## =========================================================

n_intensity <- significance_summary %>%
  filter(Term == "Intensity") %>%
  pull(n_significant_FDR_0.05)

n_frequency <- significance_summary %>%
  filter(Term == "Frequency") %>%
  pull(n_significant_FDR_0.05)

n_if <- significance_summary %>%
  filter(Term == "Intensity:Frequency") %>%
  pull(n_significant_FDR_0.05)

title_intensity <- paste0("Intensity (", n_intensity, " pathways)")
title_frequency <- paste0("Frequency (", n_frequency, " pathways)")
title_if <- paste0("Intensity × Frequency (", n_if, " pathways)")

## =========================================================
## 4. Panel A: Intensity
## =========================================================

p_path_intensity <- ggplot(plot_cat_eff_I, aes(x = mean_diff, y = Label, color = Direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = mean_diff - se_diff, xmax = mean_diff + se_diff),
                width = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Dry" = "#E7782F", "Wet" = "#1994E5", "Neutral" = "grey50")) +
  labs(x = "Mean change in pathway abundance (Dry − Wet, CLR)",
       y = "Functional category", color = NULL, title = title_intensity) +
  facet_grid(Class_L1 ~ ., scales = "free", space = "free", switch = "y") +
  theme_bw() +
  theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        strip.placement = "outside",
        strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.90),
        legend.justification = c(1, 1),
        legend.background = element_rect(color = "black", linewidth = 0.3),
        legend.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

## =========================================================
## 5. Panel B: Frequency
## =========================================================

p_path_frequency <- ggplot(plot_cat_eff_F, aes(x = mean_diff, y = Label, color = Direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = mean_diff - se_diff, xmax = mean_diff + se_diff),
                width = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Low" = "#E5D019", "High" = "#B38AF1", "Neutral" = "grey50")) +
  labs(x = "Mean change in pathway abundance (Low − High, CLR)",
       y = "Functional category", color = NULL, title = title_frequency) +
  facet_grid(Class_L1 ~ ., scales = "free", space = "free", switch = "y") +
  theme_bw() +
  theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        strip.placement = "outside",
        strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25, 0.90),
        legend.justification = c(1, 1),
        legend.background = element_rect(color = "black", linewidth = 0.3),
        legend.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

## =========================================================
## 6. Panel C: Intensity × Frequency
## =========================================================

bg_freq <- plot_cat_IF %>%
  distinct(Class_L1, Frequency) %>%
  mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
         fill_bg = ifelse(Frequency == "Low", "Low", "High"))

p_path_if <- ggplot(plot_cat_IF, aes(x = mean_diff, y = Label, color = Direction)) +
  geom_rect(data = bg_freq, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_bg),
            inherit.aes = FALSE, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = mean_diff - se_diff, xmax = mean_diff + se_diff),
                width = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Dry" = "#E7782F", "Wet" = "#1994E5", "Neutral" = "grey50")) +
  scale_fill_manual(values = c("Low" = "#E5D019", "High" = "#B38AF1"), guide = "none") +
  labs(x = "Mean change in pathway abundance (Dry − Wet, CLR)",
       y = "Functional category", color = NULL, title = title_if) +
  facet_grid(Class_L1 ~ Frequency, scales = "free_y", space = "free_y", switch = "y") +
  theme_bw() +
  theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        strip.placement = "outside",
        strip.background.y = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
        strip.text.x = element_text(face = "bold", size = 11),
        strip.background.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25, 0.90),
        legend.justification = c(1, 1),
        legend.background = element_rect(color = "black", linewidth = 0.3),
        legend.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

## =========================================================
## 7. Combine and save figure
## =========================================================

p_pathway <- ggarrange(p_path_intensity, p_path_frequency, p_path_if,
                       ncol = 1, labels = c("(a)", "(b)", "(c)"),
                       heights = c(2.1, 1.4, 1.6))
p_pathway

# - Figure S5
ggsave(filename = "./results/figures/picrust2_bacterial_predicted_pathways.pdf",
       plot = p_pathway, width = 12, height = 14, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/picrust2_bacterial_predicted_pathways.png",
       plot = p_pathway, width = 12, height = 14, units = "in", dpi = 300)
