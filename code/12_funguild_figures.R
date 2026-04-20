### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: FUNGuild trophic mode figures ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

## =========================================================
## 1. Load saved tables
## =========================================================

trophic_rel <- read.csv("./results/tables/funguild_trophic_relative_abundance_long.csv",
                        header = TRUE, check.names = FALSE)

guild_stats <- read.csv("./results/tables/funguild_trophic_mode_beta_models.csv",
                        header = TRUE, check.names = FALSE)

## =========================================================
## 2. Restore factor levels
## =========================================================

trophic_mode_levels <- c("Saprotroph", "Symbiotroph", "Pathotroph",
                         "Multi-trophic", "Unassigned")

ecological_mode_levels <- c("Saprotroph", "Symbiotroph", "Pathotroph", "Multi-trophic")

trophic_rel <- trophic_rel %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         trophicMode_clean = factor(trophicMode_clean, levels = trophic_mode_levels))

## =========================================================
## 3. Main stacked bar plot across treatments
## =========================================================

plot_I <- trophic_rel %>%
  group_by(Level = Intensity, trophicMode_clean) %>%
  summarise(RelAbundance = mean(RelAbundance), .groups = "drop") %>%
  mutate(Panel = "Intensity")

plot_F <- trophic_rel %>%
  group_by(Level = Frequency, trophicMode_clean) %>%
  summarise(RelAbundance = mean(RelAbundance), .groups = "drop") %>%
  mutate(Panel = "Frequency")

plot_C <- trophic_rel %>%
  group_by(Level = Composition, trophicMode_clean) %>%
  summarise(RelAbundance = mean(RelAbundance), .groups = "drop") %>%
  mutate(Panel = "Composition")

plot_main <- bind_rows(plot_I, plot_F, plot_C) %>%
  mutate(Panel = factor(Panel, levels = c("Intensity", "Frequency", "Composition")),
         Level = factor(Level, levels = c("Dry", "Wet", "Low", "High", "Conspecific", "Heterospecific")),
         trophicMode_clean = factor(trophicMode_clean, levels = trophic_mode_levels))

guild_cols <- c("Saprotroph"    = "#F28E7F",
                "Symbiotroph"   = "#C2A83E",
                "Pathotroph"    = "#67C5B8",
                "Multi-trophic" = "#F4EFE6",
                "Unassigned"    = "#D9D9D9")

p_guild_main <- ggplot(plot_main, aes(x = Level, y = RelAbundance, fill = trophicMode_clean)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = guild_cols) +
  # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance", fill = "Fungal trophic mode") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(face = "bold", size = 9))
p_guild_main

## =========================================================
## 4. Intensity response plot with significance labels
## =========================================================

plot_intensity <- trophic_rel %>%
  filter(!is.na(trophicMode_clean), trophicMode_clean %in% ecological_mode_levels) %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         trophicMode_clean = factor(trophicMode_clean, levels = ecological_mode_levels))

p_lab <- guild_stats %>%
  filter(Term == "Intensity") %>%
  mutate(label = case_when(p < 0.001 ~ "***", p < 0.01  ~ "**", p < 0.05  ~ "*",
                           TRUE      ~ "")) %>%
  dplyr::select(Guild, label) %>%
  rename(trophicMode_clean = Guild)

ann_y <- plot_intensity %>%
  group_by(trophicMode_clean) %>%
  summarise(y = max(RelAbundance, na.rm = TRUE) * 0.7, .groups = "drop")

ann_df <- left_join(p_lab, ann_y, by = "trophicMode_clean") %>%
  mutate(x = 1.5, trophicMode_clean = factor(trophicMode_clean, levels = ecological_mode_levels))

intensity_cols <- c("Dry" = "#E7782F", "Wet" = "#1994E5")

p_intensity <- ggplot(plot_intensity, aes(x = Intensity, y = RelAbundance, color = Intensity)) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.15) +
  stat_summary(fun = mean, geom = "point", size = 2.8) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               width = 0.10, linewidth = 0.8) +
  geom_text(data = ann_df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  facet_wrap(~ trophicMode_clean, scales = "free_y", nrow = 2) +
  scale_color_manual(values = intensity_cols) +
  scale_y_continuous(breaks = pretty_breaks(4), expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Intensity", y = "Relative abundance") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 11),
        axis.title.x = element_text(face = "bold", color = "black", size = 11),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", linewidth = 0.2),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing.x = unit(0.9, "lines"))
p_intensity

## =========================================================
## 5. Combine and save figure
## =========================================================

p_func <- ggarrange(p_guild_main, p_intensity, nrow = 1, labels = c("(a)", "(b)"))
p_func

# Save - Figure S6
ggsave(filename = "./results/figures/fungal_trophic_modes.pdf", plot = p_func,
       width = 10, height = 6, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/fungal_trophic_modes.png", plot = p_func,
       width = 10, height = 6, units = "in", dpi = 300)
