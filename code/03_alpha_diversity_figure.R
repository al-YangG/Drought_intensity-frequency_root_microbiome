### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Alpha diversity visualization ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(ggplot2)
library(ggpubr)
library(tibble)
library(grid)

## =========================================================
## 1. Load data
## =========================================================

bac <- read.csv("./results/tables/alpha_diversity_bacteria.csv", header = TRUE, check.names = FALSE)
fun <- read.csv("./results/tables/alpha_diversity_fungi.csv", header = TRUE, check.names = FALSE)
anova_all <- read.csv("./results/tables/alpha_diversity_anova_results.csv", header = TRUE, check.names = FALSE)

## =========================================================
## 2. Prepare data
## =========================================================

prepare_alpha <- function(df) {
  df %>%
    mutate(Intensity   = factor(Intensity, levels = c("Dry", "Wet")),
           Frequency   = factor(Frequency, levels = c("Low", "High")),
           Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
           Pot         = factor(Pot),
           Observed_log = log1p(Observed))
}

bac <- prepare_alpha(bac) %>% mutate(Microbes = "Bacteria")
fun <- prepare_alpha(fun) %>% mutate(Microbes = "Fungi")

bac_long <- bac %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "Metric", values_to = "Value")

fun_long <- fun %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "Metric", values_to = "Value")

alpha_all <- bind_rows(bac_long, fun_long) %>%
  filter(Metric %in% c("Observed", "Shannon")) %>%
  mutate(Metric = dplyr::recode(Metric, "Observed" = "Richness"),
         Microbes = factor(Microbes, levels = c("Bacteria", "Fungi")),
         Metric = factor(Metric, levels = c("Richness", "Shannon")))

options(contrasts = c("contr.sum", "contr.poly"))

## =========================================================
## 3. Refit only models needed for letters in the figure
## =========================================================

models_bac <- list(Richness = lmer(Observed_log ~ Intensity * Frequency * Composition + (1 | Pot), data = bac),
                   Shannon  = lmer(Shannon      ~ Intensity * Frequency * Composition + (1 | Pot), data = bac))

models_fun <- list(Richness = lmer(Observed_log ~ Intensity * Frequency * Composition + (1 | Pot), data = fun),
                   Shannon  = lmer(Shannon      ~ Intensity * Frequency * Composition + (1 | Pot), data = fun))

## =========================================================
## 4. Significance label helpers
## =========================================================

p_to_stars <- function(p) {
  case_when(is.na(p)      ~ "",
            p < 0.001     ~ "***",
            p < 0.01      ~ "**",
            p < 0.05      ~ "*",
            TRUE          ~ "")
}

factor_order <- c("Intensity", "Frequency", "Composition")

sig_labels <- anova_all %>%
  filter(Metric %in% c("Richness", "Shannon"), Term %in% factor_order) %>%
  transmute(Microbes = factor(Microbe, levels = c("Bacteria", "Fungi")),
            Metric = factor(Metric, levels = c("Richness", "Shannon")),
            Factor = factor(Term, levels = factor_order),
            label = p_to_stars(`Pr(>Chisq)`))

sig_if <- anova_all %>%
  filter(Metric %in% c("Richness", "Shannon"), Term == "Intensity:Frequency") %>%
  transmute(Microbes = factor(Microbe, levels = c("Bacteria", "Fungi")),
            Metric = factor(Metric, levels = c("Richness", "Shannon")),
            label = p_to_stars(`Pr(>Chisq)`))

## =========================================================
## 5. Main-effect plot data
## =========================================================

alpha_main <- alpha_all %>%
  mutate(Intensity   = factor(Intensity,   levels = c("Dry", "Wet")),
         Frequency   = factor(Frequency,   levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific"))) %>%
  pivot_longer(cols = c(Intensity, Frequency, Composition),
               names_to = "Factor", values_to = "Treatment") %>%
  mutate(Factor = factor(Factor, levels = factor_order))

sig_labels <- sig_labels %>%
  left_join(alpha_main %>%
              group_by(Microbes, Metric, Factor) %>%
              summarise(y = max(Value, na.rm = TRUE) * 0.95, .groups = "drop"),
            by = c("Microbes", "Metric", "Factor"))

## =========================================================
## 6. Letters for main effects
## =========================================================

get_letters <- function(model, factor_name, microbe, metric) {
  emm <- emmeans(model, as.formula(paste("~", factor_name)))
  cld_df <- multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame()
  
  cld_df %>%
    dplyr::select(Treatment = all_of(factor_name), .group) %>%
    mutate(Letters  = gsub(" ", "", .group),
           Factor   = factor_name,
           Microbes = microbe,
           Metric   = metric) %>%
    dplyr::select(Microbes, Metric, Factor, Treatment, Letters)
}

letters_df <- bind_rows(
  get_letters(models_bac$Richness, "Intensity",   "Bacteria", "Richness"),
  get_letters(models_bac$Richness, "Frequency",   "Bacteria", "Richness"),
  get_letters(models_bac$Richness, "Composition", "Bacteria", "Richness"),
  get_letters(models_bac$Shannon,  "Intensity",   "Bacteria", "Shannon"),
  get_letters(models_bac$Shannon,  "Frequency",   "Bacteria", "Shannon"),
  get_letters(models_bac$Shannon,  "Composition", "Bacteria", "Shannon"),
  get_letters(models_fun$Richness, "Intensity",   "Fungi",    "Richness"),
  get_letters(models_fun$Richness, "Frequency",   "Fungi",    "Richness"),
  get_letters(models_fun$Richness, "Composition", "Fungi",    "Richness"),
  get_letters(models_fun$Shannon,  "Intensity",   "Fungi",    "Shannon"),
  get_letters(models_fun$Shannon,  "Frequency",   "Fungi",    "Shannon"),
  get_letters(models_fun$Shannon,  "Composition", "Fungi",    "Shannon")) %>%
  mutate(Microbes  = factor(Microbes, levels = c("Bacteria", "Fungi")),
         Metric    = factor(Metric, levels = c("Richness", "Shannon")),
         Factor    = factor(Factor, levels = factor_order),
         Treatment = as.character(Treatment))

letter_pos <- alpha_main %>%
  group_by(Microbes, Metric, Factor, Treatment) %>%
  summarise(y = max(Value, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  mutate(Treatment = as.character(Treatment))

letters_df <- letters_df %>%
  left_join(letter_pos, by = c("Microbes", "Metric", "Factor", "Treatment"))

## =========================================================
## 7. Main-effect panel
## =========================================================

c1 <- c("Dry" = "#E7782F", "Wet" = "#1994E5")
c2 <- c("Low" = "#E5D019", "High" = "#B38AF1")
c3 <- c("Conspecific" = "#D95F8D", "Heterospecific" = "#33A27F")

p1 <- ggplot(alpha_main, aes(x = Treatment, y = Value, color = Treatment)) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.15) +
  stat_summary(fun = mean, geom = "point", size = 2.8) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.10, linewidth = 0.8) +
  geom_text(data = sig_labels, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_grid(Microbes + Metric ~ Factor, scales = "free", switch = "y") +
  scale_color_manual(values = c(c1, c2, c3)) +
  labs(x = "Treatment level", y = "Alpha diversity") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 11),
        axis.title.x = element_text(face = "bold", color = "black", size = 11),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", linewidth = 0.2),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing.x = unit(0.9, "lines"))
p1

## =========================================================
## 8. Interaction panel: Intensity × Frequency
## =========================================================

alpha_if <- alpha_all %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Interaction = factor("Intensity × Frequency", levels = "Intensity × Frequency"))

sig_if <- sig_if %>%
  left_join(alpha_if %>%
              group_by(Microbes, Metric) %>%
              summarise(y = max(Value, na.rm = TRUE) * 0.95, .groups = "drop"),
            by = c("Microbes", "Metric"))

get_if_letters <- function(model, microbe, metric) {
  emm <- emmeans(model, ~ Intensity | Frequency)
  cld_df <- multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame()
  
  cld_df %>%
    dplyr::select(Intensity, Frequency, .group) %>%
    mutate(Letters = gsub(" ", "", .group),
           Microbes = microbe, Metric = metric,
           Interaction = "Intensity × Frequency") %>%
    dplyr::select(Microbes, Metric, Interaction, Frequency, Intensity, Letters)
}

letters_if <- bind_rows(get_if_letters(models_bac$Richness, "Bacteria", "Richness"),
                        get_if_letters(models_bac$Shannon,  "Bacteria", "Shannon"),
                        get_if_letters(models_fun$Richness, "Fungi",    "Richness"),
                        get_if_letters(models_fun$Shannon,  "Fungi",    "Shannon")) %>%
  mutate(Microbes    = factor(Microbes, levels = c("Bacteria", "Fungi")),
         Metric      = factor(Metric, levels = c("Richness", "Shannon")),
         Frequency   = factor(Frequency, levels = c("Low", "High")),
         Intensity   = factor(Intensity, levels = c("Dry", "Wet")),
         Interaction = factor(Interaction, levels = "Intensity × Frequency"))

letter_if_pos <- alpha_if %>%
  group_by(Microbes, Metric, Frequency, Intensity) %>%
  summarise(
    y = max(Value, na.rm = TRUE) + 0.06 * diff(range(Value, na.rm = TRUE)),
    .groups = "drop"
  )

letters_if <- letters_if %>%
  left_join(letter_if_pos, by = c("Microbes", "Metric", "Frequency", "Intensity"))

pd <- position_dodge(width = 0.35)

p2 <- ggplot(alpha_if, aes(x = Frequency, y = Value, color = Intensity, group = Intensity)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.22, dodge.width = 0.35),
              size = 0.8, alpha = 0.15) +
  stat_summary(fun = mean, geom = "point", position = pd, size = 2.6) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = pd, width = 0.10, linewidth = 0.7) +
  geom_text(data = sig_if, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE,
            fontface = "bold", size = 4) +
  geom_text(data = letters_if, aes(x = Frequency, y = y, label = Letters, group = Intensity),
            position = pd, inherit.aes = FALSE, fontface = "bold", size = 3.5,
            color = "black", show.legend = FALSE) +
  facet_grid(Microbes + Metric ~ Interaction, scales = "free_y", switch = "y") +
  scale_color_manual(values = c("Wet" = "#1994E5", "Dry" = "#E7782F")) +
  labs(x = "Frequency", y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 11),
        axis.title.x = element_text(face = "bold", color = "black", size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.61, 0.15),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "gray95"),
        legend.box.background = element_blank(),
        legend.text = element_text(face = "bold", color = "black", size = 9),
        legend.title = element_text(face = "bold", color = "black", size = 9),
        panel.grid = element_blank(),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        strip.background.x = element_rect(fill = "grey95", linewidth = 0.2),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing.x = unit(0.9, "lines"))
p2

## =========================================================
## 9. Combine and save figure - Figure 1
## =========================================================

p <- ggarrange(p1, p2, nrow = 1, widths = c(2.6, 1))

ggsave(filename = "./results/figures/alpha_diversity_figure.pdf", plot = p,
       width = 12, height = 8, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/alpha_diversity_figure.png", plot = p,
       width = 12, height = 8, dpi = 300)

