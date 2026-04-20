### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Alpha diversity mixed-effects models ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(r2glmm)

## =========================================================
## 1. Load data
## =========================================================

bac <- read.csv("./results/tables/alpha_diversity_bacteria.csv", header = TRUE, check.names = FALSE)
fun <- read.csv("./results/tables/alpha_diversity_fungi.csv", header = TRUE, check.names = FALSE)

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

bac <- prepare_alpha(bac)
fun <- prepare_alpha(fun)

options(contrasts = c("contr.sum", "contr.poly"))

## =========================================================
## 3. Fit final models
## Richness is modeled as log-transformed Observed richness
## =========================================================

fit_alpha_models <- function(df) {
  list(Richness = lmer(Observed_log ~ Intensity * Frequency * Composition + (1 | Pot), data = df),
       Shannon  = lmer(Shannon      ~ Intensity * Frequency * Composition + (1 | Pot), data = df),
       Simpson  = lmer(Simpson      ~ Intensity * Frequency * Composition + (1 | Pot), data = df))
}

models_bac <- fit_alpha_models(bac)
models_fun <- fit_alpha_models(fun)

## =========================================================
## 4. Extract ANOVA tables
## =========================================================

extract_anova <- function(model, microbe, metric) {
  as.data.frame(car::Anova(model, type = 3)) %>%
    tibble::rownames_to_column("Term") %>%
    mutate(Microbe = microbe, Metric  = metric) %>%
    select(Microbe, Metric, Term, Chisq, Df, `Pr(>Chisq)`)
}

anova_all <- bind_rows(extract_anova(models_bac$Richness, "Bacteria", "Richness"),
                       extract_anova(models_bac$Shannon,  "Bacteria", "Shannon"),
                       extract_anova(models_bac$Simpson,  "Bacteria", "Simpson"),
                       extract_anova(models_fun$Richness, "Fungi",    "Richness"),
                       extract_anova(models_fun$Shannon,  "Fungi",    "Shannon"),
                       extract_anova(models_fun$Simpson,  "Fungi",    "Simpson")) %>%
  arrange(Microbe, Metric)

# Save - Table S4
write.csv(anova_all, "./results/tables/alpha_diversity_anova_results.csv", row.names = FALSE)

## =========================================================
## 5. Extract variance explained
## =========================================================

extract_r2 <- function(model, microbe, metric) {
  as.data.frame(r2beta(model, method = "nsj", partial = TRUE)) %>%
    mutate(Microbe = microbe, Metric  = metric) %>%
    filter(Effect != "Model") %>%
    select(Microbe, Metric, Effect, Rsq, lower.CL, upper.CL)
}

variance_all <- bind_rows(extract_r2(models_bac$Richness, "Bacteria", "Richness"),
                          extract_r2(models_bac$Shannon,  "Bacteria", "Shannon"),
                          extract_r2(models_bac$Simpson,  "Bacteria", "Simpson"),
                          extract_r2(models_fun$Richness, "Fungi",    "Richness"),
                          extract_r2(models_fun$Shannon,  "Fungi",    "Shannon"),
                          extract_r2(models_fun$Simpson,  "Fungi",    "Simpson")) %>%
  arrange(Microbe, Metric)

write.csv(variance_all, "./results/tables/alpha_diversity_variance_explained.csv", row.names = FALSE)
