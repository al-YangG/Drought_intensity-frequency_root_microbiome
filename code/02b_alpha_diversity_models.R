### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Alpha diversity mixed-effects models ##

## Soil Moisture ##

# Change to working directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(r2glmm)
library(performance)

## =========================================================
## 1. Load alpha-diversity data
## =========================================================

bac <- read.csv("./results/tables/alpha_diversity_bacteria.csv", header = TRUE, check.names = FALSE)

fun <- read.csv("./results/tables/alpha_diversity_fungi.csv", header = TRUE, check.names = FALSE)

## =========================================================
## 2. Create treatment–soil-moisture lookup table
## =========================================================

# Mean VWC during treatment periods, calculated from the logger data
soil_moisture_key <- tibble(Treatment = c("Dry-High", "Dry-Low", "Wet-High", "Wet-Low"),
                            SoilMoisture = c(16.35, 18.92, 31.39, 31.36))

soil_moisture_key


## =========================================================
## 3. Prepare data and add soil moisture
## =========================================================

prepare_alpha <- function(df) {
  df %>%
    mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
           Frequency = factor(Frequency, levels = c("Low", "High")),
           Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
      Pot = factor(Pot),
      
      # Combine water amount and frequency
      Treatment = paste(Intensity, Frequency, sep = "-"),
      
      # Richness transformation
      Observed_log = log1p(Observed)) %>%
    left_join(soil_moisture_key, by = "Treatment") %>%
    mutate(Treatment = factor(Treatment, levels = c("Dry-High", "Dry-Low", "Wet-High", "Wet-Low")),
      
      # Optional centering improves interpretation
      SoilMoisture_c = SoilMoisture - mean(SoilMoisture, na.rm = TRUE))
}

bac <- prepare_alpha(bac)
fun <- prepare_alpha(fun)

# Confirm that the join worked
bac %>%
  distinct(Intensity, Frequency, Treatment, SoilMoisture) %>%
  arrange(Treatment)

fun %>%
  distinct(Intensity, Frequency, Treatment, SoilMoisture) %>%
  arrange(Treatment)



options(contrasts = c("contr.sum", "contr.poly"))

## =========================================================
## 4. Fit treatment and soil-moisture models
## =========================================================

fit_comparison_models <- function(df) {
  
  treatment_models <- list(
    Richness = lmer(Observed_log ~ Intensity * Frequency * Composition + (1 | Pot), data = df, REML = FALSE),
    Shannon = lmer(Shannon ~ Intensity * Frequency * Composition + (1 | Pot), data = df, REML = FALSE),
    Simpson = lmer(Simpson ~ Intensity * Frequency * Composition + (1 | Pot), data = df, REML = FALSE))
  
  moisture_models <- list(
    Richness = lmer(Observed_log ~ SoilMoisture_c * Composition + (1 | Pot), data = df, REML = FALSE),
    Shannon = lmer(Shannon ~ SoilMoisture_c * Composition + (1 | Pot), data = df, REML = FALSE),
    Simpson = lmer(Simpson ~ SoilMoisture_c * Composition + (1 | Pot), data = df, REML = FALSE))
  
  list(Treatment = treatment_models, SoilMoisture = moisture_models)
}

models_bac_compare <- fit_comparison_models(bac)
models_fun_compare <- fit_comparison_models(fun)



## =========================================================
## 5. Extract soil-moisture ANOVA results
## =========================================================

extract_moisture_anova <- function(model, microbe, metric) {
  as.data.frame(car::Anova(model, type = 3)) %>%
    tibble::rownames_to_column("Term") %>%
    mutate(Microbe = microbe, Metric = metric) %>%
    select(Microbe, Metric, Term, Chisq, Df, `Pr(>Chisq)`)
}

moisture_anova <- bind_rows(
  extract_moisture_anova(models_bac_compare$SoilMoisture$Richness, "Bacteria", "Richness"),
  extract_moisture_anova(models_bac_compare$SoilMoisture$Shannon, "Bacteria", "Shannon"),
  extract_moisture_anova(models_bac_compare$SoilMoisture$Simpson, "Bacteria", "Simpson"),
  extract_moisture_anova(models_fun_compare$SoilMoisture$Richness, "Fungi", "Richness"),
  extract_moisture_anova(models_fun_compare$SoilMoisture$Shannon, "Fungi", "Shannon"),
  extract_moisture_anova(models_fun_compare$SoilMoisture$Simpson, "Fungi", "Simpson"))

moisture_anova

write.csv(moisture_anova, "./results/tables/alpha_diversity_soil_moisture_anova.csv", row.names = FALSE)


## =========================================================
## 6. Compare model performance
## =========================================================

compare_one_metric <- function(treatment_model,
                               moisture_model,
                               microbe,
                               metric) {
  
  # Ensure that paired models use the same observations
  if (nobs(treatment_model) != nobs(moisture_model)) {
    stop(paste(microbe, metric, "models were fitted with different numbers of observations."))
  }
  
  # Extract Nakagawa R2
  r2_treatment <- performance::r2_nakagawa(treatment_model)
  r2_moisture  <- performance::r2_nakagawa(moisture_model)
  
  marginal_r2_treatment <- unname(r2_treatment$R2_marginal)
  
  marginal_r2_moisture <- unname(r2_moisture$R2_marginal)
  
  conditional_r2_treatment <- unname(r2_treatment$R2_conditional)
  
  conditional_r2_moisture <- unname(r2_moisture$R2_conditional)
  
  # Number of observations
  n_treatment <- nobs(treatment_model)
  n_moisture  <- nobs(moisture_model)
  
  # Number of fixed-effect predictors, excluding the intercept
  p_treatment <- length(lme4::fixef(treatment_model)) - 1
  
  p_moisture <- length(lme4::fixef(moisture_model)) - 1
  
  # Conventional adjustment applied to marginal Nakagawa R2
  adjusted_marginal_r2_treatment <-
    1 - (1 - marginal_r2_treatment) *
    ((n_treatment - 1) /
       (n_treatment - p_treatment - 1))
  
  adjusted_marginal_r2_moisture <-
    1 - (1 - marginal_r2_moisture) *
    ((n_moisture - 1) /
       (n_moisture - p_moisture - 1))
  
  # Extract AIC
  aic_treatment <- AIC(treatment_model)
  aic_moisture  <- AIC(moisture_model)
  
  tibble(Microbe = microbe, Metric = metric, N = n_treatment,
         Fixed_Effect_Parameters_Treatment = p_treatment,
         Fixed_Effect_Parameters_SoilMoisture = p_moisture,
         AIC_Treatment = aic_treatment, AIC_SoilMoisture = aic_moisture,
    
    # Negative values favor the soil-moisture model
    Delta_AIC_Moisture_vs_Treatment = aic_moisture - aic_treatment,
    
    Better_Supported_Model_AIC = case_when(aic_moisture + 2 < aic_treatment ~ "Soil moisture",
                                           aic_treatment + 2 < aic_moisture ~ "Treatment",
                                           TRUE ~ "Similar support"),
    
    Marginal_R2_Treatment = marginal_r2_treatment,
    
    Marginal_R2_SoilMoisture = marginal_r2_moisture,
    
    Adjusted_Marginal_R2_Treatment = adjusted_marginal_r2_treatment,
    
    Adjusted_Marginal_R2_SoilMoisture = adjusted_marginal_r2_moisture,
    
    # Positive values favor the soil-moisture model
    Delta_Adjusted_Marginal_R2_Moisture_vs_Treatment =
      adjusted_marginal_r2_moisture - adjusted_marginal_r2_treatment,
    
    Better_Model_Adjusted_R2 = case_when(adjusted_marginal_r2_moisture >
                                           adjusted_marginal_r2_treatment ~ "Soil moisture",
                                         adjusted_marginal_r2_treatment >
                                           adjusted_marginal_r2_moisture ~ "Treatment",
                                         TRUE ~ "Equal"),
    
    Conditional_R2_Treatment = conditional_r2_treatment,
    Conditional_R2_SoilMoisture = conditional_r2_moisture)
}


model_comparison <- bind_rows(compare_one_metric(models_bac_compare$Treatment$Richness,
                                                 models_bac_compare$SoilMoisture$Richness,
                                                 "Bacteria", "Richness"),
                              compare_one_metric(models_bac_compare$Treatment$Shannon,
                                                 models_bac_compare$SoilMoisture$Shannon,
                                                 "Bacteria", "Shannon"),
                              compare_one_metric(models_bac_compare$Treatment$Simpson,
                                                 models_bac_compare$SoilMoisture$Simpson,
                                                 "Bacteria", "Simpson"),
                              compare_one_metric(models_fun_compare$Treatment$Richness,
                                                 models_fun_compare$SoilMoisture$Richness,
                                                 "Fungi", "Richness"),
                              compare_one_metric(models_fun_compare$Treatment$Shannon,
                                                 models_fun_compare$SoilMoisture$Shannon,
                                                 "Fungi", "Shannon"),
                              compare_one_metric(models_fun_compare$Treatment$Simpson,
                                                 models_fun_compare$SoilMoisture$Simpson,
                                                 "Fungi", "Simpson"))

model_comparison

# Save the comparison table
write.csv(model_comparison, "./results/tables/alpha_diversity_soil_moisture_model_comparison.csv", row.names = FALSE)

