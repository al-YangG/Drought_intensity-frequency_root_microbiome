### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Beta diversity, PCoA, and PERMANOVA ##

## Soil moisture ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(tibble)

## =========================================================
## 1. Helper functions
## =========================================================

prepare_metadata <- function(ps) {
  data.frame(sample_data(ps)) %>%
    rownames_to_column("SampleID") %>%
    mutate(Intensity   = factor(Intensity, levels = c("Dry", "Wet")),
           Frequency   = factor(Frequency, levels = c("Low", "High")),
           Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
           Pot         = factor(Pot))
}

extract_pcoa_df <- function(ps_rel, dataset_name, flip_axis1 = FALSE) {
  dist_bc <- phyloseq::distance(ps_rel, method = "bray")
  ord <- cmdscale(dist_bc, k = 2, eig = TRUE)
  
  eig <- ord$eig
  var_expl <- 100 * eig / sum(eig[eig > 0], na.rm = TRUE)
  
  ord_df <- as.data.frame(ord$points)
  colnames(ord_df) <- c("PCoA1", "PCoA2")
  
  if (flip_axis1) {
    ord_df$PCoA1 <- -ord_df$PCoA1
  }
  
  meta <- prepare_metadata(ps_rel)
  
  ord_df <- ord_df %>%
    rownames_to_column("SampleID") %>%
    left_join(meta, by = "SampleID") %>%
    mutate(Dataset = dataset_name)
  
  list(ord_df = ord_df, dist_bc = dist_bc, ord = ord,
       var_expl = c(PCoA1 = round(var_expl[1], 1), PCoA2 = round(var_expl[2], 1)))
}

run_permanova <- function(dist_obj, meta_df, strata = NULL, permutations = 999) {
  if (is.null(strata)) {
    adonis2(dist_obj ~ Intensity * Frequency * Composition,
            data = meta_df, permutations = permutations, by = "term")
  } else {
    adonis2(dist_obj ~ Intensity * Frequency * Composition,
            data = meta_df, permutations = permutations, by = "term", strata = strata)
  }
}

format_adonis_table <- function(adonis_obj) {
  out <- as.data.frame(adonis_obj)
  out$Term <- rownames(out)
  out[, c("Term", setdiff(names(out), "Term"))]
}

run_betadisper_tests <- function(dist_obj, meta_df, group_var, permutations = 999) {
  group <- meta_df[[group_var]]
  bd <- betadisper(dist_obj, group)
  
  list(model = bd, anova = anova(bd), permutest = permutest(bd, permutations = permutations))
}

format_anova_table <- function(anova_obj) {
  out <- data.frame(anova_obj)
  out$Term <- rownames(out)
  out[, c("Term", setdiff(names(out), "Term"))]
}

format_permutest_table <- function(permutest_obj) {
  out <- as.data.frame(permutest_obj$tab)
  out$Term <- rownames(out)
  out[, c("Term", setdiff(names(out), "Term"))]
}

## =========================================================
## 2. Load phyloseq objects
## =========================================================

ps_Bac <- readRDS("./data/phyloseq/ps_Bac_filtered.rds")
ps_Fun <- readRDS("./data/phyloseq/ps_Fun.rds")

## =========================================================
## 3. Transform to relative abundance
## =========================================================

ps_Bac_rel <- transform_sample_counts(ps_Bac, function(x) x / sum(x))
ps_Fun_rel <- transform_sample_counts(ps_Fun, function(x) x / sum(x))

## =========================================================
## 4. PCoA coordinates
## =========================================================

# Fungal PCoA1 is flipped for visual consistency only
bac_pcoa <- extract_pcoa_df(ps_Bac_rel, dataset_name = "Bacteria", flip_axis1 = FALSE)
fun_pcoa <- extract_pcoa_df(ps_Fun_rel, dataset_name = "Fungi", flip_axis1 = TRUE)

write.csv(bac_pcoa$ord_df, "./results/tables/pcoa_coordinates_bacteria.csv", row.names = FALSE)
write.csv(fun_pcoa$ord_df, "./results/tables/pcoa_coordinates_fungi.csv", row.names = FALSE)

## =========================================================
## 5. Sample-level PERMANOVA
## =========================================================

set.seed(331)

meta_bac <- prepare_metadata(ps_Bac_rel)
meta_fun <- prepare_metadata(ps_Fun_rel)

perm_bac <- run_permanova(bac_pcoa$dist_bc, meta_bac, permutations = 999)
perm_fun <- run_permanova(fun_pcoa$dist_bc, meta_fun, permutations = 999)

perm_bac_strata <- run_permanova(bac_pcoa$dist_bc, meta_bac, strata = meta_bac$Pot, permutations = 999)
perm_fun_strata <- run_permanova(fun_pcoa$dist_bc, meta_fun, strata = meta_fun$Pot,permutations = 999)

# Save - Table S5
write.csv(format_adonis_table(perm_bac), "./results/tables/permanova_bacteria_sample_level.csv", row.names = FALSE)
write.csv(format_adonis_table(perm_fun), "./results/tables/permanova_fungi_sample_level.csv", row.names = FALSE)

write.csv(format_adonis_table(perm_bac_strata), "./results/tables/permanova_bacteria_sample_level_stratified_by_pot.csv", row.names = FALSE)
write.csv(format_adonis_table(perm_fun_strata), "./results/tables/permanova_fungi_sample_level_stratified_by_pot.csv", row.names = FALSE)


## =========================================================
## 5a. Add treatment-level soil moisture to metadata
## =========================================================

soil_moisture_key <- tibble(Treatment = c("Dry-High", "Dry-Low", "Wet-High", "Wet-Low"),
                            SoilMoisture = c(16.35, 18.92, 31.39, 31.36))

add_soil_moisture <- function(meta_df) {
  
  meta_df %>%
    mutate(Treatment = paste(Intensity, Frequency, sep = "-")) %>%
    left_join(soil_moisture_key, by = "Treatment") %>%
    mutate(Treatment = factor(Treatment, levels = c("Dry-High", "Dry-Low", "Wet-High", "Wet-Low")),
           # Center using the values assigned to the samples
           SoilMoisture_c = SoilMoisture - mean(SoilMoisture, na.rm = TRUE))
}

meta_bac_moisture <- add_soil_moisture(meta_bac)
meta_fun_moisture <- add_soil_moisture(meta_fun)

# Check the treatment mapping
meta_bac_moisture %>%
  distinct(Intensity, Frequency, Treatment, SoilMoisture) %>%
  arrange(Treatment)

meta_fun_moisture %>%
  distinct(Intensity, Frequency, Treatment, SoilMoisture) %>%
  arrange(Treatment)

stopifnot(!anyNA(meta_bac_moisture$SoilMoisture))
stopifnot(!anyNA(meta_fun_moisture$SoilMoisture))


stopifnot(identical(rownames(meta_bac_moisture), labels(bac_pcoa$dist_bc)))

stopifnot(identical(rownames(meta_fun_moisture), labels(fun_pcoa$dist_bc)))


align_metadata <- function(meta_df, dist_obj) {
  
  sample_order <- labels(dist_obj)
  
  meta_aligned <- meta_df %>%
    filter(SampleID %in% sample_order) %>%
    arrange(match(SampleID, sample_order))
  
  stopifnot(identical(meta_aligned$SampleID, sample_order))
  
  rownames(meta_aligned) <- meta_aligned$SampleID
  
  meta_aligned
}

meta_bac_moisture <- align_metadata(meta_bac_moisture, bac_pcoa$dist_bc)

meta_fun_moisture <- align_metadata(meta_fun_moisture, fun_pcoa$dist_bc)



## =========================================================
## 5b. Soil-moisture PERMANOVA
## =========================================================

run_moisture_permanova <- function(dist_obj, meta_df, permutations = 999) {
  
  vegan::adonis2(dist_obj ~ SoilMoisture_c * Composition, data = meta_df,
                 permutations = permutations, by = "terms")
}

set.seed(331)

perm_bac_moisture <- run_moisture_permanova(bac_pcoa$dist_bc, meta_bac_moisture)

set.seed(331)

perm_fun_moisture <- run_moisture_permanova(fun_pcoa$dist_bc, meta_fun_moisture)

perm_bac_moisture
perm_fun_moisture


write.csv(format_adonis_table(perm_bac_moisture),
          "./results/tables/permanova_bacteria_soil_moisture.csv", row.names = FALSE)

write.csv(format_adonis_table(perm_fun_moisture),
          "./results/tables/permanova_fungi_soil_moisture.csv", row.names = FALSE)

extract_permanova_effects <- function(adonis_object, microbe, model_name) {
  
  as.data.frame(adonis_object) %>%
    rownames_to_column("Term") %>%
    filter(!Term %in% c("Residual", "Total")) %>%
    transmute(Microbe = microbe, Model = model_name,
              Term, Df, SumOfSqs, R2, F = F, p_value = `Pr(>F)`)
}

permanova_comparison <- bind_rows(
  extract_permanova_effects(perm_bac,"Bacteria", "Factorial watering"),
  extract_permanova_effects(perm_bac_moisture, "Bacteria", "Soil moisture"),
  extract_permanova_effects(perm_fun, "Fungi", "Factorial watering"),
  extract_permanova_effects(perm_fun_moisture, "Fungi", "Soil moisture"))

permanova_comparison

write.csv(permanova_comparison, row.names = FALSE,
          "./results/tables/permanova_factorial_vs_soil_moisture.csv")



## Adjusted R2 comparison using dbRDA

extract_adjusted_r2 <- function(dist_obj, meta_df, microbe) {
  
  factorial_model <- vegan::dbrda(dist_obj ~ Intensity * Frequency * Composition, data = meta_df)
  
  moisture_model <- vegan::dbrda(dist_obj ~ SoilMoisture_c * Composition, data = meta_df)
  
  factorial_r2 <- vegan::RsquareAdj(factorial_model)
  moisture_r2  <- vegan::RsquareAdj(moisture_model)
  
  tibble(Microbe = microbe, R2_Factorial = unname(factorial_r2$r.squared),
         Adjusted_R2_Factorial = unname(factorial_r2$adj.r.squared),
         R2_SoilMoisture = unname(moisture_r2$r.squared),
         Adjusted_R2_SoilMoisture = unname(moisture_r2$adj.r.squared),
         Delta_Adjusted_R2_Moisture_vs_Factorial = unname(moisture_r2$adj.r.squared - factorial_r2$adj.r.squared))
}

multivariate_r2_comparison <- bind_rows(
  extract_adjusted_r2(bac_pcoa$dist_bc, meta_bac_moisture, "Bacteria"),
  extract_adjusted_r2(fun_pcoa$dist_bc, meta_fun_moisture, "Fungi"))

multivariate_r2_comparison

write.csv(multivariate_r2_comparison, row.names = FALSE,
          "./results/tables/multivariate_adjusted_r2_comparison.csv")
