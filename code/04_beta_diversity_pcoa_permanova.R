### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Beta diversity, PCoA, and PERMANOVA ##
## Date: April 2026 | Author: Yang ##

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

perm_bac_strata <- run_permanova(bac_pcoa$dist_bc, meta_bac, strata = meta_bac$Pot,
                                 permutations = 999)

perm_fun_strata <- run_permanova(fun_pcoa$dist_bc, meta_fun, strata = meta_fun$Pot,
                                 permutations = 999)

# Save - Table S5
write.csv(format_adonis_table(perm_bac), "./results/tables/permanova_bacteria_sample_level.csv", row.names = FALSE)
write.csv(format_adonis_table(perm_fun), "./results/tables/permanova_fungi_sample_level.csv", row.names = FALSE)

write.csv(format_adonis_table(perm_bac_strata), "./results/tables/permanova_bacteria_sample_level_stratified_by_pot.csv", row.names = FALSE)
write.csv(format_adonis_table(perm_fun_strata), "./results/tables/permanova_fungi_sample_level_stratified_by_pot.csv", row.names = FALSE)

## =========================================================
## 6. Multivariate dispersion tests
## =========================================================

dispersion_vars <- c("Intensity", "Frequency", "Composition")

# Bacteria
for (v in dispersion_vars) {
  res <- run_betadisper_tests(bac_pcoa$dist_bc, meta_bac, v, permutations = 999)
  
  write.csv(format_anova_table(res$anova),
            file = paste0("./results/tables/beta_dispersion_bacteria_", tolower(v), "_anova.csv"),
            row.names = FALSE)
  
  write.csv(format_permutest_table(res$permutest),
            file = paste0("./results/tables/beta_dispersion_bacteria_", tolower(v), "_permutest.csv"),
            row.names = FALSE)
}

# Fungi
for (v in dispersion_vars) {
  res <- run_betadisper_tests(fun_pcoa$dist_bc, meta_fun, v, permutations = 999)
  
  write.csv(format_anova_table(res$anova),
            file = paste0("./results/tables/beta_dispersion_fungi_", tolower(v), "_anova.csv"),
            row.names = FALSE)
  
  write.csv(format_permutest_table(res$permutest),
            file = paste0("./results/tables/beta_dispersion_fungi_", tolower(v), "_permutest.csv"),
            row.names = FALSE)
}

## =========================================================
## 7. Save variance explained by ordination axes
## =========================================================

pcoa_variance <- bind_rows(data.frame(Dataset = "Bacteria", Axis = c("PCoA1", "PCoA2"),
                                      Variance_explained_percent = unname(bac_pcoa$var_expl)),
                           data.frame(Dataset = "Fungi", Axis = c("PCoA1", "PCoA2"),
                                      Variance_explained_percent = unname(fun_pcoa$var_expl)))

write.csv(pcoa_variance, "./results/tables/pcoa_variance_explained.csv", row.names = FALSE)
