### Project: Drought intensity and frequency (DIF) - Microbes ###
## Indicator species analysis (ASVs) ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

library(phyloseq)
library(indicspecies)
library(dplyr)
library(stringr)
library(permute)

# =========================================================
# 1. Load data
# =========================================================

ps_Bac <- readRDS("./data/phyloseq/ps_Bac_filtered.rds")
ps_Fun <- readRDS("./data/phyloseq/ps_Fun.rds")

# =========================================================
# 2. Indicator species analysis function
# =========================================================

run_indval <- function(ps, grouping_var, min_occurrence = 3,
                       seed = 39, nperm = 999) {
  
  asv_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) asv_mat <- t(asv_mat)
  
  meta <- data.frame(sample_data(ps))
  meta <- meta[rownames(asv_mat), , drop = FALSE]
  
  group <- factor(meta[[grouping_var]])
  
  # filter rare taxa
  keep_taxa <- colSums(asv_mat > 0) >= min_occurrence
  asv_filt <- asv_mat[, keep_taxa, drop = FALSE]
  
  set.seed(seed)
  res <- multipatt(x = asv_filt, cluster = group, func = "IndVal.g",
                   control = how(nperm = nperm))
  
  tab <- data.frame(ASV = rownames(res$sign), res$sign, check.names = FALSE) %>%
    mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
    arrange(p.value)
  
  tax_df <- as.data.frame(tax_table(ps))
  tax_df$ASV <- rownames(tax_df)
  
  left_join(tab, tax_df, by = "ASV")
}

# =========================================================
# 3. Helper: assign indicator group labels
# =========================================================

prep_indval_sig <- function(indval_tab, factor_name) {
  
  out <- indval_tab %>%
    filter(p_adj < 0.05)
  
  if (factor_name == "Intensity") {
    out <- out %>%
      mutate(Indicator = case_when(s.Wet == 1 ~ "Wet", s.Dry == 1 ~ "Dry",
                                   TRUE ~ NA_character_))
  }
  
  if (factor_name == "Frequency") {
    out <- out %>%
      mutate(Indicator = case_when(s.Low == 1 ~ "Low", s.High == 1 ~ "High",
                                   TRUE ~ NA_character_))
  }
  
  if (factor_name == "Composition") {
    out <- out %>%
      mutate(Indicator = case_when(s.Conspecific == 1 ~ "Conspecific",
                                   s.Heterospecific == 1 ~ "Heterospecific",
                                   TRUE ~ NA_character_))
  }
  
  out %>%
    dplyr::select(ASV, Indicator, stat, p.value, p_adj,
                  Phylum, Class, Order, Family, Genus) %>%
    arrange(desc(stat))
}

# =========================================================
# 4. Run indicator analysis
# =========================================================

# --- Bacteria ---
bac_indval_intensity   <- run_indval(ps_Bac, "Intensity")
bac_indval_frequency   <- run_indval(ps_Bac, "Frequency")
bac_indval_composition <- run_indval(ps_Bac, "Composition")

# --- Fungi ---
fun_indval_intensity   <- run_indval(ps_Fun, "Intensity")
fun_indval_frequency   <- run_indval(ps_Fun, "Frequency")
fun_indval_composition <- run_indval(ps_Fun, "Composition")

# =========================================================
# 5. Extract significant indicators
# =========================================================

# Bacteria
bac_sig_intensity   <- prep_indval_sig(bac_indval_intensity, "Intensity")
bac_sig_frequency   <- prep_indval_sig(bac_indval_frequency, "Frequency")
bac_sig_composition <- prep_indval_sig(bac_indval_composition, "Composition")

# Fungi
fun_sig_intensity   <- prep_indval_sig(fun_indval_intensity, "Intensity")
fun_sig_frequency   <- prep_indval_sig(fun_indval_frequency, "Frequency")
fun_sig_composition <- prep_indval_sig(fun_indval_composition, "Composition")

# =========================================================
# 6. Summary: number of significant indicators
# =========================================================

cat("Bacteria indicators:\n")
cat("  Intensity:", nrow(bac_sig_intensity), "\n")
cat("  Frequency:", nrow(bac_sig_frequency), "\n")
cat("  Composition:", nrow(bac_sig_composition), "\n\n")

cat("Fungi indicators:\n")
cat("  Intensity:", nrow(fun_sig_intensity), "\n")
cat("  Frequency:", nrow(fun_sig_frequency), "\n")
cat("  Composition:", nrow(fun_sig_composition), "\n\n")

# NOTE:
# No significant indicator ASVs were detected for Composition (FDR < 0.05).
# Composition results are retained for completeness.

# =========================================================
# 7. Save outputs - Table S10
# =========================================================

# Bacteria
write.csv(bac_sig_intensity, "./results/tables/indicator_asvs_bacteria_intensity.csv",
          row.names = FALSE)
write.csv(bac_sig_frequency, "./results/tables/indicator_asvs_bacteria_frequency.csv",
          row.names = FALSE)

# Fungi
write.csv(fun_sig_intensity, "./results/tables/indicator_asvs_fungi_intensity.csv",
          row.names = FALSE)
write.csv(fun_sig_frequency, "./results/tables/indicator_asvs_fungi_frequency.csv",
          row.names = FALSE)

