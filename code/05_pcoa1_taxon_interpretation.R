### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: PCoA1 interpretation via pot-level taxon correlations ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(stringr)

## =========================================================
## 1. Helper functions
## =========================================================

# Ensure OTU table is numeric without losing taxa names
make_numeric_ps <- function(ps) {
  mat <- as(otu_table(ps), "matrix")
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  otu_table(ps) <- otu_table(mat, taxa_are_rows = taxa_are_rows(ps))
  ps
}

# Build one row of metadata per pot
make_pot_metadata_df <- function(ps) {
  md <- data.frame(sample_data(ps))
  stopifnot(all(c("Pot", "Intensity", "Frequency", "Composition") %in% colnames(md)))
  
  md$Pot <- as.character(md$Pot)
  
  pot_md <- md %>%
    group_by(Pot) %>%
    summarise(Intensity   = first(Intensity),
              Frequency   = first(Frequency),
              Composition = first(Composition),
              .groups = "drop") %>%
    as.data.frame()
  
  rownames(pot_md) <- pot_md$Pot
  pot_md$Pot <- NULL
  
  pot_md$Intensity   <- factor(pot_md$Intensity, levels = c("Dry", "Wet"))
  pot_md$Frequency   <- factor(pot_md$Frequency, levels = c("Low", "High"))
  pot_md$Composition <- factor(pot_md$Composition, levels = c("Conspecific", "Heterospecific"))
  
  pot_md
}

# Pot-level PCoA extraction
extract_pcoa_axes_pot <- function(ps, dataset_name, n_axes = 2, flip_axis1 = FALSE) {
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- make_numeric_ps(ps)
  
  pot_md <- make_pot_metadata_df(ps)
  
  # Merge counts by pot
  ps_pot <- merge_samples(ps, "Pot")
  
  # Restore clean pot metadata
  pots <- sample_names(ps_pot)
  pot_md <- pot_md[pots, , drop = FALSE]
  stopifnot(identical(rownames(pot_md), pots))
  
  sample_data(ps_pot) <- sample_data(pot_md)
  
  # Relative abundance after merging
  ps_pot <- transform_sample_counts(ps_pot, function(x) x / sum(x))
  
  # Bray-Curtis distance + PCoA
  dist_bc <- phyloseq::distance(ps_pot, method = "bray")
  ord <- cmdscale(dist_bc, k = n_axes, eig = TRUE)
  
  axes_df <- as.data.frame(ord$points[, seq_len(n_axes), drop = FALSE])
  colnames(axes_df) <- paste0("PCoA", seq_len(n_axes))
  
  if (flip_axis1) {
    axes_df$PCoA1 <- -axes_df$PCoA1
  }
  
  eig_vals <- ord$eig
  var_exp <- 100 * eig_vals / sum(eig_vals[eig_vals > 0], na.rm = TRUE)
  
  axes_df$Pot <- rownames(axes_df)
  axes_df <- axes_df %>%
    mutate(Intensity   = pot_md[Pot, "Intensity"],
           Frequency   = pot_md[Pot, "Frequency"],
           Composition = pot_md[Pot, "Composition"],
           Dataset     = dataset_name) %>%
    relocate(Dataset, Pot, Intensity, Frequency, Composition)
  
  list(ps_pot = ps_pot, dist_bc = dist_bc, ord = ord,
       axes_df = axes_df, var_exp = var_exp[seq_len(n_axes)])
}

# Correlate pot-level taxon relative abundances with PCoA1
interpret_pcoa1_ends <- function(ps, axes_df, axis_col = "PCoA1", rho_cutoff = 0.7,
                                 low_label = "Low_PCoA1", high_label = "High_PCoA1") {
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- make_numeric_ps(ps)
  
  # Ensure Pot exists
  if (!"Pot" %in% colnames(data.frame(sample_data(ps)))) {
    sample_data(ps)$Pot <- paste0("Pot-", str_extract(sample_names(ps), "^[0-9]+"))
  }
  
  # Relative abundance before merging
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_pot <- merge_samples(ps_rel, "Pot")
  
  otu <- as.data.frame(otu_table(ps_pot))
  if (!taxa_are_rows(ps_pot)) {
    otu <- t(otu)
  }
  otu <- as.data.frame(otu)
  
  stopifnot("Pot" %in% colnames(axes_df))
  common_pots <- intersect(colnames(otu), axes_df$Pot)
  
  otu2 <- otu[, common_pots, drop = FALSE]
  meta2 <- axes_df %>%
    filter(Pot %in% common_pots) %>%
    as.data.frame()
  
  rownames(meta2) <- meta2$Pot
  meta2 <- meta2[common_pots, , drop = FALSE]
  
  stopifnot(all(colnames(otu2) == rownames(meta2)))
  
  taxa_cor <- apply(otu2, 1, function(x) {
    cor(x, meta2[[axis_col]], method = "spearman", use = "pairwise.complete.obs")
  })
  
  tax <- as.data.frame(tax_table(ps_pot)) %>%
    rownames_to_column("ASV")
  
  cor_df <- data.frame(ASV = names(taxa_cor), rho = as.numeric(taxa_cor),
                       stringsAsFactors = FALSE) %>%
    left_join(tax, by = "ASV")
  
  strong_taxa <- cor_df %>%
    filter(!is.na(rho), abs(rho) > rho_cutoff)
  
  pos_taxa <- strong_taxa %>%
    filter(rho > 0) %>%
    mutate(PCoA1_end = high_label)
  
  neg_taxa <- strong_taxa %>%
    filter(rho < 0) %>%
    mutate(PCoA1_end = low_label)
  
  strong_both <- bind_rows(pos_taxa, neg_taxa)
  
  list(ps_pot = ps_pot, otu = otu2, meta = meta2, cor_df = cor_df,
       strong_taxa = strong_taxa, strong_both = strong_both)
}

## =========================================================
## 2. Load phyloseq objects
## =========================================================

ps_Bac <- readRDS("./data/phyloseq/ps_Bac_filtered.rds")
ps_Fun <- readRDS("./data/phyloseq/ps_Fun.rds")

## =========================================================
## 3. Pot-level PCoA
## =========================================================

# Fungal PCoA1 is flipped for consistency with the main ordination figure
bac_pcoa <- extract_pcoa_axes_pot(ps = ps_Bac, dataset_name = "Bacteria",
                                  n_axes = 2, flip_axis1 = FALSE)

fun_pcoa <- extract_pcoa_axes_pot(ps = ps_Fun, dataset_name = "Fungi",
                                  n_axes = 2, flip_axis1 = TRUE)

write.csv(bac_pcoa$axes_df, "./results/tables/pcoa_axes_bacteria_pot_level.csv",
          row.names = FALSE)

write.csv(fun_pcoa$axes_df, "./results/tables/pcoa_axes_fungi_pot_level.csv",
          row.names = FALSE)

pcoa_variance_pot <- bind_rows(data.frame(Dataset = "Bacteria", Axis = c("PCoA1", "PCoA2"),
                                          Variance_explained_percent = round(bac_pcoa$var_exp[1:2], 3)),
                               data.frame(Dataset = "Fungi", Axis = c("PCoA1", "PCoA2"),
                                          Variance_explained_percent = round(fun_pcoa$var_exp[1:2], 3)))

write.csv(pcoa_variance_pot, "./results/tables/pcoa_variance_explained_pot_level.csv",
          row.names = FALSE)

## =========================================================
## 4. PCoA1 interpretation via taxon correlations
## =========================================================

bac_pcoa1_interp <- interpret_pcoa1_ends(ps = ps_Bac, axes_df = bac_pcoa$axes_df,
                                         axis_col = "PCoA1", rho_cutoff = 0.7)

fun_pcoa1_interp <- interpret_pcoa1_ends(ps = ps_Fun, axes_df = fun_pcoa$axes_df,
                                         axis_col = "PCoA1", rho_cutoff = 0.6)

## =========================================================
## 5. Save full and filtered correlation tables
## =========================================================

write.csv(bac_pcoa1_interp$cor_df, "./results/tables/pcoa1_taxon_correlations_bacteria.csv",
          row.names = FALSE)
write.csv(fun_pcoa1_interp$cor_df, "./results/tables/pcoa1_taxon_correlations_fungi.csv",
          row.names = FALSE)

write.csv(bac_pcoa1_interp$strong_both, "./results/tables/pcoa1_associated_asvs_bacteria.csv",
          row.names = FALSE)
write.csv(fun_pcoa1_interp$strong_both, "./results/tables/pcoa1_associated_asvs_fungi.csv",
          row.names = FALSE)

## Combined table for downstream figure summaries
tab_bac <- bac_pcoa1_interp$strong_both %>%
  mutate(Domain = "Bacteria")

tab_fun <- fun_pcoa1_interp$strong_both %>%
  mutate(Domain = "Fungi")

tab_pcoa1_asvs <- bind_rows(tab_bac, tab_fun) %>%
  dplyr::select(Domain, PCoA1_end, ASV, rho, everything()) %>%
  arrange(Domain, PCoA1_end, desc(abs(rho)))

# Save - Table S6
write.csv(tab_pcoa1_asvs, "./results/tables/pcoa1_associated_asvs_bacteria_fungi.csv",
          row.names = FALSE)

## =========================================================
## 6. Save simple summaries
## =========================================================

summary_counts <- bind_rows(bac_pcoa1_interp$strong_both %>%
                              count(PCoA1_end, name = "n_asvs") %>%
                              mutate(Dataset = "Bacteria"),
                            fun_pcoa1_interp$strong_both %>%
                              count(PCoA1_end, name = "n_asvs") %>%
                              mutate(Dataset = "Fungi")) %>%
  dplyr::select(Dataset, PCoA1_end, n_asvs)

write.csv(summary_counts, "./results/tables/pcoa1_associated_asv_counts_summary.csv",
          row.names = FALSE)
