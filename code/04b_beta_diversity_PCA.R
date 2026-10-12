### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Beta diversity, PCA, and PERMANOVA ##

# Change to working directory

setwd("~/Desktop/DIF-Microbes/Data availability")

# Relative-abundance CLR-PCA / Aitchison sensitivity analysis for bacterial and fungal ASV tables.
#
# Sample-level analysis matching the structure of the Bray-Curtis PCoA.

library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)

## =========================================================
## 1. Helper functions
## =========================================================

prepare_metadata <- function(ps) {
  
  data.frame(sample_data(ps), check.names = FALSE) %>%
    rownames_to_column("SampleID") %>%
    mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
           Frequency = factor(Frequency, levels = c("Low", "High")),
           Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
           Pot = factor(Pot))
}


# Extract the sample-by-ASV count matrix and align it
# with the sample metadata.

extract_sample_data <- function(ps) {
  
  meta <- prepare_metadata(ps)
  
  counts <- as(otu_table(ps), "matrix")
  
  if (taxa_are_rows(ps)) {
    counts <- t(counts)
  }
  
  counts <- counts[meta$SampleID, , drop = FALSE]
  
  if (!identical(rownames(counts), meta$SampleID)) {
    stop("Sample order differs between the ASV table and metadata.")
  }
  
  list(counts = counts, metadata = meta)
}


# Convert raw counts to sample-wise relative abundances,
# add a small pseudocount to accommodate zeros, and apply
# the centred log-ratio transformation.

clr_transform_relative <- function(count_matrix, pseudocount = 1e-6) {
  
  count_matrix <- as.matrix(count_matrix)
  
  if (any(!is.finite(count_matrix))) {
    stop("The count matrix contains non-finite values.")
  }
  
  if (any(count_matrix < 0)) {
    stop("Counts must be non-negative.")
  }
  
  if (!is.numeric(pseudocount) ||
      length(pseudocount) != 1L ||
      !is.finite(pseudocount) ||
      pseudocount <= 0) {
    stop("The pseudocount must be a single positive finite number.")
  }
  
  sample_totals <- rowSums(count_matrix)
  
  if (any(sample_totals <= 0)) {
    zero_samples <- rownames(count_matrix)[sample_totals <= 0]
    
    stop(paste0("Samples with zero total abundance cannot be transformed: ",
                paste(zero_samples, collapse = ", ")))
  }
  
  # Convert counts to relative abundances.
  
  relative_abundance <- sweep(count_matrix, 1, sample_totals,FUN = "/")
  
  # Add the pseudocount and re-close each sample to sum to one.
  #
  # Re-closing is not required mathematically for CLR because
  # the CLR transformation is scale invariant, but it makes the
  # zero-replacement procedure explicit.
  
  relative_adjusted <- relative_abundance + pseudocount
  
  relative_adjusted <- sweep(relative_adjusted, 1, rowSums(relative_adjusted), FUN = "/")
  
  # Apply the CLR transformation.
  
  log_abundance <- log(relative_adjusted)
  
  clr <- sweep(log_abundance, 1, rowMeans(log_abundance), FUN = "-")
  
  list(relative_abundance = relative_abundance, relative_adjusted = relative_adjusted, clr = clr)
}


run_clr_pca <- function(ps, dataset_name, pseudocount = 1e-6, flip_axis1 = FALSE) {
  
  sample_data_list <- extract_sample_data(ps)
  
  transformed <- clr_transform_relative(sample_data_list$counts, pseudocount = pseudocount)
  
  clr <- transformed$clr
  
  # Remove invariant ASVs before PCA. Invariant ASVs do not
  # contribute to distances between samples.
  
  asv_sd <- apply(clr, 2, sd)
  keep <- is.finite(asv_sd) & asv_sd > 0
  
  if (sum(keep) < 2) {
    stop(paste0(dataset_name, ": fewer than two non-invariant ASVs remain for PCA."))
  }
  
  clr_pca <- clr[, keep, drop = FALSE]
  
  pca <- prcomp(clr_pca, center = TRUE, scale. = FALSE)
  
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  
  colnames(scores) <- c("PC1", "PC2")
  
  if (flip_axis1) {
    scores$PC1 <- -scores$PC1
  }
  
  scores <- scores %>%
    rownames_to_column("SampleID") %>%
    left_join(sample_data_list$metadata, by = "SampleID") %>%
    mutate(Dataset = dataset_name)
  
  list(scores = scores, relative_abundance = transformed$relative_abundance,
       relative_adjusted = transformed$relative_adjusted, clr = clr,
       distance = dist(clr, method = "euclidean"), pca = pca,
       variance = c(PC1 = variance[1], PC2 = variance[2]),
       metadata = sample_data_list$metadata, pseudocount = pseudocount)
}


run_permanova <- function(dist_obj, meta_df, strata = NULL, permutations = 999) {
  
  if (is.null(strata)) {
    
    adonis2(dist_obj ~ Intensity * Frequency * Composition, data = meta_df,
            permutations = permutations, by = "term")
    
  } else {
    
    adonis2(dist_obj ~ Intensity * Frequency * Composition, data = meta_df,
            permutations = permutations, by = "term", strata = strata)
  }
}


format_adonis_table <- function(adonis_obj) {
  
  out <- as.data.frame(adonis_obj)
  out$Term <- rownames(out)
  
  out[, c("Term", setdiff(names(out), "Term"))]
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

dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Relative-abundance CLR transformation and PCA
## =========================================================

# Axis signs are arbitrary. Flip an axis only if needed for
# visual comparison with the Bray-Curtis PCoA. Axis flipping
# does not affect distances or statistical tests.

bac_clr <- run_clr_pca(ps = ps_Bac, dataset_name = "Bacteria", pseudocount = 1e-6, flip_axis1 = FALSE)

fun_clr <- run_clr_pca(ps = ps_Fun, dataset_name = "Fungi", pseudocount = 1e-6, flip_axis1 = FALSE)


## =========================================================
## 4. Sample-level PERMANOVA using Aitchison distance
## =========================================================

# Unrestricted analyses

set.seed(807)

perm_bac_clr <- run_permanova(bac_clr$distance, bac_clr$metadata, permutations = 999)

set.seed(807)

perm_fun_clr <- run_permanova(fun_clr$distance, fun_clr$metadata, permutations = 999)


write.csv(format_adonis_table(perm_bac_clr), row.names = FALSE,
          paste0("./results/tables/", "permanova_aitchison_bacteria_sample_level.csv"))

write.csv(format_adonis_table(perm_fun_clr), row.names = FALSE,
          paste0("./results/tables/", "permanova_aitchison_fungi_sample_level.csv"))


## =========================================================
## 5. Save variance explained and analysis settings
## =========================================================

pca_variance <- bind_rows(data.frame(Dataset = "Bacteria", Axis = c("PC1", "PC2"),
                                     Variance_explained_percent = round(unname(bac_clr$variance), 1)),
                          data.frame(Dataset = "Fungi", Axis = c("PC1", "PC2"),
                                     Variance_explained_percent = round(unname(fun_clr$variance), 1)))

write.csv(pca_variance, paste0("./results/tables/", "clr_pca_variance_explained_sample_level.csv"), row.names = FALSE)


analysis_settings <- data.frame(Setting = c("analysis_level", "subsample_aggregation",
                                            "input_abundance", "zero_replacement",
                                            "post_replacement_closure", "ordination",
                                            "distance", "permutations", "permanova_tests"),
                                Value = c("Sample", "None",
                                          "Relative abundance calculated from raw ASV counts",
                                          "pseudocount = 1e-6 added to relative abundances",
                                          "Each adjusted sample re-closed to sum to one",
                                          "PCA of CLR-transformed relative abundances",
                                          "Euclidean distance on CLR values (Aitchison distance)",
                                          "999", "Sequential tests (by = term)"),
                                stringsAsFactors = FALSE)

write.csv(analysis_settings, "./results/tables/clr_pca_analysis_settings.csv", row.names = FALSE)
