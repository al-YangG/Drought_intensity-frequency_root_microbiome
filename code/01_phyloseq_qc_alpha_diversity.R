### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Phyloseq construction, sequencing-depth QC, and alpha diversity ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

## =========================================================
## 1. Helper functions
## =========================================================

# Remove leading X only if followed by a digit
fix_ids <- function(x) sub("^X(?=[0-9])", "", x, perl = TRUE)

# Basic sequencing-depth summary
qc_summary <- function(ps, dataset, cutoffs = c(500, 1000, 2000, 5000, 10000)) {
  r <- sample_sums(ps)
  
  out <- data.frame(Dataset = dataset, Samples = nsamples(ps), ASVs = ntaxa(ps),
                    Total_reads = sum(r), Min = min(r), Median = median(r),
                    Mean = round(mean(r), 1), Max = max(r), stringsAsFactors = FALSE)
  
  for (c in cutoffs) {
    out[[paste0("n<", c)]] <- sum(r < c)
  }
  
  out
}

# Convert phyloseq object to samples x taxa matrix
otu_mat <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# Build long dataframe of rarefaction curves
compute_rare_df_ps <- function(ps, dataset, step = 200) {
  otu <- otu_mat(ps)
  depths <- rowSums(otu)
  
  map_dfr(seq_len(nrow(otu)), function(i) {
    d <- depths[i]
    if (d < step) return(NULL)
    
    grid <- seq(step, d, by = step)
    
    exp_rich <- suppressWarnings(as.numeric(vegan::rarefy(otu[i, , drop = FALSE], sample = grid)))
    
    tibble(Dataset = dataset, SampleID = rownames(otu)[i],
           Depth = grid, ObservedExp = exp_rich)
  })
}

# Count retained and dropped samples at candidate rarefaction depths
keep_table <- function(ps, depths) {
  s <- sample_sums(ps)
  data.frame(depth = depths,
             n_keep = sapply(depths, function(d) sum(s >= d)),
             n_drop = sapply(depths, function(d) sum(s < d)))
}

## =========================================================
## 2. Build phyloseq objects
## =========================================================

## Input data
# ASV table: rows = ASVs, columns = samples
# Taxonomy table: rows = ASVs, columns = taxonomic information
# Metadata table: rows = samples, first column = SampleID

# ---- Bacteria ----
asv_bac <- read.csv("./data/asv_table_Bac.csv", header = TRUE, check.names = FALSE, row.names = 1)
tax_bac <- read.csv("./data/taxonomy_Bac.csv", header = TRUE, check.names = FALSE, row.names = 1)
meta_bac <- read.csv("./data/metadata_Bac.csv", header = TRUE, check.names = FALSE, row.names = 1)

asv_bac[] <- lapply(asv_bac, as.numeric)

OTU_bac  <- otu_table(as.matrix(asv_bac), taxa_are_rows = TRUE)
TAX_bac  <- tax_table(as.matrix(tax_bac))
SAMP_bac <- sample_data(meta_bac)

ps_Bac <- phyloseq(OTU_bac, TAX_bac, SAMP_bac)
ps_Bac <- prune_taxa(taxa_sums(ps_Bac) > 0, ps_Bac)

saveRDS(ps_Bac, "./data/phyloseq/ps_Bac.rds")

# ---- Fungi ----
asv_fun <- read.csv("./data/asv_table_Fun.csv", header = TRUE, check.names = FALSE, row.names = 1)
tax_fun <- read.csv("./data/taxonomy_Fun.csv", header = TRUE, check.names = FALSE, row.names = 1)
meta_fun <- read.csv("./data/metadata_Fun.csv", header = TRUE, check.names = FALSE, row.names = 1)

asv_fun[] <- lapply(asv_fun, as.numeric)

OTU_fun  <- otu_table(as.matrix(asv_fun), taxa_are_rows = TRUE)
TAX_fun  <- tax_table(as.matrix(tax_fun))
SAMP_fun <- sample_data(meta_fun)

ps_Fun <- phyloseq(OTU_fun, TAX_fun, SAMP_fun)
ps_Fun <- prune_taxa(taxa_sums(ps_Fun) > 0, ps_Fun)

saveRDS(ps_Fun, "./data/phyloseq/ps_Fun.rds")

## =========================================================
## 3. Sequencing-depth QC
## =========================================================

# Initial summary
summary_stats <- bind_rows(qc_summary(ps_Bac, "Bacteria"),
                           qc_summary(ps_Fun, "Fungi"))

# Identify bacterial samples with very low depth
low_depth_bac <- sample_sums(ps_Bac)[sample_sums(ps_Bac) < 1000]
low_depth_bac

# Remove low-depth bacterial sample(s)
ps_Bac <- prune_samples(sample_sums(ps_Bac) >= 1000, ps_Bac)
saveRDS(ps_Bac, "./data/phyloseq/ps_Bac_filtered.rds")

# Updated summary after filtering
summary_stats_filtered <- bind_rows(qc_summary(ps_Bac, "Bacteria"),
                                    qc_summary(ps_Fun, "Fungi"))

# Save - Table S2
write.csv(summary_stats_filtered, "./results/tables/qc_summary_bac_fun.csv", row.names = FALSE)

## =========================================================
## 4. Rarefaction assessment
## =========================================================

ps_list <- list(Bacteria = ps_Bac, Fungi = ps_Fun)

# Build rarefaction dataframe
rare_df <- imap_dfr(ps_list, ~compute_rare_df_ps(.x, .y, step = 200))

# Median rarefaction curve per dataset
median_df <- rare_df %>%
  group_by(Dataset) %>%
  group_modify(~{
    upper <- as.numeric(quantile(.x$Depth, 0.95, na.rm = TRUE))
    upper <- max(upper, 1000)
    
    grid <- unique(round(exp(seq(log(1000), log(upper), length.out = 120))))
    by_samp <- group_split(.x, .x$SampleID)
    
    interp <- lapply(by_samp, function(d) {
      d <- arrange(d, Depth)
      if (nrow(d) >= 2) {
        out <- approx(x = d$Depth, y = d$ObservedExp, xout = grid,
                      method = "linear", rule = 2)
        tibble(Depth = out$x, ObservedExp = out$y)
      } else {
        tibble(Depth = grid,
               ObservedExp = ifelse(grid <= max(d$Depth), d$ObservedExp[1], NA_real_))
      }
    })
    
    bind_rows(interp) %>%
      group_by(Depth) %>%
      summarise(MedianObserved = median(ObservedExp, na.rm = TRUE), .groups = "drop")
  }) %>%
  ungroup()

# Keep plotting order
order_vec <- names(ps_list)
rare_df$Dataset <- factor(rare_df$Dataset, levels = order_vec)
median_df$Dataset <- factor(median_df$Dataset, levels = order_vec)

# Chosen rarefaction depths
depth_lines <- tibble(Dataset = factor(c("Bacteria", "Fungi"), levels = order_vec),
                      RareDepth = c(5000, 3000))

# Plot
p_rarefaction <- ggplot(rare_df, aes(Depth, ObservedExp, group = SampleID)) +
  geom_path(alpha = 0.35, linewidth = 0.3, color = "grey40", na.rm = TRUE) +
  geom_smooth(data = median_df, aes(x = Depth, y = MedianObserved),
              inherit.aes = FALSE, method = "loess", span = 0.25, se = FALSE,
              linewidth = 1.2, color = "#2C7FB8") +
  geom_vline(data = depth_lines, aes(xintercept = RareDepth), inherit.aes = FALSE,
             linetype = "dashed", linewidth = 0.6, color = "red3") +
  geom_text(data = depth_lines, aes(x = RareDepth, y = Inf, label = paste0("Depth = ", RareDepth)),
            inherit.aes = FALSE, angle = 90, vjust = 1.2, hjust = 1.1, size = 3, color = "red3") +
  facet_wrap(~Dataset, scales = "free_y", ncol = 2) +
  labs(x = "Sequencing depth (reads)", y = "Observed ASVs (expected)") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

p_rarefaction

# Save - Figure S2 
ggsave(filename = "./results/figures/rarefaction_curves.pdf", plot = p_rarefaction,
       width = 10, height = 5, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/rarefaction_curves.png", plot = p_rarefaction,
       width = 10, height = 5, dpi = 300)


## =========================================================
## 5. Alpha diversity analysis
## =========================================================

# Dataset-specific rarefaction depths
# Bacteria: 5000 reads
# Fungi:    3000 reads

set.seed(34)

ps_Bac_rar <- rarefy_even_depth(ps_Bac, sample.size = 5000, rngseed = 123,
                                replace = FALSE, verbose = TRUE)

ps_Fun_rar <- rarefy_even_depth(ps_Fun, sample.size = 3000, rngseed = 123,
                                replace = FALSE, verbose = TRUE)

# Fix sample IDs in phyloseq objects
sample_names(ps_Bac_rar) <- fix_ids(sample_names(ps_Bac_rar))
sample_names(ps_Fun_rar) <- fix_ids(sample_names(ps_Fun_rar))

sd_b <- data.frame(sample_data(ps_Bac_rar))
sd_f <- data.frame(sample_data(ps_Fun_rar))

rownames(sd_b) <- fix_ids(rownames(sd_b))
rownames(sd_f) <- fix_ids(rownames(sd_f))

sample_data(ps_Bac_rar) <- sample_data(sd_b)
sample_data(ps_Fun_rar) <- sample_data(sd_f)

# Sanity checks
stopifnot(identical(sample_names(ps_Bac_rar), rownames(data.frame(sample_data(ps_Bac_rar)))))
stopifnot(identical(sample_names(ps_Fun_rar), rownames(data.frame(sample_data(ps_Fun_rar)))))

# Calculate alpha diversity
alpha_bac <- estimate_richness(ps_Bac_rar, measures = c("Observed", "Shannon", "Simpson")) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = fix_ids(SampleID))

alpha_fun <- estimate_richness(ps_Fun_rar, measures = c("Observed", "Shannon", "Simpson")) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = fix_ids(SampleID))

# Join metadata
meta_bac_join <- data.frame(sample_data(ps_Bac_rar)) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = fix_ids(SampleID))

meta_fun_join <- data.frame(sample_data(ps_Fun_rar)) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = fix_ids(SampleID))

alpha_bac <- left_join(alpha_bac, meta_bac_join, by = "SampleID")
alpha_fun <- left_join(alpha_fun, meta_fun_join, by = "SampleID")

# Final checks
stopifnot(nrow(alpha_bac) == nsamples(ps_Bac_rar))
stopifnot(nrow(alpha_fun) == nsamples(ps_Fun_rar))
stopifnot(sum(is.na(alpha_bac$Composition)) == 0)
stopifnot(sum(is.na(alpha_fun$Composition)) == 0)

# Save sample-level alpha diversity
write.csv(alpha_bac, "./results/tables/alpha_diversity_bacteria.csv", row.names = FALSE)
write.csv(alpha_fun, "./results/tables/alpha_diversity_fungi.csv", row.names = FALSE)

# Pot-level summaries
alpha_bac_pot <- alpha_bac %>%
  mutate(Pot = as.character(Pot)) %>%
  group_by(Pot) %>%
  summarise(Intensity = first(Intensity),
            Frequency = first(Frequency),
            Composition = first(Composition),
            Observed_mean = mean(Observed, na.rm = TRUE),
            Shannon_mean = mean(Shannon, na.rm = TRUE),
            Simpson_mean = mean(Simpson, na.rm = TRUE),
            n_samples = n(),
            .groups = "drop")

alpha_fun_pot <- alpha_fun %>%
  mutate(Pot = as.character(Pot)) %>%
  group_by(Pot) %>%
  summarise(Intensity = first(Intensity),
            Frequency = first(Frequency),
            Composition = first(Composition),
            Observed_mean = mean(Observed, na.rm = TRUE),
            Shannon_mean = mean(Shannon, na.rm = TRUE),
            Simpson_mean = mean(Simpson, na.rm = TRUE),
            n_samples = n(),
            .groups = "drop")

write.csv(alpha_bac_pot, "./results/tables/alpha_diversity_bacteria_pot_level.csv", row.names = FALSE)
write.csv(alpha_fun_pot, "./results/tables/alpha_diversity_fungi_pot_level.csv", row.names = FALSE)
