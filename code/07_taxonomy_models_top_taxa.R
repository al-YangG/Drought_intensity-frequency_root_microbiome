### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Taxonomy models and top-taxa data preparation ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(phyloseq)
library(dplyr)
library(purrr)
library(tidyr)
library(glmmTMB)
library(car)
library(tibble)

## =========================================================
## 1. Helper functions
## =========================================================

clean_tax_table <- function(ps) {
  tax <- as.data.frame(tax_table(ps))
  
  tax[] <- lapply(tax, function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "" | x == "NA"] <- "Unclassified"
    x
  })
  
  tax_table(ps) <- tax_table(as.matrix(tax))
  ps
}

top_taxa_table <- function(df, tax_rank, top_n) {
  df %>%
    filter(.data[[tax_rank]] != "Unclassified") %>%
    group_by(.data[[tax_rank]]) %>%
    summarise(mean_abundance = mean(Abundance),
              total_abundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(mean_abundance)) %>%
    slice_head(n = top_n)
}

make_heatmap_data <- function(df, tax_rank, top_taxa) {
  df_sub <- df %>%
    filter(.data[[tax_rank]] %in% top_taxa)
  
  df_intensity <- df_sub %>%
    group_by(Level = Intensity, .data[[tax_rank]]) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Effect = "Intensity")
  
  df_frequency <- df_sub %>%
    group_by(Level = Frequency, .data[[tax_rank]]) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Effect = "Frequency")
  
  df_composition <- df_sub %>%
    group_by(Level = Composition, .data[[tax_rank]]) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Effect = "Composition")
  
  bind_rows(df_intensity, df_frequency, df_composition)
}

prepare_barplot_data <- function(df, tax_rank, top_taxa) {
  df_plot <- df %>%
    mutate(Taxon = ifelse(.data[[tax_rank]] %in% top_taxa, .data[[tax_rank]], "Other"))
  
  df_intensity <- df_plot %>%
    group_by(Level = Intensity, Taxon) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Effect = "Intensity")
  
  df_frequency <- df_plot %>%
    group_by(Level = Frequency, Taxon) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Effect = "Frequency")
  
  df_composition <- df_plot %>%
    group_by(Level = Composition, Taxon) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Effect = "Composition")
  
  bind_rows(df_intensity, df_frequency, df_composition) %>%
    group_by(Effect, Level) %>%
    mutate(mean_abundance = mean_abundance / sum(mean_abundance)) %>%
    ungroup() %>%
    mutate(
      Effect = factor(Effect, levels = c("Intensity", "Frequency", "Composition")),
      Level = factor(Level, levels = c("Dry", "Wet", "Low", "High", "Conspecific", "Heterospecific"))
    )
}

beta_transform <- function(x) {
  n <- length(x)
  (x * (n - 1) + 0.5) / n
}

run_taxon_models <- function(df, taxon_col) {
  taxa <- unique(df[[taxon_col]])
  
  map_dfr(taxa, function(taxon) {
    dat <- df %>% filter(.data[[taxon_col]] == taxon)
    
    mod <- glmmTMB(Abundance_beta ~ Intensity * Frequency * Composition + (1 | Pot),
                   data = dat, family = beta_family(link = "logit"))
    
    an <- car::Anova(mod, type = 3)
    out <- as.data.frame(an)
    out$Term <- rownames(out)
    out$Taxon <- taxon
    rownames(out) <- NULL
    out
  })
}

order_heatmap_taxa <- function(df, tax_rank) {
  df %>%
    group_by(.data[[tax_rank]]) %>%
    summarise(total = mean(mean_abundance), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(.data[[tax_rank]])
}

order_barplot_taxa <- function(df) {
  df %>%
    group_by(Taxon) %>%
    summarise(total = mean(mean_abundance), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(Taxon)
}

## =========================================================
## 2. Load phyloseq objects
## =========================================================

ps_Bac <- readRDS("./data/phyloseq/ps_Bac_filtered.rds")
ps_Fun <- readRDS("./data/phyloseq/ps_Fun.rds")

## =========================================================
## 3. Clean taxonomy and agglomerate
## =========================================================

ps_Bac <- clean_tax_table(ps_Bac)
ps_Fun <- clean_tax_table(ps_Fun)

ps_Bac_phylum <- tax_glom(ps_Bac, taxrank = "Phylum", NArm = FALSE)
ps_Fun_phylum <- tax_glom(ps_Fun, taxrank = "Phylum", NArm = FALSE)

ps_Bac_genus <- tax_glom(ps_Bac, taxrank = "Genus", NArm = FALSE)
ps_Fun_genus <- tax_glom(ps_Fun, taxrank = "Genus", NArm = FALSE)

## =========================================================
## 4. Relative abundance and long-format tables
## =========================================================

ps_Bac_phylum_rel <- transform_sample_counts(ps_Bac_phylum, function(x) x / sum(x))
ps_Fun_phylum_rel <- transform_sample_counts(ps_Fun_phylum, function(x) x / sum(x))

ps_Bac_genus_rel <- transform_sample_counts(ps_Bac_genus, function(x) x / sum(x))
ps_Fun_genus_rel <- transform_sample_counts(ps_Fun_genus, function(x) x / sum(x))

bac_phylum_df <- psmelt(ps_Bac_phylum_rel)
fun_phylum_df <- psmelt(ps_Fun_phylum_rel)

bac_genus_df <- psmelt(ps_Bac_genus_rel)
fun_genus_df <- psmelt(ps_Fun_genus_rel)

# Optional label cleanup used in plotting
fun_genus_df$Genus <- gsub("Stachybotryaceae_gen_Incertae_sedis",
                           "Stachybotryaceae sp.", fun_genus_df$Genus)

write.csv(bac_phylum_df, "./results/tables/taxonomy_bacteria_phylum_long.csv", row.names = FALSE)
write.csv(fun_phylum_df, "./results/tables/taxonomy_fungi_phylum_long.csv", row.names = FALSE)
write.csv(bac_genus_df, "./results/tables/taxonomy_bacteria_genus_long.csv", row.names = FALSE)
write.csv(fun_genus_df, "./results/tables/taxonomy_fungi_genus_long.csv", row.names = FALSE)

## =========================================================
## 5. Identify top taxa
## =========================================================

# Top 5 for heatmaps and beta-regression
top5_bac_phyla <- top_taxa_table(bac_phylum_df, "Phylum", 5)$Phylum
top5_fun_phyla <- top_taxa_table(fun_phylum_df, "Phylum", 5)$Phylum
top5_bac_genera <- top_taxa_table(bac_genus_df, "Genus", 5)$Genus
top5_fun_genera <- top_taxa_table(fun_genus_df, "Genus", 5)$Genus

# Top 10 for barplots
top10_bac_phyla <- top_taxa_table(bac_phylum_df, "Phylum", 10)$Phylum
top10_fun_phyla <- top_taxa_table(fun_phylum_df, "Phylum", 10)$Phylum
top10_bac_genera <- top_taxa_table(bac_genus_df, "Genus", 10)$Genus
top10_fun_genera <- top_taxa_table(fun_genus_df, "Genus", 10)$Genus

top_taxa_summary <- bind_rows(
  data.frame(Dataset = "Bacteria", Rank = "Phylum", Purpose = "Heatmap", TopN = 5, Taxon = top5_bac_phyla),
  data.frame(Dataset = "Fungi",    Rank = "Phylum", Purpose = "Heatmap", TopN = 5, Taxon = top5_fun_phyla),
  data.frame(Dataset = "Bacteria", Rank = "Genus",  Purpose = "Heatmap", TopN = 5, Taxon = top5_bac_genera),
  data.frame(Dataset = "Fungi",    Rank = "Genus",  Purpose = "Heatmap", TopN = 5, Taxon = top5_fun_genera),
  data.frame(Dataset = "Bacteria", Rank = "Phylum", Purpose = "Barplot", TopN = 10, Taxon = top10_bac_phyla),
  data.frame(Dataset = "Fungi",    Rank = "Phylum", Purpose = "Barplot", TopN = 10, Taxon = top10_fun_phyla),
  data.frame(Dataset = "Bacteria", Rank = "Genus",  Purpose = "Barplot", TopN = 10, Taxon = top10_bac_genera),
  data.frame(Dataset = "Fungi",    Rank = "Genus",  Purpose = "Barplot", TopN = 10, Taxon = top10_fun_genera))

write.csv(top_taxa_summary, "./results/tables/top_taxa_summary.csv", row.names = FALSE)

## =========================================================
## 6. Prepare heatmap input tables
## =========================================================

effect_levels <- c("Intensity", "Frequency", "Composition")
level_order <- c("Dry", "Wet", "Low", "High", "Conspecific", "Heterospecific")

bac_phylum_heat <- make_heatmap_data(bac_phylum_df, "Phylum", top5_bac_phyla)
fun_phylum_heat <- make_heatmap_data(fun_phylum_df, "Phylum", top5_fun_phyla)
bac_genus_heat  <- make_heatmap_data(bac_genus_df,  "Genus",  top5_bac_genera)
fun_genus_heat  <- make_heatmap_data(fun_genus_df,  "Genus",  top5_fun_genera)

bac_phylum_order <- order_heatmap_taxa(bac_phylum_heat, "Phylum")
fun_phylum_order <- order_heatmap_taxa(fun_phylum_heat, "Phylum")
bac_genus_order  <- order_heatmap_taxa(bac_genus_heat,  "Genus")
fun_genus_order  <- order_heatmap_taxa(fun_genus_heat,  "Genus")

bac_phylum_heat <- bac_phylum_heat %>%
  mutate(Phylum = factor(Phylum, levels = rev(bac_phylum_order)),
         Effect = factor(Effect, levels = effect_levels),
         Level  = factor(Level, levels = level_order))

fun_phylum_heat <- fun_phylum_heat %>%
  mutate(Phylum = factor(Phylum, levels = rev(fun_phylum_order)),
         Effect = factor(Effect, levels = effect_levels),
         Level  = factor(Level, levels = level_order))

bac_genus_heat <- bac_genus_heat %>%
  mutate(Genus  = factor(Genus, levels = rev(bac_genus_order)),
         Effect = factor(Effect, levels = effect_levels),
         Level  = factor(Level, levels = level_order))

fun_genus_heat <- fun_genus_heat %>%
  mutate(Genus  = factor(Genus, levels = rev(fun_genus_order)),
         Effect = factor(Effect, levels = effect_levels),
         Level  = factor(Level, levels = level_order))

write.csv(bac_phylum_heat, "./results/tables/taxonomy_heatmap_bacteria_phylum.csv", row.names = FALSE)
write.csv(fun_phylum_heat, "./results/tables/taxonomy_heatmap_fungi_phylum.csv", row.names = FALSE)
write.csv(bac_genus_heat,  "./results/tables/taxonomy_heatmap_bacteria_genus.csv", row.names = FALSE)
write.csv(fun_genus_heat,  "./results/tables/taxonomy_heatmap_fungi_genus.csv", row.names = FALSE)

## =========================================================
## 7. Prepare barplot input tables
## =========================================================

bac_phylum_bar <- prepare_barplot_data(bac_phylum_df, "Phylum", top10_bac_phyla)
fun_phylum_bar <- prepare_barplot_data(fun_phylum_df, "Phylum", top10_fun_phyla)
bac_genus_bar  <- prepare_barplot_data(bac_genus_df, "Genus", top10_bac_genera)
fun_genus_bar  <- prepare_barplot_data(fun_genus_df, "Genus", top10_fun_genera)

bac_phylum_bar$Taxon <- factor(bac_phylum_bar$Taxon, levels = order_barplot_taxa(bac_phylum_bar))
fun_phylum_bar$Taxon <- factor(fun_phylum_bar$Taxon, levels = order_barplot_taxa(fun_phylum_bar))
bac_genus_bar$Taxon  <- factor(bac_genus_bar$Taxon,  levels = order_barplot_taxa(bac_genus_bar))
fun_genus_bar$Taxon  <- factor(fun_genus_bar$Taxon,  levels = order_barplot_taxa(fun_genus_bar))

write.csv(bac_phylum_bar, "./results/tables/taxonomy_barplot_bacteria_phylum.csv", row.names = FALSE)
write.csv(fun_phylum_bar, "./results/tables/taxonomy_barplot_fungi_phylum.csv", row.names = FALSE)
write.csv(bac_genus_bar,  "./results/tables/taxonomy_barplot_bacteria_genus.csv", row.names = FALSE)
write.csv(fun_genus_bar,  "./results/tables/taxonomy_barplot_fungi_genus.csv", row.names = FALSE)

## =========================================================
## 8. Beta-regression models for top 5 taxa
## =========================================================

options(contrasts = c("contr.sum", "contr.poly"))

# Bacterial phyla
bac_phylum_test <- bac_phylum_df %>%
  filter(Phylum %in% top5_bac_phyla) %>%
  mutate(Abundance_beta = beta_transform(Abundance))

bac_phylum_results <- run_taxon_models(bac_phylum_test, "Phylum") %>%
  group_by(Term) %>%
  mutate(p_adj = p.adjust(`Pr(>Chisq)`, method = "BH")) %>%
  ungroup() %>%
  mutate(Rank = "Bacterial phylum")

# Fungal phyla
fun_phylum_test <- fun_phylum_df %>%
  filter(Phylum %in% top5_fun_phyla) %>%
  mutate(Abundance_beta = beta_transform(Abundance))

fun_phylum_results <- run_taxon_models(fun_phylum_test, "Phylum") %>%
  group_by(Term) %>%
  mutate(p_adj = p.adjust(`Pr(>Chisq)`, method = "BH")) %>%
  ungroup() %>%
  mutate(Rank = "Fungal phylum")

# Bacterial genera
bac_genus_test <- bac_genus_df %>%
  filter(Genus %in% top5_bac_genera) %>%
  mutate(Abundance_beta = beta_transform(Abundance))

bac_genus_results <- run_taxon_models(bac_genus_test, "Genus") %>%
  group_by(Term) %>%
  mutate(p_adj = p.adjust(`Pr(>Chisq)`, method = "BH")) %>%
  ungroup() %>%
  mutate(Rank = "Bacterial genus")

# Fungal genera
fun_genus_test <- fun_genus_df %>%
  filter(Genus %in% top5_fun_genera) %>%
  mutate(Abundance_beta = beta_transform(Abundance))

fun_genus_results <- run_taxon_models(fun_genus_test, "Genus") %>%
  group_by(Term) %>%
  mutate(p_adj = p.adjust(`Pr(>Chisq)`, method = "BH")) %>%
  ungroup() %>%
  mutate(Rank = "Fungal genus")

taxa_results <- bind_rows(bac_phylum_results, fun_phylum_results,
                          bac_genus_results, fun_genus_results) %>%
  dplyr::select(Rank, Taxon, Term, everything())

# Save - Table S7
write.csv(taxa_results, "./results/tables/taxonomy_top5_taxa_beta_models.csv", row.names = FALSE)
