### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: FUNGuild trophic mode annotation and models ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(tidyr)
library(stringr)
library(phyloseq)
library(glmmTMB)
library(car)
library(tibble)

## =========================================================
## 1. Load input data
## =========================================================

funguild <- read.csv("./data/funguild_assignments.csv", stringsAsFactors = FALSE)

ps_Fun <- readRDS("./data/phyloseq/ps_Fun.rds")

## =========================================================
## 2. Clean FUNGuild annotations
## =========================================================

funguild_clean <- funguild %>%
  mutate(trophicMode = str_trim(trophicMode),
         guild = str_trim(guild)) %>%
  filter(!is.na(guild), guild != "") %>%
  filter(confidenceRanking %in% c("Probable", "Highly Probable")) %>%
  mutate(trophicMode_clean = case_when(is.na(trophicMode) ~ NA_character_,
                                       str_detect(trophicMode, "-") ~ "Multi-trophic",
                                       TRUE ~ trophicMode),
         guild_clean = str_split(guild, "-", simplify = TRUE)[, 1],
         guild_clean = str_remove_all(guild_clean, "\\|"),
         guild_clean = str_trim(guild_clean)) %>%
  filter(!str_detect(guild_clean, "Undefined")) %>%
  dplyr::select(OTU, trophicMode, trophicMode_clean, guild, guild_clean, confidenceRanking) %>%
  distinct(OTU, .keep_all = TRUE)

write.csv(funguild_clean, "./results/tables/funguild_cleaned.csv", row.names = FALSE)

## =========================================================
## 3. Match FUNGuild annotations to fungal ASVs
## =========================================================

asv_phyloseq <- taxa_names(ps_Fun)
funguild_asv <- funguild_clean$OTU

funguild_match_summary <- data.frame(n_asvs_phyloseq = length(asv_phyloseq),
                                     n_asvs_funguild_clean = length(funguild_asv),
                                     n_matched_asvs = sum(asv_phyloseq %in% funguild_asv),
                                     prop_matched_asvs = sum(asv_phyloseq %in% funguild_asv) / length(asv_phyloseq))

write.csv(funguild_match_summary, "./results/tables/funguild_match_summary.csv", row.names = FALSE)

funguild_matched <- funguild_clean %>%
  filter(OTU %in% asv_phyloseq)

## =========================================================
## 4. Attach FUNGuild annotations to fungal phyloseq object
## =========================================================

funguild_df <- funguild_matched %>%
  dplyr::select(OTU, trophicMode_clean, guild_clean)

tax <- as.data.frame(tax_table(ps_Fun))
tax$ASV <- rownames(tax)

tax2 <- left_join(tax, funguild_df, by = c("ASV" = "OTU"))

rownames(tax2) <- tax2$ASV
tax2$ASV <- NULL

tax_table(ps_Fun) <- tax_table(as.matrix(tax2))

saveRDS(ps_Fun, "./data/phyloseq/ps_Fun_guild.rds")

guild_annotation_summary <- bind_rows(
  as.data.frame(table(tax_table(ps_Fun)[, "trophicMode_clean"], useNA = "always")) %>%
    mutate(Type = "trophicMode_clean") %>%
    rename(Category = Var1, Count = Freq),
  as.data.frame(table(tax_table(ps_Fun)[, "guild_clean"], useNA = "always")) %>%
    mutate(Type = "guild_clean") %>%
    rename(Category = Var1, Count = Freq))

write.csv(guild_annotation_summary, "./results/tables/funguild_annotation_summary.csv",
          row.names = FALSE)

## =========================================================
## 5. Build sample-level trophic mode abundance tables
## =========================================================

otu <- as(otu_table(ps_Fun), "matrix")
if (!taxa_are_rows(ps_Fun)) {
  otu <- t(otu)
}
otu <- as.data.frame(otu)
otu$ASV <- rownames(otu)

tax <- as.data.frame(tax_table(ps_Fun))
tax$ASV <- rownames(tax)

meta <- data.frame(sample_data(ps_Fun))
meta$SampleID <- rownames(meta)

otu_trophic <- otu %>%
  left_join(tax %>%
              dplyr::select(ASV, trophicMode_clean), by = "ASV") %>%
  mutate(trophicMode_clean = ifelse(is.na(trophicMode_clean), "Unassigned", trophicMode_clean))

trophic_counts <- otu_trophic %>%
  pivot_longer(cols = -c(ASV, trophicMode_clean), names_to = "SampleID",
               values_to = "Abundance") %>%
  group_by(SampleID, trophicMode_clean) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

trophic_mode_levels <- c("Saprotroph", "Symbiotroph", "Pathotroph",
                         "Multi-trophic", "Unassigned")

trophic_counts$trophicMode_clean <- factor(trophic_counts$trophicMode_clean,
                                           levels = trophic_mode_levels)

trophic_rel <- trophic_counts %>%
  group_by(SampleID) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

trophic_counts <- trophic_counts %>%
  left_join(meta, by = "SampleID")

trophic_rel <- trophic_rel %>%
  left_join(meta, by = "SampleID")

trophic_counts_wide <- trophic_counts %>%
  dplyr::select(SampleID, trophicMode_clean, Abundance) %>%
  pivot_wider(names_from = trophicMode_clean, values_from = Abundance,
              values_fill = 0) %>%
  left_join(meta, by = "SampleID")

trophic_rel_wide <- trophic_rel %>%
  dplyr::select(SampleID, trophicMode_clean, RelAbundance) %>%
  pivot_wider(names_from = trophicMode_clean, values_from = RelAbundance,
              values_fill = 0) %>%
  left_join(meta, by = "SampleID")

write.csv(trophic_counts, "./results/tables/funguild_trophic_counts_long.csv",
          row.names = FALSE)

write.csv(trophic_rel, "./results/tables/funguild_trophic_relative_abundance_long.csv",
          row.names = FALSE)

write.csv(trophic_counts_wide, "./results/tables/funguild_trophic_counts_wide.csv",
          row.names = FALSE)

write.csv(trophic_rel_wide, "./results/tables/funguild_trophic_relative_abundance_wide.csv",
          row.names = FALSE)

## =========================================================
## 6. Beta-regression models for trophic modes
## =========================================================

options(contrasts = c("contr.sum", "contr.poly"))

guild_test <- trophic_rel %>%
  filter(trophicMode_clean != "Unassigned") %>%
  dplyr::select(SampleID, trophicMode_clean, RelAbundance,
                Intensity, Frequency, Composition, Pot) %>%
  group_by(SampleID, Intensity, Frequency, Composition, Pot) %>%
  complete(trophicMode_clean, fill = list(RelAbundance = 0)) %>%
  ungroup() %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         trophicMode_clean = factor(trophicMode_clean,
                                    levels = c("Saprotroph", "Symbiotroph", "Pathotroph", "Multi-trophic")))

fit_guild_model <- function(df) {
  n_df <- nrow(df)
  
  df <- df %>%
    mutate(RelAbundance_beta = (RelAbundance * (n_df - 1) + 0.5) / n_df)
  
  glmmTMB(RelAbundance_beta ~ Intensity * Frequency * Composition + (1 | Pot),
          data = df, family = beta_family(link = "logit"))
}

model_list <- guild_test %>%
  split(.$trophicMode_clean) %>%
  lapply(fit_guild_model)

anova_list <- lapply(model_list, function(m) car::Anova(m, type = 3))

guild_stats <- lapply(names(anova_list), function(g) {
  tab <- as.data.frame(anova_list[[g]])
  tab$Term <- rownames(tab)
  tab$Guild <- g
  tab
}) %>%
  bind_rows() %>%
  filter(Term != "(Intercept)") %>%
  rename(p = `Pr(>Chisq)`) %>%
  dplyr::select(Guild, Term, Chisq, Df, p)

# Save - Table S9
write.csv(guild_stats, "./results/tables/funguild_trophic_mode_beta_models.csv",
          row.names = FALSE)

