### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: PICRUSt2 pathway models and data preparation ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lme4)
library(car)
library(compositions)
library(stringr)

## =========================================================
## 1. Load input data
## =========================================================

NSTI <- read.table("./data/picrust2/picrust2_NSTI.tsv", header = TRUE, sep = "\t")

pathways <- read.table("./data/picrust2/picrust2_pathways.tsv", header = TRUE,
                       sep = "\t", row.names = 1, check.names = FALSE)

metadata <- read.csv("./data/picrust2/Metadata_B.csv", header = TRUE, check.names = FALSE)

pathway_info <- read.csv("./data/picrust2/metacyc-pwy_info.csv", header = TRUE, check.names = FALSE)

## =========================================================
## 2. NSTI summary
## =========================================================

nsti_summary <- data.frame(n_samples = length(NSTI$metadata_NSTI),
                           mean_nsti = mean(NSTI$metadata_NSTI, na.rm = TRUE),
                           median_nsti = median(NSTI$metadata_NSTI, na.rm = TRUE),
                           min_nsti = min(NSTI$metadata_NSTI, na.rm = TRUE),
                           max_nsti = max(NSTI$metadata_NSTI, na.rm = TRUE),
                           n_below_0.5 = sum(NSTI$metadata_NSTI < 0.5, na.rm = TRUE),
                           n_below_1 = sum(NSTI$metadata_NSTI < 1, na.rm = TRUE))

write.csv(nsti_summary, "./results/tables/picrust2_nsti_summary.csv", row.names = FALSE)

## =========================================================
## 3. Prepare pathway table
## =========================================================

# Transpose pathway table so rows = samples and columns = pathways
pathways_t <- pathways %>%
  t() %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("SampleID")

# Merge with metadata
pathways_meta <- pathways_t %>%
  left_join(metadata, by = "SampleID")

# Pathway columns correspond to original row names of pathway table
pathway_cols <- rownames(pathways)

# Check that all pathway columns are present
stopifnot(all(pathway_cols %in% colnames(pathways_meta)))

# Normalize to relative abundance
pathways_meta[pathway_cols] <- pathways_meta[pathway_cols] /
  rowSums(pathways_meta[pathway_cols], na.rm = TRUE)

# Filter rare pathways
keep <- colMeans(pathways_meta[, pathway_cols], na.rm = TRUE) > 0.0001
pathway_cols_filt <- pathway_cols[keep]

pathway_filter_summary <- data.frame(n_pathways_total = length(pathway_cols),
                                     n_pathways_retained = length(pathway_cols_filt),
                                     mean_abundance_threshold = 0.0001)

write.csv(pathway_filter_summary, "./results/tables/picrust2_pathway_filter_summary.csv",
          row.names = FALSE)

## =========================================================
## 4. Prepare metadata and CLR-transformed pathway table
## =========================================================

pathways_meta <- pathways_meta %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         Pot = factor(Pot))

options(contrasts = c("contr.sum", "contr.poly"))

pathway_mat <- as.matrix(pathways_meta[, pathway_cols_filt])

# Add small pseudocount to avoid log(0)
pseudocount <- 1e-6
pathway_clr <- clr(pathway_mat + pseudocount)

pathways_meta_clr <- bind_cols(pathways_meta %>%
                                 dplyr::select(SampleID, Pot, Intensity, Frequency, Composition),
                               as.data.frame(pathway_clr, check.names = FALSE))

write.csv(pathways_meta_clr, "./results/tables/picrust2_pathways_clr_table.csv", row.names = FALSE)

## =========================================================
## 5. Mixed-effects models for all retained pathways
## =========================================================

fit_pathway_model <- function(pathway_id, clr_table) {
  dat <- clr_table %>%
    dplyr::select(Intensity, Frequency, Composition, Pot, all_of(pathway_id)) %>%
    rename(PathwayAbundance = all_of(pathway_id)) %>%
    na.omit()
  
  mod <- lmer(PathwayAbundance ~ Intensity * Frequency * Composition + (1 | Pot),
              data = dat, REML = FALSE)
  
  aov_tab <- Anova(mod, type = 3)
  
  data.frame(PathwayID = pathway_id,
             Intensity_chisq = aov_tab["Intensity", "Chisq"],
             Intensity_p = aov_tab["Intensity", "Pr(>Chisq)"],
             Frequency_chisq = aov_tab["Frequency", "Chisq"],
             Frequency_p = aov_tab["Frequency", "Pr(>Chisq)"],
             Composition_chisq = aov_tab["Composition", "Chisq"],
             Composition_p = aov_tab["Composition", "Pr(>Chisq)"],
             IF_chisq = aov_tab["Intensity:Frequency", "Chisq"],
             IF_p = aov_tab["Intensity:Frequency", "Pr(>Chisq)"],
             IC_chisq = aov_tab["Intensity:Composition", "Chisq"],
             IC_p = aov_tab["Intensity:Composition", "Pr(>Chisq)"],
             FC_chisq = aov_tab["Frequency:Composition", "Chisq"],
             FC_p = aov_tab["Frequency:Composition", "Pr(>Chisq)"],
             IFC_chisq = aov_tab["Intensity:Frequency:Composition", "Chisq"],
             IFC_p = aov_tab["Intensity:Frequency:Composition", "Pr(>Chisq)"],
             singular = isSingular(mod))
}

pathway_stats_clr <- map_dfr(pathway_cols_filt, fit_pathway_model, clr_table = pathways_meta_clr)

## =========================================================
## 6. FDR correction and pathway names
## =========================================================

pathway_stats_clr <- pathway_stats_clr %>%
  mutate(Intensity_FDR = p.adjust(Intensity_p, method = "BH"),
         Frequency_FDR = p.adjust(Frequency_p, method = "BH"),
         Composition_FDR = p.adjust(Composition_p, method = "BH"),
         IF_FDR = p.adjust(IF_p, method = "BH"),
         IC_FDR = p.adjust(IC_p, method = "BH"),
         FC_FDR = p.adjust(FC_p, method = "BH"),
         IFC_FDR = p.adjust(IFC_p, method = "BH"))

pathway_info <- pathway_info %>%
  distinct(PathwayID, .keep_all = TRUE)

pathway_stats_clr <- pathway_stats_clr %>%
  left_join(pathway_info, by = "PathwayID") %>%
  relocate(PathwayName, .after = PathwayID)

# Save - Table S8
write.csv(pathway_stats_clr, "./results/tables/picrust2_pathway_model_results_clr.csv",
          row.names = FALSE)

## =========================================================
## 7. Rule-based pathway classification
## =========================================================

pathway_info_classified <- pathway_info %>%
  mutate(PathwayName_lc = str_to_lower(PathwayName),
    
    Class_L2 = case_when(
      str_detect(PathwayName_lc, "arginine|ornithine|citrulline|lysine|threonine|methionine|cysteine|serine|glycine|alanine|valine|leucine|isoleucine|tryptophan|tyrosine|phenylalanine|histidine|proline|aspartate|glutamate|glutamine|amino acid|ectoine|cyanophycin|norspermidine|polyamine") ~ "Amino acid metabolism",
      
      str_detect(PathwayName_lc, "cobalamin|thiamine|biotin|riboflavin|folate|folic acid|pyridoxal|pyridoxine|nad|fad|coenzyme a|coenzyme m|pantothenate|tetrahydrofolate|molybdopterin|heme|siroheme|cofactor|vitamin|quinone|flavin|preq0|queuosine|mycothiol|ergothioneine|cob\\(|menaquinol|demethylmenaquinol|naphthoate|naphthoquinol") ~ "Cofactor / vitamin metabolism",
      
      str_detect(PathwayName_lc, "purine|pyrimidine|adenosine|guanosine|inosine|uridine|thymidine|nucleotide|ump biosynthesis|atp|gtp|ctp|utp|dtdp") ~ "Nucleotide metabolism",
      
      str_detect(PathwayName_lc, "fatty acid|lipid|phospholipid|triacylglycerol|glycerolipid|glycerophospholipid|beta-oxidation|diacylglycerol|vaccenate|gondoate|ceramide|cholesterol|oleate|stearate|petroselinate|dodecenoate|farnesol|phosphatidate|phospholipase|sitosterol|phytol") ~ "Lipid metabolism",
      
      str_detect(PathwayName_lc, "starch|glycogen|sucrose|trehalose|maltose|galactose|arabinose|mannose|rhamnose|fucose|xylose|glucose|fructose|glucuronate|galacturonate|glucarate|galactarate|ketogluconate|ascorbate|stachyose|carbohydrate") ~ "Carbohydrate metabolism",
      
      str_detect(PathwayName_lc, "chitin|cellulose|glucan|xylan|pectin|polysaccharide|mannan|polymer|anhydromuropeptide") ~ "Polymer degradation",
      
      str_detect(PathwayName_lc, "aromatic|phenol|phenyl|benzoate|benzoyl|catechol|protocatechuate|cinnamate|gallate|mandelate|hydroxyphenylacetate|hydroxycinnamate|lignin|chorismate|salicylate|toluene|cresol") ~ "Aromatic compound metabolism",
      
      str_detect(PathwayName_lc, "nitrogen|nitrate|nitrite|ammonia|ammonium|urea|ureide|allantoin|denitrification|nitrification|creatinine|ethanolamine") ~ "Nitrogen metabolism",
      
      str_detect(PathwayName_lc, "sulfur|sulphur|sulfate|sulfite|sulfide") ~ "Sulfur metabolism",
      
      str_detect(PathwayName_lc, "phosphate|phosphonate|phosphorus|glyphosate") ~ "Phosphorus metabolism",
      
      str_detect(PathwayName_lc, "fermentation|methanogenesis|glycolysis|tca|tricarboxylic|electron transfer|respiration|reductive acetyl-coa|acetyl-coa|pyruvate|gluconeogenesis|glyoxylate|calvin|rubisco|methylcitrate|gaba shunt|entner-doudoroff|formaldehyde assimilation|formaldehyde oxidation|r[u]?mp|3-hydroxypropanoate cycle|ethylmalonyl-coa|ketogenesis|c1 compounds oxidation|oxidation to co2") ~ "Energy / central carbon metabolism",
      
      str_detect(PathwayName_lc, "peptidoglycan|lipopolysaccharide|cell wall|teichoic acid|teichuronic acid|o-antigen|muramoyl|acetylglucosamine|heptose|colanic acid|legionaminate|mannuronate|octulosonate|building blocks biosynthesis") ~ "Cell envelope / glycan metabolism",
      
      str_detect(PathwayName_lc, "isoprene|hopanoid|chlorophyll|bioluminescence|germacrene|flexixanthin|terpen|carotenoid|taxadiene") ~ "Secondary metabolite / isoprenoid metabolism",
      
      str_detect(PathwayName_lc, "dichloro|chlorobenzene|xenobiotic|octane oxidation|pollutant") ~ "Xenobiotic / pollutant degradation",
      
      str_detect(PathwayName_lc, "propanediol|butanediol|methylglyoxal|hexitol|mannitol|inositol|diol|dopamine|serotonin|adrenaline|noradrenaline|small-molecule|2-oxobutanoate") ~ "Organic acid / small-molecule metabolism",
      
      str_detect(PathwayName_lc, "trna|rrna|ribosome|ppgpp|dna|rna processing") ~ "Cellular processes",
      
      str_detect(PathwayName_lc, "detoxification|oxidative stress|repair|response") ~ "Stress / detoxification",
      
      TRUE ~ "Other")) %>%
  mutate(Class_L1 = case_when(Class_L2 %in% c("Amino acid metabolism",
                                              "Cofactor / vitamin metabolism",
                                              "Nucleotide metabolism",
                                              "Lipid metabolism",
                                              "Cell envelope / glycan metabolism",
                                              "Secondary metabolite / isoprenoid metabolism"
                                              ) ~ "Biosynthesis / cellular components",
      
                              Class_L2 %in% c("Carbohydrate metabolism",
                                              "Polymer degradation",
                                              "Aromatic compound metabolism",
                                              "Nitrogen metabolism",
                                              "Phosphorus metabolism",
                                              "Sulfur metabolism",
                                              "Xenobiotic / pollutant degradation",
                                              "Organic acid / small-molecule metabolism"
                                              ) ~ "Substrate utilization / nutrient cycling",
                              Class_L2 %in% c("Energy / central carbon metabolism") ~ "Energy metabolism",
                              Class_L2 %in% c("Cellular processes") ~ "Cellular processes",
                              Class_L2 %in% c("Stress / detoxification") ~ "Stress response",
                              TRUE ~ "Other")) %>%
  dplyr::select(-PathwayName_lc)

write.csv(pathway_info_classified, "./results/tables/picrust2_pathway_classification.csv",
          row.names = FALSE)

## =========================================================
## 8. Join classification to model results
## =========================================================

pathway_stats_clr_classified <- pathway_stats_clr %>%
  dplyr::select(-matches("^Category$|^ShortName$"), everything()) %>%
  left_join(pathway_info_classified %>%
              dplyr::select(PathwayID, Class_L1, Class_L2), by = "PathwayID")

write.csv(pathway_stats_clr_classified,
          "./results/tables/picrust2_pathway_model_results_clr_classified.csv",
          row.names = FALSE)

## =========================================================
## 9. Summary of significant pathways by term
## =========================================================

significance_summary <- data.frame(
  Term = c("Intensity", "Frequency", "Composition", "Intensity:Frequency",
           "Intensity:Composition", "Frequency:Composition", "Intensity:Frequency:Composition"),
  n_significant_FDR_0.05 = c(
    sum(pathway_stats_clr_classified$Intensity_FDR < 0.05, na.rm = TRUE),
    sum(pathway_stats_clr_classified$Frequency_FDR < 0.05, na.rm = TRUE),
    sum(pathway_stats_clr_classified$Composition_FDR < 0.05, na.rm = TRUE),
    sum(pathway_stats_clr_classified$IF_FDR < 0.05, na.rm = TRUE),
    sum(pathway_stats_clr_classified$IC_FDR < 0.05, na.rm = TRUE),
    sum(pathway_stats_clr_classified$FC_FDR < 0.05, na.rm = TRUE),
    sum(pathway_stats_clr_classified$IFC_FDR < 0.05, na.rm = TRUE)))

write.csv(significance_summary, "./results/tables/picrust2_pathway_significance_summary.csv",
          row.names = FALSE)

singular_summary <- data.frame(n_models = nrow(pathway_stats_clr_classified),
                               n_singular = sum(pathway_stats_clr_classified$singular, na.rm = TRUE))

write.csv(singular_summary, "./results/tables/picrust2_pathway_singular_summary.csv",
          row.names = FALSE)

## =========================================================
## 10. Effect-size tables for figure generation
## =========================================================

# ---------- Intensity ----------
pathway_effects_I <- pathways_meta_clr %>%
  dplyr::select(Intensity, all_of(pathway_cols_filt)) %>%
  pivot_longer(cols = all_of(pathway_cols_filt),
               names_to = "PathwayID", values_to = "CLR") %>%
  group_by(PathwayID, Intensity) %>%
  summarise(mean_CLR = mean(CLR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Intensity, values_from = mean_CLR) %>%
  mutate(diff_Dry_Wet = Dry - Wet,
         Direction = case_when(diff_Dry_Wet > 0 ~ "Higher in Dry",
                               diff_Dry_Wet < 0 ~ "Higher in Wet",
                               TRUE ~ "No change"))

sig_I_cat <- pathway_stats_clr_classified %>%
  filter(Intensity_FDR < 0.05) %>%
  dplyr::select(PathwayID, PathwayName, Class_L1, Class_L2, Intensity_FDR) %>%
  left_join(pathway_effects_I, by = "PathwayID") %>%
  filter(Direction != "No change")

write.csv(sig_I_cat, "./results/tables/picrust2_effects_intensity_significant.csv",
          row.names = FALSE)

plot_cat_eff_I <- sig_I_cat %>%
  group_by(Class_L1, Class_L2) %>%
  summarise(mean_diff = mean(diff_Dry_Wet, na.rm = TRUE),
            sd_diff = sd(diff_Dry_Wet, na.rm = TRUE),
            n = sum(!is.na(diff_Dry_Wet)),
            se_diff = sd_diff / sqrt(n),
            .groups = "drop") %>%
  mutate(Class_L1 = factor(Class_L1, levels = c("Biosynthesis / cellular components",
                                                "Substrate utilization / nutrient cycling",
                                                "Energy metabolism",
                                                "Cellular processes",
                                                "Stress response",
                                                "Other")),
         Direction = case_when(mean_diff > 0 ~ "Dry", mean_diff < 0 ~ "Wet",
                               TRUE ~ "Neutral")) %>%
  arrange(Class_L1, mean_diff) %>%
  mutate(Label = paste0(Class_L2, " (n=", n, ")"),
         Label = factor(Label, levels = unique(Label)))

write.csv(plot_cat_eff_I, "./results/tables/picrust2_effect_plot_intensity.csv",
          row.names = FALSE)

# ---------- Frequency ----------
pathway_effects_F <- pathways_meta_clr %>%
  dplyr::select(Frequency, all_of(pathway_cols_filt)) %>%
  pivot_longer(cols = all_of(pathway_cols_filt), names_to = "PathwayID", values_to = "CLR") %>%
  group_by(PathwayID, Frequency) %>%
  summarise(mean_CLR = mean(CLR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Frequency, values_from = mean_CLR) %>%
  mutate(diff_Low_High = Low - High)

sig_F_cat <- pathway_stats_clr_classified %>%
  filter(Frequency_FDR < 0.05) %>%
  dplyr::select(PathwayID, PathwayName, Class_L1, Class_L2, Frequency_FDR) %>%
  left_join(pathway_effects_F, by = "PathwayID")

write.csv(sig_F_cat, "./results/tables/picrust2_effects_frequency_significant.csv",
          row.names = FALSE)

plot_cat_eff_F <- sig_F_cat %>%
  group_by(Class_L1, Class_L2) %>%
  summarise(mean_diff = mean(diff_Low_High, na.rm = TRUE), sd_diff = sd(diff_Low_High, na.rm = TRUE),
            n = sum(!is.na(diff_Low_High)), se_diff = sd_diff / sqrt(n), .groups = "drop") %>%
  mutate(Class_L1 = factor(Class_L1,
                           levels = c("Biosynthesis / cellular components",
                                      "Substrate utilization / nutrient cycling",
                                      "Energy metabolism",
                                      "Cellular processes",
                                      "Stress response",
                                      "Other")),
         Direction = case_when(mean_diff > 0 ~ "Low", mean_diff < 0 ~ "High",
                               TRUE ~ "Neutral")) %>%
  arrange(Class_L1, mean_diff) %>%
  mutate(Label = paste0(Class_L2, " (n=", n, ")"),
         Label = factor(Label, levels = unique(Label)))

write.csv(plot_cat_eff_F, "./results/tables/picrust2_effect_plot_frequency.csv",
          row.names = FALSE)

# ---------- Intensity × Frequency ----------
sig_IF_ids <- pathway_stats_clr_classified %>%
  filter(IF_FDR < 0.05) %>%
  pull(PathwayID)

pathway_effects_IF <- pathways_meta_clr %>%
  dplyr::select(Intensity, Frequency, all_of(sig_IF_ids)) %>%
  pivot_longer(cols = all_of(sig_IF_ids),
               names_to = "PathwayID", values_to = "CLR") %>%
  group_by(PathwayID, Frequency, Intensity) %>%
  summarise(mean_CLR = mean(CLR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Intensity, values_from = mean_CLR) %>%
  mutate(diff_Dry_Wet = Dry - Wet)

sig_IF_cat <- pathway_stats_clr_classified %>%
  filter(IF_FDR < 0.05) %>%
  dplyr::select(PathwayID, PathwayName, Class_L1, Class_L2, IF_FDR) %>%
  left_join(pathway_effects_IF, by = "PathwayID")

write.csv(sig_IF_cat, "./results/tables/picrust2_effects_intensity_frequency_significant.csv",
          row.names = FALSE)

plot_cat_IF <- sig_IF_cat %>%
  group_by(Class_L1, Class_L2, Frequency) %>%
  summarise(mean_diff = mean(diff_Dry_Wet, na.rm = TRUE),
            sd_diff = sd(diff_Dry_Wet, na.rm = TRUE),
            n = sum(!is.na(diff_Dry_Wet)),
            se_diff = ifelse(n > 1, sd_diff / sqrt(n), 0),
            .groups = "drop") %>%
  mutate(Class_L1 = factor(Class_L1,
                           levels = c("Biosynthesis / cellular components",
                                      "Substrate utilization / nutrient cycling",
                                      "Energy metabolism",
                                      "Cellular processes",
                                      "Stress response",
                                      "Other")),
         Direction = case_when(mean_diff > 0 ~ "Dry", mean_diff < 0 ~ "Wet",
                               TRUE ~ "Neutral"),
         Frequency = factor(Frequency, levels = c("Low", "High"))) %>%
  arrange(Class_L1, mean_diff)

cat_order_IF <- sig_IF_cat %>%
  distinct(PathwayID, Class_L1, Class_L2) %>%
  count(Class_L1, Class_L2, name = "n_pathways") %>%
  left_join(plot_cat_IF %>%
              group_by(Class_L1, Class_L2) %>%
              summarise(overall_mean = mean(mean_diff, na.rm = TRUE), .groups = "drop"),
            by = c("Class_L1", "Class_L2")) %>%
  arrange(Class_L1, overall_mean) %>%
  mutate(Label = paste0(Class_L2, " (n=", n_pathways, ")"))

plot_cat_IF <- plot_cat_IF %>%
  left_join(cat_order_IF %>%
              dplyr::select(Class_L1, Class_L2, Label),
            by = c("Class_L1", "Class_L2")) %>%
  mutate(Label = factor(Label, levels = unique(cat_order_IF$Label)))

write.csv(plot_cat_IF, "./results/tables/picrust2_effect_plot_intensity_frequency.csv",
          row.names = FALSE)
