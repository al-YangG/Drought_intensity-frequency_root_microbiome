### Project: Drought intensity and frequency (DIF) - Microbes ###
## Indicator ASV figures ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

library(phyloseq)
library(dplyr)
library(stringr)
library(forcats)
library(tidytext)
library(ggplot2)
library(ggpubr)
library(grid)

# =========================================================
# 1. Load data
# =========================================================

ps_Bac <- readRDS("./data/phyloseq/ps_Bac_filtered.rds")
ps_Fun <- readRDS("./data/phyloseq/ps_Fun.rds")

bac_sig_intensity <- read.csv("./results/tables/indicator_asvs_bacteria_intensity.csv",
                              header = TRUE, check.names = FALSE)

bac_sig_frequency <- read.csv("./results/tables/indicator_asvs_bacteria_frequency.csv",
                              header = TRUE, check.names = FALSE)

fun_sig_intensity <- read.csv("./results/tables/indicator_asvs_fungi_intensity.csv",
                              header = TRUE, check.names = FALSE)

fun_sig_frequency <- read.csv("./results/tables/indicator_asvs_fungi_frequency.csv",
                              header = TRUE, check.names = FALSE)

# =========================================================
# 2. Phylum-level summary for bubble plot
# =========================================================

make_phylum_summary <- function(df, kingdom, factor_name) {
  if (nrow(df) == 0) {
    return(data.frame(Kingdom = character(), Factor = character(), Indicator = character(),
                      Phylum = character(), n_ASVs = integer()))
  }
  
  df %>%
    mutate(Kingdom = kingdom, Factor = factor_name,
           Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum)) %>%
    count(Kingdom, Factor, Indicator, Phylum, name = "n_ASVs")
}

plot_bac_intensity <- make_phylum_summary(bac_sig_intensity, "Bacteria", "Intensity")
plot_bac_frequency <- make_phylum_summary(bac_sig_frequency, "Bacteria", "Frequency")

plot_fun_intensity <- make_phylum_summary(fun_sig_intensity, "Fungi", "Intensity")
plot_fun_frequency <- make_phylum_summary(fun_sig_frequency, "Fungi", "Frequency")

plot_indval_combined <- bind_rows(plot_bac_intensity, plot_bac_frequency,
                                  plot_fun_intensity, plot_fun_frequency) %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Bacteria", "Fungi")),
         Factor = factor(Factor, levels = c("Intensity", "Frequency")),
         Indicator = factor(Indicator,
                            levels = c("Dry", "Wet", "Low", "High")))

# Keep top 10 phyla within each kingdom based on total number of indicator ASVs
top_phyla_by_kingdom <- plot_indval_combined %>%
  group_by(Kingdom, Phylum) %>%
  summarise(indicator_ASVs = sum(n_ASVs), .groups = "drop_last") %>%
  arrange(desc(indicator_ASVs), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

plot_indval_top <- plot_indval_combined %>%
  semi_join(top_phyla_by_kingdom, by = c("Kingdom", "Phylum"))

# =========================================================
# 3. Total ASVs per phylum in original datasets
# =========================================================

tax_bac <- as.data.frame(tax_table(ps_Bac))
tax_bac$ASV <- rownames(tax_bac)

tax_fun <- as.data.frame(tax_table(ps_Fun))
tax_fun$ASV <- rownames(tax_fun)

phylum_total_bac <- tax_bac %>%
  mutate(Kingdom = "Bacteria",
         Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum)) %>%
  count(Kingdom, Phylum, name = "total_ASVs")

phylum_total_fun <- tax_fun %>%
  mutate(Kingdom = "Fungi",
         Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum)) %>%
  count(Kingdom, Phylum, name = "total_ASVs")

phylum_total_all <- bind_rows(phylum_total_bac, phylum_total_fun)

phylum_indicator_total <- plot_indval_combined %>%
  group_by(Kingdom, Phylum) %>%
  summarise(indicator_ASVs = sum(n_ASVs), .groups = "drop")

phylum_labels <- phylum_indicator_total %>%
  left_join(phylum_total_all, by = c("Kingdom", "Phylum")) %>%
  mutate(Phylum_label = paste0(Phylum, "\n(", indicator_ASVs, "/", total_ASVs, ")"))

plot_indval_top <- plot_indval_top %>%
  left_join(phylum_labels %>%
              dplyr::select(Kingdom, Phylum, indicator_ASVs, total_ASVs, Phylum_label),
            by = c("Kingdom", "Phylum"))

label_order <- phylum_labels %>%
  semi_join(top_phyla_by_kingdom, by = c("Kingdom", "Phylum")) %>%
  group_by(Kingdom) %>%
  arrange(indicator_ASVs, .by_group = TRUE) %>%
  ungroup() %>%
  pull(Phylum_label)

plot_indval_top <- plot_indval_top %>%
  mutate(Phylum_label = factor(Phylum_label, levels = label_order))

# =========================================================
# 4. Bubble plot
# =========================================================

indicator_colors <- c("Dry" = "#E7782F", "Wet" = "#1994E5",
                      "Low" = "#E5D019", "High" = "#B38AF1")

p_indval_bac_fun <- ggplot(plot_indval_top,
                           aes(x = Indicator, y = Phylum_label, size = n_ASVs, fill = Indicator)) +
  geom_point(shape = 21, color = "black", alpha = 0.6) +
  geom_text(aes(label = n_ASVs), size = 3, fontface = "bold") +
  facet_grid(Kingdom ~ Factor, scales = "free", space = "free") +
  scale_size_continuous(range = c(3, 13)) +
  scale_fill_manual(values = indicator_colors) +
  labs(x = NULL, y = "Phylum (indicator / total ASVs)") +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(face = "bold", size = 9, color = "black"),
        axis.title.y = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.spacing = unit(1.5, "lines"))
p_indval_bac_fun

# =========================================================
# 5. Helper for taxon labels
# =========================================================

make_tax_label_clean <- function(df) {
  bad_terms <- c("", NA, "uncultured", "Uncultured", "uncultured_bacterium", "uncultured_fungus",
                 "Incertae_sedis", "incertae_sedis", "gen_Incertae_sedis",
                 "family_Incertae_sedis", "order_Incertae_sedis", "class_Incertae_sedis")
  
  df %>%
    mutate(Genus_clean  = ifelse(Genus  %in% bad_terms | str_detect(Genus,  "uncultured|Incertae|incertae"), NA, Genus),
           Family_clean = ifelse(Family %in% bad_terms | str_detect(Family, "uncultured|Incertae|incertae"), NA, Family),
           Order_clean  = ifelse(Order  %in% bad_terms | str_detect(Order,  "uncultured|Incertae|incertae"), NA, Order),
           Class_clean  = ifelse(Class  %in% bad_terms | str_detect(Class,  "uncultured|Incertae|incertae"), NA, Class),
           Phylum_clean = ifelse(Phylum %in% bad_terms | str_detect(Phylum, "uncultured|Incertae|incertae"), NA, Phylum),
           Taxon_base = case_when(!is.na(Genus_clean)  ~ paste0("g_", Genus_clean),
                                  !is.na(Family_clean) ~ paste0("f_", Family_clean),
                                  !is.na(Order_clean)  ~ paste0("o_", Order_clean),
                                  !is.na(Class_clean)  ~ paste0("c_", Class_clean),
                                  !is.na(Phylum_clean) ~ paste0("p_", Phylum_clean),
                                  TRUE ~ "Unclassified")) %>%
    dplyr::select(-ends_with("_clean"))
}

# =========================================================
# 6. Top 10 indicator ASVs for Intensity
# =========================================================

top_bac_intensity <- bac_sig_intensity %>%
  make_tax_label_clean() %>%
  mutate(Kingdom = "Bacteria") %>%
  group_by(Kingdom, Indicator) %>%
  slice_max(order_by = stat, n = 10, with_ties = FALSE) %>%
  ungroup()

top_fun_intensity <- fun_sig_intensity %>%
  make_tax_label_clean() %>%
  mutate(Kingdom = "Fungi") %>%
  group_by(Kingdom, Indicator) %>%
  slice_max(order_by = stat, n = 10, with_ties = FALSE) %>%
  ungroup()

top_indval_intensity <- bind_rows(top_bac_intensity, top_fun_intensity) %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Bacteria", "Fungi")),
         Indicator = factor(Indicator, levels = c("Dry", "Wet")),
         Taxon_label = paste0(Taxon_base, " [", Indicator, "]"),
         Taxon_id = paste0(Taxon_label, "___", row_number())) %>%
  mutate(Taxon_id = reorder_within(Taxon_id, stat, Kingdom))

indicator_colors_intensity <- c("Dry" = "#E7782F", "Wet" = "#1994E5")

p_top_indval_intensity <- ggplot(top_indval_intensity,
                                 aes(x = stat, y = Taxon_id, color = Indicator)) +
  geom_point(size = 3.5) +
  facet_grid(Kingdom ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered(labels = function(x) gsub("___.*$", "", x)) +
  scale_color_manual(values = indicator_colors_intensity) +
  labs(x = "Indicator value (IndVal statistic)", y = "Taxon [indicator group]",
       color = "Indicator") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9, color = "black", face = "bold"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 11),
        legend.title = element_blank(),
        legend.position = c(0.35, 0.95),
        legend.text = element_text(size = 9, face = "bold"),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.spacing = unit(1.2, "lines"))
p_top_indval_intensity

# =========================================================
# 7. Top 10 indicator ASVs for Frequency
# =========================================================

top_bac_frequency <- bac_sig_frequency %>%
  make_tax_label_clean() %>%
  mutate(Kingdom = "Bacteria") %>%
  group_by(Kingdom, Indicator) %>%
  slice_max(order_by = stat, n = 10, with_ties = FALSE) %>%
  ungroup()

top_fun_frequency <- fun_sig_frequency %>%
  make_tax_label_clean() %>%
  mutate(Kingdom = "Fungi") %>%
  group_by(Kingdom, Indicator) %>%
  slice_max(order_by = stat, n = 10, with_ties = FALSE) %>%
  ungroup()

top_indval_frequency <- bind_rows(top_bac_frequency, top_fun_frequency) %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Bacteria", "Fungi")),
         Indicator = factor(Indicator, levels = c("Low", "High")),
         Taxon_label = paste0(Taxon_base, " [", Indicator, "]"),
         Taxon_id = paste0(Taxon_label, "___", row_number())) %>%
  mutate(Taxon_id = reorder_within(Taxon_id, stat, Kingdom))

indicator_colors_frequency <- c("Low" = "#E5D019", "High" = "#B38AF1")

p_top_indval_frequency <- ggplot(top_indval_frequency,
                                 aes(x = stat, y = Taxon_id, color = Indicator)) +
  geom_point(size = 3.5) +
  facet_grid(Kingdom ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered(labels = function(x) gsub("___.*$", "", x)) +
  scale_color_manual(values = indicator_colors_frequency) +
  labs(x = "Indicator value (IndVal statistic)",
       y = "Taxon [indicator group]", color = "Indicator") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9, color = "black", face = "bold"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 11),
        legend.title = element_blank(),
        legend.position = c(0.35, 0.95),
        legend.text = element_text(size = 9, face = "bold"),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.spacing = unit(1.2, "lines"))
p_top_indval_frequency

# =========================================================
# 8. Combine and save figure
# =========================================================

p_ind <- ggarrange(p_indval_bac_fun, p_top_indval_intensity,
                   p_top_indval_frequency, nrow = 1,
                   labels = c("(a)", "(b)", "(c)"))
p_ind

# Save - Figure 4
ggsave(filename = "./results/figures/indicator_asvs.pdf", plot = p_ind,
       width = 18, height = 10, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/indicator_asvs.png", plot = p_ind,
       width = 18, height = 10, dpi = 600)

