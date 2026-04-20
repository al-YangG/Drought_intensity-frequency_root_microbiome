### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: SEM correlation and screening analyses ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(Hmisc)
library(reshape2)
library(ggplot2)

# =========================================================
# 1. Load data
# =========================================================

data <- read.csv("./data/SEM data-plant+microbes.csv", header = TRUE, check.names = FALSE)

data <- data %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))

# =========================================================
# 2. Define variable groups
# =========================================================

fitness_vars <- c("SPs", "RPs", "TBs")
biomass_vars <- c("RB", "TB", "R_S")
root_trait_vars <- c("RD", "RDMC", "SRL", "SRA", "RTD")
root_chem_vars <- c("RCC", "RNC", "Root_CN")

micro_div_vars <- c("Bac_Richness", "Bac_Shannon", "Fun_Richness", "Fun_Shannon")
micro_comp_vars <- c("Bac_PCoA1", "Bac_PCoA2", "Fun_PCoA1", "Fun_PCoA2")

ordered_vars <- c(fitness_vars, biomass_vars, root_trait_vars, root_chem_vars,
                  micro_div_vars, micro_comp_vars)

plant_groups <- list("Plant fitness"   = fitness_vars,
                     "Plant biomass"   = biomass_vars,
                     "Root morphology" = root_trait_vars,
                     "Root chemistry"  = root_chem_vars)

micro_groups <- list("Microbial diversity"   = micro_div_vars,
                     "Microbial composition" = micro_comp_vars)

plant_vars <- unlist(plant_groups, use.names = FALSE)
micro_vars <- unlist(micro_groups, use.names = FALSE)

# =========================================================
# 3. Full correlation matrix table
# =========================================================

cor_data_full <- data %>%
  dplyr::select(all_of(ordered_vars)) %>%
  mutate(across(everything(), as.numeric))

cor_res_full <- rcorr(as.matrix(cor_data_full), type = "pearson")
cor_mat_full <- cor_res_full$r
p_mat_full <- cor_res_full$P

cor_mat_round <- round(cor_mat_full, 2)

sig_mat <- matrix("", nrow = nrow(p_mat_full), ncol = ncol(p_mat_full))
sig_mat[p_mat_full < 0.05] <- "*"
sig_mat[p_mat_full < 0.01] <- "**"
sig_mat[p_mat_full < 0.001] <- "***"

cor_display <- matrix(paste0(cor_mat_round, sig_mat),
                      nrow = nrow(cor_mat_full),
                      dimnames = dimnames(cor_mat_full))

cor_table_wide <- as.data.frame(cor_display)

write.csv(cor_table_wide, "./results/tables/sem_correlation_matrix_table.csv", row.names = TRUE)

# =========================================================
# 4. Plant–microbe correlation heatmap
# =========================================================

all_vars <- intersect(c(plant_vars, micro_vars), colnames(data))

cor_data_pm <- data %>%
  dplyr::select(all_of(all_vars)) %>%
  mutate(across(everything(), as.numeric))

cor_res_pm <- rcorr(as.matrix(cor_data_pm), type = "pearson")
cor_mat_pm <- cor_res_pm$r
p_mat_pm <- cor_res_pm$P

plant_vars2 <- intersect(plant_vars, colnames(cor_mat_pm))
micro_vars2 <- intersect(micro_vars, colnames(cor_mat_pm))

cor_pm <- cor_mat_pm[plant_vars2, micro_vars2, drop = FALSE]
p_pm <- p_mat_pm[plant_vars2, micro_vars2, drop = FALSE]

cor_df <- melt(cor_pm)
p_df <- melt(p_pm)

colnames(cor_df) <- c("Plant_var", "Micro_var", "Correlation")
colnames(p_df)   <- c("Plant_var", "Micro_var", "P_value")

plot_df <- left_join(cor_df, p_df, by = c("Plant_var", "Micro_var"))

plot_df <- plot_df %>%
  mutate(Plant_group = case_when(Plant_var %in% plant_groups[["Plant fitness"]]   ~ "Plant fitness",
                                 Plant_var %in% plant_groups[["Plant biomass"]]   ~ "Plant biomass",
                                 Plant_var %in% plant_groups[["Root morphology"]] ~ "Root morphology",
                                 Plant_var %in% plant_groups[["Root chemistry"]]  ~ "Root chemistry"),
         Micro_group = case_when(Micro_var %in% micro_groups[["Microbial diversity"]]   ~ "Microbial diversity",
                                 Micro_var %in% micro_groups[["Microbial composition"]] ~ "Microbial composition"),
         sig = case_when(P_value < 0.001 ~ "***", P_value < 0.01  ~ "**", P_value < 0.05  ~ "*",
                         TRUE ~ ""))

plot_df$Plant_group <- factor(plot_df$Plant_group, levels = names(plant_groups))
plot_df$Micro_group <- factor(plot_df$Micro_group, levels = names(micro_groups))
plot_df$Plant_var <- factor(plot_df$Plant_var, levels = rev(plant_vars2))
plot_df$Micro_var <- factor(plot_df$Micro_var, levels = micro_vars2)

p_corr_heatmap <- ggplot(plot_df, aes(x = Micro_var, y = Plant_var, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(P_value < 0.05, sprintf("%.2f", Correlation), "")),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  facet_grid(Plant_group ~ Micro_group, scales = "free", space = "free") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
        strip.text = element_text(face = "bold", size = 11),
        axis.title = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"))
p_corr_heatmap

# Save - Figure S7
ggsave(filename = "./results/figures/microbial_plant_correlation_heatmap.pdf", plot = p_corr_heatmap,
       width = 6, height = 8, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/microbial_plant_correlation_heatmap.png", plot = p_corr_heatmap,
       width = 6, height = 8, units = "in", dpi = 300)

# =========================================================
# 5. Screening plant traits against microbial predictors
# =========================================================

plant_response_vars <- c("SPs", "RPs", "TBs", "RB", "TB", "R_S",
                         "RD", "RDMC", "SRL", "SRA", "RTD", "RCC", "RNC", "Root_CN")

micro_predictors <- c("Bac_Richness", "Bac_Shannon", "Fun_Richness", "Fun_Shannon",
                      "Bac_PCoA1", "Bac_PCoA2", "Fun_PCoA1", "Fun_PCoA2")

options(contrasts = c("contr.sum", "contr.poly"))

fit_trait_model <- function(response, predictor, data) {
  form <- as.formula(paste(response, "~ Intensity * Frequency * Composition +", predictor))
  
  mod <- lm(form, data = data)
  
  anova_tab <- car::Anova(mod, type = 3)
  anova_df <- as.data.frame(anova_tab)
  anova_df$term <- rownames(anova_df)
  
  pred_row <- anova_df %>%
    filter(term == predictor) %>%
    transmute(Df = Df, F_value = `F value`, P_value = `Pr(>F)`)
  
  coef_tab <- broom::tidy(mod) %>%
    filter(term == predictor) %>%
    transmute(Estimate = estimate, Std_error = std.error, t_value = statistic)
  
  fit_tab <- tibble(Adj_R2 = summary(mod)$adj.r.squared, AIC = AIC(mod), N = nobs(mod))
  
  bind_cols(pred_row, coef_tab, fit_tab)
}

all_results <- expand.grid(Response = plant_response_vars,
                           Predictor = micro_predictors,
                           stringsAsFactors = FALSE) %>%
  mutate(result = map2(Response, Predictor, ~ fit_trait_model(.x, .y, data))) %>%
  unnest(result) %>%
  mutate(P_adj_BH = p.adjust(P_value, method = "BH")) %>%
  arrange(P_value)

sig_results <- all_results %>%
  filter(P_value < 0.05) %>%
  arrange(P_value)


write.csv(sig_results, "./results/tables/plant_microbe_screening_significant.csv", row.names = FALSE)

