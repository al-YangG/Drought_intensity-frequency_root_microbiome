### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: SEM follow-up plots ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

library(dplyr)
library(ggplot2)
library(ggpubr)

# =========================================================
# 1. Load data
# =========================================================

data <- read.csv("./data/SEM data-plant+microbes.csv", header = TRUE, check.names = FALSE)

data_sem <- data %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Intensity_num = ifelse(Intensity == "Dry", -0.5, 0.5),
         Frequency_num = ifelse(Frequency == "Low", -0.5, 0.5),
         SRL_z = scale(SRL)[, 1],
         RD_z = scale(RD)[, 1],
         TBs_z = scale(TBs)[, 1],
         Bac_PCoA1_z = scale(Bac_PCoA1)[, 1]) %>%
  dplyr::select(Intensity, Frequency, Intensity_num, Frequency_num,
                SRL, RD, TBs, Bac_PCoA1, SRL_z, RD_z, TBs_z, Bac_PCoA1_z) %>%
  na.omit()

# =========================================================
# 2. Helper functions
# =========================================================

treatment_cols <- c("Dry" = "#E7782F", "Wet" = "#1994E5")

get_lm_stats <- function(df, yvar, xvar = "Bac_PCoA1") {
  form <- as.formula(paste(yvar, "~", xvar))
  mod <- lm(form, data = df)
  sm <- summary(mod)
  
  slope <- coef(mod)[2]
  pval <- sm$coefficients[2, 4]
  
  p_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", pval)))
  
  ltype <- ifelse(pval < 0.05, "solid", "dashed")
  
  data.frame(slope = slope, pval = pval, p_label = p_label, ltype = ltype)
}

get_lm_label <- function(df, yvar, xvar = "Bac_PCoA1_c") {
  form <- as.formula(paste(yvar, "~", xvar))
  mod <- lm(form, data = df)
  sm <- summary(mod)
  
  slope <- coef(mod)[2]
  pval <- sm$coefficients[2, 4]
  
  p_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", pval)))
  
  paste0("\u03b2 = ", sprintf("%.2f", slope), ", ", p_label)
}

# =========================================================
# 3. Raw-scale plots by intensity
# =========================================================

x_min <- min(data_sem$Bac_PCoA1, na.rm = TRUE)
x_max <- max(data_sem$Bac_PCoA1, na.rm = TRUE)
x_limits <- c(x_min, x_max)

# ---------- SRL ----------
srl_stats <- data_sem %>%
  group_by(Intensity) %>%
  do(get_lm_stats(., "SRL")) %>%
  ungroup() %>%
  mutate(label = paste0(Intensity, ": \u03b2 = ", sprintf("%.2f", slope), ", ", p_label))

y_min_srl <- min(data_sem$SRL, na.rm = TRUE)
y_max_srl <- max(data_sem$SRL, na.rm = TRUE)

srl_stats <- srl_stats %>%
  mutate(x = x_min + 0.05 * (x_max - x_min),
         y = c(y_max_srl - 0.05 * (y_max_srl - y_min_srl),
               y_max_srl - 0.12 * (y_max_srl - y_min_srl)))

p_srl_plot <- ggplot(data_sem, aes(x = Bac_PCoA1, y = SRL, color = Intensity)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(aes(linetype = Intensity), method = "lm", se = FALSE,
              linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = treatment_cols) +
  scale_linetype_manual(values = setNames(srl_stats$ltype, srl_stats$Intensity)) +
  coord_cartesian(xlim = x_limits) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  labs(x = "Bacterial composition (PCoA1)", y = "Specific root length (SRL)") +
  geom_text(data = srl_stats, aes(x = x, y = y, label = label, color = Intensity),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5,
            fontface = "bold", show.legend = FALSE)
p_srl_plot

# ---------- RD ----------
rd_stats <- data_sem %>%
  group_by(Intensity) %>%
  do(get_lm_stats(., "RD")) %>%
  ungroup() %>%
  mutate(label = paste0(Intensity, ": \u03b2 = ", sprintf("%.2f", slope), ", ", p_label))

y_min_rd <- min(data_sem$RD, na.rm = TRUE)
y_max_rd <- max(data_sem$RD, na.rm = TRUE)

rd_stats <- rd_stats %>%
  mutate(x = x_min + 0.05 * (x_max - x_min),
         y = c(y_max_rd - 0.05 * (y_max_rd - y_min_rd),
               y_max_rd - 0.12 * (y_max_rd - y_min_rd)))

p_rd_plot <- ggplot(data_sem, aes(x = Bac_PCoA1, y = RD, color = Intensity)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(aes(linetype = Intensity), method = "lm", se = FALSE,
              linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = treatment_cols) +
  scale_linetype_manual(values = setNames(rd_stats$ltype, rd_stats$Intensity)) +
  coord_cartesian(xlim = x_limits) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  labs(x = "Bacterial composition (PCoA1)", y = "Root diameter (RD)") +
  geom_text(data = rd_stats, aes(x = x, y = y, label = label, color = Intensity),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5,
            fontface = "bold", show.legend = FALSE)
p_rd_plot

# =========================================================
# 4. Centered plots consistent with SEM logic
# =========================================================

data_centered <- data_sem %>%
  group_by(Intensity) %>%
  mutate(Bac_PCoA1_c = Bac_PCoA1 - mean(Bac_PCoA1, na.rm = TRUE),
         SRL_c = SRL - mean(SRL, na.rm = TRUE),
         RD_c = RD - mean(RD, na.rm = TRUE)) %>%
  ungroup()

# ---------- Centered SRL ----------
srl_label <- get_lm_label(data_centered, "SRL_c", "Bac_PCoA1_c")

p_srl_cent <- ggplot(data_centered, aes(x = Bac_PCoA1_c, y = SRL_c, color = Intensity)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4, linetype = "dotted") +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4, linetype = "dotted") +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  scale_color_manual(values = treatment_cols) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10)) +
  labs(x = "Bacterial composition (centered)", y = "SRL (centered)") +
  annotate("text", x = Inf, y = Inf, label = srl_label, hjust = 1.1, vjust = 1.5,
           size = 3.5, fontface = "bold")
p_srl_cent

# ---------- Centered RD ----------
rd_label <- get_lm_label(data_centered, "RD_c", "Bac_PCoA1_c")

p_rd_cent <- ggplot(data_centered, aes(x = Bac_PCoA1_c, y = RD_c, color = Intensity)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4, linetype = "dotted") +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4, linetype = "dotted") +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  scale_color_manual(values = treatment_cols) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10)) +
  labs(x = "Bacterial composition (centered)", y = "RD (centered)") +
  annotate("text", x = Inf, y = Inf, label = rd_label, hjust = 1.1, vjust = 1.5,
           size = 3.5, fontface = "bold")
p_rd_cent

# =========================================================
# 5. Combine and save figure
# =========================================================

final_plot <- ggarrange(p_srl_plot, p_rd_plot, p_srl_cent, p_rd_cent,
                        nrow = 2, ncol = 2, common.legend = TRUE,
                        legend = "bottom", font.label = list(size = 12),
                        labels = c("(a)", "(b)", "(c)", "(d)"))

# Save - Figure S8
ggsave(filename = "./results/figures/srl_rd_vs_bac_pcoa1.pdf", plot = final_plot,
       width = 10, height = 8, units = "in", device = cairo_pdf)

ggsave(filename = "./results/figures/png/srl_rd_vs_bac_pcoa1.png", plot = final_plot,
       width = 10, height = 8, units = "in", dpi = 300, bg = "white")

