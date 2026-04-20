### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Piecewise SEM ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

library(dplyr)
library(piecewiseSEM)

# =========================================================
# 1. Load data
# =========================================================

data <- read.csv("./data/SEM data-plant+microbes.csv", header = TRUE, check.names = FALSE)

# =========================================================
# 2. Prepare SEM data
# =========================================================

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
# 3. Fit component models
# =========================================================

m_bac <- lm(Bac_PCoA1_z ~ Intensity_num * Frequency_num, data = data_sem)

m_srl <- lm(SRL_z ~ Intensity_num + Bac_PCoA1_z, data = data_sem)

m_rd <- lm(RD_z ~ Intensity_num + Bac_PCoA1_z, data = data_sem)

m_tbs <- lm(TBs_z ~ Frequency_num + Bac_PCoA1_z, data = data_sem)

# =========================================================
# 4. Fit piecewise SEM
# =========================================================

sem1 <- psem(m_bac, m_srl, m_rd, m_tbs, SRL_z %~~% RD_z)
summary(sem1)

# =========================================================
# 5. SEM summaries
# =========================================================

# Standardized path coefficients
sem_coefs <- coefs(sem1, standardize = "scale")
sem_coefs

# R2 for each component model
sem_r2 <- rsquared(sem1)
sem_r2

# Directed separation / global fit table
sem_summary <- summary(sem1)
sem_dsep <- sem_summary$dTable
sem_dsep
