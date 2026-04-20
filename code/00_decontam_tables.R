### Project: Drought intensity and frequency (DIF) - Microbes ###
## Script: Decontamination of ASV tables ##
## Date: April 2026 | Author: Yang ##

# Change to work directory
setwd("~/Desktop/DIF-Microbes/Data availability")

# Load required packages
library(decontam)

## =========================================================
## 1. Helper function for decontamination
## =========================================================

run_decontam <- function(asv_file, tax_file, blank_ids, dataset_name,
                         output_asv_file, output_tax_file, output_contam_file,
                         threshold = 0.1) {
  
  ## Read ASV table
  count_tab <- read.table(asv_file, sep = "\t", header = TRUE, row.names = 1,
                          check.names = FALSE, comment.char = "", quote = "")
  
  ## Check that all blanks are present
  stopifnot(all(blank_ids %in% colnames(count_tab)))
  
  ## Identify negative controls
  neg_vec <- colnames(count_tab) %in% blank_ids
  cat(dataset_name, "- number of blanks detected:", sum(neg_vec), "\n")
  
  ## decontam expects samples as rows and ASVs as columns
  seqtab <- t(as.matrix(count_tab))
  
  ## Run decontam prevalence method
  contam_df <- isContaminant(seqtab, method = "prevalence", neg = neg_vec,
                             threshold = threshold)
  
  ## Extract contaminant ASVs
  contam_asvs <- rownames(contam_df)[contam_df$contaminant]
  
  cat(dataset_name, "- number of contaminants detected:", length(contam_asvs), "\n")
  
  ## Save contaminant summary
  contam_out <- data.frame(ASV = rownames(contam_df),
                           contaminant = contam_df$contaminant,
                           p = contam_df$p,
                           stringsAsFactors = FALSE)
  
  write.csv(contam_out, output_contam_file, row.names = FALSE)
  
  ## Remove contaminants
  count_tab_nocontam <- count_tab[!rownames(count_tab) %in% contam_asvs, , drop = FALSE]
  
  ## Remove negative controls from final analysis table
  count_tab_final <- count_tab_nocontam[, !(colnames(count_tab_nocontam) %in% blank_ids), drop = FALSE]
  
  ## Save filtered ASV table
  write.csv(count_tab_final, output_asv_file, quote = FALSE)
  
  ## Read taxonomy table
  tax_raw <- read.table(tax_file, sep = "\t", header = TRUE, quote = "",
                        comment.char = "", check.names = FALSE)
  
  ## Use Feature.ID as rownames
  rownames(tax_raw) <- tax_raw$Feature.ID
  
  ## Keep only ASVs retained in final ASV table and preserve order
  tax_filtered <- tax_raw[rownames(count_tab_final), , drop = FALSE]
  
  ## Save filtered taxonomy table
  write.table(tax_filtered, file = output_tax_file, sep = "\t",
              quote = FALSE, col.names = NA)
  
  ## Final consistency check
  stopifnot(identical(rownames(count_tab_final), rownames(tax_filtered)))
  
  cat(dataset_name, "- final ASVs retained:", nrow(count_tab_final), "\n")
  cat(dataset_name, "- final samples retained:", ncol(count_tab_final), "\n\n")
}

## =========================================================
## 2. Fungal decontamination
## =========================================================

run_decontam(asv_file           = "./data_raw/asv-table-F.tsv",
             tax_file           = "./data_raw/taxonomy_F.tsv",
             blank_ids          = c("NCE1F", "NPCR1F"),
             dataset_name       = "Fungi",
             output_asv_file    = "./data/asv_table_Fun.csv",
             output_tax_file    = "./data/taxonomy_Fun.csv",
             output_contam_file = "./results/tables/decontam_fungi_summary.csv",
             threshold          = 0.1)

## =========================================================
## 3. Bacterial decontamination
## =========================================================

run_decontam(asv_file           = "./data_raw/asv-table-B.tsv",
             tax_file           = "./data_raw/taxonomy_B.tsv",
             blank_ids          = c("NCE1B", "NPCR2B"),
             dataset_name       = "Bacteria",
             output_asv_file    = "./data/asv_table_Bac.csv",
             output_tax_file    = "./data/taxonomy_Bac.csv",
             output_contam_file = "./results/tables/decontam_bacteria_summary.csv",
             threshold          = 0.1)

