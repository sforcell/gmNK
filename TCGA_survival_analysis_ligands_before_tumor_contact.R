library(tidyverse)
library(xml2)
library(survival)
library(survminer)

# Ligand gene list
ligand_genes <- c(
  "THPO",        # Ligand of CD110 (MPL)
  "PDCD1LG2",    # Alternative ligand of CD279 (PD-1), in addition to CD274
  "IGSF8",       # Interactor of CD81 (EWI-2)
  "CD48",        # Ligand of CD244 (2B4)
  "PECAM1",      # Ligand of CD38
  "COL1A1",      # Collagen, ligand of CD29 (ITGB1)
  "FN1",         # Fibronectin, ligand of CD29 and CD44
  "VCAM1",       # Ligand of CD29
  "LAMA1",       # Laminin, ligand of CD29
  "CD58",        # Ligand of CD2
  "HAS2",        # Synthesizes hyaluronic acid, ligand of CD44
  "SPP1",        # Osteopontin, alternative ligand of CD44
  "PPIA",        # Cyclophilin A, ligand of CD147 (BSG)
  "S100A9",      # Another ligand/interactor of CD147 (BSG)
  "CD97",        # Ligand of CD55
  "ICAM1",       # Ligand of CD43 (SPN)
  "SELE",        # E-selectin, ligand of CD43 (SPN)
  "CD99"         # Homophilic ligand of CD99
)

# Function to extract clinical data
extract_clinical <- function(xml_file) {
  doc <- read_xml(xml_file)
  ns <- xml_ns(doc)
  
  patient_id <- xml_text(xml_find_first(doc, ".//shared:bcr_patient_barcode", ns))
  if (length(patient_id) == 0 || is.na(patient_id)) return(NULL)
  
  vital_status <- tolower(xml_text(xml_find_first(doc, ".//clin_shared:vital_status", ns)))
  days_to_death <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_death", ns)))
  days_to_last_followup <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_last_followup", ns)))
  days_to_recurrence <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_first_recurrence", ns)))
  
  # OS: time and status
  if (vital_status == "dead" && !is.na(days_to_death)) {
    os_time <- days_to_death
    os_status <- 1
  } else if (!is.na(days_to_last_followup)) {
    os_time <- days_to_last_followup
    os_status <- 0
  } else {
    return(NULL)
  }
  
  # EFS and DFS are not used here
  
  tibble(
    patient_id = substr(patient_id, 1, 12),
    OS_time = os_time,
    OS_status = os_status
  )
}

# Function to extract expression for all genes of interest
extract_expression_multiple <- function(tsv_file, mapping_df, genes) {
  file_id <- basename(dirname(tsv_file))
  barcode_row <- mapping_df %>% filter(file_id == !!file_id)
  if(nrow(barcode_row) == 0) return(NULL)
  
  data <- suppressMessages(read_tsv(tsv_file, comment = "#"))
  
  # Keep only genes of interest and average duplicates
  data_filtered <- data %>%
    filter(gene_name %in% genes) %>%
    group_by(gene_name) %>%
    summarise(tpm = mean(tpm_unstranded, na.rm = TRUE), .groups = "drop")

  if(nrow(data_filtered) == 0) return(NULL)

  tibble(
    patient_id = barcode_row$barcode[1],
    !!!set_names(data_filtered$tpm, data_filtered$gene_name)
  )

}

# Base paths
base_path <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
mapping_file <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/fileid_barcode_association.csv"
plot_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/Survivalplots_Ligands_OS_before_tumor_contact"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

mapping_df <- read_csv(mapping_file, col_types = cols())
tumor_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)

# Loop over tumors
for (tumor_path in tumor_dirs) {
  tumor_name <- basename(tumor_path)
  cat("\n### OS analysis for", tumor_name, "###\n")
  
  expression_files <- list.files(
    file.path(tumor_path, "Transcriptome_Profiling", "Gene_Expression_Quantification"),
    pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE
  )
  
  clinical_files <- list.files(
    file.path(tumor_path, "Clinical", "Clinical_Supplement"),
    pattern = "\\.xml$", recursive = TRUE, full.names = TRUE
  )
  
  if (length(expression_files) == 0 || length(clinical_files) == 0) {
    cat("Missing data for", tumor_name, "- skipping.\n")
    next
  }
  
  clinical_data <- bind_rows(lapply(clinical_files, extract_clinical)) %>%
    filter(!is.na(OS_time) & !is.na(OS_status))
  
  expression_list <- lapply(expression_files, extract_expression_multiple, mapping_df = mapping_df, genes = ligand_genes)
  expression_list <- expression_list[!sapply(expression_list, is.null)]
  expression_data <- bind_rows(expression_list)
  
  merged_data <- inner_join(clinical_data, expression_data, by = "patient_id")
  cat("  → Valid patients:", nrow(merged_data), "\n")
  
  if (nrow(merged_data) > 0) {
    tumor_plot_dir <- file.path(plot_dir, tumor_name)
    if (!dir.exists(tumor_plot_dir)) dir.create(tumor_plot_dir, recursive = TRUE)
    
    for (gene in ligand_genes) {
      if (!(gene %in% colnames(merged_data))) {
        cat("  → Gene", gene, "not found in expression data, skipping.\n")
        next
      }
      
      gene_data <- merged_data %>% filter(!is.na(.data[[gene]]))
      if (nrow(gene_data) == 0) {
        cat("  → No data for gene", gene, "\n")
        next
      }
      
      q20 <- quantile(gene_data[[gene]], probs = 0.20, na.rm = TRUE)
      q80 <- quantile(gene_data[[gene]], probs = 0.80, na.rm = TRUE)
      
      gene_data <- gene_data %>%
        mutate(Gene_group = case_when(
          .data[[gene]] <= q20 ~ "Low",
          .data[[gene]] >= q80 ~ "High",
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(Gene_group))
      
      if (nrow(gene_data) < 10) {
        cat("  → Too few patients (", nrow(gene_data), ") for gene", gene, "- skipping.\n")
        next
      }
      
      os_surv <- Surv(time = gene_data$OS_time, event = gene_data$OS_status)
      os_fit <- survfit(os_surv ~ Gene_group, data = gene_data)
      
      os_plot <- ggsurvplot(os_fit, data = gene_data, pval = TRUE, risk.table = TRUE,
                            legend.labs = c("Low", "High"),
                            palette = c("#1f77b4", "#d62728"),
                            title = paste("Overall Survival in", tumor_name, "-", gene),
                            legend.title = paste(gene, "expression group"),
                            ggtheme = theme_minimal() + theme(
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 12),
                              axis.title = element_text(size = 14),
                              axis.text = element_text(size = 12)
                            ))
      
      ggsave(filename = file.path(tumor_plot_dir, paste0("OS_", gene, "_", tumor_name, ".pdf")),
             plot = os_plot$plot, width = 5, height = 5, dpi = 300)
    }
  } else {
    cat("  → Not enough data for OS analysis.\n")
  }
}

