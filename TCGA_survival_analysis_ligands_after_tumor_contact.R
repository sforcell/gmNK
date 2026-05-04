library(tidyverse)
library(xml2)
library(survival)
library(survminer)

# Gene signature list
ligand_genes <- c(
  "ADA",         # Ligand of CD26 (DPP4), adenosine deaminase
  "CXCL10",      # Substrate and functional interactor of CD26
  "CD38",        # Interactor of PECAM1 (CD31)
  "CD177",       # Ligand of CD31 (PECAM1)
  "ENG",         # Endoglin, interacts with CD31 (PECAM1)
  "VEGFA",       # Induces signaling via PECAM1/CD31
  "SPN",         # Sialophorin, known interactor of CD69
  "CD69L",       # Putative ligand of CD69 (poorly characterized)
  "PTPRC",       # CD45, parental isoform of CD45RB
  "IL6ST",       # Functional interactor with CD45RB via signaling
  "CD53",        # Tetraspanin, modulates interaction with CD2, CD4, CD5
  "CD2",         # Interactor of CD53 in tetraspanin domains
  "CD300A",      # Ligand/interactor of CD148 (PTPRJ)
  "PTPRJ",       # CD148, interacts with CD300 family receptors
  "TNFRSF1A",    # Signaling interactor modulated by CD148
  "HMGB1",       # Immunomodulatory interactor, involved in CD69-related pathways
  "LGALS9",      # Galectin-9, co-interactor in pathways involving CD45RB
  "SELL",        # L-selectin, interacts with CD45RB in lymphocyte migration
  "CD44",        # Involved in pathways with CD26 and CD45
  "NT5E"         # CD73, functionally cooperates with CD26 (DPP4) in extracellular ATP degradation
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
  
  # EFS and DFS can remain unused for now
  
  tibble(
    patient_id = substr(patient_id, 1, 12),
    OS_time = os_time,
    OS_status = os_status
  )
}

# Updated function to extract expression for multiple genes
extract_expression_multiple <- function(tsv_file, mapping_df, genes) {
  file_id <- basename(dirname(tsv_file))
  barcode_row <- mapping_df %>% filter(file_id == !!file_id)
  if(nrow(barcode_row) == 0) return(NULL)
  
  data <- suppressMessages(read_tsv(tsv_file, comment = "#"))
  
  # Filter only genes of interest
  data_filtered <- data %>% filter(gene_name %in% genes)
  
  if(nrow(data_filtered) == 0) return(NULL)
  
  tibble(
    patient_id = barcode_row$barcode[1],
    !!!set_names(data_filtered$tpm_unstranded, data_filtered$gene_name)
  )
}

# Base paths
base_path <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
mapping_file <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/fileid_barcode_association.csv"
plot_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/Survivalplots_Ligands_OS_after_tumor_contact"

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

