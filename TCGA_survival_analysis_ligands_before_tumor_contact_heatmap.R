library(tidyverse)
library(survival)
library(survminer)
library(xml2)

# Lista geni
ligand_genes <- c(
  "THPO",        # Ligando di CD110 (MPL)
  "PDCD1LG2",    # Ligando alternativo di CD279 (PD-1), oltre a CD274
  "IGSF8",       # Interagente di CD81 (EWI-2)
  "CD48",        # Ligando di CD244 (2B4)
  "PECAM1",      # Ligando di CD38
  "COL1A1",      # Collagene, ligando di CD29 (ITGB1)
  "FN1",         # Fibronectina, ligando di CD29 e CD44
  "VCAM1",       # Ligando di CD29
  "LAMA1",       # Laminina, ligando di CD29
  "CD58",        # Ligando di CD2
  "HAS2",        # Sintetizza acido ialuronico, ligando di CD44
  "SPP1",        # Osteopontina, ligando alternativo di CD44
  "PPIA",        # Cyclophilin A, ligando di CD147 (BSG)
  "S100A9",      # Altro ligando/interagente di CD147 (BSG)
  "CD97",        # Ligando di CD55
  "ICAM1",       # Ligando di CD43 (SPN)
  "SELE",        # E-selectina, ligando di CD43 (SPN)
  "CD99"         # Ligando omofilico di CD99
)

# Percorsi
base_path <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
mapping_file <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/fileid_barcode_association.csv"

# Estrai dati clinici
extract_clinical <- function(xml_file) {
  doc <- read_xml(xml_file)
  ns <- xml_ns(doc)
  patient_id <- xml_text(xml_find_first(doc, ".//shared:bcr_patient_barcode", ns))
  if (length(patient_id) == 0 || is.na(patient_id)) return(NULL)
  vital_status <- tolower(xml_text(xml_find_first(doc, ".//clin_shared:vital_status", ns)))
  days_to_death <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_death", ns)))
  days_to_last_followup <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_last_followup", ns)))
  
  if (vital_status == "dead" && !is.na(days_to_death)) {
    os_time <- days_to_death
    os_status <- 1
  } else if (!is.na(days_to_last_followup)) {
    os_time <- days_to_last_followup
    os_status <- 0
  } else return(NULL)
  
  tibble(patient_id = substr(patient_id, 1, 12), OS_time = os_time, OS_status = os_status)
}

# Estrai espressione genica
extract_expression_multiple <- function(tsv_file, mapping_df, genes) {
  file_id <- basename(dirname(tsv_file))
  barcode_row <- mapping_df %>% filter(file_id == !!file_id)
  if(nrow(barcode_row) == 0) return(NULL)
  
  data <- suppressMessages(read_tsv(tsv_file, comment = "#"))
  
  data_filtered <- data %>%
    filter(gene_name %in% genes) %>%
    group_by(gene_name) %>%
    summarise(tpm_unstranded = mean(tpm_unstranded, na.rm = TRUE), .groups = "drop")
  
  if(nrow(data_filtered) == 0) return(NULL)
  
  tibble(patient_id = barcode_row$barcode[1], !!!set_names(data_filtered$tpm_unstranded, data_filtered$gene_name))
}

# Analisi
mapping_df <- read_csv(mapping_file, col_types = cols())
tumor_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
results_list <- list()

for (tumor_path in tumor_dirs) {
  tumor_name <- basename(tumor_path)
  cat("\nAnalizzando:", tumor_name, "\n")
  
  expression_files <- list.files(
    file.path(tumor_path, "Transcriptome_Profiling", "Gene_Expression_Quantification"),
    pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  
  clinical_files <- list.files(
    file.path(tumor_path, "Clinical", "Clinical_Supplement"),
    pattern = "\\.xml$", recursive = TRUE, full.names = TRUE)
  
  if (length(expression_files) == 0 || length(clinical_files) == 0) {
    cat("  → Dati mancanti, salto.\n")
    next
  }
  
  clinical_data <- bind_rows(lapply(clinical_files, extract_clinical)) %>%
    filter(!is.na(OS_time) & !is.na(OS_status))
  
  expression_list <- lapply(expression_files, extract_expression_multiple, mapping_df = mapping_df, genes = ligand_genes)
  expression_data <- bind_rows(expression_list[!sapply(expression_list, is.null)])
  
  merged_data <- inner_join(clinical_data, expression_data, by = "patient_id")
  if (nrow(merged_data) < 20) {
    cat("  → Troppi pochi pazienti validi, salto.\n")
    next
  }

  for (gene in ligand_genes) {
    if (!(gene %in% colnames(merged_data))) next
    gene_data <- merged_data %>% filter(!is.na(.data[[gene]]))
    if (nrow(gene_data) < 20) next
    
    q20 <- quantile(gene_data[[gene]], 0.20, na.rm = TRUE)
    q80 <- quantile(gene_data[[gene]], 0.80, na.rm = TRUE)
    
    gene_data <- gene_data %>%
      mutate(Gene_group = case_when(
        .data[[gene]] <= q20 ~ "Low",
        .data[[gene]] >= q80 ~ "High",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(Gene_group))
    
    if (nrow(gene_data) < 10) next
    
    surv_obj <- Surv(gene_data$OS_time, gene_data$OS_status)
    model <- try(coxph(surv_obj ~ Gene_group, data = gene_data), silent = TRUE)
    if (inherits(model, "try-error")) next
    
    model_sum <- summary(model)
    hr <- model_sum$coefficients[,"exp(coef)"][1]
    pval <- model_sum$coefficients[,"Pr(>|z|)"][1]
    
    results_list[[length(results_list)+1]] <- tibble(
      tumor = tumor_name,
      gene = gene,
      HR = hr,
      logHR = log10(hr),
      p_value = pval
    )
  }
}

# Combina risultati
results_df <- bind_rows(results_list) %>%
  mutate(significance = ifelse(p_value < 0.05, "*", ""))

  # Heatmap
  heatmap_plot <- ggplot(results_df, aes(x = gene, y = tumor, fill = logHR)) +
    geom_tile(color = "white") +
    geom_text(aes(label = significance), size = 5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "log10(HR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Survival trend (High vs Low expression)",
       x = "Ligands", y = "Tumors")

# Salva o mostra
ggsave("Survival_Heatmap_Ligands_OS_before_tumor_contact.pdf", heatmap_plot, width = 10, height = 8)
print(heatmap_plot)

