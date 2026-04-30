# conda activate TCGA_env

library(tidyverse)
library(xml2)
library(survival)
library(survminer)

# 1) File di mapping
mapping_file <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/fileid_barcode_association.csv"
mapping_df <- read_csv(mapping_file, col_types = cols())

# Funzione per estrarre dati clinici
extract_clinical <- function(xml_file) {
  doc <- read_xml(xml_file)
  ns <- xml_ns(doc)

  patient_id <- xml_text(xml_find_first(doc, ".//shared:bcr_patient_barcode", ns))
  if (length(patient_id) == 0 || is.na(patient_id)) return(NULL)

  vital_status <- tolower(xml_text(xml_find_first(doc, ".//clin_shared:vital_status", ns)))
  days_to_death <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_death", ns)))
  days_to_last_followup <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_last_followup", ns)))
  days_to_recurrence <- as.numeric(xml_text(xml_find_first(doc, ".//clin_shared:days_to_first_recurrence", ns)))

  # OS: tempo e status
  if (vital_status == "dead" && !is.na(days_to_death)) {
    os_time <- days_to_death
    os_status <- 1
  } else if (!is.na(days_to_last_followup)) {
    os_time <- days_to_last_followup
    os_status <- 0
  } else {
    return(NULL)
  }

  # EFS: primo evento tra recidiva o morte
  possible_events_efs <- c(days_to_recurrence, days_to_death)
  possible_events_efs <- possible_events_efs[!is.na(possible_events_efs)]

  if (length(possible_events_efs) == 0 && !is.na(days_to_last_followup)) {
    efs_time <- days_to_last_followup
    efs_status <- 0
  } else if (length(possible_events_efs) > 0) {
    efs_time <- min(possible_events_efs)
    efs_status <- 1
  } else {
    efs_time <- NA
    efs_status <- NA
  }

  # DFS: solo recidiva è evento, morte è censura
  if (!is.na(days_to_recurrence)) {
    dfs_time <- days_to_recurrence
    dfs_status <- 1
  } else if (!is.na(days_to_death)) {
    dfs_time <- days_to_death
    dfs_status <- 0
  } else if (!is.na(days_to_last_followup)) {
    dfs_time <- days_to_last_followup
    dfs_status <- 0
  } else {
    dfs_time <- NA
    dfs_status <- NA
  }

  tibble(
    patient_id = substr(patient_id, 1, 12),
    OS_time = os_time,
    OS_status = os_status,
    EFS_time = efs_time,
    EFS_status = efs_status,
    DFS_time = dfs_time,
    DFS_status = dfs_status
  )
}

# Funzione per estrarre TPM HLAE
extract_expression <- function(tsv_file, mapping_df) {
  file_id <- basename(dirname(tsv_file))
  barcode_row <- mapping_df %>% filter(file_id == !!file_id)
  if(nrow(barcode_row) == 0) return(NULL)
  data <- suppressMessages(read_tsv(tsv_file, comment = "#"))
  hlae_row <- data %>% filter(gene_name == "HLA-E")
  if(nrow(hlae_row) == 0) return(NULL)
  tibble(
    patient_id = barcode_row$barcode[1],
    HLAE_tpm = hlae_row$tpm_unstranded[1]
  )
}

# 2) Cartella con tutti i tumori
base_path <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
tumor_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)

# Directory per salvare i plot
plot_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/Survivalplots_HLAE"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# 3) Cicla su ogni tumore
for (tumor_path in tumor_dirs) {
  tumor_name <- basename(tumor_path)
  cat("\n### Analisi per", tumor_name, "###\n")
  
  expression_files <- list.files(
    file.path(tumor_path, "Transcriptome_Profiling", "Gene_Expression_Quantification"),
    pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE
  )
  
  clinical_files <- list.files(
    file.path(tumor_path, "Clinical", "Clinical_Supplement"),
    pattern = "\\.xml$", recursive = TRUE, full.names = TRUE
  )
  
  if (length(expression_files) == 0 || length(clinical_files) == 0) {
    cat("Dati mancanti per", tumor_name, "- salto.\n")
    next
  }
  
  clinical_data <- bind_rows(lapply(clinical_files, extract_clinical)) %>%
    filter(!is.na(OS_time) & !is.na(OS_status))
  
  expression_list <- lapply(expression_files, extract_expression, mapping_df = mapping_df)
  expression_list <- expression_list[!sapply(expression_list, is.null)]
  expression_data <- bind_rows(expression_list)
  
  merged_data <- inner_join(clinical_data, expression_data, by = "patient_id")
  cat("  → Pazienti validi:", nrow(merged_data), "\n")
  
  if(nrow(merged_data) > 0){
    # Calcolo quantili 20% e 80% e definizione gruppi
    merged_data <- merged_data %>%
      mutate(
        quantile_20 = quantile(HLAE_tpm, 0.20),
        quantile_80 = quantile(HLAE_tpm, 0.80),
        HLAE_group = case_when(
          HLAE_tpm < quantile_20 ~ "Low",
          HLAE_tpm > quantile_80 ~ "High",
          TRUE ~ "Intermediate"
        )
      ) %>%
      select(-quantile_20, -quantile_80)

    ## OS plot
    os_surv <- Surv(time = merged_data$OS_time, event = merged_data$OS_status)
    os_fit <- survfit(os_surv ~ HLAE_group, data = merged_data)
    os_plot <- ggsurvplot(os_fit, data = merged_data, pval = TRUE, risk.table = TRUE,
                          legend.labs = c("Low", "Intermediate", "High"),
                          palette = c("#1f77b4", "#7f7f7f", "#d62728"),
                          title = paste("Overall Survival in", tumor_name),
                          legend.title = "HLA-E group",
                          ggtheme = theme_minimal() + theme(
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 12),
                              axis.title = element_text(size = 14),
                              axis.text = element_text(size = 12)
                          ))
    ggsave(filename = file.path(plot_dir, paste0("HLAE_OS_", tumor_name, ".pdf")),
           plot = os_plot$plot, width = 4, height = 4, dpi = 300)

    ## EFS plot
    efs_data <- merged_data %>% filter(!is.na(EFS_time) & !is.na(EFS_status))
    if (nrow(efs_data) > 0) {
      efs_surv <- Surv(time = efs_data$EFS_time, event = efs_data$EFS_status)
      efs_fit <- survfit(efs_surv ~ HLAE_group, data = efs_data)
      efs_plot <- ggsurvplot(efs_fit, data = efs_data, pval = TRUE, risk.table = TRUE,
                            legend.labs = c("Low", "Intermediate", "High"),
                            palette = c("darkgreen", "gray", "orange"),
                            title = paste("Event-Free Survival in", tumor_name),
                            legend.title = "HLA-E group",
                            ggtheme = theme_minimal() + theme(
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 12),
                              axis.title = element_text(size = 14),
                              axis.text = element_text(size = 12)
                            ))
      ggsave(filename = file.path(plot_dir, paste0("HLAE_EFS_", tumor_name, ".pdf")),
             plot = efs_plot$plot, width = 4, height = 4, dpi = 300)
    } else {
      cat("  → Nessun dato sufficiente per EFS.\n")
    }

    ## DFS plot
    dfs_data <- merged_data %>% filter(!is.na(DFS_time) & !is.na(DFS_status))
    if (nrow(dfs_data) > 0) {
      dfs_surv <- Surv(time = dfs_data$DFS_time, event = dfs_data$DFS_status)
      dfs_fit <- survfit(dfs_surv ~ HLAE_group, data = dfs_data)
      dfs_plot <- ggsurvplot(dfs_fit, data = dfs_data, pval = TRUE, risk.table = TRUE,
                            legend.labs = c("Low", "Intermediate", "High"),
                            palette = c("purple", "gray", "brown"),
                            title = paste("Disease-Free Survival in", tumor_name),
                            legend.title = "HLA-E group",
                            ggtheme = theme_minimal() + theme(
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 12),
                              axis.title = element_text(size = 14),
                              axis.text = element_text(size = 12)
                            ))
      ggsave(filename = file.path(plot_dir, paste0("HLAE_DFS_", tumor_name, ".pdf")),
             plot = dfs_plot$plot, width = 4, height = 4, dpi = 300)
    } else {
      cat("  → Nessun dato sufficiente per DFS.\n")
    }

  } else {
    cat("  → Nessun dato sufficiente per analisi.\n")
  }
}

