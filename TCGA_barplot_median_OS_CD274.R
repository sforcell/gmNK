# conda activate TCGA_env

# Definition:
# Median OS (overall survival) is the time point at which 50% of patients have died.
# If fewer than 50% of patients are deceased, the median is not reached.

library(tidyverse)
library(xml2)
library(survival)

# === 1. Load mapping file ===
mapping_file <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/fileid_barcode_association.csv"
mapping_df <- read_csv(mapping_file, col_types = cols())

# === 2. Functions ===
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

extract_expression <- function(tsv_file, mapping_df, gene_symbol = "CD274") {
  file_id <- basename(dirname(tsv_file))
  barcode_row <- mapping_df %>% filter(file_id == !!file_id)
  if (nrow(barcode_row) == 0) return(NULL)

  data <- suppressMessages(read_tsv(tsv_file, comment = "#"))
  gene_row <- data %>% filter(gene_name == gene_symbol)
  if (nrow(gene_row) == 0) return(NULL)

  tibble(patient_id = barcode_row$barcode[1], TPM = gene_row$tpm_unstranded[1])
}

get_median_OS <- function(data) {
  surv_obj <- Surv(data$OS_time, data$OS_status)
  fit <- survfit(surv_obj ~ 1)
  median_val <- summary(fit)$table["median"]
  if (is.na(median_val)) median_val <- NA  # Mediana non raggiunta
  return(as.numeric(median_val))
}

# === 3. Directories ===
base_path <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
plot_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/Survivalplots_CD274"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

tumor_dirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
summary_df <- tibble()

# === 4. Tumor loop ===
for (tumor_path in tumor_dirs) {
  tumor_name <- basename(tumor_path)
  cat("\nAnalysis for:", tumor_name, "\n")

  expression_files <- list.files(file.path(tumor_path, "Transcriptome_Profiling", "Gene_Expression_Quantification"),
                                 pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

  clinical_files <- list.files(file.path(tumor_path, "Clinical", "Clinical_Supplement"),
                               pattern = "\\.xml$", recursive = TRUE, full.names = TRUE)

  if (length(expression_files) == 0 || length(clinical_files) == 0) {
    cat("Missing data, skipping.\n")
    next
  }

  clinical_data <- bind_rows(lapply(clinical_files, extract_clinical)) %>%
    filter(!is.na(OS_time))

  expression_data <- bind_rows(
    lapply(expression_files, extract_expression, mapping_df = mapping_df, gene_symbol = "CD274")
  ) %>% filter(!is.na(TPM))

  merged_data <- inner_join(clinical_data, expression_data, by = "patient_id")

  if (nrow(merged_data) < 10) {
    cat("Insufficient data (n <", nrow(merged_data), "), skipping.\n")
    next
  }

  merged_data <- merged_data %>%
    mutate(group = ifelse(TPM >= median(TPM, na.rm = TRUE), "High", "Low"))

  summary <- merged_data %>%
    group_by(group) %>%
    group_modify(~ {
      tibble(median_OS = get_median_OS(.x))
    }) %>%
    mutate(tumor = tumor_name) %>%
    select(tumor, group, median_OS)

  print(summary)
  summary_df <- bind_rows(summary_df, summary)
}

# === 5. Final plot ===
if (nrow(summary_df) > 0) {
  p <- ggplot(summary_df, aes(x = tumor, y = median_OS, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_text(aes(label = ifelse(is.na(median_OS), "NR", round(median_OS))),
              vjust = -0.5, position = position_dodge(width = 0.8), size = 3) +
    scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#d62728")) +
    labs(title = "Kaplan-Meier Estimated Median OS by CD274 Expression",
         x = "Tumor Type",
         y = "Median OS (days)",
         fill = "Expression Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(plot_dir, "Kaplan_Meier_Median_OS_CD274_barplot.pdf"), p, width = 10, height = 5)
} else {
  cat("No data available for plotting.\n")
}

