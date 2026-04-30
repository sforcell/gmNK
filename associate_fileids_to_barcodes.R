# conda activate TCGA_env

library(TCGAbiolinks)

base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"

cat("Estrazione file_id da cartelle locali...\n")
all_files <- list.files(base_dir, pattern = "\\.tsv(\\.gz)?$", 
                        full.names = TRUE, recursive = TRUE)

# Estrai file_id dalla directory genitore del file
file_ids <- basename(dirname(all_files))
cat("Totale file_id trovati:", length(file_ids), "\n")

# Prendi tutti i tumori trovati (cartelle principali)
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

df_association <- data.frame()

for(tumor in tumors) {
  cat("Query metadati per tumore:", tumor, "\n")
  
  query <- GDCquery(
    project = tumor,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  results <- getResults(query)
  
  # Filtra solo i file_id che hai scaricato
  results_filtered <- results[results$file_id %in% file_ids, ]
  
  df_association <- rbind(df_association, 
                          data.frame(file_id = results_filtered$file_id,
                                     barcode = results_filtered$cases.submitter_id,
                                     tumor = tumor,
                                     stringsAsFactors = FALSE))
}

cat("Associazioni trovate:\n")
print(head(df_association))

write.csv(df_association, file = "fileid_barcode_association.csv", 
          row.names = FALSE)

cat("Tabella di associazione salvata in 'fileid_barcode_association.csv'\n")

