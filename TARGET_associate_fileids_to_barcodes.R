# conda activate TCGA_env

library(TCGAbiolinks)

base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"

cat("Getting file_id from local folder...\n")
all_files <- list.files(base_dir, pattern = "\\.tsv(\\.gz)?$", 
                        full.names = TRUE, recursive = TRUE)

# Exstract file_id from the parent directory of the file
file_ids <- basename(dirname(all_files))
cat("Totale file_id trovati:", length(file_ids), "\n")

# Take all the found tumors
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

df_association <- data.frame()

for(tumor in tumors) {
  cat("Query metadata for tumors:", tumor, "\n")
  
  query <- GDCquery(
    project = tumor,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  results <- getResults(query)
  
  # Filter only the downloaded file_id
  results_filtered <- results[results$file_id %in% file_ids, ]
  
  df_association <- rbind(df_association, 
                          data.frame(file_id = results_filtered$file_id,
                                     barcode = results_filtered$cases.submitter_id,
                                     tumor = tumor,
                                     stringsAsFactors = FALSE))
}

cat("Obtained associations:\n")
print(head(df_association))

write.csv(df_association, file = "fileid_barcode_association.csv", 
          row.names = FALSE)

cat("Table with associations saved in 'fileid_barcode_association.csv'\n")
