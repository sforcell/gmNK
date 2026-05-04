# conda activate TCGA_env

library(TCGAbiolinks)

base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"

cat("Extracting file_id from local directories...\n")
all_files <- list.files(base_dir, pattern = "\\.tsv(\\.gz)?$", 
                        full.names = TRUE, recursive = TRUE)

# Extract file_id from parent directory of each file
file_ids <- basename(dirname(all_files))
cat("Total file_id found:", length(file_ids), "\n")

# Get all tumor directories
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

df_association <- data.frame()

for(tumor in tumors) {
  cat("Querying metadata for tumor:", tumor, "\n")
  
  query <- GDCquery(
    project = tumor,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  results <- getResults(query)
  
   # Filter only downloaded file_ids
  results_filtered <- results[results$file_id %in% file_ids, ]
  
  df_association <- rbind(df_association, 
                          data.frame(file_id = results_filtered$file_id,
                                     barcode = results_filtered$cases.submitter_id,
                                     tumor = tumor,
                                     stringsAsFactors = FALSE))
}

cat("Associations found:\n")
print(head(df_association))

write.csv(df_association, file = "fileid_barcode_association.csv", 
          row.names = FALSE)

cat("Association table saved as 'fileid_barcode_association.csv'\n")

