# conda activate TCGA_env

library(TCGAbiolinks)

tumors <- c(
  "TCGA-LGG", "TCGA-THYM", "TCGA-KIRC", "TCGA-PAAD", "TCGA-SKCM", 
  "TCGA-TGCT", "TCGA-BRCA", "TCGA-SARC", "TCGA-PCPG", "TCGA-UVM",
  "TCGA-ACC", "TCGA-HNSC", "TCGA-COAD", "TCGA-LAML", "TCGA-UCS", 
  "TCGA-GBM", "TCGA-LIHC", "TCGA-OV", "TCGA-KIRP", "TCGA-KICH", 
  "TCGA-READ", "TCGA-THCA", "TCGA-LUAD", "TCGA-MESO", "TCGA-CHOL", 
  "TCGA-PRAD", "TCGA-STAD", "TCGA-LUSC", "TCGA-DLBC", "TCGA-CESC", 
  "TCGA-ESCA", "TCGA-BLCA", "TCGA-UCEC"
)

# Project base folder
base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis"
data_dir <- file.path(base_dir, "GDCdata")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Function to attempt download with retry using API method
download_with_retry <- function(query, tumor_name, max_attempts = 3) {
  for (i in seq_len(max_attempts)) {
    cat(sprintf("Attempt %d for %s\n", i, tumor_name))
    tryCatch({
      GDCdownload(query,
                  method = "api",
                  directory = data_dir,
                  files.per.chunk = 1)  # download one file at a time for stability
      cat("Download successful for:", tumor_name, "\n")
      return(TRUE)
    }, error = function(e) {
      cat("Error at attempt", i, "for", tumor_name, ":", e$message, "\n")
      if (i == max_attempts) {
        cat("Definitively failed:", tumor_name, "\n")
      }
    })
  }
  return(FALSE)
}

# Loop over all tumors
for (tumor in tumors) {
  cat("Starting download for:", tumor, "\n")
  
  query <- GDCquery(
    project = tumor,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  success <- download_with_retry(query, tumor)
  
  # Organize files into tumor-specific subfolders
  if (success) {
    tumor_dir <- file.path(data_dir, tumor)
    dir.create(tumor_dir, showWarnings = FALSE, recursive = TRUE)
    
    files_downloaded <- getResults(query, cols = "file_id")
    for (fid in files_downloaded) {
      src <- file.path(data_dir, fid)
      if (file.exists(src)) {
        file.rename(src, file.path(tumor_dir, fid))
      }
    }
  }
}

