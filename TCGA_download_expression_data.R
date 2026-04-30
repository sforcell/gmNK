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

# Cartella base del progetto
base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis"
data_dir <- file.path(base_dir, "GDCdata")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Funzione per tentare il download con retry e metodo API
download_with_retry <- function(query, tumor_name, max_attempts = 3) {
  for (i in seq_len(max_attempts)) {
    cat(sprintf("Tentativo %d per %s\n", i, tumor_name))
    tryCatch({
      GDCdownload(query,
                  method = "api",
                  directory = data_dir,
                  files.per.chunk = 1)  # scarica uno alla volta per stabilità
      cat("Download riuscito per:", tumor_name, "\n")
      return(TRUE)
    }, error = function(e) {
      cat("Errore al tentativo", i, "per", tumor_name, ":", e$message, "\n")
      if (i == max_attempts) {
        cat("Fallito definitivamente:", tumor_name, "\n")
      }
    })
  }
  return(FALSE)
}

# Loop su tutti i tumori
for (tumor in tumors) {
  cat("Inizio download per:", tumor, "\n")
  
  query <- GDCquery(
    project = tumor,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  success <- download_with_retry(query, tumor)
  
  # Organizza i file in sottocartelle per tumore
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

