# conda activate TCGA_env

# Tumori pediatrici TARGET
tumors <- c(
  "TARGET-NBL",
  "TARGET-RT",
  "TARGET-WT",
  "TARGET-CCSK",
  "TARGET-ALL-P3",
  "TARGET-ALL-P2",
  "TARGET-ALL-P1",
  "TARGET-AML",
  "TARGET-OS"
)

# Cartella base del progetto
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis"
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
  
  # Tenta con "STAR - Counts", poi fallback a "HTSeq - Counts"
  workflows <- c("STAR - Counts", "HTSeq - Counts")
  query <- NULL
  
  for (workflow in workflows) {
    tryCatch({
      query <- GDCquery(
        project = tumor,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = workflow
      )
      if (nrow(query$results[[1]]) > 0) {
        cat("Trovato workflow valido per", tumor, ":", workflow, "\n")
        break
      }
    }, error = function(e) {
      # workflow non disponibile, passa al prossimo
    })
  }
  
  if (is.null(query)) {
    cat("Nessun dato RNA-seq disponibile per", tumor, "\n")
    next
  }
  
  # Download dei dati
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

