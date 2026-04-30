# conda activate TCGA_env

# Librerie
library(TCGAbiolinks)

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

# Directory di salvataggio
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"

# Funzione con tentativi multipli
download_with_retry <- function(query, tumor_name, max_attempts = 3) {
  for (i in seq_len(max_attempts)) {
    cat(sprintf("Tentativo %d per %s (clinical)\n", i, tumor_name))
    tryCatch({
      GDCdownload(query, method = "api", directory = base_dir)
      cat("Download riuscito per:", tumor_name, "(clinical)\n")
      return(TRUE)
    }, error = function(e) {
      cat("Errore al tentativo", i, "per", tumor_name, "(clinical):", e$message, "\n")
      if (i == max_attempts) {
        cat("Fallito definitivamente:", tumor_name, "(clinical)\n")
      }
    })
  }
  return(FALSE)
}

# Loop per ciascun tumore TARGET
for (tumor in tumors) {
  cat("\nInizio download dati clinici per:", tumor, "\n")
  
  # Tentativo: prima con "Clinical Supplement", poi fallback a "Clinical"
  query_clin <- tryCatch({
    GDCquery(
      project = tumor,
      data.category = "Clinical",
      data.type = "Clinical Supplement"
    )
  }, error = function(e) {
    cat("'Clinical Supplement' non trovato per", tumor, ", provo con 'Clinical'\n")
    tryCatch({
      GDCquery(
        project = tumor,
        data.category = "Clinical",
        data.type = "Clinical"
      )
    }, error = function(e2) {
      cat("Nessun dato clinico disponibile per", tumor, "\n")
      return(NULL)
    })
  })
  
  if (!is.null(query_clin)) {
    success_clin <- download_with_retry(query_clin, tumor)
    
    if (success_clin) {
      cat("Download completato per:", tumor, "\n")
    }
  }
}

