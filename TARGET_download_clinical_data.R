# conda activate TCGA_env

# Libreries
library(TCGAbiolinks)

# Pediatric TARGET tumors
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

# Output directory
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"

# Function with retry mechanism
download_with_retry <- function(query, tumor_name, max_attempts = 3) {
  for (i in seq_len(max_attempts)) {
    cat(sprintf("Attempt %d for %s (clinical)\n", i, tumor_name))
    tryCatch({
      GDCdownload(query, method = "api", directory = base_dir)
      cat("Download successful for:", tumor_name, "(clinical)\n")
      return(TRUE)
    }, error = function(e) {
      cat("Error at attempt", i, "per", tumor_name, "(clinical):", e$message, "\n")
      if (i == max_attempts) {
        cat("Final failure:", tumor_name, "(clinical)\n")
      }
    })
  }
  return(FALSE)
}

# Loop over TARGET tumors
for (tumor in tumors) {
  cat("\nStarting clinical data download for:", tumor, "\n")
  
  # Tentativo: prima con "Clinical Supplement", poi fallback a "Clinical"
  query_clin <- tryCatch({
    GDCquery(
      project = tumor,
      data.category = "Clinical",
      data.type = "Clinical Supplement"
    )
  }, error = function(e) {
    cat("'Clinical Supplement' not available for", tumor, ", switching to 'Clinical'\n")
    tryCatch({
      GDCquery(
        project = tumor,
        data.category = "Clinical",
        data.type = "Clinical"
      )
    }, error = function(e2) {
      cat("No clinical data available for", tumor, "\n")
      return(NULL)
    })
  })
  
  if (!is.null(query_clin)) {
    success_clin <- download_with_retry(query_clin, tumor)
    
    if (success_clin) {
      cat("Download completed for:", tumor, "\n")
    }
  }
}

