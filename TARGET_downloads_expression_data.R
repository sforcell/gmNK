# conda activate TCGA_env

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

# Project base directory
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis"
data_dir <- file.path(base_dir, "GDCdata")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Retry download function (API mode, chunked download for stability)
download_with_retry <- function(query, tumor_name, max_attempts = 3) {
  for (i in seq_len(max_attempts)) {
    cat(sprintf("Attempt %d for %s\n", i, tumor_name))
    tryCatch({
      GDCdownload(query,
                  method = "api",
                  directory = data_dir,
                  files.per.chunk = 1)  # scarica uno alla volta per stabilità
      cat("Download successful for:", tumor_name, "\n")
      return(TRUE)
    }, error = function(e) {
      cat("Error at attempt", i, "per", tumor_name, ":", e$message, "\n")
      if (i == max_attempts) {
        cat("Final failure:", tumor_name, "\n")
      }
    })
  }
  return(FALSE)
}

# Loop over all tumors
for (tumor in tumors) {
  cat("Starting download for:", tumor, "\n")
  
  # Try STAR first, then fallback to HTSeq
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
        cat("Valid workflow found for", tumor, ":", workflow, "\n")
        break
      }
    }, error = function(e) {
      # workflow not available, try nextmo
    })
  }
  
  if (is.null(query)) {
    cat("No RNA-seq data available for:", tumor, "\n")
    next
  }
  
  # Download data
  success <- download_with_retry(query, tumor)
  
  # Organize files into tumor-specific folders
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

