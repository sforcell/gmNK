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

base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"

# Function: download with retry mechanism
download_with_retry <- function(query, tumor_name, max_attempts = 3) {
  for (i in seq_len(max_attempts)) {
    cat(sprintf("Attempt %d for %s (clinical data)\n", i, tumor_name))
    tryCatch({
      GDCdownload(query, method = "api", directory = base_dir)
      cat("Download successful for:", tumor_name, "(clinical data)\n")
      return(TRUE)
    }, error = function(e) {
      cat("Error at attempt", i, "for", tumor_name, "(clinical data):", e$message, "\n")
      if (i == max_attempts) {
        cat("Final failure:", tumor_name, "(clinical data)\n")
      }
    })
  }
  return(FALSE)
}

# Loop over tumors
for (tumor in tumors) {
  cat("Starting clinical data download for:", tumor, "\n")
  
  # Query specifica solo per Clinical Supplement
  query_clin <- GDCquery(
    project = tumor,
    data.category = "Clinical",
    data.type = "Clinical Supplement"
  )
  
  success_clin <- download_with_retry(query_clin, tumor)
  
  if (success_clin) {
    cat("Download completed for:", tumor, "\n\n")
  }
}

