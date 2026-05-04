# conda activate TCGA_env

# Libreries
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Directories
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"
output_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/correlation_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

target_genes <- c("CD274", "HLA-E")
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

# Function: Read and filter protein-coding genes
read_sample_file <- function(file_path) {
  df <- fread(file_path, data.table = FALSE)
  df <- df[!grepl("^N_", df$gene_id), ]
  df <- df[df$gene_type == "protein_coding", ]
  df <- df[, c("gene_name", "tpm_unstranded")]
  df$tpm_unstranded <- as.numeric(df$tpm_unstranded)
  
  df <- df %>%
    group_by(gene_name) %>%
    summarise(tpm = mean(tpm_unstranded, na.rm = TRUE)) %>%
    ungroup()
  
  return(df)
}

# Function: GO enrichment analysis + plots
run_go_analysis <- function(gene_symbols, out_file_prefix) {
  if (length(gene_symbols) == 0 || all(is.na(gene_symbols))) {
    cat("No valid genes for GO analysis:", out_file_prefix, "\n")
    return(NULL)
  }
  
  entrez_ids <- tryCatch({
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error during SYMBOL -> ENTREZID conversion:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    cat("No valid ENTREZ IDs found for:", out_file_prefix, "\n")
    return(NULL)
  }

  ego <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    readable = TRUE
  )
  
  go_df <- as.data.frame(ego)
  if (nrow(go_df) == 0) {
    cat("No significant GO terms found for:", out_file_prefix, "\n")
    return(NULL)
  }
  
  write.csv(go_df, paste0(out_file_prefix, "_GO_BP.csv"), row.names = FALSE)
  
  ego@result$Description <- gsub(" ", "_", ego@result$Description)

  # Barplot
  pdf(paste0(out_file_prefix, "_GO_BP_barplot.pdf"), width = 10, height = 8)
  print(barplot(ego, showCategory = 10, title = "GO Biological Process") +
          theme(axis.text.y = element_text(size = 10, hjust = 1)))
  dev.off()
  
  # Dotplot
  pdf(paste0(out_file_prefix, "_GO_BP_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(ego, showCategory = 10, title = "GO Biological Process") +
          theme(axis.text.y = element_text(size = 10, hjust = 1)))
  dev.off()
}

# Loop over tumor types
for (tumor in tumors) {
  cat("\nTumore:", tumor, "\n")
  tumor_path <- file.path(base_dir, tumor)
  
  files <- list.files(tumor_path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  
  if (length(files) == 0) {
    cat("No TPM files found for:", tumor, "\n")
    next
  }
  
  expr_list <- list()
  
  for (f in files) {
    df <- tryCatch(read_sample_file(f), error = function(e) NULL)
    if (is.null(df)) next
    
    sample_name <- tools::file_path_sans_ext(basename(f))
    sample_name <- sub("\\.tsv$", "", sample_name)
    sample_name <- sub("\\.gz$", "", sample_name)
    
    expr_list[[sample_name]] <- df
  }
  
  if (length(expr_list) < 3) {
    cat("Not enough data for:", tumor, "\n")
    next
  }
  
  all_genes <- unique(unlist(lapply(expr_list, function(df) df$gene_name)))
  
  expr_mat <- matrix(NA, nrow = length(all_genes), ncol = length(expr_list),
                     dimnames = list(all_genes, names(expr_list)))
  
  for (i in seq_along(expr_list)) {
    df <- expr_list[[i]]
    expr_mat[df$gene_name, i] <- df$tpm
  }
  
  expr_mat[is.na(expr_mat)] <- 0
  
  genes_found <- intersect(target_genes, rownames(expr_mat))
  if (length(genes_found) == 0) {
    cat("Target genes not found in:", tumor, "\n")
    next
  }
  
  for (target in genes_found) {
    cat("Correlation analysis with gene:", target, "\n")
    
    target_expr <- expr_mat[target, ]
    
    cor_vals <- apply(expr_mat, 1, function(gene_expr) {
      suppressWarnings(cor(gene_expr, target_expr, method = "pearson"))
    })
    
    p_vals <- apply(expr_mat, 1, function(gene_expr) {
      suppressWarnings(cor.test(gene_expr, target_expr, method = "pearson")$p.value)
    })
    
    cor_df <- data.frame(
      gene = rownames(expr_mat),
      cor = cor_vals,
      pval = p_vals,
      stringsAsFactors = FALSE
    )
    
    cor_df$adj_pval <- p.adjust(cor_df$pval, method = "BH")
    
    pos_cor <- subset(cor_df, cor > 0.3 & adj_pval < 0.05)
    neg_cor <- subset(cor_df, cor < -0.3 & adj_pval < 0.05)
    
    prefix <- file.path(output_dir, paste0(tumor, "_", target))
    
    write.csv(pos_cor, paste0(prefix, "_positive.csv"), row.names = FALSE)
    write.csv(neg_cor, paste0(prefix, "_negative.csv"), row.names = FALSE)
    
    if (nrow(pos_cor) > 0 && any(!is.na(pos_cor$gene))) {
      cat("Running GO analysis for positive correlations (n =", nrow(pos_cor), ")\n")
      run_go_analysis(pos_cor$gene, paste0(prefix, "_positive"))
    } else {
      cat("No significant positively correlated genes for GO analysis\n")
    }
    
    if (nrow(neg_cor) > 0 && any(!is.na(neg_cor$gene))) {
      cat("Running GO analysis for negative correlations (n =", nrow(neg_cor), ")\n")
      run_go_analysis(neg_cor$gene, paste0(prefix, "_negative"))
    } else {
      cat("No significant negatively correlated genes for GO analysis\n")
    }
    
    cat("Saved:", nrow(pos_cor), "positive and", nrow(neg_cor), "negative + GO + plots (if available)\n")
  }
}

