# Libreries
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# Directories
base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
output_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/correlation_scatterplots_signatures_immunosuppression_infiltration"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Target genes and signatures
gene1 <- "CD274"
gene2 <- "HLA-E"

infiltration_genes <- c("KLRK1", "NCAM1", "KLRD1", "KLRC1", "KLRC2", "KLRB1", "GZMB", "GZMA", "GNLY", 
                        "PRF1", "FCGR3A", "NKG7", "XCL1", "XCL2", "TYROBP", "IL2RB", "ZAP70", 
                        "SH2D1B", "CD160", "CD244")

immunosuppression_genes <- c("PDCD1", "TIGIT", "HAVCR2", "CD96", "KLRC1", "LAG3", "CISH", "SOCS1", 
                             "SOCS3", "IL10RA", "TGFBR2", "EOMES", "FOXP3", "CEACAM1", "BTLA", 
                             "ZFP36", "PRDM1", "IKZF2")

immune_genes <- c("CD8A", "CD8B", "CXCL9", "CXCL10", "IFNG", "GZMB", "PRF1")

# Function to read TPM from a single file
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

# Function to compute a gene signature score
get_signature_score <- function(expr_mat, genes) {
  genes_present <- intersect(genes, rownames(expr_mat))
  if (length(genes_present) < 3) return(rep(NA, ncol(expr_mat)))
  colMeans(expr_mat[genes_present, , drop = FALSE])
}

# List to store correlations
correlation_list <- list()

# Loop over tumors
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

for (tumor in tumors) {
  cat("\nTumor:", tumor, "\n")
  tumor_path <- file.path(base_dir, tumor)
  files <- list.files(tumor_path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) next

  expr_list <- list()
  for (f in files) {
    df <- tryCatch(read_sample_file(f), error = function(e) NULL)
    if (is.null(df)) next
    sample_name <- tools::file_path_sans_ext(basename(f))
    sample_name <- sub("\\.tsv$", "", sample_name)
    sample_name <- sub("\\.gz$", "", sample_name)
    expr_list[[sample_name]] <- df
  }

  if (length(expr_list) < 3) next

  all_genes <- unique(unlist(lapply(expr_list, function(df) df$gene_name)))
  expr_mat <- matrix(NA, nrow = length(all_genes), ncol = length(expr_list),
                     dimnames = list(all_genes, names(expr_list)))
  for (i in seq_along(expr_list)) {
    df <- expr_list[[i]]
    expr_mat[df$gene_name, i] <- df$tpm
  }
  expr_mat[is.na(expr_mat)] <- 0

  # Signature scores
  infiltration_score <- get_signature_score(expr_mat, infiltration_genes)
  immunosuppression_score <- get_signature_score(expr_mat, immunosuppression_genes)
  immune_score <- get_signature_score(expr_mat, immune_genes)

  # Classification: hot / intermediate / cold
  immune_class <- rep("intermediate", length(immune_score))
  quant <- quantile(immune_score, probs = c(0.2, 0.8), na.rm = TRUE)
  immune_class[immune_score <= quant[1]] <- "cold"
  immune_class[immune_score >= quant[2]] <- "hot"
  names(immune_class) <- colnames(expr_mat)

  for (gene in c(gene1, gene2)) {
    if (!(gene %in% rownames(expr_mat))) next
    gene_expr <- expr_mat[gene, ]

    for (type in c("Infiltration", "Immunosuppression")) {
      score <- if (type == "Infiltration") infiltration_score else immunosuppression_score
      if (all(is.na(score))) next

      df <- data.frame(Sample = names(gene_expr), Gene = gene_expr, Signature = score, Class = immune_class)
      cor_val <- cor(df$Gene, df$Signature, use = "complete.obs")
      p_val <- cor.test(df$Gene, df$Signature)$p.value

      # Store correlation for heatmap
      correlation_list[[length(correlation_list) + 1]] <- data.frame(
        Tumor = tumor,
        Comparison = paste(gene, "vs", type),
        Correlation = cor_val,
        P_value = p_val
      )

      # Scatterplot
      p <- ggplot(df, aes(x = Signature, y = Gene, color = Class)) +
        geom_point(size = 2.5) +
        geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
        scale_color_manual(values = c("cold" = "#1f77b4", "intermediate" = "gray70", "hot" = "#d62728")) +
        scale_y_continuous(trans = "log10") +
        scale_x_continuous(trans = "log10") +
        labs(
          title = paste(gene, "vs", type, "in", tumor),
          subtitle = paste("Pearson r =", round(cor_val, 3), "/ p =", signif(p_val, 3)),
          x = paste(type, "score (log10)"),
          y = paste(gene, "TPM (log10)"),
          color = "Immune Class"
        ) +
        theme_minimal(base_size = 14)

      filename <- file.path(output_dir, paste0(tumor, "_", gene, "_vs_", type, ".pdf"))
      ggsave(filename, plot = p, width = 7, height = 5)
      cat("Plot saved:", filename, "\n")
    }
  }
}

# Generate final heatmap
if (length(correlation_list) > 0) {
  cor_df <- do.call(rbind, correlation_list)

  # Correlation matrix
  cor_matrix <- reshape2::dcast(cor_df, Comparison ~ Tumor, value.var = "Correlation")
  rownames(cor_matrix) <- cor_matrix$Comparison
  cor_matrix <- cor_matrix[, -1]

  # P-value matrix
  pval_matrix <- reshape2::dcast(cor_df, Comparison ~ Tumor, value.var = "P_value")
  rownames(pval_matrix) <- pval_matrix$Comparison
  pval_matrix <- pval_matrix[, -1]

  # Mask non-significant values
  cor_matrix_masked <- cor_matrix
  cor_matrix_masked[pval_matrix > 0.05] <- NA

  # Remove rows/columns with all NA
  cor_matrix_masked <- cor_matrix_masked[rowSums(!is.na(cor_matrix_masked)) > 0, ]
  cor_matrix_masked <- cor_matrix_masked[, colSums(!is.na(cor_matrix_masked)) > 0]

  # Plot heatmap (no text labels)
  pdf(file.path(output_dir, "TCGA_heatmap_correlations.pdf"), width = 10, height = 2.5)
  pheatmap(cor_matrix_masked,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize = 12,
           main = "TCGA database",
           na_col = "white")
  dev.off()

  cat("Heatmap saved in:", file.path(output_dir, "TCGA_heatmap_correlations.pdf"), "\n")
}


