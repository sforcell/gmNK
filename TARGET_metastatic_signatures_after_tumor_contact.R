# Librerie
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# Directory
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"
output_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/ligand_vs_metastatic_after_tumor_contact"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Signature genes
ligand_genes <- c(
  "ADA",         # Ligando di CD26 (DPP4), adenosina deaminasi
  "CXCL10",      # Substrato e interagente funzionale di CD26
  "CD38",        # Interagente di PECAM1 (CD31)
  "CD177",       # Ligando di CD31 (PECAM1)
  "ENG",         # Endoglina, interagente con CD31 (PECAM1)
  "VEGFA",       # Induce segnalazione tramite PECAM1/CD31
  "SPN",         # Sialoforina, interagente noto di CD69
  "CD69L",       # Ligando putativo di CD69 (non ben caratterizzato)
  "PTPRC",       # CD45, isoforma parentale di CD45RB
  "IL6ST",       # Interagente funzionale con CD45RB via segnalazione
  "CD53",        # Tetraspanina, modula l'interazione con CD2, CD4, CD5
  "CD2",         # Interagente di CD53 nei domini tetraspaninici
  "CD300A",      # Ligando/interagente di CD148 (PTPRJ)
  "PTPRJ",       # CD148, interagisce con recettori della famiglia CD300
  "TNFRSF1A",    # Interagente segnaletico modulato da CD148
  "HMGB1",       # Interagente immunomodulatore, coinvolto in pathway con CD69
  "LGALS9",      # Galectina-9, co-interagente con altri recettori in pathway con CD45RB
  "SELL",        # L-selectina, interagente con CD45RB in migrazione linfocitaria
  "CD44",        # Coinvolto in pathway con CD26 e CD45
  "NT5E"         # CD73, coopera funzionalmente con CD26 (DPP4) nella degradazione dell’ATP extracellulare
)

metastatic_genes <- c("SNAI1", "TWIST1", "ZEB1", "VIM", "MMP9", "CXCR4", "CD44",
                      "ITGA5", "ITGB1", "LOX", "TNC", "IL11", "TGFB1", "SPARC", "NEDD9")

# Funzione per leggere TPM
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

# Score medio di una signature
get_signature_score <- function(expr_mat, genes) {
  genes_present <- intersect(genes, rownames(expr_mat))
  if (length(genes_present) < 3) return(rep(NA, ncol(expr_mat)))
  colMeans(expr_mat[genes_present, , drop = FALSE])
}

# Lista per correlazioni
correlation_list <- list()

# Loop tumori
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

for (tumor in tumors) {
  cat("\nTumore:", tumor, "\n")
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
  ligand_score <- get_signature_score(expr_mat, ligand_genes)
  metastatic_score <- get_signature_score(expr_mat, metastatic_genes)

  if (all(is.na(ligand_score)) || all(is.na(metastatic_score))) next

  df <- data.frame(Sample = colnames(expr_mat),
                   Ligand = ligand_score,
                   Metastatic = metastatic_score)

  cor_val <- cor(df$Ligand, df$Metastatic, use = "complete.obs")
  p_val <- cor.test(df$Ligand, df$Metastatic)$p.value

  correlation_list[[length(correlation_list) + 1]] <- data.frame(
    Tumor = tumor,
    Comparison = "Ligand vs Metastatic",
    Correlation = cor_val,
    P_value = p_val
  )

  # Scatterplot
  p <- ggplot(df, aes(x = Ligand, y = Metastatic)) +
    geom_point(size = 2.5, alpha = 0.7, color = "#2c3e50") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
      title = paste("Ligand vs Metastatic in", tumor),
      subtitle = paste("Pearson r =", round(cor_val, 3), "/ p =", signif(p_val, 3)),
      x = "Ligand Signature (log10)",
      y = "Metastatic Signature (log10)"
    ) +
    theme_minimal(base_size = 14)

  filename <- file.path(output_dir, paste0(tumor, "_ligand_vs_metastatic.pdf"))
  ggsave(filename, plot = p, width = 7, height = 5)
  cat("Plot salvato:", filename, "\n")
}

# Heatmap finale
if (length(correlation_list) > 0) {
  cor_df <- do.call(rbind, correlation_list)

  cor_matrix <- reshape2::dcast(cor_df, Comparison ~ Tumor, value.var = "Correlation")
  rownames(cor_matrix) <- cor_matrix$Comparison
  cor_matrix <- cor_matrix[, -1]

  pval_matrix <- reshape2::dcast(cor_df, Comparison ~ Tumor, value.var = "P_value")
  rownames(pval_matrix) <- pval_matrix$Comparison
  pval_matrix <- pval_matrix[, -1]

  cor_matrix_masked <- cor_matrix
  cor_matrix_masked[pval_matrix > 0.05] <- NA

  cor_matrix_masked <- cor_matrix_masked[rowSums(!is.na(cor_matrix_masked)) > 0, ]
  cor_matrix_masked <- cor_matrix_masked[, colSums(!is.na(cor_matrix_masked)) > 0]

  pdf(file.path(output_dir, "ligand_vs_metastatic_correlation_heatmap.pdf"), width = 10, height = 2.5)
  pheatmap(cor_matrix_masked,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize = 12,
           main = "Ligand vs Metastatic Signature (TCGA)",
           na_col = "white")
  dev.off()

  cat("Heatmap salvata in:", file.path(output_dir, "ligand_vs_metastatic_correlation_heatmap.pdf"), "\n")
}

