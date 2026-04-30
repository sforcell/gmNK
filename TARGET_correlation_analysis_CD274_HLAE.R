# conda activate TCGA_env

# Librerie
library(data.table)
library(dplyr)
library(ggplot2)

# Directory
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"
plot_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/correlation_plots_CD274_HLAE"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Geni target
gene1 <- "CD274"
gene2 <- "HLA-E"

# Funzione: Lettura TPM da un file
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

# Funzione: Scatterplot
plot_gene_correlation <- function(expr_mat, gene1, gene2, tumor_name, plot_dir) {
  if (!(gene1 %in% rownames(expr_mat)) || !(gene2 %in% rownames(expr_mat))) {
    cat("Geni non trovati in", tumor_name, "\n")
    return(NULL)
  }

  df <- data.frame(
    Sample = colnames(expr_mat),
    CD274 = expr_mat[gene1, ],
    HLAE = expr_mat[gene2, ]
  )

  # Ordina i sample in base a CD274
  df <- df[order(df$CD274), ]
  df$Sample <- factor(df$Sample, levels = df$Sample)

  # Calcola correlazione
  cor_val <- cor(df$CD274, df$HLAE, method = "pearson")
  p_val <- cor.test(df$CD274, df$HLAE)$p.value

  # Plot (cerchi CD274 dietro, triangoli HLA-E davanti)
  p <- ggplot(df, aes(x = Sample)) +
    # CD274 sullo sfondo
    geom_point(aes(y = CD274, color = "CD274"), size = 2.5) +
    # HLA-E in primo piano
    geom_point(aes(y = HLAE * (max(df$CD274) / max(df$HLAE)), color = "HLA-E"), size = 2.5, shape = 17) +
    scale_y_continuous(
      name = "TPM CD274",
      sec.axis = sec_axis(~ . * (max(df$HLAE) / max(df$CD274)), name = "TPM HLA-E")
    ) +
    scale_color_manual(values = c("CD274" = "steelblue", "HLA-E" = "#e7423c")) +
    labs(
      title = paste("Correlation Plot: CD274 vs HLA-E in", tumor_name),
      subtitle = paste("Pearson r =", round(cor_val, 3), "/ p =", signif(p_val, 3)),
      x = "Samples (ordered by CD274 expression)",
      color = "Gene"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.line.y.left = element_line(color = "black"),
      axis.line.y.right = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15, face = "plain"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold")
    )

  # Salvataggio
  ggsave(
    filename = file.path(plot_dir, paste0(tumor_name, "_CD274_vs_HLAE_dotplot.pdf")),
    plot = p,
    width = 9, height = 5
  )

  cat("Plot finale salvato per", tumor_name, "\n")
}

# Loop su tumori
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

for (tumor in tumors) {
  cat("\nTumore:", tumor, "\n")
  tumor_path <- file.path(base_dir, tumor)
  files <- list.files(tumor_path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE, recursive = TRUE)

  if (length(files) == 0) {
    cat("Nessun file TPM trovato per:", tumor, "\n")
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
    cat("Dati insufficienti per:", tumor, "\n")
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

  # Plotta solo se entrambi i geni sono presenti
  plot_gene_correlation(expr_mat, gene1, gene2, tumor, plot_dir)
}

