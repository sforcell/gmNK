# Librerie
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Directory
base_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/GDCdata"
output_dir <- "/home/sergio/Scrivania/gmNK/TCGA_analysis/boxplot_expression_CD274_HLAE"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Geni target
target_genes <- c("CD274", "HLA-E")

# Lettura TPM da un file campione
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

# Raccoglie espressioni
expression_data <- data.frame(
  Tumor = character(), Sample = character(), Gene = character(), Expression = numeric(),
  stringsAsFactors = FALSE
)

# Loop tumori
tumors <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
for (tumor in tumors) {
  cat("Tumore:", tumor, "\n")
  tumor_path <- file.path(base_dir, tumor)
  files <- list.files(tumor_path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) next

  for (f in files) {
    df_sample <- tryCatch(read_sample_file(f), error = function(e) NULL)
    if (is.null(df_sample)) next

    sample_name <- tools::file_path_sans_ext(basename(f))
    sample_name <- sub("\\.tsv$", "", sample_name)
    sample_name <- sub("\\.gz$", "", sample_name)

    for (gene in target_genes) {
      if (gene %in% df_sample$gene_name) {
        expr_val <- df_sample$tpm[df_sample$gene_name == gene]
        expression_data <- rbind(expression_data,
          data.frame(Tumor = tumor, Sample = sample_name, Gene = gene, Expression = expr_val,
                     stringsAsFactors = FALSE)
        )
      }
    }
  }
}

if (nrow(expression_data) == 0) stop("Nessun dato raccolto.")

# Palette colore estesa (Pastel + Viridis soft style)
tumor_levels_all <- unique(expression_data$Tumor)
n_tumors <- length(tumor_levels_all)
palette_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_tumors)
names(palette_colors) <- tumor_levels_all

# Funzione di plot
plot_gene_boxplot <- function(data, gene_name) {
  df_gene <- data %>% filter(Gene == gene_name)

  # Ordina tumori per mediana
  tumor_medians <- df_gene %>%
    group_by(Tumor) %>%
    summarise(Median = median(Expression, na.rm = TRUE)) %>%
    arrange(Median)
  df_gene$Tumor <- factor(df_gene$Tumor, levels = tumor_medians$Tumor)

  # Boxplot
  p <- ggplot(df_gene, aes(x = Tumor, y = Expression, fill = Tumor)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
    scale_y_log10() +
    scale_fill_manual(values = palette_colors[tumor_medians$Tumor]) +
    labs(title = paste(gene_name, " expression across tumors"),
         x = "Tumor", y = "TPM (log10)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

    filename <- file.path(output_dir, paste0("Boxplot_", gene_name, ".pdf"))
    ggsave(filename, plot = p, width = 9, height = 5.5)
    cat("Plot salvato:", filename, "\n")
  }

# Crea plot per ogni gene
for (gene in target_genes) {
  plot_gene_boxplot(expression_data, gene)
}

