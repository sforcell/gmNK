# conda activate TCGA_env

# Libreries
library(data.table)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

# Directories
base_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/GDCdata"
output_dir <- "/home/sergio/Scrivania/gmNK/TARGET_analysis/immunotype_results_ligands_after_tumor_contact"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Target genes and ligands
target_genes <- c("CD274", "HLA-E")

# Genes upregulated after tumor contact
ligand_genes <- c(
  # CD26 ligands
  "ADA",       # Adenosine deaminase (direct ligand)
  "CXCL10",    # Functional substrate cleaved by CD26
  "CXCL12",    # SDF-1, CD26 substrate
  "FN1",       # Fibronectin, binds CD26
  "LGALS3",    # Galectin-3, functional interaction
  "NT5E",      # CD73, metabolic cooperation (ATP → adenosine)

  # CD31 ligands
  "PECAM1",    # CD31 (homophilic interaction)
  "CD177",     # Direct ligand of CD31
  "ENG",       # Endoglin, interacts with PECAM1 in angiogenesis
  "ITGAV",     # Integrin αv (functional/indirect interaction)
  "ITGB3",     # Integrin β3

  # CD69 ligands
  "S100A8",    # Binds CD69 in tumor microenvironment
  "S100A9",    # Together with S100A8
  "SPN",       # Sialophorin/CD43, adhesion and regulation

  # CD45RB ligands
  "LGALS1",    # Galectin-1, modulates CD45RB
  "LGALS9",    # Galectin-9, immunoregulatory interaction
  "SELL",      # L-selectin, involved in migration

  # CD53 ligands
  "CD2",       # Involved in CD53 microdomains
  "HLA-A",     # MHC-I functional interaction
  "HLA-B",
  "HLA-C",

  # CD148 (PTPRJ) ligands
  "THBS1",     # Thrombospondin-1, direct extracellular ligand
  "CD300A"     # Reported interaction with PTPRJ
)

#ordered_genes <- c(target_genes, ligand_genes)
ordered_genes <- c(ligand_genes)

# Function to read a single sample file
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

# Loop over tumor types
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

  if (length(expr_list) < 5) next

  # Build expression matrix
  all_genes <- unique(unlist(lapply(expr_list, function(df) df$gene_name)))
  expr_mat <- matrix(0, nrow = length(all_genes), ncol = length(expr_list),
                     dimnames = list(all_genes, names(expr_list)))
  for (i in seq_along(expr_list)) {
    df <- expr_list[[i]]
    expr_mat[df$gene_name, i] <- df$tpm
  }

  # Immune score
  immune_genes <- c("CD8A", "CD8B", "CXCL9", "CXCL10", "IFNG", "GZMB", "PRF1")
  immune_found <- intersect(immune_genes, rownames(expr_mat))
  if (length(immune_found) < 3) next
  immune_score <- colMeans(expr_mat[immune_found, , drop = FALSE])
  quant <- quantile(immune_score, probs = c(0.2, 0.8))
  immune_class <- rep("intermediate", length(immune_score))
  immune_class[immune_score <= quant[1]] <- "cold"
  immune_class[immune_score >= quant[2]] <- "hot"
  names(immune_class) <- names(immune_score)

  # Select genes
  found_genes <- ordered_genes[ordered_genes %in% rownames(expr_mat)]
  if (length(found_genes) < 2) next
  selected_expr <- expr_mat[found_genes, , drop = FALSE]

  # Metadata
  meta <- data.frame(
    sample = colnames(selected_expr),
    tumor = tumor,
    immune_score = immune_score[colnames(selected_expr)],
    immune_class = immune_class[colnames(selected_expr)],
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample

  meta_filtered <- meta[meta$immune_class %in% c("hot", "cold"), ]
  if (nrow(meta_filtered) < 2) next
  expr_filtered <- selected_expr[, meta_filtered$sample, drop = FALSE]

  # Log transform and Z-score scaling
  heat_data_log <- log2(expr_filtered + 1)
  heat_data_log[is.na(heat_data_log)] <- 0
  heat_data_log[is.infinite(heat_data_log)] <- 0
  heat_data_scaled <- t(scale(t(heat_data_log)))
  heat_data_scaled[is.na(heat_data_scaled)] <- 0

  # Ordering
  meta_filtered <- meta_filtered[order(factor(meta_filtered$immune_class, levels = c("cold", "hot"))), ]
  heat_data_scaled <- heat_data_scaled[, meta_filtered$sample]

  anno_col <- data.frame(Types = meta_filtered$immune_class)
  rownames(anno_col) <- rownames(meta_filtered)

  ann_colors <- list(
    Types = c(cold = "#1f77b4", hot = "#d62728")
  )
  
  # Row annotation (ligands only)
  row_annotation <- data.frame(
    Groups = factor(
      c(
        #rep("Targets", length(target_genes[target_genes %in% found_genes])),
        rep("Ligands", length(ligand_genes[ligand_genes %in% found_genes]))
      ),
      #levels = c("Targets", "Ligands")
      levels = c("Ligands")
    )
  )
  rownames(row_annotation) <- found_genes

  row_colors <- list(
    Groups = c(
      #Targets = "#1b9e77",  # verde acqua
      Ligands = "#FF7F50"    # viola
    )
  )

  # Color scale
  max_abs <- max(abs(heat_data_scaled), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  legend_breaks <- pretty(c(-max_abs, max_abs), n = 5)

  # Heatmap
  total_width <- 10
  num_cols <- ncol(heat_data_scaled)
  cellwidth <- total_width / num_cols

  pheatmap(
    heat_data_scaled,
    annotation_col = anno_col,
    annotation_row = row_annotation,
    annotation_colors = c(ann_colors, row_colors),
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 12,
    main = paste("Z-score expression -", tumor),
    filename = file.path(output_dir, paste0("heatmap_", tumor, "_hot_vs_cold_and_ligand_target_genes.pdf")),
    width = total_width,
    height = 4 + 0.1 * length(found_genes),
    legend = TRUE,
    breaks = breaks,
    legend_breaks = legend_breaks,
    border_color = "grey",
    gaps_row = c(length(target_genes[target_genes %in% found_genes])),
    gaps_row_color = "white",
    cellwidth = cellwidth*30,
    cellheight = 10,
    annotation_names_row = FALSE
  )

  # Boxplots for ligand genes
  for (gene in ligand_genes) {
    if (!(gene %in% rownames(heat_data_log))) next
    plot_df <- data.frame(
      Expression = heat_data_log[gene, meta_filtered$sample],
      ImmuneClass = meta_filtered$immune_class
    )

    p <- ggplot(plot_df, aes(x = ImmuneClass, y = Expression, fill = ImmuneClass)) +
      geom_boxplot(outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
      scale_fill_manual(values = c("cold" = "#1f77b4", "hot" = "#d62728")) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        comparisons = list(c("cold", "hot")),
        tip.length = 0.01
      ) +
      labs(
        title = paste(gene, "-", tumor),
        x = NULL, y = "log2(TPM + 1)"
      ) +
      theme_minimal(base_size = 18) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
      )

    ggsave(
      filename = file.path(output_dir, paste0("boxplot_", gene, "_", tumor, "_hot_vs_cold_and_target_genes.pdf")),
      plot = p, width = 5, height = 5
    )
  }
}  

  # Split of genes
  row_split <- c(
    rep("Targets", length(target_genes[target_genes %in% found_genes])),
    rep("Ligands", length(ligand_genes[ligand_genes %in% found_genes]))
  )
  names(row_split) <- found_genes
  
  max_abs <- max(abs(heat_data_scaled), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  legend_breaks <- pretty(c(-max_abs, max_abs), n = 5)

  # Heatmap
  pheatmap(
    heat_data_scaled,
    annotation_col = anno_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    fontsize = 16,
    fontsize_row = 14,
    fontsize_col = 12,
    main = paste("Z-score immunity score -", tumor),
    filename = file.path(output_dir, paste0("heatmap_", tumor, "_hot_vs_cold_targets_and_ligands.pdf")),
    width = 18,
    height = 4 + 0.2 * length(found_genes),
    legend = TRUE,
    breaks = breaks,
    legend_breaks = legend_breaks,
    border_color = "grey",
    gaps_row = c(length(target_genes[target_genes %in% found_genes]))
  )

  # Boxplots for target and ligand genes
  for (gene in found_genes) {
    plot_df <- data.frame(
      Expression = heat_data_log[gene, meta_filtered$sample],
      ImmuneClass = meta_filtered$immune_class
    )

    p <- ggplot(plot_df, aes(x = ImmuneClass, y = Expression, fill = ImmuneClass)) +
      geom_boxplot(outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
      scale_fill_manual(values = c("cold" = "#1f77b4", "hot" = "#d62728")) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        comparisons = list(c("cold", "hot")),
        tip.length = 0.01
      ) +
      labs(
        title = paste(gene, "-", tumor),
        x = NULL, y = "log2(TPM + 1)"
      ) +
      theme_minimal(base_size = 18) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 16),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
      )

    ggsave(
      filename = file.path(output_dir, paste0("boxplot_", gene, "_", tumor, "_hot_vs_cold.pdf")),
      plot = p, width = 4, height = 5
    )
  }
}

