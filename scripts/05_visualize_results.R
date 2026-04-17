library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(grid)

#
## Load data
#

# Load DE results, normalised counts, and PCA results
de_results <- readRDS("data/processed/de_results.rds")
normalised_counts <- readRDS("data/processed/normalised_counts.rds")
pca_results <- readRDS("data/processed/pca_results.RDS")

# List to store all visualisation plots
de_plots <- list()

#
## Prepare normalised counts matrix for significant genes
#

# Extract gene IDs from significant DE results (sig means padj < 0.05)
sig_gene_ids <- de_results$sig$gene_id

# Get normalized counts for significant genes
sig_normalised_counts <- normalised_counts |>
  dplyr::filter(gene_id %in% sig_gene_ids) |>
  dplyr::select(-gene_description) |>
  dplyr::mutate(
    gene_symbol  = dplyr::coalesce(gene_symbol, gene_id),
    .row_mean    = rowMeans(dplyr::pick(where(is.numeric)))
  ) |>
  dplyr::select(-gene_id) |>
  dplyr::slice_max(
    order_by   = .row_mean,
    n          = 1,
    by         = gene_symbol,
    with_ties  = FALSE       # guarantees exactly one row per symbol
  ) |>
  dplyr::select(-.row_mean) |>
  tibble::column_to_rownames("gene_symbol") |>
  as.matrix()


# Scale each row: subtract mean and divide by SD.
# The 2 transpositions are required because, by default, scale applies to the columns.
sig_normalised_counts_scaled <- t(scale(t(sig_normalised_counts)))   # scale rows, not columns

# Inspect your actual range first
range(sig_normalised_counts_scaled[1:200, ], na.rm = TRUE)
quantile(sig_normalised_counts_scaled[1:200, ], c(0.01, 0.99), na.rm = TRUE)


# # Clip extreme z-scores to avoid outliers dominating the color scale
# sig_normalised_counts_scaled[sig_normalised_counts_scaled >  4] <-  4
# sig_normalised_counts_scaled[sig_normalised_counts_scaled < -4] <- -4


# Confirm range after clipping
range(sig_normalised_counts_scaled, na.rm = TRUE)

#
## Heatmap | Top 200 DE genes
#

set.seed(42)

de_plots$ht <- ComplexHeatmap::Heatmap(
  sig_normalised_counts_scaled[1:200, ],
  name         = "Exprs (z-score)",
  column_title = "Cancer vs Healthy | Top 200 DE genes",
  # Hierarchical clustering of rows and columns
  cluster_columns = TRUE,
  cluster_rows    = TRUE,
  # Split rows into 2 clusters via K-means
  row_km        = 3,
  row_title     = c("A", "B", "C"),
  row_title_rot = 90,
  row_gap       = grid::unit(2, "mm"),
  # Split columns into 3 clusters via K-means (cancer vs healthy x 2)
  column_km  = 3,
  column_gap = grid::unit(2, "mm"),
  border = "grey",
  na_col = "white",
  # Color scale based on clipped range [-4, 40]
  col = circlize::colorRamp2(c(-3, 0, 3), c("skyblue3", "white", "forestgreen")),
  # row_names_gp      = grid::gpar(fontsize = 2),
  show_row_names      = FALSE,
  rect_gp             = grid::gpar(col = "grey", lwd = 0.1),
  show_column_names   = FALSE
)

# Draw the heatmap
ComplexHeatmap::draw(de_plots$ht, heatmap_legend_side = "right")


# Volcano plot

# # Add a column with differential expression status and add gene symbol to the results
sig_res_annot <-
  de_results$sig |>
  dplyr::mutate(diffexpressed = case_when(
    log2FoldChange > 5 & padj < 0.0001 ~ 'upregulated',
    log2FoldChange < -5 & padj < 0.0001 ~ 'downregulated',
    TRUE ~ 'not_de')) |>
  dplyr::arrange(padj, log2FoldChange)

# 
# set.seed(42)
# sig_res_annot_plot <- dplyr::bind_rows(
#   dplyr::filter(sig_res_annot, diffexpressed != "not_de"),
#   dplyr::filter(sig_res_annot, diffexpressed == "not_de") |>
#     dplyr::slice_sample(n = 10000)
# ) 


# # Create a volcano plot using ggplot2
# de_plots$volcano_plot <-
#   ggplot(data = sig_res_annot, aes(
#     x = log2FoldChange,
#     y = -log10(padj),
#     col = diffexpressed))+
#   geom_point(size = 0.6) +
#   # geom_text_repel(data = filter(sig_res_annot,
#   #                               ((abs(log2FoldChange) > log2(8)) & (padj < -log10(0.05)))),
#   #                 aes(label = gene_symbol), size = 2.5, max.overlaps = Inf) +
#   ggtitle("DE genes | Cancer vs Healthy") +
#   geom_vline(xintercept = c(-5, 5), col = "black", linetype = 'dashed', linewidth = 0.2) +
#   geom_hline(yintercept = -log10(0.0001), col = "black", linetype = 'dashed', linewidth = 0.2) +
#   theme(plot.title = element_text(size = rel(1.25), hjust = 0.5),
#         axis.title = element_text(size = rel(1))) +
#   scale_color_manual(values = c("upregulated" = "red",
#                                 "downregulated" = "blue",
#                                 "not_de" = "grey")) +
#   labs(color = 'DE genes') +
#   xlim(-25, 25) +   # Caution: This hides some genes
#   ylim(0, 310) +  # Caution: This hides some genes
#   theme_light()


# Build the four-level significance factor
sig_res_annot_plot <- sig_res_annot |>
  dplyr::mutate(
    significance = dplyr::case_when(
      padj < 0.0001 & abs(log2FoldChange) > 5 ~ "padj < 0.0001 & |LFC| > 5",
      padj < 0.0001                            ~ "padj < 0.0001",
      abs(log2FoldChange) > 5                  ~ "|LFC| > 5",
      TRUE                                     ~ "Not significant"
    ),
    significance = factor(
      significance,
      levels = c(
        "padj < 0.0001 & |LFC| > 5",
        "padj < 0.0001",
        "|LFC| > 5",
        "Not significant"
      )
    )
  )

# Downsample the not-significant bulk
set.seed(42)
sig_res_annot_plot <- dplyr::bind_rows(
  dplyr::filter(sig_res_annot_plot, significance != "Not significant"),
  dplyr::filter(sig_res_annot_plot, significance == "Not significant") |>
    dplyr::slice_sample(n = 10000)
)

# Top genes to label: most significant from the top tier
top_genes <- sig_res_annot_plot |>
  dplyr::filter(significance == "padj < 0.0001 & |LFC| > 5") |>
  dplyr::slice_min(padj, n = 20)

# p-value at the FDR threshold (for the horizontal line)
pval_fdr_threshold <- sig_res_annot |>
  dplyr::filter(!is.na(padj), !is.na(pvalue)) |>
  dplyr::filter(padj < 0.0001) |>
  dplyr::slice_max(pvalue, n = 1) |>
  dplyr::pull(pvalue)

sig_colours <- c(
  "padj < 0.0001 & |LFC| > 5" = "tomato",
  "padj < 0.0001"              = "darkorange",
  "|LFC| > 5"                  = "steelblue",
  "Not significant"            = "grey70"
)

de_plots$volcano_plot <-
  ggplot(
    sig_res_annot_plot,
    aes(x = log2FoldChange, y = -log10(pvalue), colour = significance)
  ) +
    geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(
    xintercept = c(-5, 5),
    linetype   = "dashed",
    linewidth  = 0.3,
    colour     = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(pval_fdr_threshold),
    linetype   = "dashed",
    linewidth  = 0.3,
    colour     = "grey40"
  ) +
  # ggrepel::geom_text_repel(
  #   data           = top_genes,
  #   aes(label      = gene_symbol),
  #   size           = 2.8,
  #   max.overlaps   = 30,
  #   colour         = "black",
  #   segment.colour = "grey60",
  #   segment.size   = 0.2,
  #   show.legend    = FALSE
  # ) +
  scale_colour_manual(values = sig_colours) +
  scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x       = "log\u2082 fold change (cancer / healthy)",
    y       = "-log\u2081\u2080(adjusted p-value)",
    colour  = NULL,
    title   = "Differential expression: cancer vs healthy",
    caption = paste0(
      "Dashed lines: |LFC| = 5 and p-value at padj 0.0001 threshold. ",
      "Top 20 significant genes labelled."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    plot.title       = element_text(size = rel(1.1), hjust = 0.5),
    plot.caption     = element_text(size = rel(0.75), colour = "grey50"),
    panel.grid.minor = element_blank()
  )

# Print the volcano plot
de_plots$volcano_plot
