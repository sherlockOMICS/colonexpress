library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Load data
counts <- readRDS(here::here("data/processed/merged/counts.rds"))

# Format count matrix
gene_counts <- counts |>
  dplyr::select(samples, gene_id, counts) |>
  tidyr::pivot_wider(names_from = samples, values_from = counts) |>
  tibble::column_to_rownames("gene_id") |>
  as.matrix() |>
  na.omit()

# Keep genes with at least 10 counts in at least 3 samples to remove lowly expressed genes
keep <- rowSums(gene_counts >= 10) >= 3
message("Keeping ", sum(keep), " of ", length(keep), " genes after filtering")
gene_counts <- gene_counts[keep, ]

# Extract unique sample annotations and align order with matrix columns
samples <- counts |>
  dplyr::select(samples, tissue, disease_status) |>
  dplyr::distinct() |>
  dplyr::filter(samples %in% colnames(gene_counts)) |>
  dplyr::slice(match(colnames(gene_counts), samples)) |>
  dplyr::mutate(
    # Set factor levels: DESeq2 uses the first level as baseline
    disease_status = factor(disease_status, levels = c("healthy", "cancer"))
  )

# Ensure sample order in metadata matches matrix column order
stopifnot(
  "Column names must match sample IDs" =
    identical(colnames(gene_counts), dplyr::pull(samples, samples))
)

# Create a list to store all QC results
qc <- list()

# VST normalization
qc$vst <- DESeq2::vst(gene_counts, blind = TRUE)

# Run PCA on transposed VST matrix (samples as rows, genes as columns)
qc$pca_vst <- prcomp(t(qc$vst))

# Extract PC coordinates for each sample
qc$components <- qc$pca_vst[["x"]]
qc$components <- tibble::as_tibble(qc$components, rownames = "samples")

# Extract sample annotations to color PCA points
sample_metadata <- counts |>
  dplyr::select(samples, tissue, disease_status) |>
  dplyr::distinct() |>
  dplyr::mutate(
    disease_status = factor(disease_status, levels = c("healthy", "cancer"))
    )

# Join PC coordinates with sample annotations
qc$components_annot <- dplyr::left_join(qc$components,
                                        sample_metadata,
                                        by = "samples") |>
  dplyr::relocate(tissue, disease_status, .after = samples)

# Calculate % variance explained by each PC
qc$pca_percent_var <- round(qc$pca_vst$sdev^2 / sum(qc$pca_vst$sdev^2) * 100)

# Compute Mahalanobis distance on PC1 and PC2 to detect multivariate outliers
pc_scores <- qc$components_annot |>
  dplyr::select(samples, PC1, PC2)

# Compute covariance matrix and centroid of the PC1/PC2 space
cov_mat <- cov(pc_scores[, c("PC1", "PC2")])
center  <- colMeans(pc_scores[, c("PC1", "PC2")])

# Flag samples exceeding the chi-squared threshold (p < 0.05, df = 2)
pc_scores <- pc_scores |>
  dplyr::mutate(
    mahal_dist = mahalanobis(cbind(PC1, PC2), center, cov_mat),
    is_outlier = mahal_dist > qchisq(0.95, df = 2)
  )

# Join Mahalanobis distance and outlier flag back to the annotated components
plot_data <- qc$components_annot |>
  dplyr::left_join(
    dplyr::select(pc_scores, samples, mahal_dist, is_outlier),
    by = "samples"
  )

# Report number of detected outliers
cat("According to the Mahalanobis distance cutoff, there are ", sum(plot_data$is_outlier), " outliers.")

#
## PCA plots
#

# Color by disease_status (healthy vs cancer)
qc$pca_disease_status <-
  ggplot2::ggplot(qc$components_annot, ggplot2::aes(x = PC1, y = PC2, color = disease_status)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA gene expression | Colored by disease status",
    x = paste0("PC1 (", qc$pca_percent_var[1], "% variance)"),
    y = paste0("PC2 (", qc$pca_percent_var[2], "% variance)"),
    color  = "Disease status"
  ) +
  ggplot2::theme_minimal()

# Color by tissue subtype
qc$pca_tissue <-
  ggplot2::ggplot(qc$components_annot, ggplot2::aes(x = PC1, y = PC2, color = tissue)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA gene expression | Colored by tissue",
    x = paste0("PC1 (", qc$pca_percent_var[1], "% variance)"),
    y = paste0("PC2 (", qc$pca_percent_var[2], "% variance)"),
    color  = "Tissue"
  ) +
  ggplot2::theme_minimal()

# Color by tissue and shape by disease_status to visualise both variables simultaneously
qc$pca_tissue_disease_status <-
  ggplot2::ggplot(qc$components_annot, ggplot2::aes(x = PC1, y = PC2,
                                                    color = tissue,
                                                    shape = disease_status)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA gene expression | Colored by tissue | Shape by disease status",
    x = paste0("PC1 (", qc$pca_percent_var[1], "% variance)"),
    y = paste0("PC2 (", qc$pca_percent_var[2], "% variance)"),
    color  = "Tissue",
    shape  = "Disease status"
  ) +
  ggplot2::theme_minimal()

#
## Sample-to-sample distance heatmap
#

# Calculate Euclidean distances between samples (t() transposes so samples are rows)
qc$sample_dist_matrix <- as.matrix(dist(t(qc$vst), method = "euclidean"))

# Define a green color palette for the distance heatmap (reversed: dark = similar)
qc$colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greens")))(255)

# Build annotation data frame for heatmap sidebars (rows and columns)
heatmap_annot <- samples |>
  tibble::column_to_rownames("samples") |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.factor)) |>
  dplyr::rename(
    "Disease Status" = disease_status,
    "Tissue"         = tissue
  )
  
# Plot distance heatmap with hierarchical clustering and annotation bars
qc$dist_clustering <- pheatmap::pheatmap(
  qc$sample_dist_matrix,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  col            = qc$colors,
  annotation_row = heatmap_annot,
  annotation_col = heatmap_annot,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames  = FALSE,
  show_colnames  = FALSE,
  fontsize       = 8
)

#
## Sample correlation heatmap
#

# Compute pairwise Pearson correlation between all samples
qc$sample_corr <- cor(qc$vst)

# Plot correlation heatmap with hierarchical clustering and annotation bars
qc$corr_clustering <- pheatmap::pheatmap(
  qc$sample_corr,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  annotation_row = heatmap_annot,
  annotation_col = heatmap_annot,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames  = FALSE,
  show_colnames  = FALSE,
  fontsize       = 8
)

# Save PCA plots
ggplot2::ggsave(
  filename = here::here("output/pca_disease_status.pdf"),
  plot     = qc$pca_disease_status,
  width    = 8,
  height   = 6
)

ggplot2::ggsave(
  filename = here::here("output/pca_tissue.pdf"),
  plot     = qc$pca_tissue,
  width    = 8,
  height   = 6
)

ggplot2::ggsave(
  filename = here::here("output/pca_tissue_disease_status.pdf"),
  plot     = qc$pca_tissue_disease_status,
  width    = 8,
  height   = 6
)

# Save heatmaps
pdf(here::here("output/dist_clustering.pdf"), width = 8, height = 6)
print(qc$dist_clustering)
dev.off()

pdf(here::here("output/corr_clustering.pdf"), width = 8, height = 6)
print(qc$corr_clustering)
dev.off()



# Export PCA object to disc
data_processed_path <- here::here("data/processed/pca_results.RDS")
if (!file.exists(data_processed_path)) saveRDS(qc[2:5], file = data_processed_path)




