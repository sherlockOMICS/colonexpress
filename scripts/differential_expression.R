# Load libraries
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)

# Load data
counts <- readRDS("data/processed/merged/counts.rds")

#
## Annotate genes upfront
#

# Remove Ensembl version suffixes from gene IDs
counts <- counts |>
  dplyr::mutate(gene_id = stringr::str_remove(gene_id, "\\.\\d+$"))

# Fetch gene annotations for all genes in the dataset
gene_annotations <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = unique(counts$gene_id),
    columns  = c("SYMBOL", "GENENAME"),
    keytype  = "ENSEMBL"
  ) |>
  tibble::as_tibble() |>
  dplyr::rename(
    gene_id          = ENSEMBL,
    gene_symbol      = SYMBOL,
    gene_description = GENENAME
  ) |>
  # Keep one row per gene (AnnotationDbi can return duplicates)
  dplyr::distinct(gene_id, .keep_all = TRUE)

#
## Format count matrix
#

gene_counts <-
  counts |>
  dplyr::select(samples, gene_id, counts) |>
  tidyr::pivot_wider(names_from = samples, values_from = counts) |>
  tibble::column_to_rownames("gene_id") |>
  as.matrix() |>
  na.omit()

# Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(gene_counts >= 10) >= 3
message("Keeping ", sum(keep), " of ", length(keep), " genes after filtering")
gene_counts <- gene_counts[keep, ]

#
## Sample metadata
#

samples <-
  counts |>
  dplyr::select(samples, tissue, disease_status) |>
  dplyr::distinct() |>
  dplyr::filter(samples %in% colnames(gene_counts)) |>
  dplyr::slice(match(colnames(gene_counts), samples)) |>
  dplyr::mutate(
    # Set factor levels: DESeq2 uses the first level as baseline
    disease_status = factor(disease_status, levels = c("healthy", "cancer"))
  )

stopifnot(
  "Column names must match sample IDs" =
    identical(colnames(gene_counts), dplyr::pull(samples, samples))
)

#
## DESeq2 analysis
#

de_deseq <- list()

# DE Step 1: Create DESeqDataSet object
# Uncomment if needed

# de_deseq$dds <- DESeq2::DESeqDataSetFromMatrix(
#   countData = gene_counts,
#   colData   = samples,
#   design    = ~ disease_status
# )


# DE Step 2: Run DESeq analysis
# Uncomment if needed: this step takes around 25 minutes in rise-cafe VM

# de_deseq$dds <- DESeq2::DESeq(de_deseq$dds)


#
## Normalised counts
#

normalised_counts <-
  DESeq2::counts(de_deseq$dds, normalized = TRUE) |>
  tibble::as_tibble(rownames = "gene_id") |>
  dplyr::left_join(gene_annotations, by = "gene_id") |>
  dplyr::relocate(gene_id, gene_symbol, gene_description)


#
## Differential expression results
#

# Find names of estimated coefficients
DESeq2::resultsNames(de_deseq$dds)

# Extract DE results
# The results function applies Benjamini-Hochberg correction by default
de_deseq$results$not_shrunk <- DESeq2::results(
  de_deseq$dds,
  name = "disease_status_cancer_vs_healthy"
)


# Apply log fold change shrinkage for more reliable estimates
# Requires the apeglm package
de_deseq$results$shrunk <- DESeq2::lfcShrink(
  de_deseq$dds,
  coef = "disease_status_cancer_vs_healthy",
  type = "apeglm"
)

# Summarise results
DESeq2::summary(de_deseq$results$shrunk)

#
## Select significant results with annotations
#

de_deseq$results$sig <-
  de_deseq$results$shrunk |>
  tibble::as_tibble(rownames = "gene_id") |>
  dplyr::filter(!is.na(padj), padj < 0.05) |>
  dplyr::arrange(padj, log2FoldChange) |>
  dplyr::left_join(gene_annotations, by = "gene_id") |>
  dplyr::mutate(
    gene_symbol      = dplyr::coalesce(gene_symbol, gene_id),
    gene_description = dplyr::coalesce(gene_description, gene_id)
  ) |>
  dplyr::relocate(gene_id, gene_symbol, gene_description)


#
## Save files to RDS objects
#

# Save dds object (it's large and takes a long time to compute)
if (!file.exists("data/processed/dds.rds")) saveRDS(de_deseq$dds, "data/processed/dds.rds")

# Save the normalised counts
if (!file.exists("data/processed/normalised_counts.rds")) saveRDS(normalised_counts, "data/processed/normalised_counts.rds")

# Save the significant results
if (!file.exists("data/processed/de_results.rds")) saveRDS(de_deseq$results, "data/processed/de_results.rds")

