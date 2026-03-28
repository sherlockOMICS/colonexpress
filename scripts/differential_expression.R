# Load libraries
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)

# Load data
counts <- readRDS("data/processed/merged/counts.rds")

#
## Format data
#

# Pivot gene counts and transform to matrix
gene_counts <-
  counts |>
  dplyr::select(samples, gene_id, counts) |>
  tidyr::pivot_wider(
    names_from = samples,
    values_from = counts
  ) |>
  tibble::column_to_rownames("gene_id") |>
  as.matrix() |>
  na.omit()    # Remove genes that are not present in both datasets

# Pre-filter low count genes
# Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(gene_counts >= 10) >= 3
message(
  "Keeping ", sum(keep), " of ", length(keep),
  
  " genes after filtering low counts"
)
gene_counts <- gene_counts[keep, ]

# Samples data
samples <-
  counts |>
  dplyr::select(samples, tissue, disease_status) |>
  dplyr::distinct() |>
  dplyr::filter(samples %in% colnames(gene_counts))

# Ensure samples are ordered to match gene_counts columns
samples <-
  samples |>
  dplyr::slice(match(colnames(gene_counts), samples))

# Verify alignment
stopifnot(
  "Column names must match sample IDs" =
    identical(colnames(gene_counts), dplyr::pull(samples, samples))
)

# Set factor levels: DESeq2 uses the first level as baseline
samples <-
  samples |>
  dplyr::mutate(
    disease_status = factor(disease_status, levels = c("healthy", "cancer"))
  )

#
## Differential Expression Analysis
#

# Create list to store analysis objects
de_deseq <- list()

# DE Step 1: Create DESeqDataSet object
# Uncomment if needed
# 
#   de_deseq$dds <- DESeq2::DESeqDataSetFromMatrix(
#     countData = gene_counts,
#     colData = samples,
#     design = ~ disease_status
#   )

# DE Step 2: Run DESeq analysis
# Uncomment if needed: this step takes over 15 minutes
  # de_deseq$dds <- DESeq2::DESeq(de_deseq$dds)

# Save dds object (it's large and takes a long time to compute)
if (!file.exists("data/processed/dds.rds")) saveRDS(de_deseq$dds, "data/processed/dds.rds")

#
## Inspect the dds object 
#

# Check the design formula
DESeq2::design(de_deseq$dds)

# Check the sample info
SummarizedExperiment::colData(de_deseq$dds)

# Display first rows of raw counts
head(DESeq2::counts(de_deseq$dds))

# Display first rows of normalised counts
head(DESeq2::counts(de_deseq$dds, normalized = TRUE))

# Convert normalised counts to tibble
normalised_counts <- DESeq2::counts(de_deseq$dds, normalized = TRUE) |>
  tibble::as_tibble(rownames = "ensembl_gene_id")

# Save the normalised counts
if (!file.exists("data/processed/normalised_counts.rds")) saveRDS(normalised_counts, "data/processed/normalised_counts.rds")

head(normalised_counts)

#
## Extract Differential Expression results
#

# Find names of estimated coefficients
DESeq2::resultsNames(de_deseq$dds)

# Extract DE results
# The results function applies Benjamini-Hochberg correction by default
de_deseq$results$not_shrunk <- DESeq2::results(
  de_deseq$dds,
  name = "disease_status_cancer_vs_healthy"
)

# Summarise results
DESeq2::summary(de_deseq$results$not_shrunk)

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
## Select significant DE results
#

# Extract significant results (padj < 0.05) and convert to tibble
de_deseq$results$sig <-
  de_deseq$results$shrunk |>
  tibble::as_tibble(rownames = "ensembl_gene_id") |>
  dplyr::filter(!is.na(padj), padj < 0.05) |>
  dplyr::arrange(padj, log2FoldChange)

# View top results
head(de_deseq$results$sig)


#
## Add gene annotation to results
#

# Get the Ensembl IDs from the results
# Strip Ensembl version suffix from Ensembl ids
ensembl_ids <- 
  de_deseq$results$sig |>
  dplyr::pull(ensembl_gene_id) |>
  stringr::str_remove("\\.\\d+$")

# Get the gene annotations from the database
gene_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
) |>
  tibble::as_tibble() |>
  dplyr::rename(
    ensembl_gene_id = ENSEMBL,
    gene_symbol = SYMBOL,
    gene_description = GENENAME
  )

# Join annotations to the results
de_deseq$results$sig <-
  de_deseq$results$sig |>
  dplyr::mutate(ensembl_gene_id = stringr::str_remove(ensembl_gene_id, "\\.\\d+$")) |>
  dplyr::left_join(gene_annotations, by = "ensembl_gene_id") |>
  dplyr::mutate(
    gene_symbol = dplyr::coalesce(gene_symbol, ensembl_gene_id),
    gene_description = dplyr::coalesce(gene_description, ensembl_gene_id)
  ) |>
  dplyr::relocate(ensembl_gene_id, gene_symbol, gene_description)

# Save the results
if (!file.exists("data/processed/de_results.rds")) saveRDS(de_deseq$results, "data/processed/de_results.rds")


