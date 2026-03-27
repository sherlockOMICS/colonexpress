library(SummarizedExperiment)
library(tibble)
library(tidyr)
library(dplyr)
library(zeallot)
library(here)

#
## TCGA
#

# Process TCGA data and assign only final tables to the global environment
c(tcga_counts, tcga_tpm) %<-% local({

  # Load TCGA colon (COAD) and rectum (READ) data into an isolated environment to avoid exposing 'data' object
  tcga_coad <- local({
    load(here::here("data/processed/TCGA/TCGA-COAD.rds"))
    data
  })
  tcga_read <- local({
    load(here::here("data/processed/TCGA/TCGA-READ.rds"))
    data
  })
  
  # Group datasets by tissue name for iteration
  tcga_list <- list(colon = tcga_coad, rectum = tcga_read)
  
  # Initialize empty lists to store counts and TPM tables per tissue
  counts_list <- list()
  tpm_list    <- list()
  
  # Loop over each tissue and extract counts and TPM assays
  for (name in names(tcga_list)) {
    counts_list[[name]] <- SummarizedExperiment::assay(tcga_list[[name]], "unstranded") |>
      tibble::as_tibble(rownames = "gene_id") |>
      tidyr::pivot_longer(-1, names_to = "samples", values_to = "counts") |>
      dplyr::mutate(tissue = as.factor(name), disease_status = as.factor("cancer"))
    
    tpm_list[[name]] <- SummarizedExperiment::assay(tcga_list[[name]], "tpm_unstrand") |>
      tibble::as_tibble(rownames = "gene_id") |>
      tidyr::pivot_longer(-1, names_to = "samples", values_to = "tpm") |>
      dplyr::mutate(tissue = as.factor(name), disease_status = as.factor("cancer"))
  }
  
  # Return combined tables for all tissues
  list(dplyr::bind_rows(counts_list),
       dplyr::bind_rows(tpm_list))
})

#
## GTEx
#

# Process GTEx data and assign only final tables to the global environment
c(gtex_counts, gtex_tpm) %<-% local({
  
  # Load GTEx data
  gtex <- readRDS(here::here("data/processed/GTEx/gtex.RDS"))
  
  # Initialize empty lists to store counts and TPM tables per tissue
  counts_list <- list()
  tpm_list    <- list()

  # Loop over each GTEx element and process counts (reads_) and TPMs (tpm_) separately
  for (name in names(gtex)) {
    if (grepl("^reads_", name)) {
      tissue <- gsub("reads_", "", name)
      counts_list[[name]] <- gtex[[name]] |>
        dplyr::rename(gene_id = Name, gene_name = Description) |>
        tidyr::pivot_longer(-c(1,2), names_to = "samples", values_to = "counts") |>
        dplyr::mutate(tissue = as.factor(tissue), disease_status = as.factor("healthy"))
      
    } else if (grepl("^tpm_", name)) {
      tissue <- gsub("tpm_", "", name)
      tpm_list[[name]] <- gtex[[name]] |>
        dplyr::rename(gene_id = Name, gene_name = Description) |>
        tidyr::pivot_longer(-c(1,2), names_to = "samples", values_to = "tpm") |>
        dplyr::mutate(tissue = as.factor(tissue), disease_status = as.factor("healthy"))
    }
  }
  
  # Return combined tables for all tissues
  list(dplyr::bind_rows(counts_list),
       dplyr::bind_rows(tpm_list))
})

# Combine TCGA and GTEx counts and TPM into single tables
counts <- dplyr::bind_rows(tcga_counts, gtex_counts |> dplyr::select(-gene_name))
tpm    <- dplyr::bind_rows(tcga_tpm, gtex_tpm |> dplyr::select(-gene_name))

# Create merged directory if it doesn't exist
if (!dir.exists(here::here("data/processed/merged"))) {
  dir.create(here::here("data/processed/merged"), recursive = TRUE, showWarnings = FALSE)
}

# Save merged counts and TPM tables
counts_path <- here::here("data/processed/merged/counts.rds")
tpm_path <- here::here("data/processed/merged/tpm.rds")

if (!file.exists(counts_path)) saveRDS(counts, counts_path)
if (!file.exists(tpm_path)) saveRDS(tpm, tpm_path)

