library(TCGAbiolinks)
library(SummarizedExperiment)
library(here)

#
## TCGA
#

### TCGA data downloaded in 13 March 2026 

# Create directories
if(!dir.exists(here::here("data/raw/TCGA"))) {
  dir.create(here::here("data/raw/TCGA"), recursive = TRUE, showWarnings = FALSE)
}

if(!dir.exists(here::here("data/processed/TCGA"))) {
  dir.create(here::here("data/processed/TCGA"), recursive = TRUE, showWarnings = FALSE)
}

# Path
tcga_data_path <- here::here("data/raw/TCGA")

# Cancer types for colorectal cancer
cancer_types <- c("TCGA-COAD", "TCGA-READ")

# Lists to store results
tcga_query   <- list()
tcga_data_se <- list()

# Loop through cancer types
for (cancer_type in cancer_types) {
  
  message("Processing: ", cancer_type)
  
  cancer_type2 <- gsub("-", "_", cancer_type)
  
  # Query
  tcga_query[[cancer_type2]] <- TCGAbiolinks::GDCquery(
    project = cancer_type,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download
  TCGAbiolinks::GDCdownload(
    query = tcga_query[[cancer_type2]],
    method = "api",
    directory = tcga_data_path,
    files.per.chunk = 10
  )
  
  # Prepare SummarizedExperiment
  tcga_data_se[[cancer_type2]] <- TCGAbiolinks::GDCprepare(
    query = tcga_query[[cancer_type2]],
    directory = tcga_data_path,
    save = TRUE,
    save.filename = here::here(
      "data/processed/TCGA",
      paste0(cancer_type, ".rds")
    )
  )
}

# Check tissue of origin
lapply(tcga_data_se, function(x) {
  table(SummarizedExperiment::colData(x)$tissue_or_organ_of_origin, useNA = "ifany")
})

# Save lists
saveRDS(tcga_query, here::here("data/raw/TCGA/tcga_data_query.rds"))
saveRDS(tcga_data_se, here::here("data/processed/TCGA/tcga_data_se.rds"))

#
## GTEx
#

### GTEx data downloaded in 13 March 2026 
### GTEx files download from: https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression

# Create directories
   # Note: No need to create counts and tpm directories because the file names are different:
   # the counts are named "_reads_" and the others are named "_tpm_". This is why the download
   # was failing.

if(!dir.exists(here::here("data/raw/GTEx"))) {
  dir.create(here::here("data/raw/GTEx"), recursive = TRUE, showWarnings = FALSE)
}

if(!dir.exists(here::here("data/processed/GTEx"))) {
  dir.create(here::here("data/processed/GTEx"), recursive = TRUE, showWarnings = FALSE)
}


# Base URL for the data to download
gtex_base_url1 <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/counts-by-tissue/"
gtex_base_url2 <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/tpms-by-tissue/"

# File names to download
gtex_file_names <- c("gene_reads_v11_colon_sigmoid.gct.gz",
                     "gene_reads_v11_colon_transverse.gct.gz", 
                     "gene_reads_v11_colon_transverse_mixed_cell.gct.gz",
                     "gene_reads_v11_colon_transverse_mucosa.gct.gz",
                     "gene_reads_v11_colon_transverse_muscularis.gct.gz",
                     "gene_tpm_v11_colon_sigmoid.gct.gz",
                     "gene_tpm_v11_colon_transverse.gct.gz",
                     "gene_tpm_v11_colon_transverse_mixed_cell.gct.gz",
                     "gene_tpm_v11_colon_transverse_mucosa.gct.gz",
                     "gene_tpm_v11_colon_transverse_muscularis.gct.gz"
)

# Loop over the file names and download each file
for (file_name in gtex_file_names[1:5]) {
  download.file(paste0(gtex_base_url1, file_name),
                destfile = here("data/raw/GTEx/", file_name))
}

for (file_name in gtex_file_names[6:10]) {
  download.file(paste0(gtex_base_url2, file_name),
                destfile = here("data/raw/GTEx/", file_name))
}


### Read GTEx data and save in a RDS file

## Create list to save all gtex data
gtex_data <- list()

## Loop over the file names and read each file
for (file_name in gtex_file_names) {
  
  ## Parse the file name to get the tissue name
  tissue_name <- gsub(pattern = "^gene_(.*)_v11_(.*)\\.gct\\.gz$", 
                      replacement = "\\1_\\2", 
                      file_name)
  
  # Create the file path
  gtex_file_path <- here("data/raw/GTEx", file_name)
  
  # Read the file and save it to a list (skip the first 2 lines that are not relevant data)
  gtex_data[[tissue_name]] <- readr::read_delim(gtex_file_path, col_names = TRUE, 
                                                delim = "\t", skip = 2, 
                                                comment = "",
                                                quote = "",
                                                name_repair = "unique")
}

## Save the GTEx data object
saveRDS(gtex_data, file = here("data/processed/GTEx/gtex.RDS"))



