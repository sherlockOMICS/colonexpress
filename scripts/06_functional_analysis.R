library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# Enrichment analysis (ORA)

de_results <- readRDS("data/processed/de_results.rds")
# normalised_counts <- readRDS("data/processed/normalised_counts.rds")

# ensembl2symbol <- normalised_counts |>
#   dplyr::select(gene_id, gene_symbol) |>
#   dplyr::distinct() |>
#   dplyr::mutate(
#     gene_name = dplyr::coalesce(gene_symbol, gene_id)
#     )

# Create a list to save the enrichment analysis results
fun_enrich <- list()

# Prepare list of significant DE genes in descending Log2FoldChange
fun_enrich$de_genes_fc <-
  de_results$sig |>
  dplyr::select(gene_id, log2FoldChange) |>
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Up regulated genes
fun_enrich$de_genes_up <-
  de_results$sig |>
  dplyr::select(gene_id, log2FoldChange) |>
  filter(log2FoldChange > 2) |>
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Down regulated genes
fun_enrich$de_genes_down <-
  de_results$sig |>
  dplyr::select(gene_id, log2FoldChange) |>
  filter(log2FoldChange < -2) |>
  dplyr::arrange(dplyr::desc(log2FoldChange))


#
## Functional enrichment analysis
#

# Run GO enrichment analysis using the enrichGO function
fun_enrich$ego_all <- clusterProfiler::enrichGO(
  gene = fun_enrich$de_genes_fc$gene_id,        # Genes of interest
  universe = rownames(de_results$shrunk),       # Background gene set
  OrgDb = org.Hs.eg.db,                         # Annotation database
  keyType = 'ENSEMBL',                          # Key type for gene identifiers
  readable = TRUE,                              # Convert gene IDs to gene names
  ont = "BP",                                   # Ontology: can be "BP", "MF", "CC", or "ALL"
  pvalueCutoff = 0.05,                          # P-value cutoff for significance
  qvalueCutoff = 0.10                           # Q-value cutoff for significance
)

# Visualize the enriched GO terms
fun_enrich$dotplot <- 
  enrichplot::dotplot(fun_enrich$ego, showCategory = 20, title = "GO BP | Enrichment barplot")

fun_enrich$heatplot <- 
  enrichplot::heatplot(fun_enrich$ego, showCategory = 10, 
                       foldChange = fun_enrich$de_genes_fc$log2FoldChange) +
  ggplot2::ggtitle("GO BP | Enrichment heatplot")

fun_enrich$emapplot <- 
  enrichplot::emapplot(enrichplot::pairwise_termsim(fun_enrich$ego), showCategory = 15, layout = "nicely")

fun_enrich$cnetplot <- 
  enrichplot::cnetplot(fun_enrich$ego, 
                       # categorySize= "pvalue", 
                       showCategory = 5, 
                       layout = "nicely", foldChange = fun_enrich$de_genes_fc$log2FoldChange)

# fun_enrich$treeplot <-
#   enrichplot::treeplot(
#     enrichplot::pairwise_termsim(fun_enrich$ego),
#     showCategory = 20,
#     nCluster = 3,
#     offset = rel(2)
#   ) +
#   ggplot2::ggtitle("GO BP | Enrichment treeplot") +
#   ggplot2::theme(text = element_text(size = 8))

# Combine the enrichment plots into panels from a single figure
fun_enrich$dotplot
fun_enrich$heatplot
fun_enrich$emapplot
fun_enrich$cnetplot



#
## Up regulated genes
#

# Run GO enrichment analysis using the enrichGO function
fun_enrich$ego_up <- clusterProfiler::enrichGO(
  gene = fun_enrich$de_genes_up$gene_id,        # Genes of interest
  universe = rownames(de_results$shrunk),       # Background gene set
  OrgDb = org.Hs.eg.db,                         # Annotation database
  keyType = 'ENSEMBL',                          # Key type for gene identifiers
  readable = TRUE,                              # Convert gene IDs to gene names
  ont = "BP",                                   # Ontology: can be "BP", "MF", "CC", or "ALL"
  pvalueCutoff = 0.05,                          # P-value cutoff for significance
  qvalueCutoff = 0.10                           # Q-value cutoff for significance
)



# Visualize the enriched GO terms
fun_enrich$dotplot_up <- 
  enrichplot::dotplot(fun_enrich$ego_up, showCategory = 20, title = "GO BP | Enrichment barplot")

fun_enrich$heatplot_up <- 
  enrichplot::heatplot(fun_enrich$ego_up, showCategory = 10, 
                       foldChange = fun_enrich$de_genes_fc$log2FoldChange) +
  ggplot2::ggtitle("GO BP | Enrichment heatplot")

fun_enrich$emapplot_up <- 
  enrichplot::emapplot(enrichplot::pairwise_termsim(fun_enrich$ego_up), showCategory = 15, layout = "nicely")

fun_enrich$cnetplot_up <- 
  enrichplot::cnetplot(fun_enrich$ego_up, 
                       categorySizeBy = "p.adjust", 
                       showCategory = 7, 
                       layout = "nicely", foldChange = fun_enrich$de_genes_up$log2FoldChange)

# Combine the enrichment plots into panels from a single figure
fun_enrich$dotplot_up
fun_enrich$heatplot_up
fun_enrich$emapplot_up
fun_enrich$cnetplot_up


#
## Down regulated genes
#






#
## PCA loadings
#

