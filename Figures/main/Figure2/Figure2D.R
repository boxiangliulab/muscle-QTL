# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(igraph)
library(dplyr)

# Set working directory for input/output files
setwd("/Users/wenjingwang/NUS/Liulab/muscleQTL/differential_expression/Plot/")

# Load and prepare candidate hub gene data
candidate_hub <- read.delim("candidate_hub.txt")
from_genes <- candidate_hub$fromNode
to_genes <- candidate_hub$toNode
gene_union <- union(from_genes, to_genes)

# Convert gene identifiers from ENSEMBL to ENTREZ
gene_entrez <- bitr(gene_union, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_entrez$ENTREZID

# Perform KEGG and GO enrichment analyses
kegg_result <- enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "kegg")
go_result <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "ENTREZID")

# Plot enrichment results using dotplot
pdf("KEGG_GO_Enrichment_Plots.pdf")
dotplot(kegg_result) + ggtitle("KEGG Enrichment Analysis")
dotplot(go_result) + ggtitle("GO Enrichment Analysis")
dev.off()

# Generate network graph from gene co-expression data
data <- read.delim("gene_network_data.txt")
network <- graph_from_data_frame(data, directed = FALSE)
V(network)$label <- V(network)$name
V(network)$size <- 25
V(network)$color <- "skyblue"
E(network)$width <- E(network)$weight * 10
E(network)$color <- "grey"

# Define various network layouts
layouts <- list(
  fr = layout_with_fr(network),
  circle = layout_in_circle(network),
  kk = layout_with_kk(network)
)

# Plot and save network graphs using different layouts
pdf("Gene_Network_Layouts.pdf")
for (name in names(layouts)) {
  plot(network, layout = layouts[[name]],
       vertex.label.cex = 0.8, vertex.label.color = "black",
       edge.width = E(network)$width, edge.color = E(network)$color,
       main = paste("Gene Co-expression Network (", name, " layout)", sep = ""))
}
dev.off()

# Retrieve and display gene information for a set of ENTREZ IDs
entrez_ids_example <- c("4541", "4536", "4540", "4509")
gene_info <- select(org.Hs.eg.db, keys = entrez_ids_example, columns = c("SYMBOL", "GENENAME"), keytype = "ENTREZID")

# Print gene information
print(gene_info)
