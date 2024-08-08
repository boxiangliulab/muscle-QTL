# Muscle bulk - 04. Enrichment
# Load necessary libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ReactomePA)

# Ensure the IDs in up_down1 do not contain spaces
up_down1$ID <- str_replace_all(up_down1$ID, " ", "")

# Separate up- and down-regulated genes
up_genes <- up_down1[grepl("Up", up_down1$gene_type),]
down_genes <- up_down1[grepl("Down", up_down1$gene_type),]

up_test <- up_genes[, c(1, 2)]
down_test <- down_genes[, c(1, 2)]

# Set databases
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

# Convert gene symbols to ENTREZ IDs for down-regulated genes
gene_down <- bitr(down_genes$ID, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

# GO enrichment analysis for down-regulated genes
GO_down <- enrichGO(gene_down$ENTREZID,
                    OrgDb = GO_database,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# Plot GO enrichment heatmap
enrichplot::heatplot(GO_down, showCategory = 20)

# KEGG enrichment analysis for down-regulated genes
KEGG_down <- enrichKEGG(gene_down$ENTREZID,
                        organism = KEGG_database,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

# GSEA analysis for up-regulated genes
colnames(up_test)[1] <- "SYMBOL"
info_merge <- merge(up_test, gene_down, by = "SYMBOL")
GSEA_input <- info_merge$logFC
names(GSEA_input) <- info_merge$ENTREZID
GSEA_input <- sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.8)

# Plot GSEA enrichment plot
enrichplot::gseaplot2(GSEA_KEGG, geneSetID = "hsa00010", title = "GSEA Enrichment Plot")

# Barplot and dotplot for GO enrichment results
barplot(GO_down, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
dotplot(GO_down, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")

# Barplot for KEGG enrichment results
barplot(KEGG_down, showCategory = 40, title = "KEGG Pathway")

# Reactome pathway enrichment analysis for down-regulated genes
bgGene <- anno$gene
bgGene <- bitr(bgGene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

reactome_enrich_down <- enrichPathway(gene = gene_down$ENTREZID,
                                      organism = 'human',
                                      universe = bgGene$ENTREZID,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.05,
                                      readable = TRUE)

# Plot Reactome enrichment results
dotplot(reactome_enrich_down, showCategory = 20, title = "Reactome Pathway Enrichment")
