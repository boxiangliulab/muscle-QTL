# Install and load necessary packages
if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  install.packages("ComplexUpset")
}
library(ComplexUpset)
library(ggplot2)

# Assuming eGenes_preINpostFDR50, eGenes_postINpreFDR50, sig_lead_pre, and sig_lead_post
# are already loaded in your environment

# Extract 'gene' columns from each dataframe
genes_preINpostFDR50 <- eGenes_preINpostFDR50$gene
genes_postINpreFDR50 <- eGenes_postINpreFDR50$gene
genes_sig_lead_pre <- sig_lead_pre$gene
genes_sig_lead_post <- sig_lead_post$gene

# Create a dataframe to hold these gene sets
all_genes <- unique(c(genes_preINpostFDR50, genes_postINpreFDR50, genes_sig_lead_pre, genes_sig_lead_post))
gene_data <- data.frame(
  gene = all_genes,
  Sig_Lead_Pre = as.numeric(all_genes %in% genes_sig_lead_pre),
  Sig_Lead_Post = as.numeric(all_genes %in% genes_sig_lead_post),
  Pre_in_Post = as.numeric(all_genes %in% genes_preINpostFDR50),
  Post_in_Pre = as.numeric(all_genes %in% genes_postINpreFDR50)
)

# Generate an UpSet plot
upset_plot <- upset(
  gene_data,
  sets = c("Post_in_Pre", "Pre_in_Post", "Sig_Lead_Post", "Sig_Lead_Pre"),
  order.by = "freq",
  keep.order = TRUE,  # Maintain the order of the sets as provided
  sets.bar.color = c("#a6cee3", "#1f78b4", "#fb9a99", "#e31a1c"),
  main.bar.color = "black"
)

# Save the plot
ggsave("~/Desktop/upset_plot.pdf", plot = upset_plot, width = 10, height = 6)
