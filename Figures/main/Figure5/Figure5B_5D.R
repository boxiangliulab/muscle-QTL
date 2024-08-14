# Load necessary libraries
library(ggplot2)
library(gcookbook)
library(dplyr)

# Inputs
target_rs_id <- "target_rs_id"  # Target SNP ID
gene_data_file <- "gene_data_file"  # Path to the gene data file
ld_matrix_file <- "ld_matrix_file"  # Path to the LD matrix file
output_file_path <- "path_to_output"  # Output file path

# Read gene data and LD matrix
eQTL_data <- read.table(gene_data_file, header = TRUE, stringsAsFactors = FALSE)
LD_matrix <- read.table(ld_matrix_file, header = TRUE, stringsAsFactors = FALSE)

# Process eQTL data
eQTL_data$position <- as.numeric(eQTL_data$position)
eQTL_data$distance <- abs(eQTL_data$position - eQTL_data$position[eQTL_data$rs_id == target_rs_id])
sorted_data <- eQTL_data[order(eQTL_data$distance, eQTL_data$pval_nominal_sqtl), ]
top_300_rs_id <- head(sorted_data$rs_id, 300)

# Write top 300 rs_id to a file
write.table(top_300_rs_id, file.path(output_file_path, "top_300_rs_id.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# Process LD data
num <- which(LD_matrix[, 1] == target_rs_id)
LD <- LD_matrix[num, ]

# Update LD score in GWAS data based on the LD matrix
eQTL_data$LDscore <- 0
for (i in 2:ncol(LD)) {
  num <- which(eQTL_data$rs_id == colnames(LD)[i])
  if (length(num) > 0) {
    eQTL_data$LDscore[num] <- LD[1, i]
  }
}

# Set colors based on LD score
color_palette <- c("#EE4000", "#EEAD0E", "#33cc00", "#00ccff", "#6959CD")
for (i in 1:nrow(eQTL_data)) {
  LD_index <- findInterval(eQTL_data$LDscore[i], c(0, 0.2, 0.4, 0.6, 0.8, 1))
  eQTL_data$color[i] <- color_palette[LD_index]
}

# Plot
p1 <- ggplot(eQTL_data, aes(x = position / 1000000, y = -log(pval_nominal_sqtl, base = 10), fill = color)) +
  geom_point(size = 4, shape = 21) +
  theme_classic() +
  scale_fill_manual(values = color_palette) +
  labs(x = NULL, y = "-log10(P-value)", title = "GWAS vs LD Score") +
  theme(panel.border = element_rect(color = "black", linewidth = 1, fill = NA))

# Save plot
ggsave(file.path(output_file_path, "gene_ld_gwas_plot.pdf"), p1, width = 5, height = 3)

