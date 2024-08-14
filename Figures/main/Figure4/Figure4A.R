# Load necessary libraries
library(circlize)

shared_df <- data.frame(gene = shared_genes, group = "Shared")
pre_specific_df <- data.frame(gene = pre_specific_genes, group = "Pre-specific")
post_specific_df <- data.frame(gene = post_specific_genes, group = "Post-specific")

# Combine data frames
gene_data <- rbind(shared_df, pre_specific_df, post_specific_df)

# Merge with additional details from original data frames
gene_details <- rbind(Pre_sGene, Post_sGene)
colnames(gene_details)[1] <- "gene"  # Ensure the gene column name matches
plot_data <- merge(gene_data, gene_details, by = "gene", all.x = TRUE)

# Calculate the average position for visualization
plot_data$average <- rowMeans(plot_data[, c("V3", "V4")], na.rm = TRUE)

# Prepare data for circlize plotting
plot_data$chr <- factor(plot_data$chr, levels = paste0("chr", 1:22))

# Initialize circos plot
circos.clear()
circos.par("track.height" = 0.1, gap.degree = 3, start.degree = 90,
           clock.wise = TRUE, track.margin = c(0, 0.02), cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0('chr', 1:22))

# Define colors for groups
group_colors <- c("Shared" = "black", "Pre-specific" = "#F87660", "Post-specific" = "#619CFF")

# Add tracks for each group
for(group_name in names(group_colors)) {
  group_data <- subset(plot_data, group == group_name)
  circos.trackPlotRegion(factors = group_data$chr, ylim = c(0, 1), track.height = 0.05, 
                         panel.fun = function(x, y) {
                           sector.index <- CELL_META$sector.index
                           chr_data <- subset(group_data, chr == sector.index)
                           for (pos in chr_data$average) {
                             circos.lines(c(pos, pos), c(0, 1), col = group_colors[group_name])
                           }
                         })
}

# Saving the circos plot
circos.save("sQTL_circular_plot.png")

