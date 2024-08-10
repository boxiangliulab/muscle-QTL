library(circlize)
library(dplyr)

# Assume Pre_eGene and Post_eGene dataframes are already loaded and have a 'family' column
# Extract shared and specific eQTL families
eQTL_shared <- intersect(Pre_eGene$family, Post_eGene$family)
Pre_eQTL_specific <- setdiff(Pre_eGene$family, Post_eGene$family)
Post_eQTL_specific <- setdiff(Post_eGene$family, Pre_eGene$family)

# Convert vectors to data frames and set column names
eQTL_shared <- data.frame(family = eQTL_shared)
Pre_eQTL_specific <- data.frame(family = Pre_eQTL_specific)
Post_eQTL_specific <- data.frame(family = Post_eQTL_specific)

# Assuming 'gencode' is preloaded and contains genomic annotations
gencode$V5 <- gsub("\\..*", "", gencode$V5)
colnames(gencode)[5] <- "family"

# Merge with gencode to get full annotation details
Pre_eQTL_full <- merge(Pre_eQTL_specific, gencode, by = "family")
Post_eQTL_full <- merge(Post_eQTL_specific, gencode, by = "family")
eQTL_shared_full <- merge(data.frame(family = eQTL_shared), gencode, by = "family")

# Calculate average genomic position for plotting (assuming V2 and V3 are positions)
Pre_eQTL_full$average <- round(rowMeans(Pre_eQTL_full[, c("V2", "V3")], na.rm = TRUE))
Post_eQTL_full$average <- round(rowMeans(Post_eQTL_full[, c("V2", "V3")], na.rm = TRUE))
eQTL_shared_full$average <- round(rowMeans(eQTL_shared_full[, c("V2", "V3")], na.rm = TRUE))

# Combine all for plotting
eQTL_data <- bind_rows(
  Pre_eQTL_full %>% mutate(group = "Pre-specific"),
  Post_eQTL_full %>% mutate(group = "Post-specific"),
  eQTL_shared_full %>% mutate(group = "Shared")
)

# Setting the chromosomal data correctly assuming chromosome data is in a column `chr`
eQTL_data$chr <- factor(eQTL_data$chr, levels = paste0("chr", 1:22))

# Initialize circos plot
circos.clear()
circos.par("track.height" = 0.1, "gap.degree" = 3, "start.degree" = 90, "clock.wise" = TRUE, "track.margin" = c(0, 0.02))

# Initialize with human ideogram (hg38), assuming it has been properly set up
circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0("chr", 1:22))

# Define colors for each group
colors <- c("Pre-specific" = "#619CFF", "Post-specific" = "#F87660", "Shared" = "black")

# Create tracks for each group
for(group in unique(eQTL_data$group)) {
  chr_data <- subset(eQTL_data, group == group)
  circos.trackPlotRegion(factors = chr_data$chr, ylim = c(0, 1), track.height = 0.05, panel.fun = function(x, y) {
    sector.index <- CELL_META$sector.index
    chr_data_sub <- chr_data[chr_data$chr == sector.index, ]
    
    for(pos in chr_data_sub$average) {
      circos.lines(c(pos, pos), c(0, 1), col = colors[group])
    }
  })
}

# Print or save the plot
# circos.save("eQTL_circos_plot.png")
