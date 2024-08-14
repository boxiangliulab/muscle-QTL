library(ggplot2)
library(readr)

diff_quant_data_file <- paste(out, "DiffQuantData.txt", sep = "")
volcano_plot_file <- paste(out, "volcano.pdf", sep = "")

data <- read_delim(diff_quant_data_file, delim = "\t")

cols <- c("No_Sig" = "#b4b4d8", "Up" = "#e94234", "Down" = "#269846")

# Plot
g <- ggplot(data, aes(x = log2(FC), y = -log10(pvalue))) +
  geom_point(aes(color = Sig), size = 4, alpha = 0.8) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = 6, size = 0.6, color = "black") +
  geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linetype = 6, size = 0.6, color = "black") +
  xlab("log2(FoldChange)") +
  ylab("-log10(pvalue)") +
  theme_classic()

# Save plot
pdf(volcano_plot_file, width = 6, height = 5)
print(g)
dev.off()
