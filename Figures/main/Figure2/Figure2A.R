# Set the working directory
setwd("/home/project/11003054/e1101919/muscle_QTL/RNAseq/06.DEG")

# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(ggrepel)

# Load the data
data <- readRDS("all_info_for_DEGplot.rds")

# Clean and prepare the data
data <- data %>%
  mutate(
    chrom = factor(V1, levels = paste0("chr", c(1:22, "X", "Y", "M"))),
    chr = as.numeric(gsub("chr", "", chrom)),  # Convert chromosome names to numeric
    chr = factor(chr, levels = c(1:22, "X", "Y", "M")),  # Ensure chromosome order
    logFC = as.numeric(logFC),  # Convert logFC to numeric if not already
    direction = case_when(  # Determine the direction of change
      FDR > 0.05 ~ 'None',
      abs(logFC) < 0.3 ~ 'None',
      logFC >= 0.3 ~ 'Up',
      TRUE ~ 'Down'
    )
  ) %>%
  arrange(chrom, logFC)

# Plotting differential expression results
P1 <- ggplot(de_result, aes(x = chrom, y = logFC)) +
  geom_jitter(data = filter(de_result, direction == "None"), aes(color = chrom), width = 0.5, size = 1.5) +
  geom_jitter(data = filter(de_result, direction == "Up"), color = 'red', width = 0.5, size = 1.5) +
  geom_jitter(data = filter(de_result, direction == "Down"), color = 'blue', width = 0.5, size = 1.5) +
  scale_color_manual(values = c(rep(c("grey30", "grey70"), 12), "grey30")) +
  geom_hline(yintercept = c(0.3, -0.3), linetype = "dashed", color = c("red", "blue")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = paste0("chr", c(1:22, "X", "Y", "M"))) +
  labs(y = "logFC") +
  geom_text_repel(
    data = filter(de_result, direction == "Up"), 
    nudge_y = 0.5, 
    point.padding = NA, 
    aes(label = gene), 
    max.overlaps = 4
  ) +
  geom_text_repel(
    data = filter(de_result, direction == "Down"), 
    nudge_y = -0.5, 
    point.padding = NA, 
    aes(label = gene), 
    max.overlaps = 3.5
  )

# Display the plot
print(P1)

# Save the plot to a file
ggsave('Mirrow_Manhattan_v4.pdf', P1, width = 12, height = 6)
