# Load necessary libraries
library(ggplot2)

# Load the lipidomics data
load("~/OneDrive - National University of Singapore/SAMS2_bulk/github/lipidomics_data.RData")

# Run PCA on the lipidomics data excluding non-numeric columns and the 'Time' column
# Typically, you'd select only the relevant lipid measurements columns, here assumed all except the first
pca_results <- prcomp(lipidomics_data[, -which(names(lipidomics_data) %in% c("Time"))], scale. = TRUE)

# Extract scores
scores <- data.frame(pca_results$x)

# Add the 'Time' column back for coloring
scores$Time <- lipidomics_data$Time

# Convert 'Time' into a factor and set levels and colors
scores$Time <- factor(scores$Time, levels = c("Pre", "Post"))
colors <- c("Pre" = "green", "Post" = "red")

# Create the PCA plot with PC2 on the x-axis and PC1 on the y-axis
p <- ggplot(scores, aes(x = PC2, y = PC1, color = Time)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust point size and transparency
  scale_color_manual(values = colors) +
  labs(x = "PC2", y = "PC1", title = "PCA of Lipids Data") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove legend title

# Print the plot
print(p)

# Optionally, save the plot
ggsave("PCA_Lipids_Plot.pdf", plot = p, width = 8, height = 6, units = "in")

