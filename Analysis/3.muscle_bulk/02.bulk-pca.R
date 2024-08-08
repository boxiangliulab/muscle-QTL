# Clinical Data - 02. pre-process
# Load necessary libraries
library(ggplot2)

# Load TPM data
load("SAMS2_muscle_TPM.RData")  # This loads the TPM data into the variable 'tpm'

# Perform PCA on the TPM data
pca_result <- prcomp(t(tpm), scale. = TRUE)
scores <- as.data.frame(pca_result$x)
scores$group <- ifelse(grepl("Pre", rownames(scores)), "Pre", "Post")

# Summarize PCA results
pca_summary <- summary(pca_result)
print(pca_summary$importance[, 1:2])  # Print the importance of the first two principal components

# Plot PCA results
ggplot(scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +  # Draw points
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.5) +  # Add confidence ellipse
  theme_minimal() +
  ggtitle("PCA of Pre vs Post Groups with 95% CI") +
  xlab(paste("Principal Component 1 (", round(pca_summary$importance[2, 1] * 100, 1), "%)", sep = "")) +
  ylab(paste("Principal Component 2 (", round(pca_summary$importance[2, 2] * 100, 1), "%)", sep = "")) +
  scale_color_manual(values = c("Pre" = "#619CFF", "Post" = "#F87660"))
