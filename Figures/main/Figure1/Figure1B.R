# 'change' dataframe is from Analysis 1. clinical data
matrix_data <- as.data.frame(change[, 4:26])

matrix_data$ID = rownames(matrix_data)

library(reshape2)
melted_data <- melt(matrix_data, id.vars = "ID")

library(ggplot2)
p <- ggplot(melted_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_text(aes(label = round(value, 2)), position = position_dodge(width = 0.75), vjust = -0.5, check_overlap = TRUE) + p
  labs(x = "Trait", y = "Post-Pre Difference", title = "Boxplot of Trait Differences") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Print the plot
print(p)

# Save the plot to a file
ggsave("trait_boxplot.pdf", p, width = 10, height = 6, dpi = 300)

