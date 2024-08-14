# Load the ggplot2 package
library(ggplot2)

# Create the data
data <- data.frame(
  group = c("Pre", "Post"),
  pi1 = c(1, 0.91)
)

data$group <- factor(data$group, levels = c("Pre", "Post"))

P1 <- ggplot(data, aes(x = group, y = pi1, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  # Add borders and adjust bar width
  scale_fill_manual(values = c("Pre" = "skyblue", "Post" = "salmon")) +  # Customize bar colors
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = scales::percent_format(scale = 1)) +  # Adjust y-axis
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),  # Enhance axis titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold plot title
    legend.position = "none"  # Remove the legend
  ) +
  labs(
    x = "Group",  # Custom x-axis label
    y = "Pi1 value",  # Custom y-axis label
    title = "Comparison of Pi1 Values Between Groups"  # Custom and descriptive title
  )

print(P1)

ggsave("Pi1_comparison_plot.pdf", plot = P1, width = 5, height = 4)
