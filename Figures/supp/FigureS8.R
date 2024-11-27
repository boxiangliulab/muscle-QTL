library(tidyverse)
library(ggplot2)

SAMS2_count_TPM <- read.csv("SAMS2_count_TPM.csv", row.names = 1)

SAMS2_long <- SAMS2_count_TPM %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "expression") %>%
  mutate(
    condition = ifelse(grepl("pre", sample), "Pre", "Post"),
    sample = gsub("_pre|_post", "", sample),
    gene = rownames(SAMS2_count_TPM)
  ) %>%
  select(sample, gene, condition, expression)

genes_of_interest <- c("MT-ND2", "MT-ND5", "MT-ND6", "MT-ATP8")
data_filtered <- SAMS2_long %>%
  filter(gene %in% genes_of_interest)

# Convert expression to log2(TPM+1) if not already transformed
data_filtered <- data_filtered %>%
  mutate(expression = log2(expression + 1))

p <- ggplot(data_filtered, aes(x = condition, y = expression, group = sample)) +
  geom_line(aes(group = sample), alpha = 0.3) +  # Draw lines connecting paired samples
  geom_point(size = 1.5, shape = 21, fill = "white") +
  facet_wrap(~ gene, scales = "free_y") +  # Separate plot for each gene
  labs(x = NULL, y = "log2(TPM + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text angle

# Adding p-values
# You can manually add them using annotate() if you have a predefined list
p <- p + annotate("text", x = 1.5, y = Inf, label = "Paired t-test p = X.XXe-XX", vjust = 1.5)

# Print the plot
print(p)
