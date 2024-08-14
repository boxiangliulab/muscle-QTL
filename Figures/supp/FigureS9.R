## Figure S9A
# Load necessary library
library(ggplot2)

# Prepare data
P1 <- ggplot(data, aes(x = Number_of_PC)) +
  geom_line(aes(y = PRE, color = "PRE"), size = 0.5) +
  geom_point(aes(y = PRE), color = "darkred", size = 4) +
  geom_line(aes(y = POST, color = "POST"), size = 0.5) +
  geom_point(aes(y = POST), color = "darkblue", size = 4) +
  labs(
    x = "Number of hidden factors",
    y = "Number of Genes",
    title = "Analysis of Gene Number by Hidden Factors"
  ) +
  scale_color_manual(values = c("PRE" = "darkred", "POST" = "darkblue")) +
  scale_x_continuous(breaks = 1:10) +
  geom_vline(xintercept = 1:2, linetype = "dashed", color = "black", size = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# Display the plot
print(P1)

# Save the plot to a file
ggsave("sIntron_PC_correction.pdf", P1, width = 7, height = 5)

## Figure S9B
# Histogram for 'joint_filtered' dataset
eqtl_pval_check <- ggplot(joint_filtered, aes(x = pval)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "grey") +
  labs(x = "p-value", y = "Frequency") +
  scale_x_continuous(limits = c(0, 1)) +
  ylim(0, 15000) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggsave("../../../muscleQTL/eQTL mapping/eqtl_pval_check.pdf", plot = eqtl_pval_check, width = 8, height = 6)

# Histogram for 'test_for_pval_check' dataset
test_pval_check_plot <- ggplot(test_for_pval_check, aes(x = pval)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "grey") +
  labs(x = "p-value", y = "Frequency") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# This line saves the plot to a specified directory; adjust the path as needed.
ggsave("../eQTL mapping/eqtl_test_pval_check.pdf", plot = test_pval_check_plot, width = 8, height = 6)

## Figure S9C
# QQ plot for 'post_pval_check' dataset
qqplot <- qq(post_pval_check$pval)  # Assuming 'post_pval_check' is correctly formatted

# Save the QQ plot
ggsave("../eQTL mapping/eqtl_qq_check.pdf", plot = qqplot, width = 8, height = 6)

# Additional QQ plots for demonstration (ensure these are valid calls)
qq(Post_rasqual_pval_check$p_value)  # Assumes 'Post_rasqual_pval_check' is a defined data frame
qq(Post_full_variants.txt$V19)        # Assumes 'Post_full_variants.txt' is read and V19 is the column of interest

## Figure S9D 
# Function to plot density of significant eQTLs by distance
plot_eqtl_vs_dist <- function(ras_filtered) {
  sig_eqtl <- ras_filtered[ras_filtered$pval < 1e-4, ]
  
  ggplot(sig_eqtl, aes(x = dist)) +
    geom_density(color = 'black', fill = 'black') +
    xlab('Distance to TSS') +
    ylab('Frequency')
}

## Figure S9E
plot_pvalue_vs_position <- function(ras_filtered) {
  ras_subset <- ras_filtered[abs(ras_filtered$dist) < 2e5, ]
  idx <- sample(nrow(ras_subset), size = min(5e3, nrow(ras_subset)))
  ras_subset <- ras_subset[idx, ]
  
  ggplot(ras_subset, aes(x = dist, y = -log10(pval), size = (10 * abs(pi - 0.5))^2)) +
    geom_point(alpha = 0.2) +
    scale_size_continuous(name = expression('|'*pi*'-0.5|'), 
                           breaks = (10 * c(0, 0.1, 0.2, 0.3, 0.4, 0.5))^2, 
                           labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    scale_y_continuous(trans = 'sqrt') +
    xlab('Distance to TSS') +
    ylab(expression(-log[10]*'(p-value)')) +
    theme(legend.position = c(0.05, 0.99), legend.justification = c('left', 'top'))
}
