## Load necessary library
library(ggplot2)
library(qqman)

## Figure S11A - Line chart similar to the previous Figure S9A
P1 <- ggplot(data, aes(x = Number_of_PC)) +
  geom_line(aes(y = PRE, color = "PRE"), size = 0.5) +
  geom_point(aes(y = PRE), color = "darkred", size = 4) +
  geom_line(aes(y = POST, color = "POST"), size = 0.5) +
  geom_point(aes(y = POST), color = "darkblue", size = 4) +
  labs(x = "Number of hidden factors", y = "Number of Genes", title = "Analysis of Gene Number by Hidden Factors") +
  scale_color_manual(values = c("PRE" = "darkred", "POST" = "darkblue")) +
  scale_x_continuous(breaks = 1:10) +
  geom_vline(xintercept = 1:2, linetype = "dashed", color = "black", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), legend.title = element_blank(), legend.text = element_text(size = 12),
        legend.position = "top", panel.border = element_rect(color = "black", fill = NA, size = 1))
ggsave("line_chart_analysis.pdf", plot = P1, width = 7, height = 5)

## Figure S11B - Histogram of p-value distribution
eqtl_pval_check <- ggplot(joint_filtered, aes(x = pval)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "grey") +
  labs(x = "p-value", y = "Frequency") +
  scale_x_continuous(limits = c(0, 1)) +  
  ylim(0, 15000) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
ggsave("pval_distribution_histogram.pdf", plot = eqtl_pval_check, width = 8, height = 6)

## Figure S11C - QQ plot
qqplot_data <- qqunif(joint_filtered$pval, plot.it = FALSE)
x <- as.data.frame(qqplot_data)
p <- ggplot(x, aes(x = x, y = y)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  xlab('Expected p-values') + 
  ylab('Observed p-values')
ggsave("qq_plot_pvalues.pdf", plot = p, width = 8, height = 6)

## Figure S11D - Density plot for distance to splicing donor/acceptor
plot_eqtl_vs_dist <- function(ras_filtered) {
  sig_eqtl <- ras_filtered[ras_filtered$pval < 1e-4, ]
  ggplot(sig_eqtl, aes(x = dist_to_splicing)) +
    geom_density(color = 'black', fill = 'black') +
    xlab('Distance to Splicing Donor/Acceptor') +
    ylab('Frequency')
}

## Figure S11E - Plot for fraction distance within intron
plot_fraction_distance <- function(ras_filtered) {
  ras_subset <- ras_filtered[abs(ras_filtered$intron_fraction_dist) < 1, ]
  p <- ggplot(ras_subset, aes(x = intron_fraction_dist, y = -log10(pval))) +
    geom_point(alpha = 0.5) +
    xlab('Fraction Distance Within Intron') +
    ylab(expression(-log[10]*'(p-value)'))
  return(p)
}
