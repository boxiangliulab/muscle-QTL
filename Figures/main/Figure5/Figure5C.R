library(ggplot2)

data <- data.frame(
  Population = c("SAMS2_EAS", "1KG_EAS", "1KG_SAS", "1KG_AMR", "1KG_EUR", "1KG_AFR"),
  AlleleFrequency = c(0.23, 0.22, 0.1, 0.03, 0.02, 0.01),
  stringsAsFactors = FALSE
)

# Adding a reference line value
threshold <- 0.05

# Create the plot
ggplot(data, aes(x = Population, y = AlleleFrequency, fill = Population)) +
  geom_bar(stat = "identity", color = "black") + 
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  scale_fill_grey(start = 0.8, end = 0.2) + 
  labs(title = "Alt Allele Frequency of rs129437434",
       x = NULL, # Removing the x-axis label
       y = "Alt Allele Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 
