library(dplyr)
library(tidyr)
library(ggplot2)

# Define the count data
sig_coloc_count <- data.frame(
  Trait = c("T2D_sakaue", "T2D_spracklen", "BMI", "Body Weight", 
            "Coronary Artery Disease", "Glucose", "HbA1c", "LDL", 
            "Total Cholesterol", "Triglycerides", "HDL", "Height"),
  A = c(0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 8),
  B = c(0, 3, 1, 2, 1, 0, 0, 0, 1, 1, 0, 6),
  C = c(1, 2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 5),
  D = c(0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 3)
)

# Define the proportion data
coloc_proportion <- data.frame(
  Trait = c("T2D_sakaue", "T2D_spracklen", "BMI", "Body Weight", 
            "Coronary Artery Disease", "Glucose", "HbA1c", "LDL", 
            "Total Cholesterol", "Triglycerides", "HDL", "Height"),
  A = c(0, 0.02, 0.05, 0.03, 0, 0, 0.1, 0, 0, 0, 0, 0.05),
  B = c(0, 0.06, 0.04, 0.05, 0.07, 0, 0, 0, 0.07, 0.1, 0, 0.04),
  C = c(0.03, 0.03, 0, 0.03, 0.08, 0, 0, 0, 0, 0, 0, 0.02),
  D = c(0, 0.01, 0, 0.02, 0.04, 0, 0, 0, 0, 0, 0, 0.01)
)

# Merge count and proportion data by Trait
data <- left_join(sig_coloc_count, coloc_proportion, by = "Trait", suffix = c("_count", "_proportion"))

# Reshape data to long format for plotting
data_long <- pivot_longer(data, cols = -Trait, names_to = "variable", values_to = "value") %>%
  mutate(
    group = sub(".*_", "", variable),  # Extract group (count/proportion)
    variable = sub("_.*", "", variable)  # Extract variable name
  ) %>%
  pivot_wider(names_from = group, values_from = value)  # Convert to wider format

# Generate the plot
ggplot(data_long, aes(x = variable, y = Trait, size = proportion, color = count)) +
  geom_point() +
  scale_color_gradient(low = "blue3", high = "brown3", name = "Sig Coloc Count") +
  labs(size = "Coloc Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
