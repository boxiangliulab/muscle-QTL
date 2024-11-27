library(ggplot2)
library(dplyr)
library(tidyr)

merged_data_with_baseline$recessive <- ifelse(merged_data_with_baseline$Code_numeric == 1, 1, 0)
merged_data_with_baseline$dominant <- ifelse(merged_data_with_baseline$Code_numeric == 0, 1, 0)
merged_data_with_baseline$additive <- ifelse(merged_data_with_baseline$Code_numeric == 2, 1, 0)
merged_data_with_baseline$additive <- as.numeric(as.character(merged_data_with_baseline$additive))

# additive
additive_model <- lm(weight_percent_change ~ Code_numeric + age + no.of.gym, 
                     data = merged_data_with_baseline)
summary(additive_model)

# recessive
recessive_model <- lm(weight_percent_change ~ recessive + age + no.of.gym, 
                      data = merged_data_with_baseline)
summary(recessive_model)

# dominant
dominant_model <- lm(weight_percent_change ~ dominant + age + no.of.gym, 
                     data = merged_data_with_baseline)
summary(dominant_model)

AIC_comparison <- data.frame(
  Model = c("Additive", "Dominant", "Recessive"),
  AIC = c(AIC(additive_model), AIC(dominant_model), AIC(recessive_model)),
  R_squared = c(summary(additive_model)$r.squared,
                summary(dominant_model)$r.squared,
                summary(recessive_model)$r.squared)
)
print(AIC_comparison)

detailed_stats <- merged_data_with_baseline %>%
  group_by(Code) %>%
  summarise(
    n = n(),
    mean_change = mean(weight_percent_change, na.rm = TRUE),
    sd_change = sd(weight_percent_change, na.rm = TRUE),
    se_change = sd_change / sqrt(n),
    ci_lower = mean_change - 1.96 * se_change,
    ci_upper = mean_change + 1.96 * se_change
  )
print(detailed_stats)

ggplot(detailed_stats, aes(x = Code, y = mean_change)) +
  geom_bar(stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(x = "Genotype", y = "Mean Weight Change (%)",
       title = "Mean Weight Change by Genotype",
       subtitle = "Error bars represent 95% CI") +
  theme_minimal()

# 4.2 violin plot with boxplot
ggplot(merged_data_with_baseline, aes(x = Code, y = weight_percent_change)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7) +
  labs(x = "Genotype", y = "Weight Change (%)",
       title = "Distribution of Weight Change by Genotype") +
  theme_minimal()

