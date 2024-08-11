# Load necessary libraries
library(dplyr)
library(ggplot2)

# Function to calculate fat and lean percentages
calculate_percentages <- function(df) {
  df %>%
    mutate(
      rarm_fat_perc = rarm_fat_post / rarm_totalmass_post * 100,
      rarm_lean_perc = rarm_lean_post / rarm_totalmass_post * 100,
      larm_fat_perc = larm_fat_post / larm_totalmass_post * 100,
      larm_lean_perc = larm_lean_post / larm_totalmass_post * 100,
      trunk_fat_perc = trunk_fat_post / trunk_totalmass_post * 100,
      trunk_lean_perc = trunk_lean_post / trunk_totalmass_post * 100,
      rleg_fat_perc = rleg_fat_post / rleg_totalmass_post * 100,
      rleg_lean_perc = rleg_lean_post / rleg_totalmass_post * 100,
      lleg_fat_perc = lleg_fat_post / lleg_totalmass_post * 100,
      lleg_lean_perc = lleg_lean_post / lleg_totalmass_post * 100,
      twb_fat_perc = twb_fat_post / twb_totalmass_post * 100,
      twb_lean_perc = twb_lean_post / twb_totalmass_post * 100,
      st_fat_perc = st_fat_post / st_totalmass_post * 100,
      st_lean_perc = st_lean_post / st_totalmass_post * 100,
      head_fat_perc = head_fat_post / head_totalmass_post * 100,
      head_lean_perc = head_lean_post / head_totalmass_post * 100
    )
}

# Apply the function to calculate percentages for pre and post datasets
pre_perc <- calculate_percentages(pre)
post_perc <- calculate_percentages(post)

# Function to calculate percentage changes
calculate_change <- function(pre, post, column) {
  (post[[column]] - pre[[column]]) / pre[[column]] * 100
}

# Function to remove outliers
remove_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.25)
  Q3 <- quantile(df[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  df %>% filter(df[[column]] >= lower_bound & df[[column]] <= upper_bound)
}

# List of metrics to be analyzed
metrics <- c(
  "liver.fat", "vat", "sat", "dsat", "ssat", "skeletal.muscle.fat",
  "rarm_fat_perc", "rarm_lean_perc", "larm_fat_perc", "larm_lean_perc",
  "trunk_fat_perc", "trunk_lean_perc", "rleg_fat_perc", "rleg_lean_perc",
  "lleg_fat_perc", "lleg_lean_perc", "twb_fat_perc", "twb_lean_perc",
  "st_fat_perc", "st_lean_perc", "head_fat_perc", "head_lean_perc"
)

# Calculate mean changes and standard errors for each metric after removing outliers
changes_stats <- t(sapply(metrics, function(metric) {
  pre_clean <- remove_outliers(pre_perc, metric)
  post_clean <- remove_outliers(post_perc, metric)
  mean_change <- calculate_change(pre_clean, post_clean, metric)
  se_change <- sd(post_clean[[metric]]) / sqrt(nrow(post_clean))
  c(mean_change, se_change)
}))

# Convert changes stats to a data frame
change_df <- data.frame(
  Metric = metrics,
  Mean_Change = changes_stats[,1],
  SE_Change = changes_stats[,2]
)

# Calculate 95% confidence intervals
confidence_level <- 0.95
t_value <- qt(1 - (1 - confidence_level) / 2, df = nrow(pre) - 1)
change_df <- change_df %>%
  mutate(
    CI_lower = Mean_Change - t_value * SE_Change,
    CI_upper = Mean_Change + t_value * SE_Change,
    Significant = ifelse(CI_lower > 0 | CI_upper < 0, "Significant", "Not Significant")
  )

# Sort by mean change percentage
change_df <- change_df %>% arrange(Mean_Change)

# Create a forest plot of body composition changes
ggplot(change_df, aes(x = Mean_Change, y = reorder(Metric, Mean_Change))) +
  geom_point() +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Mean Change (%)", y = "", title = "Forest Plot of Body Composition Changes") +
  theme_minimal()
