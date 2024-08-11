library(dplyr)
library(tidyr)
library(ggplot2)

# Define the categories for the different measurements
categories_df <- data.frame(
  Trait = c("d14.waist.circ", "d14.hip.circ", "d14.waist.hip", "arm1_pre", "d14.arm.s.span", 
            "wd1l_pre", "wd1r_pre", "kd1l_pre", "kd1r_pre",
            "subs1_pre", "suprailiac_pre", "paraumb_pre", "skinsum_pre", "si1_pre", "pu1_pre", 
            "d14.biceps", "d14.triceps",
            "l1_l4_area_post", "l1_area_post", "l2_area_post", "l3_area_post", "l4_area_post", 
            "fneck_area_post", "lthip_area_post", "ud_area_post", "mid_area_post", "ut_area_post", 
            "r_u_area_post", "larm_area_post", "rarm_area_post", "lleg_area_post", "rleg_area_post", 
            "st_area_post", "head_area_post", "twb_area_post", "lr_area_post", "rr_area_post", 
            "ts_area_post", "ls_area_post", "pelvis_area_post"),
  Category = c(rep("BaseMeasurements", 9), rep("SkinfoldThickness", 8), rep("DEXASkeleton", 23))
)

# Convert 'Trait' column to character to ensure correct merging
categories_df$Trait <- as.character(categories_df$Trait)

# Example data frame for percent changes (this should be your actual data)
change_percent_long <- data.frame(
  Trait = sample(categories_df$Trait, 100, replace = TRUE),  # example Trait names
  PercentChange = rnorm(100)  # example PercentChanges
)

# Merge percent change data with categories
change_percent_long <- left_join(change_percent_long, categories_df, by = "Trait")

# Define custom color scheme for the plot
colors <- c("BaseMeasurements" = "#ccebc5", "SkinfoldThickness" = "#fbb4ae", "DEXASkeleton" = "#b3cde3")

# Create a boxplot displaying percent changes by category with customized aesthetics
ggplot(change_percent_long, aes(x = PercentChange, y = fct_reorder(Trait, PercentChange, .fun = median), fill = Category)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  scale_fill_manual(values = colors) +  # Apply custom colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Percentage Change in Body Measurement Traits by Category",
       x = "Percent Change (%)",
       y = "Trait") +
  theme(legend.position = "bottom")  # Display legend at the bottom

# Calculate and plot differences using detailed data handling (Example placeholders)
# Variables of interest from the hypothetical data
variables_of_interest <- names(categories_df$Trait)

# Assuming you have pre and post data frames named 'final_pre_data_clinical' and 'final_post_data_clinical'
pre_selected <- final_pre_data_clinical[variables_of_interest] %>% mutate_all(as.numeric)
post_selected <- final_post_data_clinical[variables_of_interest] %>% mutate_all(as.numeric)

# Compute percent changes
change_percent <- 100 * (post_selected - pre_selected) / pre_selected

# Convert to long format for easier plotting
change_percent_long <- pivot_longer(change_percent, cols = everything(), names_to = "Trait", values_to = "PercentChange")

# Plot percent changes with a forest plot showing confidence intervals (CIs)
# Note: This example assumes CIs are calculated and included in your data.
ggplot(change_percent_long, aes(x = PercentChange, y = fct_reorder(Trait, PercentChange, .fun = median), fill = Trait)) +
  geom_point() +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +  # Example CIs
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Percentage Difference and 95% CI Between Post and Pre Measurements",
       x = "Percentage Difference (%)", y = "Trait") +
  theme_minimal()
