library(VIM)
library(mice)
library(dplyr)
library(ggplot2)

# Function to perform imputation and plot the data
plot_levels_over_time <- function(data, variable_prefix, time_order, output_filename) {
  # Perform mice imputation
  imp <- mice(data, m = 5, maxit = 50, meth = 'pmm', seed = 123)
  completed_data <- complete(imp)
  
  # Transform data to long format
  df_long <- completed_data %>%
    pivot_longer(cols = starts_with(variable_prefix), names_to = "Time", values_to = "Level")
  
  # Set the factor levels for Time based on provided order
  df_long$Time <- factor(df_long$Time, levels = time_order)
  
  # Calculate summary statistics
  summary_data <- df_long %>%
    group_by(Condition, Time) %>%
    summarise(mean_level = mean(Level, na.rm = TRUE),
              std_err = sd(Level, na.rm = TRUE) / sqrt(n()),
              .groups = 'drop')
  
  # Define shape mapping
  shape_mapping <- c("Post" = 17, "Pre" = 16)
  
  # Create the plot
  p <- ggplot(summary_data, aes(x = Time, y = mean_level, color = Condition, group = Condition, shape = Condition)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_level, ymax = mean_level + std_err), width = 0.2, color = "gray50", size = 0.5) +
    geom_point(size = 4) +
    labs(x = "Time", y = paste("Mean", variable_prefix, sep = " "), title = paste("Mean", variable_prefix, "Levels Over Time by Condition")) +
    scale_color_manual(values = c("black", "black")) +
    scale_shape_manual(values = shape_mapping) +
    scale_x_discrete(labels = time_order) +
    theme_minimal() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.ticks = element_line(size = 0.5), axis.line = element_line(size = 0.5))
  
  # Display the plot
  print(p)
  
  # Save the plot to a file
  ggsave(output_filename, plot = p, width = 6, height = 4, units = "in")
}

# Variables of interest and their respective time orders
time_order_insulin <- c("insulin.0h","insulin.1.5h","insulin.2h")
time_order_glucose <- c("glucose.0h", "glucose.0.5h", "glucose.1h", "glucose.1.5h", "glucose.2h", "glucose.4h")

# Assuming 'insulin_time' contains both insulin and glucose data or adjust accordingly
plot_levels_over_time(insulin_time, "insulin", time_order_insulin, "Insulin level over time.pdf")
plot_levels_over_time(insulin_time, "glucose", time_order_glucose, "Glucose level over time.pdf")

