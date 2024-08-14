library(ggplot2)

# Part A: Histogram
# Load exercise session data
exercise_data <- read.csv("SAMS2_exercise_data.csv")  # Extract the gym_visits variable

# Create histogram
p1 <- ggplot(exercise_data, aes(x = session_count)) +
  geom_histogram(binwidth = 5, fill = "gray", color = "black") +
  geom_vline(xintercept = 48, linetype = "dashed", color = "red") +  # Add dashed line at x = 48
  labs(title = "The number of attended structured exercise sessions",
       x = "Number of Sessions",
       y = "Frequency") +
  theme_minimal()

# Part B: Error Bar Plot
# Load diet change data
diet_data <- read.csv("SAMS2_diet_data.csv")  # Extract the diet data 

# Create error bar plot
p2 <- ggplot(diet_data, aes(x = component, y = percentage_change)) +
  geom_point() +
  geom_errorbar(aes(ymin = percentage_change - 10, ymax = percentage_change + 10), width = 0.2) +  # Adjust error range
  coord_flip() +
  labs(title = "Dietary Intake Percentage Change (%)",
       x = NULL,
       y = "Percentage Change (%)") +
  theme_minimal()

# Print the plots
print(p1)
print(p2)

