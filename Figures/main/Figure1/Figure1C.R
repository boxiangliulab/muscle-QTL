library(dplyr)
library(tidyr)
library(ggplot2)

# Convert selected variables to numeric type because data might be in factor or character form.
# Variables_of_interest here means the fat and lean-related variables.
pre_selected <- pre_imputed[variables_of_interest] %>% mutate(across(everything(), ~as.numeric(as.character(.))))
post_selected <- post_imputed[variables_of_interest] %>% mutate(across(everything(), ~as.numeric(as.character(.))))

# Add a time column to distinguish between 'Pre' and 'Post' data
pre_selected$Time <- 'Pre'
post_selected$Time <- 'Post'

# Combine the data for easier manipulation and plotting
data <- bind_rows(pre_selected, post_selected)

# Convert the wide format data to long format for easier plotting
data_long <- pivot_longer(data, cols = variables_of_interest, names_to = "Variable", values_to = "Value")

# Calculate the percentage change for each variable
change_percent <- (post_selected - pre_selected) / pre_selected * 100
change_percent_long <- pivot_longer(change_percent, cols = variables_of_interest, names_to = "Variable", values_to = "PercentChange")

# Plotting the percentage change as a boxplot
ggplot(change_percent_long, aes(x = Variable, y = PercentChange, fill = Variable)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Percentage Change of Fat-Related Variables", x = "Variable", y = "Percent Change (%)")

# Dot plot for individual percentage changes
ggplot(change_percent_long, aes(x = Variable, y = PercentChange, color = Variable)) +
  geom_point(alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Dot Plot of Percentage Change for Each Subject", x = "Variable", y = "Percent Change (%)")

# Calculate the median of percentage changes by variable and sort by median
median_changes <- change_percent_long %>%
  group_by(Variable) %>%
  summarize(Median = median(PercentChange, na.rm = TRUE)) %>%
  arrange(Median)

# Print the sorted median changes
print(median_changes)

# Update factor levels to sort the variables in the plots based on the median change
change_percent_long$Variable <- factor(change_percent_long$Variable, levels = median_changes$Variable)

# Plotting a sorted boxplot of percentage changes
ggplot(change_percent_long, aes(x = Variable, y = PercentChange, fill = Variable)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-100, 100)) +  # Limiting y-axis to -100% to 100% for better visualization
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Sorted Boxplot of Percentage Change of Fat-Related Variables", x = "Variable", y = "Percent Change (%)")

