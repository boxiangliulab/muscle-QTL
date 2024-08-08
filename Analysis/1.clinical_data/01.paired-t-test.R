# Clinical Data - 01. paired t-test
# Load necessary packages
library(dplyr)
library(tidyr)
library(broom)
library(stats)

# Ensure that the rows of Pre and Post correspond 
pre_data <- pre_imputed # backup
post_data <- post_imputed

# Select the clinical data rows
pre_data <- pre_data[,6:340]
post_data <- post_data[,6:340]

# Calculate the mean and standard deviation for each variable
pre_stats <- data.frame(
  mean = apply(pre_data, 2, mean),
  sd = apply(pre_data, 2, sd)
)

post_stats <- data.frame(
  mean = apply(post_data, 2, mean),
  sd = apply(post_data, 2, sd)
)

# Paired t-test
p_values <- mapply(function(x, y) t.test(x, y, paired = TRUE)$p.value,
                   pre_data, post_data)


# Adjust the p-value with BH
q_values <- p.adjust(p_values, method = "BH")

# Organize the outputs
results <- data.frame(
  variable = colnames(pre_data),
  pre_mean = pre_stats$mean,
  pre_sd = pre_stats$sd,
  post_mean = post_stats$mean,
  post_sd = post_stats$sd,
  q_value = q_values,
  significant = q_values < 0.05
)

# Keep two decimal points and organize the format
results$pre_mean <- round(results$pre_mean, 2)
results$pre_sd <- round(results$pre_sd, 2)
results$post_mean <- round(results$post_mean, 2)
results$post_sd <- round(results$post_sd, 2)
results$q_value <- q_values

results <- results %>%
  mutate(Pre = paste0(results$pre_mean, " ± ", results$pre_sd))

results <- results %>%
  mutate(Post = paste0(results$post_mean, " ± ", results$post_sd))

results <- as.data.frame(results)

results$Pre <- round(results$Pre, 2)
results$Post <- round(results$Post, 2)
