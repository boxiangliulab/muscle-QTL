library(ggplot2)
library(dplyr)

# Extract 'myocyte' column from each dataset and create a combined data frame with an intervention label
pre_data <- pre_celltype_fraction_bmind %>%
  select(myocyte) %>%
  mutate(Intervention = 'Pre')

post_data <- post_celltype_fraction_bmind %>%
  select(myocyte) %>%
  mutate(Intervention = 'Post')

# Combine the data frames
combined_data <- bind_rows(pre_data, post_data)

# Generate the violin plot
ggplot(combined_data, aes(x = Intervention, y = myocyte, fill = Intervention)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("Pre" = "red", "Post" = "blue")) +
  labs(title = "Violin Plot of Myocyte Fraction in Skeletal Muscle RNA-seq",
       subtitle = "Pre- and Post-intervention",
       x = NULL,
       y = "Myocyte Fraction") +
  theme_minimal() +
  theme(legend.position = "none") 

