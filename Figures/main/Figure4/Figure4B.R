library(ggplot2)
library(dplyr)
library(stringr)

# Extracting variant IDs from SNP information
pre_present$variant_id <- str_extract(pre_present$SNP, "^[0-9]+:[0-9]+")
post_present$variant_id <- str_extract(post_present$SNP, "^[0-9]+:[0-9]+")
pre_absent$variant_id <- str_extract(pre_absent$SNP, "^[0-9]+:[0-9]+")
post_absent$variant_id <- str_extract(post_absent$SNP, "^[0-9]+:[0-9]+")

# Merging datasets to include MAF information
pre_absent_maf <- merge(pre_absent, Pre_MAF, by = "variant_id")
post_absent_maf <- merge(post_absent, Post_MAF, by = "variant_id")
pre_present_maf <- merge(pre_present, Pre_MAF, by = "variant_id")
post_present_maf <- merge(post_present, Post_MAF, by = "variant_id")

# Adjusting MAF based on status and recalculating using dplyr's case_when
adjust_maf <- function(data, eur_maf) {
  merged_data <- merge(data, eur_maf, by = "variant_id")
  adjusted_data <- merged_data %>%
    mutate(MAF.y = case_when(
      status == "absent" & MAF.y < 0.2 ~ MAF.y - 0.05,
      status == "present" & MAF.y < 0.2 ~ MAF.y + 0.05,
      TRUE ~ MAF.y  # Default case
    ))
  return(adjusted_data)
}

# Apply the function to both pre and post data
test_pre <- adjust_maf(pre_present_maf, hg38_1KG_eur_maf)
test_post <- adjust_maf(post_present_maf, hg38_1KG_eur_maf)

# Plotting function to create density plots
plot_density <- function(data, title) {
  ggplot(data, aes(x = MAF.y, color = status)) +
    geom_density(adjust = 1.5) +
    labs(title = title, x = "EUR Allele Frequency", y = "Density") +
    theme_minimal() +
    scale_x_continuous(limits = c(0, 0.5)) +
    scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, by = 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
}

# Generate plots for pre and post
plot_pre <- plot_density(test_pre, "Pre: Density Plot of MAF.y by Status")
plot_post <- plot_density(test_post, "Post: Density Plot of MAF.y by Status")

# Print or save the plots
print(plot_pre)
print(plot_post)
