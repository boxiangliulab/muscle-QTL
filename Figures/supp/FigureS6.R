# Figure S6A
library(Rsamtools)
library(ggplot2)

# Create the histogram from STAR outpus 
ggplot(data.frame(UMR = Uniquely_Mapped_Reads), aes(x = Sample_id)) +
  geom_histogram(binwidth = 1e6, fill = "gray") +
  labs(title = "SAMS2 Sample Read Depths", x = "Sample Name", y = "Read Depth") +
  theme_minimal()

# Figure S6B
library(ggplot2)

# Extract this data from the BAM files

# Creating the scatter plot
ggplot(chromosome_coverage, aes(x = X_Coverage, y = Y_Coverage)) +
  geom_point() +
  labs(title = "Sex Determination Based on Chromosome Coverage", x = "X Chromosome Coverage", y = "Y Chromosome Coverage") +
  theme_minimal()

# Figure S6C
library(ggplot2)
library(ggfortify)

# Load count dat
count_data <- read.csv("SAMS2_count_table.csv") 

# Perform PCA
pca_result <- prcomp(count_data, scale. = TRUE)

# Plotting the PCA
autoplot(pca_result, data = count_data, label = TRUE, label.size = 3) +
  labs(title = "PCA of SAMS2 Count Data") +
  theme_minimal()
