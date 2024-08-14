# Figure S5A
library(Rsamtools)
library(ggplot2)

# List all BAM files
bam_files <- list.files(path = "/home/project/11003054/e1101919/muscle_QTL/RNAseq/01.wasp_Mapping/combined_sorted_bam", pattern = "*.bam", full.names = TRUE)

# Function to get read depth from a BAM file
get_read_depth <- function(bam_file) {
  idxstats <- idxstatsBam(bam_file)
  sum(idxstats$reads)
}

# Apply the function to each BAM file
read_depths <- sapply(bam_files, get_read_depth)

# Create a histogram
ggplot(data.frame(Read_Depth = read_depths), aes(x = Read_Depth)) +
  geom_histogram(binwidth = 1e6, fill = "gray") +
  labs(title = "SAMS2 Sample Read Depths", x = "Sample Name", y = "Read Depth") +
  theme_minimal()

# Figure S5B
library(ggplot2)

# Extract this data from the BAM files

# Creating the scatter plot
ggplot(chromosome_coverage, aes(x = X_Coverage, y = Y_Coverage)) +
  geom_point() +
  labs(title = "Sex Determination Based on Chromosome Coverage", x = "X Chromosome Coverage", y = "Y Chromosome Coverage") +
  theme_minimal()

# Figure S5C
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
