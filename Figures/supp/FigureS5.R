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

# Example data for X and Y chromosome coverage
# Normally you'd extract this data from the BAM files or a summary file
chromosome_coverage <- data.frame(
  Sample = paste("Sample", 1:100),
  X_Coverage = runif(100, 2, 8),
  Y_Coverage = runif(100, 1.5, 3.5)
)

# Creating the scatter plot
ggplot(chromosome_coverage, aes(x = X_Coverage, y = Y_Coverage)) +
  geom_point() +
  labs(title = "Sex Determination Based on Chromosome Coverage", x = "X Chromosome Coverage", y = "Y Chromosome Coverage") +
  theme_minimal()

# Figure S5C
library(ggplot2)
library(ggfortify)

# Load your count data (this should be a matrix or data frame of counts)
count_data <- read.csv("/path/to/SAMS2_count_table.csv")  # Update with the actual path

# Perform PCA
pca_result <- prcomp(count_data, scale. = TRUE)

# Plotting the PCA
autoplot(pca_result, data = count_data, label = TRUE, label.size = 3) +
  labs(title = "PCA of SAMS2 Count Data") +
  theme_minimal()
