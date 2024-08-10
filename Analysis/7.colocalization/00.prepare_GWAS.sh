# 7. colocalization - 00.prepare_GWAS.sh

library(data.table)

# Set the working directory (assuming your .txt files are located in this directory)
setwd("/home/project/11003054/share/data/T2D_related_GWAS/processed_files")

# Get all .txt files in the current directory
file_list <- list.files(pattern = "\\.txt$")

# Specify the output directory
output_dir <- "/home/project/11003054/share/data/T2D_related_GWAS/processed_files/threshold_1e-7"

# Check and create the output directory if it does not exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Loop through each file
for (file_name in file_list) {
    # Read data using fread for faster reading and automatic formatting
    data <- fread(file_name)

    # Filter rows where Pvalue is less than 1e-7
    filtered_data <- data[Pvalue < 1e-7]

    # Define the output file path
    output_file_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file_name)), "_1e_7.txt"))

    # Write the filtered data to a new file
    fwrite(filtered_data, file = output_file_path, sep = "\t", quote = FALSE)
}

