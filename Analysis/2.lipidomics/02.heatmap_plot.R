# Lipidomics - 02. heatmap_plot
# Load necessary packages
library(readr)
library(dplyr)
library(pheatmap)

# Load data
load("~/OneDrive - National University of Singapore/SAMS2_bulk/github/lipidomics_data.RData")

# Define file paths
original_data_file <- paste(f[i], "data_original.csv", sep = "")
diff_quant_data_file <- paste(out, "3.1DiffQuantData.txt", sep = "")
heatmap_plot_file <- paste(out, "04.heatmap.pdf", sep = "")

# Read data
mdata <- read_delim(original_data_file, delim = ",") %>%
  dplyr::rename(Name = ...1)
data <- read_delim(diff_quant_data_file, delim = "\t")

# Extract group names
grpp <- unlist(strsplit(f[i], split = "-vs-"))
indexx <- c(grep(grpp[1], names(mdata)), grep(grpp[2], names(mdata)))

# Filter significant data
hdata <- data %>%
  filter(Sig != "No_Sig") %>%
  select(Name) %>%
  inner_join(mdata, by = "Name") %>%
  select(Name, indexx) %>%
  as.data.frame() %>%
  column_to_rownames(var = "Name")

# Class labels
classLabel <- data.frame(t(hdata[1, ]))
classLabel[, 1] <- row.names(classLabel)
names(classLabel) <- "class"
classLabel$class <- sapply(classLabel$class, function(x) strsplit(x, "_")[[1]][1])

# Annotation colors
ann_colors <- list(class = c(Post = "#e6194b", Pre = "#3cb44b"))

# Convert data to numeric
hdata <- apply(hdata, 2, as.numeric)

# Plot heatmap
pheatmap(hdata, annotation = classLabel, scale = "row", fontsize = 7,
         annotation_colors = ann_colors, fontsize_row = 6, fontsize_col = 5,
         color = colorRampPalette(c("blue", "white", "red"))(255), show_colnames = TRUE,
         show_rownames = TRUE, cluster_cols = TRUE, angle_col = "45",
         filename = heatmap_plot_file, width = 10, height = 7, border_color = "grey")

# Save processed heatmap data
write.table(hdata, file = paste(out, "03.heatmap.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
