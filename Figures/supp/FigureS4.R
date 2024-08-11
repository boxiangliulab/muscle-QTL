# Plot correlation heatmap using corrplot
pdf("Cor_heatmap.pdf")
corrplot(cor_matrix, method = "color", col = rev(col2(200)),
         addtl.cex = 0.6, insig = "label_sig",
         tl.col = "black", mar = c(0, 0, 1, 0), hclust.method = "complete", order = "hclust")
dev.off()

# Define annotations for rows (or columns)
annotation_row <- data.frame(
  dataClass = c("Basic Information", "Basic Information", "Basic Information",
                "Measurements", "Measurements", "Measurements", "Measurements",
                "Measurements", "Fat Distribution", "Fat Distribution",
                "Fat Distribution", "Fat Distribution", "Fat Distribution",
                "Fat Distribution", "Fat Distribution", "Metabolism",
                "Metabolism", "Metabolism", "Metabolism", "Metabolism",
                "Metabolism", "Metabolism", "Energy Expenditure",
                "Energy Expenditure", "Energy Expenditure", "Energy Expenditure")
)

# Define column labels
label <- c("Age", "Birth Weight", "Gym Visits", "Weight", "BMI", "Waist Circ",
           "Hip Circ", "Waist Hip Ratio", "Total Whole Body Fat",
           "Total Whole Body Lean", "% Liver Fat", "SAT", "VAT", "dSAT",
           "sSAT", "Fasting Insulin", "Fasting Glucose", "ISI", "IMCR",
           "Carbohydrate Expenditure", "Fat Expenditure", "Cholesterol",
           "Total Energy Expenditure", "Accelerometer", "Resting Metabolic Rate",
           "Respiratory Quotient")

# Assign row names to annotation data frame
rownames(annotation_row) <- label

# Verify annotation data frame and labels match
stopifnot(all(rownames(annotation_row) == colnames(cor_matrix)))

# Define custom colors for annotations
annotation_colors <- list(
  dataClass = c(
    "Basic Information" = "#fbb4ae",
    "Measurements" = "#b3cde3",
    "Fat Distribution" = "#ccebc5",
    "Metabolism" = "#decbe4",
    "Energy Expenditure" = "#fed9a6"
  )
)

# Assign labels to correlation matrix
colnames(cor_matrix) <- label
rownames(cor_matrix) <- label

# Create heatmap using pheatmap
pheatmap(
  cor_matrix,
  border_color = NA,
  color = rev(col2(200)),
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE,
  display_numbers = TRUE,
  main = "Changes in Clinical Information",
  angle_col = 45,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
  )
