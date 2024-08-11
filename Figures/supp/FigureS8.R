# Figure S8A
library(ggplot2)

ggplot(plink_check_sex, aes(x = F_index)) +
  geom_histogram(bins = 10, fill = "blue", color = "black") +  # You can adjust the number of bins as needed
  labs(title = "Gender Determination Based on F Index", x = "F Index (plink --check-sex)", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the plot title

# Figure S8B
# Load required libraries
library(ggplot2)

# Read PCA data file
dat <- read.table("PCA_merge2_pruned.pca.eigenvec", header=FALSE)
colnames(dat)[1] <- "FID"  # Change column header of first column to FID
colnames(dat)[2] <- "IID"  # Change column header of second column to IID

# Read population data
pop <- read.table("racefile.txt", header=TRUE)
colnames(pop)[3] <- "group"  # Rename the third column to "group"

# Merge PCA data with population data
data <- merge(dat, pop, by=c("IID", "FID"), all=FALSE)

# Ensure 'group' is treated as a factor
data$group <- as.factor(data$group)

# Identify indices for different population groups
OWN <- which(data$group == "OWN")
EUR <- which(data$group == "EUR")
EAS <- which(data$group == "EAS")
AMR <- which(data$group == "AMR")
AFR <- which(data$group == "AFR")
SAS <- which(data$group == "SAS")

# Start a PDF device to save the plot
pdf("pca-ancestry-plot2.pdf")

# Create an empty plot setting the limits
plot(0, 0, pch="", xlim=c(-0.02, 0.04), ylim=c(-0.05, 0.05), xlab="Principal Component 1", ylab="Principal Component 2")

# Add points for each group with different colors and point characters
points(data$V3[EAS], data$V4[EAS], pch=20, col="#66c2a5", cex=0.1)
points(data$V3[AMR], data$V4[AMR], pch=20, col="#fc8d62", cex=0.1)
points(data$V3[AFR], data$V4[AFR], pch=20, col="#8da0cb", cex=0.1)
points(data$V3[EUR], data$V4[EUR], pch=20, col="#e78ac3", cex=0.1)
points(data$V3[SAS], data$V4[SAS], pch=20, col="#a6d854", cex=0.1)
points(data$V3[OWN], data$V4[OWN], pch="+", col="BLACK", cex=0.6)

# Add grid lines for reference
abline(v=-0.01, col="gray32", lty=2)
abline(h=-0.02, col="gray32", lty=2)

# Add a legend to the plot
legend("topright", pch=c(20, 20, 20, 20, 20, 3), labels=c("EUR", "AMR", "AFR", "EAS", "SAS", "OWN"), col=c("#e78ac3", "#fc8d62", "#8da0cb", "#66c2a5", "#a6d854", "BLACK"), bty="o", cex=1)

# Close the PDF device
dev.off()
