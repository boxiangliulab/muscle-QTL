# muscle bulk - 07.WGCNA.R
# Load required packages
library("dynamicTreeCut")
library("fastcluster")
library("WGCNA")
options(stringsAsFactors = FALSE)
allowWGCNAThreads(18)

# Read expression data
expo <- read.csv("selected_TPM_bulk.txt", sep = '\t', row.names = 1)

# Check the data structure
dim(expo)
names(expo)
datExpr <- as.data.frame(t(expo))

# Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

# Remove samples and genes with missing data if necessary
if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
        printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
    }
    if (sum(!gsg$goodSamples) > 0) {
        printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
    }
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Filter out lowly expressed genes
meanFPKM <- 0.5
n <- nrow(datExpr)
datExpr[n + 1, ] <- apply(datExpr[1:n, ], 2, mean)
datExpr <- datExpr[1:n, datExpr[n + 1, ] > meanFPKM]

# Transpose the data back for further analysis
filtered_fpkm <- t(datExpr)
filtered_fpkm <- data.frame(rownames(filtered_fpkm), filtered_fpkm)
names(filtered_fpkm)[1] <- "sample"
head(filtered_fpkm)
write.table(filtered_fpkm, file = "FPKM_filter.xls", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')

# Cluster samples
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file = "1_sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Read trait data
traitData <- read.csv("trits.txt", sep = '\t', row.names = 1)
alltraits <- traitData

# Align trait data with expression data
fpkmSamples <- rownames(datExpr)
traitsamples <- rownames(alltraits)
traitRows <- match(fpkmSamples, traitsamples)
dataTraits <- alltraits[traitRows, ]
rownames(dataTraits)
collectGarbage()

# Re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(dataTraits, signed = FALSE)
pdf(file = "2_sample dendrogram and trait heatmap.pdf", width = 15, height = 12)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(dataTraits), main = "Sample dendrogram and trait heatmap")
dev.off()

# Network construction
enableWGCNAThreads()
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf(file = "3_Scale independence.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = 0.9, col = "red")
abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
dev.off()

# Choose the appropriate soft-thresholding power and construct a network
softPower <- sft$powerEstimate
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf(file = "4_Gene clustering on TOM-based dissimilarity.pdf", width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# Define module minimum size and identify modules using dynamic tree cut
minModulesize <- 10
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModulesize)
dynamicColors <- labels2colors(dynamicMods)
pdf(file = "5_Dynamic Tree cut.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Save processed data
save(datExpr, dataTraits, file = "fpkm_and_traits_forAnalysis.RData")

# Process and export gene significance and module membership
geneInfo0 <- data.frame(probes = rownames(datExpr), moduleColor = moduleColors)
for (Tra in 1:ncol(dataTraits)) {
    oldNames <- names(geneInfo0)
    geneInfo0 <- data.frame(geneInfo0, geneTraitSignificance[, Tra], GSPvalue[, Tra])
    names(geneInfo0) <- c(oldNames, names(geneTraitSignificance)[Tra], names(GSPvalue)[Tra])
}
write.table(geneInfo0, file = "10_GS_and_MM.xls", sep = "\t", row.names = FALSE)
