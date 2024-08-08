# Clinical Data - 01. pre-process
# Load necessary packages
library(statmod)

# Read expression data and annotation files
expression1 <- read.table("outprefix.txt", header = TRUE, sep = "\t")
expression1 <- expression1[, -c(2:6)]
anno <- read.table("gencode.txt", sep = "\t", header = TRUE)  # Gene annotation file

# Convert ENSG to gene symbol
colnames(expression1)[1] <- "id"
expression <- merge(anno, expression1, by = "id") #30847  108 
expression <- expression[, -c(1, 3:6)]
expression_res <- aggregate(expression[, 2:ncol(expression)], by = list(expression[, 1]), FUN = max)
rownames(expression_res) <- expression_res[, 1]
expression_res <- expression_res[, -1]
expression_res <- expression_res[-which(rowSums(expression_res) == 0), ] #25248  108 # Remove genes not expressed in any sample

# Filter genes with average read count below 10 and with zero counts in more than 20% of samples
filtered_genes <- rowMeans(expression_res) >= 10 & rowSums(expression_res == 0) <= ncol(expression_res) * 0.2
expression_res <- expression_res[filtered_genes, ] # 15658 108

# Read gene lengths
gene_lengths <- read.table("gene_lengths.txt", header = TRUE, sep = "\t")
colnames(gene_lengths) <- c("gene_id", "length")

# Filter genes with average read count below 10 and with zero counts in more than 20% of samples
filtered_genes <- rowMeans(expression_res) >= 10 & rowSums(expression_res == 0) <= ncol(expression_res) * 0.2
expression_res <- expression_res[filtered_genes, ]

# Merge with gene lengths
expression_res <- merge(expression_res, gene_lengths, by.x = "row.names", by.y = "gene_id")
rownames(expression_res) <- expression_res$Row.names
expression_res <- expression_res[, -1]

# Calculate FPKM
counts <- expression_res[, -ncol(expression_res)]
lengths <- expression_res[, ncol(expression_res)]
fpkm <- t(t(counts) / lengths) * 1e9 / colSums(counts)

# Convert FPKM to TPM
tpm <- t(t(fpkm) / colSums(fpkm)) * 1e6

save(expression_res, file = "SAMS2_muscle_count_table.RData")
save(tpm, file = "SAMS2_muscle_TPM.RData"

)
