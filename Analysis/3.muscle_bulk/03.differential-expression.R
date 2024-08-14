# 3.muscle bulk - 03. differential expression
library(edgeR)
library(statmod)

# Load saved RData files
load("SAMS2_muscle_count_table.RData")  # loads expression_res
load("SAMS2_muscle_TPM.RData")          # loads tpm
load("sample_info.RData")               

group <- factor(ifelse(grepl("Pre", colnames(expression_res)), "Pre", "Post"))
group <- relevel(group, "Pre")  # Set the control group to 'Pre'

sample_info <- sample_info[match(colnames(expression_res), sample_info$ID), ]

design <- model.matrix(~ group + Age + factor(ID))

dge <- DGEList(counts = expression_res, group = group)
dge <- calcNormFactors(dge, method = 'TMM')  # Normalize library sizes
dge <- estimateDisp(dge, design, robust = TRUE)

# Fit the model and perform quasi-likelihood F-test
fit <- glmQLFit(dge, design, robust = TRUE)
lt <- glmQLFTest(fit, contrast = c(-1, 1, rep(0, ncol(design)-2)))  # Adjust contrast to test the effect of 'group'

# Get differentially expressed genes
tempDEG <- topTags(lt, n = Inf)
tempDEG <- as.data.frame(tempDEG)
nrDEG <- na.omit(tempDEG)

# Classify gene expression changes
nrDEG$gene_type <- ifelse(nrDEG$PValue > 0.05, 'Not',
                          ifelse(nrDEG$logFC > 1, 'Up',
                                 ifelse(nrDEG$logFC < -1, 'Down', 'Not')))
table(nrDEG$gene_type)

# Extract up- and down-regulated genes and save to file
up_down <- nrDEG[grep("Up|Down", nrDEG$gene_type), ]
up_down1 <- cbind(rownames(up_down), up_down)
colnames(up_down1)[1] <- "ID"
write.table(up_down1, "DEG.csv", sep = ",", row.names = FALSE, col.names = TRUE)

