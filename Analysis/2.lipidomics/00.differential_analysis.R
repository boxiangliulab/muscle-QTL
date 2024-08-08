# Lipidomics - 00. differential_analysis
# Load necessary packages
library(MetaboAnalystR)
library(dplyr)
library(readr)

# Load data
load("~/OneDrive - National University of Singapore/SAMS2_bulk/github/lipidomics_data.RData")

# Create a MetaboAnalystR object
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "~/OneDrive - National University of Singapore/SAMS2_bulk/github/lipidomics_data.txt", "rowu", "disc")

# Add data to MetaboAnalystR object
mSet$dataSet$preData <- Pre_metabolite_54ind
mSet$dataSet$postData <- Post_metabolite_54ind

# Normalize the data
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "NULL")

# Perform differential analysis
mSet <- Ttests.Anal(mSet)
mSet <- SAM.Anal(mSet)

# Get differential analysis results
ttest_results <- mSet$analSet$tt$sig.mat
sam_results <- mSet$analSet$sam$sig.mat

# Save differential analysis results
write.csv(ttest_results, "~/OneDrive - National University of Singapore/SAMS2_bulk/github/ttest_results.csv")
write.csv(sam_results, "~/OneDrive - National University of Singapore/SAMS2_bulk/github/sam_results.csv")

# Input modilfed-file paths
vip_file <- paste(f[i], "oplsda_vip.csv", sep = "")
fcp_file <- paste(f[i], "volcano.csv", sep = "")
output_file <- paste(out, "DiffQuantData.txt", sep = "")

# Read data
vip <- read_delim(vip_file, delim = ",")
fcp <- read_delim(fcp_file, delim = ",")

# Merge and process data
data <- vip %>%
  inner_join(fcp, by = c("...1" = "...1")) %>%
  dplyr::rename(VIP = V1, Name = ...1, pvalue = raw.pval) %>%
  select(Name, VIP, FC, pvalue) %>%
  mutate(t.test_qvalue = p.adjust(pvalue, method = "fdr"))

# Define significance
data$Sig <- "No_Sig"
data$Sig[data$FC >= 1.2 & data$t.test_qvalue < 0.05] <- "Up"
data$Sig[data$FC <= 1/1.2 & data$t.test_qvalue < 0.05] <- "Down"

# Output data
write_delim(data, output_file, delim = "\t")

# Output significance count
significance_count <- as.data.frame(table(data$Sig)) %>%
  t() %>%
  as.data.frame()
write_delim(significance_count, paste(out, "DiffNum.txt", sep = ""), delim = "\t", col_names = FALSE)