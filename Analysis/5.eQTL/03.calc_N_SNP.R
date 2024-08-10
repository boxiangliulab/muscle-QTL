# 5.eQTL mapping - 03.calc_N_SNP.R

library(dplyr)
library(devtools)
library(rasqualTools)

for (i in 1:22) {
  gene_data <- read.table(paste("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/12.rasqual/gene_data/chr", i, ".txt", sep = ""), header = TRUE)
  snp_coords <- read.table(paste("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/12.rasqual/snp_coords/snp_coords_input_byCHR/chr", i, ".txt", sep = ""), header = TRUE)
  
  snp_counts <- countSnpsOverlapingExons(gene_data, snp_coords, cis_window = 5e5)
  selected_data <- dplyr::select(snp_counts, gene_id, feature_snp_count, cis_snp_count)
  
  output_file <- paste("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/12.rasqual/input_filt/rasqual.input.chr", i, ".filt.txt", sep = "")
  write.table(selected_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

