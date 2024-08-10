# 5.eQTL mapping - 07.treeQTL.R

library(data.table)
library(stringr)
library(cowplot)
library(dplyr)
library(TreeQTL)
library(qvalue)

# Define significance levels for the analysis
level1 = 0.05
level2 = 0.05

# Encapsulate the file reading and analysis process into a function
process_QTL <- function(input_file, output_prefix) {
    # Read the data from the specified input file
    dt <- fread(input_file)

    # Generate gene and SNP maps based on unique identifiers
    gene_map <- get_gene_map(unique(dt$fid))
    snp_map <- get_snp_map(dt$sid)
    names(snp_map) <- c("family", "chr", "pos")
    names(gene_map) <- c("family", "chr", "s1", "s2")

    # Display the first few rows of the gene and SNP maps
    print("Gene map:")
    head(gene_map)

    print("SNP map:")
    head(snp_map)

    # Prepare mQTL data for TreeQTL analysis
    meqtl = dt[, .(SNP = sid, gene = fid, beta = -1, `t-stat` = -1, `p-value` = pval, FDR = -1)]
    # Run TreeQTL analysis
    temp = treeQTL(meqtl, snp_map, gene_map, level1 = level1, level2 = level2, eSNP = FALSE)

    # Extract and order results
    eGenes = temp[[1]]
    eAssocs = temp[[2]]
    setorder(eGenes, fam_p)
    setorder(eAssocs, BBFDR)

    # Print and write the number of eGenes and eAssociations
    print("Number of eGenes:")
    print(dim(eGenes))
    write.table(eGenes, paste0(output_prefix, "/eGene_test.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    print("Number of eAssociations:")
    print(dim(eAssocs))
    write.table(eAssocs, paste0(output_prefix, "/eAssocs_test.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    # Print the first few rows of eGenes and eAssociations
    print("eGene:")
    print(head(eGenes))

    print("eAssocs:")
    print(head(eAssocs))
}

# Process pre and post experiment data separately
process_QTL("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/pre_rasqual_output_full.txt", "/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/TreeQTL/Pre")
process_QTL("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/post_rasqual_output_full.txt", "/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/TreeQTL/Post")

