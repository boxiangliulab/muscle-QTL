# 7. colocalization - 04.colocalization.R

library("coloc")
library(dplyr)
library(data.table)

# Set the base directories for GWAS and inter loci
gwas_dir <- "/home/project/11003054/share/data/T2D_related_GWAS/processed_files/"
base_inter_loci_dir <- "/home/project/11003054/e1101919/muscle_QTL/RNAseq/"

# Function to extract gene names from filenames
extract_genename <- function(filename) {
    gsub(".txt", "", basename(filename))
}

# Iterate over each analysis type
for (analysis_type in c("eQTL", "sQTL")) {
    # Adjust directories based on the analysis type
    inter_loci_dir <- paste0(base_inter_loci_dir, "11.coloc/inter_loci/", analysis_type, "/")
    nominal_snp_dir <- paste0(base_inter_loci_dir, "11.coloc/inter_loci/", analysis_type, "/")
    output_base_dir <- paste0(base_inter_loci_dir, "11.coloc/coloc_results/", analysis_type, "/")

    for (celltype in c("Pre", "Post")) {
        output_dir <- paste0(output_base_dir, celltype, "/")

        # List all GWAS files
        gwas_files <- list.files(gwas_dir, full.names = TRUE, pattern = "\\.txt$")

        for (gwas_file in gwas_files) {
            genename <- extract_genename(gwas_file)
            GWAS_pre <- fread(gwas_file, sep="\t", header=TRUE)

            # Set file paths for gene and SNP data
            gene_file_path <- paste0(inter_loci_dir, celltype, "/", genename, "_1e_7.txt/gene.txt")
            snp_file_path <- paste0(nominal_snp_dir, celltype, "_nominal_snp/", genename, "_1e_7.txt")

            if (!file.exists(gene_file_path) || !file.exists(snp_file_path)) {
                next  # Skip if any file doesn't exist
            }

            loci <- read.table(gene_file_path, sep="\t")

            # Load MAF data
            MAF <- fread("/home/project/11003054/e1101919/muscle_QTL/RNAseq/00.data/raw_genotype/raw_genotype_hg19_organized/impute/imputed_results/merged/calc_vcf_AS/AS_maf0.05/vcf_maf/SAMs2_calc_maf_freq.frq")
            MAF$variant_id <- sapply(strsplit(MAF$SNP, ":"), function(x) paste(x[1], x[2], sep = ":"))
            MAF$variant_id <- gsub("chr", "", MAF$variant_id)

            print("Finished preparation.")

            if (nrow(loci) > 0) {
                test_rs <- as.data.frame(matrix(NA, nrow(loci), 8))
                colnames(test_rs) <- c("variant_id", "gene", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

                for (j in 1:nrow(loci)) {
                    snp_data <- fread(paste(snp_file_path, "/", loci[j, 1], ".txt", sep = ""), sep = "\t")
                    snp_data <- snp_data[, .(variant_id, pval_nominal, beta)]
                    snp_data <- snp_data[order(pval_nominal)]

                    input <- merge(snp_data, MAF, by = "variant_id")
                    input <- merge(input, GWAS_pre, by = "variant_id", suffixes = c("_sqtl", "_gwas"))
                    input <- input[order(input$pval_nominal_gwas),]

                    if (nrow(input) > 10) {
                        result <- coloc.abf(dataset1 = list(pvalues = input$pval_nominal_gwas, type = "quant", beta = input$beta_gwas, varbeta = input$se_gwas^2, N = 80000, snp = input$variant_id),
                                            dataset2 = list(pvalues = input$pval_nominal_sqtl, type = "quant", N = 54, snp = input$variant_id), MAF = input$MAF)
                        test_rs[j, 3:8] <- t(as.data.frame(result$summary))[1, 1:6]
                    }
                }

                # Save results
                output_file_path <- paste0(output_dir, genename, "_", celltype, "_coloc.txt")
                if (nrow(test_rs) > 0) {
                    write.table(test_rs, output_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
                }
            }
        }
    }
}
print("Script executed successfully.")

