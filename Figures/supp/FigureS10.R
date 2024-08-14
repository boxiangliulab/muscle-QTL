library(vcfR)
library(ggplot2)

# Function to extract genotype data from a VCF file for a specific SNP
extract_genotypes <- function(vcf_path, snp_id) {
    vcf <- read.vcfR(vcf_path)

    snp_index <- which(vcf@fix[, "ID"] == snp_id)
    if (length(snp_index) == 0) {
        stop("SNP ID not found in the VCF file.")
    }

    genotypes <- extract.gt(vcf, element = "GT", snp.index = snp_index)

    readable_genotypes <- apply(genotypes, 2, function(g) {
        gsub("0", vcf@ref[snp_index], g)
        gsub("1", vcf@alt[snp_index, 1], g) # Assuming biallelic SNP for simplicity
    })

    # Return a dataframe with Sample IDs and their genotypes
    data.frame(SampleID = names(readable_genotypes), Genotype = readable_genotypes, stringsAsFactors = FALSE)
}




# Function to plot gene expression by genotype
plot_gene_expression <- function(expr_data, gene_id, genotype_data) {
    gene_expr <- as.data.frame(t(expr_data[gene_id, ]), stringsAsFactors = FALSE)
    colnames(gene_expr) <- "Expression"
    gene_expr$SampleID <- rownames(gene_expr)

    merged_data <- merge(gene_expr, genotype_data, by = "SampleID")

    merged_data$Genotype <- factor(merged_data$Genotype)

    ggplot(merged_data, aes(x = Genotype, y = Expression)) +
        geom_boxplot(aes(fill = Genotype)) +
        labs(title = paste("Expression of", gene_id, "by Genotype"), x = "Genotype", y = "Expression") +
        theme_minimal()
}

# Examples
vcf_path <- "SAMS2.vcf"
snp_id <- "rs12437434"
gene_id <- "ENSG00000187630"
expr_data <- corrected_TPM_SAMS2_bulk

genotype_data <- extract_genotypes(vcf_path, snp_id)
expression_plot <- plot_gene_expression(expr_data, gene_id, genotype_data)
print(expression_plot)
