library(vcfR)
library(ggplot2)

# Function to extract genotype data from a VCF file for a specific SNP
extract_genotypes <- function(vcf_path, snp_id) {
    # Load VCF file
    vcf <- read.vcfR(vcf_path)

    # Extract specific SNP by ID
    snp_index <- which(vcf@fix[, "ID"] == snp_id)
    if (length(snp_index) == 0) {
        stop("SNP ID not found in the VCF file.")
    }

    # Extract genotype information
    genotypes <- extract.gt(vcf, element = "GT", snp.index = snp_index)

    # Parse genotypes to a more readable format (e.g., AA, AT, TT)
    readable_genotypes <- apply(genotypes, 2, function(g) {
        gsub("0", vcf@ref[snp_index], g)
        gsub("1", vcf@alt[snp_index, 1], g)
    })

    data.frame(SampleID = names(readable_genotypes), Genotype = readable_genotypes, stringsAsFactors = FALSE)
}

plot_gene_splicing <- function(splicing_ratio_data, gene_id, genotype_data) {
    if (!gene_id %in% rownames(splicing_ratio_data)) {
        stop("Gene ID not found in the splicing data.")
    }

    # Extract splicing data for the gene
    gene_splicing <- as.data.frame(t(splicing_ratio_data[gene_id, ]), stringsAsFactors = FALSE)
    colnames(gene_splicing) <- "Spling"
    gene_expr$SampleID <- rownames(gene_splicing)

    # Merge splicing data with genotype data
    merged_data <- merge(gene_splicing, genotype_data, by = "SampleID")

    # Convert genotype to a factor for plotting
    merged_data$Genotype <- factor(merged_data$Genotype)

    # Plotting using ggplot2
    ggplot(merged_data, aes(x = Genotype, y = Splicing)) +
        geom_boxplot(aes(fill = Genotype)) +
        labs(title = paste("Splicing of", gene_id, "by Genotype"), x = "Genotype", y = "Splicing") +
        theme_minimal()
}

