# 7. colocalization - 02.inter_GWAS_QTLs.sh

args <- commandArgs(TRUE)
library(data.table)

# Determine cell type and analysis type (eQTL or sQTL) from command line arguments
celltype <- args[1]
analysis_type <- args[2]  # "eQTL" or "sQTL"

# Set base directory paths depending on the analysis type
if (analysis_type == "eQTL") {
    data_dir <- "/home/project/11003054/e1101919/muscle_QTL/RNAseq/13.ras_coloc/eQTL_data/"
    output_base_dir <- "/home/project/11003054/e1101919/muscle_QTL/RNAseq/13.ras_coloc/inter_loci/"
} else if (analysis_type == "sQTL") {
    data_dir <- "/home/project/11003054/e1101919/muscle_QTL/RNAseq/11.coloc/all_sQTL/nominal_overstandard/"
    output_base_dir <- "/home/project/11003054/e1101919/muscle_QTL/RNAseq/11.coloc/inter_loci/sQTL/"
}

gwas_dir <- "/home/project/11003054/share/data/T2D_related_GWAS/prepare_required_columns/threshold_1e-7/"
file_list <- list.files(gwas_dir, pattern="\\.txt$", full.names=TRUE)

for (file in file_list) {
    GWAS_pre <- fread(file, sep="\t", header=TRUE)
    GWAS_pre <- as.data.frame(GWAS_pre)
    GWAS_pre <- GWAS_pre[, c(2, 3, 4, 5, 1)]
    colnames(GWAS_pre) <- c("variant_id", "beta", "se", "Pvalue", "rsID")

    for (i in 1:22) {
        if (analysis_type == "eQTL") {
            nominal_file <- paste(data_dir, celltype, "/chr", i, ".txt", sep="")
        } else {
            nominal_file <- paste(data_dir, celltype, "/new/newchr", i, ".txt", sep="")
        }

        if (file.exists(nominal_file)) {
            nominal <- fread(nominal_file, sep="\t")
            nominal <- as.data.frame(nominal)
            nominal$variant_id <- sub("(:[^:]+){2}$", "", nominal$sid)

            if (analysis_type == "sQTL") {
                colnames(nominal) <- c("intron_clu", "chr", "start", "end", "strand", "distance", "V7", "sid", "chr_num", "V10", "V11", "V12", "V13", "V14")
                nominal$variant_id <- gsub("chr", "", nominal$variant_id)
            }

            print("Read in nominal")

            common_snp <- merge(nominal, GWAS_pre, by = "variant_id", all = FALSE)
            common_snp <- na.omit(common_snp)

            if (nrow(common_snp) > 0) {
                output_dir <- paste(output_base_dir, celltype, "/", basename(file), sep = "")
                if (!dir.exists(output_dir)) {
                    dir.create(output_dir, recursive = TRUE)
                }
                output_file <- paste(output_dir, "/", i, "loci.txt", sep = "")
                write.table(common_snp, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
            }
        }
    }
}
print("Script executed successfully.")

