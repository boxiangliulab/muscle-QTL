# 4. sQTL mapping - 5.ethnicity_specific_QTL.R

# Load necessary libraries
library(ggplot2)

# Remove the string "chr" from column V1 in the AF_data dataframe
AF_data$V1 <- gsub("chr", "", AF_data$V1)

# Combine columns V1 to V4 into a new column named 'combined'
AF_data$combined <- paste(AF_data$V1, AF_data$V2, AF_data$V3, AF_data$V4, sep = ":")

# Calculate the Minor Allele Frequency (MAF) for each SNP
AF_data$MAF <- pmin(AF_data$V5, 1 - AF_data$V5)

# Create a combined column in the European frequency data
eur_frq$combined <- paste(eur_frq$chromosome.start, eur_frq$A2, eur_frq$A1, sep = ":")

# Merge SNP datasets with allele frequency data for both pre and post datasets
pre_sig_af <- merge(common_pre_sig, all_af, by = "SNP")
pre_gtex_af <- merge(common_pre_gtex, all_af, by = "SNP")
pre_not_af <- merge(pre_not_in_gtex, all_af, by = "SNP")

post_sig_af <- merge(common_post_sig, all_af, by = "SNP")
post_gtex_af <- merge(common_post_gtex, all_af, by = "SNP")
post_not_af <- merge(post_not_in_gtex, all_af, by = "SNP")

# Assign data to different groups for analysis clarity
pre_sig_af$group <- 'Pre-significant AF'
pre_gtex_af$group <- 'Pre-GTEx AF'
pre_not_af$group <- 'Pre-not in GTEx AF'

post_sig_af$group <- 'Post-significant AF'
post_gtex_af$group <- 'Post-GTEx AF'
post_not_af$group <- 'Post-not in GTEx AF'

# Select specific columns to create subsets for plotting
pre_sig_af_sub <- pre_sig_af[, c("eur_maf", "group")]
pre_gtex_af_sub <- pre_gtex_af[, c("eur_maf", "group")]
pre_not_af_sub <- pre_not_af[, c("eur_maf", "group")]

post_sig_af_sub <- post_sig_af[, c("eur_maf", "group")]
post_gtex_af_sub <- post_gtex_af[, c("eur_maf", "group")]
post_not_af_sub <- post_not_af[, c("eur_maf", "group")]

# Combine pre and post datasets to form complete datasets for analysis
combined_df_pre <- rbind(pre_sig_af_sub, pre_gtex_af_sub, pre_not_af_sub)
combined_df_post <- rbind(post_sig_af_sub, post_gtex_af_sub, post_not_af_sub)

