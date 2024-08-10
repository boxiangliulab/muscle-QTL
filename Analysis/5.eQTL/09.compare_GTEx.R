# 5. eQTL mapping - 09.compare_GTEx.R

# 加载必要的库
library(data.table)
library(dplyr)
library(qvalue)

# 读取数据
all <- fread("/home/project/11003054/e1101919/muscle_QTL/RNAseq/16.ethnicity_specific/eQTL/GTEx_EUR_sm/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Muscle_Skeletal.allpairs_formatted.txt")
print("GTEx data loaded:")
head(all)

pre <- fread("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/TreeQTL/Pre/eAssocs_test.txt")
print("Pre data loaded:")
head(pre)

post <- fread("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/TreeQTL/Post/eAssocs_test.txt")
print("Post data loaded:")
head(post)

pre_filter <- fread("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/pre_filter.txt")

post_filter <- fread("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/post_filter.txt")

colnames(pre_filter)[1] <- "gene"
colnames(pre_filter)[2] <- "SNP"
colnames(post_filter)[1] <- "gene"
colnames(post_filter)[2] <- "SNP"

pre_combined <- inner_join(pre_filter, pre, by = c("gene", "SNP"))

post_combined <- inner_join(post_filter, post, by = c("gene", "SNP"))

pre_combined <- mutate(pre_combined, gene_SNP = paste(gene, SNP, sep = "_"))
post_combined <- mutate(post_combined, gene_SNP = paste(gene, SNP, sep = "_"))
all <- mutate(all, gene_SNP = paste(gene, SNP, sep = "_"))

pre_inALL <- filter(all, gene_SNP %in% pre_combined$gene_SNP)
post_inALL <- filter(all, gene_SNP %in% post_combined$gene_SNP)

pre_inALL <- select(pre_inALL, -gene_SNP)
post_inALL <- select(post_inALL, -gene_SNP)

fwrite(pre_inALL, "/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/compare_gtex/pre_inALL.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

fwrite(post_inALL, "/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/compare_gtex/post_inALL.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

pre_combined$pi <- pre_combined$pi - 0.5
post_combined$pi <- post_combined$pi - 0.5

result_pre <- qvalue(pre_inALL$pval_nominal)
result_post <- qvalue(post_inALL$pval_nominal)

print("pre_pi0:")
print(result_pre$pi0)

print("post_pi0:")
print(result_post$pi0)

merged_pre <- inner_join(pre_inALL, pre_combined, by = c("gene", "SNP"))

print("merged_pre:")
head(merged_pre)

merged_pre <- merged_pre%>%
  mutate(direction_consistent = (slope > 0 & pi > 0) | (slope < 0 & pi < 0))

merged_post <- inner_join(post_inALL, post_combined, by = c("gene", "SNP"))

print("merged_post:")
head(merged_post)

merged_post <- merged_post %>%
  mutate(direction_consistent = (slope > 0 & pi > 0) | (slope < 0 & pi < 0))

fwrite(merged_pre, "/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/compare_gtex/pre_inALL_direction.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

fwrite(merged_post, "/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/compare_gtex/post_inALL_direction.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


print("pre_direction_gtex:")
table(merged_pre$direction_consistent)

print("post_direction_gtex:")
table(merged_post$direction_consistent)

