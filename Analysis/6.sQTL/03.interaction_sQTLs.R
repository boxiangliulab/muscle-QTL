library(data.table)
library(dplyr)

adjust_genotypes <- function(genotypes_df, expression_df) {
  sample_names <- colnames(expression_df)
  
  new_genotypes <- data.frame(matrix(nrow = nrow(genotypes_df), ncol = length(sample_names)))
  rownames(new_genotypes) <- rownames(genotypes_df)
  colnames(new_genotypes) <- sample_names
  
  genotype_cols <- colnames(genotypes_df)
  
  for (sample in sample_names) {
    if (startsWith(sample, "Pre_")) {
      id <- substr(sample, 5, nchar(sample))
    } else if (startsWith(sample, "Post_")) {
      id <- substr(sample, 6, nchar(sample))
    } else {
      warning(paste("Unexpected sample name format:", sample))
      next
    }
    
    possible_ids <- c(paste0("ID_", id), id, paste0("id", id))
    matching_col <- genotype_cols[genotype_cols %in% possible_ids]
    
    if (length(matching_col) == 1) {
      new_genotypes[, sample] <- genotypes_df[, matching_col]
    } else {
      warning(paste("No matching column found for sample:", sample))
    }
  }
  
  return(new_genotypes)
}

genotypes_encoded_adj <- genotypes_encoded_adj[, colnames(SAM2_count_table_final_adj)]

if (!all(colnames(genotypes_encoded_adj) == colnames(SAM2_count_table_final_adj))) {
  stop("Column names (sample order) in genotypes_encoded_adj and SAM2_count_table_final_adj do not match!")
}

sQTL_both_significant <- eQTL_summary_statistics_adj %>%
  group_by(eGene, SNP) %>%
  filter(n_distinct(condition) == 2) %>%
  ungroup()

prepare_data <- function(gene, snp) {
  splicing_ratio <- as.numeric(SAM2_count_table_final_adj[gene, ])
  genotype <- as.numeric(genotypes_encoded_adj[snp, ])
  
  data <- data.frame(
    sample_id = colnames(SAM2_count_table_final_adj),
    splicing_ratio = splicing_ratio,
    genotype = genotype
  )
  
  data <- merge(data, lifestyle, by.x = "sample_id", by.y = "sample_id", all.x = TRUE)
  data <- merge(data, t(cov3_adj), by.x = "sample_id", by.y = "row.names", all.x = TRUE)
  
  return(data)
}

analyze_sQTL <- function(gene, snp) {
  data <- prepare_data(gene, snp)
  data$condition <- factor(data$condition, levels = c("pre", "post"))
  
  print(paste("Analyzing:", gene, snp))
  print(table(data$condition, data$genotype))
  
  model <- lm(splicing_ratio ~ genotype + condition + genotype:condition + PC1 + PC2 + PC3, data = data)
  summary <- summary(model)
  
  interaction_term <- "genotype:conditionpost"
  if (interaction_term %in% rownames(summary$coefficients)) {
    interaction_p <- summary$coefficients[interaction_term, "Pr(>|t|)"]
    interaction_effect <- summary$coefficients[interaction_term, "Estimate"]
  } else {
    warning(paste("Interaction term not found for", gene, snp))
    interaction_p <- NA
    interaction_effect <- NA
  }
  
  return(data.frame(gene = gene, snp = snp, interaction_p = interaction_p, interaction_effect = interaction_effect))
}

results <- lapply(1:nrow(sQTL_both_significant), function(i) {
  tryCatch({
    analyze_sQTL(sQTL_both_significant$eGene[i], sQTL_both_significant$SNP[i])
  }, error = function(e) {
    message(paste("Error in analyzing", sQTL_both_significant$eGene[i], sQTL_both_significant$SNP[i], ":", e$message))
    return(data.frame(gene = sQTL_both_significant$eGene[i],
                      snp = sQTL_both_significant$SNP[i],
                      interaction_p = NA,
                      interaction_effect = NA))
  })
})

results_df <- do.call(rbind, results)
results_df$fdr <- p.adjust(results_df$interaction_p, method = "fdr")
effect_size_threshold <- median(abs(results_df$interaction_effect), na.rm = TRUE)
results_df$is_interaction_sQTL <- results_df$fdr < 0.05 & abs(results_df$interaction_effect) > effect_size_threshold

write.csv(results_df, "sQTL_interaction_results.csv", row.names = FALSE)
interaction_sQTLs <- results_df[results_df$is_interaction_sQTL, ]
print(interaction_sQTLs)

num_interaction_sQTLs <- sum(results_df$is_interaction_sQTL, na.rm = TRUE)
print(paste("Number of interaction sQTLs:", num_interaction_sQTLs))

