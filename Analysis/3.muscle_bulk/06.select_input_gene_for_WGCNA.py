# muscle bulk - 06.select_input_gene_for_WGCNA.py
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys

# Receive file names from command line
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the data
df = pd.read_csv(input_file, sep="\t", index_col=0)

# Retrieve columns for Post and Pre conditions
post_columns = [col for col in df.columns if col.startswith('Post')]
pre_columns = [col for col in df.columns if col.startswith('Pre')]

# Filter genes with at least 2/3 of FPKM values greater than 1
threshold_post = len(post_columns) * 2 / 3
threshold_pre = len(pre_columns) * 2 / 3

post_filtered = df[post_columns].apply(lambda x: (x > 1).sum(), axis=1) > threshold_post
pre_filtered = df[pre_columns].apply(lambda x: (x > 1).sum(), axis=1) > threshold_pre

# Combine filtering results from both conditions
filtered_genes = df[post_filtered & pre_filtered]

# Calculate mean values for both groups
mean_post = filtered_genes[post_columns].mean(axis=1)
mean_pre = filtered_genes[pre_columns].mean(axis=1)

# Calculate the coefficient of variation (standard deviation / mean) between the two groups
cv_between_groups = np.abs(mean_post - mean_pre) / ((mean_post + mean_pre) / 2)

# Select genes where the coefficient of variation is greater than 0.2
final_filtered_genes = filtered_genes[cv_between_groups > 0.2]

# Output to new file
final_filtered_genes.to_csv(output_file, sep='\t')

print(f"Processed data has been saved to {output_file}")
