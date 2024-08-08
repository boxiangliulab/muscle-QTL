# Genotype - 01.after imputation qc

#!/bin/bash

# Directory containing the VCF files
input_dir="/home/project/11003054/e1101919/SAMS2/genotype/SAMS2/imputation_results"
output_dir="/home/project/11003054/e1101919/SAMS2/genotype/SAMS2/imputation_results/filtered"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over all .vcf.gz files in the input directory
for file in "$input_dir"/*.vcf.gz; do
  # Extract the file name without extension
  filename=$(basename "$file" .vcf.gz)

  # Calculate HWE and save intermediate file
  bcftools +fill-tags "$file" -- -t HWE | \
  bcftools view -Oz -o "$output_dir/${filename}_with_hwe.vcf.gz"

  # Filter the VCF file based on MAF and HWE
  bcftools view -m2 -M2 -v snps "$output_dir/${filename}_with_hwe.vcf.gz" | \ # only keep biallelic snps
  bcftools filter -Oz -e 'INFO/MAF < 0.05 | INFO/HWE < 1e-6 | INFO/R2 < 0.8' -o "$output_dir/${filename}_biallelic_filtered.vcf.gz"

  # Index the filtered VCF file
  bcftools index -t "$output_dir/${filename}_biallelic_filtered.vcf.gz"
  
  # Remove intermediate file
  rm "$output_dir/${filename}_with_hwe.vcf.gz"
done

echo "Filtering complete. Filtered files are located in $output_dir."

