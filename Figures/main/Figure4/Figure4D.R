#!/bin/bash

# Define the path to the LDBlockShow binary
LDBlockShow_path="/home/project/11003054/e1101919/software/LDBlockShow/bin/LDBlockShow"

# Specify the input VCF file for the ANK1 gene
vcf_input="/home/project/11003054/e1101919/muscle_QTL/RNAseq/07.splicing/ggsashimi/input_bams_ANK1_pre.tsv"

# Define output parameters
output_base="/home/project/11003054/e1101919/muscle_QTL/RNAseq/17.susie"
output_name="ANK1_susie"

# Define the genomic region of interest
region="chr8:41663473:41667473"

# Specify the input GWAS summary statistics file
gwas_input="/home/project/11003054/e1101919/muscle_QTL/RNAseq/17.susie/pip.value"

# Name of the file containing specific SNP names to highlight
snp_names_file="/home/project/11003054/e1101919/muscle_QTL/RNAseq/17.susie/target.txt"

# Run LDBlockShow to visualize linkage disequilibrium
$LDBlockShow_path \
  -InVCF $vcf_input \
  -OutPut $output_base/$output_name \
  -Region $region \
  -OutPng \
  -SeleVar 2 \
  -InGWAS $gwas_input \
  -ShowGWASSpeSNP \
  -SpeSNPName $snp_names_file

# Check if there is an unknown argument error, adjust command accordingly
if [ $? -ne 0 ]; then
  echo "Error encountered with LDBlockShow. Checking command options."

  # Attempt to run without the problematic parameter
  $LDBlockShow_path \
    -InVCF $vcf_input \
    -OutPut $output_base/$output_name \
    -Region $region \
    -OutPng \
    -SeleVar 2 \
    -InGWAS $gwas_input \
    -ShowGWASSpeSNP \
    -SpeSNPName $snp_names_file
fi

echo "LDBlockShow operation completed."
