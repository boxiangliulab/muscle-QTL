# muscle bulk - 08.rsem_isoform_quantification.sh
#!/bin/bash

# Set directories and parameters
input_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/01.wasp_Mapping/aligned.toTranscriptome.out.bam/merged_bam"
converted_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/01.wasp_Mapping/aligned.toTranscriptome.out.bam/merged_bam/rsem_bam"
annotation_gtf="/home/project/11003054/e1101919/muscle_QTL/RNAseq/18.isoform_quantification/rsem/reference/reference_name"
output_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/18.isoform_quantification/rsem/isoform_expression"
num_threads_convert=64
num_threads_rsem=8

# Create output directories if they do not exist
mkdir -p "$converted_dir"
mkdir -p "$output_dir"

# Convert BAM files for RSEM
echo "Starting BAM file conversion for RSEM..."
for bam_file in "$input_dir"/*.bam; do
    filename=$(basename "$bam_file" .bam)
    output_file="$converted_dir/${filename}_rsem.bam"
    convert-sam-for-rsem "$bam_file" "$output_file" -p $num_threads_convert
    echo "Converted $bam_file to $output_file"
done
echo "All BAM files have been converted to RSEM format."

# Calculate expression using RSEM
echo "Starting RSEM calculation..."
for bam_file in "$converted_dir"/*.bam; do
  out_prefix=$(basename "$bam_file" .bam)
  output_path="$output_dir/$out_prefix"
  rsem-calculate-expression --paired-end --no-bam-output --alignments -p $num_threads_rsem "$bam_file" "$annotation_gtf" "$output_path"
done
echo "RSEM expression calculation completed."
