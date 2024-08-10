# 7. colocalization - 01.prepare_input_QTLs.sh
#!/bin/bash

# Define directories
input_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/nominal_pass/all/phenotype_test_cov2"
# input_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/nominal_pass/all/phenotype_test_cov2"

output_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/11.coloc/all_eQTL/nominal_overstandard"
# output_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/11.coloc/all_sQTL/nominal_overstandard"


# Array of time points
timepoints=("Pre" "Post")

# Loop through each time point
for timepoint in "${timepoints[@]}"; do
    # Ensure the output directory for each timepoint exists
    mkdir -p "$output_dir/$timepoint"

    # Process each chromosome
    for i in $(seq 1 22); do
        echo "Processing ${timepoint} chromosome ${i}"

        # Path to the input file
        input_file="${input_dir}/${timepoint}_nominals_chr${i}.txt"
        # Intermediate output file path
        intermediate_output_file="${output_dir}/${timepoint}/chr${i}.txt"
        # Final output file path
        final_output_file="${output_dir}/${timepoint}/newchr${i}.txt"

        # Filter records with p-value less than 1e-3 and write to intermediate file
        awk '$12 < 1e-3' "$input_file" > "$intermediate_output_file"

        # Remove 'chr' prefix from the 9th column and write to final output file
        awk '{ gsub(/chr/,"", $9); print }' "$intermediate_output_file" > "$final_output_file"
    done
done

echo "All files have been processed."

