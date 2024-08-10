#7. colocalization - 03.build_GWAS_QTL_data.sh
#!/bin/bash

# Define the analysis types and their respective paths
declare -A dir_paths
dir_paths[eQTL]="/home/project/11003054/e1101919/muscle_QTL/RNAseq/13.ras_coloc/inter_loci"
dir_paths[sQTL]="/home/project/11003054/e1101919/muscle_QTL/RNAseq/11.coloc/inter_loci/sQTL"

declare -A data_paths
data_paths[eQTL]="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/joint_output2"
data_paths[sQTL]="/home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/nominal_pass_intron_clu/All_for_summary_statistics/phenotype_test_cov1"

# Define the patterns for gene file discovery
declare -A gene_patterns
gene_patterns[eQTL]="*/gene.txt"
gene_patterns[sQTL]="*/gene.txt"

# Loop through both eQTL and sQTL analysis types
for analysis in "${!dir_paths[@]}"; do
    echo "Processing ${analysis} data..."

    for celltype in "Pre" "Post"; do
        for gene_file in ${dir_paths[$analysis]}/${celltype}/${gene_patterns[$analysis]}; do
            while IFS= read -r gene; do
                echo "Processing gene: ${gene} for celltype: ${celltype}"

                # Set up output directory
                output_dir="${dir_paths[$analysis]}/${celltype}_nominal_snp/$(basename "$(dirname "$gene_file")")"
                mkdir -p "${output_dir}"

                # Define file patterns and processing based on analysis type
                if [[ "${analysis}" == "eQTL" ]]; then
                    grep -E "\b${gene}_" ${data_paths[$analysis]}/${celltype}/*/*.txt > "${output_dir}/${gene}.txt"
                elif [[ "${analysis}" == "sQTL" ]]; then
                    grep -F "${gene}" ${data_paths[$analysis]}/${celltype}_nominals*.txt | \
                    awk '{print $1 "\t" $8 "\t" $9 "\t" $12 "\t" $13}' > "${output_dir}/${gene}.txt"
                fi

            done < "${gene_file}"
        done
    done
done

echo "Script execution completed."

