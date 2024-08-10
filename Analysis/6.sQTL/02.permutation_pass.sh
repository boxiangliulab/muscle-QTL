# 6.sQTL mapping - 02.permutation_pass.sh

dir=("Pre" "Post")
chromosomes=$(seq 1 22)

output_base_dir="/home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/permutation_pass"

# Generate
for timepoint in "${dir[@]}"; do
    for i in $chromosomes; do
        for cov_type in "genotype_test" "phenotype_test"; do
            for cov_folder in /home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/cov/${cov_type}/cov*; do
                cov_name=$(basename "$cov_folder")
                output_dir="${output_base_dir}/${cov_type}_${cov_name}/${timepoint}"
                mkdir -p "$output_dir"
                {
                    /home/project/11003054/e1101919/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools cis \
                        --vcf /home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/00.data/vcf/${timepoint}/reheader/${timepoint}_merge.vcf.gz \
                        --bed /home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/07.splicing/prepare_phenotype/update_geneInfo/final/sorted/${timepoint}_chr${i}.bed.gz \
                        --cov "$cov_folder" \
                        --permute 1000 --normal \
                        --out "${output_dir}/${timepoint}_perm_chr${i}.txt"
                }
            done
        done
    done
done


# Organize
for cov_type in "genotype_test" "phenotype_test"; do
    for cov_folder in /home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/cov/phenotype_test/${cov_type}/cov*; do
        cov_name=$(basename "$cov_folder")
        output_dir="${output_base_dir}/${cov_type}_${cov_name}"
        for timepoint in "${dir[@]}"; do
            cat "${output_dir}/${timepoint}"/*_perm_chr*.txt > "${output_dir}/${timepoint}_full_variants.txt"
            gzip "${output_dir}/${timepoint}_full_variants.txt"
        done
    done
done

# Run R script for FDR
for timepoint in "${dir[@]}"; do
    for cov_type in "genotype_test" "phenotype_test"; do
        for cov_folder in /home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/cov/phenotype_test/${cov_type}/cov*; do
            cov_name=$(basename "$cov_folder")
            output_dir="${output_base_dir}/${cov_type}_${cov_name}"
            Rscript /home/project/11003054/e1101919/software/QTLtools_1.2_Ubuntu16.04_x86_64/script/runFDR_cis.R \
                "${output_dir}/${timepoint}_full_variants.txt.gz" \
		0.05 \
                "${output_dir}/${timepoint}_permutations_afterFDR.txt"
	    	
        done
    done
done
