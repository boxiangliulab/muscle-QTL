# 6.sQTL mapping - 01.nominal_pass.sh

dir=("Pre" "Post")

output_base_dir="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/08.sQTL/nominal"

#Generate
for timepoint in "${dir[@]}"; do
    for i in $(seq 1 22); do
            for cov_folder in /home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/08.sQTL/cov/phenotype_test/cov*; do
                cov_name=$(basename "$cov_folder")
                output_dir="${output_base_dir}/${cov_type}_${cov_name}"
                mkdir -p "$output_dir"
                {
                    /home/project/11003054/e1101919/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools cis \
                        --vcf /home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/00.data/vcf/${timepoint}/${timepoint}_merge.vcf.gz \
                        --bed /home/project/11003054/e1101919/muscle_QTL/RNAseq/07.splicing/prepare_phenotype/update_geneInfo/final/sorted/${timepoint}_chr${i}.bed.gz \
                        --cov "$cov_folder" \
                        --nominal 1 --normal \
                        --out "${output_dir}/${timepoint}_nominals_chr${i}.txt"
                }
            done
        done
    done

# Organize
for cov_folder in /home/project/11003054/e1101919/muscle_QTL/RNAseq/08.sQTL/cov/${cov_type}/cov*; do
        cov_name=$(basename "$cov_folder")
        output_dir="${output_base_dir}/${cov_type}_${cov_name}"
	rm "${output_dir}"/*gz
        for timepoint in "${dir[@]}"; do
                cat "${output_dir}"/${timepoint}_nominals_chr*.txt | awk '{if($14==1) print $0}' > "${output_dir}/${timepoint}_top_variants.txt"
                gzip "${output_dir}/${timepoint}_top_variants.txt"
	done
done
