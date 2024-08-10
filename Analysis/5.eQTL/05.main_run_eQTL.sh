# 5.eQTL mapping - 05.main_run_eQTL.sh

in_dir=/home/project/11003054/e1101919/muscle_QTL/RNAseq/12.rasqual/input_filt/new/Pre
expr_dir=/home/project/11003054/e1101919/muscle_QTL/RNAseq/12.rasqual/rasqualtools/Pre/new
geno_dir=/home/project/11003054/e1101919/muscle_QTL/RNAseq/12.rasqual/vcf_by_chrom/Pre
out_dir=/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/joint_output1/Pre
log_dir=/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/joint_output1/log
cov_dir=/home/project/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/cov/phenotype_test

for i in {1..22}; do
    (
        echo INFO - chr$i
        n_genes=$(wc -l "$in_dir/rasqual.input.chr$i.filt.txt" | cut -d" " -f1)
        echo INFO - $n_genes genes.

        if [[ ! -d "$out_dir/chr$i/" ]]; then mkdir -p "$out_dir/chr$i/"; fi
        if [[ ! -d "$log_dir/Pre/chr$i/" ]]; then mkdir -p "$log_dir/Pre/chr$i/"; fi

        n=0
        for j in $(seq $n_genes); do
                n=$((n+1))
                if [[ $n -gt 500 ]]; then wait; n=0; fi
                bash /home/project/11003054/e1101919/muscle_QTL/RNAseq/12.rasqual/rasqual/rasqual_pre_test2.sh \
                    $in_dir/rasqual.input.chr$i.filt.txt \
                    $j \
                    $expr_dir/sams2_pre.expression.bin \
                    $expr_dir/sams2_pre.size_factors.bin \
                    $geno_dir/chr$i_modified.vcf.gz \
                    joint \
                    $out_dir/chr$i \
                    $cov_dir/pre_phenotype_cov1.bin &> $log_dir/Pre/chr$i/rasqual.chr$i.$j.log &
        done
        wait
    )
done

