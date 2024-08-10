# 5.eQTL mapping - 00.genotype_add_AS_counts.sh

#!/bin/bash
module load anaconda3
conda activate muscle
module load gsl
export RASQUALDIR=/home/users/nus/e1101919/scratch/software/rasqual

# 定义时间点列表
timepoints=("Pre" "Post")

# 遍历时间点列表
for timepoint in "${timepoints[@]}"
do
    # 构建输入文件名和输出文件名
    timepoint_bam="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/12.rasqual/${timepoint}_bam.txt"
    timepoint_SAMS2_merge="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/12.rasqual/${timepoint}_SAMS2_merge.vcf.gz"
    timepoint_with_AS="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/12.rasqual/${timepoint}_with_AS.vcf.gz"

    # 执行你的命令，替换$timepoint变量
    bash /home/users/nus/e1101919/scratch/software/rasqual/src/ASVCF/createASVCF.sh paired_end "$timepoint_bam" "$timepoint_SAMS2_merge" "$timepoint_with_AS" rna
done

