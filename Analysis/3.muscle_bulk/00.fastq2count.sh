#!/bin/bash
# Program: This is a shell script for calculate read count.
# Date: 20230225
# Input: Paired-end raw fastq/fastq.gz data
# Output: Read count matrix

# Muscle bulk - 00.fastq2count
# set working directory path
work_dir=/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq

# output prefix
output=SAMS2

# set genome information
genome=$work_dir/db/GRCh38.primary_assembly.genome.fa
gff=$work_dir/db/gencode.v43.annotation.gff3
gtf=$work_dir/db/gencode.v43.annotation.gtf

index=$work_dir/db/genome

# set samples info
sampleInfo=$work_dir/00.data/samples.txt

# set thread number
thread=30

# build index
cd $work_dir/db
# hisat2-build -p ${thread} ${genome} ${index}
gffread ${gff} -T -o ${gtf}
mkdir -p $work_dir/db/STARindex
STAR --runThreadN $thread --runMode genomeGenerate --genomeDir ./STARindex2 --genomeFastaFiles $genome --sjdbGTFfile $gtf --sjdbOverhang 99

cat ${sampleInfo} | while read group sample fq1 fq2
do
	# filter raw fastq data
	cd ${work_dir}/00.data/01.clean_data
	fastp -i ${fq1} -o ./${sample}_1.clean.fastq.gz \
		-I ${fq2} -O ./${sample}_2.clean.fastq.gz \
		--json=./${sample}.json --html=${sample}.html --report_title="${sample} fastp report" \
		--thread=${thread}

	cd ${work_dir}/01.Mapping
	# mapping clean read, default is non strand-specific library, paired-end read
#	hisat2 --new-summary -p ${thread} \
#		-x ${index} \
#		-1 ../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
#		-2 ../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
#		-S ${sample}.sam \
#		1> ${sample}.log 2>&1
#	samtools sort -@ ${thread} -o ${sample}.sort.bam ${sample}.sam
#	samtools index ${sample}.sort.bam
#	rm ${sample}.sam
	 STAR --twopassMode Basic \
               --quantMode TranscriptomeSAM GeneCounts \
               --runThreadN $thread --genomeDir ../db/STARindex \
               --outSAMtype BAM SortedByCoordinate \
               --sjdbOverhang 99 \
               --outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA \
               --outSJfilterReads Unique \
               --outSAMmultNmax 1 \
               --outFileNamePrefix $sample \
               --outSAMmapqUnique 60 \
               --readFilesCommand gunzip -c \
#              --outSAMattributes NH HI AS nM NM MD jM jI MC ch vA vW vG \
#              --waspOutputMode SAMtag --varVCFfile $vcfFile \
               --readFilesIn ../00.data/01.clean_data/${sample}_1.clean.fastq.gz ../00.data/01.clean_data/${sample}_2.clean.fastq.gz
#done
done

cd ${work_dir}/02.Quantification
cat ${sampleInfo} | while read group sample fq1 fq2
do
	ln -s ../01.Mapping/$sample.sort.bam ./$sample
done

bamfile=`cut -f2 $sampleInfo`
# count fragments (read pair) using featureCounts, default is non strand-specific library, paired-end read. Multi-mapping reads will not be counted.
featureCounts -t exon -g gene_id -s 0 -p --countReadPairs -T $thread -a $gtf -o ${output}.txt $bamfile
grep -v "^# Program:featureCounts" ${output}.txt | cut -f1,2,3,4,5,6 > ${output}.geneInfo.txt
grep -v "^# Program:featureCounts" ${output}.txt | cut -f2,3,4,5,6 --complement > ${output}.read_count.txt
rm $bamfile
