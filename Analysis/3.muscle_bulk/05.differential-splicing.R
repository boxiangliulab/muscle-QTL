#!/bin/bash

# Set paths and filenames
junction_files="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/07.splicing/script/test_juncfiles.txt"
output_prefix="SAMS2"
min_intron_length=50
max_intron_length=500000
exonbed="/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/07.splicing/leafcutter_ds/exonbed_hg38.txt"
num_threads=4
min_samples_per_intron=4
counts_file="${output_prefix}_perind_numers.counts.gz"
group_file="leafcutter_ds_group.txt"

# Step 1: Run LeafCutter Cluster
echo "Running LeafCutter Cluster..."
python /home/users/nus/e1101919/scratch/software/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $junction_files -m $min_intron_length -o $output_prefix -l $max_intron_length

# Step 2: Check exon file content
echo "Checking exon file..."
head $exonbed

# Step 3: Run LeafCutter differential splicing analysis
echo "Running LeafCutter DS..."
/home/users/nus/e1101919/scratch/software/leafcutter/scripts/leafcutter_ds.R --num_threads $num_threads \
  --exon_file=$exonbed --min_samples_per_intron=$min_samples_per_intron \
  $counts_file $group_file

echo "LeafCutter analysis complete."

