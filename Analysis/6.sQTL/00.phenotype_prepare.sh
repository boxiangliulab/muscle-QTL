# 6.sQTL mapping - 00.phenotype_prepare.sh

# Define directories for "Pre" and "Post" analysis
dirs=("Pre" "Post")

# Iterate over each directory to process BAM files
for dir in "${dirs[@]}"; do
    # Create a file to list junction files for leafcutter
    junc_file_list="/home/project/11003054/e1101919/muscle_QTL/RNAseq/01.wasp_Mapping/combined_sorted_bam/${dir}_juncfiles.txt"
    > $junc_file_list  # Clear the contents of the file to ensure it's empty

    # Process each BAM file in the directory
    for bamfile in /home/project/11003054/e1101919/muscle_QTL/RNAseq/01.wasp_Mapping/combined_sorted_bam/${dir}/*.bam; do
        # Extract junctions from BAM using regtools
        junction_file="${bamfile}.junc"
        regtools junctions extract -a 8 -m 50 -M 500000 -s FR $bamfile -o $junction_file
        echo $junction_file >> $junc_file_list
    done

    # Use leafcutter to cluster junctions
    python /home/project/11003054/e1101919/software/leafcutter/clustering/leafcutter_cluster_regtools.py \
        -j $junc_file_list -m 50 -o "${dir}_aida" -l 500000
done

