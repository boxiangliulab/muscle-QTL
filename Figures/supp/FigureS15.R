/home/project/11003054/e1101919/software/ggsashimi.py \  # Path to the ggsashimi.py script
-b /home/project/11003054/e1101919/muscle_QTL/RNAseq/07.splicing/ggsashimi/input_bams_ANK1_pre.tsv \ 
-c chr8:41652753-41664801 \  # Genomic coordinates for the plot
-g /home/project/11003054/e1101919/muscle_QTL/RNAseq/db/gencode.v44.primary_assembly.annotation.gtf \  # Path to the GTF file with gene annotations
-M 10 \  # Maximum number of reads to show in each pileup
-C 3 \  # Coverage track height
-O 3 \  # Offset between different sashimi plots (different BAM files)
-A mean \  # Method for calculating the average coverage (mean or median)
--alpha 1 \  # Alpha transparency for the sashimi arcs
-F pdf \  # Output file format
-R 350 \  # DPI for the output file
--base-size=16 \  # Base font size for plot text
--height=3 \  # Height of the plot in inches
--width=18 \  # Width of the plot in inches
--ann-height=4 \  # Height of the annotation tracks
-P /home/project/11003054/e1101919/software/ggsashimi/examples/palette.txt  # Path to the color palette file

