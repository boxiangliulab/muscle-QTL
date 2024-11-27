exonbed=/data/projects/11003054/e1101919/muscle_QTL/RNAseq/07.splicing/leafcutter_ds/exonbed_hg38.txt
head $exonbed
/data/projects/11003054/e1101919/software/leafcutter/scripts/leafcutter_ds.R --num_threads 4 \
--exon_file=$exonbed  --min_samples_per_intron=4 \
SAMS2_perind_numers.counts.gz leafcutter_ds_group.txt
