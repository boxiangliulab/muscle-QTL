# 5.eQTL mapping - 06.rasqual_filter.R

library(data.table)
library(cowplot)
library(R.utils)

pre <- fread("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/pre_rasqual_output_full.txt")
post <- fread("/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/post_rasqual_output_full.txt")

phi_ub = 0.75; phi_lb = 0.25
delta_ub = 0.1
r2_fsnp_lb = r2_rsnp_lb = 0.9

pre_filter <- pre[phi < phi_ub & phi > phi_lb & 
	delta < delta_ub & r2_fsnp > r2_fsnp_lb & 
	r2_rsnp > r2_rsnp_lb, ]

post_filter <- post[phi < phi_ub & phi > phi_lb &
        delta < delta_ub & r2_fsnp > r2_fsnp_lb &
        r2_rsnp > r2_rsnp_lb, ]

write.table(pre_filter,"/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/pre_filter.txt",col.names=T,row.names=F,sep="\t",quote=F)

write.table(post_filter,"/home/users/nus/e1101919/scratch/muscle_QTL/RNAseq/rasqual_output/check/rasqual_filter/post_filter.txt",col.names=T,row.names=F,sep="\t",quote=F)
