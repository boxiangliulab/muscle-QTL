# 5. eQTL mapping - 08.treeQTL_multi_tissue.R

library(data.table)
library(stringr)
library(cowplot)
library(dplyr)
library(TreeQTL)

rasqual_fn_pre='/data/projects/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/rasqual_output/check/rasqual_filter/pre_filter_all_association.txt'
rasqual_fn_post='/data/projects/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/rasqual_output/check/rasqual_filter/post_filter_all_association.txt'

out_dir_pre='/data/projects/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/rasqual_output/TreeQTL/multiTreeQTL_eGene_2try/pre'
out_dir_post='/data/projects/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/rasqual_output/TreeQTL/multiTreeQTL_eGene_2try/post'

out_dir='/data/projects/11003054/e1101919/muscle_QTL/RNAseq/09.eQTL/rasqual_output/TreeQTL/multiTreeQTL_eGene_2try'

joint_pre=fread(rasqual_fn_pre, select = c(1,2,10,11,12), col.names=c('fid','sid','log10qval','chisq','pi'))
joint_post=fread(rasqual_fn_post, select = c(1,2,10,11,12), col.names=c('fid','sid','log10qval','chisq','pi'))

get_gene_map=function(gene_ids){
        gene_ids=unique(gene_ids)
        gencode=read.table(text = paste(readLines("/home/project/11003054/e1101919/muscle_QTL/RNAseq/db/gencode.v43.primary_assembly.annotation.gtf.gz")[which(!grepl("^##", readLines("/home/project/11003054/e1101919/muscle_QTL/RNAseq/db/gencode.v43.primary_assembly.annotation.gtf.gz")))[1]:length(readLines("/home/project/11003054/e1101919/muscle_QTL/RNAseq/db/gencode.v43.primary_assembly.annotation.gtf.gz"))], collapse = "\n"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        gencode$geneid <- sub('.*gene_id (.*?);.*', '\\1', gencode$V9)
        gene_map <- gencode[gencode$V3 == 'gene', ]
        gene_map <- gene_map[,c(10,1,4,5)]
        colnames(gene_map)[1] <- "geneid"
        gene_map$geneid <- sub("\\..*", "", gene_map$geneid)
        gene_map=gene_map[gene_map$geneid%in%gene_ids,]
        return(gene_map)
}

get_snp_map=function(snps){
        snp_map=data.table(snpid=unique(snps))
        snp_map[,chr:=str_split_fixed(snpid,'_',5)[,1]]
        snp_map[,chr:=paste0('chr',chr)]
        snp_map$chr <- sub("^(.*?):.*", "\\1", snp_map$chr)
        snp_map[,pos:=as.integer(str_split_fixed(snpid,'_',5)[,2])]
        snp_map$pos <- as.integer(sub("^[^:]+:(\\d+):.*", "\\1", snp_map$snpid))
        return(snp_map)
}

joint_pre[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
joint_post[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]

setorder(joint_pre,pval)
setorder(joint_post,pval)

meqtl_pre=joint_pre[,list(SNP=sid,gene=fid,beta=pi,`t-stat`=-1,`p-value`=pval,FDR=-1)]
meqtl_post=joint_post[,list(SNP=sid,gene=fid,beta=pi,`t-stat`=-1,`p-value`=pval,FDR=-1)]

gene_map=get_gene_map(unique(c(joint_pre$fid, joint_post$fid)))
snp_map=get_snp_map(unique(c(joint_pre$sid, joint_post$sid)))

MTtreeQTL=function(meqtl_glucose, meqtl_galactose, snp_map, gene_map, level1=0.05,level2=0.05,level3=0.05,eSNP=FALSE, out_dir,cis_dist=1e6) {
	fwrite(meqtl_pre,sprintf('%s/mEQTL_pre_out_cis.txt',out_dir),sep='\t')
	fwrite(meqtl_post,sprintf('%s/mEQTL_post_out_cis.txt',out_dir),sep='\t')
	snps_by_tissue=data.frame(snp_name=snp_map$snpid, glucose=rep(1,nrow(snp_map)), galactose=rep(1,nrow(snp_map)))
	genes_by_tissue=data.frame(gene_name=gene_map$geneid, glucose=rep(1,nrow(gene_map)), galactose=rep(1,nrow(gene_map)))
	if (eSNP){
		print('INFO - level1 is eSNP')
		print('INFO - calculating number of tests per SNP...')
		n_tests_per_SNP=get_n_tests_per_SNP(snp_map,gene_map,nearby=TRUE,dist=1e+06)
		print('INFO - getting eSNPs...')
		eSNPs=get_eSNPs_multi_tissue(genes_by_tissue=genes_by_tissue, snps_by_tissue=snps_by_tissue, n_tests_per_SNP=n_tests_per_SNP, m_eqtl_out_dir=out_dir, tissue_names=c("post","pre"), level1 = 0.05, level2 = 0.05, level3 = 0.05)
		return(eSNPs)
		} else {
			print('INFO - level1 is eGene')
			print('INFO - getting eGenes and eAssociations...')
			eGenes=get_eGenes_multi_tissue(genes_by_tissue = genes_by_tissue, snps_by_tissue = snps_by_tissue, snp_map = snp_map, gene_map = gene_map, nearby = TRUE, dist = 1e+06, m_eqtl_out_dir=out_dir,  tissue_names=c("post","pre"), level1 = 0.05, level2 = 0.05, level3 = 0.05)
			return(eGenes)
	}
}

eGenes=MTtreeQTL(meqtl_pre, meqtl_post, snp_map, gene_map, level1=0.05,level2=0.05,level3=0.05, eSNP=FALSE,out_dir,cis_dist=1e6)
fwrite(eGenes,sprintf('%s/eGenesMT.txt',out_dir),sep='\t')
table(eGenes[,2:3])


