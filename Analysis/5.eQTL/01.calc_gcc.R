# 5.eQTL mapping - 01.calc_gcc.R

library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

# Set file path
genome_file <- "/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa"
annotation_file <- "/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf"

# Read in the genome
genome <- FaFile(genome_file)
open(genome)

# Read in the annotation file
gtf <- import(annotation_file, format = "gtf")

# Calculategene or exon
option <- "gene"  # 可以根据需要改变为 "exon"

if (option == "gene") {
  genes <- gtf[gtf$type == "gene"]
  gene_ids <- sapply(elementMetadata(genes)$gene_id, function(x) unlist(strsplit(x, split = ";"))[1])
  
  for (i in seq_along(genes)) {
    seq <- getSeq(genome, genes[i])
    gcc <- letterFrequency(seq, letters = c("G", "C"), as.prob = TRUE)
    cat(sprintf("%s\t%f\n", gene_ids[i], sum(gcc)))
  }
  
} else if (option == "exon") {
  genes <- gtf[gtf$type == "exon"]
  gene_ids <- sapply(elementMetadata(genes)$gene_id, function(x) unlist(strsplit(x, split = ";"))[1])
  
  current_gene <- NULL
  seqs <- DNAStringSet()
  
  for (i in seq_along(genes)) {
    gene_id <- gene_ids[i]
    if (!is.null(current_gene) && current_gene != gene_id) {
      gcc <- letterFrequency(reduce(seqs), letters = c("G", "C"), as.prob = TRUE)
      cat(sprintf("%s\t%f\n", current_gene, sum(gcc)))
      seqs <- DNAStringSet()
    }
    
    current_gene <- gene_id
    seqs <- append(seqs, getSeq(genome, genes[i]))
  }
  
  if (length(seqs) > 0) {
    gcc <- letterFrequency(reduce(seqs), letters = c("G", "C"), as.prob = TRUE)
    cat(sprintf("%s\t%f\n", current_gene, sum(gcc)))
  }
  
} else {
  stop("3rd option must be either gene or exon")
}

close(genome)

