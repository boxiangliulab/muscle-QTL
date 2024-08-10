# 5.eQTL mapping - 02.peer.R

library(peer)
library(data.table)

expr <- read.delim("/home/project/11003054/e1101919/muscle_QTL/RNAseq/10.peer/Pre_TPM.txt")
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setNk(model,15)
PEER_getNk(model)
PEER_update(model)
factors = as.data.frame(t(PEER_getX(model)))
colnames(factors) <- colnames(expr)
write.table(t(factors),"/home/project/11003054/e1101919/muscle_QTL/RNAseq/10.peer/new_pre_peer_covariates_15.txt",quote=F,row.names=T,sep="\t",col.names=T)

expr <- read.delim("/home/project/11003054/e1101919/muscle_QTL/RNAseq/10.peer/Post_TPM.txt")
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setNk(model,15)
PEER_getNk(model)
PEER_update(model)
factors = as.data.frame(t(PEER_getX(model)))
colnames(factors) <- colnames(expr)
write.table(t(factors),"/home/project/11003054/e1101919/muscle_QTL/RNAseq/10.peer/new_post_peer_covariates_15.txt",quote=F,row.names=T,sep="\t",col.names=T)

