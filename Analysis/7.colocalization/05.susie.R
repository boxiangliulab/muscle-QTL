# 7. colocalization - 05.susie.R

# library(vcfR)
library(susieR)
# install.packages("remotes")
# remotes::install_github("etnite/bwardr")
library(bwardr)

args = commandArgs(trailingOnly=TRUE)
# args[1] = '/ebs1/users/liufei/project/susie_pipeline/data/processed/vcf_lead_SNPs_2Mbwindow/CD16+_Monocyte.0.vcf'
# args[2] = '/ebs1/users/liufei/project/susie_pipeline/data/processed/plink_matrix/CD16+_Monocyte.0.matrix.ld'
# args[3] = '/ebs1/users/liufei/project/susie_pipeline/data/processed/phenotype/CD16+_Monocyte.0.y.txt'

# print(args[1])
# print(args[2])
# print(args[3])


vcf<-read.csv(args[1],sep='\t',comment.char='#',header=FALSE, row.names=3)
X = gt2num(as.matrix(t(vcf)))$genomat

# R <- read.csv(args[2],header=FALSE, sep=' ')
# R <- as.matrix(R[c(1:length(R)-1),c(1:length(R)-1)])


y <- read.csv(args[2], row.names = 1, header=TRUE, sep='\t')

X <- scale(X,center = TRUE,scale = FALSE)
ss  <- univariate_regression(X,y$X0)

dat <- compute_suff_stat(X,y$X0,standardize = FALSE)
R   <- cov2cor(dat$XtX)

result = susie_rss(bhat = ss$betahat,shat = ss$sebetahat,R = R,n=length(y$X0),var_y = var(y$X0),L = 10,estimate_residual_variance = TRUE,coverage = 0.90)
# print(result$sets$cs_index)
a = data.frame(result$sets$cs$L1)
write.table(a, args[3], row.names=FALSE, col.names=FALSE,sep='\t')

