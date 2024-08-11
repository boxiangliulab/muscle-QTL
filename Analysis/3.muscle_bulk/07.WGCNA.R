# install.packages("BiocManager") 
# BiocManager::install("WGCNA")

################################################################WGCNA 

library("dynamicTreeCut")
library("fastcluster")
library("WGCNA")
options(stringsAsFactors = FALSE)
allowWGCNAThreads(18)
# expro=read.csv("corrected_TPM_bulk.txt",sep = '\t', row.names = 1)
expo = read.csv("selected_TPM_bulk.txt",sep='\t',row.names=1)

## Check the data structure：
dim(expo)
names(expo)
datExpr = as.data.frame(t(expo))


##check missing value
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

if(!gsg$allOK)
{
	#Optionally,print the gene and sample names that were removed;
	if(sum(!gsg$goodGenes) >0)
		printFlush(paste("Removing genes:",paste(names(datExpr)[!gsg$goodGenes],collapse = ", ")))
	if(sum(!gsg$goodSamples)>0)
		printFlush(paste("Removing samples:",paste(rownames(datExpr)[!gsg$goodSamples],collapse = ", ")))
		#remove the offending genes and samples from the data
	datExpr = datExpr[gsg$goodSamples,gsg$goodGenes]
}

## filter low expression genes
meanFPKM=0.5
n=nrow(datExpr)
datExpr[n+1,]=apply(datExpr[c(1:nrow(datExpr)),],2,mean)
datExpr=datExpr[1:n,datExpr[n+1,] > meanFPKM]

filtered_fpkm=t(datExpr)

filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
write.table(filtered_fpkm,file="FPKM_filter.xls",row.names=F, col.names=T, quote=FALSE, sep='\t')

############################sample cluster
sampleTree = hclust(dist(datExpr),method = "average")
#Plot the sample tree: Open a graphic output window of size 12 by 9 inches
#The user should change the dimensions if the window is too large or too small

#sizeGrWindow(12,9)
pdf(file = "1_sampleClustering.pdf",width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main ="sample clusterint to detect outliners", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main =2)


###Plot a line to show the cut
## abline(h = 15, col = "red") 

dev.off()
traitData = read.csv("trits.txt",sep = '\t',row.names = 1) 
dim(traitData)
names(traitData)
 
#remove colnames that hold information we do not need
alltraits = traitData
dim(alltraits)
names(alltraits)

#form a data fram analogous to expression data that will hold the clinical traites
fpkmSamples = rownames(datExpr)
traitsamples= rownames(alltraits)
traitRows = match(fpkmSamples,traitsamples)
dataTraits = alltraits[traitRows, ]
rownames(dataTraits)
collectGarbage()

#Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(dataTraits, signed = FALSE)
#sizeGrWindow(12,12)
pdf(file = "2_sample dendrogram and trait heatmap.pdf",width = 15, height = 12)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(dataTraits), main = "Sample dendrogram and trait heatmap")
dev.off()
save(datExpr,file = "fpkm_forAnalysis.RData")
save(dataTraits,file = "trait_forAnalysis.RData")

#########################################network constr######
enableWGCNAThreads()

#choose a set of soft-thresholding power
power = c(1:30)

#call the network topology analysis funciton
sft = pickSoftThreshold(datExpr, powerVector = power,verbose = 5)

#plot the result
# sizeGrWindow(9,5)
pdf(file = "3_Scale independence.pdf",width = 9, height = 5)
par(mfrow = c(1,2))
cexl = 0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=power,cex=cexl,col="red");
#this line corresponds to using an R^2 cut-off of h 
abline(h=0.8,col="red")
#Mean connectivity as a funciton of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cexl,col="red")
dev.off()

####choose the softPower
softPower = sft$powerEstimate
adjacency = adjacency(datExpr, power = softPower)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
#####Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);


dissTOM = 1-TOM

#call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM),method="average");
#plot the resulting clustering tree(dendrogram)

# sizeGrWindow(12,9)
pdf(file = "4_Gene clustering on TOM-based dissimilarity.pdf",width = 12, height = 9)
plot(geneTree, xlab = "", sub="",main = "Gene clustering on TOM-based dissimilarity", labels=FALSE, hang = 0.04)
dev.off()

#we like large modules, so we set the mininum module size relatively high:
minModulesize=10
#module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModulesize);
table(dynamicMods)

#convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
pdf(file = "5_Dynamic Tree cut.pdf",width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors,"Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,  main = "Gene dendrogram and module colors")
dev.off()

# calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs)
#Cluster module eigengenes
METree = hclust(as.dist(MEDiss),method = "average")
#Plot the result
# sizeGrWindow(7,6)
pdf(file = "6_Clustering of module eigengenes.pdf",width = 7, height = 6)
plot(METree, main = "clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.25
#plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

#call an automatic merging funciton
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres,verbose = 3)
#the merged module colors
mergedColors = merge$colors
#Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# sizeGrWindow(12,9)
pdf(file = "7_merged dynamic.pdf",width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut","Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
 
#Rename to moduleColors
moduleColors = mergedColors
#construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, TOM, dissTOM, moduleColors, geneTree, sft, file="networkConstruction-stepbystep.RData")
##########################################relate module to external clinical triats ##################
#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleTraitCor = cor(MEs, dataTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# sizeGrWindow(10,6)
pdf(file = "8_Module-trait relationships.pdf",width = 5, height = 6)
# will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor,2),"\n(", signif(moduleTraitPvalue,1),")",sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6,8.5,3,3))

#Dispaly the correlation value within a heatmap plot			   
labeledHeatmap(Matrix = moduleTraitCor, 
                xLabels = names(dataTraits), yLabels = names(MEs), cex.lab = 0.9,  yColorWidth=0.01, xColorWidth = 0.03,
              ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()	
###################Define variable weight contatining all column of datTraits
 
##MM and GS
 
#names(colors) of the modules
modNames = substring(names(MEs),3)

geneModuleMembership = as.data.frame(cor(datExpr,MEs,use = "p"))
MMpvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMpvalue) = paste("p,MM", modNames, sep = "")

#names of those trait
traitNames = names(dataTraits)

geneTraitSignificance = as.data.frame(cor(datExpr,dataTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep = "")
names(GSPvalue) = paste("p.GS.", traitNames, sep = "")
 names(datExpr)
probes = names(datExpr)

################################export GS and MM

geneInfo0 = data.frame(probes = probes, moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
	oldNames = names(geneInfo0)
	geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra], GSPvalue[,Tra])
	names(geneInfo0) = c(oldNames, names(geneTraitSignificance)[Tra], names(GSPvalue)[Tra])
}
for(mod in 1:ncol(geneModuleMembership))
{
	oldNames = names(geneInfo0)
	geneInfo0 = data.frame(geneInfo0,geneModuleMembership[,mod], MMpvalue[,mod])
	names(geneInfo0) = c(oldNames, names(geneModuleMembership)[mod], names(MMpvalue)[mod])
}
geneorder = order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneorder, ]
write.table(geneInfo,file="10_GS_and_MM.xls",sep="\t",row.names=F)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
#set diagonal to NA for a nicer plot
diag(plotTOM) = NA

#Call the plot function 
#sizeGrWindow(9,9)
pdf(file="12_Network heatmap plot_all gene",width = 9, height = 9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot. all genes")
dev.off()

nSelect = 400
#For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select,select]
selectTree = hclust(as.dist(selectTOM),method = "average")
selectColors = moduleColors[select]

#Open a graphical window
#sizeGrWindow(9,9)
#Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
plotDiss = selectTOM^7
diag(plotDiss) = NA

pdf(file="13_Network heatmap plot_selected genes.pdf",width = 9, height = 9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot. selected genes")
dev.off()
###############################exporting to cytoscape all one by one########################################


#select each module
for (mod in 1:nrow(table(moduleColors)))
{
	modules = names(table(moduleColors))[mod]
	#Select module probes 
	probes = names(datExpr)
	inModule = (moduleColors == modules)
	modProbes = probes[inModule]
	modGenes = modProbes
	#Select the corrsponding Topologyical overlap 
	modTOM = TOM[inModule,inModule]
	
	dimnames(modTOM) = list(modProbes, modProbes)
	#Export the network into edge and node list file cytoscape can read
	cyt = exportNetworkToCytoscape(modTOM,
									edgeFile = paste("CytoscapeInput-edges-",modules,".txt",sep = ""),
									nodeFile = paste("CytoscapeInput-nodes-",modules,".txt",sep = ""),
									weighted = TRUE, threshold = 0.02, nodeNames =  modProbes,
									altNodeNames = modGenes,
									nodeAttr = moduleColors[inModule])
}
