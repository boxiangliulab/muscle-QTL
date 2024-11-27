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
