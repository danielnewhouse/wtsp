source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
install.packages("WGCNA") 
library("WGCNA")
library("DESeq2")
#install.packages("flashClust")
library("flashClust")
#install.packages("pheatmap")
library("pheatmap")
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

################################
# DEseq2 to get vsd for WGCNA  #
################################


countData <- read.table ("star_htseq_all_geneIDs.txt", header=TRUE)
head(countData)

wtspDesign = data.frame (
  row.names = colnames(countData),
  morph = c("Tan","White","White","White","White","Tan","Tan","Tan","White","White","White","Tan","Tan","Tan","Tan","White","Tan","White","White","White","White","Tan","Tan","Tan","Tan","Tan","White","White","White","White","White","Tan"),
  nest = c("TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","WxT","WxT","WxT","WxT","TxW","TxW","WxT","WxT","WxT","TxW","WxT","WxT","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW"),
  sex = c("F","F","F","F","M","M","M","F","F","F","F","M","M","F","M","F","F","M","F","M","F","F","F","F","M","M","M","M","F","M","F","U"))

ddsAll<- DESeqDataSetFromMatrix (countData = countData,
                                 colData = wtspDesign,
                                 design = ~ nested.nest + sex + morph + parentalcare)

ddsAll <- DESeq(ddsAll)

vsd=getVarianceStabilizedData(ddsAll)
head(vsd)

write.csv(as.data.frame(vsd), file = "WTSP_WGCNA_genes_greaterthan5.csv")

###### Genes used for WGCNA were first filtered prior to DEseq2 - only keeping genes with >5 counts.
###### Output vsd of those genes and again only kept genes >5 to use as input for WGCNA.



################################
#            WGCNA             #
################################


#-----Load expression data
dis0= read.csv("WTSP_WGCNA_genes_greaterthan5_orderedNestType.csv")
head(dis0)
names(dis0)
dis=dis0[c(1:31)]
dim(dis)
names(dis)
rownames(dis) <- dis$X
datExpr0= as.data.frame(t(dis[, -c(1)]))
names(datExpr0)= dis$X
rownames(datExpr0)=names(dis)[-c(1)]
dim(datExpr0)

gsg=goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts


#-----Load trait data
# removed Age from traits. it's corrolated with one module and that module is not correlated with any other trait, i.e. age doesn't influence other traits.

traitData= read.csv("WTSP.traits.noindividual_outliersremoved_condensed_orderednesttype_noage.csv")
#traitData= read.csv("WTSP.traits.noindividual_outliersremoved_condensed_orderednesttype_age_nestsize.csv")
dim(traitData)
head(traitData)
names(traitData)

rownames(traitData) <- traitData$Sample
traitData$Sample <- NULL
datTraits= traitData


#######   #################    ################   #######    
#                 Call sample outliers
#######   #################    ################   #######   

#-----Sample dendrogram and traits
A=adjacency(t(datExpr0),type="signed")
#-----Calculate whole network connectivity
k=as.numeric(apply(A,2,sum))-1
#-----Standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
#-----Convert traits to colors
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlier=outlierColor,traitColors)

#-----Plot the sample dendrogram
quartz()
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")

# no outliers here. Note, two samples were identified as outliers previously and removed - 25896, 25898

#-----Remove outlying samples 
#remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
#datExpr0=datExpr0[!remove.samples,]
#datTraits=datTraits[!remove.samples,]
#A=adjacency(t(datExpr0),type="distance")
#k=as.numeric(apply(A,2,sum))-1
#Z.k=scale(k)

save(datExpr0, datTraits, file="SamplesAndTraits2.RData")

#######   #################    ################   #######    
#                     Choose soft threshold
#######   #################    ################   #######     
options(stringsAsFactors = FALSE)
lnames= load(file="SamplesAndTraits2.RData")
lnames
dim(datExpr0)
dim(datTraits)
powers= c(seq(1,10,by=0.5), seq(from =8, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose =5,networkType="signed") #call network topology analysis function

quartz()
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

softPower=12
### see manual for more info here. should pick the first number as the curve falttens out in the graph generated above (Scale Free Topology)


#######   #################    ################   #######    
#                    Construct network
#######   #################    ################   #######     

adjacency=adjacency(datExpr0, power=softPower, type="signed") 
TOM= TOMsimilarity(adjacency, TOMType="signed")
dissTOM= 1-TOM

geneTree= flashClust(as.dist(dissTOM), method="average")

quartz()
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#######   #################    ################   #######    
#                    Make modules
#######   #################    ################   ####### 

minModuleSize=30
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

dynamicColors= labels2colors(dynamicMods)

quartz()
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#-----Merge modules whose expression profiles are very similar
MEList= moduleEigengenes(datExpr0, colors= dynamicColors)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

quartz()
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres= 0.20
abline(h=MEDissThres, col="red")
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

# MinSize 30, MEDisThres-> 0.10=27 modules, 0.20=26 modules, 0.30=22 modules, 0.35=19 modules
# MinSize 50, MEDisThres-> 0.10=19 modules, 0.20=19 modules, 0.30=17 modules, 0.35=16 modules

mergedColors= merge$colors
mergedMEs= merge$newMEs

quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)


moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs


save(MEs, moduleLabels, moduleColors, geneTree, file= "SamplesAndColors_thresh12merge20_signed2.RData")


#######   #################    ################   #######    
#                Relate modules to traits
#######   #################    ################   ####### 

datt=datExpr0

#-----Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
#-----Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#-----Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#---------------------Module-trait heatmap with only sig rows shown

# module-trait correlations #1 WITH ONLY SIGNIFICANT CORRELATIONS
MEcor <- as.data.frame(matrix("NA", nrow = 26, ncol = 6), row.names = rownames(moduleTraitPvalue))
colnames(MEcor) <- colnames(moduleTraitCor)

# only show significant modules
for (i in 1:26){
  for (j in 1:6){
    MEcor[i,j] <- as.character(ifelse(moduleTraitPvalue[i,j] < 0.05, signif(moduleTraitCor[i,j], 2), as.character("-")))
  }
}

quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 20, 3, 3));
moduleTraitCor[moduleTraitCor > -0.36 & moduleTraitCor < 0.36]<-0     #make all non-significant modules white
contrasting4 = colorRampPalette(rev(c("chocolate1","chocolate1","white","white", "steelblue1","steelblue")))(100)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = contrasting4,
               textMatrix = MEcor, #textMatrix for all pvalues, MEcor for only sig
               setStdMargins = FALSE,
               cex.text = 0.8,  #1.5 to make text in cells bigger
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))






######## output module membership
datME=moduleEigengenes(datt,mergedColors)$eigengenes
datKME=signedKME(datt, datME, outputColumnName="MM.")
genes=names(datt)
geneInfo0 = data.frame(gene=genes,moduleColor=moduleColors)
color=data.frame(geneInfo0,datKME) #these are from your original WGCNA analysis 
head(color)
write.csv(as.data.frame(color), file = "WTSP_Blood_WGCNA_GeneModule_membership_outliers_removed.csv")

#### MM pvalues
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
pvals=data.frame(geneModuleMembership,MMPvalue)
head(pvals)
write.csv(as.data.frame(pvals), file = "WTSP_Blood_WGCNA_GeneModule_membership_outliers_removed_pvalues.csv")





#---------------------Gene significance by Module membership 

#scatterplots of trait-module correlations

whichTrait="PairType" #Replace this with the trait of interest

quartz()
nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(datTraits[,whichTrait]);
names(selTrait) = whichTrait
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
par(mfrow=c(2,3))
counter=0
for(module in modNames[1:length(modNames)]){
  counter=counter+1
  if (counter>6) {
    quartz()
    par(mfrow=c(2,3))
    counter=1
  }
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste(module,"module membership"),
                     ylab = paste("GS for", whichTrait),
                     col = module,mgp=c(2.3,1,0))
}
######--------------------end--------------------#######


#---------------------Eigengene heatmap
which.module="tan" #replace with module of interest
datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
#       nrgcols=30,rlabels=F,rcols=which.module,
#      main=which.module, cex.main=2)
#par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", names.arg=c(row.names(datt)), cex.names=0.5, cex.main=2,
        ylab="eigengene expression",xlab="sample")
######--------------------end--------------------#######






#######   #################    ################   #######    
#             Gene expression within modules
#######   #################    ################   ####### 



#---------------------Heatmap for top-ModuleMembership genes in a module
# based on MM score, not gene significance score. These top genes are the "hubs"



vsd=read.csv("WTSP_WGCNA_genes_greaterthan5_orderedMorph.csv") #place original vsd file  or individual module vsd stats
names(vsd)

a.vsd=vsd[c(2:31)] #Columns with vsd
row.names(a.vsd)=vsd$X
head(a.vsd)
names(a.vsd)

allkME =as.data.frame(signedKME(t(a.vsd), MEs))

gg=read.table("WTSP_Blood_WGCNA_GeneModule_membership_outliers_removed.csv", sep="\t", header = T)  # place module membership file or individual module
head(gg)

whichModule="lightcyan"
top=183

modcol=paste("kME",whichModule,sep="")
sorted=a.vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

wgcnaDataWTSP <- read.table ("WTSP_WGCNA_genes_greaterthan5_orderedMorph.txt", header=TRUE) #tab-delim formatted vst file
wtspDesign = data.frame (
  row.names = colnames(wgcnaDataWTSP),
  #pairtype = c("BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","BiParental","FemaleBiased","FemaleBiased","FemaleBiased","FemaleBiased","FemaleBiased","FemaleBiased","FemaleBiased","FemaleBiased","FemaleBiased"))
  morph= c("White","White","White","White","White","White","White","White","White","White","White","White","White","White","White","White","White","Tan","Tan","Tan","Tan","Tan","Tan","Tan","Tan","Tan","Tan","Tan","Tan","Tan"))

annotation_colors = list(
  #parentalcare = c(BiParental = "black", FemaleBiased = "gray"))
  morph = c(White = "grey94", Tan = "tan"))

pheatmap(hubs,scale="row",cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = F, show_colnames = T, annotation_colors = annotation_colors, annotation_col = wtspDesign, col=contrasting,border_color=NA,  main=paste(whichModule,"module",sep=" "))
######--------------------end--------------------#######







#---------------------Extract genes and VSD by module
head(vsd)

col="lightgreen"

cands=names(datt[moduleColors==col])
c.vsd=vsd[vsd$X %in% cands,]
head(c.vsd)
length(c.vsd[,1])
write.csv(c.vsd,paste("vsd_",col,".csv",sep=""),quote=F, row.names=F)
######--------------------end--------------------#######














########### Get Pvalues associated with Module Membership scores. 
# correlation between MM and Gene significance 
# i.e. finding genes highly correlated with module AND with external trait

# Define variable time containing the time_h column of datTrait
pc = as.data.frame(datTraits$ParentalCare);
names(pc) = "pc"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, pc, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
module = "lightgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Parental Care",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 


## repeat to get scatterplots for top 2 each of positive/negative correlations
module = "grey60"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for time",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 



### output GS with MM scores
annot = read.csv(file = "WTSP_WGCNA_genes_greaterthan5_orderedMorph.csv")
dim(annot)
names(annot)
probes = names(datExpr0)
names(probes)
probes2annot = match(probes, annot)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for trait of interest
modOrder = order(-abs(cor(MEs, pc, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
} 
write.csv(geneInfo0, file = " geneInfo0.csv", row.names = FALSE) 










############## VisAnt
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 12);
# Read in the annotation file
annot = read.csv(file = "WTSP_Blood_WGCNA_GeneModule_membership_outliers_removed.csv");
# Select module
module = "pink";
# Select module probes
probes = names(datExpr0)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)


# ntop is number of hub genes to restrict network to. i.e. Top 30 genes based on MM or IM (should be the same)
# within visant, you'll plot the ### strongest connections based on TOM scores from the interactions of these 30 genes
# if module is small, probably don't have to filter





nTop = 300;
IMConn = softConnectivity(datExpr0[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top300.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0 )



topGenesblue<-data.frame(IMConn,modProbes)[order(-IMConn),]
topGenesblue<-topGenesblue[1:10,]
names(topGenesblue)<-c("IMConnectivity","Genes")
topGenesblue




