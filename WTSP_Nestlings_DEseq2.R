############# Load DEseq2 and data ##########


source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
biocLite("pcaExplorer")
library("pcaExplorer")


countData <- read.table ("star_htseq_all_geneIDs_no896_898_greaterthan5.txt", header=TRUE)
head(countData)


############# Pairwise DE ##########

wtspDesign = data.frame (
  row.names = colnames(countData),
  morph = c("Tan","White","White","White","White","Tan","Tan","Tan","White","White","Tan","Tan","Tan","Tan","Tan","White","White","White","White","White","Tan","White","White","White","Tan","Tan","White","Tan","White","White"),
  sex = c("F","F","F","F","M","M","M","F","M","F","F","F","F","M","M","M","M","F","M","F","M","F","F","F","M","M","F","F","F","M"),
  sexPC = c("F_TxW","F_TxW","F_TxW","F_TxW","M_TxW","M_TxW","M_TxW","F_TxW","M_TxW","F_TxW","F_TxW","F_TxW","F_TxW","M_TxW","M_TxW","M_TxW","M_TxW","F_TxW","M_TxW","F_TxW","M_TxW","F_WxT","F_WxT","F_WxT","M_WxT","M_WxT","F_WxT","F_WxT","F_WxT","M_WxT"),
  pairtype = c("TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","TxW","WxT","WxT","WxT","WxT","WxT","WxT","WxT","WxT","WxT"),
  whichnest = c("TxW_1","TxW_1","TxW_1","TxW_1","TxW_2","TxW_2","TxW_2","TxW_2","TxW_3","TxW_4","TxW_4","TxW_4","TxW_4","TxW_5","TxW_5","TxW_5","TxW_5","TxW_6","TxW_6","TxW_6","TxW_7","WxT_1","WxT_1","WxT_1","WxT_1","WxT_2","WxT_2","WxT_2","WxT_3","WxT_3"),
  nested.nest = c("1","1","1","1","2","2","2","2","3","4","4","4","4","5","5","5","5","6","6","6","7","1","1","1","1","2","2","2","3","3"),
  morphnest = c("Tan_TxW","White_TxW","White_TxW","White_TxW","White_TxW","Tan_TxW","Tan_TxW","Tan_TxW","White_TxW","White_TxW","Tan_TxW","Tan_TxW","Tan_TxW","Tan_TxW","Tan_TxW","White_TxW","White_TxW","White_TxW","White_TxW","White_TxW","Tan_TxW","White_WxT","White_WxT","White_WxT","Tan_WxT","Tan_WxT","White_WxT","Tan_WxT","White_WxT","White_WxT"),
  age = c("six","six","six","six","six","six","six","six","six","six","six","seven","seven","six","six","six","six","six","six","six","seven","seven","six","seven","seven","six","five","five","five","six"))
wtspDesign




#MORPH x NEST differences with nested nest
# White FB vs White BP
ddsmorphnest<- DESeqDataSetFromMatrix (countData = countData,
                                       colData = wtspDesign,
                                       design = ~ nested.nest + sex + morphnest)

ddsmorphnest <- DESeq(ddsmorphnest)
ddsmorphnest <- ddsmorphnest[which(mcols(ddsmorphnest)$betaConv),]
resmorphnest <- results(ddsmorphnest, contrast=c("morphnest","White_FB","White_BP"))
head(resmorphnest)
sum(resmorphnest$padj < 0.1, na.rm=TRUE) ## 622 genes
sum(resmorphnest$padj < 0.05, na.rm=TRUE) ## 344 genes
write.csv(as.data.frame(resmorphnest), file = "WTSP_Blood_White_FBvsBP.csv")


# Tan vs White in FB nests
resMorphNest_Contrast2 <- results(ddsmorphnest, contrast=c("morphnest","Tan_FB","White_FB"))
head(resMorphNest_Contrast2)
sum(resMorphNest_Contrast2$padj < 0.1, na.rm=TRUE) ## 40
sum(resMorphNest_Contrast2$padj < 0.05, na.rm=TRUE) ## 18
write.csv(as.data.frame(resMorphNest_Contrast2), file = "WTSP_Blood_FBnests_TanvsWhite.csv")


# Tan FB vs Tan BP nests
resMorphNest_Contrast3 <- results(ddsmorphnest, contrast=c("morphnest","Tan_FB","Tan_BP"))
head(resMorphNest_Contrast3)
sum(resMorphNest_Contrast3$padj < 0.1, na.rm=TRUE) ## 154
sum(resMorphNest_Contrast3$padj < 0.05, na.rm=TRUE) ## 74
write.csv(as.data.frame(resMorphNest_Contrast3), file = "WTSP_Blood_Tan_FBvsBP.csv")


# Tan vs White in BP nests
resMorphNest_Contrast4 <- results(ddsmorphnest, contrast=c("morphnest","Tan_BP","White_BP"))
head(resMorphNest_Contrast4)
sum(resMorphNest_Contrast4$padj < 0.1, na.rm=TRUE) ## 2
sum(resMorphNest_Contrast4$padj < 0.05, na.rm=TRUE) ## 2
write.csv(as.data.frame(resMorphNest_Contrast4), file = "WTSP_Blood_BPnests_TanvsWhite_test.csv")


vst <- varianceStabilizingTransformation (ddsmorphnest, blind=TRUE)
pcaExplorer(dds = ddsmorphnest, rlt = vst)



# MORPH controlling for sex
ddsMorph <- DESeqDataSetFromMatrix (countData = countData,
                                    colData = wtspDesign,
                                    design = ~sex + morph)
ddsMorph <- DESeq(ddsMorph)
ddsMorphClean <- ddsMorph[which(mcols(ddsMorph)$betaConv),]
resMorph <- results (ddsMorph) #change to ddsMorphClean to remove row that did not converge This is 'LOC102067071' which is one of two genes DE in BP T vs W above
head(resMorph)
sum(resMorph$padj < 0.1, na.rm=TRUE) ## 92
sum(resMorph$padj < 0.05, na.rm=TRUE) ## 58
write.csv(as.data.frame(resMorph), file = "WTSP_Blood_WvsT_controllingSex_outliersremoved_test.csv")
rld <- rlogTransformation (ddsMorphClean, blind=TRUE)
pcaExplorer(dds = ddsMorphClean, rlt = rld)



#PARENTAL CARE differences with nested nest
ddsNestedPC<- DESeqDataSetFromMatrix (countData = countData,
                                      colData = wtspDesign,
                                      design = ~ nested.nest + sex + morph + parentalcare)

ddsNestedPC <- DESeq(ddsNestedPC)
ddsNestedPC <- ddsNestedPC[which(mcols(ddsNestedPC)$betaConv),]
resnested <- results (ddsNestedPC)
head(resnested)
sum(resnested$padj < 0.1, na.rm=TRUE) ## 881
sum(resnested$padj < 0.05, na.rm=TRUE) ## 554
write.csv(as.data.frame(resnested), file = "WTSP_ParentalCare_controlling_nest_sex_morph.csv")
vst <- varianceStabilizingTransformation (ddsNestedPC, blind=TRUE)
pcaExplorer(dds = ddsNestedPC, rlt = vst)
plotPCA(vst, intgroup=c("whichnest"))


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plotCounts(ddsNestedPC,"NR3C1",intgroup = "parentalcare", normalized = TRUE,transform=TRUE,xlab="group")

#example:
d <- plotCounts(ddsNestedPC,"CHD4",intgroup = c("parentalcare"), returnData = TRUE, xlab="group")
quartz()
ggplot(d, aes(x = pairtype, y = count, colour = pairtype)) + ggtitle("CHD4", subtitle = "") + xlab("Group") + ylab("Normalized Counts") + scale_y_continuous() + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 5, width = 0.15) +  scale_colour_manual(values=cbPalette) + theme_classic() + theme(legend.position = "none",axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), plot.title = element_text(hjust = 0.5))

