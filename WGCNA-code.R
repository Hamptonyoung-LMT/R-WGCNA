install.packages("WGCNA")
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db","impute","preprocessCore"))
dir.create("helix-WGCNA-lesson3")
setwd("helix-WGCNA-lesson3")
fileUrl <-
"https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiv
er-Data.zip"
download.file(fileUrl,destfile="FemaleLiver-Data.zip")
unzip("FemaleLiver-Data.zip")
list.files()
## [1] "ClinicalTraits.csv" "FemaleLiver-Data.zip" "GeneAnnotation.csv"
## [4] "LiverFemale3600.csv"

###########################芯片数据载入#############################
library(WGCNA)
library(RColorBrewer)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
femData <- read.csv("LiverFemale3600.csv",header=T,sep=",",check.names=F)
head(femData[,1:13])
## substanceBXH gene_symbol LocusLinkID ProteomeID cytogeneticLoc CHROMOSOME
## 1 MMT00000044 1700007N18Rik 69339 286025 0 16
## 2 MMT00000046 Mast2 17776 157466 0 4
## 3 MMT00000051 Ankrd32 105377 321939 0 13
## 4 MMT00000076 0 383154 0 0 16
## 5 MMT00000080 Ldb2 16826 157383 0 5
## 6 MMT00000102 Rdhs 216453 0 10_70.0_cM 10
## StartPosition EndPosition F2_2 F2_3 F2_14 F2_15 F2_19
## 1 50911260 50912491 -0.01810 0.0642 6.44e-05 -0.05800 0.04830
## 2 115215318 115372404 -0.07730 -0.0297 1.12e-01 -0.05890 0.04430
## 3 74940309 74982847 -0.02260 0.0617 -1.29e-01 0.08710 -0.11500
## 4 49345114 49477048 -0.00924 -0.1450 2.87e-02 -0.04390 0.00425
## 5 43546124 43613704 -0.04870 0.0582 -4.83e-02 -0.03710 0.02510
## 6 1337265 1347607 0.17600 -0.1890 -6.50e-02 -0.00846 -0.00574
########################过滤过多缺失值的样品和基因####################
femData <- femData[femData$gene_symbol!="0",]
datExpr0 <- femData[,-c(1:8)]
rownames(datExpr0) <- femData$substanceBXH
datExpr0 <- t(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if(!gsg$allOK){
if(sum(!gsg$goodGenes)>0){
printFlush(paste("Removing genes:",paste(names(datExpr0)[!gsg$goodGenes],collapse=",")))
}
if (sum(!gsg$goodSamples)>0){
printFlush(paste("Removing samples:",paste(rownames(datExpr0)[!gsg$goodSamples],collapse=",")))
}
datExpr0 <- datExpr0[gsg$goodSamples,gsg$goodGenes]
}

#########################对样本进行聚类分析############################
###接下来，我们对样本进行聚类(与稍后对基因进行聚类形成对比)，以查看是否有明显的异常值。

sampleTree <- hclust(dist(datExpr0),method="average")
pdf("01.samples_cluster_tree.pdf",width=25,height=8)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="")
cutHeight <- 15
abline(h=cutHeight,col="red")
dev.off()
####更改了cutoff值
clust <- cutreeStatic(sampleTree,cutHeight=cutHeight,minSize=10)
clust
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples,]

############################读取样本表型数据###############################
traitData <- read.csv("ClinicalTraits.csv",header=T,sep=",",check.names=F)
head(traitData)
####################自己修改的代码，为保证表型数据与count数据一致#########
########################################################################
#############Loading clinical trait data###############################
traitData = BLCA_clinic
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(12, 11,2)];
allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.

###################################################
#将datexpr与临床恕不不一样的提出
#保留408个样本
femaleSamples = rownames(datExpr);
allTraits_newID = allTraits$newID
write.csv(femaleSamples,file="femaleSamples2.csv")
write.csv(allTraits_newID,file="allTraits_newID.csv")
delect_ID<-c("TCGA-BT-A2LA-11","TCGA-GC-A3WC-11","TCGA-CU-A0YR-11","TCGA-BT-A2LB-11","TCGA-CU-A0YN-11","TCGA-BL-A13J-11","TCGA-BT-A20U-11","TCGA-K4-A5RI-11","TCGA-BT-A20N-11","TCGA-GC-A6I3-11",         "TCGA-BT-A20W-11","TCGA-BT-A20Q-11","TCGA-K4-A54R-11","TCGA-GD-A3OQ-11","TCGA-GD-A2C5-11",
"TCGA-GD-A3OP-11","TCGA-K4-A3WV-11","TCGA-GC-A3BM-11","TCGA-BT-A20R-11")
delect_sort <-  sort(match(delect_ID,rownames(datExpr)),decreasing = T)
datExpr <- datExpr[-delect_sort,]
####################################################
femaleSamples = rownames(datExpr);

traitRows = match(femaleSamples, allTraits$newID);

datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 11];
collectGarbage();

#######################绘制样本聚类与表型展示图#############################
allTraits <- traitData[,c(11:15,17:30)]
rownames(allTraits) <- traitData[,2]
femaleSamples <- rownames(datExpr)
datTraits <- allTraits[femaleSamples,]
sampleTree2 <- hclust(dist(datExpr),method="average")
traitColors <- numbers2colors(datTraits,colors=blueWhiteRed(50))
pdf("01.cluster_trait.pdf",width=25,height=10)
plotDendroAndColors(sampleTree2,traitColors,
groupLabels=names(datTraits),
main="Sample dendrogram and trait heatmap")
dev.off()
save(datExpr, datTraits_1, file = "FemaleLiver-01-dataInput.RData")

####################picksoftthreshold函数筛选软阈值###############################
powers <- c(c(1:10),seq(from=12,to=30,by=2))
networkType <- "signed"
sft <- pickSoftThreshold(datExpr,
powerVector=powers,
networkType="signed",
verbose=5)
# networkType = "unsigned", adjacency = |cor|^power;
# networkType = "signed", adjacency = ((1+cor)/2)^power;
# networkType = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise;
# networkType = "distance", adjacency = (1-(dist/max(dist))^2)^power.
#####################绘图展示软阈值######################
pdf(file="02.soft_threshold.pdf",width=8,height=4.5)
par(mfrow=c(1,2))
cex1 <- 0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft threshold (power)",
ylab="Scale free topology model fit, signed R^2",
type="n",
main=paste("Scale independence"))
text(sft$fitIndices[,1],
-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],
xlab="Soft threshold (power)",ylab="Mean connectivity",
type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")
dev.off()

#软阈值的检验
#软阈值的检验
beta1=9 #根据自己的检验结果来选择
Connectivity=softConnectivity(datExpr,power=beta1,corOptions = "use = 'p'",
                              type = "unsigned",)
pdf("02.scalefree.pdf",16,8)
par(mfrow=c(1,2))
hist(Connectivity)
scaleFreePlot(Connectivity,main='Check scale free topology\n')
dev.off()
#下面代码有待商榷
ADJ1=abs(cor(datExpr,use='p'))^9
k=as.vector(apply(ADJ1,2,sum,na.rm=T))
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main='Check scale free topology\n')
dev.off()
enableWGCNAThreads()
####################################################################
##################################################################
########### One-step network construction and module detection#############
softPower <- 5
networkType <- "signed"
net <- blockwiseModules(datExpr,power=softPower,
                        networkType=networkType,numericLabels=TRUE,
                        mergeCutHeight=0.25,minModuleSize=30,
                        maxBlockSize=30000,saveTOMs=TRUE,
                        saveTOMFileBase="WGCNA_TOM",verbose=5)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
pdf(file="03.auto_modules_color.pdf",width=12,height=4)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels=FALSE,addGuide=TRUE)
dev.off()
save(datExpr,sft,softPower,networkType,net,moduleColors,
     file="WGCNA_2.RData")

##########################################################
###计算模块的ME
##########################################################

MEs0 <- moduleEigengenes(datExpr[,net$goodGenes],
                         moduleColors[net$blockGenes[[1]]])$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames(datExpr[,net$goodGenes])
text <- cbind(rownames(MEs),MEs)
colnames(text)[1] <- "samples"
write.table(text,file="05.module_eigengenes.xls",
            quote=F,sep="\t",row.names=F)


##########################################################
###############基于ME绘制模块聚类树#######################
names(MEs) <- substring(names(MEs),3)
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss),method="average")
pdf(file="05.modules_cluster_tree.pdf",width=7,height=5)
plot(METree,main="Clustering of module eigengenes",xlab="",sub="")
dev.off()

############################################################
############基于ME绘制热图##########################
moduleCor <- corAndPvalue(MEs,use="p")
rowLabels <- paste("ME",names(MEs),sep="")
textMatrix <- paste(signif(moduleCor$cor,2),
                    "\n(",signif(moduleCor$p,1),")",sep="")
dim(textMatrix) <- dim(moduleCor$cor)
pdf(file="05.modules_relationships.pdf",12,8)
par(mar=c(10,10,1,2))
labeledHeatmap(Matrix=moduleCor$cor,
               textMatrix=textMatrix,
               xLabels=rowLabels,
               yLabels=rowLabels,
               xSymbols=names(MEs),
               ySymbols=names(MEs),
               colorLabels=TRUE,
               colors=blueWhiteRed(50),
               setStdMargins=FALSE,
               xLabelsAngle=90,zlim=c(-1,1))
dev.off()

####保存模块间相关信息####
text <- paste("cor=",round(moduleCor$cor,4),
              ";p-value=",round(moduleCor$p,4),sep="")
dim(text) <- dim(moduleCor$cor)
rownames(text) <- rowLabels
colnames(text) <- rowLabels
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="05.modules_relationships.xls",
            quote=F,sep="\t",row.names=F)

###########################################################
#########绘制模块ME和节点基因表达热图##################
dir.create("05.expression_ME")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  dir <- "05.expression_ME/"
  pdf(file=paste(dir,"05.expression_ME_",which.module,".pdf",sep=""),
      30,10)
  ME <- MEs[,which.module]
  ME <- t(as.matrix(MEs[,which.module]))
  colnames(ME) <- rownames(datExpr[,net$goodGenes])
  layout(matrix(c(1,2)),heights=c(1.5,3))
  par(mar=c(0.3,9,3,5))
  plotMat(t(scale(datExpr[,net$goodGenes][,moduleColors[net$blockGenes[[1]]]==which.module])),
          nrgcols=30,rlabels=F,rcols=which.module,
          main=paste(which.module),cex.main=1)
  par(mar=c(5,4,0,1))
  barplot(ME,col=which.module,main="",cex.names=0.5,cex.axis=1,
          ylab="module eigengene",las=3)
  dev.off()
}

########################################################################
###########将模块与样本关联############################
##cor计算相关性
sample_cor <- cor(t(datExpr[,net$goodGenes]),t(datExpr[,net$goodGenes]),
                  use='pairwise.complete.obs')
moduleSampleCor <- cor(MEs,sample_cor,use="p")
nSamples <- nrow(datExpr[,net$goodGenes])
moduleSamplePvalue <- corPvalueStudent(moduleSampleCor,nSamples)
textMatrix <- paste(signif(moduleSampleCor,2),
                    "\n(",signif(moduleSamplePvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleSampleCor)
rowLabels <- paste("ME",names(MEs),sep="")
pdf(file="06.modules_samples_relationships.pdf",30,6)
par(mar=c(5,12,1,1))
labeledHeatmap(Matrix=moduleSampleCor,
               xLabels=colnames(sample_cor),
               yLabels=rowLabels,
               ySymbols=names(MEs),
               colorLabels=TRUE,
               colors=blueWhiteRed(50),
               setStdMargins=FALSE,
               xLabelsAngle=90,zlim=c(-1,1),
               cex.lab.x = 0.5
               )
dev.off()
#############保存模块与样本关联度信息#####################
text <- paste("cor=",round(moduleSampleCor,4),
              ";p-value=",round(moduleSamplePvalue,4),sep="")
dim(text) <- dim(moduleSampleCor)
rownames(text) <- rownames(moduleSampleCor)
colnames(text) <- colnames(moduleSampleCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="06.modules_samples_relationships.xls",
            quote=F,sep="\t",row.names=F)


##############将模块与表型关联##############################3
moduleTraitCor <- cor(MEs,datTraits_1,use="p")
nSamples <- nrow(datExpr[,net$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste(signif(moduleTraitCor,2),
                    "\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- paste("ME",names(MEs),sep="")
pdf(file="06.modules_traits_relationships.pdf",15,8)
par(mar=c(7,12,1,2))
labeledHeatmap(Matrix=moduleTraitCor,
               textMatrix=textMatrix,
               xLabels=colnames(datTraits_1),
               yLabels=rowLabels,
               ySymbols=names(MEs),
               colorLabels=TRUE,
               colors=blueWhiteRed(50),
               setStdMargins=FALSE,
               xLabelsAngle=90,zlim=c(-1,1))
dev.off()
#############保存模块与表型关联度信息#################
text <- paste("cor=",round(moduleTraitCor,4),
              ";p-value=",round(moduleTraitPvalue,4),sep="")
dim(text) <- dim(moduleTraitCor)
rownames(text) <- rownames(moduleTraitCor)
colnames(text) <- colnames(moduleTraitCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="06.modules_traits_relationships.xls",
            quote=F,sep="\t",row.names=F)


####################################################################
#################将节点基因与样本关联###################
modNames <- names(MEs)
geneModuleMembership <- cor(datExpr[,net$goodGenes],MEs,use="p")
nSamples <- nrow(datExpr[,net$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership,nSamples)
colnames(geneModuleMembership) <- paste("MM",modNames,sep="")
colnames(MMPvalue) <- paste("p.MM",modNames,sep="")
text <- paste("cor=",round(geneModuleMembership,4),
              ";p-value=",round(MMPvalue,4),sep="")
dim(text) <- dim(geneModuleMembership)
rownames(text) <- rownames(geneModuleMembership)
colnames(text) <- colnames(geneModuleMembership)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="07.genes_module_membership.xls",
            quote=F,sep="\t",row.names=F)


###########################################################
############将节点基因与表型(deadORlive)关联###################
modNames <- names(MEs)
deadORlive <- as.data.frame(datTraits_1$deadORlive)
names(deadORlive) <- "deadORlive"
geneTraitSignificance <- as.data.frame(cor(datExpr[,net$goodGenes],deadORlive,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(deadORlive),sep="")
names(GSPvalue) <- paste("p.GS.",names(deadORlive),sep="")
head(geneTraitSignificance)
head(GSPvalue)
#####绘制各模块节点基因与表型(deadORlive)关联图###############

dir.create("07.MM_vs_deadORlive")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColors[net$blockGenes[[1]]]==which.module
  dir <- "07.MM_vs_deadORlive/"
  pdf(file=paste(dir,"07.",which.module,"_MM_vs_deadORlive.pdf",sep=""),6,6)
  verboseScatterplot(geneModuleMembership[moduleGenes,column],
                     geneTraitSignificance[moduleGenes,1],
                     xlab=paste("Module membership (MM) in",which.module,"module"),
                     ylab="Gene significance for deadORlive",
                     main=paste("Module membership vs. gene significance\n"),
                     col=which.module)
  dev.off()
}
##保存模块节点基因与表型关联度结果#####
text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
write.table(text,file="07.genes_trait_significance.xls",
            quote=F,sep="\t",row.names=F)

######################导出共表达网络##############################
#载入TOM矩阵
load(file="WGCNA_TOM-block.1.RData")
ATOM <- as.matrix(TOM)
TOM1 <- ATOM[1:round((nrow(ATOM)/2)),1:round((nrow(ATOM)/2))]
TOM2 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),1:round((nrow(ATOM)/2))]
TOM3 <- ATOM[1:round((nrow(ATOM)/2)),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
TOM4 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
###提取各模块节点基因的TOM矩阵
dir.create("08.module_result")
setwd("08.module_result")
for(i in 1:(ncol(MEs)-1)) {
  module <- labels2colors(i)
  inModule <- moduleColors[net$blockGenes[[1]]]==module
  genename <- colnames(datExpr[,net$goodGenes])
  modGenes <- genename[inModule]
  modTOM1 <- TOM1[inModule[1:round((nrow(ATOM)/2))],
                  inModule[1:round((nrow(ATOM)/2))]]
  modTOM2 <- TOM2[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)],
                  inModule[1:round((nrow(ATOM)/2))]]
  modTOM3 <- TOM3[inModule[1:round((nrow(ATOM)/2))],
                  inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
  modTOM4 <- TOM4[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)],
                  inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
  modTOM <- rbind(cbind(modTOM1,modTOM3),cbind(modTOM2,modTOM4))
 


####保存节点基因对cytoscape输入为文件
cyt1 <- exportNetworkToCytoscape(modTOM,
                                 edgeFile=paste("CytoscapeInput-edges-",paste(module,collapse="-"),".txt",sep=""),
                                 nodeFile=paste("CytoscapeInput-nodes-",paste(module,collapse="-"),".txt",sep=""),
                                 threshold=0.02,
                                 nodeNames=modGenes,altNodeNames=modGenes,
                                 nodeAttr=moduleColors[net$blockGenes[[1]]][inModule])
####保存各模块节点基因和hub节点信息
IMConn <- softConnectivity(datExpr[,net$goodGenes][,modGenes])
out <- cbind(modGenes,IMConn)
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(module,"-module-gene.txt",sep=""),
            sep="\t",quote=F,row.names=F)
nTop <- 0.05*length(modGenes)
top <- (rank(-IMConn) <= nTop)
out <- cbind(modGenes[top],IMConn[top])
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(module,"-5%hubgene.txt",sep=""),
            sep="\t",quote=F,row.names=F) }


