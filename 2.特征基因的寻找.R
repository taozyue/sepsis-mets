#BiocManager::install('limma',ask = F,update = F)
library(limma)
# 没有先安装
#BiocManager::install('sva',ask = F,update = F)
library(sva)

# 运行以下代码去批次
mergeFile="merge.preNorm.txt"            #合并后的文件名称
normalizeFile="merge.normalzie.txt"      #矫正后的文件名称
# 如果不想用第二个，也可以就两个去批次,实际上只需训练集PR和MS合并
   files=c('GSE154918.txt',"GSE28750.txt",'GSE98895.txt')       #输入文件名称
#获取交集基因
geneList=list()
for(i in 1:length(files)){
  fileName=files[i]
  rt=read.table(fileName, header=T, sep="\t", check.names=F)
  header=unlist(strsplit(fileName, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
# 共同基因
intersectGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  fileName=files[i]
  header=unlist(strsplit(fileName, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(fileName, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对数值大的数据取log2
  #qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  #LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  #if(LogC){
  #rt[rt<0]=0
  #rt=log2(rt+1)}
  #rt=normalizeBetweenArrays(rt)
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(header[1],ncol(rt)))
}

allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file=mergeFile, sep="\t", quote=F, col.names=F)

#对数据进行批次矫正，输出矫正后的结果
normalizeTab=ComBat(allTab, batchType, par.prior=TRUE)
normalizeTab=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTab, file=normalizeFile, sep="\t", quote=F, col.names=F)

#合并前的PCA###################
#install.packages('ggplot2')
library(ggplot2)        #引用包
#读取输入文件,提取数据
rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project,group)
#group=group$group
#group=read.csv('pca.csv')

#write.csv(PCA,file = 'pca.csv',row.names = T)

PCA2=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], PC3=pcaPredict[,3],Type=Project)
#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图
#先绘制基础散点图：
p <- ggplot(data = PCA, aes(x = PC1,y = PC2))+
  geom_point(size = 2,
             aes(col = group, shape = Type)) #颜色区分实验/对照组，散点形状区分队列地区；
p

p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type,shape = group)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
ggsave(p,file = 'PCA_pre.pdf',width = 5.1,height = 3.8)

#合并后PCA###################
#读取输入文件,提取数据
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project,group)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图

p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type,shape = group)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(p)
ggsave(p,file = 'PCA_norm.pdf',width = 5.1,height = 3.8)

####################################WGCNA（双病联合版）################
#######################################################
#BiocManager::install('WGCNA',update=F,ask = F)
library(WGCNA)
options(stringsAsFactors = FALSE)
exprSet = read.table("merge.normalzie.txt",header = T,check.names = F,row.names = 1)
## 针对se1和MS作WGCNA
A=grep('GSE154918',colnames(exprSet))
B=grep('GSE98895',colnames(exprSet))
exprSet=exprSet[,c(A,B)]

load('se1.Rdata')
load('MS.Rdata')

group_list1 <- as.character(group_list1)
group_list3 <- as.character(group_list3)

#样本信息
group_list=c(group_list1,group_list3)

anno=data.frame(row.names=colnames(exprSet),group=group_list)


#dim(femData)
###选择最高变的5000基因，若容量大可以不选
#WGCNA_matrix = t(exprSet[order(apply(exprSet,1,mad), decreasing = T)[1:5000],])
WGCNA_matrix = t(exprSet[order(apply(exprSet,1,mad), decreasing = T)[1:nrow(exprSet)],])
datExpr0 <- WGCNA_matrix  ## top mad genes
datExpr0 <- as.data.frame(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # 返回TRUE则继续
# （可选）如果存在太多的缺失值
if (!gsg$allOK)
{
  # 把含有缺失值的基因或样本打印出来
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 去掉那些缺失值
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## 样本过滤前
  sampleTree = hclust(dist(datExpr0), method = "average")
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)

## 根据图片挑选cutheight,80可变
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust) # 0代表切除的，1代表保留的
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

## 更新anno
anno=anno[rownames(datExpr),,drop=F]

## 不调样本则用下面的代码
# datExpr=datExpr0

## 样本过滤后
sampleTree = hclust(dist(datExpr), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers (after)", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

## 制作更适合WGCNA的表型矩阵
anno$sepsis=ifelse(anno$group=="Seps_P",1,0)
anno$MS=ifelse(anno$group=='MS',1,0)
#anno$CT_blood=ifelse(anno$group=='CT',1,0)
#anno$MS_blood=ifelse(anno$group=='MS',1,0)

datTraits =anno[,-1]
datExpr=datExpr[rownames(datTraits),]
sampleNames = rownames(datExpr)
# 能全部对上
traitRows = match(sampleNames, rownames(datTraits))  

###power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,blockSize = 100000)

dev.off()
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
# 发现不合适就自定义
softPower=6
adjacency = adjacency(datExpr, power = softPower)

net = blockwiseModules(datExpr, power = softPower,
                       TOMType = "unsigned", minModuleSize = 60,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
# 显示模块数量以及各自包含的基因数目
# 0表示未分入任何模块的基因
# 1是最大的模块，往后依次降序排列，分别对应各自模块的基因
table(net$colors)


mergedColors = labels2colors(net$colors)

##手动保存
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]


nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 用color labels重新计算MEs（Module Eigengenes:模块的第一主成分）
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") #（这是重点）计算ME和表型相关性
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# 设置热图上的文字（两行数字：第一行是模块与各种表型的相关系数；
# 第二行是p值）
# signif 取有效数字
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 然后对moduleTraitCor画热图
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# 把各个module的名字提取出来（从第三个字符开始），用于一会重命名
modNames = substring(names(MEs), 3)
# 得到矩阵
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# 矩阵t检验
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# 修改列名
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


identical(rownames(datTraits),rownames(datExpr))

## 蓝色模块基因和脓毒症
# 先将感兴趣的表型提取出来，用于计算矩阵，注意更换名字
se = as.data.frame(datTraits$sepsis)
names(se) = "sepsis"
# 得到矩阵
geneTraitSignificance = as.data.frame(cor(datExpr, se, use = "p"))
# 矩阵t检验
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# 修改列名
names(geneTraitSignificance) = paste("GS.", names(se), sep="")
names(GSPvalue) = paste("p.GS.", names(se), sep="")  

### 相关性##########
##################你要的###################
module = "brown"
column = match(module, modNames) #找到目标模块所在列
moduleGenes = moduleColors==module #找到模块基因所在行
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for sepsis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


####模块基因和代谢综合征（MS）
# 先将感兴趣的表型提取出来，用于计算矩阵
MS = as.data.frame(datTraits$MS)
names(MS) = "MS"
# 得到矩阵
geneTraitSignificance = as.data.frame(cor(datExpr, MS, use = "p"))
# 矩阵t检验
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# 修改列名
names(geneTraitSignificance) = paste("GS.", names(MS), sep="")
names(GSPvalue) = paste("p.GS.", names(MS), sep="")  

### 相关性##########
##################你要的###################
module = "brown"
column = match(module, modNames) #找到目标模块所在列
moduleGenes = moduleColors==module #找到模块基因所在行
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for MS",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



# 选择导出模块
module = "brown"
# 选择模块中基因/探针
probes = names(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]
#modprobes可以后续分析
modProbes
write.table(modProbes,file ='brown_gene.txt',row.names = F,col.names = F,quote=F)
