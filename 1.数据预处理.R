
################### 第一个数据集fpkm#############################
#############################################################
# BiocManager::install('GEOquery',update = F,ask = F)
library(GEOquery)
gset=getGEO('GSE154918',getGPL = F,destdir = '.')
pdata=pData(gset[[1]])
pdata1=pdata[,c(1,43)]
group_list=pdata1$`status:ch1`

#挑选样本
Hlty=grep(('Hlty'),pdata1$`status:ch1`)
Seps_P=grep(('Seps_P'),pdata1$`status:ch1`)
pdata2=pdata1[c(Hlty,Seps_P),]
table(pdata2$`status:ch1`)

## 读取数据
se1=read.table('GSE154918_Schughart_Sepsis_200320.txt',header = T)
rownames(se1)=se1$gene_symbol
colnames(se1)
sematrix=se1[,-c(1:9)]
max(sematrix)

# 分组信息
group_list=pdata2$`status:ch1`
group_list
group_list1=factor(group_list,levels=c("Hlty","Seps_P"))
table(group_list1)
## 调整exprSet列名使其与pdata2的title列一致
ss=intersect(pdata2$title,colnames(sematrix))
exprSet1=sematrix[,ss]

exprSet1=as.data.frame(exprSet1)

#BiocManager::install('limma',update=F,ask=F)
library(limma)
design=model.matrix(~ group_list1)
fit=lmFit(exprSet1,design)
fit=eBayes(fit) 
allDiff_se1=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=0.05) 

#显著的上下调
allDiff_se1_up=allDiff_se1[allDiff_se1$logFC>0.5,]
allDiff_se1_down=allDiff_se1[allDiff_se1$logFC< -0.5,]
save(exprSet1,group_list1,allDiff_se1_up,allDiff_se1_down,file='se1.Rdata')
write.table(exprSet1,file='GSE154918.txt',quote = F,sep = '\t',col.names = NA)
load('se1.Rdata')
################################################
##################第二个数据集count##################
################################################
library(GEOquery)
library(dplyr)
library(tidyverse)
pdata=pData(gset[[1]])

##读取数据
exprSet2=read.table('GSE28750.txt',row.names = NULL,header = T,sep='\t',check.names = F)
rownames(exprSet2)=exprSet2$ID
exprSet2=exprSet2[,-1]

max(exprSet2)
## 赋值给exprSet2
exprSet2=log2(fpkm+1)
max(exprSet2)


# 分组信息
group_list=read.table('GSE28750_sample.txt',header = F,row.names = 1)
group_list2=group_list$V2
group_list2=factor(group_list2,levels=c("Normal",'sepsis'))
table(group_list2)


## 调整exprSet列名使其与group_list的title列一致
ss=intersect(rownames(group_list),colnames(exprSet2))
exprSet2=exprSet2[,ss]

exprSet2=as.data.frame(exprSet2)

library(limma)
design=model.matrix(~ group_list2)
fit=lmFit(exprSet2,design)
fit=eBayes(fit) 
#
allDiff_se2=topTable(fit,adjust='fdr',coef=2,number=Inf) 

allDiff_se2_up=allDiff_se2[allDiff_se2$logFC>0.5,]
allDiff_se2_down=allDiff_se2[allDiff_se2$logFC< -0.5,]

save(exprSet2,group_list2,allDiff_se2_up,allDiff_se2_down,file='se2.Rdata')
write.table(exprSet2,file='GSE28750.txt',quote = F,sep = '\t',col.names = NA)
load('se2.Rdata')

###############################################
#####代谢综合征################################
###############################################
library(GEOquery)
gset=getGEO('GSE98895',getGPL = F,destdir = '.')
#我们不要他处理好的数据集
#exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

x <- read.ilmn(files="GSE98895_non-normalized.txt",
               expr="Signal",
               probeid='PROBE_ID',
               other.columns="Detection Pval")

boxplot(log2(x$E),range=0,ylab="log2 intensity")
## Reading file GSE16997_raw.txt ... ...
# 背景校正和标准化
y <- neqc(x,detection.p="Detection Pval")
boxplot(log2(y$E),range=0,ylab="log2 intensity")
dim(y)
## 过滤
expressed <- rowSums(y$other$`Detection Pval` < 0.05) >= 3 ;table(expressed)
y <- y[expressed,]
dim(y)

#获取表达矩阵
exprSet=as.data.frame(y$E)

#install.packages(stringr)
colnames(exprSet)=stringr::str_remove(colnames(exprSet),'\\.AVG_')
max(exprSet)
dim(exprSet)

pdata$description2=stringr::str_replace(pdata$description,pattern=' ',replacement = '_')

## 调整顺序(必须pdata放前面才能调整顺序)
ss=intersect(pdata$description2,colnames(exprSet))
exprSet=exprSet[,ss]
colnames(exprSet)=rownames(pdata)

group_list=c(rep('MS',20),rep('CT',20))
group_list3=factor(group_list,levels = c('CT','MS'))


## ID转换
#install.packages('data.table')
library(data.table)
anno = fread("GPL6947-13512.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,14)] 
#!!!!可变！
anno = anno[!anno$Symbol== "",] 
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 

## 赋值给exprSet3
exprSet3=exprSet

design=model.matrix(~ group_list3)
fit=lmFit(exprSet3,design)
fit=eBayes(fit) 
#tumor处需修改case
allDiff_MS=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=0.05) 

allDiff_MS_up=allDiff_MS[allDiff_MS$logFC>0, ]
allDiff_MS_down=allDiff_MS[allDiff_MS$logFC< 0, ]

save(exprSet3,group_list3,allDiff_MS_up,allDiff_MS_down,file ='MS.Rdata')
write.table(exprSet3,file='GSE98895.txt',quote = F,sep = '\t',col.names = NA)
load('MS.Rdata')

load("se1.Rdata")
load("se2.Rdata")
###交集差异基因############################
#############
## 上调
dev.off()
#install.packages('VennDiagram')
library(VennDiagram)
venn.plot_up <- venn.diagram(
  x = list (
    sepsis_up  = rownames(allDiff_se1_up),
    MS_up=rownames(allDiff_MS_up)
  ),  
  cat.col=c("#e8c559","#ea9c9d"),
  fill = c("#e8c559","#ea9c9d"),
  filename = NULL
)

#将venn.plot通过grid.draw画到pdf文件中,手动保存
grid.draw(venn.plot_up)



## 下调
dev.off()
library(VennDiagram)
venn.plot_down <- venn.diagram(
  x = list (
    sepsis_down  = rownames(allDiff_se1_down),
    MS_down=rownames(allDiff_MS_down)
  ),  
  cat.col=c("#e8c559","#ea9c9d"),
  fill = c("#e8c559","#ea9c9d"),
  filename = NULL
)

#将venn.plot通过grid.draw画到pdf文件中，手动保存
grid.draw(venn.plot_down)



up=intersect(rownames(allDiff_se1_up),rownames(allDiff_MS_up))
write.table(up,file ='up.txt',quote = F,row.names = F,col.names = F)
down=intersect(rownames(allDiff_se1_down),rownames(allDiff_MS_down))
write.table(down,file ='down.txt',quote = F,row.names = F,col.names = F)

################veen 绘制
#文件名
A="Sepsis_down"
B="MS_down"
#构建一个列表
geneList=list()
geneNames=rownames(allDiff_se1_down)            
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[A]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("1",uniqLength,sep=" "))
geneNames=rownames(allDiff_MS_down)
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[B]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))

mycol <- distinctColorPalette(100)

pdf(file="DOWN.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()
