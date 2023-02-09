############CIBERSORT免疫浸润
#############################################
#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")          #引用包
#运行CIBERSORT，得到免疫细胞含量结果
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "merge.normalzie.txt", perm=100, QN=F)


#############脓毒症中可视化CIBERSORT的结果###########################
# 每包先安装
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("vioplot")


#引用包
library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"      #免疫细胞浸润文件

pFilter=0.05          #免疫细胞浸润结果的过滤条件

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
#根据P过滤，这些不适合分析
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

## 读取表达矩阵
gene=read.csv('genes.csv')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]

#按照分型对样品分组
load('se1.Rdata')
load('se2.Rdata')
load('MS.Rdata')

#### 先看脓毒症中
colnames(rt)
a=grep('GSE154918',colnames(rt))
b=grep('GSE28750',colnames(rt))
rt_pr=rt[,c(a,b)]
group1=ifelse(group_list1=='Hlty','Normal','sepsis')
group2=ifelse(group_list2=='Normal','Normal','sepsis' )
group=c(group1,group2)

anno=data.frame(row.names = colnames(rt_pr),group=group)
dev.off()
# 由于前面过滤了P，要取交集
ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_pr=rt_pr[,ss]

# 更新group
group=anno$group

lowName=row.names(anno)[anno[,1]=='Normal']
highName=row.names(anno)[anno[,1]=='sepsis']

#提取不同分型的免疫细胞含量
lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)

#绘制小提琴图
outTab=data.frame()
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n",      
     cex.main=2,cex.lab=1.4,cex.sub=1.2)

#对每个免疫细胞循环，绘制vioplot，低表达用蓝色表示，高表达用红色表示
bioCol=c("#0072b5","#bc3c29","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
  if(sd(rt[1:lowNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
    rt[(lowNum+1),i]=0.00001
  }
  lowData=rt[1:lowNum,i]
  highData=rt[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Normal", "sepsis"),
       lwd=4.5,bty="n",cex=1,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1.3,srt = 30,pos=2)


#输出免疫细胞和p值表格文件
write.table(outTab,file="diff_result_pr.txt",sep="\t",row.names=F,quote=F)


### 脓毒症免疫细胞和基因相关性#######
##############################
##############################
#没包先安装
#install.packages('Hmisc')

library(Hmisc)

#引用包
library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"      #免疫细胞浸润文件

pFilter=0.05          #免疫细胞浸润结果的过滤条件

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
#根据P过滤，这些不适合分析
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
immune

## 读取表达矩阵
gene=read.csv('genes.csv')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]

#按照分型对样品分组
load('se1.Rdata')
load('se2.Rdata')
load('MS.Rdata')

#### 先看脓毒症中相关性
colnames(rt)
a=grep('GSE154918',colnames(rt))
b=grep('GSE28750',colnames(rt))
rt_pr=rt[,c(a,b)]
group1=ifelse(group_list1=='Hlty','Normal','sepsis')
group2=ifelse(group_list2=='Normal','Normal','sepsis' )
group=c(group1,group2)

anno=data.frame(row.names = colnames(rt_pr),group=group)

# 取交集
ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_pr=rt_pr[,ss]

#相关性检验
nc =cbind(immune,t(rt_pr))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
m=as.data.frame(m)
#!!!!!!!!!基因要换
m =dplyr::filter(m,m$STOM != 'NaN')
p = rcorr(nc)$P[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
p =p[rownames(m),]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


library(pheatmap)

## 手动保存
pheatmap(t(m),
         display_numbers =t(tmp),
         angle_col =45,
         color = colorRampPalette(c("#0072b5", "white", "#bc3c29"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)

######################################################
##################MS的免疫结果###########################

# 每包先安装
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("vioplot")


#引用包
library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"      #免疫细胞浸润文件

pFilter=0.05          #免疫细胞浸润结果的过滤条件

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
#根据P过滤，这些不适合分析
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

## 读取表达矩阵
gene=read.csv('genes.csv')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]

#按照分型对样品分组
load('se1.Rdata')
load('se2.Rdata')
load('MS.Rdata')

#### 看代谢综合征中
a=grep('GSE98895',colnames(rt))
rt_MS=rt[,a]
anno=data.frame(row.names = colnames(rt_MS),group=group_list3)

# 由于前面过滤了P，要取交集
ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_MS=rt_MS[,ss]

# 更新group
group=anno$group

lowName=row.names(anno)[anno[,1]=='CT']
highName=row.names(anno)[anno[,1]=='MS']

#提取不同分型的免疫细胞含量
lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)

#绘制小提琴图
outTab=data.frame()
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n",     
     cex.main=2,cex.lab=1.4,cex.sub=1.2)

#对每个免疫细胞循环，绘制vioplot，低表达用蓝色表示，高表达用红色表示
bioCol=c("#20854e","#e18727","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
  if(sd(rt[1:lowNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
    rt[(lowNum+1),i]=0.00001
  }
  lowData=rt[1:lowNum,i]
  highData=rt[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("CT_blood", "MS_blood"),
       lwd=3.0,bty="n",cex=1,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1.3,srt = 30,pos=2)


#输出免疫细胞和p值表格文件
write.table(outTab,file="diff_result_MS.txt",sep="\t",row.names=F,quote=F)


### 免疫细胞和基因相关性#######
##############################
##############################
#没包先安装
#install.packages('Hmisc')

library(Hmisc)

#引用包
library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"      #免疫细胞浸润文件

pFilter=0.05          #免疫细胞浸润结果的过滤条件

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
#根据P过滤，这些不适合分析
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


## 读取表达矩阵
gene=read.csv('genes.csv')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]



#### 看代谢综合征中相关性
a=grep('GSE98895',colnames(rt))
rt_MS=rt[,a]
anno=data.frame(row.names = colnames(rt_MS),group=group_list3)

# 取交集
ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_MS=rt_MS[,ss]


#相关性检验
nc =cbind(immune,t(rt_MS))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
m = as.data.frame(m)
#!!!!!!!!!基因要换！！！！！！！！！
m =dplyr::filter(m,m$STOM != 'NaN')
p = rcorr(nc)$P[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
p =p[rownames(m),]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


library(pheatmap)

## 手动保存
pheatmap(t(m),
         display_numbers =t(tmp),
         angle_col =45,
         color = colorRampPalette(c("#20854e", "white", "#e18727"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)






##########################################
#####脓毒症的代谢情况########################
#########################################
BiocManager::install('GSVA',update = F,ask = F)
#BiocManager::install('Biobase',update = F,ask = F)
load('GS.Rdata')

library(GSVA)
library(Biobase)

data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

uni_matrix=data

list= list


gsva_matrix<- gsva(as.matrix(uni_matrix), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

##前5个字符
rownames(gsva_matrix)=sub('^......','',rownames(gsva_matrix))

save(gsva_matrix,file = 'gsva_array.Rdata')
load('gsva_array.Rdata')

metabolism=gsva_matrix[c('HALLMARK_HYPOXIA',
                         'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
                         'HALLMARK_ADIPOGENESIS',
                         'HALLMARK_XENOBIOTIC_METABOLISM',
                         'HALLMARK_FATTY_ACID_METABOLISM',
                         'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                         'HALLMARK_GLYCOLYSIS',
                         'HALLMARK_HEME_METABOLISM',
                         'HALLMARK_BILE_ACID_METABOLISM'),]
                         
## 读取表达矩阵
gene=read.csv('genes.csv')
gene=gene$x
rt=data[gene,]


#按照分型对样品分组
load('se1.Rdata')
load('se2.Rdata')
load('MS.Rdata')

#### 先看脓毒症中
a=grep('GSE154918',colnames(rt))
b=grep('GSE28750',colnames(rt))
rt_pr=rt[,c(a,b)]
group1=ifelse(group_list1=='Hlty','Normal','sepsis')
group2=ifelse(group_list2=='Normal','Normal','sepsis' )
group=c(group1,group2)

anno=data.frame(row.names = colnames(rt_pr),group=group)

# 要取交集
metabolism_pr=metabolism[,rownames(anno)]
metabolism_pr=t(metabolism_pr)

library(Hmisc)
#相关性检验
nc =cbind(metabolism_pr,t(rt_pr))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(metabolism_pr),(ncol(nc)-length(gene)+1):ncol(nc)]
rownames(m)=stringr::str_remove(rownames(m),pattern = 'HALLMARK_')
p = rcorr(nc)$P[1:ncol(metabolism_pr),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


library(pheatmap)


## 手动保存
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#0072b5", "white", "#bc3c29"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)


##########################################
#####MS代谢情况########################
#########################################

load('gsva_array.Rdata')

metabolism=gsva_matrix[c('HALLMARK_HYPOXIA',
                         'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
                         'HALLMARK_ADIPOGENESIS',
                         'HALLMARK_XENOBIOTIC_METABOLISM',
                         'HALLMARK_FATTY_ACID_METABOLISM',
                         'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                         'HALLMARK_GLYCOLYSIS',
                         'HALLMARK_HEME_METABOLISM',
                         'HALLMARK_BILE_ACID_METABOLISM'),]

## 读取表达矩阵
gene=read.csv('genes.csv')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)
rt=data[gene,]


#### 看代谢综合征代谢相关性
a=grep('GSE98895',colnames(rt))
rt_MS=rt[,a]
anno=data.frame(row.names = colnames(rt_MS),group=group_list3)


#要取交集
metabolism_MS=metabolism[,rownames(anno)]
metabolism_MS=t(metabolism_MS)

#相关性检验
nc =cbind(metabolism_MS,t(rt_MS))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(metabolism_MS),(ncol(nc)-length(gene)+1):ncol(nc)]
rownames(m)=stringr::str_remove(rownames(m),pattern = 'HALLMARK_')
p = rcorr(nc)$P[1:ncol(metabolism_MS),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


library(pheatmap)


## 手动保存
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#20854E", "white", "#E18727"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)

