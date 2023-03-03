#加载R包
setwd("C:/生信分析最新思路总结/2.脓毒症和代谢（双疾病案例分析）(待完成)/5 单细胞数据处理")
getwd()
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(clustree)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(CellChat)

library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)

dir_name=c('N1','N2','S1','S2')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 250)
}

#修改名称
names(datalist)=dir_name
#添加分组信息
datalist[[1]]@meta.data$group <- 'normal'
datalist[[2]]@meta.data$group <- 'normal'
datalist[[3]]@meta.data$group <- 'sepsis'
datalist[[4]]@meta.data$group <- 'sepsis'
gc()

#计算线粒体 核糖体比例

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  datalist[[i]] <- sce
  rm(sce)
  gc()
}

sce <- merge(datalist[[1]],y=datalist[2:length(datalist)]) #merge 

#细胞数统计
raw_meta=sce@meta.data
raw_count <- table(raw_meta$orig.ident)
raw_count 
sum(raw_count)#27013
# 小提琴图查看metrics分布情况
cols <- c("#606f8a","#e8c559","#ea9c9d","#005496")
pearplot_befor<-VlnPlot(sce,group.by ='orig.ident', 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        pt.size = 0, 
                        ncol = 3,cols = cols)
pearplot_befor

#featurescatter
Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA',feature2 = 'nCount_RNA',group.by = 'orig.ident',cols = cols)
Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt',feature2 = 'nCount_RNA',group.by = 'orig.ident',cols = cols)
Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt',feature2 = 'nFeature_RNA',group.by = 'orig.ident',cols = cols)
Feature_ber1=Feature_ber1+theme(legend.position = 'none')
Feature_ber2=Feature_ber2+theme(legend.position = 'none')
Feature_ber<-ggarrange(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths = c(1,1,1.2))
Feature_ber

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 100 & 
              nFeature_RNA < 4000 & 
              percent.mt < 25 &
              nCount_RNA > 100 )
})

#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_meta=sce@meta.data
clean_count <- table(clean_meta$orig.ident)
clean_count
sum(clean_count)#26455
pearplot_after <- VlnPlot(sce,group.by ='orig.ident', 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          pt.size = 0, 
                          ncol = 3,cols = cols)
pearplot_after
save(datalist,file = 'datalist.RData')
gc()
load("datalist.RData")
#归一化
datalist <- lapply(datalist,function(x) {
  NormalizeData(x)
})
#寻找高变异度基因
datalist <- lapply(datalist, function(x) {
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
datalist
#整合成一个对象,寻找锚点
sce.anchors <- FindIntegrationAnchors(object.list = datalist, dims = 1:20)
#根据锚点来整合
sce  <- IntegrateData(anchorset = sce.anchors, dims = 1:20)
DefaultAssay(sce) <- "integrated" #更改默认数组

#保存sce文件
save(sce,file = 'sce.RData')


#对整合后的数据进行尺度变换
#ScaleData
all.genes <- rownames(sce[["RNA"]]@data) #针对所有基因
length(all.genes)
sce <- ScaleData(sce, features = all.genes)
#####
#Normalizing the data
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
#Identification of highly variable features (feature selection) 高变基因
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
#PCA降维，选择合适的拐点
sce <- RunPCA(sce, features = VariableFeatures(sce))
sce <- RunTSNE(sce, features = VariableFeatures(sce))
sce <- RunUMAP(sce, features = VariableFeatures(sce))
dimplot1 <- DimPlot(sce, reduction = "pca",group.by = 'orig.ident') 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1

figs1=ggarrange(Feature_ber,pearplot_befor,pearplot_after,sc_pca,
                nrow = 4,labels = c('A','B','C','D'),ncol = 1,widths = c(1,1,1,2))
ggsave(filename = 'FigS1.pdf',plot = figs1,he=12,width = 10)
ggsave(filename = 'FigS1.jpg',plot = figs1,he=12,width = 10,dpi = 300)


Dims <- 20
#降维
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)
sce <- RunUMAP(sce, 
               dims=1:Dims, 
               reduction="pca")

figs1=ggarrange(Feature_ber,pearplot_befor,pearplot_after,sc_pca,
                nrow = 4,labels = c('A','B','C','D'),ncol = 1,widths = c(1,1,1,2))
ggsave(filename = 'FigS1.pdf',plot = figs1,he=12,width = 10)
ggsave(filename = 'FigS1.jpg',plot = figs1,he=12,width = 10,dpi = 300)


library(clustree)
sce <- FindNeighbors(sce, dims = 1:Dims)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1,.1))
)
colnames(sce@meta.data)
clustree(sce@meta.data, prefix = "integrated_snn_res.")

pdf('clust.snn_res.pdf',he=15,wi=15)
clustree(sce@meta.data, prefix = "integrated_snn_res.")
dev.off()

colnames(sce@meta.data)

#保存sce文件
save(sce,file = 'sce.RData')



#聚类分析
Resolution <- 0.7
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#clustree(sce)
colnames(sce@meta.data) 
pdf('tsne_multi_samples_combined.pdf',width = 11,height = 6)
p1 <- DimPlot( sce , reduction = "tsne", group.by = "group")
p2 <- DimPlot(sce , reduction = "tsne", label = TRUE)
p1 + p2
dev.off()

#保存sce文件
save(sce,file = 'sce.RData')
load("sce.RData")

  ###3.细胞类型鉴定（SingleR）
  ### 安装SingleR   BiocManager::install('SingleR')
  dir.create('cell_identify')
  # 
  library(SingleR)
  refdata <- HumanPrimaryCellAtlasData()
  
  testdata <- GetAssayData(sce, slot="data")
  clusters <- sce@meta.data$integrated_snn_res.0.7
  cellpred <- SingleR(test = testdata, ref = refdata,
                      labels =refdata$label.main,
                      method = "cluster", clusters = clusters, 
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
  write.csv(celltype,"cell_identify/celltype_singleR.csv",row.names = F)
  sce@meta.data$celltype = "NA"
  for(i in 1:nrow(celltype)){
    sce@meta.data[which(sce@meta.data$integrated_snn_res.0.7 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  
  p2 = DimPlot(sce, group.by="celltype", label=T, label.size=4.5, reduction='tsne')+ggsci::scale_color_lancet()
  p2
  
  ggsave("cell_identify/tsne_celltype.pdf", p2, width=7 ,height=6)
  
  markers <- c("CD2","CD69","CD3D",   # T_cells
               "GNLY","IL7R",  #NK_cell
               "CD79A",   #B_cell
               "CD14","CD74","LYZ"# Monocyte
               )
                
               genes <- list("T_cells" = c("CD2","CD69","CD3D"),
                             "Monocyte"=c("CD14","CD74","LYZ"),
                             "NK_cell" = c("GNLY","IL7R","CD160","NKG7"),
                             "B cells" = c("CD79A","CD37"))

                             
  pdf(file = "05.ann_cluster_marker.pdf",width =20,height = 7)
  do_DotPlot(sample = sce,features = genes,dot.scale = 12,colors.use = c("yellow","red"),legend.length = 50,
             legend.framewidth = 2, font.size =12)
  dev.off()
  

  
  
  
  
  
  
  

  load("sce.RData")
 DefaultAssay(sce)="RNA"
  load('rf_lasso_genes.Rdata')
VlnPlot(sce,features = genes,group.by = 'celltype',pt.size = 0)

VlnPlot(sce,features = genes,group.by = 'group',pt.size = 0)

save(sce,file = 'sce.RData')
load("sce.RData")
#############
df <- as.data.frame(table(sce$celltype, sce$orig.ident))
colnames(df) <- c("celltype","sample","value")
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill", width = 0.9) + 
  ggsci::scale_color_lancet() +
  scale_y_continuous(labels = scales::percent) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90), axis.title = element_blank())
ggsave('celltype_proportion_1.pdf', p, width = 6, height = 4)
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill", width = 0.9) + 
  ggsci::scale_color_lancet() +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90), axis.title = element_blank())
ggsave('celltype_proportion_2.pdf', p, width = 6, height = 3.5)

############################################
###### 细胞分类结果可视化
dir.create("Figures")

### 全局展示细胞类型
## umap图
colnames(sce@meta.data)
p1 = DimPlot(sce, reduction = "tsne", group.by = "celltype", label = T)
p2 = DimPlot(sce, reduction = "umap", group.by = "celltype", label = T)
pc = p1 + p2
pc
ggsave("Figures/celltype_main_manual.pdf", pc, width = 15, height = 6.5)
# 色板调整颜色
p1 = p1 + ggsci::scale_color_npg()
p2 = p2 + ggsci::scale_color_npg()
pc = p1 + p2
ggsave("Figures/celltype_main_manual_c1.pdf", pc, width = 15, height = 6.5)
# 手工调整颜色
#https://colorhunt.co/palettes
mycol = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
          "#7579E7","#F08080","#1E90FF","#F037A5","#52734D")
p1 = p1 + scale_color_manual(values = mycol)
p2 = p2 + scale_color_manual(values = mycol)
pc = p1 + p2
p1
ggsave("Figures/celltype_main_manual_c2.pdf", pc, width = 15, height = 6.5)








#############################################
DimPlot(sce, reduction = "tsne",group.by = "celltype",
        label = T,label.box = T,label.size = 3,repel = T) + theme_bw() +
  labs( x= "tSNE 1",y= "tSNE 2",title = "cell type") +
  theme(panel.grid=element_blank(), # 去网格线
        plot.title = element_text(size = 15,color="black",hjust = 0.5),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black'),
        axis.ticks = element_line(color = 'black', lineend = 'round'),
        legend.position = 'bottom',
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black'),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
  scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))








########################
## 重要热图
table(sce@meta.data$celltype)
# 加入tissue_cell列
tissue_cell=paste0(sce@meta.data$group,"_",sce@meta.data$celltype)
sce <- AddMetaData(sce,tissue_cell,'tissue_cell')
table(sce@meta.data$tissue_cell)

##生成一个参考矩阵，并不重要，后面都要替换里面的数字
ap <- AverageExpression(sce,features = genes,group.by ='tissue_cell' ,slot = 'data')[[1]]
## 循环
rowgroup=genes
colgroup=colnames(ap)
## 循环的逻辑是：先取每种细胞，再取每个基因，计算该细胞种类中阳性(表达>0)的细胞比例占所有细胞个数的百分比
## 下面的循环的巧妙之处在于提取自己sobj，然后在sobj中巧妙地提取S4对象中的RNA矩阵
## 稍微等几十秒
for( i in rowgroup){
  for(j in colgroup){
    sobj <- subset(sce,tissue_cell==j)
    ap[i,j] <- as.numeric(table(sobj@assays$RNA@data[i,]>0)[2])/length(sobj@assays$RNA@data[i,])*100
  }}
## NA是0，替换
ap[is.na(ap)] <- 0
colnames(ap)
## 热图
### ap:average proportion
#ap=ap[,c('CT_DC','PR_DC','CT_Keratinocytes','PR_Keratinocytes','CT_Monocyte','PR_Monocyte','CT_NK_cell','CT_NK_cell','CT_T_cells','PR_T_cells')]
ap=ap[,c('normal_B_cell','sepsis_B_cell','normal_Monocyte','sepsis_Monocyte'
         ,'sepsis_Neutrophils','normal_NK_cell','sepsis_NK_cell'
         ,'normal_Platelets','sepsis_Platelets','normal_T_cells','sepsis_T_cells')]

p=pheatmap::pheatmap(ap,display_numbers = T,  
                   color = colorRampPalette(c(rep("white",1), rep("firebrick3",1)))(100),
                   cluster_rows = F,
                   cluster_cols = F,angle_col = 45,main = 'Proportion')
p
ggsave(filename = 'pheatmap_proportion.pdf',plot = p,width = 6,height = 5)


#############表达量热图##########################
#################################################
sce@meta.data$tissue_cell=paste0(sce@meta.data$group,'_',sce@meta.data$celltype)
library(Seurat)

ae=AverageExpression(sce,assays = 'RNA',group.by = 'tissue_cell',features = genes)
ae=as.data.frame(ae$RNA)
ae=log2(ae+1)
ae=na.omit(ae)
ae=ae[,c('normal_B_cell','sepsis_B_cell','normal_Monocyte','sepsis_Monocyte'
         ,'sepsis_Neutrophils','normal_NK_cell','sepsis_NK_cell'
         ,'normal_Platelets','sepsis_Platelets','normal_T_cells','sepsis_T_cells')]
p=pheatmap::pheatmap(ae,display_numbers = T,  
                   color = colorRampPalette(c(rep("white",1), rep("firebrick3",1)))(100),
                   cluster_rows = F,
                   cluster_cols = F,angle_col = 45 ,main='Expression')
p
ggsave(filename = 'pheatmap_expression.pdf',plot = p,width = 6,height = 5)
genes
sce$celltype
Idents(sce)="celltype"
p2=FeaturePlot(sce, features = c("CASP4", "CFLAR"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T,
            split.by = 'group',label = T) 
p2
ggsave(filename = 'CASP4-CFLAR.pdf',plot = p2,width = 10,height = 5)


p3=FeaturePlot(sce, features = c("CFLAR"), cols =c("lightgrey", "orange"),pt.size = .9, order=T,
               split.by = 'group',label = T) 
p3
ggsave(filename = 'CFLAR.pdf',plot = p3,width = 8,height = 4)

load('GS.Rdata')
sce <-AddModuleScore(sce, features= list,name = names(list))
names(x = sce[[]])

## 选择感兴趣的代谢通路
##提取细胞子集(Monocyte)
sce$celltype
Cells.sub <- subset(sce@meta.data, celltype=='Monocyte')
scRNAsub <- subset(sce, cells=row.names(Cells.sub))
Idents(scRNAsub)=scRNAsub$group
Monocyte_normal <- subset(scRNAsub, idents = c("normal"))
Monocyte_sepsis <- subset(scRNAsub, idents = c("sepsis"))
scRNAsub$M5937_HALLMARK_GLYCOLYSIS35

VlnPlot(scRNAsub, features="M5937_HALLMARK_GLYCOLYSIS35", group.by = "group", pt.size = 0,) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Glycolysi score", x="Monocyte") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
wilcox.test(Monocyte_normal$M5937_HALLMARK_GLYCOLYSIS35, Monocyte_sepsis$M5937_HALLMARK_GLYCOLYSIS35, alternative = "two.sided") 
#p-value < 2.2e-16

## 选择感兴趣的代谢通路
##提取细胞子集（表皮）
table(sce$celltype)
Cells.sub <- subset(sce@meta.data, celltype=='NK_cell')
scRNAsub <- subset(sce, cells=row.names(Cells.sub))
Idents(scRNAsub)=scRNAsub$group
NK_normal <- subset(scRNAsub, idents = c("normal"))
NK_sepsis <- subset(scRNAsub, idents = c("sepsis"))

VlnPlot(scRNAsub, features="M5937_HALLMARK_GLYCOLYSIS35", group.by = "group", pt.size = 0,) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Glycolysi score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
wilcox.test(NK_normal$M5937_HALLMARK_GLYCOLYSIS35, NK_sepsis$M5937_HALLMARK_GLYCOLYSIS35, alternative = "two.sided") #p-value < 2.2e-16

sce$Glycolysi_score=sce$M5937_HALLMARK_GLYCOLYSIS35
p=FeaturePlot(sce, features = c("CFLAR", "Glycolysi_score"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T,
               split.by = 'group',label = T) 

p
ggsave(filename = 'Glycolysi_score_CFLAR.pdf',plot = p,width = 14,height = 8)
