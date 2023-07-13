library(Seurat)
library(ggplot2)
library(data.table)
library(hdf5r)
library(dplyr)
library("ggthemes")
library("RColorBrewer")
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)

setwd("C:\\addIH\\Fig1")

B<-readRDS("B.rds")
setwd("C:\\addIH\\Fig4")
B <- FindVariableFeatures(B, selection.method = "vst", nfeatures = 2000) 
al.genes<-row.names(B) 
B <- ScaleData(B, features = al.genes) 
B <- RunPCA(B, features = VariableFeatures(object = B), npcs = 50)   
B <- FindNeighbors(B, reduction = "pca", dims = 1:30) 
B <- FindClusters(B, resolution = 0.2)
B <- RunTSNE(B,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(B, reduction = "tsne",label="T",label.size=3.5,pt.size = 1)
B <- RunUMAP(B,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(B, reduction = "umap",label="T",label.size=3.5,pt.size = 1)

B.markers <- FindAllMarkers(B, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- B.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(B, features = top10$gene) + NoLegend()
dev.off()

saveRDS(B,"B.rds")

B<-readRDS("B.rds")
###############################
#平均表达谱函数AverageExpression
AverageExp<-AverageExpression(B,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pdf("corr.pdf",width=20,height=20)
pheatmap(coorda$r)
dev.off()
##################################
##堆叠小提琴图
#读入seurat处理后的rds文件，以下是R代码
library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#配色方案
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

#MHC I
StackedVlnPlot(B, c("Tapbp","H2-Q1","H2-Q2","H2-Q4","H2-Q6","H2-Q7","H2-Q10"), pt.size=0, cols=my36colors)
#MHC II
StackedVlnPlot(B, c('H2-Ab1', 'H2-Aa', 'H2-Eb1', 'H2-Eb2', 'H2-Oa', 'H2-Ob',
                    'H2-DMa', 'H2-DMb1', 'H2-DMb2'), pt.size=0, cols=my36colors)
#B细胞发育
StackedVlnPlot(B, c("Cd19","Cd34","Cd38"), pt.size=0, cols=my36colors)

#LZ score
#StackedVlnPlot(B, c("Aicda"), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c("Cd83","Nfkbia", "Bcl2a1a", "Egr1","Il4i1"), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c("Aicda","Ccnd2","Ccr6", "S1pr3", "Myc"), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c("Gpr183","Cd69","Cd40", "Cd86","Egr2"), pt.size=0, cols=my36colors)
#DZ score
StackedVlnPlot(B, c('Cdc20', 'Cxcr4', 'Ccnb2', 'Ccnb1'),pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Hmmr', 'Polh', 'Lmo4'), pt.size=0, cols=my36colors)
#生发中心
StackedVlnPlot(B, c('Ccl19', 'Cxcl13', 'Ccr7', 'Cxcr5', 'Sell', 'Lamp3'), pt.size=0, cols=my36colors)


#clusters
StackedVlnPlot(B, c('Tsc22d3', 'Ccr7', 'Klf2', 'Cd83'), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Cd52', 'Serpinb1a', 'S1pr4', 'Ms4a4c'), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Plac8', 'H2-K1', 'H2-Q7', 'Ly6a'), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Aicda', 'S1pr2', 'Lipc', 'Rgs1'), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Birc5', 'Ube2c', 'Cdc20', 'Stmn1'), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Ly6c1', 'Ly6c2', 'Xbp1', 'Jchain'), pt.size=0, cols=my36colors)
StackedVlnPlot(B, c('Hsp90ab1', 'Eif5a', 'C1qbp', 'Mif'), pt.size=0, cols=my36colors)

#IG-specific
StackedVlnPlot(B, c('Pax5','Rel','Mcl1','Mycbp2'))

#Plasma B
StackedVlnPlot(B, c('Ly6c1', 'Ly6c2', 'Xbp1', 'Jchain'), pt.size=0, cols=my36colors)
#APB
StackedVlnPlot(B, c('Ccr7', 'Cxcr4', 'Cxcr5'), pt.size=0, cols=my36colors)

#memory B 
rownames(B)[grep("Hla",rownames(B))]
#StackedVlnPlot(B, c('IgA1','Ig', 'IgA2', 'IgG1','IgG2', 'IgG3', 'IgM'))
StackedVlnPlot(B, c('Ccr6','Cd79a','Ms4a1','Cd19','Casp3'))
StackedVlnPlot(B, c('Cotl1','Cd82','Aim2','Bank1','Cd80'))
StackedVlnPlot(B, c('Ms4a1','Gpr183','Blk','Scimp'))

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7442689/
StackedVlnPlot(B, c('Cd38','Cd27','Zeb2','S1pr2','Ighd'))

StackedVlnPlot(B, c('Plp2','Plac8','Maml2','Plxnc1'))
StackedVlnPlot(B, c('P2ry10','Celf2','Rnase6'))

#human
#GC:0,2
pdf("GC.pdf")
StackedVlnPlot(BC, c("BCL6","AICDA"), pt.size=0, cols=my36colors)
dev.off()

pdf("GC.light.pdf")
StackedVlnPlot(BC, c("CD83","CD86"), pt.size=0, cols=my36colors)
dev.off()
#Bfh
pdf("GC.dark.pdf")
StackedVlnPlot(BC, c("CXCR4"), pt.size=0, cols=my36colors)
dev.off()
#Follicular:3
pdf("Follicular.pdf")
StackedVlnPlot(BC, c("CD22","IGHD","FCER2","CR2"), pt.size=0, cols=my36colors)
dev.off()
#plasma:0
pdf("Plasma.pdf")
StackedVlnPlot(BC, c("SDC1","PRDM1","CD38"), pt.size=0, cols=my36colors)
dev.off()
#memory:3
pdf("Memory.pdf")
StackedVlnPlot(BC, c("IGHA1","MS4A1","CD40","CD80"), pt.size=0, cols=my36colors)
dev.off()



0,1,7,11;t1
0,1;t2
0,1,7,11;t3

#B.bak<-B
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)

B.markers <- FindAllMarkers(B, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- B.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.cell.pdf",width=25,height=30)
DoHeatmap(B, features = top10$gene) + NoLegend()
dev.off()


#################### 挑选样品分析比例
Idents(B) <- "stim"
B1<-subset(B,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(B) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")

names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(B1) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B1)
B1 <- RenameIdents(B1, new.cluster.ids)
#B1@active.ident=factor(B1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

pdf("DimPlot.tsne.pdf",width=8,height=5)
DimPlot(B1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.tsne.jpg",width=800,height=500)
DimPlot(B1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.tsne.eps",width=8,height=5)
DimPlot(B1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

pdf("DimPlot.umap.pdf",width=10,height=5)
DimPlot(B1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.umap.jpg",width=1000,height=500)
DimPlot(B1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.umap.eps",width=10,height=5)
DimPlot(B1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

B1$batch<-factor(B1$batch,levels = c("WT","AS","AAA","IH","IG","AG"))
pdf("DimPlot.umap-batch.pdf",width=18,height=3)
DimPlot(B1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.umap-batch.jpg",width=1800,height=300)
DimPlot(B1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.umap-batch.eps",width=18,height=3)
DimPlot(B1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

pdf("DimPlot.tsne-batch.pdf",width=18,height=3)
DimPlot(B1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.tsne-batch.jpg",width=1800,height=300)
DimPlot(B1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.tsne-batch.eps",width=18,height=3)
DimPlot(B1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()


### percent
df1<-table(B1@active.ident,B1$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#df1$Type<-factor(df1$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("percent_barplot.pdf",width=5,height=6)
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()
jpeg("percent_barplot.jpg")
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

setEPS()
postscript("percent_barplot.eps",width=5,height=6)
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

pdf("meta-percent.pdf",width=13,height=4)
#jpeg("metafig1.jpg",width=800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
# scale_colour_manual(breaks = df1$Type, values = colors) +
#  scale_colour_manual(values = colors, 0.5) +
  theme_bw() +  theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("metafig-percent.jpg",width=1800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
#  scale_colour_manual(breaks = df1$Type, values = colors) +
#  scale_colour_manual(values = colors, 0.5) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

setEPS()
postscript("metafig-percent.eps",width=13,height=4)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
#  scale_colour_manual(breaks = df1$Type, values = colors) +
#  scale_colour_manual(values = colors, 0.5) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

##########################################
#Fig2 E
df<-table(B1@active.ident,B1$stim)

df<-df[,c("AAA1","AAA2","AG2","AG4","ApoE4","ApoE5","IG2","IG4","Ldlr2","Ldlr3","IH2","IH4")]

df<-prop.table(df,2)

AS1<-data.frame(apply(df[,c(5,9)],1,mean))
AS2<-data.frame(apply(df[,c(6,10)],1,mean))
colnames(AS1)<-"AS1"
colnames(AS2)<-"AS2"
df<-cbind(df[,1:ncol(df)],AS1=AS1$AS1,AS2=AS2$AS2)
df<-df[,c(-5,-6,-9,-10)]

write.table(df,"percent.txt",sep="\t",quote=F)
#perl deal_plot.pl percent.txt > percent.trans.txt
per<-read.table("percent.trans.txt",header=T,sep="\t")

per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
pdf("time.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y = percent, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Percent") +
  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = percent, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Percent") +
  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = percent, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Percent") +
  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

### bar
pdf("time-barplot.pdf",width=18,height=3)
ggplot(data = per,aes(x = batch, y = percent,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-barplot.jpeg",width=900,height=300)
ggplot(data = per,aes(x = batch, y = percent,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-barplot.eps",width=18,height=3)
ggplot(data = per,aes(x = batch, y = percent,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

############ F heatmap and featureplot ####
setwd("C:\\addIH\\Fig4")
B<-readRDS("B.rds")
setwd("new")

new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)

Idents(B) <- "stim"
B1<-subset(B,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(B1) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B1)
B1 <- RenameIdents(B1, new.cluster.ids)

library(pheatmap)
library(ComplexHeatmap)

B.markers<-readRDS("B.markers.rds")

top20 <- B.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

cts <- GetAssayData(B1, slot = 'scale.data')
group = data.frame(group =B1@active.ident,row.names = names(B1@active.ident))
#saveRDS(cts,"cts.rds")
#saveRDS(group,group.rds)

#new_cluster <- sort(scRNA1@active.ident)
#ctstop <- as.matrix(cts[top20$gene, colnames(cts)])
ctstop <- as.matrix(cts[top20$gene, colnames(B1)])

#ctstop_bak<-ctstop

#select 200 genes 
types<-c("Memory B","Antigen-Presenting B","Inflammatory B","DZ B","LZ B","IG-specific B","Plasma B")
cells<-vector()
for (i in 1:length(types)){
set.seed(1234)
index<-which(group$group==types[i])
select<-sample(index,35)
ids<-rownames(group)[select]
cells<-c(cells,ids)
}
tmp2<-ctstop[,cells]
gindex<-which(rownames(group) %in% cells)
group2<-data.frame(row.names=rownames(group)[gindex],group=group[gindex,])
tmp2[tmp2>6] <- 6
tmp2[tmp2< -5] <- -5

pdf("heatmap.pdf")
pheatmap(tmp2,cluster_row=F,cluster_col=F,annotation_col=group2, color = colorRampPalette(colors = c("#136BA5", "#ffffff","#BF172A"))(100),show_rownames=T,show_colnames=F,fontsize_row = 4,fontsize_col = 8,border=F)
dev.off()
setEPS()
postscript("heatmap.eps")
pheatmap(tmp2,cluster_row=F,cluster_col=F,annotation_col=group2, color = colorRampPalette(colors = c("#136BA5", "#ffffff","#BF172A"))(1000),show_rownames=T,show_colnames=F,fontsize_row = 4,fontsize_col = 8,border=F)
dev.off()

types<-c("Memory B","Antigen-Presenting B","Inflammatory B","DZ B","LZ B","IG-specific B","Plasma B")
for (i in 1:length(types)){
	index<-which(top20$cluster==types[i])
	select<-index[1:2]
	features <- top20$gene[select]
	out<-paste("FeaturePlot",types[i],"pdf",sep=".")
	pdf(out,width=10,height=4)
	p<-FeaturePlot(B1, features = features,reduction="tsne")
	print(p)
	dev.off()
	out1<-paste("FeaturePlot",types[i],"eps",sep=".")
	setEPS()
	postscript(out1,width=10,height=4)
	p1<-FeaturePlot(B1, features = features,reduction="tsne")
	print(p1)
	dev.off()
}


####### 细胞的增殖能力打分（D）；T细胞活化打分（C）；免疫细胞受体重排（A）###
#################
setwd("C:\\addIH\\Fig4\\new")
B<-readRDS("B1.rds")

library(clusterProfiler)
library(biomaRt)
library("org.Mm.eg.db")  

mmu_kegg <- clusterProfiler::download_KEGG("mmu")
names(mmu_kegg)
head(mmu_kegg$KEGGPATHID2NAME)
head(mmu_kegg$KEGGPATHID2EXTID)

PATH2ID <- mmu_kegg$KEGGPATHID2EXTID
PATH2NAME <- mmu_kegg$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")

### 增殖 ###
g<-c("Angpt1","Angpt2","Csf1","Csf1r","Efna5","Egfr","Epha2","Fgf1","Fgfr3","Vegfd","Flt1","Flt3","Flt4","Gnb4","Gngt2","Hgf","Igf1","Igf1r","Kdr","Kit","Lat","Kitl","Mras","Pdgfa","Pdgfb","Pdgfrb","Pla2g4a","Rgl1","Tek","Rapgef5","Tgfa","Vegfa","Rasal2","Rasgrp4","Rasgrp3","Rasa4","Ralb","Gng11","Plce1","Calml4","Pla1a")
### 活化 
g<-c("Ctla4","Cd28","Cd3d","Cd3g","Cd4","Cd8a","Cd8b1","Csf2","Fos","Ptpn6","Ifng","Il10","Il4","Jun","Lat","Nfkbia","Pdcd1","Tnf","Mapk13","Icos")
#### 重排
re<- c("mmu04662")
pathway<-PATH_ID_NAME[PATH_ID_NAME$KEGGID %in% re,] ### all pathway gene
SYMBOL<- bitr(pathway$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(B@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(B@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B@meta.data, B@reductions$tsne@cell.embeddings,Type=B@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

###early
Idents(B) <- "stim"
early<-subset(B,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT2","WT3","IH2"))

Idents(B) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)

Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)

### 增殖 ###
g<-c("Angpt1","Angpt2","Csf1","Csf1r","Efna5","Egfr","Epha2","Fgf1","Fgfr3","Vegfd","Flt1","Flt3","Flt4","Gnb4","Gngt2","Hgf","Igf1","Igf1r","Kdr","Kit","Lat","Kitl","Mras","Pdgfa","Pdgfb","Pdgfrb","Pla2g4a","Rgl1","Tek","Rapgef5","Tgfa","Vegfa","Rasal2","Rasgrp4","Rasgrp3","Rasa4","Ralb","Gng11","Plce1","Calml4","Pla1a")
### 活化 
g<-c("Ctla4","Cd28","Cd3d","Cd3g","Cd4","Cd8a","Cd8b1","Csf2","Fos","Ptpn6","Ifng","Il10","Il4","Jun","Lat","Nfkbia","Pdcd1","Tnf","Mapk13","Icos")
#### 重排
re<- c("mmu04662")
pathway<-PATH_ID_NAME[PATH_ID_NAME$KEGGID %in% re,] ### all pathway gene
SYMBOL<- bitr(pathway$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(early@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(early@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(early@meta.data, early@reductions$tsne@cell.embeddings,Type=early@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

###late
#laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
Idents(B) <- "stim"
late<-subset(B,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT2","WT3","IH4"))

Idents(B) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B","IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)

Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("Memory B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B","DZ B", "Antigen-Presenting B","LZ B","LZ B","IG-specific B","Plasma B")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)

### 增殖 ###
g<-c("Angpt1","Angpt2","Csf1","Csf1r","Efna5","Egfr","Epha2","Fgf1","Fgfr3","Vegfd","Flt1","Flt3","Flt4","Gnb4","Gngt2","Hgf","Igf1","Igf1r","Kdr","Kit","Lat","Kitl","Mras","Pdgfa","Pdgfb","Pdgfrb","Pla2g4a","Rgl1","Tek","Rapgef5","Tgfa","Vegfa","Rasal2","Rasgrp4","Rasgrp3","Rasa4","Ralb","Gng11","Plce1","Calml4","Pla1a")
### 活化 
g<-c("Ctla4","Cd28","Cd3d","Cd3g","Cd4","Cd8a","Cd8b1","Csf2","Fos","Ptpn6","Ifng","Il10","Il4","Jun","Lat","Nfkbia","Pdcd1","Tnf","Mapk13","Icos")
#### 重排
re<- c("mmu04662")
pathway<-PATH_ID_NAME[PATH_ID_NAME$KEGGID %in% re,] ### all pathway gene
SYMBOL<- bitr(pathway$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(late@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(late@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(late@meta.data, late@reductions$tsne@cell.embeddings,Type=late@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))


########################
#### 1 热图 ####
library(ComplexHeatmap)
library(scales)
show_col(hue_pal()(11))
hue_pal()(11)
### all
df2<-aggregate(Sig.score~Type+batch,data=data,mean)
df2$batch<-factor(df2$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#df2$Type<-factor(df2$Type,levels=c("Stromal cells","Lymphocyte","Myeloid cells"))
dat<-reshape2::dcast(df2,Type~batch)
rownames(dat)<-dat$Type
dat<-dat[,-1]
dat[is.na(dat)]<-0

col_annotation <- colSums(dat)
row_annotation <- rowSums(dat)

col = HeatmapAnnotation(barplot = anno_barplot(col_annotation, height = unit(3, "cm"),box_width = 0.9, outline = FALSE, gp = gpar(fill = c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00"))),border = F)
row = HeatmapAnnotation(barplot = anno_barplot(row_annotation, height = unit(3, "cm"),box_width = 0.9, outline = FALSE, gp = gpar(fill = c("green", "blue", "orange"))), which = "row",border = F)

pdf("heatmap.all.pdf",width=8,height=6)
Heatmap(dat,name = "Score",cluster_rows = F,cluster_columns=F,top_annotation = col,right_annotation = row,na_col = "grey",row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 16),column_names_rot = 0,row_names_side = "left")
dev.off()
setEPS()
postscript("heatmap.all.eps",width=8,height=6)
Heatmap(dat,name = "Score",cluster_rows = F,cluster_columns=F,top_annotation = col,right_annotation = row,na_col = "grey",row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 16),column_names_rot = 0,row_names_side = "left")
dev.off()


########## iTALK ########
identity = data.frame(cellID = names(B@active.ident),cellType =B@active.ident) 
identity$type<-type
write.table(identity,"B-cellType.txt",quote=F,sep="\t",row.names=F)

setwd("C:\\addIH\\new")
scRNA<-readRDS("scRNA1.rds")
Idents(scRNA) <- "seurat_clusters"
new.cluster.ids <- c("MSC", "SMC","SMC", "Macrophage", "Macrophage", "MSC", 
                       "NK", "EC", "Monocyte", "T cell", "B cell",
                       "Macrophage", "DC", "T cell", "Macrophage", "Macrophage",
                       "MSC", "MSC", "T cell", "B cell", "Macrophage",
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell","SMC")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

identity = data.frame(cellID = names(scRNA@active.ident),cellType =scRNA@active.ident) 
type<-ifelse((identity$cellType == "Macrophage" | identity$cellType == "Monocyte" | identity$cellType == "DC"),"Myeloid cells",ifelse((identity$cellType == "NK" | identity$cellType == "T cell" |identity$cellType == "B cell" | identity$cellType == "pDC"),"Lymphocyte","Stromal cells"))
identity$type<-type

newtype<-read.table("B-cellType.txt",header = T,sep="\t")
newtype$batch<-scRNA$batch

scRNA$newtype<-newtype$newtype

##### early
Idents(scRNA) <- "stim"
early<-subset(scRNA,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","IH2"))

Idents(scRNA) <- "seurat_clusters"
new.cluster.ids <- c("MSC", "SMC","SMC", "Macrophage", "Macrophage", "MSC", 
                       "NK", "EC", "Monocyte", "T cell", "B cell",
                       "Macrophage", "DC", "T cell", "Macrophage", "Macrophage",
                       "MSC", "MSC", "T cell", "B cell", "Macrophage",
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell","SMC")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("MSC", "SMC","SMC", "Macrophage", "Macrophage", "MSC", 
                       "NK", "EC", "Monocyte", "T cell", "B cell",
                       "Macrophage", "DC", "T cell", "Macrophage", "Macrophage",
                       "MSC", "MSC", "T cell", "B cell", "Macrophage",
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell","SMC")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)

identity = data.frame(group =early@active.ident, row.names = names(early@active.ident)) # 
type<-ifelse((identity$group == "Macrophage" | identity$group == "Monocyte" | identity$group == "DC"),"Myeloid cells",ifelse((identity$group == "NK" | identity$group == "T cell" |identity$group == "B cell" | identity$group == "pDC"),"Lymphocyte","Stromal cells"))
identity$type<-type
early$type<-identity$type

group<-paste(early$newtype,early$batch,early$type,sep='-')
early$group<-group

iTalk_data <- as.data.frame(t(early@assays$RNA@counts))
iTalk_data$cell_type <- early$group
iTalk_data$compare_group <- early$batch

### late
Idents(scRNA) <- "stim"
late<-subset(scRNA,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT2","IH4"))

Idents(scRNA) <- "seurat_clusters"
new.cluster.ids <- c("MSC", "SMC","SMC", "Macrophage", "Macrophage", "MSC", 
                       "NK", "EC", "Monocyte", "T cell", "B cell",
                       "Macrophage", "DC", "T cell", "Macrophage", "Macrophage",
                       "MSC", "MSC", "T cell", "B cell", "Macrophage",
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell","SMC")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("MSC", "SMC","SMC", "Macrophage", "Macrophage", "MSC", 
                       "NK", "EC", "Monocyte", "T cell", "B cell",
                       "Macrophage", "DC", "T cell", "Macrophage", "Macrophage",
                       "MSC", "MSC", "T cell", "B cell", "Macrophage",
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell","SMC")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)

identity = data.frame(group =late@active.ident, row.names = names(late@active.ident)) # 
type<-ifelse((identity$group == "Macrophage" | identity$group == "Monocyte" | identity$group == "DC"),"Myeloid cells",ifelse((identity$group == "NK" | identity$group == "T cell" |identity$group == "B cell" | identity$group == "pDC"),"Lymphocyte","Stromal cells"))
identity$type<-type

late$type<-identity$type
group<-paste(late$newtype,late$batch,late$type,sep='-')
late$group<-group

iTalk_data <- as.data.frame(t(late@assays$RNA@counts))
iTalk_data$cell_type <- late$group
iTalk_data$compare_group <- late$batch


#### iTALK ###
library(biomaRt)
library(iTALK)
library(homologene)
genes<-colnames(iTalk_data)
###人鼠同源基因转化
#genes<-read.table("genes.txt")
human=useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse=useMart("ensembl",dataset="mmusculus_gene_ensembl")

genesV2=getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=genes, mart=mouse, attributesL=c("hgnc_symbol"), martL=human, uniqueRows=TRUE)

genesV2_dup=genesV2[!duplicated(genesV2[,2]),]
genesV2_dup=genesV2_dup[!duplicated(genesV2_dup[,1]),]

type<-c("cell_type","cell_type")
group<-c("compare_group","compare_group")

genesV2_dup<-rbind(genesV2_dup,type,group)
##########
aa<-which(colnames(iTalk_data) %in% genesV2_dup[,1]) 
bb<-iTalk_data[,aa]
bbn<-colnames(bb)
gsn<-vector()

for (i in 1:length(bbn)){
	index<-which(unique(genesV2_dup[,1]) %in% bbn[i])
	gsn<-c(gsn,genesV2_dup[index,2])
}
colnames(bb)<-gsn

iTalk_data<-bb

my10colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282',"#f7aa5d")
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)

iTalk_res1 <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res1 <- rbind(iTalk_res1, res_cat)
}
#write.table(iTalk_res1,"early.iTalk_res1.xls",sep="\t",quote=F,row.names = F)
write.table(iTalk_res1,"late.iTalk_res1.xls",sep="\t",quote=F,row.names = F)



#################### GO #####################
library(clusterProfiler)
#library(biomaRt)
library("org.Mm.eg.db")  

logFCfilter=0.25
adjPvalFilter=0.25
B.markers <- FindAllMarkers(object = B,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
							   
B.markers<-readRDS("B.markers.rds")
group_g <- B.markers[,c(6,7)]
colnames(group_g) <- c("group","gene")
#library(clusterProfiler)
# Convert gene ID into entrez genes
tmp <- bitr(group_g$gene, fromType="ALIAS", 
            toType="ENTREZID", 
            OrgDb="org.Mm.eg.db")
de_gene_clusters=merge(tmp,group_g,by.x='ALIAS',by.y='gene')
table(de_gene_clusters$group)
# Run full GO enrichment test
formula_res <- compareCluster(
  ENTREZID~group, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.10,
  qvalueCutoff  = 0.05
)
# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)
# Plot both analysis results
jpeg('compared_GO_term_DE_cluster.jpg',width = 1200,height = 900)
dotplot(formula_res, showCategory=10)
dev.off()
#jpeg('compared_GO_term_DE_cluster_simplified.jpg',width = 1000,height = 600)
pdf("compared_GO_term_DE_cluster_simplified.pdf",width=18,height=20)
p<-dotplot(lineage1_ego, showCategory=15)
p+theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
dev.off()	

p<-dotplot(lineage1_ego, showCategory=15)
p+theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
ggsave("compared_GO_term_DE_cluster_simplified.eps", width = 25, height = 20, units = "cm")

write.table(formula_res@compareClusterResult,file="GO.txt",sep="\t",quote=F,row.names = F)
write.table(lineage1_ego@compareClusterResult,file="GO_simplified.txt",sep="\t",quote=F,row.names = F)

########### zscore B细胞相关
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
Idents(B) <- "stim"
B1<-subset(B,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(B) <- "seurat_clusters"
new.cluster.ids <- c("Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", 
                     "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B",
                     "DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B",
                     "IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

Idents(B1) <- "seurat_clusters"
new.cluster.ids <- c("Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", 
                     "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B",
                     "DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B",
                     "IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B1)
B1 <- RenameIdents(B1, new.cluster.ids)
#B1@active.ident=factor(B1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


k<-B1
g<-c("Ccr5","Cxcr3","Csf2","Igsf6","Il2ra","Cd38","Cd40","Cd5","Ms4a1")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=500,height=300)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=5,height=3)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=5,height=3)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

data<-cbind(data,Sig.score)
#data <- data %>% group_by(batch) %>% mutate(med = median(Sig.score))
#median.batch.sig<-aggregate(Sig.score~batch,data=data,median)

data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
jpeg("batch-score2.jpg",width=800,height=600)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+
#	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("batch-score2.pdf",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("batch-score2.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score2.jpg",width=800,height=700)
ggplot(data,aes(x=B1@active.ident,y=Sig.score,fill=B1@active.ident,color=B1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
ggplot(data,aes(x=B1@active.ident,y=Sig.score,fill=B1@active.ident,color=B1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
ggplot(data,aes(x=B1@active.ident,y=Sig.score,fill=B1@active.ident,color=B1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### early
#data.bak<-data
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}
data<-data[earlyindex,]
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
jpeg("batch-score2-early.jpg",width=800,height=600)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+
#	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("batch-score2-early.pdf",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("batch-score2-early.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score2-early.jpg",width=800,height=700)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### late
data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}
data<-data[lateindex,]
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
jpeg("batch-score2-late.jpg",width=800,height=600)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+
#	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("batch-score2-late.pdf",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("batch-score2-late.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score2-late.jpg",width=800,height=700)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",		"CD16 hi Monocyte","ISG hi Monocyte"))

pdf("time-zscore-Brela.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-Brela.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-Brela.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

### bar ##
pdf("time-zscore-Brela.pdf",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-Brela.jpeg",width=900,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-Brela.eps",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()


########### zscore 控制淋巴细胞转运的基因
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
Idents(B) <- "stim"
B1<-subset(B,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(B) <- "seurat_clusters"
new.cluster.ids <- c("Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", 
                     "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B",
                     "DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B",
                     "IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B)
B <- RenameIdents(B, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

Idents(B1) <- "seurat_clusters"
new.cluster.ids <- c("Antigen-Presenting B", "Antigen-Presenting B", "Antigen-Presenting B", 
                     "Antigen-Presenting B", "Antigen-Presenting B", "Inflammatory B", "Antigen-Presenting B",
                     "DZ B", "Antigen-Presenting B", "Inflammatory B","LZ B","LZ B",
                     "IG-specific B","Plasma B","Inflammatory B")
names(new.cluster.ids) <- levels(B1)
B1 <- RenameIdents(B1, new.cluster.ids)
#B1@active.ident=factor(B1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

k<-B1
g<-c('Glycam1', 'Fut7', 'Gcnt1', 'Chst4', 'B3gnt3', 'Ccl21a')
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=500,height=300)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=5,height=3)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=5,height=3)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

data<-cbind(data,Sig.score)
#data <- data %>% group_by(batch) %>% mutate(med = median(Sig.score))
#median.batch.sig<-aggregate(Sig.score~batch,data=data,median)

data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
jpeg("batch-score.jpg",width=800,height=600)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+
#	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("batch-score.pdf",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("batch-score.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score.jpg",width=800,height=700)
ggplot(data,aes(x=B1@active.ident,y=Sig.score,fill=B1@active.ident,color=B1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score.pdf",width=10,height=10)
ggplot(data,aes(x=B1@active.ident,y=Sig.score,fill=B1@active.ident,color=B1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
ggplot(data,aes(x=B1@active.ident,y=Sig.score,fill=B1@active.ident,color=B1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### early
#data.bak<-data
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}
data<-data[earlyindex,]
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
jpeg("batch-score-early.jpg",width=800,height=600)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+
#	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("batch-score-early.pdf",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("batch-score-early.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score-early.jpg",width=800,height=700)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### late
data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}
data<-data[lateindex,]
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
jpeg("batch-score-late.jpg",width=800,height=600)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+
#	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("batch-score-late.pdf",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("batch-score-late.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score-late.jpg",width=800,height=700)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Antigen-Presenting B","Inflammatory B"),c("Antigen-Presenting B","DZ B"),c("Antigen-Presenting B","LZ B"),c("Antigen-Presenting B","IG-specific B"),c("Antigen-Presenting B","Plasma B"),c("Inflammatory B","DZ B"),c("Inflammatory B","LZ B"),c("Inflammatory B","IG-specific B"),c("Inflammatory B","Plasma B"),c("DZ B","LZ B"),c("DZ B","IG-specific B"),c("DZ B","Plasma B"),c("LZ B","IG-specific B"),c("LZ B","Plasma B"),c("IG-specific B","Plasma B")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",		"CD16 hi Monocyte","ISG hi Monocyte"))

pdf("time-zscore-Recruit.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-Recruit.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-Recruit.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

### bar ##
pdf("time-zscore-Recruit.pdf",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-Recruit.jpeg",width=900,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-Recruit.eps",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()


############# Fig4 K
library(clusterProfiler)
#library(biomaRt)
library("org.Mm.eg.db")  

mmu_kegg <- clusterProfiler::download_KEGG("mmu")
names(mmu_kegg)
head(mmu_kegg$KEGGPATHID2NAME)
head(mmu_kegg$KEGGPATHID2EXTID)

PATH2ID <- mmu_kegg$KEGGPATHID2EXTID
PATH2NAME <- mmu_kegg$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")

#innate immune
mmu04650<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04650",]
SYMBOL<- bitr(mmu04650$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-B1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

### early
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}

### late
data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata-immune.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-immune2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-B.immune.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-B.immune.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-B.immune.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#抗原提呈加工
mmu04612<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04612",]
SYMBOL<- bitr(mmu04612$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-B1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

### early
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}

### late
#data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata-APC.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-APC2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-B.APC.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-B.APC.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-B.APC.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#浆细胞产生抗体
mmu04672<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04672",]
SYMBOL<- bitr(mmu04672$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-B1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

### early
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}

### late
#data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata-plasma.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-plasma2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-B.plasma.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-B.plasma.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-B.plasma.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#细胞生长、增殖
mmu04014<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04014",]
SYMBOL04014<- bitr(mmu04014$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

mmu04310<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04310",]
SYMBOL04310<- bitr(mmu04310$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

proliferation<-rbind(SYMBOL04014,SYMBOL04310)
SYMBOL<-unique(proliferation$SYMBOL)

k<-B1
g<-SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

### early
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}

### late
#data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata-proliferation.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-proliferation2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-B.proliferation.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-B.proliferation.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-B.proliferation.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#Chemokine signaling pathway
mmu04062<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04062",]
SYMBOL<- bitr(mmu04062$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-B1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

### early
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}

### late
#data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata-Chemokine.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-Chemokine2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-B.Chemokine.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-B.Chemokine.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-B.Chemokine.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#Cytokine signaling pathway
mmu04061<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04061",]
SYMBOL<- bitr(mmu04061$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-B1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(B1@meta.data, B1@reductions$tsne@cell.embeddings,Type=B1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

### early
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}

### late
#data<-data.bak
laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
lateindex<-vector()
for (i in 1:length(laten)){
	aa<-which(data$stim==laten[i])
	lateindex<-c(lateindex,aa)
}

################ Fig2 H
early<-data[earlyindex,]
late<-data[lateindex,]

earlydata<-aggregate(Sig.score~batch+Type,data=early,mean)
latedata<-aggregate(Sig.score~batch+Type,data=late,mean)

group<-rep("early",nrow(earlydata))
earlydata<-cbind(earlydata,group)

group<-rep("late",nrow(latedata))
latedata<-cbind(latedata,group)

newdata<-rbind(earlydata,latedata)
write.table(newdata,"newdata-Cytokine.txt",sep="\t",quote=F,row.names=F)
#perl sort.sig.time.pl newdata-Cytokine.txt > newdata-Cytokine2.txt

per<-read.table("newdata-Cytokine2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-B.Cytokine.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-B.Cytokine.jpeg",width=900,height=300)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
setEPS()
postscript("time-zscore-B.Cytokine.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()
