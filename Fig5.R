library(Seurat)
library(ggplot2)
library(data.table)
library(hdf5r)
library(dplyr)
library("ggthemes")
library("RColorCrewer")
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)

setwd("C:\\addIH\\Fig1")
combined<-readRDS("fig1.rds")
combined<-subset(combined,idents = 24,invert=T)
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)

#cDC<-readRDS("DC.rds")
#pDC<-readRDS("pDC.rds")
#C<-merge(cDC,pDC)

C<-subset(combined,idents=c("cDC","pDC"))
setwd("C:\\addIH\\Fig5")
C <- FindVariableFeatures(C, selection.method = "vst", nfeatures = 2000) 
al.genes<-row.names(C) 
C <- ScaleData(C, features = al.genes) 
C <- RunPCA(C, features = VariableFeatures(object = C), npcs = 50)   
C <- FindNeighbors(C, reduction = "pca", dims = 1:30) 
C <- FindClusters(C, resolution = 0.2)
C <- RunTSNE(C,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(C, reduction = "tsne",label="T",label.size=3.5,pt.size = 1)
C <- RunUMAP(C,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(C, reduction = "umap",label="T",label.size=3.5,pt.size = 1)

C.markers <- FindAllMarkers(C, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
saveRDS(C.markers,"C.markers.rds")
top10 <- C.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(C, features = top10$gene) + NoLegend()
dev.off()

###############################
#平均表达谱函数AverageExpression
AverageExp<-AverageExpression(C,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pdf("corr.pdf",width=20,height=20)
pheatmap(coorda$r)
dev.off()

saveRDS(C,"C.rds")

C<-readRDS("C.rds")
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
my36colors <- c('#E5D2DD', '#53A85F', '#F1CC72', '#F3C1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AC3282', '#23452F', '#CD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DECA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0CE', '#C53E2C',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

#markers 
#pDC
#StackedVlnPlot(C,  c('Klk1', 'Cox6a2', 'Siglech', 'Cd7', 'Ccr9', "Rnase6", "Irf8", "Bst2", "Cd8b1", "Tcf4", "Dnajc7", "D13Ertd608e"), pt.size=1, cols=my36colors)
StackedVlnPlot(C,  c('Klk1', 'Cox6a2', 'Siglech', 'Cd7', 'Ccr9'), pt.size=1, cols=my36colors)
StackedVlnPlot(C,  c("Rnase6", "Irf8", "Bst2", "Cd8b1", "Tcf4"), pt.size=1, cols=my36colors)
StackedVlnPlot(C,  c("Dnajc7", "D13Ertd608e"), pt.size=1, cols=my36colors)

#CD209+ DC
StackedVlnPlot(C,  c('Cd209a', 'Mgl2', 'Clec10a', 'Ckb'), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c('H2afz', "Napsa","Syngr2"), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c("Cst3", "Fcgrt", "Stmn1", "H2-DMb2", "Tuba1a"), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c('Xcr1', 'Irf8', 'Naaa', 'Ifi30', 'Ifitm1', 'Cd74'), pt.size=0, cols=my36colors)

#CCR7+ DC
#StackedVlnPlot(C,  c('Ccr7', 'Cacnb3', 'Ramp3', 'Ccl22', 'Ccl17', "Serpinb6b", "ll4i1", "Cd209a", "Tmem123", "Mgl2", "Ccl5", "Fscn1"), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c('Ccr7', 'Cacnb3', 'Ramp3', 'Ccl22'), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c('Ccl17', "Serpinb6b", "Il4i1", "Cd209a"), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c("Tmem123", "Mgl2", "Ccl5", "Fscn1"), pt.size=0, cols=my36colors)

#XCR1+ DC
#StackedVlnPlot(C,  c('Sep3', 'Xcr1', 'Wdfy4', 'Naga', 'Ckb', "Naaa", "Pldb1", "Cst3", "Ppt1", "lrf8", "Eef1b2", "ld2"), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c('Xcr1', 'Wdfy4', 'Naga'), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c('Ckb', "Naaa","Cst3"), pt.size=0, cols=my36colors)
StackedVlnPlot(C,  c("Ppt1", "Irf8", "Eef1b2", "Id2"), pt.size=0, cols=my36colors)


#C.bak<-C
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(C)
C <- RenameIdents(C, new.cluster.ids)

C.markers <- FindAllMarkers(C, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- C.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.cell.pdf",width=25,height=30)
DoHeatmap(C, features = top10$gene) + NoLegend()
dev.off()


#################### 挑选样品分析比例
Idents(C) <- "stim"
C1<-subset(C,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(C) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")

names(new.cluster.ids) <- levels(C)
C <- RenameIdents(C, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(C1) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")

names(new.cluster.ids) <- levels(C1)
C1 <- RenameIdents(C1, new.cluster.ids)
#C1@active.ident=factor(C1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))

pdf("DimPlot.tsne.pdf",width=8,height=5)
DimPlot(C1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.tsne.jpg",width=800,height=500)
DimPlot(C1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.tsne.eps",width=8,height=5)
DimPlot(C1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

pdf("DimPlot.umap.pdf",width=10,height=5)
DimPlot(C1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.umap.jpg",width=1000,height=500)
DimPlot(C1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.umap.eps",width=10,height=5)
DimPlot(C1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

C1$batch<-factor(C1$batch,levels = c("WT","AS","AAA","IH","IG","AG"))
pdf("DimPlot.umap-batch.pdf",width=18,height=3)
DimPlot(C1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.umap-batch.jpg",width=1800,height=300)
DimPlot(C1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.umap-batch.eps",width=18,height=3)
DimPlot(C1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

pdf("DimPlot.tsne-batch.pdf",width=18,height=3)
DimPlot(C1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.tsne-batch.jpg",width=1800,height=300)
DimPlot(C1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.tsne-batch.eps",width=18,height=3)
DimPlot(C1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()


### percent
df1<-table(C1@active.ident,C1$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Catch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#df1$Type<-factor(df1$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
mycolors<-c("#6495ED","#FFA500","#228C22","#FF4500")

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
ggplot(data = df1, aes(x = Catch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Catch", y = "Percent") +
# scale_colour_manual(breaks = df1$Type, values = colors) +
#  scale_colour_manual(values = colors, 0.5) +
  theme_bw() +  theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("metafig-percent.jpg",width=1800,height=300)
ggplot(data = df1, aes(x = Catch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Catch", y = "Percent") +
#  scale_colour_manual(breaks = df1$Type, values = colors) +
#  scale_colour_manual(values = colors, 0.5) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

setEPS()
postscript("metafig-percent.eps",width=13,height=4)
ggplot(data = df1, aes(x = Catch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Catch", y = "Percent") +
#  scale_colour_manual(breaks = df1$Type, values = colors) +
#  scale_colour_manual(values = colors, 0.5) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

##########################################
#Fig2 E
df<-table(C1@active.ident,C1$stim)

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
pdf("time-barplot.pdf",width=8,height=3)
ggplot(data = per,aes(x = batch, y = percent,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-barplot.jpeg",width=600,height=300)
ggplot(data = per,aes(x = batch, y = percent,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-barplot.eps",width=8,height=3)
ggplot(data = per,aes(x = batch, y = percent,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()


#################### GO #####################
library(clusterProfiler)
#library(biomaRt)
library("org.Mm.eg.db")  

logFCfilter=0.25
adjPvalFilter=0.25
C.markers <- FindAllMarkers(object = C,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
							   
C.markers<-readRDS("C.markers.cells.rds")
group_g <- C.markers[,c(6,7)]
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

####### 细胞的增殖能力打分（D）；T细胞活化打分（C）；免疫细胞受体重排（A）###
#################
setwd("C:\\addIH\\Fig5")
C<-readRDS("C.rds")
Idents(C) <- "stim"
C1<-subset(C,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(C) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")

names(new.cluster.ids) <- levels(C)
C <- RenameIdents(C, new.cluster.ids)

Idents(C1) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")

names(new.cluster.ids) <- levels(C1)
C1 <- RenameIdents(C1, new.cluster.ids)

########## iTALK ########
identity = data.frame(cellID = names(C1@active.ident),cellType =C1@active.ident) 
identity$type<-type
write.table(identity,"C-cellType.txt",quote=F,sep="\t",row.names=F)

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

newtype<-read.table("all.newtype-C.txt",header = T,sep="\t")
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
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell")
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











################################################
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

### 抗原提呈 ###
g<-c("Ciita","Calr","Cd4","Cd8a","Cd8b1","Ctsb","Ctsl","Ctss","H2-Aa","H2-Ab1","H2-Eb1","H2-DMa","H2-DMb1","H2-DMb2","H2-Oa","H2-Ob","Hspa8","Hspa1b","Hsp90ab1","Hsp90aa1","Ifng","Cd74","Klrc1","Klrd1","Lgmn","Hspa1a","Tnf")
gindex<-which(g %in% rownames(C1@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(C1@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(C1@meta.data, C1@reductions$tsne@cell.embeddings,Type=C1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

data<-cbind(data,Sig.score)
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

###early
Idents(C1) <- "stim"
early<-subset(C1,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT2","WT3","IH2"))

Idents(C1) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(C1)
C1 <- RenameIdents(C1, new.cluster.ids)

Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)

### 抗原提呈 ###
g<-c("Ciita","Calr","Cd4","Cd8a","Cd8b1","Ctsb","Ctsl","Ctss","H2-Aa","H2-Ab1","H2-Eb1","H2-DMa","H2-DMb1","H2-DMb2","H2-Oa","H2-Ob","Hspa8","Hspa1b","Hsp90ab1","Hsp90aa1","Ifng","Cd74","Klrc1","Klrd1","Lgmn","Hspa1a","Tnf")

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
Idents(C1) <- "stim"
late<-subset(C1,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT2","WT3","IH4"))

Idents(C1) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(C1)
C1 <- RenameIdents(C1, new.cluster.ids)

Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)

### 抗原提呈 ###
g<-c("Ciita","Calr","Cd4","Cd8a","Cd8b1","Ctsb","Ctsl","Ctss","H2-Aa","H2-Ab1","H2-Eb1","H2-DMa","H2-DMb1","H2-DMb2","H2-Oa","H2-Ob","Hspa8","Hspa1b","Hsp90ab1","Hsp90aa1","Ifng","Cd74","Klrc1","Klrd1","Lgmn","Hspa1a","Tnf")
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



########### zscore 与三级淋巴组织形成有关的细胞趋化因子
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
Idents(C) <- "stim"
C1<-subset(C,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(C) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(C)
C <- RenameIdents(C, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

Idents(C1) <- "seurat_clusters"
new.cluster.ids <- c("CD209a+ DC","CD209a+ DC","CD209a+ DC","pDC","Ccr7+ DC","Ccr7+ DC","CD209a+ DC","CD209a+ DC","Ccr7+ DC","CD209a+ DC","pDC","pDC")
names(new.cluster.ids) <- levels(C1)
C1 <- RenameIdents(C1, new.cluster.ids)
#C1@active.ident=factor(C1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


k<-C1
g<-c("Ccl2","Ccl3","Ccl4","Ccl5","Ccl8","Ccl18","Ccl19","Ccl21","Cxcl9","Cxcl10","Cxcl11","Cxcl13","Ccr7","Cxcr5","Sell","Lamp3")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(C1@meta.data, C1@reductions$tsne@cell.embeddings,Type=C1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=550,height=300)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=6,height=3.5)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=6,height=3.5)
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
ggplot(data,aes(x=C1@active.ident,y=Sig.score,fill=C1@active.ident,color=C1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
ggplot(data,aes(x=C1@active.ident,y=Sig.score,fill=C1@active.ident,color=C1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
ggplot(data,aes(x=C1@active.ident,y=Sig.score,fill=C1@active.ident,color=C1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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

pdf("time-zscore-Crela.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-Crela.jpeg",width=900,height=300)
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
postscript("time-zscore-Crela.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

pdf("time-zscore-Crela.barplot.pdf",width=8,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-Crela.barplot.jpeg",width=800,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-Crela.barplot.eps",width=8,height=3)
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
Idents(C) <- "stim"
C1<-subset(C,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT2","WT3","IH2","IH4"))

Idents(C) <- "seurat_clusters"
new.cluster.ids <- c("Antigen-Presenting C", "Antigen-Presenting C", "Antigen-Presenting C", 
                     "Antigen-Presenting C", "Antigen-Presenting C", "Inflammatory C", "Antigen-Presenting C",
                     "DZ C", "Antigen-Presenting C", "Inflammatory C","LZ C","LZ C",
                     "IG-specific C","Plasma C","Inflammatory C")
names(new.cluster.ids) <- levels(C)
C <- RenameIdents(C, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

Idents(C1) <- "seurat_clusters"
new.cluster.ids <- c("Antigen-Presenting C", "Antigen-Presenting C", "Antigen-Presenting C", 
                     "Antigen-Presenting C", "Antigen-Presenting C", "Inflammatory C", "Antigen-Presenting C",
                     "DZ C", "Antigen-Presenting C", "Inflammatory C","LZ C","LZ C",
                     "IG-specific C","Plasma C","Inflammatory C")
names(new.cluster.ids) <- levels(C1)
C1 <- RenameIdents(C1, new.cluster.ids)
#C1@active.ident=factor(C1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

k<-C1
g<-c("Glycam1","Fut7","Gcnt1","Chst4","B3gnt3","Ccl21a")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(C1@meta.data, C1@reductions$tsne@cell.embeddings,Type=C1@active.ident,k)
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
pdf("sig.score.pdf",width=6,height=3.5)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=6,height=3.5)
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
ggplot(data,aes(x=C1@active.ident,y=Sig.score,fill=C1@active.ident,color=C1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
ggplot(data,aes(x=C1@active.ident,y=Sig.score,fill=C1@active.ident,color=C1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
ggplot(data,aes(x=C1@active.ident,y=Sig.score,fill=C1@active.ident,color=C1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("CD209a+ DC","pDC"),c("CD209a+ DC","Ccr7+ DC"),c("pDC","Ccr7+ DC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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

pdf("time-zscore-Recruit.barplot.pdf",width=8,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-Recruit.barplot.jpeg",width=800,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-Recruit.barplot.eps",width=8,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

