library(Seurat)
library(ggplot2)
library(data.table)
library(hdf5r)
library(dplyr)
library("ggthemes")
library("RColorBrewer")

setwd("C:\\addIH\\Fig3")

TNK<-subset(combined,idents=c("T cell","NK"))

TNK <- FindVariableFeatures(TNK, selection.method = "vst", nfeatures = 2000) 
al.genes<-row.names(TNK) 
TNK <- ScaleData(TNK, features = al.genes) 
TNK <- RunPCA(TNK, features = VariableFeatures(object = TNK), npcs = 50)   
TNK <- FindNeighbors(TNK, reduction = "pca", dims = 1:30) 
TNK <- FindClusters(TNK, resolution = 0.6)
TNK <- RunTSNE(TNK,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(TNK, reduction = "tsne",label="T",label.size=3.5,pt.size = 1)
TNK <- RunUMAP(TNK,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(TNK, reduction = "umap",label="T",label.size=3.5,pt.size = 1)

T.markers <- FindAllMarkers(TNK, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- T.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(T, features = top10$gene) + NoLegend()
dev.off()

saveRDS(TNK,"TNK.rds")

TNK<-readRDS("TNK.rds")
###############################
#平均表达谱函数AverageExpression
AverageExp<-AverageExpression(TNK,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pdf("corr.pdf",width=20,height=20)
pheatmap(coorda$r)
dev.off()
##################################
###区分T 和 NK
StackedVlnPlot(TNK, c('Cd3d', 'Cd3g', 'Cd4', 'Cd8a', 'Cd8b1', 'Ncr1'), pt.size=0, cols=my36colors)
#Cd8+ CTL
StackedVlnPlot(TNK, c('Gzmk', 'Gzmb', 'Epsti1', 'Ly6a'), pt.size=0, cols=my36colors)#Gzmk hi CD8+ CTL
StackedVlnPlot(TNK, c('Spp1', 'Ifng', 'Ccl3', 'Xcl1', 'Ccl4'), pt.size=0, cols=my36colors)   #Ifng hi CD8+ CTL
#CTLs
StackedVlnPlot(TNK, c('Pdcd1','Nkg7', 'Klrc1', 'Klrd1', 'Klrk1'), pt.size=0, cols=my36colors)
StackedVlnPlot(TNK, c('Tcf7', 'Sh2d1a', 'Cd40lg', 'Cd69', 'Slamf6'), pt.size=0, cols=my36colors)
#Th17 cell
StackedVlnPlot(TNK, c('Il17a', 'Il17f', 'Il23r', 'Blk', 'Tcf12'), pt.size=0, cols=my36colors)    



new.cluster.ids <- c("Cd8+ T","Cd4+Cd8+ T","Cd4+Cd8+ T","Cd3+ T","Cd8+CTL","Cd8+CTL","Cd4+Cd8+ T","Cd3+ T","Cd8+CTL","Cd8+CTL","Cd8+CTL","Cd8+CTL","Cd3+ T","Cd4+ T","Cd8+ T","Cd3+ T","Cd3+ T","Cd8+CTL","NK","Cd4+Cd8+ T","NK","Cd4+Cd8+ T","NKT","Cd8+CTL","Th17","Cd8+CTL","Cd3+ T","NK","CTL","Cd4+Cd8+ T","Cd3+ T","Cd4+Cd8+ T")
names(new.cluster.ids) <- levels(TNK)
TNK <- RenameIdents(TNK, new.cluster.ids)

T.markers <- FindAllMarkers(TNK, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- T.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(TNK, features = top10$gene) + NoLegend()
dev.off()


#################### 挑选样品分析比例
Idents(TNK) <- "stim"
TNK1<-subset(TNK,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE2","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT1","WT2","WT3","IH2","IH4"))

Idents(TNK) <- "seurat_clusters"
new.cluster.ids <- c("Cd8+ T","Cd4+Cd8+ T","Cd4+Cd8+ T","Cd3+ T","Cd8+CTL","Cd8+CTL","Cd4+Cd8+ T","Cd3+ T","Cd8+CTL","Cd8+CTL","Cd8+CTL","Cd8+CTL","Cd3+ T","Cd4+ T","Cd8+ T","Cd3+ T","Cd3+ T","Cd8+CTL","NK","Cd4+Cd8+ T","NK","Cd4+Cd8+ T","NKT","Cd8+CTL","Th17","Cd8+CTL","Cd3+ T","NK","CTL","Cd4+Cd8+ T","Cd3+ T","Cd4+Cd8+ T")

names(new.cluster.ids) <- levels(TNK)
TNK <- RenameIdents(TNK, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(TNK1) <- "seurat_clusters"
new.cluster.ids <- c("Cd8+ T","Cd4+Cd8+ T","Cd4+Cd8+ T","Cd3+ T","Cd8+CTL","Cd8+CTL","Cd4+Cd8+ T","Cd3+ T","Cd8+CTL","Cd8+CTL","Cd8+CTL","Cd8+CTL","Cd3+ T","Cd4+ T","Cd8+ T","Cd3+ T","Cd3+ T","Cd8+CTL","NK","Cd4+Cd8+ T","NK","Cd4+Cd8+ T","NKT","Cd8+CTL","Th17","Cd8+CTL","Cd3+ T","NK","CTL","Cd4+Cd8+ T","Cd3+ T","Cd4+Cd8+ T")

names(new.cluster.ids) <- levels(TNK1)
TNK1 <- RenameIdents(TNK1, new.cluster.ids)
#TNK1@active.ident=factor(TNK1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))

pdf("DimPlot.tsne.pdf",width=8,height=5)
DimPlot(TNK1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.tsne.jpg",width=800,height=500)
DimPlot(TNK1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.tsne.eps",width=8,height=5)
DimPlot(TNK1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

pdf("DimPlot.umap.pdf",width=10,height=5)
DimPlot(TNK1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.umap.jpg",width=1000,height=500)
DimPlot(TNK1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.umap.eps",width=10,height=5)
DimPlot(TNK1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

TNK1$batch<-factor(TNK1$batch,levels = c("WT","AS","AAA","IH","IG","AG"))
pdf("DimPlot.umap-batch.pdf",width=18,height=3)
DimPlot(TNK1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.umap-batch.jpg",width=1800,height=300)
DimPlot(TNK1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.umap-batch.eps",width=18,height=3)
DimPlot(TNK1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

pdf("DimPlot.tsne-batch.pdf",width=18,height=3)
DimPlot(TNK1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.tsne-batch.jpg",width=1800,height=300)
DimPlot(TNK1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.tsne-batch.eps",width=18,height=3)
DimPlot(TNK1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()


### percent
df1<-table(TNK1@active.ident,TNK1$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#df1$Type<-factor(df1$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("metafig-percent.pdf",width=13,height=4)
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
df<-table(TNK1@active.ident,TNK1$stim)

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

#################### GO #####################
#T.markers<-readRDS("T.markers.rds")
logFCfilter=0.25
adjPvalFilter=0.25
T.markers <- FindAllMarkers(object = T,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
							   

group_g <- T.markers[,c(6,7)]
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

########### zscore T细胞相关
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
Idents(TNK) <- "stim"
TNK1<-subset(TNK,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE2","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT1","WT2","WT3","IH2","IH4"))

Idents(TNK) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(TNK)
TNK <- RenameIdents(TNK, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(TNK1) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(TNK1)
TNK1 <- RenameIdents(TNK1, new.cluster.ids)
TNK1@active.ident=factor(TNK1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))


k<-TNK1
g<-c("Cxcl13","Cd200","Fbln7","Icos","Sgpp2","Sh2d1a","Tigit","Pdcd1","Ccr5","Cxcr3","Csf2","Igsf6","Il2ra","Cd38","Cd40","Cd5","Ms4a1")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=800,height=300)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=8,height=4)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=8,height=4)
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
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score.pdf",width=10,height=10)
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### early
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
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Cd8+ T","Cd4+Cd8+ T"),c("Cd8+ T","Cd3+ T"),c("Cd8+ T","Cd8+CTL"),c("Cd8+ T","Cd4+ T"),c("Cd8+ T","NK"),c("Cd8+ T","NKT"),c("Cd8+ T","Th17"),c("Cd8+ T","CTL"),c("Cd4+Cd8+ T","Cd3+ T"),c("Cd4+Cd8+ T","Cd8+CTL"),c("Cd4+Cd8+ T","Cd4+ T"),c("Cd4+Cd8+ T","NK"),c("Cd4+Cd8+ T","NKT"),c("Cd4+Cd8+ T","Th17"),c("Cd4+Cd8+ T","CTL"),c("Cd3+ T","Cd8+CTL"),c("Cd3+ T","Cd4+ T"),c("Cd3+ T","NK"),c("Cd3+ T","NKT"),c("Cd3+ T","Th17"),c("Cd3+ T","CTL"),c("Cd8+CTL","Cd4+ T"),c("Cd8+CTL","NK"),c("Cd8+CTL","NKT"),c("Cd8+CTL","Th17"),c("Cd8+CTL","CTL"),c("Cd4+ T","NK"),c("Cd4+ T","NKT"),c("Cd4+ T","Th17"),c("Cd4+ T","CTL"),c("NK","NKT"),c("NK","Th17"),c("NK","CTL"),c("NKT","Th17"),c("NKT","CTL"),c("Th17","CTL")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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

per<-read.table("newdata.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-Trela.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-Trela.jpeg",width=900,height=300)
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
postscript("time-zscore-Trela.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

########### Fig3-K TFH
k<-TNK1
g<-c("Cxcr5","Sh2d1a","Il21","Icos","Tcf","Asap1","Nt5e","Rgs10","Tnfsf8","Cxcr4","Cxcr5","Bcl6","S100a11")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=800,height=300)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=8,height=4)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=8,height=4)
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
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score.pdf",width=10,height=10)
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

########### Fig3-K CD3_T细胞向TFH细胞分化的调控基因
k<-TNK1
g<-c("Asap1","Tcf7","Rgs10","Tbc1d4","Gm10275","Folr4","Zfp36l1","Cxcr5","Ddx24","Rps12-ps3","Rps20","Erdr1","Tnfsf8","Cdk2ap2","Kcnn4","Rpl10a","D198wg1357e","Ivnslabp","Bcas2","Znf512b","Tagap1","Ptma","Nop10","Tnfrsf26","Alyref","Rpl36a","Npm1","Dennd2d","Lsm7","Il16","1500012F01Rik","Trp53","Rpl15-ps3","Clec2i","Rpl10","Samp","Srsf2","Inpp4b","Rps20","Emg1","Phb","Hnmpr","Nap1l1","Cops6","Prdx6","Lsm3","Adk","Mrps6","Cdk4")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=800,height=300)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=8,height=4)
ggplot(cor1,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=8,height=4)
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
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
ggplot(data,aes(x=TNK1@active.ident,y=Sig.score,fill=TNK1@active.ident,color=TNK1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()


############# Fig2 L
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

#T细胞抗原提呈和加工能力
mmu04612<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04612",]
SYMBOL<- bitr(mmu04612$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-TNK1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
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
write.table(newdata,"newdata-APC.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-APC2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-T.APC.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-T.APC.jpeg",width=900,height=300)
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
postscript("time-zscore-T.APC.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#T细胞增殖活化
mmu04660<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04660",]
SYMBOL<- bitr(mmu04660$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-TNK1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
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
pdf("time-zscore-T.proliferation.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-T.proliferation.jpeg",width=900,height=300)
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
postscript("time-zscore-T.proliferation.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#T细胞效应因子得分
#mmu04660<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04660",]
#SYMBOL<- bitr(mmu04660$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-TNK1
g<-c("Gzmb","Tnf","Bnlhe40","Cd101","Il2","Prf1","Elf4","Ifng","Tbx21","Hck","Cd274 ","Ifngr1")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
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
write.table(newdata,"newdata-effector_factor.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-effector_factor2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-T.effector_factor.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-T.effector_factor.jpeg",width=900,height=300)
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
postscript("time-zscore-T.effector_factor.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#T细胞竭耗
#mmu04660<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04660",]
#SYMBOL<- bitr(mmu04660$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-TNK1
g<-c("Pdcd1","Tox","Lag3","Tigit","Pdcd1","Cd274","Pdcd1lg2","Slc5a11","Havcr2","Lgals9","Ctla4")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
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
write.table(newdata,"newdata-exhaustion.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-exhaustion2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-T.exhaustion.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-T.exhaustion.jpeg",width=900,height=300)
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
postscript("time-zscore-T.exhaustion.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

#细胞凋亡
mmu04660<-PATH_ID_NAME[PATH_ID_NAME$KEGGID=="mmu04660",]
SYMBOL<- bitr(mmu04660$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)

k<-TNK1
g<-SYMBOL$SYMBOL
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(TNK1@meta.data, TNK1@reductions$tsne@cell.embeddings,Type=TNK1@active.ident,k)
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
write.table(newdata,"newdata-apoptotic.txt",sep="\t",quote=F,row.names=F)

per<-read.table("newdata-apoptotic2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-T.apoptotic.pdf",width=10,height=3)
ggplot(data = per, aes(x = group, y =Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
#  scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("time-zscore-T.apoptotic.jpeg",width=900,height=300)
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
postscript("time-zscore-T.apoptotic.eps",width=10,height=3)
ggplot(data = per, aes(x = group, y = Sig.score, group= batch, colour=batch)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Time", y = "Sig.score") +
 # scale_fill_manual(values = my36colors, 0.5) +
  theme_bw() + 
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

################################# new ############
TNK<-readRDS("TNK.rds")

TNK <- FindClusters(TNK, resolution = 1)
TNK <- RunTSNE(TNK,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(TNK, reduction = "tsne",label="T",label.size=3.5,pt.size = 1)
TNK <- RunUMAP(TNK,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(TNK, reduction = "umap",label="T",label.size=3.5,pt.size = 1)

T.markers <- FindAllMarkers(TNK, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- T.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(T, features = top10$gene) + NoLegend()
dev.off()

saveRDS(TNK,"TNK.rds")

