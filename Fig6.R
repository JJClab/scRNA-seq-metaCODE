library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
library(ggplot2)

scRNA<-readRDS("scRNA.rds")
new.cluster.ids <- c("MSC", "SMC","SMC", "Macrophage", "Macrophage", "MSC", 
                       "NK", "EC", "Monocyte", "T cell", "B cell",
                       "Macrophage", "DC", "T cell", "Macrophage", "Macrophage",
                       "MSC", "MSC", "T cell", "B cell", "Macrophage",
                       "pDC", "Macrophage", "Neuron", "Macrophage", "SMC", "B cell","SMC")

names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)


EC<-subset(scRNA,idents="EC")
saveRDS(EC,"EC.rds")

EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000) 
al.genes<-row.names(EC) 
EC <- ScaleData(EC, features = al.genes) 
EC <- RunPCA(EC, features = VariableFeatures(object = EC), npcs = 50)   
EC <- FindNeighbors(EC, reduction = "pca", dims = 1:30) 
EC <- FindClusters(EC, resolution = 0.8)
EC <- RunTSNE(EC,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(EC, reduction = "tsne",label="T",label.size=3.5,pt.size = 1)
EC <- RunUMAP(EC,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(EC, reduction = "umap",label="T",label.size=3.5,pt.size = 1)

#8,9,12,14,17,20,22,
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")
names(new.cluster.ids) <- levels(EC)
EC <- RenameIdents(EC, new.cluster.ids)

EC.markers <- FindAllMarkers(EC, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- EC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.cell.pdf",width=25,height=30)
DoHeatmap(EC, features = c('Glycam1', 'Fut7', 'Gcnt1', 'Chst4', 'B3gnt3', 'Ccl21a')) + NoLegend()
dev.off()

#################### 挑选样品分析比例
Idents(EC) <- "stim"
EC1<-subset(EC,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE2","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7","VG","WT1","WT2","IH2","IH4"))

Idents(EC) <- "seurat_clusters"
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")

names(new.cluster.ids) <- levels(EC)
EC <- RenameIdents(EC, new.cluster.ids)
#M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(EC1) <- "seurat_clusters"
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")

names(new.cluster.ids) <- levels(EC1)
EC1 <- RenameIdents(EC1, new.cluster.ids)
#B1@active.ident=factor(B1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

pdf("DimPlot.tsne.pdf",width=8,height=5)
DimPlot(EC1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.tsne.jpg",width=800,height=500)
DimPlot(EC1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.tsne.eps",width=8,height=5)
DimPlot(EC1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

pdf("DimPlot.umap.pdf",width=10,height=5)
DimPlot(EC1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.umap.jpg",width=1000,height=500)
DimPlot(EC1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.umap.eps",width=10,height=5)
DimPlot(EC1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

EC1$batch<-factor(EC1$batch,levels = c("WT","AS","AAA","IH","IG","AG"))
pdf("DimPlot.umap-batch.pdf",width=18,height=3)
DimPlot(EC1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.umap-batch.jpg",width=1800,height=300)
DimPlot(EC1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.umap-batch.eps",width=18,height=3)
DimPlot(EC1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

pdf("DimPlot.tsne-batch.pdf",width=18,height=3)
DimPlot(EC1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.tsne-batch.jpg",width=1800,height=300)
DimPlot(EC1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.tsne-batch.eps",width=18,height=3)
DimPlot(EC1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

### percent
df1<-table(EC1@active.ident,EC1$batch)
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

### percent early
Idents(EC) <- "stim"
early<-subset(EC,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","IH2"))
Idents(EC) <- "seurat_clusters"
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")
names(new.cluster.ids) <- levels(EC)
EC <- RenameIdents(EC, new.cluster.ids)
#EC@active.ident=factor(EC@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)
#early@active.ident=factor(early@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


df1<-table(early@active.ident,early$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#df1$Type<-factor(df1$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("meta-percent-early.barplot.pdf",width=5,height=6)
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

jpeg("metafig-percent-early.barplot.jpg")
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

setEPS()
postscript("metafig-percent-early.barplot.eps",width=5,height=6)
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

### percent late
Idents(EC) <- "stim"
late<-subset(EC,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","IH4"))
Idents(EC) <- "seurat_clusters"
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")
names(new.cluster.ids) <- levels(EC)
EC <- RenameIdents(EC, new.cluster.ids)
#EC@active.ident=factor(EC@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Vascular EC","Lymphatic EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Vascular EC","Lymphatic EC","Vascular EC","Lymphatic EC")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)
#late@active.ident=factor(late@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


df1<-table(late@active.ident,late$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
df1$Type<-factor(df1$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("meta-percent-late.barplot.pdf",width=5,height=6)
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

jpeg("metafig-percent-late.barplot.jpg")
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()

setEPS()
postscript("metafig-percent-late.barplot.eps",width=5,height=6)
ggplot(data = df1,aes(x = Batch, y = Percent,fill=Type)) + 
  geom_bar(stat = 'identity') + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=12,face = "bold"))
dev.off()


#Fig2 E
df<-table(EC1@active.ident,EC1$stim)

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


########### zscore 与招募有关的细胞趋化因子
k<-EC
g<-c('Glycam1', 'Fut7', 'Gcnt1', 'Chst4', 'B3gnt3', 'Ccl21a')
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(EC@meta.data, EC@reductions$tsne@cell.embeddings,Type=EC@active.ident,k)
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

### early
data.bak<-data
earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
earlyindex<-vector()
for (i in 1:length(earlyn)){
	aa<-which(data$stim==earlyn[i])
	earlyindex<-c(earlyindex,aa)
}
data<-data[earlyindex,]
data$batch<-factor(data$batch,levels=c("WT","AS","AAA","IH","IG","AG"))

jpeg("sig.score-Recruit-early.jpg",width=1000,height=300)
p<-ggplot(data,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	facet_grid(.~batch, scale='free')
	print(p)
dev.off()
pdf("sig.score-Recruit-early.pdf",width=10,height=3)
p<-ggplot(data,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(.~batch, scale='free')
	print(p)
dev.off()
setEPS()
postscript("sig.score-Recruit-early.eps",width=10,height=3)
p<-ggplot(data,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))+facet_grid(.~batch, scale='free')
	print(p)
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

jpeg("sig.score-Recruit-late.jpg",width=1000,height=300)
p<-ggplot(data,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	facet_grid(.~batch, scale='free')
	print(p)
dev.off()
pdf("sig.score-Recruit-late.pdf",width=10,height=3)
p<-ggplot(data,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(.~batch, scale='free')
	print(p)
dev.off()
setEPS()
postscript("sig.score-Recruit-late.eps",width=10,height=3)
p<-ggplot(data,aes(x=tSNE_1,y=tSNE_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))+facet_grid(.~batch, scale='free')
	print(p)
dev.off()

########### zscore 与三级淋巴组织形成有关的细胞趋化因子
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)

k<-EC1
g<-c("Ccl2","Ccl3","Ccl4","Ccl5","Ccl8","Ccl18","Ccl19","Ccl21",
"Cxcl9","Cxcl10","Cxcl11","Cxcl13","Ccr7","Cxcr5","Sell","Lamp3")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(EC1@meta.data, EC1@reductions$tsne@cell.embeddings,Type=EC1@active.ident,k)
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
ggplot(data,aes(x=EC1@active.ident,y=Sig.score,fill=EC1@active.ident,color=EC1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score.pdf",width=10,height=10)
ggplot(data,aes(x=EC1@active.ident,y=Sig.score,fill=EC1@active.ident,color=EC1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
ggplot(data,aes(x=EC1@active.ident,y=Sig.score,fill=EC1@active.ident,color=EC1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### early
data.bak<-data
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
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Lymphatic EC","Vascular EC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

################ Fig2 H
data<-data.bak
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
#perl sort.sig.time.pl newdata.txt > newdata2.txt

per<-read.table("newdata2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
#per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

pdf("time-zscore-TLS.pdf",width=5,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-TLS.jpeg",width=500,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-TLS.eps",width=5,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()


