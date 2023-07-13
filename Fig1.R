library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
library(ggplot2)

####################################################################################
scRNA<-readRDS("scRNA.rds")
Macrophage<-subset(scRNA,idents="Macrophage")
saveRDS(Macrophage,"Macrophage.rds")
NK<-subset(scRNA,idents="NK")
saveRDS(NK,"NK.rds")
Monocyte<-subset(scRNA,idents="Monocyte")
saveRDS(Monocyte,"Monocyte.rds")
T<-subset(scRNA,idents="T cell")
saveRDS(T,"T.rds")
B<-subset(scRNA,idents="B cell")
saveRDS(B,"B.rds")
DC<-subset(scRNA,idents="DC")
saveRDS(DC,"DC.rds")
pDC<-subset(scRNA,idents="pDC")
saveRDS(pDC,"pDC.rds")
###########################################
combined<-subset(scRNA,idents=c("Macrophage","NK","Monocyte","T cell","B cell","DC","pDC"))
pc.num=1:30
aa <- RunTSNE(combined, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.6)
saveRDS(aa,"fig1.rds")

combined<-readRDS("fig1.rds")
scRNA.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(scRNA.markers,"combined.markers.rds")
#write.table (scRNA.markers, file ="tsne0.8.csv", sep =",", row.names = F, col.names =TRUE, quote =TRUE)
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(combined, features = top10$gene) + NoLegend()
dev.off()

AverageExp<-AverageExpression(combined,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pdf("corr.pdf",width=20,height=20)
pheatmap(coorda$r)
dev.off()

combined.bak<-combined

saveRDS(combined.bak,"fig1.rds")
combined<-readRDS("fig1.rds")
combined<-subset(combined,idents = 24,invert=T)
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)

combined@active.ident=factor(combined@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))

################################
Bgroup = data.frame(group =B@active.ident,cellID = names(B@active.ident))
write.table(Bgroup,"Bgroup.xls",sep="\t",quote=F,row.names=F)
cDCgroup = data.frame(group =cDC@active.ident,cellID = names(cDC@active.ident))
write.table(cDCgroup,"cDCgroup.xls",sep="\t",quote=F,row.names=F)
Macgroup = data.frame(group =Macrophage@active.ident,cellID = names(Macrophage@active.ident))
write.table(Macgroup,"Macgroup.xls",sep="\t",quote=F,row.names=F)
Mongroup = data.frame(group =Monocyte@active.ident,cellID = names(Monocyte@active.ident))
write.table(Mongroup,"Mongroup.xls",sep="\t",quote=F,row.names=F)
NKgroup = data.frame(group =NK@active.ident,cellID = names(NK@active.ident))
write.table(NKgroup,"NKgroup.xls",sep="\t",quote=F,row.names=F)
pDCgroup = data.frame(group =pDC@active.ident,cellID = names(pDC@active.ident))
write.table(pDCgroup,"pDCgroup.xls",sep="\t",quote=F,row.names=F)
Tgroup = data.frame(group =T@active.ident,cellID = names(T@active.ident))
write.table(Tgroup,"Tgroup.xls",sep="\t",quote=F,row.names=F)
write.table(combined.bak@meta.data,"allT.metadata.xls",sep="\t",quote=F)

#cat Bgroup.xls cDCgroup.xls Macgroup.xls Mongroup.xls NKgroup.xls  Tgroup.xls pDCgroup.xls |perl -ne 'chomp;@a=split/\t/;print "$a[1]\t$a[0]\n"' > newTgroup.xls

newg<-read.table("newTgroup.xls",header=TRUE,sep="\t")
combined$group<-newg$group
Idents(combined)<-"group"
################################

### 查看默认配色 ### 
show_col(hue_pal()(11))
hue_pal()(11)
cells<-c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC")
colors<-c("#AEA200","#64B200","#00C1A7","#00BADE","#00A6FF","#B385FF","#EF67EB")

#################### 挑选样品分析比例
Idents(combined) <- "stim"
combined1<-subset(combined,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE2","ApoE3","ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7","Ldlr8","VG","WT1","WT2","WT3","IH2","IH4"))

Idents(combined) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined@active.ident=factor(combined@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(combined1) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined1)
combined1 <- RenameIdents(combined1, new.cluster.ids)
combined1@active.ident=factor(combined1@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))

pdf("DimPlot.tsne.pdf",width=8,height=5)
DimPlot(combined1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0,cols=colors)
dev.off()
jpeg("DimPlot.tsne.jpg",width=800,height=500)
DimPlot(combined1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0,cols=colors)
dev.off()
setEPS()
postscript("DimPlot.tsne.eps",width=8,height=5)
DimPlot(combined1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0,cols=colors)
dev.off()

pdf("DimPlot.umap.pdf",width=8,height=5)
DimPlot(combined1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0,cols=colors)
dev.off()
jpeg("DimPlot.umap.jpg",width=800,height=500)
DimPlot(combined1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0,cols=colors)
dev.off()
setEPS()
postscript("DimPlot.umap.eps",width=8,height=5)
DimPlot(combined1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0,cols=colors)
dev.off()

combined1$batch<-factor(combined1$batch,levels = c("WT","AS","AAA","IH","IG","AG"))
pdf("DimPlot.umap-batch.pdf",width=18,height=3)
DimPlot(combined1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch",cols=colors)
dev.off()
jpeg("DimPlot.umap-batch.jpg",width=1800,height=300)
DimPlot(combined1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch",cols=colors)
dev.off()
setEPS()
postscript("DimPlot.umap-batch.eps",width=18,height=3)
DimPlot(combined1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch",cols=colors)
dev.off()

pdf("DimPlot.tsne-batch.pdf",width=18,height=3)
DimPlot(combined1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch",cols=colors)
dev.off()
jpeg("DimPlot.tsne-batch.jpg",width=1800,height=300)
DimPlot(combined1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch",cols=colors)
dev.off()
setEPS()
postscript("DimPlot.tsne-batch.eps",width=18,height=3)
DimPlot(combined1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch",cols=colors)
dev.off()


### percent
df1<-table(combined1@active.ident,combined1$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
df1$Type<-factor(df1$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("meta-percent.pdf",width=13,height=4)
#jpeg("metafig1.jpg",width=800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() +  theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("metafig-percent.jpg",width=1800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
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
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

### percent early
Idents(combined) <- "stim"
early<-subset(combined,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2"))
Idents(combined) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined@active.ident=factor(combined@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)
early@active.ident=factor(early@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


df1<-table(early@active.ident,early$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
df1$Type<-factor(df1$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("meta-percent-early.pdf",width=13,height=4)
#jpeg("metafig1.jpg",width=800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() +  theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("metafig-percent-early.jpg",width=1800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

setEPS()
postscript("metafig-percent-early.eps",width=13,height=4)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

### percent late
Idents(combined) <- "stim"
late<-subset(combined,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4"))
Idents(combined) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined@active.ident=factor(combined@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)
late@active.ident=factor(late@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


df1<-table(late@active.ident,late$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
df1$Type<-factor(df1$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("meta-percent-late.pdf",width=13,height=4)
#jpeg("metafig1.jpg",width=800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() +  theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

jpeg("metafig-percent-late.jpg",width=1800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()

setEPS()
postscript("metafig-percent-late.eps",width=13,height=4)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
  scale_colour_manual(breaks = df1$Type, values = colors) +
  theme_bw() + theme(legend.position="none",strip.text.x = element_text(size = 14)) +
  facet_grid(.~Type, scale='free')
# facet_grid(Type~., scale='free')
dev.off()


###### feature plot
pdf("FeaturePlot1.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Ccl2","Ctsc"),reduction="umap")
dev.off()
jpeg("FeaturePlot1.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Ccl2","Ctsc"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot1.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Ccl2","Ctsc"),reduction="umap")
dev.off()
pdf("FeaturePlot1.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Ccl2","Ctsc"),reduction="tsne")
dev.off()
jpeg("FeaturePlot1.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Ccl2","Ctsc"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot1.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Ccl2","Ctsc"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot2.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Wfdc21","G0s2"),reduction="umap")
dev.off()
jpeg("FeaturePlot2.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Wfdc21","G0s2"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot2.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Wfdc21","G0s2"),reduction="umap")
dev.off()
pdf("FeaturePlot2.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Wfdc21","G0s2"),reduction="tsne")
dev.off()
jpeg("FeaturePlot2.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Wfdc21","G0s2"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot2.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Wfdc21","G0s2"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot3.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Ms4a4b","Cd8b1"),reduction="umap")
dev.off()
jpeg("FeaturePlot3.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Ms4a4b","Cd8b1"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot3.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Ms4a4b","Cd8b1"),reduction="umap")
dev.off()
pdf("FeaturePlot3.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Ms4a4b","Cd8b1"),reduction="tsne")
dev.off()
jpeg("FeaturePlot3.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Ms4a4b","Cd8b1"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot3.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Ms4a4b","Cd8b1"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot4.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Cd79a","Cd79b"),reduction="umap")
dev.off()
jpeg("FeaturePlot4.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Cd79a","Cd79b"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot4.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Cd79a","Cd79b"),reduction="umap")
dev.off()
pdf("FeaturePlot4.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Cd79a","Cd79b"),reduction="tsne")
dev.off()
jpeg("FeaturePlot4.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Cd79a","Cd79b"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot4.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Cd79a","Cd79b"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot5.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("S100a9","S100a8"),reduction="umap")
dev.off()
jpeg("FeaturePlot5.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("S100a9","S100a8"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot5.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("S100a9","S100a8"),reduction="umap")
dev.off()
pdf("FeaturePlot5.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("S100a9","S100a8"),reduction="tsne")
dev.off()
jpeg("FeaturePlot5.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("S100a9","S100a8"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot5.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("S100a9","S100a8"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot6.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Naaa","Wdfy4"),reduction="umap")
dev.off()
jpeg("FeaturePlot6.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Naaa","Wdfy4"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot6.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Naaa","Wdfy4"),reduction="umap")
dev.off()
pdf("FeaturePlot6.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Naaa","Wdfy4"),reduction="tsne")
dev.off()
jpeg("FeaturePlot6.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Naaa","Wdfy4"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot6.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Naaa","Wdfy4"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot7.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Nkg7","Klre1"),reduction="umap")
dev.off()
jpeg("FeaturePlot7.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Nkg7","Klre1"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot7.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Nkg7","Klre1"),reduction="umap")
dev.off()
pdf("FeaturePlot7.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Nkg7","Klre1"),reduction="tsne")
dev.off()
jpeg("FeaturePlot7.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Nkg7","Klre1"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot7.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Nkg7","Klre1"),reduction="tsne")
dev.off()

##############
pdf("FeaturePlot8.umap.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Klk1","Siglech"),reduction="umap")
dev.off()
jpeg("FeaturePlot8.umap.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Klk1","Siglech"),reduction="umap")
dev.off()
setEPS()
postscript("FeaturePlot8.umap.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Klk1","Siglech"),reduction="umap")
dev.off()
pdf("FeaturePlot8.tsne.pdf",width=8,height=4)
FeaturePlot(combined1, features = c("Klk1","Siglech"),reduction="tsne")
dev.off()
jpeg("FeaturePlot8.tsne.jpg",width=800,height=400)
FeaturePlot(combined1, features = c("Klk1","Siglech"),reduction="tsne")
dev.off()
setEPS()
postscript("FeaturePlot8.tsne.eps",width=8,height=4)
FeaturePlot(combined1, features = c("Klk1","Siglech"),reduction="tsne")
dev.off()

###############饼图 WT
WT<-df1[df1$Batch=="WT",]
WT$Type<-factor(WT$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
label <- paste(WT$Type, paste('(', round(WT$Percent * 100, 1), '%)', sep = ''), sep = '')
label<-label[c(1,5,2,3,4,6,7)]
pdf("WT.pie.pdf",width=6,height=6)
p <- ggplot(data = WT, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
jpeg("WT.pie.jpg",width=500,height=500)
p <- ggplot(data = WT, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
setEPS()
postscript("WT.pie.eps",width=6,height=6)
p <- ggplot(data = WT, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()

###############饼图 AS
AS<-df1[df1$Batch=="AS",]
AS$Type<-factor(AS$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
label <- paste(AS$Type, paste('(', round(AS$Percent * 100, 1), '%)', sep = ''), sep = '')
label<-label[c(1,5,2,3,4,6,7)]
pdf("AS.pie.pdf",width=6,height=6)
p <- ggplot(data = AS, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
jpeg("AS.pie.jpg",width=500,height=500)
p <- ggplot(data = AS, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
setEPS()
postscript("AS.pie.eps",width=6,height=6)
p <- ggplot(data = AS, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()

###############饼图 AAA
AAA<-df1[df1$Batch=="AAA",]
AAA$Type<-factor(AAA$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
label <- paste(AAA$Type, paste('(', round(AAA$Percent * 100, 1), '%)', sep = ''), sep = '')
label<-label[c(1,5,2,3,4,6,7)]
pdf("AAA.pie.pdf",width=6,height=6)
p <- ggplot(data = AAA, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
jpeg("AAA.pie.jpg",width=500,height=500)
p <- ggplot(data = AAA, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
setEPS()
postscript("AAA.pie.eps",width=6,height=6)
p <- ggplot(data = AAA, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()

###############饼图 IH
IH<-df1[df1$Batch=="IH",]
IH$Type<-factor(IH$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
label <- paste(IH$Type, paste('(', round(IH$Percent * 100, 1), '%)', sep = ''), sep = '')
label<-label[c(1,5,2,3,4,6,7)]
pdf("IH.pie.pdf",width=6,height=6)
p <- ggplot(data = IH, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
jpeg("IH.pie.jpg",width=500,height=500)
p <- ggplot(data = IH, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
setEPS()
postscript("IH.pie.eps",width=6,height=6)
p <- ggplot(data = IH, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()

###############饼图 IG
IG<-df1[df1$Batch=="IG",]
IG$Type<-factor(IG$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
label <- paste(IG$Type, paste('(', round(IG$Percent * 100, 1), '%)', sep = ''), sep = '')
label<-label[c(1,5,2,3,4,6,7)]
pdf("IG.pie.pdf",width=6,height=6)
p <- ggplot(data = IG, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
jpeg("IG.pie.jpg",width=500,height=500)
p <- ggplot(data = IG, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
setEPS()
postscript("IG.pie.eps",width=6,height=6)
p <- ggplot(data = IG, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()

###############饼图 AG
AG<-df1[df1$Batch=="AG",]
AG$Type<-factor(AG$Type,levels=c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))
label <- paste(AG$Type, paste('(', round(AG$Percent * 100, 1), '%)', sep = ''), sep = '')
label<-label[c(1,5,2,3,4,6,7)]
pdf("AG.pie.pdf",width=6,height=6)
p <- ggplot(data = AG, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
jpeg("AG.pie.jpg",width=500,height=500)
p <- ggplot(data = AG, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()
setEPS()
postscript("AG.pie.eps",width=6,height=6)
p <- ggplot(data = AG, mapping = aes(x = 'Percent', y = Percent, fill = Type)) + geom_bar(stat = 'identity', position = 'stack', colour="black") + scale_fill_manual(values = colors,labels=label)
p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(panel.background = element_blank(),panel.grid=element_blank(),axis.text = element_blank(),axis.ticks = element_blank()) 
dev.off()

########### zscore 控制淋巴细胞转运的基因
#select early
Idents(combined) <- "stim"
early<-subset(combined,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2"))

Idents(combined) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined@active.ident=factor(combined@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)
early@active.ident=factor(early@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))

##### zscore
k<-early
g<-c('Glycam1', 'Fut7', 'Gcnt1', 'Chst4', 'B3gnt3', 'Ccl21a')
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(early@meta.data, early@reductions$tsne@cell.embeddings,k)
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
postscript("sig.score2.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score.jpg",width=800,height=700)
ggplot(data,aes(x=early@active.ident,y=Sig.score,fill=early@active.ident,color=early@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Macrophage","NK"),c("Macrophage","Monocyte"),c("Macrophage","T"),c("Macrophage","cell"),c("Macrophage","B"),c("Macrophage","cell"),c("Macrophage","cDC"),c("Macrophage","pDC"),c("NK","Monocyte"),c("NK","T"),c("NK","cell"),c("NK","B"),c("NK","cell"),c("NK","cDC"),c("NK","pDC"),c("Monocyte","T"),c("Monocyte","cell"),c("Monocyte","B"),c("Monocyte","cell"),c("Monocyte","cDC"),c("Monocyte","pDC"),c("T","cell"),c("T","B"),c("T","cell"),c("T","cDC"),c("T","pDC"),c("cell","B"),c("cell","cell"),c("cell","cDC"),c("cell","pDC"),c("B","cell"),c("B","cDC"),c("B","pDC"),c("cell","cDC"),c("cell","pDC"),c("cDC","pDC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score.pdf",width=10,height=10)
ggplot(data,aes(x=early@active.ident,y=Sig.score,fill=early@active.ident,color=early@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Macrophage","NK"),c("Macrophage","Monocyte"),c("Macrophage","T"),c("Macrophage","cell"),c("Macrophage","B"),c("Macrophage","cell"),c("Macrophage","cDC"),c("Macrophage","pDC"),c("NK","Monocyte"),c("NK","T"),c("NK","cell"),c("NK","B"),c("NK","cell"),c("NK","cDC"),c("NK","pDC"),c("Monocyte","T"),c("Monocyte","cell"),c("Monocyte","B"),c("Monocyte","cell"),c("Monocyte","cDC"),c("Monocyte","pDC"),c("T","cell"),c("T","B"),c("T","cell"),c("T","cDC"),c("T","pDC"),c("cell","B"),c("cell","cell"),c("cell","cDC"),c("cell","pDC"),c("B","cell"),c("B","cDC"),c("B","pDC"),c("cell","cDC"),c("cell","pDC"),c("cDC","pDC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
ggplot(data,aes(x=early@active.ident,y=Sig.score,fill=early@active.ident,color=early@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Macrophage","NK"),c("Macrophage","Monocyte"),c("Macrophage","T"),c("Macrophage","cell"),c("Macrophage","B"),c("Macrophage","cell"),c("Macrophage","cDC"),c("Macrophage","pDC"),c("NK","Monocyte"),c("NK","T"),c("NK","cell"),c("NK","B"),c("NK","cell"),c("NK","cDC"),c("NK","pDC"),c("Monocyte","T"),c("Monocyte","cell"),c("Monocyte","B"),c("Monocyte","cell"),c("Monocyte","cDC"),c("Monocyte","pDC"),c("T","cell"),c("T","B"),c("T","cell"),c("T","cDC"),c("T","pDC"),c("cell","B"),c("cell","cell"),c("cell","cDC"),c("cell","pDC"),c("B","cell"),c("B","cDC"),c("B","pDC"),c("cell","cDC"),c("cell","pDC"),c("cDC","pDC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

cs<-levels(combined)
cellType<-data.frame(Type=early@active.ident)
data1<-cbind(data,cellType)

for (i in 1:length(cs)){
	plotscore<-data1[data1$Type==cs[i],]
	out<-paste(cs[i],"score2","jpg",sep=".")
	jpeg(out,width=800,height=700)
	p<-ggplot(plotscore,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(title=cs[i], x="",y="Zscore",size=4)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=18),plot.title = element_text(hjust = 0.5)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
	print(p)
	dev.off()
	
	out<-paste(cs[i],"score2","pdf",sep=".")
	pdf(out,width=10,height=8)
	p<-ggplot(plotscore,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(title=cs[i], x="",y="Zscore",size=4)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=18),plot.title = element_text(hjust = 0.5)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
	print(p)
	dev.off()
	
	out<-paste(cs[i],"score2","eps",sep=".")
	setEPS()
	postscript(out,width=10,height=8)
	p<-ggplot(plotscore,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(title=cs[i], x="",y="Zscore",size=4)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=18),plot.title = element_text(hjust = 0.5)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
	print(p)
	dev.off()
}



##########################select late
Idents(combined) <- "stim"
late<-subset(combined,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4"))

Idents(combined) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined@active.ident=factor(combined@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))


Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("Macrophage", "Macrophage","Monocyte", "T cell", "B cell", "NK", 
                       "Macrophage", "Macrophage", "NK", "T cell", "Macrophage",
                       "Macrophage", "cDC", "T cell", "Macrophage", "B cell",
                       "Macrophage", "cDC", "pDC", "Macrophage", "Macrophage",
                       "Monocyte", "T cell", "B cell")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)
late@active.ident=factor(late@active.ident,levels = c("Macrophage","NK","Monocyte","T cell","B cell","cDC","pDC"))

##### zscore
k<-late
g<-c('Glycam1', 'Fut7', 'Gcnt1', 'Chst4', 'B3gnt3', 'Ccl21a')
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(late@meta.data, late@reductions$tsne@cell.embeddings,k)
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
postscript("sig.score2.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score.jpg",width=800,height=700)
ggplot(data,aes(x=late@active.ident,y=Sig.score,fill=late@active.ident,color=late@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Macrophage","NK"),c("Macrophage","Monocyte"),c("Macrophage","T"),c("Macrophage","cell"),c("Macrophage","B"),c("Macrophage","cell"),c("Macrophage","cDC"),c("Macrophage","pDC"),c("NK","Monocyte"),c("NK","T"),c("NK","cell"),c("NK","B"),c("NK","cell"),c("NK","cDC"),c("NK","pDC"),c("Monocyte","T"),c("Monocyte","cell"),c("Monocyte","B"),c("Monocyte","cell"),c("Monocyte","cDC"),c("Monocyte","pDC"),c("T","cell"),c("T","B"),c("T","cell"),c("T","cDC"),c("T","pDC"),c("cell","B"),c("cell","cell"),c("cell","cDC"),c("cell","pDC"),c("B","cell"),c("B","cDC"),c("B","pDC"),c("cell","cDC"),c("cell","pDC"),c("cDC","pDC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score.pdf",width=10,height=10)
ggplot(data,aes(x=late@active.ident,y=Sig.score,fill=late@active.ident,color=late@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Macrophage","NK"),c("Macrophage","Monocyte"),c("Macrophage","T"),c("Macrophage","cell"),c("Macrophage","B"),c("Macrophage","cell"),c("Macrophage","cDC"),c("Macrophage","pDC"),c("NK","Monocyte"),c("NK","T"),c("NK","cell"),c("NK","B"),c("NK","cell"),c("NK","cDC"),c("NK","pDC"),c("Monocyte","T"),c("Monocyte","cell"),c("Monocyte","B"),c("Monocyte","cell"),c("Monocyte","cDC"),c("Monocyte","pDC"),c("T","cell"),c("T","B"),c("T","cell"),c("T","cDC"),c("T","pDC"),c("cell","B"),c("cell","cell"),c("cell","cDC"),c("cell","pDC"),c("B","cell"),c("B","cDC"),c("B","pDC"),c("cell","cDC"),c("cell","pDC"),c("cDC","pDC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
ggplot(data,aes(x=late@active.ident,y=Sig.score,fill=late@active.ident,color=late@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("Macrophage","NK"),c("Macrophage","Monocyte"),c("Macrophage","T"),c("Macrophage","cell"),c("Macrophage","B"),c("Macrophage","cell"),c("Macrophage","cDC"),c("Macrophage","pDC"),c("NK","Monocyte"),c("NK","T"),c("NK","cell"),c("NK","B"),c("NK","cell"),c("NK","cDC"),c("NK","pDC"),c("Monocyte","T"),c("Monocyte","cell"),c("Monocyte","B"),c("Monocyte","cell"),c("Monocyte","cDC"),c("Monocyte","pDC"),c("T","cell"),c("T","B"),c("T","cell"),c("T","cDC"),c("T","pDC"),c("cell","B"),c("cell","cell"),c("cell","cDC"),c("cell","pDC"),c("B","cell"),c("B","cDC"),c("B","pDC"),c("cell","cDC"),c("cell","pDC"),c("cDC","pDC")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

cs<-levels(combined)
cellType<-data.frame(Type=late@active.ident)
data1<-cbind(data,cellType)

for (i in 1:length(cs)){
	plotscore<-data1[data1$Type==cs[i],]
	out<-paste(cs[i],"score2","jpg",sep=".")
	jpeg(out,width=800,height=700)
	p<-ggplot(plotscore,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(title=cs[i], x="",y="Zscore",size=4)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=18),plot.title = element_text(hjust = 0.5)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
	print(p)
	dev.off()
	
	out<-paste(cs[i],"score2","pdf",sep=".")
	pdf(out,width=10,height=8)
	p<-ggplot(plotscore,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(title=cs[i], x="",y="Zscore",size=4)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=18),plot.title = element_text(hjust = 0.5)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
	print(p)
	dev.off()
	
	out<-paste(cs[i],"score2","eps",sep=".")
	setEPS()
	postscript(out,width=10,height=8)
	p<-ggplot(plotscore,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(title=cs[i], x="",y="Zscore",size=4)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=18),plot.title = element_text(hjust = 0.5)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
	print(p)
	dev.off()
}


