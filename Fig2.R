library(Seurat)
library(ggplot2)
library(data.table)
library(hdf5r)
library(dplyr)
library("ggthemes")
library("RColorBrewer")

setwd("C:\\addIH\\Fig1")
#Mac<-readRDS("Macrophage.rds")
#Mon<-readRDS("Monocyte.rds")
#MacMon<-merge(Mac,Mon)

MacMon<-subset(combined,idents=c("Macrophage","Monocyte"))

MacMon <- FindVariableFeatures(MacMon, selection.method = "vst", nfeatures = 2000) 
al.genes<-row.names(MacMon) 
MacMon <- ScaleData(MacMon, features = al.genes) 
MacMon <- RunPCA(MacMon, features = VariableFeatures(object = MacMon), npcs = 50)   
MacMon <- FindNeighbors(MacMon, reduction = "pca", dims = 1:30) 
MacMon <- FindClusters(MacMon, resolution = 0.5)
MacMon <- RunTSNE(MacMon,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(MacMon, reduction = "tsne",label="T",label.size=3.5,pt.size = 1)
MacMon <- RunUMAP(MacMon,Reduction = "pca",dims = 1:20, check_duplicates = FALSE)
DimPlot(MacMon, reduction = "umap",label="T",label.size=3.5,pt.size = 1)

M.markers <- FindAllMarkers(MacMon, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top10 <- M.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DoHeatmap.pdf",width=25,height=30)
DoHeatmap(M, features = top10$gene) + NoLegend()
dev.off()

saveRDS(MacMon,"MacMon.rds")

setwd("C:\\addIH\\Fig2")
M<-readRDS("MacMon.rds")
###############################
#平均表达谱函数AverageExpression
AverageExp<-AverageExpression(M,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pdf("corr.pdf",width=20,height=20)
pheatmap(coorda$r)
dev.off()
##################################
#TREM2 hi Macrophage
StackedVlnPlot(M,  c('Trem2', 'Cd9', 'Lgals3', 'Spp1', 'Aldoa', "Ctsd"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c('Ms4a7', 'C1qa', 'C1qb', 'Cx3cr1', 'Cd81', "Ly86"), pt.size=0, cols=my36colors)
#RELMa+ Resident-like Macrophage
StackedVlnPlot(M,  c('Fcrls', 'Fgfr1', 'Retnla', 'Ear2', 'Clec4b1', "Mgl2"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c('Lyz1', 'Pdlim1', 'Lpl', 'Cd81', "Fn1"), pt.size=0, cols=my36colors)
#FOLR2+ Resident-like Macrophage视为抗炎M2
StackedVlnPlot(M,  c('Ccl8', 'Folr2', 'Cbr2', 'Gas6', 'C4b', "Mrc1"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c('Clec10a', 'Selenop', 'Pltp', 'F13a1', 'Fcgrt', "Pf4"), pt.size=0, cols=my36colors)
#Inflammation Macrophage
StackedVlnPlot(M,  c('Il1b', 'Il1a', 'Cebpb', 'Egr1', 'Ccl4', "Cxcl2"), pt.size=0, cols=my36colors)
#iNOS+ Inflammation Macrophage/NMES1 Macrophage
StackedVlnPlot(M,  c('Cxcl9', 'AA467197', 'Ly6i', 'Ass1', 'Upp1', "Gbp2"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c('Ly6a', 'Slamf8', 'Fam26f', 'Prdx5', 'Cxcl16', "Cxcl10"), pt.size=0, cols=my36colors)
#ISG hi Inflammation Macrophage 
StackedVlnPlot(M,  c('Ccr2', 'Lyz2', 'Ms4a6c', 'Plac8', 'Clec4a3', "Pld4"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c('Ifitm6', 'Cybb', 'Clec4a1', 'Cx3cr1', 'Ifitm3', "Chil3"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c('Ifit3b', 'Ifit3', 'Ifit1', 'Cxcl10', 'Ifi205'), pt.size=0, cols=my36colors)
#CD16 hi Monocyte
StackedVlnPlot(M,  c('Mmp9', 'Cxcr2', 'Il1r2', 'Hdc', 'Msrb1', "Il1b"), pt.size=0, cols=my36colors)
StackedVlnPlot(M,  c( 'Cxcl2', 'Ccrl2', 'Ccl3', 'Cd274', "Acod1"), pt.size=0, cols=my36colors)

new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M)
M <- RenameIdents(M, new.cluster.ids)

#################### 挑选样品分析比例
Idents(M) <- "stim"
M1<-subset(M,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE2","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT1","WT2","WT3","IH2","IH4"))

Idents(M) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M)
M <- RenameIdents(M, new.cluster.ids)
M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage",
"Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

Idents(M1) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)
M1@active.ident=factor(M1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

pdf("DimPlot.tsne.pdf",width=8,height=5)
DimPlot(M1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.tsne.jpg",width=800,height=500)
DimPlot(M1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.tsne.eps",width=8,height=5)
DimPlot(M1, reduction = "tsne",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

pdf("DimPlot.umap.pdf",width=10,height=5)
DimPlot(M1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
jpeg("DimPlot.umap.jpg",width=1000,height=500)
DimPlot(M1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()
setEPS()
postscript("DimPlot.umap.eps",width=10,height=5)
DimPlot(M1, reduction = "umap",label=T,label.size=3.5,pt.size = 2.0)
dev.off()

M1$batch<-factor(M1$batch,levels = c("WT","AS","AAA","IH","IG","AG"))
pdf("DimPlot.umap-batch.pdf",width=18,height=3)
DimPlot(M1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.umap-batch.jpg",width=1800,height=300)
DimPlot(M1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.umap-batch.eps",width=18,height=3)
DimPlot(M1, reduction = "umap",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

pdf("DimPlot.tsne-batch.pdf",width=18,height=3)
DimPlot(M1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
jpeg("DimPlot.tsne-batch.jpg",width=1800,height=300)
DimPlot(M1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()
setEPS()
postscript("DimPlot.tsne-batch.eps",width=18,height=3)
DimPlot(M1, reduction = "tsne",label=T,label.size=3.5,pt.size = 1,split.by="batch")
dev.off()

###### heatmap top200 ###
library(pheatmap)
library(ComplexHeatmap)

M.markers <- FindAllMarkers(M1, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)
top20 <- M.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

cts <- GetAssayData(M1, slot = 'scale.data')
group = data.frame(group =M1@active.ident,row.names = names(M1@active.ident))

#new_cluster <- sort(scRNA1@active.ident)
#ctstop <- as.matrix(cts[top20$gene, colnames(cts)])
ctstop <- as.matrix(cts[top20$gene, colnames(M1)])

#ctstop_bak<-ctstop

#select 200 genes 
types<-c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte")
cells<-vector()
for (i in 1:length(types)){
set.seed(1111)
index<-which(group$group==types[i])
select<-sample(index,140)
ids<-rownames(group)[select]
cells<-c(cells,ids)
}
tmp2<-ctstop[,cells]
gindex<-which(rownames(group) %in% cells)
group2<-data.frame(row.names=rownames(group)[gindex],group=group[gindex,])
tmp2[tmp2>5] <- 5
tmp2[tmp2< -4] <- -4

pdf("heatmap.pdf")
pheatmap(tmp2,cluster_row=F,cluster_col=F,annotation_col=group2, color = colorRampPalette(colors = c("#136BA5", "#ffffff","#BF172A"))(50),show_rownames=T,show_colnames=F,fontsize_row = 4,fontsize_col = 8,border=F)
dev.off()
setEPS()
postscript("heatmap.eps")
pheatmap(tmp2,cluster_row=F,cluster_col=F,annotation_col=group2, color = colorRampPalette(colors = c("#136BA5", "#ffffff","#BF172A"))(50),show_rownames=T,show_colnames=F,fontsize_row = 4,fontsize_col = 8,border=F)
dev.off()

#### G Featureplot top2 ###
types<-c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte")
for (i in 1:length(types)){
	index<-which(top20$cluster==types[i])
	select<-index[1:2]
	features <- top20$gene[select]
	out<-paste("FeaturePlot",types[i],"pdf",sep=".")
	pdf(out,width=10,height=4)
	p<-FeaturePlot(M1, features = features,reduction="tsne")
	print(p)
	dev.off()
	out1<-paste("FeaturePlot",types[i],"eps",sep=".")
	setEPS()
	postscript(out1,width=10,height=4)
	p1<-FeaturePlot(M1, features = features,reduction="tsne")
	print(p1)
	dev.off()
}

### percent
df1<-table(M1@active.ident,M1$batch)
df1<-prop.table(df1,2)
df1<-as.data.frame(df1)
colnames(df1)<-c("Type","Batch","Percent")
df1$Batch<-factor(df1$Batch,levels=c("WT","AS","AAA","IH","IG","AG"))
df1$Type<-factor(df1$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))
mycolors<-c("#6495ED","#FFA500","#228B22","#FF4500")

pdf("meta-percent.pdf",width=13,height=4)
#jpeg("metafig1.jpg",width=800,height=300)
ggplot(data = df1, aes(x = Batch, y = Percent, group= Type, colour=Type)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "Batch", y = "Percent") +
# scale_colour_manual(breaks = df1$Type, values = colors) +
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

##########################################
#Fig2 E
df<-table(M1@active.ident,M1$stim)

df<-df[,c("AAA1","AAA2","AG2","AG4","ApoE4","ApoE5","IG2","IG4","Ldlr2","Ldlr3","IH2","IH4")]

df<-prop.table(df,2)

AS1<-data.frame(apply(df[,c(5,9)],1,mean))
AS2<-data.frame(apply(df[,c(6,10)],1,mean))
colnames(AS1)<-"AS1"
colnames(AS2)<-"AS2"
df<-cbind(df[,1:ncol(df)],AS1=AS1$AS1,AS2=AS2$AS2)
df<-df[,c(-5,-6,-9,-10)]

write.table(df,"percent.txt",sep="\t",quote=F)
#perl deal_plot.pl percent.txt > percent.trans.txts
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

### bar ##
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




#################### GO #####################
#M.markers<-readRDS("M.markers.rds")
logFCfilter=0.25
adjPvalFilter=0.25
M.markers <- FindAllMarkers(object = M,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
							   

group_g <- M.markers[,c(6,7)]
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

########### zscore 与炎症因子相关基因Featureplot和小提琴图
### early
#earlyn<-c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2")
Idents(M) <- "stim"
early<-subset(M,idents=c("AAA1","AG2","ApoE4","IG2","Ldlr2","WT1","WT2","WT3","IH2"))
Idents(M) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M)
M <- RenameIdents(M, new.cluster.ids)
M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(early) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(early)
early <- RenameIdents(early, new.cluster.ids)
early@active.ident=factor(early@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))

### late
#laten<-c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4")
Idents(M) <- "stim"
late<-subset(M,idents=c("AAA2","AG4","ApoE5","IG4","Ldlr3","WT1","WT2","WT3","IH4"))
Idents(M) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M)
M <- RenameIdents(M, new.cluster.ids)
M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(late) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(late)
late <- RenameIdents(late, new.cluster.ids)
late@active.ident=factor(late@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


###### 各个batch的打分 ######
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

#entrezid to gene symbol
inf<- c("mmu04062","mmu04061","mmu04670","mmu04657","mmu04514")
batchs<-c("AAA","AG","AS","IG","IH","WT")
for(i in 1:length(inf)){
	dir.create(inf[i])
	setwd(inf[i])
#	pathway<-PATH_ID_NAME[PATH_ID_NAME$KEGGID==inf[i],]
	pathway<-PATH_ID_NAME[PATH_ID_NAME$KEGGID %in% inf,] ### all pathway gene
	SYMBOL<- bitr(pathway$ENTREZID, fromType = "ENTREZID", toType="SYMBOL",OrgDb = org.Mm.eg.db)
	
	for(j in 1:length(batchs)){
		dir.create(batchs[j])
		setwd(batchs[j])
		bb<-subset(early,batch==batchs[j])
		#bb<-subset(late,batch==batchs[j])
		nn<-table(bb$batch)
		print(nn)
		#k<-bb
		g<-SYMBOL$SYMBOL
		gindex<-which(g %in% rownames(bb@assays$RNA@data))
		genes<-g[gindex]

		s<-t(as.matrix(bb@assays$RNA@data[genes,]))
		k<-rowMeans(s)
		data <- cbind(bb@meta.data, bb@reductions$umap@cell.embeddings,Type=bb@active.ident,k)
		Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))
		
		cor1<-cbind(data,Sig.score)
#cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=500,height=300)
p<-ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
print(p)
dev.off()
pdf("sig.score.pdf",width=5,height=3)
p<-ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
		print(p)
dev.off()
setEPS()
postscript("sig.score.eps",width=5,height=3)
p<-ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()

data<-cbind(data,Sig.score)
#data <- data %>% group_by(batch) %>% mutate(med = median(Sig.score))
#median.batch.sig<-aggregate(Sig.score~batch,data=data,median)

jpeg("cells-score.jpg",width=800,height=700)
p<-ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
print(p)
dev.off()
pdf("cells-score.pdf",width=10,height=10)
p<-ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
print(p)
dev.off()
setEPS()
postscript("cells-score.eps",width=10,height=8)
p<-ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
print(p)
dev.off()

jpeg("cells-score2.jpg",width=800,height=700)
p<-ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
print(p)
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
p<-ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
print(p)
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
p<-ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
print(p)
dev.off()
		setwd("../")
	}
	setwd("../")
}

########## monocle ########


########## iTALK ########
identity = data.frame(cellID = names(M@active.ident),cellType =M@active.ident) 
#identity$type<-type
write.table(identity,"M-cellType.txt",quote=F,sep="\t",row.names=F)

#perl select-T.pl  M-cellType.txt iTALK/all.newtype.txt > all.newtype.txt

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

newtype<-read.table("all.newtype.txt",header = T,sep="\t")
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

iTalk_data <- as.data.frame(t(as.matrix(early@assays$RNA@counts)))
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


iTalk_data <- as.data.frame(t(as.matrix(late@assays$RNA@counts)))
iTalk_data$cell_type <- late$group
iTalk_data$compare_group <- late$batch



########### zscore 与三级淋巴组织形成有关的细胞趋化因子
library(reshape2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
Idents(M) <- "stim"
M1<-subset(M,idents=c("AAA2","AAA1","AG2","AG4","ApoE1","ApoE2","ApoE3",
"ApoE4","ApoE5","ApoE6","IG2","IG4","Ldlr2","Ldlr3","Ldlr4","Ldlr7",
"Ldlr8","VG","WT1","WT2","WT3","IH2","IH4"))

Idents(M) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M)
M <- RenameIdents(M, new.cluster.ids)
M@active.ident=factor(M@active.ident, levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


Idents(M1) <- "seurat_clusters"
new.cluster.ids <- c("Res-like Macrophage", "TREM2 Macrophage","Inflammation Macrophage",
					   "TREM2 Macrophage", "TREM2 Macrophage", "TREM2 Macrophage", 
                       "TREM2 Macrophage", "Inflammation Macrophage", "Res-like Macrophage", 
					   "CD16 hi Monocyte", "TREM2 Macrophage",
                       "CD16 hi Monocyte", "ISG hi Monocyte", "TREM2 Macrophage",
					   "Inflammation Macrophage", "Res-like Macrophage",
                       "Inflammation Macrophage", "TREM2 Macrophage", "ISG hi Monocyte", 
					   "ISG hi Monocyte", "TREM2 Macrophage")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)
M1@active.ident=factor(M1@active.ident,levels = c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))


k<-M1
g<-c("Ccl2","Ccl3","Ccl4","Ccl5","Ccl8","Ccl18","Ccl19","Ccl21",
"Cxcl9","Cxcl10","Cxcl11","Cxcl13","Ccr7","Cxcr5","Sell","Lamp3")
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(M1@meta.data, M1@reductions$umap@cell.embeddings,Type=M1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=500,height=300)
ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=5,height=3)
ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=5,height=3)
ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
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
postscript("batch-score.eps",width=10,height=8)
ggplot(data,aes(x=batch,y=Sig.score,fill=batch,color=batch))+geom_violin()+labs(x="model",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

jpeg("cells-score2.jpg",width=800,height=700)
ggplot(data,aes(x=M1@active.ident,y=Sig.score,fill=M1@active.ident,color=M1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
ggplot(data,aes(x=M1@active.ident,y=Sig.score,fill=M1@active.ident,color=M1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
ggplot(data,aes(x=M1@active.ident,y=Sig.score,fill=M1@active.ident,color=M1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### late
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
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

#### 与三级淋巴组织形成有关的细胞趋化因子的Featureplot
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

jpeg("sig.score-TLS-early.jpg",width=1400,height=400)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	facet_grid(.~batch, scale='free')
	print(p)
dev.off()
pdf("sig.score-TLS-early.pdf",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(.~batch, scale='free')
	print(p)
dev.off()
setEPS()
postscript("sig.score-TLS-early.eps",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
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

jpeg("sig.score-TLS-late.jpg",width=1400,height=400)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	facet_grid(.~batch, scale='free')
	print(p)
dev.off()
pdf("sig.score-TLS-late.pdf",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(.~batch, scale='free')
	print(p)
dev.off()
setEPS()
postscript("sig.score-TLS-late.eps",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))+facet_grid(.~batch, scale='free')
	print(p)
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
per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage",
		"CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-TLS.pdf",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-TLS.jpeg",width=900,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-TLS.eps",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

########### zscore 与招募有关的细胞趋化因子
k<-M1
g<-c('Glycam1', 'Fut7', 'Gcnt1', 'Chst4', 'B3gnt3', 'Ccl21a')
gindex<-which(g %in% rownames(k@assays$RNA@data))
genes<-g[gindex]

s<-t(as.matrix(k@assays$RNA@data[genes,]))
k<-rowMeans(s)
data <- cbind(M1@meta.data, M1@reductions$umap@cell.embeddings,Type=M1@active.ident,k)
Sig.score <- (data$k-min(data$k))/(max(data$k)-min(data$k))

cor1<-cbind(data,Sig.score)
cor1$batch<-factor(cor1$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
head(cor1)

jpeg("sig.score.jpg",width=500,height=300)
ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
pdf("sig.score.pdf",width=5,height=3)
ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=0.8)+
  scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(cor1$Sig.score)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+facet_grid(.~batch, scale='free')
dev.off()
setEPS()
postscript("sig.score.eps",width=5,height=3)
ggplot(cor1,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
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
ggplot(data,aes(x=M1@active.ident,y=Sig.score,fill=M1@active.ident,color=M1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2.pdf",width=10,height=10)
ggplot(data,aes(x=M1@active.ident,y=Sig.score,fill=M1@active.ident,color=M1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2.eps",width=10,height=8)
ggplot(data,aes(x=M1@active.ident,y=Sig.score,fill=M1@active.ident,color=M1@active.ident))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
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
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score2-early.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score2-early.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
#	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### late
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
	geom_boxplot(width=0.3,color="black")+
	theme_classic() + theme(axis.text.x = element_text(face="bold",size=20)) +
#	stat_compare_means(comparisons = list(c("WT","AS"),c("WT","IG"),c("WT","AG"),c("WT","AAA"),c("WT","IH"),c("AS","IG"),c("AS","AG"),c("AS","AAA"),c("AS","IH"),c("AAA","IH"),c("AAA","IG"),c("AAA","AG"),c("IH","AG"),c("IH","IG"),c("AG","IG")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6) +
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
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
pdf("cells-score-late.pdf",width=10,height=10)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()
setEPS()
postscript("cells-score-late.eps",width=10,height=8)
ggplot(data,aes(x=Type,y=Sig.score,fill=Type,color=Type))+geom_violin()+labs(x="",y="Zscore",size=0.5)+geom_boxplot(width=0.1)+theme_classic() + theme(legend.position='none',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face="bold",size=12)) +
	stat_compare_means(comparisons = list(c("TREM2 Macrophage","Inflammation Macrophage"),c("TREM2 Macrophage","Res-like Macrophage"),c("TREM2 Macrophage","CD16 hi Monocyte"),c("TREM2 Macrophage","ISG hi Monocyte"),c("Inflammation Macrophage","Res-like Macrophage"),c("Inflammation Macrophage","CD16 hi Monocyte"),c("Inflammation Macrophage","ISG hi Monocyte"),c("Res-like Macrophage","CD16 hi Monocyte"),c("Res-like Macrophage","ISG hi Monocyte"),c("CD16 hi Monocyte","ISG hi Monocyte")),aes(label=paste0(..method..,"","p=",..p.format..)),size=6)+
	stat_summary(fun = "median", geom = "crossbar", mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.5,color="black")
dev.off()

### 与招募相关基因的featurePlot
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

jpeg("sig.score-Recuit-early.jpg",width=1400,height=400)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	facet_grid(.~batch, scale='free')
	print(p)
dev.off()
pdf("sig.score-Recuit-early.pdf",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(.~batch, scale='free')
	print(p)
dev.off()
setEPS()
postscript("sig.score-Recuit-early.eps",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
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

jpeg("sig.score-Recuit-late.jpg",width=1400,height=400)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	facet_grid(.~batch, scale='free')
	print(p)
dev.off()
pdf("sig.score-Recuit-late.pdf",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(.~batch, scale='free')
	print(p)
dev.off()
setEPS()
postscript("sig.score-Recuit-late.eps",width=18,height=4)
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2,colour=Sig.score)) + geom_point(size=1)+
	scale_colour_gradient2(low="blue", mid="blue", high="red", na.value = "grey50",
                         midpoint = median(data$Sig.score)) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))+facet_grid(.~batch, scale='free')
	print(p)
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
#perl sort.sig.time.pl newdata.txt > newdata2.txt

per<-read.table("newdata2.txt",header=T,sep="\t")
per$batch<-factor(per$batch,levels=c("WT","AS","AAA","IH","IG","AG"))
per$Type<-factor(per$Type,levels=c("TREM2 Macrophage","Inflammation Macrophage","Res-like Macrophage","CD16 hi Monocyte","ISG hi Monocyte"))
pdf("time-zscore-Rec.pdf",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()

jpeg("time-zscore-Rec.jpeg",width=900,height=300)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
setEPS()
postscript("time-zscore-Rec.eps",width=10,height=3)
ggplot(data = per,aes(x = batch, y = Sig.score,fill=group)) + 
  geom_bar(stat = 'identity', position = "dodge",width=0.6) + labs(x="") +
  theme_bw() + 
  theme(legend.title = element_blank(),axis.text.x=element_text(size=9,face = "bold")) +
  facet_grid(.~Type, scale='free')
dev.off()
