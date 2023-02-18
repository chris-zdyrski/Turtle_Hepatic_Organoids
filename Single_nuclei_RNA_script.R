#Renaming the genes
library(ape)
t.gff <- read.gff("/home/debosmitak/Desktop/RA/snRNASeq/GCF_000241765.4_Chrysemys_picta_BioNano-3.0.4_genomic.gff.gz")
genes<-gsub("^.*gene=(.*);.*$", "\\1", t.gff$attributes[t.gff$seqid=="NC_023890.1" & t.gff$type=="gene"], perl=T)
gene_changed<-paste("MT-", genes, sep="")

d<-read.csv("/home/debosmitak/Desktop/RA/snRNASeq/30-723912557/01_analysis/cellranger_count/Turtle-22-Liver-Organoid/filtered_feature_bc_matrix/features.tsv.gz", sep="\t")
for(i in 1:length(genes)){
	d[,2][which(d[,2]==genes[i])]<-gene_changed[i]
}

write.table(d, "/home/debosmitak/Desktop/RA/snRNASeq/30-723912557/01_analysis/cellranger_count/Turtle-22-Liver-Organoid/filtered_feature_bc_matrix/features.tsv", sep="\t", row.names=F)


library(dplyr)
library(Seurat)
library(patchwork)

# Load the liver dataset
liver.data <- Read10X(data.dir = "/home/debosmitak/Desktop/RA/snRNASeq/30-723912557/01_analysis/cellranger_count/Turtle-22-Liver-Organoid/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
liver <- CreateSeuratObject(counts = liver.data, project = "liver3k", min.cells = 3, min.features = 200)
liver





# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT-")
#liver[["percent.mt"]]<-colSums(liver@assays$RNA@counts[gene_changed, ])*100/colSums(liver@assays$RNA@counts)


# Visualize QC metrics as a violin plot
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

liver <- subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 40)

# Visualize QC metrics as a violin plot
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

liver <- subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 40)

# Visualize QC metrics as a violin plot
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

liver <- NormalizeData(liver)
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(liver), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(liver)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)
liver <- RunPCA(liver, features = VariableFeatures(object = liver))
# Examine and visualize PCA results a few different ways
print(liver[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(liver, dims = 1:2, reduction = "pca")
DimPlot(liver, reduction = "pca")
DimHeatmap(liver, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(liver, dims = 1:15, cells = 500, balanced = TRUE)



# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
liver <- JackStraw(liver, num.replicate = 100)
liver <- ScoreJackStraw(liver, dims = 1:20)
JackStrawPlot(liver, dims = 1:15)

ElbowPlot(liver)

liver <- FindNeighbors(liver, dims = 1:10)
liver <- FindClusters(liver, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(liver), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
liver <- RunUMAP(liver, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(liver, reduction = "umap")


# find all markers of cluster 2
cluster2.markers <- FindMarkers(liver, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(liver, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
    
#cluster0.markers <- FindMarkers(liver, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#VlnPlot(liver, features = c("MT-COX3", "ENSCPBG00000025423"))
# you can plot raw counts as well
#VlnPlot(liver, features = c("DIAPH3", "VEGFA"), slot = "counts", log = TRUE)
FeaturePlot(liver, features = c("TTK", "ENSCPBG00000004637", "SYN3", "POLA1", "NDRG1", "TXN", "SORBS2", "GCG"))
    
liver.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(liver, features = top10$gene) + NoLegend()
#library(readxl)
#liver.markers_cl<-read.csv("/home/debosmitak/Desktop/RA/snRNASeq/30-723912557/Biomarkers.csv",header = T)
#clusters<-read.csv("/home/debosmitak/Desktop/RA/snRNASeq/30-723912557/clusters.csv",header = T)
#Bcells<-read_excel("/home/debosmitak/Desktop/RA/snRNASeq/Supplementary Table 1.xlsx")
#CentralHep<-read_excel("/home/debosmitak/Desktop/RA/snRNASeq/Supplementary Table 1.xlsx", sheet=2)

#common_gene<-intersect(liver.markers$gene, Bcells$...1)
#cl_0<-which(liver.markers$cluster==0)
#lm_0<-c()
#for(j in 1:length(common_gene)){
#	if(common_gene[j] %in% liver.markers$gene[cl_0]){
#		lm_0<-c(lm_0, which(liver.markers$gene[cl_0] == common_gene[j]))
#	}	
#}

#Bc_0<-c()
#for(j in 1:length(lm_0)){
#	if(liver.markers$gene[lm_0[j]] %in% Bcells$...1){
#		Bc_0<-c(Bc_0, which(Bcells$...1 == liver.markers$gene[lm_0[j]]))
#	}	
#}

#Ch_0<-c()
#for(j in 1:length(lm_0)){
#	if(liver.markers$gene[lm_0[j]] %in% CentralHep$...1){
#		Ch_0<-c(Ch_0, which(CentralHep$...1 == liver.markers$gene[lm_0[j]]))
#	}	
#}

#cor(liver.markers$avg_log2FC[cl_0][lm_0],Bcells$summary.logFC[Bc_0])
#cor(liver.markers$avg_log2FC[cl_0][lm_0],CentralHep$summary.logFC[Ch_0])

#count<-matrix(0, 8,12)
#gene<-c()


#for(i in 1:dim(top10)[1]){
#	a<-array(0,12)
#	for(j in 1:dim(clusters)[2]){
#		if(top10$gene[i] %in% clusters[,j]){
#			a[j] = which(clusters[,j]==top10$gene[i])	
#		}
#	}
#	if(sum(a)!=0){
#		index<-ifelse((i%%10)!=0,(i%/%10)+1,i%/%10)
#		index_2<-which.min(a)
#		gene<-c(gene, top10$gene[i])
#		print(index_2)
#		count[index,index_2] = count[index,index_2]+1
#	}
#	if(i%%10==0){print("This cluster is done")}
#}

#gene_cl<-matrix(0, length(gene), 12)
#cluster_gene<-array(0, length(gene))
#for(i in 1:length(gene)){
#	for(j in 1:12){
#		if(gene[i]%in%clusters[,j]){
#			gene_cl[i,j]<-which(clusters[,j]==gene[i])
#		}
#	}
#	a<-which(top10$gene==gene_cl[i])
#	cluster_gene[i] = top10$cluster[i]
#		
#}

#count<-matrix(0, 8,12)
#gene<-c()


#for(i in 1:dim(liver.markers_cl)[1]){
#	a<-array(0,12)
#	for(j in 1:dim(clusters)[2]){
#		if(liver.markers_cl$gene[i] %in% clusters[,j]){
#			a[j] = which(clusters[,j]==liver.markers_cl$gene[i])	
#		}
#	}
#	if(sum(a)!=0){
#		index<-liver.markers_cl$cluster[i]+1
#		index_2<-which.min(a)
#		count[index,index_2] = count[index,index_2]+1
#	}
#}

#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#    "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(liver)
#liver <- RenameIdents(liver, new.cluster.ids)
#DimPlot(liver, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#clusters_mmc4<-read_excel("/home/debosmitak/Desktop/RA/mmc4.xlsx", sheet=2)
#gene_mmc4<-read_excel("/home/debosmitak/Desktop/RA/mmc4.xlsx", sheet=1)
#colnames(gene_mmc4)[8]<-"GENE"


#count<-matrix(0, 8, 98)
#liver.markers$cluster<-as.numeric(liver.markers$cluster)
#gene_com<-intersect(gene_mmc4$GENE, liver.markers$gene)
#for(i in 1:length(gene_com)){
#	index<-liver.markers$cluster[which(liver.markers$gene==gene_com[i])]
#	index_2<-gene_mmc4$cluster[which(gene_mmc4$GENE==gene_com[i])]
#	count[index,index_2] = count[index,index_2]+1
#}

#cl_list<-list(cluster1=c(), cluster2=c(), cluster3=c(), cluster4=c(), cluster5=c(), cluster6=c(), cluster7=c(), cluster8=c())
#for(i in 1:8){
#	a = max(count[i,])
#	index<-which(count[i,] == a)
#	cl_list[[i]]<-clusters_mmc4$`Cell Type`[index]
#}

#new.cluster.ids <- c(cl_list[[1]][1], cl_list[[2]][1], paste(cl_list[[3]][1], cl_list[[3]][2], cl_list[[3]][3], sep="+"), cl_list[[4]][1], cl_list[[5]][1], cl_list[[6]][1], cl_list[[7]][1], cl_list[[8]][1])
#names(new.cluster.ids) <- levels(liver)
#liver <- RenameIdents(liver, new.cluster.ids)
#DimPlot(liver, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
