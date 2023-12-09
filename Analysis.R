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
liver.data <- Read10X(data.dir = "01_analysis/cellranger_count/Turtle-22-Liver-Organoid/filtered_feature_bc_matrix/")

liver <- CreateSeuratObject(counts = liver.data, project = "liver3k", min.cells = 3, min.features = 200)
liver

liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT-")
#liver[["percent.mt"]]<-colSums(liver@assays$RNA@counts[gene_changed, ])*100/colSums(liver@assays$RNA@counts)


# Visualize QC metrics as a violin plot
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

liver <- subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 40)

# Visualize QC metrics as a violin plot
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

liver <- subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 40)

# Visualize QC metrics as a violin plot
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


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

print(liver[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(liver, dims = 1:2, reduction = "pca")
DimPlot(liver, reduction = "pca")
DimHeatmap(liver, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(liver, dims = 1:15, cells = 500, balanced = TRUE)


liver <- JackStraw(liver, num.replicate = 100)
liver <- ScoreJackStraw(liver, dims = 1:20)
JackStrawPlot(liver, dims = 1:15)

ElbowPlot(liver)

liver <- FindNeighbors(liver, dims = 1:10)
liver <- FindClusters(liver, resolution = 0.5)
head(Idents(liver), 5)
liver <- RunUMAP(liver, dims = 1:10)
DimPlot(liver, reduction = "umap")

cluster2.markers <- FindMarkers(liver, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(liver, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
    
FeaturePlot(liver, features = c("DIAPH3", "ENSCPBG00000004984", "ENSCPBG00000004637", "SYN3", "POLA1", "TXN", "SORBS2", "ENSCPBG00000007381"))
    
liver.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(liver, features = top10$gene) + NoLegend()

dat1<-read_excel("clusters from andrews.xlsx", sheet=1)
dat2<-read_excel("clusters from andrews.xlsx", sheet=2)
cd_genes1<-as.vector(dat1[-c(1,11,22,32,45),3])[[1]]
cd_genes2<-as.vector(dat2[-c(1,2),2])[[1]]
DotPlot(object = liver, features = cd_genes1)
DotPlot(object = liver, features = cd_genes2)


