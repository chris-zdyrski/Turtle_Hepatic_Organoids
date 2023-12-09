# Mucin Heatmaps
####To Make Heatmaps
#see: https://www.youtube.com/watch?v=ht1r34-ifVI&t=29s
#github: sanbomics/tutorial_complex_Heatmap.Rmd

library("DESeq2")
library("tidyverse")

load("DESeqResults.RData")

#look at the object of interest G = Hatchling Tissue vs Adult Tissue
head(resGSig)

#convert it to a data frame
df <- as.data.frame(resGSig) #change object relative to dataset of interest

df$names <- rownames(resGSig)
rownames(df) <- df$names

###New Portion
map <- read.delim(file = "../../Blast/BlastKeySuperRelaxed.tsv", header = T)

keys <- map$V2 #stringtieID or other ID
values <- map$V3 #uniprot hit

#make a list where each key is assigned a value
l <- list()
for (i in 1:length(keys)){
 l[keys[i]] <- values[i]
}

#for missing values
no_values <- setdiff(rownames(df), keys)

for (i in 1:length(no_values)){
 l[no_values[i]] <- 'NA'
}

#assign a new column to the original data frame containing the gene symbol
df$symbol <- unlist(l[rownames(df)], use.names = FALSE)

df$symbol2 <- paste(df$symbol, df$names, sep = "|")

###End of New Portion

#we don't want to filter the results anymore so we will just create a new object, df.top to simplify things
df.top <- df

#get data
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
mat0 <- assay(rlog_out)
mat <- mat0[rownames(df.top), rownames(coldata)]
colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (giving a z score) then transpose
colnames(mat.scaled) <- colnames(mat)

#keep mucin genes
rows_keep <- c('MSTRG.19448','MSTRG.19447','MSTRG.19444|MUC6')

#getting log2 value for each gene we are keeping
l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"

#getting mean value for each gene we are keeping
mean <- as.matrix(df.top[rows_keep,]$baseMean)
colnames(mean) <- "AvgExpr"

#load libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between blue/white/red for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AvgExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

#generate heatmap
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), height = unit(2, "cm")))

print("working")

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, column_labels = colnames(mat.scaled), name = "Z-score", cluster_columns = T, width = unit(20, "cm"))

print("working")

h2 <- Heatmap(l2_val, row_labels = df.top[rows_keep,]$symbol2, row_names_max_width = unit(10, "cm"), cluster_rows = F, name = "logFC", top_annotation = ha, col = col_logFC, width = unit(2, "cm"))

print("working")

h3 <- Heatmap(mean, row_labels = df.top[rows_keep,]$symbol2, row_names_max_width = unit(10, "cm"), cluster_rows = F, name = "AvgExpr", col = col_AvgExpr, width = unit(2, "cm"))

print("working")

h <- h1+h2+h3

pdf("./heatmapAgeTissue_Mucin2.pdf", width = 15, height = 16.5) #change filename relative to dataset of interest
print(h)
dev.off()

# Cluster Heatmaps
## 3A
####To Make Heatmaps
#see: https://www.youtube.com/watch?v=ht1r34-ifVI&t=29s
#github: sanbomics/tutorial_complex_Heatmap.Rmd

library("tidyverse")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")

tpm <- read.csv(file = "cpm_gene_count_matrix.csv", header = T, row.names = 1) #log2 transformed TPM values

head(tpm)

#convert it to a data frame and select three columns of interest
df <- as.data.frame(tpm[,c("Cse04LiverTissue","Cpi08LiverTissue","Cpi12LiverTissue","Cse04LiverOrganoid","Cpi08LiverOrganoid","Cpi12LiverOrganoid")])

head(df)

df$names <- rownames(tpm)
rownames(df) <- df$names

rownames(df)

#keep mucin genes
rows_keep <- c("MSTRG.31703|KRT18",
"MSTRG.31701|KRT8",
"MSTRG.18138|SPP1",
"MSTRG.28494|ANXA4",
"MSTRG.2425|AQP1",
"MSTRG.5983|PIGR",
"MSTRG.20655|HNF1B",
"MSTRG.31638|TDO2",
"MSTRG.25314|BCL2",
"MSTRG.7347|FHIT",
"MSTRG.16271|PPARGC1A",
"MSTRG.29760|SOX9",
"MSTRG.13973|ONECUT1",
"MSTRG.49518|POU5F1",
"MSTRG.15232|GATA6",
"MSTRG.17377|MAML2",
"MSTRG.849|NOTCH2",
"MSTRG.7725|MAML3",
"MSTRG.5251|RUNX2",
"MSTRG.20242|FOXO3",
"MSTRG.16248|RBPJ",
"MSTRG.10209|AOX1",
"MSTRG.38367|SLC27A2",
"MSTRG.26145|ABCA1",
"MSTRG.34400|GPC6",
"MSTRG.4049|HAO1",
"MSTRG.16483|FGFR3",
"MSTRG.34658|HES1",
"MSTRG.4075|JAG1",
"MSTRG.3977|FOXA2",
"MSTRG.8235|ONECUT2",
"MSTRG.7595|PROX1")

#getting mean value for each gene we are keeping
mat <- as.matrix(df[rows_keep,c("Cse04LiverTissue","Cpi08LiverTissue","Cpi12LiverTissue","Cse04LiverOrganoid","Cpi08LiverOrganoid","Cpi12LiverOrganoid")])
colnames(mat)

col_mat <- colorRamp2(c(-20,0,10), c("blue", "white","red")) #min value is ~ -19 max should be around 10 because it is log2(counts)

h <- Heatmap(mat, cluster_rows = F, column_labels = colnames(mat), name = "log2cpm", cluster_columns = F)

pdf("./heatmapCluster3A.pdf", width = 8.5, height = 11)
print(h)
dev.off()

## 3E
####To Make Heatmaps
#see: https://www.youtube.com/watch?v=ht1r34-ifVI&t=29s
#github: sanbomics/tutorial_complex_Heatmap.Rmd

library("tidyverse")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")

tpm <- read.csv(file = "cpm_gene_count_matrix.csv", header = T, row.names = 1) #log2 transformed TPM values

head(tpm)

#convert it to a data frame and select three columns of interest
df <- as.data.frame(tpm[,c("Cse04LiverTissue","Cpi08LiverTissue","Cpi12LiverTissue","Cse04LiverOrganoid","Cpi08LiverOrganoid","Cpi12LiverOrganoid")])

head(df)

df$names <- rownames(tpm)
rownames(df) <- df$names

rownames(df)

#keep mucin genes
rows_keep <- c("MSTRG.31703|KRT18",
"MSTRG.31701|KRT8",
"MSTRG.18138|SPP1",
"MSTRG.34680|CLDN1",
"MSTRG.29760|SOX9",
"MSTRG.2425|AQP1",
"MSTRG.32146|ELF3",
"MSTRG.25314|BCL2",
"MSTRG.5983|PIGR",
"MSTRG.21579|FOXA1",
"MSTRG.15232|GATA6",
"MSTRG.16119|FGFR2",
"MSTRG.16483|FGFR3",
"MSTRG.16184|CTBP2",
"MSTRG.8235|ONECUT2",
"MSTRG.13494|WWTR1",
"MSTRG.38434|CLDN3",
"MSTRG.5869|GATA4",
"MSTRG.1379|ALCAM",
"MSTRG.1777|CXCR4",
"MSTRG.33449|TBX3",
"MSTRG.28494|ANXA4",
"MSTRG.20655|HNF1B",
"MSTRG.20242|FOXO3",
"MSTRG.13491|TM4SF4",
"MSTRG.49518|POU5F1",
"MSTRG.3977|FOXA2",
"MSTRG.13054|TF",
"MSTRG.23931|HAL",
"MSTRG.26254|HNF4A",
"MSTRG.41271|ASGR1",
"MSTRG.7595|PROX1")

#getting mean value for each gene we are keeping
mat <- as.matrix(df[rows_keep,c("Cse04LiverTissue","Cpi08LiverTissue","Cpi12LiverTissue","Cse04LiverOrganoid","Cpi08LiverOrganoid","Cpi12LiverOrganoid")])
colnames(mat)

col_mat <- colorRamp2(c(-20,0,10), c("blue", "white","red")) #min value is ~ -19 max should be around 10 because it is log2(counts)

h <- Heatmap(mat, cluster_rows = F, column_labels = colnames(mat), name = "log2cpm", cluster_columns = F)

pdf("./heatmapCluster3E.pdf", width = 8.5, height = 11)
print(h)
dev.off()
