---
title: "bulk RNA-seq analysis for the paper: Characterization of the First Turtle Organoids"
author: "Thea B. Gessler"
output: html_document
---

# Initial Check with FASTQC

Check the quality of the reads with the program fastqc. You can list as many pairs of fastq files as needed in a single command. They should be listed as pairs R1, R2, R1, R2, etc.

```
module load fastqc/0.11.7-3flwcvl

fastqc -o fastQC -t 36 -f fastq --contaminants adapters.txt R1.fastq.gz R2.fastq.gz 
```
Once the command has been run, review the output html files for quality checks. Not all quality checks need to be passed for RNA-seq data given some of its fundamental features.

# Read Trimming with Trimmomatic

Trim the reads to remove adapter sequences and low quality bases using the program trimmomatic. You will need to specify names for the output - both paired reads and unpaired reads after trimming. Enter a separate command for each library. It can be a litte tricky to correctly specify the adapter file, so if possible use one that is provided.

```
module load trimmomatic/0.39-da5npsr

trimmomatic PE -threads 36 -phred33 \
R1.fastq.gz R2.fastq.gz \
R1_pairedOutputR1.fastq R1_unpairedOutputR1.fastq R2_pairedOutputR2.fastq R2_unpairedOutputR2.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

# Second Check with FASTQC

Check the quality of the reads after trimming with the program fastqc. You can list as many pairs of fastq files as needed in a single command. They should be listed as pairs R1, R2, R1, R2, etc.

```
module load fastqc/0.11.7-3flwcvl

fastqc -o fastQC_PostTrim -t 36 -f fastq --contaminants adapters.txt \
R1_pairedOutputR1.fastq R2_pairedOutputR2.fastq
```
Once the command has been run, review the output html files for quality checks. Not all quality checks need to be passed for RNA-seq data given some of its fundamental features.

# Read Mapping with GSNAP

First generate a database with indexing files for mapping. You will specify the location and the name of this database directory. The location of the fasta file used for mapping will also need to be provided.

```
module load gcc/10.2.0-zuvaafu
module load gmap-gsnap/2021-03-08-jb3nha6

gmap_build -D </path/to/mapping/database> -d <database_name> /path/to/fasta/GCF_000241765.4_Chrysemys_picta_BioNano-3.0.4_genomic.fna
```

Next, map the reads and calculate some basic statistics with samtools. It is most efficient to map the libraries independently so multiple jobs can be run simultaneously. 

```
module load samtools/1.10-py3-xuj7ylj
module load gcc/10.2.0-zuvaafu
module load gmap-gsnap/2021-03-08-jb3nha6

fastq1=/path/to/R1_pairedOutputR1.fastq
fastq2=/path/to/R2_pairedOutputR2.fastq

gsnap -D </path/to/mapping/database> -d <database_name> -t 36 -N 1 --pairexpect=300 -A sam -o file.sam $fastq1 $fastq2

samtools view -bS file.sam | samtools sort -o file.bam
samtools flagstat file.bam > file.stat.out
rm file.sam
```

# Assembly and Annotation with StringTie

Stringtie generates initial independent assemblies that are then merged and used to calculate transcript abundances. First generate the initial assemblies for each library. This can be done with a for loop to save on scripting. The for loop is made as a .sh file which is then run using sbatch scripts.

```
#!/bin/bash
for x in /path/to/files*.bam
        do
                stringtie $x -o ./${x}.gtf -p 36 -G /path/to/gff/GCF_000241765.4_Chrysemys_picta_BioNano-3.0.4_genomic.gff --rf
        done
```

Run the for loop.

```
module load stringtie/1.3.4a-3fmlrwa

sh stringtie.sh
```

Next, merge the assemblies again using a for loop for simplicity. Specify the name of the output file (-o) and provide a list of the gtf file locations what will be included in the merge (CpiGtf.list).

```
#!/bin/bash

#stringtie --merge -G /path/to/gff/GCF_000241765.4_Chrysemys_picta_BioNano-3.0.4_genomic.gff -o outTemp.gtf gtf.list
```

Run the for loop. 

```
module load stringtie/1.3.4a-3fmlrwa

sh stringtieMerge.sh
```

Calculate abundances of the transcripts. A separate line of code will need to be run for each library. There are probably more efficient ways to do this. You will need to specify the output file name (-o), the name of the directory abundance files will be located in (-b) and the name of the previously generated gtf file (-G). The command looks like this for each library:

```
stringtie /path/to/file.bam -o outFinal.gtf -b <Directory> -p 36 -G outTemp.gtf --rf -e
```

# ERCC Spike-in Analysis

This followed a similar process as above but for the spike in reads separately. The files for the spike-in reads were provided by the company.

## Generate a database for mapping

```
module load gcc/10.2.0-zuvaafu
module load gmap-gsnap/2021-03-08-jb3nha6

gmap_build -D /path/to/mapping/dir/ERCC -d ERCC_Dir ERCC92.fa
```

## Map the reads

We map all the reads because the spike-ins are contained within the whole set. But in this case our reference is the reference file for the spike-in reads, not the genome. Again you will want to have a separate mapping script for each library.

```
module load samtools/1.10-py3-xuj7ylj
module load gcc/10.2.0-zuvaafu
module load gmap-gsnap/2021-03-08-jb3nha6

fastq1=/path/to/R1_pairedOutputR1.fastq
fastq2=/path/to/R2_pairedOutputR2.fastq

gsnap -D /path/to/mapping/directory/ERCC -d ERCC_Dir -t 36 -N 1 --pairexpect=300 -A sam -o fileERCC.sam $fastq1 $fastq2

samtools view -bS fileERCC.sam | samtools sort -o fileERCC.bam
samtools flagstat fileERCC.bam > fileERCC.stat.out
rm fileERCC.sam
```

## Assemble with Stringtie

### Initial Assembly

```
module load stringtie/1.3.4a-3fmlrwa

sh stringtie.sh
```

stringtie.sh - assembles the reads using the spike-in reference annotation file. Here we add the -e option because we don't want StringTie to discover new transcripts. We only want the known spike-in transcripts.

```
#!/bin/bash
for x in ./*.bam
        do
                stringtie $x -o ./${x}.gtf -p 36 -G ERCC92.gtf --rf -e
        done
```

### Merge

```
module load stringtie/1.3.4a-3fmlrwa

sh stringtieMerge.sh
```

stringtieMerge.sh

```
#!/bin/bash

stringtie --merge -G ERCC92.gtf -o ./mergeERCC.gtf ERCC.list
```

### Calculate Abundance

Include a command for each library.

```
module load stringtie/1.3.4a-3fmlrwa

stringtie fileERCC.bam -o finalFileERCC.gtf -b fileERCC -p 36 -G mergeERCC.gtf --rf -e
```

### Generate the counts to use for normalizaiton in differential expression analysis

Run the prepDE.py script provided by the StringTie developers. samples.txt lists the location of the final gtf files containing abundances

```
module load python/2.7.18-2ut3ogj

python prepDE.py -i samples.txt -l 150 -c
```

This will give a gene count matrix and a transcript count matrix file. These should be identical becasue we did not allow new transcript discovery.

We next join the gene count matrix with a key file derived from stringtie output to change the IDs. This key file comes from the abundance calculation output and is just an isolation of the transcript data containing only the chromosome ID (chr) and the ID generated by StringTie (gene_id).

```{R}
library("tidyverse")

geneMatrix <- read.csv("gene_count_matrix.csv")
key <- read.delim("../../ERCC/StringtieNoNewTranscripts/Cpi08LiverOrganoidERCC/key.tsv")
join1 <- left_join(geneMatrix, key, by = "gene_id")
res <- join1[,c(14,2:13)]
write.csv(res, file = "gene_count_matrix_annotated.csv", quote = F, row.names = F)
```

To run the Rscript:
```
module load gcc/7.3.0-xegsmw4
module load r/4.0.2-py3-icvulwq

Rscript --vanilla join.R
```

Finally, remove the header from gene_count_matrix_annotated.csv to give gene.csv so we can easily concatenate it to our gene count matrix from the primary assembly.

# Differential Expression Analysis with DESeq2

Note: I tested several normalization schemes: with the spike-in (ultimately used), without the spike-in, and with the spike-in but not special normalization specified. They all gave very similar results.

## prep the gene count matrix

Run prepDE.py, provided by the StringTie developers to convert the StringTie abundance output into a gene count matrix. samples.txt contains the names of the libraries in column one and the path to abunance files (gtf files) in column two.

```
module load python/2.7.18-2ut3ogj

python prepDE.py -i samples.txt -l 150 -c
```

Next concatenate the ERCC counts to this newly generated gene count matrix to include the spike-ins.

```
cat gene_count_matrix.csv gene.csv > > joint_gene_count_matrix.csv
```
## Run DESeq2 to calculate differential expression

```{R}
library("DESeq2")

#####read in data
counts <- as.matrix(read.csv("joint_gene_count_matrix.csv", row.names = 1)) # reads in the counts matrix which contains liver stuff and ercc stuff
counts <- round(counts) #make integers
head(counts)
tail(counts)
ercc <- read.delim(file = "ERCC_Control_IDs_ForDEA.txt", header = T) #names of spike-in controls
head(ercc)
coldata <- read.delim("targets.txt", header=TRUE, sep="\t", row.names=1) #classifies libraries
coldata$Age <- relevel(coldata$Age, ref="Hatchling") #hatchling is baseline so positive fold change indicates adult upregulation
coldata$SampleType <- relevel(coldata$SampleType, ref="Tissue") #tissue is baseline so positive fold change indicates organoid upregulation
coldata$Condition <- relevel(coldata$Condition, ref="B") #sets hatchling tissue as baseline, not really a necessary step
coldata

#####control vector
controlTest <- subset(counts, rownames(counts) %in% ercc$ID)

head(controlTest)

control <- if(subset(counts, rownames(counts) %in% ercc$ID)) {
TRUE
} else {
FALSE
}

head(control)
tail(control)

#####Run DESeq2
ddsCount <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~Age * SampleType) #full model
keep <- rowSums(counts(ddsCount)) >= 10
dds <- ddsCount[keep,]
dds <- estimateSizeFactors(ddsCount, type = "ratio", controlGenes=control) #direct it to use the ERCC spike-ins to estimate size factors
#Michael Love, DESeq2 creator, said to give controlGenes the names of the spike-ins in this function when using spike-ins 
#see: https://support.bioconductor.org/p/9143354/
dds <- DESeq(dds) #run DESeq2 (may need additional options)
resultsNames(dds)

#contrasts - I worked this out on paper using the notes from Meiling, so this is correct
resA <- results(dds, contrast=c(0,0,2,1)) #main effect SampleType 
resB <- results(dds, contrast=c(0,2,0,1)) #main effect Age 
resC <- results(dds, contrast=c(0,0,0,1)) #interaction coefficient
resD <- results(dds, contrast=c(0,0,1,0)) #simple effect of SampleType in Hatchlings
resE <- results(dds, contrast=c(0,0,1,1)) #simple effect of SampleType in Adults
resF <- results(dds, contrast=c(0,1,0,1)) #simple effect of Age within Organoids
resG <- results(dds, contrast=c(0,1,0,0)) #simple effect of Age within Tissues

#rearrange output
resASig <- resA[order(resA$padj),]
resBSig <- resB[order(resB$padj),]
resCSig <- resC[order(resC$padj),]
resDSig <- resD[order(resD$padj),]
resESig <- resE[order(resE$padj),]
resFSig <- resF[order(resF$padj),]
resGSig <- resG[order(resG$padj),]

#####write results
write.table(as.data.frame(resASig), file="DESeq2SampleTypeMainEffectResults.tsv", sep="\t", quote=F, row.names=T)
write.table(as.data.frame(resBSig), file="DESeq2AgeMainEffectResults.tsv", sep="\t", quote=F, row.names=T)
write.table(as.data.frame(resCSig), file="DESeq2InteractionEffectResults.tsv", sep="\t", quote=F, row.names=T)
write.table(as.data.frame(resDSig), file="DESeq2SampleTypeHatchlingsResults.tsv", sep="\t", quote=F, row.names=T)
write.table(as.data.frame(resESig), file="DESeq2SampleTypeAdultsResults.tsv", sep="\t", quote=F, row.names=T)
write.table(as.data.frame(resFSig), file="DESeq2AgeOrganoidsResults.tsv", sep="\t", quote=F, row.names=T)
write.table(as.data.frame(resGSig), file="DESeq2AgeTissuesResults.tsv", sep="\t", quote=F, row.names=T)

#####save an RData to use for heatmap generation
save(coldata, dds, resASig, resBSig, resCSig, resDSig, resESig, resFSig, resGSig, file = "DESeqResults.RData")
```

Run the Rscript.

```
module load gcc/7.3.0-xegsmw4
module load r/4.0.2-py3-icvulwq
module load r-deseq2/1.24.0-py3-jiqp45b

Rscript --vanilla deseq2.R
```

# Generation of Heatmaps

Create heatmaps  based on the differential expression results. This is an example script for all the heatmaps. For different datasets, swap out the object for the dataset of interest and change the filenames for the output accordingly.

```{R}
#### To Make Heatmaps 
# see: https://www.youtube.com/watch?v=ht1r34-ifVI&t=29s
# github: sanbomics/tutorial_complex_Heatmap.Rmd

library("DESeq2")
library("tidyverse")

load("DESeqResults.RData")

#look at the object of interest
#resFSig should be replaced for each data set.
head(resFSig)
write.table(rownames(resFSig), file = "rownamesresFSig.txt", sep = "\t", quote = F)

#convert it to a data frame
df <- as.data.frame(resFSig) #some rownames deemed invalid and pushed to NA causing problems later

df$names <- rownames(resFSig)
rownames(df) <- df$names

###New Portion - determining the IDs for genes resulting from a Uniprot blast
map <- read.delim(file = "../../Blast/BlastKeySuperRelaxed.tsv", header = T)

head(map)

keys <- map$V2 #stringtieID or other ID
values <- map$V3 #uniprot hit

head(keys)
head(values)

#make a list where each key is assigned a value
l <- list()
for (i in 1:length(keys)){
 l[keys[i]] <- values[i]
}

#for missing values
no_values <- setdiff(rownames(df), keys)

class(rownames(df))
class(keys)

head(no_values, 50)
length(no_values)

for (i in 1:length(no_values)){
 l[no_values[i]] <- 'NA'
}

#assign a new column to the original data frame containing the gene symbol
df$symbol <- unlist(l[rownames(df)], use.names = FALSE)

head(df, 12)

df$symbol2 <- paste(df$symbol, df$names, sep = "|")

head(df, 12)

###End of New Portion

#filter the top results from the data frame of differential expression results
df.top <- df[(df$padj < 0.05) & (df$baseMean > 50) & (abs(df$log2FoldChange) > 1.5),]

tail(df.top)

#order dataset by log2FoldChange column
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]

tail(df.top)

#inspect the object
head(df.top)
dim(df.top)
head(rownames(df.top))
write.table(rownames(df.top), file = "rownamesDfTop.txt", sep = "\t", quote = F)

head(rownames(coldata))

#drop NAs
df.top <- drop_na(df.top) # no clue where they came from but they are messing up later steps

#get data
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object

#inspect data
head(rlog_out)
dim(rlog_out)
head(assay(rlog_out))
mat0 <- assay(rlog_out)
head(mat0)
dim(mat0)
rownames(mat0)
write.table(rownames(mat0), file = "rownamesMat0.txt", sep = "\t", quote = F)

#which part of the code is causing the error?
mat1 <- mat0[,rownames(coldata)]
mat2 <- mat0[rownames(df.top),] #this one - due to NAs introduced
mat <- mat0[rownames(df.top), rownames(coldata)]

#mat <- assay(rlog_out)[rownames(df.top), rownames(coldata)] #sig genes x samples
colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (giving a z score) then transpose
colnames(mat.scaled) <- colnames(mat)

#keep top 25 genes and bottom 25 genes
num_keep <- 25
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))

#who is in rows_keep
write.table(rows_keep, file = "Top25Bottom25Genes.txt", sep = "\t", quote = F)

l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val) <- "logFC"

mean <- as.matrix(df.top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean) <- "AvgExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between blue/white/red for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AvgExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, column_labels = colnames(mat.scaled), name = "Z-score", cluster_columns = T, width = unit(20, "cm"))

h2 <- Heatmap(l2_val, row_labels = df.top$symbol2[rows_keep], row_names_max_width = unit(8, "cm"), cluster_rows = F, name = "logFC", top_annotation = ha, col = col_logFC, width = unit(2, "cm"))

h3 <- Heatmap(mean, row_labels = df.top$symbol2[rows_keep], row_names_max_width = unit(8, "cm"), cluster_rows = F, name = "AvgExpr", col = col_AvgExpr, width = unit(2, "cm"))

h <- h1+h2+h3

#output file names should be replaced for each data set.
tiff("./heatmapAgeOrganoid_v1.tiff", res = 300, width = 5000, height = 5500)
pdf("./heatmapAgeOrganoid_v1.pdf", width = 15, height = 16.5)
print(h)
dev.off()
```

# Generation of PCA plots

Code below generates the plot that is for the _Chrysemys picta_ liver transcriptome.

```{R}
####EdgeR transformation
library(edgeR) #calls edgeR 
x <- joint_gene_count_matrix #assign the gene count matrix (need to read this in) to object x
x <- round(x)
group <- factor(c("hatchOrganoid","hatchTissue","hatchOrganoid","hatchTissue","hatchOrganoid","hatchTissue","adultOrganoid","adultTissue","adultOrganoid","adultTissue","adultOrganoid","adultTissue")) #defines the groups, in this case the replicates
y <- DGEList(counts=x, group=group) #makes a data object DGEList used by EdgeR
keep <- filterByExpr(y) #filters out low count genes-genes are kept if there are at least 2 columns (because group size is 2) with worthwhile counts but filters without respect to which group counts fall into (default minimum count is 10, default minimum total count is 15) see details under filterByExpr in vingettes
y <- y[keep, , keep.lib.sizes=FALSE] #will recalculate the library sizes after filtering which is recommended by the developers but not required
y <- calcNormFactors(y) #normalizes library size so that a small number of highly expressed genes in a single sample don't make genes look falsely down regulated uses TMM to generate effective library size
y
tmm <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.0001) #according to numerous boards this is how you get this. Also because log = TRUE it is doing the log2 transformation at the same time
#note that this was not normalized with the ERCC spike-in

#make a PCA plot
tData <- t(tmm) #transpose data
treatment <- rownames(tData) #set treatment variable - should match names of groups
pcaData <- prcomp(tData) #get principle component data
summary(pcaData)
pcScores <- pcaData$x #isolate the values used for the pca plot
pcaDF <- data.frame(pcScores, treatment) #make a data frame

#variance
pcaData$sdev^2 / sum(pcaData$sdev^2)

#define your color and shape variables so you know what is what
colors <- c("Cpi08LiverOrganoid" = "mediumorchid4",
            "Cpi08LiverTissue" = "mediumorchid4",
            "Cpi09LiverOrganoid" = "#91bfdb",
            "Cpi09LiverTissue" = "#91bfdb",
            "Cpi11LiverOrganoid" = "#fee090",
            "Cpi11LiverTissue" = "#fee090",
            "Cpi12LiverOrganoid" = "springgreen4",
            "Cpi12LiverTissue" = "springgreen4",
            "Cpi13LiverOrganoid" = "white",
            "Cpi13LiverTissue" = "white",
            "Cpi14LiverOrganoid" = "black",
            "Cpi14LiverTissue" = "black")

shapes <- c("Cpi08LiverOrganoid" = 21,
            "Cpi08LiverTissue" = 22,
            "Cpi09LiverOrganoid" = 21,
            "Cpi09LiverTissue" = 22,
            "Cpi11LiverOrganoid" = 21,
            "Cpi11LiverTissue" = 22,
            "Cpi12LiverOrganoid" = 21,
            "Cpi12LiverTissue" = 22,
            "Cpi13LiverOrganoid" = 21,
            "Cpi13LiverTissue" = 22,
            "Cpi14LiverOrganoid" = 21,
            "Cpi14LiverTissue" = 22)

#make the plot
#start by generating a .tiff file or a .pdf file or whichever you like
#you will provide it with dimensions and resolution, these here should work fine
tiff("pcaLiverTranscriptome.tiff", width = 7, height = 5, units = 'in', res = 300)
pdf("pcaLiverTranscriptome.pdf", width = 7, height = 5)
ggplot(pcaDF, aes(PC1, PC2, shape=treatment, fill=treatment)) + #treatment is referencing the same factor variable as before
  geom_point(aes(fill=factor(treatment),shape=factor(treatment)),size=5) +
  theme(panel.background = element_rect(fill = "white", color = "black")) +
  scale_fill_manual(values=colors) + #fyi these colors are colorblind and grayscale friendly
  scale_shape_manual(values=shapes) +
  theme(axis.title = element_text(size = 16, face = "bold")) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold")) +
  xlab(label = "PC1 47.10%") +
  ylab(label = "PC2 9.32%")
dev.off() #dev.off tells R you are done providing code to be incorporated into the plot
```

Code below generates the PCA plot that includes both _Chrysemys picta_ and _Chelydra serpentina_ liver. Note this draws from a different assembly that also incorporated the _C. serpentina_ data. It was prepared identically as described above.

```{R}
#### PCA only Liver - includes Cpi and Cse
####EdgeR transformation
library(edgeR) #calls edgeR 
a <- gene_count_matrix[,c(1:10,13:16),] #select the columns for cpi liver and cse liver
a <- round(a)
groupLiver <- factor(c("hatchOrganoid","hatchTissue","hatchOrganoid","hatchTissue","hatchOrganoid","hatchTissue","adultOrganoid","adultTissue","adultOrganoid","adultTissue","adultOrganoid","adultTissue","CseOrganoid","CseTissue")) #defines the groups, in this case the replicates
b <- DGEList(counts=a, group=groupLiver) #makes a data object DGEList used by EdgeR
keepLiver <- filterByExpr(b) #filters out low count genes-genes are kept if there are at least 2 columns (because group size is 2) with worthwhile counts but filters without respect to which group counts fall into (default minimum count is 10, default minimum total count is 15) see details under filterByExpr in vingettes
b <- b[keepLiver, , keep.lib.sizes=FALSE] #will recalculate the library sizes after filtering which is recommended by the developers but not required
b <- calcNormFactors(b) #normalizes library size so that a small number of highly expressed genes in a single sample don't make genes look falsely down regulated uses TMM to generate effective library size
b
tmmLiver <- cpm(b, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.0001) #according to numerous boards this is how you get this. Also because log = TRUE it is doing the log2 transformation at the same time
#note that this was not normalized with the ERCC spike-in

#make a PCA plot
tDataLiver <- t(tmmLiver) #transpose data
treatmentLiver <- rownames(tDataLiver) #set treatment variable - should match names of groups
pcaDataLiver <- prcomp(tDataLiver) #get principle component data
summary(pcaDataLiver)
pcScoresLiver <- pcaDataLiver$x #isolate the values used for the pca plot
pcaDFLiver <- data.frame(pcScoresLiver, treatmentLiver) #make a data frame

#variance
pcaDataLiver$sdev^2 / sum(pcaDataLiver$sdev^2)

#define your color and shape variables so you know what is what
colorsLiver <- c("Cpi08LiverOrganoid" = "mediumorchid4",
               "Cpi08LiverTissue" = "mediumorchid4",
               "Cpi09LiverOrganoid" = "#91bfdb",
               "Cpi09LiverTissue" = "#91bfdb",
               "Cpi11LiverOrganoid" = "#fee090",
               "Cpi11LiverTissue" = "#fee090",
               "Cpi12LiverOrganoid" = "springgreen4",
               "Cpi12LiverTissue" = "springgreen4",
               "Cpi13LiverOrganoid" = "white",
               "Cpi13LiverTissue" = "white",
               "Cpi14LiverOrganoid" = "black",
               "Cpi14LiverTissue" = "black",
               "Cse04LiverOrganoid" = "sienna1",
               "Cse04LiverTissue" = "sienna1")

shapesLiver <- c("Cpi08LiverOrganoid" = 21,
               "Cpi08LiverTissue" = 22,
               "Cpi09LiverOrganoid" = 21,
               "Cpi09LiverTissue" = 22,
               "Cpi11LiverOrganoid" = 21,
               "Cpi11LiverTissue" = 22,
               "Cpi12LiverOrganoid" = 21,
               "Cpi12LiverTissue" = 22,
               "Cpi13LiverOrganoid" = 21,
               "Cpi13LiverTissue" = 22,
               "Cpi14LiverOrganoid" = 21,
               "Cpi14LiverTissue" = 22,
               "Cse04LiverOrganoid" = 21,
               "Cse04LiverTissue" = 22)

pdf("pcaTranscriptomeAllLiver.pdf", width = 7, height = 5)
ggplot(pcaDFLiver, aes(PC1, PC2, shape=treatmentLiver, fill=treatmentLiver)) + #treatment is referencing the same factor variable as before
  geom_point(aes(fill=factor(treatmentLiver),shape=factor(treatmentLiver)),size=5) +
  theme(panel.background = element_rect(fill = "white", color = "black")) +
  scale_fill_manual(values=colorsLiver) + #fyi these colors are colorblind and grayscale friendly
  scale_shape_manual(values=shapesLiver) +
  theme(axis.title = element_text(size = 16, face = "bold")) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold")) +
  xlab(label = "PC1 34.76%") +
  ylab(label = "PC2 17.49%")
dev.off() #dev.off tells R you are done providing code to be incorporated into the plot
```

# Generation of Venn diagrams

```{R}
#read in your gene count matrix

#isolate the organoid and tissue libraries for Cpi separately
cpiOrg2 <- joint_gene_count_matrix[,c(1,3,5,7,9,11)] #Cpi liver only assembly
cpiTis2 <- joint_gene_count_matrix[,c(2,4,6,8,10,12)] #Cpi liver only assembly

adultOrgAlt <- joint_gene_count_matrix[,c(7,9,11)] #cpi liver only assembly
adultTisAlt <- joint_gene_count_matrix[,c(8,10,12)] #cpi liver only assembly

hatchOrgAlt <- joint_gene_count_matrix[,c(1,3,5)] #cpi liver only assembly
hatchTisAlt <- joint_gene_count_matrix[,c(2,4,6)] #cpi liver only assembly

#sum the counts across rows (genes)
sumOrg2 <- rowSums(cpiOrg2)
sumTis2 <- rowSums(cpiTis2)

sumOrg7 <- rowSums(adultOrgAlt)
sumTis7 <- rowSums(adultTisAlt)

sumOrg8 <- rowSums(hatchOrgAlt)
sumTis8 <- rowSums(hatchTisAlt)

#filter out lowly expressed genes within each group
#this step is critical because otherwise we would see the same things for both
#that would be because the count matrices were collectivley constructed, not because they are the same
#I picked this more stringent threshold so we can be more confident that we are looking at genes that are actually expressed
keepOrg2 <- sumOrg2 > 100
keepTis2 <- sumTis2 > 100

keepOrg7 <- sumOrg7 > 50
keepTis7 <- sumTis7 > 50

keepOrg8 <- sumOrg8 > 50
keepTis8 <- sumTis8 > 50

#This step does the actual filtering
geneOrg2 <- cpiOrg2[keepOrg2,]
geneTis2 <- cpiTis2[keepTis2,]

geneOrg7 <- adultOrgAlt[keepOrg7,]
geneTis7 <- adultTisAlt[keepTis7,]

geneOrg8 <- hatchOrgAlt[keepOrg8,]
geneTis8 <- hatchOrgAlt[keepTis8,]

library(VennDiagram)

#now we need the row names from the filtered results so we can create the venn diagrams
namesOrg2 <- row.names(geneOrg2)
namesTis2 <- row.names(geneTis2)

namesOrg7 <- row.names(geneOrg7)
namesTis7 <- row.names(geneTis7)

namesOrg8 <- row.names(geneOrg8)
namesTis8 <- row.names(geneTis8)

venn.diagram(x=list(namesOrg2, namesTis2),
             filename = "./Output/OrganoidVsTissueOverlapCpiLiverAssembly.tiff",
             height = 10, width = 10, units = "in",
             main = "Chyrsemys picta Liver Expressed Genes Overlap",
             main.fontface = "bold",
             main.fontfamily = "sans",
             main.cex = 2.25,
             category.names = c("Organoid Genes", "Tissue Genes"),
             col=c("springgreen4","mediumorchid4"),
             fill=c("springgreen4","mediumorchid4"),
             margin=0.1,
             fontfamily = "sans",
             cat.cex = 1.25, 
             cat.fontfamily = "sans",
             cat.fontface = "bold",
             cex = 1.75,
             fontface = "bold",
             cat.pos = c(200,160))

venn.diagram(x=list(namesOrg8, namesTis8, namesOrg7, namesTis7),
             filename = "./Output/CpiAdultVsHatchlingOrganoidAndTissue.tiff",
             height = 10, width = 10, units = "in",
             main = "C. picta Orgaoind and Tissue Expressed Genes Overlap",
             main.fontface = "bold",
             main.fontfamily = "sans",
             main.cex = 2.25,
             category.names = c("Hatchling Organoid", "Hatchling Tissue", "Adult Organoid", "Adult Tissue"),
             col=c("springgreen4","mediumorchid4", "#fee090", "#91bfdb"),
             fill=c("springgreen4","mediumorchid4", "#fee090", "#91bfdb"),
             margin=0.1,
             fontfamily = "sans",
             cat.cex = 1.25, 
             cat.fontfamily = "sans",
             cat.fontface = "bold",
             cex = 1.75,
             fontface = "bold")
```

# Enrichment Analysis with PantherDB

## Data preparation

### Get TPM values for each gene from the assembly.

```{R}
library("tximport") #load packages
library("tidyverse")
cpiDir <- ("/path/to/directory/Stringtie") #this first part sets up the reading in of files
list.files(cpiDir)
samples <- c("Cpi08LiverOrganoid","Cpi08LiverTissue","Cpi09LiverOrganoid","Cpi09LiverTissue","Cpi11LiverOrganoid","Cpi11LiverTissue","Cpi12LiverOrganoid","Cpi12LiverTissue","Cpi13LiverOrganoid","Cpi13LiverTissue","Cpi14LiverOrganoid","Cpi14LiverTissue")
files <- file.path(cpiDir, samples, "t_data.ctab")
names(files) <- paste0(samples)
all(file.exists(files))
files
tmp <- read_tsv(files[1]) #this reads in the transcript file so IDs can be obtained
tx2gene <- tmp[, c("t_name","gene_id")]
write.table(tx2gene, file = "tx2gene.tsv", quote = FALSE, sep = "\t")
write.table(tmp, file = "tmp.tsv", quote = FALSE, sep = "\t")
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM") #This does the consolidating of different transcripts into gene models
names(txi)
head(txi$counts)
counts <- txi$counts
write.table(counts, file = "countsGeneTPM.tsv", quote=FALSE, sep ="\t") #Write the output
```

### Isolate transcript sequences with gffread.

```
/path/to/gffread-0.12.7.Linux_x86_64/gffread -w StringtieLiverTranscriptsMerge.fasta -g /path/to/genome/assembly/GCF_000241765.4_Chrysemys_picta_BioNano-3.0.4_genomic.fna --keep-genes CpiLiverAssembly.gtf
```

### Translate assembly transcript sequences into peptide sequences.

```
module load transdecoder/5.5.0-lr3a445

TransDecoder.LongOrfs -t StringtieLiverTranscripts.fasta --gene_trans_map sort_gene_transcript_map.tsv

echo "LongOrfs Done"

TransDecoder.Predict -t StringtieLiverTranscripts.fasta --single_best_only 

echo "Predict Done"
```

### Search translated assembly transcript sequences against panther hidden markov models.

```
module load gcc/7.3.0-xegsmw4
module load perl/5.34.0-d6utzgh
module load hmmer/3.3-py3-openmpi3-brlxcli

export PATH="/path/to/PANTHER/pantherScore2.2/lib:$PATH"

./pantherScore2.2.pl -l ./target/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/ -D B -V -i StringtieLiverTranscripts.fasta.transdecoder.pep -o TranslatedCDS.txt -n
```

### Modify output file.

```
#remove the peptide info

awk 'BEGIN {FS = "\t"} {split($1,a,".p"); {print a[1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' TranslatedCDS.txt > translatedTranscripts.txt

add headers to columns in translatedTranscripts.txt
```

### Join gene TPM values to corresponding transcript IDs and join to pather results. Filter to nonredundant set of gene best hits.

```{R}
library("tidyverse")

map <- read.delim("sort_gene_transcript_map.tsv") #map of geneID to corresponding transcriptID

transcript <- read.delim("translatedTranscripts.txt") #results of mapping transcripts to pantherIDs

gene <- read.delim("countsGeneTPM.tsv") #gene level TPM values

join1 <- distinct(left_join(gene, map, by = "GeneID"))

head(join1)

join2 <- distinct(left_join(join1, transcript, by = "TranscriptID"))

head(join2)

mygenes0 <- replace_na(join2, list(PantherID="NOHIT",GeneName="NOHIT",HmmEvalue="NOHIT",HmmBitscore="NOHIT",AlignRange="NOHIT"))

head(mygenes0)

mygenes1 <- arrange(mygenes0, desc(HmmBitscore), HmmEvalue)

mygenes2 <- !duplicated(mygenes1$GeneID)

mygenes3 <- mygenes1[mygenes2,]

mygenes4 <- !duplicated(mygenes3$PantherID, incomparables = "NOHIT")

mygenes5 <- mygenes3[mygenes4,]

mygenes <- mygenes5[,c(1,15,2:13)]

head(mygenes)

write.table(mygenes, file = "new_gene_counts_matrix.csv", sep = ",", quote = F, row.names = F)
```

### Prepare the data for searching in the pather database.

```{R}
library("tidyverse")

gene <- read.csv("new_gene_counts_matrix.csv")

a <- filter(gene[,c(1:3)], PantherID != "NOHIT")
b <- filter(gene[,c(1:2,4)], PantherID != "NOHIT")
c <- filter(gene[,c(1:2,5)], PantherID != "NOHIT")
d <- filter(gene[,c(1:2,6)], PantherID != "NOHIT")
e <- filter(gene[,c(1:2,7)], PantherID != "NOHIT")
f <- filter(gene[,c(1:2,8)], PantherID != "NOHIT")
g <- filter(gene[,c(1:2,9)], PantherID != "NOHIT")
h <- filter(gene[,c(1:2,10)], PantherID != "NOHIT")
i <- filter(gene[,c(1:2,11)], PantherID != "NOHIT")
j <- filter(gene[,c(1:2,12)], PantherID != "NOHIT")
k <- filter(gene[,c(1:2,13)], PantherID != "NOHIT")
l <- filter(gene[,c(1:2,14)], PantherID != "NOHIT")

length(a)
length(b)
length(c)
length(d)
length(e)
length(f)
length(g)
length(h)
length(i)
length(j)
length(k)
length(l)

a1 <- arrange(a, PantherID, desc(Cpi08LiverOrganoid))
b1 <- arrange(b, PantherID, desc(Cpi08LiverTissue))
c1 <- arrange(c, PantherID, desc(Cpi09LiverOrganoid))
d1 <- arrange(d, PantherID, desc(Cpi09LiverTissue))
e1 <- arrange(e, PantherID, desc(Cpi11LiverOrganoid))
f1 <- arrange(f, PantherID, desc(Cpi11LiverTissue))
g1 <- arrange(g, PantherID, desc(Cpi12LiverOrganoid))
h1 <- arrange(h, PantherID, desc(Cpi12LiverTissue))
i1 <- arrange(i, PantherID, desc(Cpi13LiverOrganoid))
j1 <- arrange(j, PantherID, desc(Cpi13LiverTissue))
k1 <- arrange(k, PantherID, desc(Cpi14LiverOrganoid))
l1 <- arrange(l, PantherID, desc(Cpi14LiverTissue))

a2 <- !duplicated(a1$PantherID, incomparables = "NOHIT")
b2 <- !duplicated(b1$PantherID, incomparables = "NOHIT")
c2 <- !duplicated(c1$PantherID, incomparables = "NOHIT")
d2 <- !duplicated(d1$PantherID, incomparables = "NOHIT") 
e2 <- !duplicated(e1$PantherID, incomparables = "NOHIT")
f2 <- !duplicated(f1$PantherID, incomparables = "NOHIT")
g2 <- !duplicated(g1$PantherID, incomparables = "NOHIT")
h2 <- !duplicated(h1$PantherID, incomparables = "NOHIT")
i2 <- !duplicated(i1$PantherID, incomparables = "NOHIT")
j2 <- !duplicated(j1$PantherID, incomparables = "NOHIT")
k2 <- !duplicated(k1$PantherID, incomparables = "NOHIT")
l2 <- !duplicated(l1$PantherID, incomparables = "NOHIT")

head(a1)
head(b1)
head(a2)
head(b2)

a3 <- a1[a2,]
b3 <- b1[b2,]
c3 <- c1[c2,]
d3 <- d1[d2,]
e3 <- e1[e2,]
f3 <- f1[f2,]
g3 <- g1[g2,]
h3 <- h1[h2,]
i3 <- i1[i2,]
j3 <- j1[j2,]
k3 <- k1[k2,]
l3 <- l1[l2,]

head(a3)
head(b3)

length(a3)
length(b3)
length(c3)
length(d3)
length(e3)
length(f3)
length(g3)
length(h3)
length(i3)
length(j3)
length(k3)
length(l3)

write.table(a3, file = "Cpi08Organoid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(b3, file = "Cpi08Tissue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(c3, file = "Cpi09Organoid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(d3, file = "Cpi09Tissue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(e3, file = "Cpi11Organoid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(f3, file = "Cpi11Tissue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(g3, file = "Cpi12Organoid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(h3, file = "Cpi12Tissue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(i3, file = "Cpi13Organoid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(j3, file = "Cpi13Tissue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(k3, file = "Cpi14Organoid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(l3, file = "Cpi14Tissue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

### Results were then manually searched on pantherdb.org using statistical enrichment test.

# Visualization of Enrichment Results with REVIGO

Prepare results for input into REVIGO by finding the common set of positively enriched terms among each set of triplicates.

```{R}
#Prepare the enrichment results for input into REVIGO and also to convert into venn diagrams

cpi08OrgBP <- read.delim("Enrichment_Data/Cpi08OrganoidBiolProcessEdit.txt")
cpi08OrgCC <- read.delim("Enrichment_Data/Cpi08OrganoidCellCompEdit.txt")
cpi08OrgMF <- read.delim("Enrichment_Data/Cpi08OrganoidMolFuncEdit.txt")

cpi08TisBP <- read.delim("Enrichment_Data/Cpi08TissueBiolProcessEdit.txt")
cpi08TisCC <- read.delim("Enrichment_Data/Cpi08TissueCellCompEdit.txt")
cpi08TisMF <- read.delim("Enrichment_Data/Cpi08TissueMolFuncEdit.txt")

cpi09OrgBP <- read.delim("Enrichment_Data/Cpi09OrganoidBiolProcessEdit.txt")
cpi09OrgCC <- read.delim("Enrichment_Data/Cpi09OrganoidCellCompEdit.txt")
cpi09OrgMF <- read.delim("Enrichment_Data/Cpi09OrganoidMolFuncEdit.txt")

cpi09TisBP <- read.delim("Enrichment_Data/Cpi09TissueBiolProcessEdit.txt")
cpi09TisCC <- read.delim("Enrichment_Data/Cpi09TissueCellCompEdit.txt")
cpi09TisMF <- read.delim("Enrichment_Data/Cpi09TissueMolFuncEdit.txt")

cpi11OrgBP <- read.delim("Enrichment_Data/Cpi11OrganoidBiolProcessEdit.txt")
cpi11OrgCC <- read.delim("Enrichment_Data/Cpi11OrganoidCellCompEdit.txt")
cpi11OrgMF <- read.delim("Enrichment_Data/Cpi11OrganoidMolFuncEdit.txt")

cpi11TisBP <- read.delim("Enrichment_Data/Cpi11TissueBiolProcessEdit.txt")
cpi11TisCC <- read.delim("Enrichment_Data/Cpi11TissueCellCompEdit.txt")
cpi11TisMF <- read.delim("Enrichment_Data/Cpi11TissueMolFuncEdit.txt")

cpi12OrgBP <- read.delim("Enrichment_Data/Cpi12OrganoidBiolProcessEdit.txt")
cpi12OrgCC <- read.delim("Enrichment_Data/Cpi12OrganoidCellCompEdit.txt")
cpi12OrgMF <- read.delim("Enrichment_Data/Cpi12OrganoidMolFuncEdit.txt")

cpi12TisBP <- read.delim("Enrichment_Data/Cpi12TissueBiolProcessEdit.txt")
cpi12TisCC <- read.delim("Enrichment_Data/Cpi12TissueCellCompEdit.txt")
cpi12TisMF <- read.delim("Enrichment_Data/Cpi12TissueMolFuncEdit.txt")

cpi13OrgBP <- read.delim("Enrichment_Data/Cpi13OrganoidBiolProcessEdit.txt")
cpi13OrgCC <- read.delim("Enrichment_Data/Cpi13OrganoidCellCompEdit.txt")
cpi13OrgMF <- read.delim("Enrichment_Data/Cpi13OrganoidMolFuncEdit.txt")

cpi13TisBP <- read.delim("Enrichment_Data/Cpi13TissueBiolProcessEdit.txt")
cpi13TisCC <- read.delim("Enrichment_Data/Cpi13TissueCellCompEdit.txt")
cpi13TisMF <- read.delim("Enrichment_Data/Cpi13TissueMolFuncEdit.txt")

cpi14OrgBP <- read.delim("Enrichment_Data/Cpi14OrganoidBiolProcessEdit.txt")
cpi14OrgCC <- read.delim("Enrichment_Data/Cpi14OrganoidCellCompEdit.txt")
cpi14OrgMF <- read.delim("Enrichment_Data/Cpi14OrganoidMolFuncEdit.txt")

cpi14TisBP <- read.delim("Enrichment_Data/Cpi14TissueBiolProcessEdit.txt")
cpi14TisCC <- read.delim("Enrichment_Data/Cpi14TissueCellCompEdit.txt")
cpi14TisMF <- read.delim("Enrichment_Data/Cpi14TissueMolFuncEdit.txt")

######################

#Make a list of each biological process file
listBP <- list(cpi08OrgBP, cpi08TisBP, cpi09OrgBP, cpi09TisBP, cpi11OrgBP, cpi11TisBP, 
               cpi12OrgBP, cpi12TisBP, cpi13OrgBP, cpi13TisBP, cpi14OrgBP, cpi14TisBP)

#make a list of each cellular component file
listCC <- list(cpi08OrgCC, cpi08TisCC, cpi09OrgCC, cpi09TisCC, cpi11OrgCC, cpi11TisCC, 
               cpi12OrgCC, cpi12TisCC, cpi13OrgCC, cpi13TisCC, cpi14OrgCC, cpi14TisCC)

#make a list of each molecular function file
listMF <- list(cpi08OrgMF, cpi08TisMF, cpi09OrgMF, cpi09TisMF, cpi11OrgMF, cpi11TisMF, 
               cpi12OrgMF, cpi12TisMF, cpi13OrgMF, cpi13TisMF, cpi14OrgMF, cpi14TisMF)

#load tidyverse
library(tidyverse)

#for loops to use a regular expression to isolate the GO ID and append it as an additional element to the sublist
for(i in 1:length(listBP)) {
  listBP[[i]]$GO <- str_extract(listBP[[i]]$PANTHER.GO.Slim.Biological.Process, "GO:[:digit:]+")
}

for(i in 1:length(listCC)) {
  listCC[[i]]$GO <- str_extract(listCC[[i]]$PANTHER.GO.Slim.Cellular.Component, "GO:[:digit:]+")
}

for(i in 1:length(listMF)) {
  listMF[[i]]$GO <- str_extract(listMF[[i]]$PANTHER.GO.Slim.Molecular.Function, "GO:[:digit:]+")
}

#now we want just the up-enriched terms and we only want the GO ID and the FDR value

#isolate the direction of enrichment the fdr and the GO ID

#filter the list based on one of the column values
listBP2 <- lapply(listBP, filter, overUnder == "+")
listCC2 <- lapply(listCC, filter, overUnder == "+")
listMF2 <- lapply(listMF, filter, overUnder == "+")

#initialize an empty list
listBP3 <- list()
listCC3 <- list()
listMF3 <- list()

#populate the list with only three of the original columns
for(i in 1:length(listBP2)) {
listBP3[[i]] <- listBP2[[i]][,c(6,5)]
}

for(i in 1:length(listCC2)) {
  listCC3[[i]] <- listCC2[[i]][,c(6,5)]
}

for(i in 1:length(listMF2)) {
  listMF3[[i]] <- listMF2[[i]][,c(6,5)]
}

#get dataframes of all the items of the list
names(listBP3) <- c("CpiO08BP", "CpiT08BP", "CpiO09BP", "CpiT09BP", "CpiO11BP", "CpiT11BP", "CpiO12BP", "CpiT12BP", "CpiO13BP", "CpiT13BP", "CpiO14BP", "CpiT14BP")
names(listCC3) <- c("CpiO08CC", "CpiT08CC", "CpiO09CC", "CpiT09CC", "CpiO11CC", "CpiT11CC", "CpiO12CC", "CpiT12CC", "CpiO13CC", "CpiT13CC", "CpiO14CC", "CpiT14CC")
names(listMF3) <- c("CpiO08MF", "CpiT08MF", "CpiO09MF", "CpiT09MF", "CpiO11MF", "CpiT11MF", "CpiO12MF", "CpiT12MF", "CpiO13MF", "CpiT13MF", "CpiO14MF", "CpiT14MF")

lapply(names(listBP3), function(i) write.table(listBP3[[i]], file=paste(i, "txt", sep="."), sep = "\t", quote = F, row.names = F))
lapply(names(listCC3), function(i) write.table(listCC3[[i]], file=paste(i, "txt", sep="."), sep = "\t", quote = F, row.names = F))
lapply(names(listMF3), function(i) write.table(listMF3[[i]], file=paste(i, "txt", sep="."), sep = "\t", quote = F, row.names = F))


#create a list with a path to all the files you want to load
#note: files are now alphabetical so it goes all organoid all tissue, rather than in matched pairs
temp <- list.files(path = "./", pattern = "*BP.txt", full.names = T)
tempCC <- list.files(path = "./", pattern = "CC.txt", full.names = T)
tempMF <- list.files(path = "./", pattern = "MF.txt", full.names = T)

#use lapply to load those files
myfiles <- lapply(temp, read.delim)
myfilesCC <- lapply(tempCC, read.delim)
myfilesMF <- lapply(tempMF, read.delim)

#biological process hatchling organoid
joinBP1 <- full_join(myfiles[[1]], myfiles[[2]], by = "GO", suffix = c(".08", ".09"))
joinBP2 <- drop_na(full_join(joinBP1, myfiles[[3]], by = "GO"))
joinBP2$meanFDR <- rowMeans(joinBP2[,c(2:4)])
bpHatchOrganoid <- joinBP2[,c(1,5)]

#biological process hatchling tissue
joinBP3 <- full_join(myfiles[[7]], myfiles[[8]], by = "GO", suffix = c(".08", ".09"))
joinBP4 <- drop_na(full_join(joinBP3, myfiles[[9]], by = "GO"))
joinBP4$meanFDR <- rowMeans(joinBP4[,c(2:4)])
bpHatchTissue <- joinBP4[,c(1,5)]

#biological process adult organoid
joinBP5 <- full_join(myfiles[[4]], myfiles[[5]], by = "GO")
joinBP6 <- drop_na(full_join(joinBP5, myfiles[[6]], by = "GO"))
joinBP6$meanFDR <- rowMeans(joinBP6[,c(2:4)])
bpAdultOrganoid <- joinBP6[,c(1,5)]

#biological process adult tissue
joinBP7 <- full_join(myfiles[[10]], myfiles[[11]], by = "GO")
joinBP8 <- drop_na(full_join(joinBP7, myfiles[[12]], by = "GO"))
joinBP8$meanFDR <- rowMeans(joinBP8[,c(2:4)])
bpAdultTissue <- joinBP8[,c(1,5)]

#cellular component hatchling organoid
joinCC1 <- full_join(myfilesCC[[1]], myfilesCC[[2]], by = "GO")
joinCC2 <- drop_na(full_join(joinCC1, myfilesCC[[3]], by = "GO"))
joinCC2$meanFDR <- rowMeans(joinCC2[,c(2:4)])
ccHatchOrganoid <- joinCC2[,c(1,5)]

#cellular component hatchling tissue
joinCC3 <- full_join(myfilesCC[[7]], myfilesCC[[8]], by = "GO")
joinCC4 <- drop_na(full_join(joinCC3, myfilesCC[[9]], by = "GO"))
joinCC4$meanFDR <- rowMeans(joinCC4[,c(2:4)])
ccHatchTissue <- joinCC4[,c(1,5)]

#cellular componenet adult organoid
joinCC5 <- full_join(myfilesCC[[4]], myfilesCC[[5]], by = "GO")
joinCC6 <- drop_na(full_join(joinCC5, myfilesCC[[6]], by = "GO"))
joinCC6$meanFDR <- rowMeans(joinCC6[,c(2:4)])
ccAdultOrganoid <- joinCC6[,c(1,5)]

#cellular component adult tissue
joinCC7 <- full_join(myfilesCC[[10]], myfilesCC[[11]], by = "GO")
joinCC8 <- drop_na(full_join(joinCC7, myfilesCC[[12]], by = "GO"))
joinCC8$meanFDR <- rowMeans(joinCC8[,c(2:4)])
ccAdultTissue <- joinCC8[,c(1,5)]

#molecular function hatchling organoid
joinMF1 <- full_join(myfilesMF[[1]], myfilesMF[[2]], by = "GO")
joinMF2 <- drop_na(full_join(joinMF1, myfilesMF[[3]], by = "GO"))
joinMF2$meanFDR <- rowMeans(joinMF2[,c(2:4)])
mfHatchOrganoid <- joinMF2[,c(1,5)]

#molecular function hatchling tissue
joinMF3 <- full_join(myfilesMF[[7]], myfilesMF[[8]], by = "GO")
joinMF4 <- drop_na(full_join(joinMF3, myfilesMF[[9]], by = "GO"))
joinMF4$meanFDR <- rowMeans(joinMF4[,c(2:4)])
mfHatchTissue <- joinMF4[,c(1,5)]

#molecular function adult organoid
joinMF5 <- full_join(myfilesMF[[4]], myfilesMF[[5]], by = "GO")
joinMF6 <- drop_na(full_join(joinMF5, myfilesMF[[6]], by = "GO"))
joinMF6$meanFDR <- rowMeans(joinMF6[,c(2:4)])
mfAdultOrganoid <- joinMF6[,c(1,5)]

#molecular function adult tissue
joinMF7 <- full_join(myfilesMF[[10]], myfilesMF[[11]], by = "GO")
joinMF8 <- drop_na(full_join(joinMF7, myfilesMF[[12]], by = "GO"))
joinMF8$meanFDR <- rowMeans(joinMF8[,c(2:4)])
mfAdultTissue <- joinMF8[,c(1,5)]

#write the output
#make a list
listAll <- list(bpHatchOrganoid,bpHatchTissue,bpAdultOrganoid,bpAdultTissue,ccHatchOrganoid,ccHatchTissue,ccAdultOrganoid,ccAdultTissue,mfHatchOrganoid,mfHatchTissue,mfAdultOrganoid,mfAdultTissue)
#designate list names
names(listAll) <- c("bpHatchOrganoid","bpHatchTissue","bpAdultOrganoid","bpAdultTissue","ccHatchOrganoid","ccHatchTissue","ccAdultOrganoid","ccAdultTissue","mfHatchOrganoid","mfHatchTissue","mfAdultOrganoid","mfAdultTissue")
#write
lapply(names(listAll), function(i) write.table(listAll[[i]], file=paste(i, "txt", sep="."), sep = "\t", quote = F, row.names = F))
```

Code for generating REVIGO plots was downloaded from the REVIGO website and run locally. Code for figures is available upon request.

# Acknowledgments

I would be remiss if I did not thank all of the people who have made excellent tutorials as well as those who have provided _kind_ answers to questions posted online.
