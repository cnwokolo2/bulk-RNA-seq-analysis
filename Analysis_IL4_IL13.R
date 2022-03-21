#Investigate transcriptional differences between IL4 organoids treated for 7 days and IL13 glands treated for 3days

library("DESeq2")
library("apeglm")
library("pasilla")
library("matrixStats")
library("reshape")
library("ggplot2")
library("ggrepel")
library("DEGreport")
library("RColorBrewer")
library("pheatmap")
library("EnhancedVolcano")
library("gplots")



setwd("/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/")
directory<-"/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/"

samplenames<-c(list.files(pattern="*.counts.txt"))
print(samplenames)
samplenames<-print(gsub("[:outputcnsx::_::.09o:]", "", samplenames))

IL13trt<-grep("13trt",list.files(directory), value=T)
IL4trt<-grep("4trt",list.files(directory), value=T)
IL13controls<-grep("m13control",list.files(directory), value=T)
IL4controls<-grep("m4control",list.files(directory), value=T)

condition<-c("IL13 treated", "IL13 treated", "IL13 treated","IL4 treated", "IL4 treated", "IL4 treated", 
             "IL13 control", "IL13 control", "IL13 control", "IL4 control", "IL4 control", "IL4 control")

Table<-data.frame(sampleNames=samplenames, filename=c( IL13trt, IL4trt, IL13controls, IL4controls), condition=condition)
print(Table)
Table$condition <- factor(Table$condition)

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = Table,
                                        directory = directory,
                                        design = ~condition)
ddsHTseq


colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('IL13 treated', 'IL4 treated', 'IL13 control', 'IL4 control'))

dds <- DESeq(ddsHTseq)
res <- results(dds)
res <- res[order(res$padj),] #order most significant to least
head(res)

#write.csv(res, file="CSV/avgexp.csv") #log2FC pvadj
sum(res$padj < 0.05, na.rm=TRUE) #11283

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
#write.csv(resdata, file="CSV/DESeq2-results-with-normalized.csv")
#write.table(as.data.frame(counts(dds),normalized=T), file = 'CSV/DESeq2_normalized_counts.txt', sep = '\t')
mcols(res, use.names = T)
#write.csv(as.data.frame(mcols(res, use.name = T)),file = "CSV/DESeq2-test-conditions.csv")

#minimize outliers
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.001,
             cleaned = results(ddsClean)$padj < 0.001)
addmargins(tab)
#write.csv(addmargins(tab), file = 'CSV/DESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
#write.csv(as.data.frame(resClean),file = 'CSV/DESeq2-replaceoutliers-results.csv')

#log normalization methods
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
#write.csv(assay(vsd), file = "CSV/DESeq2_vsd_normalized.csv")
#write.csv(assay(rld), file = "CSV/DESeq2_rlog_normalized.csv")

#visualizations
#to account for low count genes log fold changes
resultsNames(dds)
resLFC1 <- lfcShrink(dds, coef="condition_IL4.treated_vs_IL13.treated", type="apeglm")
resLFC1
#write.csv(as.data.frame(resLFC1),file = "CSV/condition_IL4.treated_vs_IL13.treated.csv")
plotMA(resLFC1, ylim=c(-6,6), alpha= 0.001, main = "IL4 Treated vs IL13 Treated RNAseq")
dev.copy(png, "./Plots/condition_IL4.treated_vs_IL13.treated<0.001.png")
dev.off()
#comparing non-shrunken values
res.noshr <- results(dds, name="condition_IL4.treated_vs_IL13.treated")
plotMA(res.noshr, ylim = c(-6, 6))

#labelling topgene on plot
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#bottomgene
bottomGene <- rownames(res)[which.max(res$padj)]
with(res[bottomGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, bottomGene, pos=2, col="dodgerblue")
})

#histogram
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#plot counts for specific genes
plotCounts(dds, gene=which.min(res$padj), intgroup="condition") #gene with smallest p-val = 0610005C13Rik
plotCounts(dds, gene=which.max(res$padj), intgroup="condition") #gene with largest p-val = Sh3gl3

opar <- par(mfrow = c(2,5), mex=0.6, mar=c(4,4,3,2) +0.3) #helps plot graphs in a 2x2 format
opar <- par(mfrow = c(2,5), mex=0.9, c(5, 4, 4, 2) + 0.2) #helps plot graphs in a 2x2 format

plotCounts(dds, gene="Gkn3", intgroup="condition")
plotCounts(dds, gene="Tff2", intgroup="condition")
plotCounts(dds, gene="Muc6", intgroup="condition")
plotCounts(dds, gene="Muc5ac", intgroup="condition")
plotCounts(dds, gene="Gif", intgroup="condition")
plotCounts(dds, gene="Chil4", intgroup="condition")
plotCounts(dds, gene="Stmn1", intgroup="condition")
plotCounts(dds, gene="Mki67", intgroup="condition")
plotCounts(dds, gene="Birc5", intgroup="condition")
plotCounts(dds, gene="Spink4", intgroup="condition")
plotCounts(dds, gene="Iqgap3", intgroup="condition")
plotCounts(dds, gene="Dmbt1", intgroup="condition")
plotCounts(dds, gene="Cd44", intgroup="condition")
plotCounts(dds, gene="Tff3", intgroup="condition")
plotCounts(dds, gene="Gkn2", intgroup="condition")

par(mfrow=c(1,1)) #closes the graphic parameter

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(0,25,50,100))

#heatmap of countmatrix
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

sampleDists <- dist(t(assay(dds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$condition, dds$sizeFactor, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap(log2(counts(dds, normalized = TRUE) [rownames(dds)%in%features,]),
        col = hmcol, scale = "row",
        Romv = TRUE, Colv = FALSE,
        dendrogram= "row",
        trace= "none",
        margin = c(4,6), cexRow=0.5, cexCol =1, keysize =1)

#volcanoplots
threshold <- c(res$padj<0.05)
length(which(threshold)) #3326
res$threshold <- threshold

#unfiltered volcanoplots
spemstem <- c("Hmgb2","Ube2c","Pclaf","Birc5","Cks2","Tff1","Gkn2","H2afz","Tubb5","Top2a","Tuba1b","Stmn1","Sptssb","Gm3776", "Nusap1",
              "Ube2s","Cdk1","Ccna2","Dpcr1", "Cenpa", "Pbk","Pttg1", "Mgst3","Smc4","Gkn1","Ran", "Mki67","Cks1b", "Rrm2", "Lgals4","Spc25",
              "Tubb4b","Smc2","Tmpo",'Tk1',"Muc5ac","Ccnb1", "Cdc20", "Ptma","Jpt1","H2afx","Arl6ip1")
genestolabel <- c(spemstem)
EnhancedVolcano(resLFC1, lab=rownames(res), x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 1.3, selectLab=genestolabel, drawConnectors = FALSE, title = 'IL4 Treated vs IL13 Treated RNAseq')
EnhancedVolcano(resLFC1, lab=rownames(res), x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 2.0, drawConnectors = FALSE, title = 'IL4 Treated vs IL13 Treated RNAseq' )

#volcanoplot of lfc2
EnhancedVolcano(resLFC1,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'IL4 Treated vs IL13 Treated RNAseq',
                pCutoff = 10e-5,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

#Special heatmap
#four conditions (consolidated)
genesofinterest <- c("Gkn3", "Muc6", "Spink4", "Mki67", "Birc5", "Stmn1", "Cd44", "Iqgap3", "Dmbt1", "Muc5ac", "Gif")
#define genes to label
genestolabel <- c(genesofinterest)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/CSV/DESeq2_rlog_normalized_consolidated.csv", header = T)
#setup matrix for heatmap
values <- data.frame(data[,c(2:5)])
values <- as.matrix(values)
rownames(values) <- data$gene
#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)
heatmap.2(subset, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "row", main = "Heatmap of genes in IL4 and IL13 treatments",
          col = bluered(100), trace = "none", density.info = "none", keysize = 0.9, add.expr = T, cexRow = 1.2)


#Special heatmap of non-consolidated date 12 samples
genesofinterest <- c("Gkn3", "Muc6", "Spink4", "Mki67", "Birc5", "Stmn1", "Cd44", "Iqgap3", "Dmbt1", "Muc5ac", "Gif", "Chil4")
#define genes to label
genestolabel <- c(genesofinterest)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/CSV/DESeq2_rlog_normalized copy.csv", header = T)
#setup matrix for heatmap
values <- data.frame(data[,c(2:13)])
values <- as.matrix(values)
rownames(values) <- data$gene
#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)
#heatmap.2: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
heatmap.2(subset, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "row", main = "Heatmap of genes in IL4 and IL13 treatments",
          col = bluered(100), trace = "none", density.info = "none", keysize = 0.9, add.expr = T, cexRow = 1.2)

#Special heatmap of consolidated trt samples
genesofinterest <- c("Gkn3", "Muc6", "Spink4", "Mki67", "Birc5", "Stmn1", "Cd44", "Iqgap3", "Dmbt1", "Muc5ac", "Gif", "Chil4")
#define genes to label
genestolabel <- c(genesofinterest)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/CSV/DESeq2_rlog_normalized_IL4_IL13.csv", header = T)
#setup matrix for heatmap
values <- data.frame(data[,c(2:3)])
values <- as.matrix(values)
rownames(values) <- data$gene
#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)
#heatmap.2: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
heatmap.2(subset, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "row", main = "Heatmap of genes in IL4 and IL13 treatments",
          col = bluered(100), trace = "none", density.info = "none", keysize = 0.9, add.expr = T, cexRow = 1.2)

#Special heatmap of non-consolidated trt samples
genesofinterest <- c("Gkn3", "Muc6", "Spink4", "Mki67", "Birc5", "Stmn1", "Cd44", "Iqgap3", "Dmbt1", "Muc5ac", "Gif", "Chil4")
#define genes to label
genestolabel <- c(genesofinterest)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/CSV/DESeq2_rlog_normalized_trt.csv", header = T)
#setup matrix for heatmap
values <- data.frame(data[,c(2:7)])
values <- as.matrix(values)
rownames(values) <- data$gene
#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)
#heatmap.2: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
heatmap.2(subset, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "row", main = "Heatmap of genes in IL4 and IL13 treatments",
          col = bluered(100), trace = "none", density.info = "none", keysize = 0.9, add.expr = T, cexRow = 1.2)


#volcanoplot
genestolabel <- c(genesofinterest)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/CSV/condition_IL4.treated_vs_IL13.treated copy.csv", header = T)

names(genestolabel)[genestolabel == 'red'] <- 'Genes of Interest'
#Log2 Fold Change cutoff modified to 2.0 (FC=1), and p.value cutoff modified to 0.05

#generate volcano plot labeled by pathway
EnhancedVolcano(data, lab=data$gene, x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 2.0, xlim=c(-15,15), ylim=c(0, 200, na.rm = TRUE), selectLab=genestolabel, drawConnectors = FALSE, title = 'IL4 Treated vs IL13 Treated RNAseq',
                labSize=4.0 , pointSize = 1.0)
EnhancedVolcano(data, lab=data$gene, x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 2.0, xlim=c(-13,13), ylim=c(0, 350, na.rm = TRUE), drawConnectors = FALSE, title = 'IL4 Treated vs IL13 Treated RNAseq',
                labSize=3.0 , pointSize = 1.0)

#converting log2FC to FC
data.transform <- data
data.transform <- data.transform[,c(1,3)]
transformed <- 2^(data.transform[,2])
data[,7] <- transformed
names(data)[7] <- "FC"
final <- data[,c(1,2,3,7,4,5,6)]
write.csv(final, file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/Mouse_IL13_IL4/CSV/DESeq2_IL4_treated_vs_IL13.treated_results_FC.csv")



