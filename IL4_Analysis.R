#Human microglia2 (some stimulated by tcells)

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

setwd("/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/IL4_organoids/data_counts")
directory<-"/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/IL4_organoids/data_counts"

samplenames<-c(list.files(pattern="*.counts.txt"))
print(samplenames)
samplenames<-print(gsub("[:outputcnsx::_::.09o:]", "", samplenames))

controls<-grep("mcontrol",list.files(directory), value=T)
IL4treat<-grep("m4",list.files(directory), value=T)
condition<-c("IL4 treated", "IL4 treated", "IL4 treated", "control", "control", "control")

Table<-data.frame(sampleNames=samplenames, filename=c(IL4treat, controls), condition=condition)
print(Table)
Table$condition <- factor(Table$condition)

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = Table,
                                        directory = directory,
                                        design = ~condition)
ddsHTseq


colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('control','IL4 treated'))

dds <- DESeq(ddsHTseq)
res <- results(dds)
res <- res[order(res$padj),] #order most significant to least
head(res)

#write.csv(res, file="avgexp.csv") #log2FC pvadj
sum(res$padj < 0.05, na.rm=TRUE) #3094

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
#write.csv(resdata, file="./CSV/DESeq2-results-with-normalized.csv")
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

#normalization methods
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
#write.csv(assay(vsd), file = "CSV/DESeq2_vsd_normalized.csv")
#write.csv(assay(rld), file = "CSV/DESeq2_rlog_normalized.csv")

#visualizations
#to account for low count genes log fold changes
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_IL4.treated_vs_control", type="apeglm")
resLFC
#write.csv(as.data.frame(resLFC),file = "./DESeq2_IL4_treated_vs_control_results.unfiltered.csv")
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "IL4 Treated vs Control RNAseq")
dev.copy(png, "./plots/MAplot_healthy_control_vs_control<0.001.png")
dev.off()


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

# genes with padj < 0.1 are colored Red automatically, alpha changes
plotMA(dds, ylim=c(-6,6), alpha= 0.001, main = "Mouse IL4 treated vs Untreated RNAseq")

#plot counts for specific genes
plotCounts(dds, gene=which.min(res$padj), intgroup="condition") #gene with smallest p-val = ACTB
plotCounts(dds, gene=which.max(res$padj), intgroup="condition") #gene with largest p-val = RP11-563J2.2

plotCounts(dds, gene="Gkn3", intgroup="condition")
plotCounts(dds, gene="Tff2", intgroup="condition")
plotCounts(dds, gene="Muc6", intgroup="condition")
plotCounts(dds, gene="Stmn1", intgroup="condition")
plotCounts(dds, gene="Mki67", intgroup="condition")
plotCounts(dds, gene="Spink4", intgroup="condition")
plotCounts(dds, gene="Iqgap3", intgroup="condition")
plotCounts(dds, gene="Birc5", intgroup="condition")
plotCounts(dds, gene="Dmbt1", intgroup="condition")
plotCounts(dds, gene="Cd44", intgroup="condition")
plotCounts(dds, gene="Tff3", intgroup="condition")
plotCounts(dds, gene="Gif", intgroup="condition")
plotCounts(dds, gene="Muc5ac", intgroup="condition")
plotCounts(dds, gene="Gkn2", intgroup="condition")
plotCounts(dds, gene="Chil4", intgroup="condition")


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)

library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(0,25,50,100))

#heatmap of countmatrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

library("RColorBrewer")
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
library("EnhancedVolcano")
spemstem <- c("Hmgb2","Ube2c","Pclaf","Birc5","Cks2","Tff1","Gkn2","H2afz","Tubb5","Top2a","Tuba1b","Stmn1","Sptssb","Gm3776", "Nusap1",
              "Ube2s","Cdk1","Ccna2","Dpcr1", "Cenpa", "Pbk","Pttg1", "Mgst3","Smc4","Gkn1","Ran", "Mki67","Cks1b", "Rrm2", "Lgals4","Spc25",
              "Tubb4b","Smc2","Tmpo",'Tk1',"Muc5ac","Ccnb1", "Cdc20", "Ptma","Jpt1","H2afx","Arl6ip1")
genestolabel <- c(spemstem)
EnhancedVolcano(resLFC, lab=rownames(res), x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 1.3, xlim=c(-5,17), ylim=c(0, max(-log10(data[["pvalue"]]), na.rm = TRUE)), selectLab=genestolabel, drawConnectors = FALSE, title = 'IL4 Treated vs Control RNAseq')
EnhancedVolcano(resLFC, lab=rownames(res), x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 2.0, xlim=c(-5,17), ylim=c(0, max(-log10(data[["pvalue"]]), na.rm = TRUE)), drawConnectors = FALSE, title = 'IL4 Treated vs Control RNAseq' )

#volcanoplot of lfc2
EnhancedVolcano(resLFC,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'IL4 Treated vs Control RNAseq',
                pCutoff = 10e-5,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

#Special heatmap
genesofinterest <- c("Gkn3", "Muc6", "Spink4", "Mki67", "Birc5", "Stmn1", "Cd44", "Iqgap3", "Dmbt1", "Muc5ac", "Gif")
spemstem <- c("Hmgb2","Ube2c","Pclaf","Birc5","Cks2","Tff1","Gkn2","H2afz","Tubb5","Top2a","Tuba1b","Stmn1","Sptssb","Gm3776", "Nusap1",
              "Ube2s","Cdk1","Ccna2","Dpcr1", "Cenpa", "Pbk","Pttg1", "Mgst3","Smc4","Gkn1","Ran", "Mki67","Cks1b", "Rrm2", "Lgals4","Spc25",
              "Tubb4b","Smc2","Tmpo",'Tk1',"Muc5ac","Ccnb1", "Cdc20", "Ptma","Jpt1","H2afx","Arl6ip1")


#define genes to label
genestolabel <- c(spemstem)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/IL4_organoids/data_counts/CSV/DESeq2_rlog_normalized.csv", header = T)


#setup matrix for heatmap
values <- data.frame(data[,c(2:7)])
values <- as.matrix(values)
rownames(values) <- data$gene

#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)

#Complex Heatmap: https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
#Complex Heatmap Manual: https://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/complex_heatmap.html
library(ComplexHeatmap)

#scale data by row
subset_scale <- t(scale(t(subset)))
Heatmap(subset, row_title = c("SPEM-STEM Genes"), 
        row_title_rot = 0, cluster_rows=FALSE, row_names_gp = grid::gpar(fontsize = 10), heatmap_legend_param = list(title = ""))


#Special heatmap
genesofinterest <- c("Gkn3", "Muc6", "Spink4", "Mki67", "Birc5", "Stmn1", "Cd44", "Iqgap3", "Dmbt1", "Muc5ac", "Gif", "Chil4")

#define genes to label
genestolabel <- c(genesofinterest)
genestolabel <- c(spemstem)

data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/IL4_organoids/data_counts/CSV/DESeq2_rlog_normalized.csv", header = T)

#setup matrix for heatmap
values <- data.frame(data[,c(2:7)])
values <- as.matrix(values)
rownames(values) <- data$gene

#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)

#generate basic heatmaps for subsets of genes
heatmap(subset, Colv = NA, Rowv = NA, scale="row")

#heatmap.2: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
library("gplots")

heatmap.2(subset, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "row", main = "Heatmap highlighting DEG in SPEM-STEM",
          col = bluered(100), trace = "none", density.info = "none", keysize = 0.9, add.expr = T, cexRow = 1.2)

#volcanoplot
genestolabel <- c(genesofinterest)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/IL4_organoids/data_counts/CSV/DESeq2_IL4_treated_vs_control_results.csv", header = T)

library(EnhancedVolcano)
names(genestolabel)[genestolabel == 'red'] <- 'Genes of Interest'
#Log2 Fold Change cutoff modified to 2.0 (FC=1), and p.value cutoff modified to 0.05

#generate volcano plot labeled by pathway
EnhancedVolcano(data, lab=data$gene, x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 2.0, xlim=c(-5,17), ylim=c(0, max(-log10(data[["pvalue"]]), na.rm = TRUE)), selectLab=genestolabel, drawConnectors = FALSE, title = 'IL4 Treated vs Control RNAseq',
                labSize=3.0 , pointSize = 1.0)
EnhancedVolcano(data, lab=data$gene, x='log2FoldChange', y="pvalue", pCutoff= 0.05, FCcutoff = 2.0, xlim=c(-5,17), ylim=c(0, max(-log10(data[["pvalue"]]), na.rm = TRUE)), drawConnectors = FALSE, title = 'IL4 Treated vs Control RNAseq',
                labSize=3.0 , pointSize = 1.0)

#converting log2FC to FC
data.transform <- data
data.transform <- data.transform[,c(1,3)]
transformed <- 2^(data.transform[,2])
data[,7] <- transformed
names(data)[7] <- "FC"
final <- data[,c(1,2,3,7,4,5,6)]
write.csv(final, file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/IL4_organoids/data_counts/CSV/DESeq2_IL4_treated_vs_control_results_FC.csv")



