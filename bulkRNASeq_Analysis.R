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

setwd("/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia2")
directory<-"/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia2"

samplenames<-c(list.files(pattern="*.counts.txt"))
print(samplenames)
samplenames<-print(gsub("[:outputcnsx::_::.09o:]", "", samplenames))

control<-grep("control",list.files(directory), value=T)
hcontrol<-grep("healthy",list.files(directory), value=T)
stim<-grep("stim",list.files(directory), value=T)
tmedia<-grep("tmedia",list.files(directory), value=T)
unstim<-grep("unactivated",list.files(directory), value=T)


condition<- c("control", "control", "hcontrol", "hcontrol", "Stim", "Stim", "Tcell_media", "Tcell_media", "Unstim", "Unstim")

Table<-data.frame(sampleNames=samplenames, filename=c(control, hcontrol, stim, tmedia, unstim), condition=condition)
print(Table)
Table$condition <- factor(Table$condition)

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = Table,
                                        directory = directory,
                                        design = ~condition)
ddsHTseq


colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c("control", "hcontrol", "Stim", "Tcell_media", "Unstim"))

dds <- DESeq(ddsHTseq)
res <- results(dds)
res <- res[order(res$padj),] #order most significant to least
head(res)

##write.csv(res, file="avgexp.csv") #log2FC pvadj
sum(res$padj < 0.05, na.rm=TRUE) #3094
resSigind = res[ which(res$padj < 0.05 & res$log2FoldChange > 0), ]
resSigrep = res[ which(res$padj < 0.05 & res$log2FoldChange < 0), ]
resSig = rbind(resSigind, resSigrep)

rownames(resSigind)

resdata <- merge(as.data.frame(resSig), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
#write.csv(resdata, file="./CSV/DESeq2-results-with-normalized.csv")
#write.table(as.data.frame(counts(dds),normalized=T), file = 'CSV/DESeq2_normalized_counts.txt', sep = '\t')
mcols(res, use.names = T)
##write.csv(as.data.frame(mcols(res, use.name = T)),file = "CSV/DESeq2-test-conditions.csv")

#minimize outliers
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.001,
             cleaned = results(ddsClean)$padj < 0.001)
addmargins(tab)
##write.csv(addmargins(tab), file = 'CSV/DESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
##write.csv(as.data.frame(resClean),file = 'CSV/DESeq2-replaceoutliers-results.csv')

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
resLFC <- lfcShrink(dds, coef="condition_hcontrol_vs_control", type="apeglm")
resLFC
#write.csv(as.data.frame(resLFC),file = "./DESeq2_hcontrol_vs_control_results.csv")
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "Healthy Control vs Control RNAseq")
#dev.copy(png, "./plots/MAplot_healthy_control_vs_control<0.001.png")
#dev.off()

resLFC2 <- lfcShrink(dds, coef="condition_Stim_vs_control", type="apeglm")
#write.csv(as.data.frame(resLFC2),file = './DESeq2_Stim_vs_Control_results.csv')
plotMA(resLFC2, ylim=c(-6,6), alpha= 0.001, main = "Stimulated vs Control RNAseq")
#dev.copy(png, "./plots/MAplot_Stimulated_vs_control<0.001.png")
#dev.off()

resLFC3 <- lfcShrink(dds, coef="condition_Tcell_media_vs_control", type="apeglm")
#write.csv(as.data.frame(resLFC3),file = './DESeq2_Tcell_media_vs_control_results.csv')
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "T cell media vs Control RNAseq")
#dev.copy(png, "./plots/MAplot_Tcellmedia_vs_control<0.001.png")
#dev.off()

resLFC4 <- lfcShrink(dds, coef="condition_Unstim_vs_control", type="apeglm")
#write.csv(as.data.frame(resLFC4),file = './DESeq2_Unstim_vs_control_results.csv')
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "Unstimulated vs Control RNAseq")
#dev.copy(png, "./plots/MAplot_unstimulated_vs_control<0.001.png")
#dev.off()
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
plotMA(dds, ylim=c(-6,6), alpha= 0.001, main = "Mouse IL13 treated vs Untreated RNAseq")

#plot counts for specific genes
plotCounts(dds, gene=which.min(res$padj), intgroup="condition") #gene with smallest p-val = ACTB
plotCounts(dds, gene=which.max(res$padj), intgroup="condition") #gene with largest p-val = RP11-563J2.2

plotCounts(dds, gene="RP11-563J2.2", intgroup="condition")
plotCounts(dds, gene="ACTB", intgroup="condition")
plotCounts(dds, gene="TNF", intgroup="condition")
plotCounts(dds, gene="IL6", intgroup="condition")
plotCounts(dds, gene="CXCR4", intgroup="condition")
plotCounts(dds, gene="CXCL12", intgroup="condition")
plotCounts(dds, gene="CD80", intgroup="condition")
plotCounts(dds, gene="CCL2", intgroup="condition")
plotCounts(dds, gene="CD86", intgroup="condition")
plotCounts(dds, gene="CXCL9", intgroup="condition")
plotCounts(dds, gene="SLC2A3", intgroup="condition")
plotCounts(dds, gene="CDCP1", intgroup="condition")


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
threshold <- c(res$log2FoldChange >2, res$padj<0.05)
length(which(threshold)) #5510
res$threshold <- threshold

library("EnhancedVolcano")
#volcanoplot of lfc2
EnhancedVolcano(resLFC2,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Stimulated vs Control RNAseq',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


resSig = c(rownames(resSigind),rownames(resSigrep))

normalized = rlog(Table, blind=FALSE)
rownames(normalized)=rownames(Table)

sigGenes.normalized = normalized[resSig,]

#Special heatmap
hemostatic <- c("P2Y12", "CXC3R1", "TMEM119", "HEXB", "GPR34", "TGFB1", "SALL1")
activation <- c("TREM2", "CD68", "CD80", "CXCR3", "CD40", "TLR4", "HLA-DQA1", "HLA-DQB1")
pro_inflammatory <- c("IL1B", "IL1RL1", "1L1RAP", "IL4", "IL13", "IL4R", "CCL2", "CCL5", "IL6", "CD86","CXCL8", 
                      "CXCL9", "CXCL10", "TNFA", "INOS", "CD32", "CD32A", "CD32B", "MMP9", "CCR4", "CCR5", "CCR3", 
                      "CCR1", "TNFRSF1A", "TNFRSF1B", "NFKBIA", "PER3", "FOS", "JUN", "HLA-DRB1")
repair <- c("MRC1", "TYROBP", "IL10", "CCL22", "IL13", "TMEM2", "CD163", "CD33")

#define genes to label
genestolabel <- c(hemostatic, activation, pro_inflammatory, repair)
data <- read.csv(file = "/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia2/CSV/DESeq2_rlog_normalized_heatmap.csv", header = T)


#setup matrix for heatmap
values <- data.frame(data[,c(2:6)])
values <- as.matrix(values)
rownames(values) <- data$gene

#select genes to plot
subset <- subset(values, rownames(values) %in% genestolabel)

#Complex Heatmap: https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
#Complex Heatmap Manual: https://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/complex_heatmap.html
library(ComplexHeatmap)


splitby <- ifelse(rownames(subset) %in% hemostatic, 'Hemostatic',
                  ifelse(rownames(subset) %in% activation, 'Activation',
                         ifelse(rownames(subset) %in% pro_inflammatory, 'Pro_inflammatory',
                                ifelse(rownames(subset) %in% repair, 'Repair',
                                       'na'))))

#scale data by row
subset_scale <- t(scale(t(subset)))
Heatmap(subset_scale, split = splitby, row_title = c("Hemostatic", "Activation", "Pro_inflammatory", 'Repair'), 
        row_title_rot = 0, cluster_rows=FALSE, row_names_gp = grid::gpar(fontsize = 10), heatmap_legend_param = list(title = ""))


