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

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
#write.csv(resdata, file="./DESeq2-results-with-normalized.csv")
##write.table(as.data.frame(counts(dds),normalized=T), file = 'CSV/DESeq2_normalized_counts.txt', sep = '\t')
mcols(res, use.names = T)
##write.csv(as.data.frame(mcols(res, use.name = T)),file = "CSV/DESeq2-test-conditions.csv")

#minimize outliers
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.001,
             cleaned = results(ddsClean)$padj < 0.001)
addmargins(tab)
##write.csv(addmargins(tab), file = 'CSV/ESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
##write.csv(as.data.frame(resClean),file = 'CSV/DESeq2-replaceoutliers-results.csv')

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

vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE) for smaller datasets <30
head(assay(vsd), 3)

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

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

plotPCA(vsd, intgroup="condition")
#plotPCA(rld, intgroup="condition")
plotPCA(ntd, intgroup="condition")

#for comparing against stimulated

setwd("/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia2")
directory<-"/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia2"

samplenames<-c(list.files(pattern="*.counts.txt"))
print(samplenames)
samplenames<-print(gsub("[:outputcnsx::_::.09o:]", "", samplenames))

new_names <- c("im-1", "im-2", "media-1", "media-2", "aivaed-1", "aivaed-2", "rl-1", "rl-2", "healhy-1", "healhy-2")
               
stim<-grep("stim",list.files(directory), value=T)
tmedia<-grep("tmedia",list.files(directory), value=T)
unstim<-grep("unactivated",list.files(directory), value=T)
control<-grep("control",list.files(directory), value=T)
hcontrol<-grep("healthy",list.files(directory), value=T)


condition<- c("Stim", "Stim", "Tcell_media", "Tcell_media", "Unstim", "Unstim", "control", "control", "hcontrol", "hcontrol")

Table<-data.frame(sampleNames=new_names, filename=c(stim, tmedia, unstim, control, hcontrol), condition=condition)
print(Table)
Table$condition <- factor(Table$condition)

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = Table,
                                        directory = directory,
                                        design = ~condition)
ddsHTseq


colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c("Stim", "Tcell_media", "Unstim", "control", "hcontrol"))

dds <- DESeq(ddsHTseq)
res <- results(dds)
res <- res[order(res$padj),] #order most significant to least
head(res)

##write.csv(res, file="avgexp.csv") #log2FC pvadj
sum(res$padj < 0.001, na.rm=TRUE) #3613

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
##write.csv(resdata, file="CSV/DESeq2-results-with-normalized.csv")
##write.table(as.data.frame(counts(dds),normalized=T), file = 'CSV/DESeq2_normalized_counts.txt', sep = '\t')
mcols(res, use.names = T)
##write.csv(as.data.frame(mcols(res, use.name = T)),file = "CSV/DESeq2-test-conditions.csv")

#minimize outliers
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.001,
             cleaned = results(ddsClean)$padj < 0.001)
addmargins(tab)
##write.csv(addmargins(tab), file = 'CSV/ESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
sum(res$padj < 0.001, na.rm=TRUE) #3613
sum(res$padj < 0.05, na.rm=TRUE) #5510
#write.csv(as.data.frame(resClean),file = 'CSV/DESeq2-replaceoutliers-results.csv')

#visualizations
#to account for low count genes log fold changes
resultsNames(dds)
resLFC5 <- lfcShrink(dds, coef="condition_control_vs_Stim", type="apeglm")
resLFC
#write.csv(as.data.frame(resLFC5),file = './DESeq2_control_vs_stim_results.csv')
plotMA(resLFC5, ylim=c(-6,6), alpha= 0.001, main = "Control vs Stim RNAseq")
#dev.copy(png, "./plots/MAplot_control_vs_stim_p<0.001.png")
#dev.off()

resLFC6 <- lfcShrink(dds, coef="condition_hcontrol_vs_Stim", type="apeglm")
#write.csv(as.data.frame(resLFC6),file = './DESeq2_hcontrol_vs_Stim.csv')
plotMA(resLFC6, ylim=c(-6,6), alpha= 0.001, main = "Healthy Control vs Stim RNAseq")
#dev.copy(png, "./plots/MAplot_hcontrol_vs_Stim_p<0.001.png")
#dev.off()

resLFC7 <- lfcShrink(dds, coef="condition_Tcell_media_vs_Stim", type="apeglm")
#write.csv(as.data.frame(resLFC7),file = './DESeq2_Tcell_media_vs_stim_results.csv')
plotMA(resLFC7, ylim=c(-6,6), alpha= 0.001, main = "T cell media vs Stim RNAseq")
#dev.copy(png, "./plots/MAplot_Tcellmedia_vs_stim_p<0.001.png")
#dev.off()

resLFC8 <- lfcShrink(dds, coef="condition_Unstim_vs_Stim", type="apeglm")
#write.csv(as.data.frame(resLFC8),file = './DESeq2_Unstim_vs_stim_results.csv')
plotMA(resLFC8, ylim=c(-6,6), alpha= 0.001, main = "Unstimulated vs Stimulated RNAseq")
#dev.copy(png, "./plots/MAplot_unstimulated_vs_stim_p<0.001.png")
#dev.off()

#volcanoplots
threshold <-  res$padj<0.05
length(which(threshold)) #5510
res$threshold <- threshold

library("EnhancedVolcano")
EnhancedVolcano(resLFC5,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Control vs Stimulated',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

EnhancedVolcano(resLFC6,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Healthy Control vs Stimulated',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

EnhancedVolcano(resLFC7,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Tcell media vs Stimulated',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

EnhancedVolcano(resLFC8,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Unstimulated vs Stimulated',
                pCutoff = 10e-6,
                FCcutoff = 1.5,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

control_stim <- read.csv("./DESeq2_control_vs_stim_results.csv", header = TRUE, stringsAsFactors = FALSE)
names(control_stim)[1] <- 'gene'
head(control_stim)
ggplot(control_stim) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text(aes(x = log2FoldChange, y = -log10(padj), label=ifelse(log2FoldChange >2,as.character(gene),'')),hjust=0,vjust=0) +
  ggtitle("Volcano Plot Control vs Stimulated") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

unstim_stim <- read.csv("./DESeq2_Unstim_vs_stim_results.csv", header = TRUE, stringsAsFactors = FALSE)
names(unstim_stim)[1] <- 'gene'
head(unstim_stim)
ggplot(unstim_stim) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text(aes(x = log2FoldChange, y = -log10(padj), label=ifelse(log2FoldChange>2,as.character(gene),'')),hjust=0,vjust=0) +
  ggtitle("Volcano Plot Unstimulated vs Stimulated") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

DEGreport::degPlot(dds = dds, res = res, n=30, xs="sizeFactor", group = "condition") # dds object is output from DESeq2

tcellmedia_stim <- read.csv("./DESeq2_Tcell_media_vs_stim_results.csv", header = TRUE, stringsAsFactors = FALSE)
names(tcellmedia_stim)[1] <- 'gene'
head(tcellmedia_stim)
ggplot(tcellmedia_stim) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text(aes(x = log2FoldChange, y = -log10(padj), label=ifelse(log2FoldChange< -2,as.character(gene),'')),hjust=0,vjust=0) +
  ggtitle("Volcano Plot T cell media vs Stimulated") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

DEGreport::degPlot(dds = dds, res = res, n=30, xs="sizeFactor", group = "condition") # dds object is output from DESeq2

unstim_stim <- read.csv("./DESeq2_Unstim_vs_stim_results.csv", header = TRUE, stringsAsFactors = FALSE)
names(unstim_stim)[1] <- 'gene'
head(unstim_stim)
ggplot(unstim_stim) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text(aes(x = log2FoldChange, y = -log10(padj), label=ifelse(log2FoldChange< -2,as.character(gene),'')),hjust=0,vjust=0) +
  ggtitle("Volcano Plot Unstimulated vs Stimulated") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 



DEGreport::degPlot(dds = dds, res = res, n=30, xs="sizeFactor", group = "condition") # dds object is output from DESeq2


#Reordered for comparing against healthy control

setwd("/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia")
directory<-"/Users/dipaololabmacbook/Desktop/BulkRNASeq21_22/microglia2"

samplenames<-c(list.files(pattern="*.counts.txt"))
print(samplenames)
samplenames<-print(gsub("[:outputcnsx::_::.09o:]", "", samplenames))

hcontrol<-grep("healthy",list.files(directory), value=T)
stim<-grep("stim",list.files(directory), value=T)
control<-grep("control",list.files(directory), value=T)
tmedia<-grep("tmedia",list.files(directory), value=T)
unstim<-grep("unactivated",list.files(directory), value=T)


condition<- c("hcontrol", "hcontrol", "Stim", "Stim", "control", "control", "Tcell_media", "Tcell_media", "Unstim", "Unstim")

Table<-data.frame(sampleNames=samplenames, filename=c(hcontrol, stim, control, tmedia, unstim), condition=condition)
print(Table)
Table$condition <- factor(Table$condition)

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = Table,
                                        directory = directory,
                                        design = ~condition)
ddsHTseq


colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c("hcontrol", "Stim" ,"control", "Tcell_media", "Unstim"))

dds <- DESeq(ddsHTseq)
res <- results(dds)
res <- res[order(res$padj),] #order most significant to least
head(res)

##write.csv(res, file="avgexp.csv") #log2FC pvadj
sum(res$padj < 0.001, na.rm=TRUE) #1787

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
##write.csv(resdata, file="CSV/DESeq2-results-with-normalized.csv")
##write.table(as.data.frame(counts(dds),normalized=T), file = 'CSV/DESeq2_normalized_counts.txt', sep = '\t')
mcols(res, use.names = T)
##write.csv(as.data.frame(mcols(res, use.name = T)),file = "CSV/DESeq2-test-conditions.csv")

#minimize outliers
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.001,
             cleaned = results(ddsClean)$padj < 0.001)
addmargins(tab)
##write.csv(addmargins(tab), file = 'CSV/ESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
##write.csv(as.data.frame(resClean),file = 'CSV/DESeq2-replaceoutliers-results.csv')

#visualizations
#to account for low count genes log fold changes
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_Stim_vs_hcontrol", type="apeglm")
resLFC
#write.csv(as.data.frame(resLFC),file = 'CSV/others_vs_hcontrol/DESeq2_Stim_vs_hcontrolresults.csv')
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "Stim vs Healthy Control RNAseq")
#dev.copy(png, "./plots/MAplot_stim_vs_healthyc_p<0.001.png")
#dev.off()

resLFC <- lfcShrink(dds, coef="condition_control_vs_hcontrol", type="apeglm")
#write.csv(as.data.frame(resLFC),file = 'CSV/others_vs_hcontrol/DESeq2_control_vs_hcontrol.csv')
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "Control vs Healthy Control RNAseq")
#dev.copy(png, "./plots/MAplot_control_vs_hcontrol_p<0.001.png")
#dev.off()

resLFC <- lfcShrink(dds, coef="condition_Tcell_media_vs_hcontrol", type="apeglm")
#write.csv(as.data.frame(resLFC),file = 'CSV/others_vs_hcontrol/DESeq2_Tcell_media_vs_hcontrol_results.csv')
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "T cell media vs Healthy Control RNAseq")
#dev.copy(png, "./plots/MAplot_Tcellmedia_vs_hcontrol_p<0.001.png")
#dev.off()

resLFC <- lfcShrink(dds, coef="condition_Unstim_vs_hcontrol", type="apeglm")
#write.csv(as.data.frame(resLFC),file = 'CSV/others_vs_hcontrol/DESeq2_Unstim_vs_hcontrol_results.csv')
plotMA(resLFC, ylim=c(-6,6), alpha= 0.001, main = "Unstimulated vs Healthy Control RNAseq")
#dev.copy(png, "./plots/MAplot_unstimulated_vs_hcontrol_p<0.001.png")
#dev.off()

#volcanoplots
Stim_hcontrol <- read.csv("./DESeq2_Stim_vs_hcontrolresults.csv", header = TRUE, stringsAsFactors = FALSE)
names(Stim_hcontrol)[1] <- 'gene'
head(Stim_hcontrol)
ggplot(Stim_hcontrol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text(aes(x = log2FoldChange, y = -log10(padj), label=ifelse(log2FoldChange< -2,as.character(gene),'')),hjust=0,vjust=0) +
  ggtitle("Volcano Plot Stimulated vs Healthy Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

DEGreport::degPlot(dds = dds, res = res, n=30, xs="sizeFactor", group = "condition") # dds object is output from DESeq2

control_hcontrol <- read.csv("./DESeq2_control_vs_hcontrol.csv", header = TRUE, stringsAsFactors = FALSE)
names(control_hcontrol)[1] <- 'gene'
head(control_hcontrol)
ggplot(control_hcontrol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text(aes(x = log2FoldChange, y = -log10(padj), label=ifelse(log2FoldChange< -2.3,as.character(gene),'')),hjust=0,vjust=0) +
  ggtitle("Volcano Plot Control vs Healthy Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


