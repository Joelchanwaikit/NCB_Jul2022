
# Count Matrix from company. Metadata generated from console 
raw_counts <- read.csv("C:\\Users\\joelc\\OneDrive\\Desktop\\Count_matrix.csv")
new_metadata <- my_metadata
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]

# Create DESeq Object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(raw_counts),colData = new_metadata,design = ~genotype)

# Normalisation of Count Data 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds,normalized = TRUE)

# Unsupervised Clustering Analyses (QC)
vsd <- vst(dds,blind=TRUE)

# Extract VST Matrix
vsd_mat <- assay(vsd)

# Compute Pairwise correlation function 
vsd_cor <- cor(vsd_mat)

# Heatmap
library(pheatmap)

pheatmap(vsd_cor)

# Principal Component Analysis 

plotPCA(vsd, intgroup = "genotype")

# DE Analysis 

dds <- DESeq(dds)

# Dispersion Estimates 
plotDispEsts(dds)


# Add Log2 Change threshold
library("ashr")
dds_res <- results(dds,contrast=c('genotype','KO','WT') ,alpha = 0.05, lfcThreshold = 0.32)
dds_res <- lfcShrink(dds,contrast=c('genotype','KO','WT') ,res = dds_res,type = "ashr")
plotMA(dds_res)

# Extract Genes 
dds_res_all <- data.frame(dds_res)
dds_res_sig <- subset(dds_res_all,padj<0.05)

# Subset Normalised counts to sig genes 
sig_norm_counts <- normalized_counts[rownames(dds_res_sig),]

# Heatmap again 
library(RColorBrewer)
heat_colors <- brewer.pal(6,"YlOrRd")
pheatmap(sig_norm_counts,color = heat_colors,cluster_rows = T,show_rownames = F, scale = "row")

# Volcano Plot
library(magrittr)
library (tidyverse)
dds_res_all <- dds_res_all %>% rownames_to_column(var="Genes") %>% mutate(threshold = padj <0.05)

ggplot(dds_res_all)+geom_point(aes(x=log2FoldChange,y=-log10(padj),color = threshold))+xlab('Log2 fold change')+ylab("-log10 adjusted p-value")

write.csv(sig_norm_counts,'C:\\Users\\joelc\\OneDrive\\Desktop\\DEResultsFinal.csv')
