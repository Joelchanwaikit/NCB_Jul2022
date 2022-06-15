library(clusterProfiler)
library(enrichplot)
library(ggplot2)

df <- read.csv("C:\\Users\\joelc\\OneDrive\\Desktop\\DEResultsforGSEA.csv", header=TRUE)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$Genes
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$Genes %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$X = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$X

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

emapplot(kk2)

cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

library("ggnewscale")
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

library(enrichplot)
barplot(kk2, showCategory=20) 

