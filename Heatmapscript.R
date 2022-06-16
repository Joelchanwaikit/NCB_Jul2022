library(pheatmap)

my_data <- read.csv("C:\\Users\\joelc\\OneDrive\\Desktop\\BreastTCGAHeatmap.csv")



my_matrix_top25 <- as.matrix(my_data[,2:26])
rownames(my_matrix_top25) = my_data$Symbol
my_df_top25 <- scale(my_matrix_top25)
my_tdf_top25 <- t(my_df_top25)

library(RColorBrewer)
heat_colors <- brewer.pal(8,"RdYlBu")
pheatmap(my_tdf_top25,main = "TCGA_BRCA (Top & Bottom 25% Z Scores)",cluster_cols=FALSE,
         cluster_rows = FALSE, show_colnames = F ,gaps_row = c(13,24),gaps_col = c(274,548),color = heat_colors)



