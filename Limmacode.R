# Generate the 3 required input files 

dfx <- read.csv(file="C:\\Users\\joelc\\OneDrive\\Desktop\\expressionmatrix.csv",row.names =1)
x <- as.matrix(dfx)
f <- read.csv("C:\\Users\\joelc\\OneDrive\\Desktop\\featuredata.csv",row.names =1)
p <- read.csv("C:\\Users\\joelc\\OneDrive\\Desktop\\microarrayphenotype.csv",row.names =1)



# Initialize Library 
library(Biobase)
library(limma)

# Create ExpressionSet object
eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))


# Create new ExpressionSet to store normalized data
eset_norm <- eset

# View the distribution of the raw data
plotDensities(eset_norm, legend = FALSE)

# Log tranform
exprs(eset_norm) <- log(exprs(eset_norm))
plotDensities(eset_norm, legend = FALSE)

# Quantile normalize
exprs(eset_norm) <- normalizeBetweenArrays(exprs(eset_norm))
plotDensities(eset_norm, legend = FALSE)

# Create new ExpressionSet to store filtered data
eset <- eset_norm

# View the normalized gene expression levels
plotDensities(eset, legend = FALSE); abline(v = 5)

# Determine the genes with mean expression level greater than 5
keep <- rowMeans(exprs(eset)) > 5
sum(keep)

# Filter the genes
eset <- eset[keep,]
plotDensities(eset, legend = FALSE)

# Plot principal components labeled by Genotype
plotMDS(eset, labels = pData(eset)[, "Genotype"], gene.selection = "common")

# Plot principal components labeled by batch
plotMDS(eset, labels = pData(eset)[, "Batch"], gene.selection = "common")

# Remove the batch effect
exprs(eset) <- removeBatchEffect(eset, batch = pData(eset)[,'Batch'])

# Plot principal components labeled by Genotype
plotMDS(eset, labels = pData(eset)[, "Genotype"], gene.selection = "common")

# Plot principal components labeled by batch
plotMDS(eset, labels = pData(eset)[, "Batch"], gene.selection = "common")

# Create design matrix for study
design <- model.matrix(~Genotype, data = pData(eset))

# Count the number of samples modeled by each coefficient
colSums(design)

# Fit the model
fit <- lmFit(eset, design)

# Calculate the t-statistics
fit <- eBayes(fit)

# Summarize results
results <- decideTests(fit[, "GenotypeWT"])
summary(results)

# Obtain the summary statistics for every gene
stats <- topTable(fit, number = nrow(fit), sort.by = "none")

# Plot a histogram of the p-values
hist(stats[, "P.Value"])

# Create a volcano plot. Highlight the top 5 genes (Not working)
volcanoplot(fit)

write.csv(stats,"C:\\Users\\joelc\\OneDrive\\Desktop\\limmaresults.csv")