library(DiffBind)
library(tidyverse)

samples <- read.csv('../resources/samplesheet_metadata.csv')
dbObj <- dba(sampleSheet=samples)

# The next step is to take the alignment files and compute count information for each of the peaks/regions in the consensus set. 
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

# To see how well the samples cluster with one another, we can draw a PCA plot using all consensus sites.
pdf('../results/diffbind/DiffBind_pca.pdf')
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

# We can also plot a correlation heatmap, to evaluate the relationship between samples.
pdf('../results/diffbind/DiffBind_corr_heatmap.pdf')
plot(dbObj)
dev.off()

# tell DiffBind which samples we want to compare to one another
print("Contrasting samples...")
dbObj <- dba.contrast(dbObj, contrast=c("Condition", "Cph_NanoDam", "control"), minMembers = 2)

# Performing the differential enrichment analysis
print("Performing differential enrichment analysis...")
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2, design="~Condition")

print("Differential enrichment analysis results:")
# Show summary
dba.show(dbObj, bContrasts=T)

# MA plots are a useful way to visualize the effect of normalization on data, as well as seeing which of the data points are being identified as differentially bound.
pdf('../results/diffbind/DiffBind_MAplot_deseq.pdf')
dba.plotMA(dbObj, method=DBA_DESEQ2)
dev.off()

# concentrations of each sample groups plotted against each other.
pdf('../results/diffbind/DiffBind_MAplot_bXY.pdf')
dba.plotMA(dbObj, bXY=TRUE)
dev.off()

# If we want to see how the reads are distributed amongst the different classes of differentially bound sites and sample groups, we can use a boxplot
pdf('../results/diffbind/DiffBind_boxplot.pdf')
pvals <- dba.plotBox(dbObj)
dev.off()

# Extracting results
res <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(res)
write.table(out, file="../results/diffbind/diffbind_results.tsv", sep="\t", quote=F, row.names=F)