## loading packages
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(clusterProfiler)
})

# Handle file names
bed_file <- snakemake@input[['bed_file']]
target <- snakemake@output[['annotated_file']]

# Set database for genome annotation
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Load contrasted peaks from MACS2
peak <- readPeakFile(bed_file)
seqlevelsStyle(peak) <- "UCSC"  # e.g. "chr2R" -> "2R"

# Define TSS regions (takes ca 30sec)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

# Peak annotation
peakAnno <- annotatePeak(peak, 
    tssRegion=c(-3000, 3000),
    TxDb=txdb)

# Save the results
peakAnno.df <- as.data.frame(peakAnno)
print(peakAnno.df[1:5,])
print(target)
# write to csv
write.csv(peakAnno.df, file = target)