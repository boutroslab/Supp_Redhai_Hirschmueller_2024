- paper: https://www.nature.com/articles/s41467-022-34270-0
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211629

Dam-ID sequencing and data analysis
- Raw sequencing reads were mapped to the D. melanogaster genome (BDGP6) using bowtie2 (version 2.2.4). 
- R package iDEAR 0.1.0 was subsequently used for establishing highly reliable profiles of Pros-binding sites in Drosophila EEs. 
- The density of Pros-binding on each genomic region was normalized to the total number of mapped reads. 
- BigWig files were generated for visualization using the Homer package. 
- Annotated oDam peaks (oDam-Pros versus oDam) with log2FoldChange more than 1 and adjusted p value <0.01 were filtered as putative Pros-binding targets.

Idea: Check if Pros is binding to Cph.