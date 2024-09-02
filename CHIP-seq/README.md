Scute binding CHIP-seq data.

- paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5430544/
- data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84283
- figures example: Fig3 in https://elifesciences.org/articles/69937

"The ChIP–seq data were 
- mapped to the Drosophila melanogaster reference genome (BDGP6) using the Bowtie software package (version 1.1.2), with three mismatches allowed. 
- Uniq mapping reads were permitted in the next analysis. 
- We used the MACS software package (version 1.4.2) to call significantly-enriched peaks with a qvalue cutoff of 0.05 relative to the input control samples. 
- To generate the bigwig files, we normalized the tag counts in each bin according to the total number of reads. 
- Input reads were processed in the same way, and their normalized signal intensity values were subtracted from the ChIP-Seq tracks."
