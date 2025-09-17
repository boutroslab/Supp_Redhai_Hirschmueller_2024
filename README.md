# Redhai & Hirschmüller _et al._, 2024

This repository contains scripts and files supporting the publication: </br>

[Redhai & Hirschmüller _et al._](https://www.biorxiv.org/content/10.1101/2024.09.08.611891v1#), **Intestinal stem cell proliferation and differentiation depends on the self-repressive zinc finger transcription factor BCL11/Cph** (2024).


## Abstract
Intestinal stem cells (ISCs) continuously renew the gut epithelium by generating regionally specialized cell types, but how they commit to distinct lineages remains poorly defined. Here, we identify a self-limiting transcriptional program, mediated by the zinc-finger transcription factor Chronophage (Cph), that drives ISC maintenance and differentiation into enteroendocrine (EE) cells in different regions of the Drosophila midgut. We show that Cph expression is transiently induced by the proneural factor scute at the onset of ISC-to-EE fate specification. Genetic and single-cell transcriptomics approaches revealed that Cph is required to intrinsically remodel the transcriptome of ISCs and sustains normal lifespan. Genome-wide chromatin profiling demonstrated that Cph directly binds and regulates the expression of key proliferation and differentiation genes, while simultaneously repressing its own expression. This autoinhibitory feedback safeguards ISCs from undergoing autophagy and cell death, thus ensuring proliferation and differentiation are faithfully executed. Our findings highlight a key mechanism that balances ISC maintenance and lineage commitment. 

## Contact
Should you encounter any issues or have any questions please contact [Nick Hirschmüller](mailto:hirschmueller.nick@gmail.com) or [Siamak Redhai](mailto:siamak.redhai@dkfz-heidelberg.de).

Raw sequencing read data along with count matrices and metadata for each sample were deposited at GEO under the Study ID: [GSE276185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276185).


## Repository structure
<ins>**scRNA-seq data**</ins> 

Analyzed by Nick Hirschmüller and Erica Valentini.

The code to reproduce all analyses is located in the `scRNAseq/analyses` folder. If you run all the scripts in order, all necessary output files will be generated to create the figures from the paper (scripts located in the `scRNAseq/figures` folder).


<ins>**NanoDam, DamID, CHIP-seq and bulk RNAseq data**</ins> 

Analyzed by Stefan Peidli.

The code to reproduce all analyses is located in the `NanoDam_DamID_CHIPseq_RNAseq` folder.



<ins>**Data explorer (ShinyApp)**</ins> 

Created by Nick Hirschmüller.

The code to recreate the app is located in the `ShinyApp` folder. 









