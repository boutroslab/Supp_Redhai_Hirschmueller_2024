# Redhai & Hirschmüller _et al._, 2025

This repository contains scripts and files supporting the publication: </br>

[Redhai & Hirschmüller _et al._](https://www.biorxiv.org/content/10.1101/2024.09.08.611891v1#), **Autoinhibitory feedback preserves intestinal stem cell maintenance and fate commitment ** (2025).


## Abstract
Intestinal stem cells (ISCs) continuously renew the gut epithelium by producing specialized cell types, yet the mechanisms that couple ISC renewal with lineage commitment remain poorly characterized. Here, we identify a self-limiting transcriptional program, mediated by the zinc-finger transcription factor *Chronophage* (*Cph*), that promotes both ISC maintenance and differentiation into enteroendocrine (EE) cells in the *Drosophila* midgut. *Cph* expression is transiently induced by the proneural factor *scute* at the onset of ISC-to-EE specification. Genetic and single-cell transcriptomic approaches revealed that *Cph* is required to reprogram ISCs and sustain normal lifespan. Cph binds to genes involved in proliferation and differentiation, and directly represses its own expression. This autoinhibitory feedback safeguards ISCs from accumulating autophagosomes and undergoing cell death, thus preserving ISC function. Our findings uncover a key regulatory mechanism that balances stem cell maintenance and differentiation, highlighting principles relevant to regenerating tissues.

## Contact
Should you encounter any issues or have any questions please contact [Nick Hirschmüller](mailto:hirschmueller.nick@gmail.com) or [Siamak Redhai](mailto:siamak.redhai@dkfz-heidelberg.de).

Raw sequencing read data along with count matrices and metadata for each sample were deposited at GEO under the Study ID: [GSE276185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276185).

## Repository structure
<ins>**scRNA-seq data**</ins> 

Analyzed by Nick Hirschmüller and Erica Valentini.

Processed Seurat objects that can be used to reproduce figures are deposited on Zenodo: https://doi.org/10.5281/zenodo.17716158

The code to reproduce all analyses is located in the `scRNAseq/analyses` folder. If you run all the scripts in order, all necessary output files will be generated to create the figures from the paper (scripts located in the `scRNAseq/figures` folder).


<ins>**NanoDam, DamID, CHIP-seq and bulk RNAseq data**</ins> 

Analyzed by Stefan Peidli.

The code to reproduce all analyses is located in the `NanoDam_DamID_CHIPseq_RNAseq` folder.



<ins>**Data explorer (ShinyApp)**</ins> 

Created by Nick Hirschmüller.

The code to recreate the app is located in the `ShinyApp` folder. 









