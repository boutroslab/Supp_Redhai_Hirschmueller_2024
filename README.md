# Redhai & Hirschmüller _et al._, 2024

This repository contains scripts and files supporting the publication: </br>

[Redhai & Hirschmüller _et al._](https://www.biorxiv.org/content/10.1101/2024.09.08.611891v1#), **Proliferation and differentiation of intestinal stem cells depends on the zinc finger transcription factor BCL11/Chronophage** (2024).


## Abstract
The  molecular  programs  that  drive  proliferation  and  differentiation  of  intestinal  stem  cells (ISCs) are essential for organismal fitness. Notch signalling regulates the binary fate decision of ISCs, favouring enterocyte commitment when Notch activity is high and enteroendocrine cell (EE) fate when activity is low. However, the gene regulatory mechanisms that underlie this process on an organ scale remain poorly understood. Here, we find that the expression of the  C2H2-type  zinc-finger  transcription  factor  Chronophage  (Cph),  homologous  to mammalian  BCL11,  increases  specifically  along  the  ISC-to-EE  lineage  when  Notch  is inactivated. We show that the expression of Cph is regulated by the Achaete-Scute Complex (AS-C) gene, scute, which directly binds to multiple sites within the Cph locus to promote its expression. Our genetic and single-cell RNA sequencing experiments demonstrate that Cph maintains  the  ISC  and  EE  populations  and  is  necessary  to  remodel  the  transcriptome  of progenitor cells with low Notch activity. By identifying and functionally validating Cph target genes, we uncover a novel role for sugar free frosting (sff) in directing proliferative and lineage commitment steps of ISCs. Our findings provide mechanistic insight into how Cph maintains intestinal epithelial homeostasis, and together with our evolutionary analysis of Cph/BCL11, suggests that it may represent a general mechanism used by adult somatic stem cells in other tissues to balance proliferation and differentiation. 


## Contact
Should you encounter any issues or have any questions please contact [Nick Hirschmüller](mailto:hirschmueller.nick@gmail.com) or [Siamak Redhai](mailto:siamak.redhai@dkfz-heidelberg.de).

Raw sequencing read data along with count matrices and metadata for each sample were deposited at GEO under the Study ID: [GSE276185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276185).


## Repository structure
<ins>**scRNA-seq data**</ins> 

Analyzed by Nick Hirschmüller and Erica Valentini.

The code to reproduce all analyses is located in the `scRNAseq/analyses` folder. If you run all the scripts in order, all necessary output files will be generated to create the figures from the paper (scripts located in the `scRNAseq/figures` folder).


<ins>**NanoDam, DamID, CHIP-seq and RNAseq data**</ins> 

Analyzed by Stefan Peidli.

The code to reproduce all analyses is located in the `NanoDam_DamID_CHIPseq_RNAseq` folder.











