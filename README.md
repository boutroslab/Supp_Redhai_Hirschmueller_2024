# Redhai & Hirschmüller _et al._, 2024

This repository contains scripts and files supporting the publication: </br>

[Redhai & Hirschmüller _et al._](https://www.biorxiv.org/content/10.1101/2024.09.08.611891v1#), **Intestinal stem cell proliferation and differentiation depends on the self-repressive zinc finger transcription factor BCL11/Cph** (2024).


## Abstract
The molecular programs that drive proliferation and differentiation of intestinal stem cells (ISCs) are essential for organismal fitness. Notch signalling regulates the binary fate decision of ISCs, favouring enterocyte commitment when Notch activity is high and enteroendocrine cell (EE) fate when activity is low. However, the temporal dynamics of transcription factor expression and its involvement in generating distinct EE lineages across different intestinal regions remain largely uncharacterised. Here, we find that the expression of the C2H2-type zinc-finger transcription factor Chronophage (Cph), homologous to mammalian BCL11, increases early along the ISC-to-EE lineage when Notch is inactivated. We show that the expression of Cph is regulated by the Achaete-Scute Complex (AS-C) gene, scute, which directly binds to multiple sites within the Cph locus to promote its expression. Our genetic and single-cell RNA sequencing experiments demonstrate that Cph maintains the ISC and EE populations across different regions and is required to remodel the transcriptome of ISCs with low Notch activity. By performing in vivo chromatin profiling, we find that Cph activates the expression of key genes involved in proliferation and lineage commitment, while simultaneously repressing its own expression. This inhibitory feedback loop prevents ISCs from undergoing autophagy-dependent cell death and ensures differentiation is faithfully executed. Our findings shed light on the mechanism by which Cph sustains intestinal epithelial homeostasis and could represent a conserved strategy for balancing proliferation and differentiation in other tissues and species.


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



<ins>**Data explorer (ShinyApp)**</ins> 

Created by Nick Hirschmüller.

The code to recreate the app is located in the `ShinyApp` folder. 









