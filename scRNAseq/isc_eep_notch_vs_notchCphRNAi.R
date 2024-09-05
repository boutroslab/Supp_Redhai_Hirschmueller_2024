.libPaths(c(
    "/g/huber/users/hirschmueller/R/4.2.0-foss-2021b",
    "/home/hirschmu/R/x86_64-pc-linux-gnu-library/4.2",
    "/g/easybuild/x86_64/Rocky/8/haswell/software/R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0",
    "/g/easybuild/x86_64/Rocky/8/haswell/software/arrow-R/6.0.0.2-foss-2021b-R-4.2.0",
    "/g/easybuild/x86_64/Rocky/8/haswell/software/R/4.2.0-foss-2021b/lib64/R/library"
))


library(Seurat)
library(tidyverse)
library(here)
library(MAST)
library(SingleCellExperiment)
library(scater)


source(here("plot_theme.R"))
source(here("helper_functions.R"))


seurat <- readRDS(here("output", "Ctrl_Notch_NotchCphRNAi_integrated_scent.rds"))
DefaultAssay(seurat) <- "RNA"
# -> remove irrelevant celltypes
seurat <- seurat[, !seurat$celltype_manual %in% c("unk", "unk2", "MT", "dEC", "daEC")]

# for this analysis we will aggregate multiple celltypes:

# we combine ISC and EEPs to progenitor population
seurat <- seurat[, seurat$high_res_annotation %in% c("ISC", "EEP")]
seurat$tmp_celltype <- "Progenitor_ISC_EEP"


sce <- as.SingleCellExperiment(seurat)

# rerun log2 transformation
sce <- scater::logNormCounts(sce)
sca <- SceToSingleCellAssay(sce, class = "SingleCellAssay")


# How to run must was extracted from here:
# https://www.sc-best-practices.org/conditions/differential_gene_expression.html#id21
# and here https://github.com/kdzimm/PseudoreplicationPaper/blob/c3059a3b361e89bde595f222757d04b89f77eb62/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

MAST_results <- lapply(setNames(unique(sce$tmp_celltype), unique(sce$tmp_celltype)), function(celltype) {
    print(celltype)
    all_conditions <- setdiff(unique(sce$perturbation), "ctrl")
    tmp_sca <- sca[, sca$tmp_celltype == celltype]

    geneFilter <- freq(tmp_sca) > 0.1
    print(str_interp("MAST will be run for ${sum(geneFilter)} genes"))

    tmp_sca <- tmp_sca[geneFilter, ]

    tmp_sca$perturbation <- factor(tmp_sca$perturbation)
    tmp_sca$perturbation <- relevel(tmp_sca$perturbation, "NotchCphRNAi") # we want to compare to NotchCphRNAi not control!

    tmp_sca$batch <- factor(substr(tmp_sca$orig.ident.unique, 1, 4))
    tmp_sca$orig.ident.unique <- factor(tmp_sca$orig.ident.unique)

    # calculate how many genes are "on", i.e. expressed.
    cdr2 <- colSums(assay(tmp_sca) > 0)
    colData(tmp_sca)$ngeneson <- scale(cdr2)


    # Differential expression analysis with random effects.
    formula <- ~ ngeneson + batch + perturbation + (1 | orig.ident.unique)
    zlmCond <- zlm(
        formula = formula,
        sca = tmp_sca,
        exprs_values = "logcounts",
        method = "glmer",
        ebayes = F,
        strictConvergence = F,
        fitArgsD = list(nAGQ = 0)
    )

    per_condition_res <- lapply(setNames(all_conditions, all_conditions), function(condition) {
        if (condition != "NotchRNAi") {
            return()
        }
        print(condition)
        summaryCond <- MAST::summary(zlmCond, doLRT = paste0("perturbation", condition))
        summaryDt <- summaryCond$datatable
        result <- merge(summaryDt[contrast == paste0("perturbation", condition) & component == "H", .(primerid, `Pr(>Chisq)`)], # hurdle P values
            summaryDt[contrast == paste0("perturbation", condition) & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
            by = "primerid"
        ) # logFC coefficients

        # do multiple testing correction
        result[, FDR := p.adjust(`Pr(>Chisq)`, "fdr")]
        result <- data.frame(result) %>%
            arrange(FDR) %>%
            tidyr::drop_na()
        return(result)
    })
    return(per_condition_res)
})
saveRDS(MAST_results, here("output", "differential_expression", "MAST", "MAST_ISC_EEP_NotchRNAi_vs_NotchCphRNAi.rds"))
