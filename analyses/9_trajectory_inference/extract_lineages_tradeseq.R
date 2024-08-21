# .libPaths(c("/g/easybuild/x86_64/CentOS/7/haswell/software/R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0",.libPaths()))
library(slingshot)
library(scran)
library(parallel)
library(tradeSeq)
library(tidyverse)
library(Seurat)

library(pbapply)
library(pbmcapply)
library(here)
gam_fit_6 <- readRDS(here("output", "tradeSeq_fit_condiments.rds"))

# the vast majority of this code was just copied straight from the source code of
# `tradeSeq::plotSmoothers`.
extract_plot_data <- function(gene) {
    counts <- assays(gam_fit_6)$counts
    models <- gam_fit_6
    nPoints <- 150


    id <- which(names(models) %in% gene)
    dm <- colData(models)$tradeSeq$dm # design matrix
    y <- unname(counts[names(models), ][id, ])
    X <- colData(models)$tradeSeq$X # linear predictor
    slingshotColData <- colData(models)$crv
    pseudotime <- slingshotColData[, grep(
        x = colnames(slingshotColData),
        pattern = "pseudotime"
    )]
    nLineages <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
    beta <- rowData(models)$tradeSeq[id, ]$beta[[1]]

    conditions <- colData(models)$tradeSeq$conditions
    nConditions <- nlevels(conditions)

    # Construct time variable based on cell assignments.
    lcol <- timeAll <- rep(0, nrow(dm))
    # Loop over all curves, i.e. all lineqges and conditions
    for (jj in seq_len(nLineages)) {
        for (kk in seq_len(nConditions)) {
            for (ii in seq_len(nrow(dm))) {
                if (dm[ii, paste0("l", jj, "_", kk)] == 1) {
                    timeAll[ii] <- dm[ii, paste0("t", jj)]
                    lcol[ii] <- paste0("Lineage ", jj, "_", levels(conditions)[kk])
                } else {
                    next
                }
            }
        }
    }

    df <- data.frame(
        "time" = timeAll,
        "gene_count" = y,
        "lineage" = as.character(lcol)
    )


    # WE WILL JUST USE LOESS REGRESSION 
    # combs <- paste0("Lineage ", seq_len(nLineages), "_")
    # combs_list <- lapply(combs, function(lin) paste0(lin, levels(conditions)))
    # combs <- do.call("c", combs_list)
    # df$lineage <- factor(df$lineage, levels = combs)
    # 
    # df_curves <- mclapply(seq_len(nLineages), function(jj) {
    #     tmp <- lapply(seq_len(nConditions), function(kk) {
    #         df <- tradeSeq:::.getPredictRangeDf(dm,
    #             lineageId = jj, conditionId = kk,
    #             nPoints = nPoints
    #         )
    #         Xdf <- tradeSeq:::predictGAM(
    #             lpmatrix = X,
    #             df = df,
    #             pseudotime = pseudotime,
    #             conditions = conditions
    #         )
    #         yhat <- c(exp(t(Xdf %*% t(beta)) + df$offset))
    # 
    #         df_curves <- data.frame(
    #             "time" = df[, paste0("t", jj)],
    #             "gene_count" = yhat,
    #             "lineage" = as.character(combs_list[[jj]][[kk]])
    #         )
    # 
    #         return(df_curves)
    #     })
    #     names(tmp) <- c("ctrl", "notch")
    #     return(tmp)
    # }, mc.cores = 1)
    # names(df_curves) <- c("Lineage1", "Lineage2", "Lineage3", "Lineage4")
    return(list(
        "df_curves" = NULL,
        "df_points" = df
    ))
}

# takes ~30 mins with 20 cores
all_genes <- rownames(gam_fit_6)
res <- pbmclapply(setNames(all_genes, all_genes), extract_plot_data, mc.cores = 20)
saveRDS(res, here("output", "plotSmoothers_data.rds"))




