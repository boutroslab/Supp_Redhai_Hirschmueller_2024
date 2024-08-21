# easy function to plot heatmap for set of genes. each column is one cell
marker_heatmap <- function(seurat, markers, celltype, group.by, cap_value = NULL, cluster_rows_ = F, to_ggplot = T) {
    library(ComplexHeatmap)
    library(circlize)
    # plot a heatmap, where each column is one cell! The cells are grouped and split according to the cluster they are from

    ordered_index <- order(seurat@meta.data[[group.by]])
    c_split <- sort(seurat@meta.data[[group.by]])

    if (!all(markers %in% rownames(seurat@assays$RNA))) {
        mising_genes <- paste(markers[!markers %in% rownames(seurat@assays$RNA)], sep = ", ")
        celltype <- celltype[markers %in% rownames(seurat@assays$RNA)]
        markers <- markers[markers %in% rownames(seurat@assays$RNA)]
        print(str_interp("Not all markers occur in the data matrix of the seurat object. Specifically, ${mising_genes} are missing. Removing it and continuing"))
    }

    # get counts and perform gene wise scaling
    cnts_scaled <- as.matrix(seurat@assays$RNA[markers, ordered_index]) %>%
        t() %>%
        scale() %>%
        t()



    if (!is.null(cap_value)) {
        message(str_interp("Zscores > |${cap_value}| are set to ${cap_value} (or -${cap_value})"))
        changed_values <- sum(cnts_scaled > cap_value) + sum(cnts_scaled < -cap_value)
        cnts_scaled[cnts_scaled > cap_value] <- cap_value
        cnts_scaled[cnts_scaled < -cap_value] <- -cap_value
        message(str_interp("This was the case for ${changed_values} values"))
        col_fun <- colorRamp2(c(-cap_value, 0, cap_value), c("blue", "white", "red"))
    } else {
        col_fun <- colorRamp2(
            breaks = c(min(cnts_scaled, na.rm = T), mean(cnts_scaled, na.rm = T), max(cnts_scaled, na.rm = T)),
            colors = c("blue", "white", "red")
        )
    }

    p <- Heatmap(cnts_scaled,
        cluster_rows = cluster_rows_,
        cluster_columns = F,
        show_row_names = T,
        show_column_names = F,
        column_split = c_split,
        row_split = celltype,
        show_row_dend = F,
        name = "Z-score",
        col = col_fun,
        row_gap = unit(4, "mm"),
        border = "black",
        column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        column_title_rot = 90,
        row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        use_raster = F,
        heatmap_legend_param = list(at = c(-cap_value, 0, cap_value))
    )
    if (to_ggplot) {
        p <- ggplotify::as.ggplot(p)
    }
    return(p)
}

# runs the "normal" seurat pipeline
run_seurat_steps <- function(seurat_object, include_leiden = F, include_norm = T, include_tsne = F, include_louvain = F) {
    if (include_norm) {
        print("Running normalization and variable feature extraction")
        seurat_object <- seurat_object %>%
            NormalizeData(., normalization.method = "LogNormalize", scale.factor = 10000) %>%
            FindVariableFeatures(., selection.method = "vst", nfeatures = 3000)
    }
    print(seurat_object)

    print("Running PCA and UMAP")
    seurat_object <- seurat_object %>%
        ScaleData(.) %>%
        RunPCA(., npcs = ifelse(ncol(seurat_object) < 50, ncol(seurat_object) - 1, 50), verbose = F) %>%
        RunUMAP(., reduction = "pca", dims = 1:20, verbose = F, n.neighbors = ifelse(ncol(seurat_object) < 30, ncol(seurat_object) - 1, 30))

    if (include_tsne) {
        print("Running Tsne")
        seurat_object <- seurat_object %>%
            RunTSNE(., reduction = "pca", dims = 1:20, verbose = F)
    }


    if (include_leiden) {
        print("Running neighborhood detection")
        seurat_object <- seurat_object %>%
            FindNeighbors(., dims = 1:20, k.param = 10, verbose = F) %>%
            FindClusters(.,
                algorithm = 4,
                resolution = c(
                    0.1, 0.3, 0.5,
                    0.6, 0.7
                ),
                method = "igraph"
            )
    }
    if (include_louvain) {
        print("Running neighborhood detection")
        seurat_object <- seurat_object %>%
            FindNeighbors(., dims = 1:20, k.param = 10, verbose = F) %>%
            FindClusters(.,
                resolution = c(
                    0.1, 0.3, 0.5,
                    0.6, 0.7, 1, 1.3, 1.5
                )
            )
    }

    return(seurat_object)
}


run_integration_steps <- function(seurat_object) {
    seurat_object_split <- SplitObject(DietSeurat(seurat_object, assays = "RNA"), split.by = "orig.ident")
    seurat_object_split <- lapply(setNames(names(seurat_object_split), names(seurat_object_split)), function(x) {
        print(x)
        tmp <- seurat_object_split[[x]]
        tmp <- run_seurat_steps(tmp)
    })
    features <- SelectIntegrationFeatures(object.list = seurat_object_split)
    min_size <- min(sapply(seurat_object_split, ncol))
    anchors <- FindIntegrationAnchors(
        object.list = seurat_object_split,
        anchor.features = features,
        dims = 1:ifelse(min_size < 30, min_size - 1, 30),
        k.score = ifelse(min_size < 30, min_size - 1, 30),
        k.filter = ifelse(min_size < 200, 50, 200)
    )
    seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:ifelse(min_size < 20, min_size - 1, 20), k.weight = ifelse(min_size < 100, min_size - 1, 100))
    DefaultAssay(seurat_integrated) <- "integrated"
    seurat_integrated <- run_seurat_steps(seurat_integrated, include_norm = F, include_louvain = T, include_tsne = F)
    return(seurat_integrated)
}



# size=0.05,
# shape=19,
# stroke=0.4
new_dimplot <- function(seurat, feature, dim_reduction = "umap", sort = F, ...) {
    embedding_df <- Embeddings(seurat, dim_reduction)
    data <- FetchData(seurat, feature)
    feature_name <- colnames(data)[1]
    if (sort) {
        data <- data %>% arrange(!!sym(feature_name))
    }

    if (n_distinct(data[[feature_name]]) > 20) {
        print("Assuming the feature of interest is numeric")
        data[[feature_name]] <- as.numeric(data[[feature]])
    }
    data_ <<- data
    data <- left_join(
        data %>% data.frame() %>% rownames_to_column("tmp_id"),
        embedding_df %>% data.frame() %>% rownames_to_column("tmp_id"),
        by = "tmp_id"
    )
    embedding_colnames <- colnames(data)[(ncol(data) - 1):ncol(data)]

    p <- ggplot(data, aes(!!sym(embedding_colnames[[1]]),
        !!sym(embedding_colnames[[2]]),
        color = !!sym(feature_name)
    )) +
        geom_point(...)
    if (n_distinct(data[[feature_name]]) < 20) {
        p <- p +
            guides(color = guide_legend(override.aes = list(size = 1.3)))
    }
    return(p)
}
