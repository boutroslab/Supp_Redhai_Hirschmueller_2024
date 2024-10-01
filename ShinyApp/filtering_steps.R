# here, we created a very reduced seurat object to save space and decrease load time
library(here)
library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
###################
# For the NotchKO #
###################
seurat <- readRDS(here("data", "ctrl_notch_integrated_scent.rds"))


# update the celltypes
seurat@meta.data$celltype_manual <- ifelse(
    seurat$high_res_annotation == "EEP",
    "EEP",
    seurat$celltype_manual
)
seurat$celltype_manual <- ifelse(seurat$celltype_manual == "unk2", "unk", seurat$celltype_manual)
seurat$celltype_manual <- factor(seurat$celltype_manual,
    levels = c("ISC", "EB", "EEP", "dEC", "daEC", "aEC", "mEC", "Copper", "LFC", "pEC", "EE", "MT", "unk")
)

# remove unnecessary data.
reductions <- c("pca", "phate", "umap")
dim_data <- data.frame(CellID = rownames(seurat@meta.data))

for (reduction in reductions) {
    dim_data <- cbind(dim_data, Embeddings(seurat, reduction)[, 1:2])
}
dim_data <- dim_data %>% select(-CellID)

seurat_small <- DietSeurat(seurat,
    layers = "data",
    assays = "RNA",
    misc = FALSE
)


# remove all genes that are expressed in less than 15 cells
seurat_small <- seurat_small[rowSums(seurat_small@assays$RNA@data > 0) >= 15, ]


# remove all counts, we will only work with normalized counts.
# this saves a lot of space
empty_matrix <- sparseMatrix(dims = c(nrow(seurat_small),ncol(seurat_small)), i={}, j={})
empty_matrix <- as(empty_matrix, "dgCMatrix")
dimnames(empty_matrix) <- dimnames(seurat_small)
seurat_small <- SetAssayData(seurat_small, slot = "counts", new.data = empty_matrix)

# only meta data that is relevant
seurat_small@meta.data <- seurat_small@meta.data[, c("celltype_manual", "perturbation")]

all(Cells(seurat_small) == rownames(dim_data))

# add the dimension data to the meta data
seurat_small@meta.data <- cbind(seurat_small@meta.data, dim_data)

# change the gene name of cph
rownames(seurat_small@assays$RNA@data) <-
    ifelse(rownames(seurat_small@assays$RNA@data) == "CG9650", "Cph", rownames(seurat_small@assays$RNA@data))
rownames(seurat_small@assays$RNA@counts) <-
    ifelse(rownames(seurat_small@assays$RNA@counts) == "CG9650", "Cph", rownames(seurat_small@assays$RNA@counts))

saveRDS(seurat_small, here("data", "filtered_seurat_object.rds"))

############################
# For the RNAi Experiments #
############################
seurat <- readRDS(here("data", "RNAi_experiments_integrated.rds"))
# remove the unkonw cells
seurat <- seurat[, seurat$celltype_manual != "unk"]
seurat@meta.data$celltype_manual <- case_when(
    seurat$high_res_annotation == "EEP" ~ "EEP",
    T ~ seurat$celltype_manual
)


seurat$celltype_manual <- factor(seurat$celltype_manual,
    levels =
        c("ISC", "EEP", "EB", "dEC", "daEC", "aEC", "mEC", "Copper", "LFC", "pEC", "EE", "MT")
)

# remove unnecessary data.
reductions <- c("pca", "phate", "umap")
dim_data <- data.frame(CellID = rownames(seurat@meta.data))

for (reduction in reductions) {
    dim_data <- cbind(dim_data, Embeddings(seurat, reduction)[, 1:2])
}
dim_data <- dim_data %>% select(-CellID)
DefaultAssay(seurat) <- "RNA"
seurat_small <- DietSeurat(seurat,
    layers = "data",
    assays = "RNA",
    misc = FALSE
)

# remove all genes that are expressed in less than 15 cells
seurat_small <- seurat_small[rowSums(seurat_small@assays$RNA@data > 0) >= 15, ]

# remove all counts, we will only work with normalized counts.
# this saves a lot of space
empty_matrix <- sparseMatrix(dims = c(nrow(seurat_small),ncol(seurat_small)), i={}, j={})
empty_matrix <- as(empty_matrix, "dgCMatrix")
dimnames(empty_matrix) <- dimnames(seurat_small)
seurat_small <- SetAssayData(seurat_small, slot = "counts", new.data = empty_matrix)


# only meta data that is relevant
seurat_small@meta.data <- seurat_small@meta.data[, c("celltype_manual", "perturbation")]

all(Cells(seurat_small) == rownames(dim_data))

# add the dimension data to the meta data
seurat_small@meta.data <- cbind(seurat_small@meta.data, dim_data)

seurat_small$perturbation <- factor(seurat_small$perturbation, levels = c("ctrl", "NotchRNAi", "NotchCphRNAi", "CphUp"))

# change the gene name of cph
rownames(seurat_small@assays$RNA@data) <-
    ifelse(rownames(seurat_small@assays$RNA@data) == "CG9650", "Cph", rownames(seurat_small@assays$RNA@data))
rownames(seurat_small@assays$RNA@counts) <-
    ifelse(rownames(seurat_small@assays$RNA@counts) == "CG9650", "Cph", rownames(seurat_small@assays$RNA@counts))

saveRDS(seurat_small, here("data", "filtered_seurat_object_RNAi.rds"))

