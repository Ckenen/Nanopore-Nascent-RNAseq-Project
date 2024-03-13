#!/usr/bin/env Rscript
library(Seurat)
library(monocle3)

obj <- readRDS("seurat_outputs/seurat_annotated_obj.rds")

data <- GetAssayData(obj, assay = "RNA", slot = "counts")
cell_metadata <- obj@meta.data
cell_metadata$celltype <- obj@active.ident
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim=50)
#plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, preprocess_method= "PCA")
#plot_cells(cds)
#plot_cells(cds, reduction_method = "UMAP", color_cells_by = "celltype")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
#plot_cells(cds, reduction_method = "UMAP", color_cells_by = "celltype")

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)





