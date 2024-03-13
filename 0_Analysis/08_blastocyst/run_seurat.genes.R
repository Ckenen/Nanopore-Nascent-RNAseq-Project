#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
# load matrix
outdir <- "results/seurat_genes"
f_expr <- paste0("results/blastocyst_counts.genes.gene_name.total.tsv")
f_meta <- paste0("results/blastocyst_counts.genes.gene_name.meta.tsv")
counts <- read.table(f_expr, sep = '\t', header = T, row.names = 1)
meta = read.table(f_meta, sep = '\t', header = T)
rownames(meta) <- colnames(counts)
strains = meta$Strain
for (i in 1:length(meta$Strain)){
    s <- strains[i]
    if (endsWith(s, "X")){
        s <- substr(s, 1, nchar(s) - 1)
    }
    if (s == "Early-blastocyst"){
        s = "Early"
    }
    else if (s == "Mid-blastocyst") {
        s = "Mid"
    }
    else if (s == "Late-blastocyst") {
        s = "Late"
    }
    strains[i] <- s
}

# create Seruat object

obj <- CreateSeuratObject(counts = counts, project = "Blastocyst", min.cells = 3, meta.data = meta)
obj$strain <- strains
obj <- subset(obj, nCount_RNA >= 4000 & nFeature_RNA > 2000)
# obj <- subset(obj, strain == "Early-blastocyst" | strain == "Mid-blastocyst" | strain == "Late-blastocyst")

# VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "strain")
# FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "strain")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 100)
# top10 <- head(VariableFeatures(obj), 10)
# plot1 <- VariableFeaturePlot(obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, max.overlaps=100)
# plot1 + plot2
# plot2

all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
# VizDimLoadings(obj, dims = 1:5, reduction = "pca")
# DimPlot(obj, reduction = "pca", group.by = "strain")
# DimHeatmap(obj, dims = 1:5, cells = 30, balanced = TRUE)

obj <- JackStraw(obj, num.replicate = 100, prop.freq=0.05)
obj <- ScoreJackStraw(obj, dims = 1:20)
# JackStrawPlot(obj, dims = 1:15)
# ElbowPlot(obj)

obj <- FindNeighbors(obj, dims = 1:5)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:5)

# UMAP

pdf(paste0(outdir, "/umap_groupby_clusters.pdf"), width=6, height=6)
DimPlot(obj, reduction = "umap", pt.size = 1, group.by="seurat_clusters", label=TRUE)
dev.off()

pdf(paste0(outdir, "/umap_groupby_strains.pdf"), width=6, height=6)
DimPlot(obj, reduction = "umap", pt.size = 1, group.by="strain")
dev.off()

# marker genes

pdf(paste0(outdir, "/seurat_markers.ICM.pdf"), width=9, height=9)
FeaturePlot(obj, features = c("Nanog", "Sox2", "Fgf4", "Lifr", "Gata4", "Gata6", "Sox17"), ncol=3)
dev.off()

pdf(paste0(outdir, "/seurat_markers.TE.pdf"), width=9, height=6)
FeaturePlot(obj, features = c("Cdx2", "Eomes", "Elf5", "Gata3", "Tfap2a"), ncol=3)
dev.off()

pdf(paste0(outdir, "/seurat_markers.pdf"), width=12, height=6)
FeaturePlot(obj, features = c("Nanog", "Sox2", "Gata4", "Sox17", "Cdx2", "Eomes", "Elf5", "Gata3"), ncol=4)
dev.off()

# annotate celltype
new.cluster.ids <- c("ICM-2", "TE-3", "EPI", "ICM-1", "PE", "TE-1", "TE-2", "ICM-3", "Unknown", "Morula")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

pdf(paste0(outdir, "/umap_groupby_celltype.pdf"), width=6, height=6)
DimPlot(obj, reduction = "umap", pt.size=1, label=TRUE)
dev.off()

# output
write.table(obj@meta.data, paste0(outdir, "/metadata.tsv"), sep = "\t")
write.table(obj@active.ident, paste0(outdir, "/active_ident.tsv"), sep = "\t")
saveRDS(obj, paste0(outdir, "/seurat.rds"))

# re-load

obj <- readRDS(paste0(outdir, "/seurat.rds"))

pdf(paste0(outdir, "/feature_s4u.pdf"), width=6, height=6)
DimPlot(obj, reduction = "umap", pt.size=1, group.by="s4U")
dev.off()

obj_t <- subset(obj, s4U == 400 & Time == 3)
pdf(paste0(outdir, "/feature_nascent_umi_ratio.pdf"), width=6, height=6)
FeaturePlot(obj_t, features = c("NascentRatio2"), pt.size = 1)
dev.off()
pdf(paste0(outdir, "/feature_tc_ratio.pdf"), width=6, height=6)
FeaturePlot(obj_t, features = c("TCRatio"), pt.size = 1)
dev.off()

# filter

obj <- readRDS(paste0(outdir, "/seurat.rds"))
obj <- subset(obj, seurat_clusters != 8 & seurat_clusters != 9)

pdf(paste0(outdir, "/umap_groupby_strains.filtered.pdf"), width=4.5, height=4.5)
DimPlot(obj, reduction="umap", pt.size=0.5, group.by="strain") & NoAxes()
dev.off()

pdf(paste0(outdir, "/umap_groupby_celltype.filtered.pdf"), width=4.5, height=4.5)
DimPlot(obj, reduction = "umap", pt.size = 1, label=TRUE, label.box=F) & NoAxes() & NoLegend()
dev.off()

pdf(paste0(outdir, "/seurat_markers.filtered.pdf"), width=14, height=4.5)
FeaturePlot(obj, features = c("Nanog", "Sox2", "Tdgf1", "Fgf4", "Gata4", "Sox17", "Gata6", "Cdx2", "Gata3", "Eomes", "Elf5", "Tfap2a"), ncol=6, pt.size=0.2) & NoAxes()
dev.off()

pdf(paste0(outdir, "/feature_s4u.filtered.pdf"), width=4.5, height=4.5)
DimPlot(obj, reduction = "umap", pt.size=0.5, group.by="s4U") & NoAxes() 
dev.off()

obj_t <- subset(obj, s4U == 400 & Time == 3)
pdf(paste0(outdir, "/feature_nascent_umi_ratio.filtered.pdf"), width=4.5, height=4.5)
FeaturePlot(obj_t, features = c("NascentRatio2"), pt.size=0.5) + NoAxes()
dev.off()
pdf(paste0(outdir, "/feature_tc_ratio.filtered.pdf"), width=4.5, height=4.5)
FeaturePlot(obj_t, features = c("TCRatio"), pt.size=0.5) + NoAxes()
dev.off()

markers <- FindAllMarkers(obj, only.pos = TRUE)
write.csv(markers, paste0(outdir, "/markers.csv"))

markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 8) %>% ungroup() -> top10
pdf(paste0(outdir, "/heatmap_markers.pdf"), width=7, height=9)
DoHeatmap(obj, features = top10$gene, size=3) # + NoLegend()
dev.off()


# FeaturePlot(obj, features = c("Folr1", "Ccdc12", "Nln", "Prrc1", "Pdgfa", "Cdk2ap1", "Tpm1", "Serpinb6c", "Ccdc28a", "Wtap", "Plac1"), ncol=4)

