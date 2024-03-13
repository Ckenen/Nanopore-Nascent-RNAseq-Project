#!/usr/bin/env Rscript
library(ggplot2)
library(Seurat)
library(dplyr)
# load matrix
outdir <- "results/seurat_transcripts"
f_expr <- paste0("results/blastocyst_counts.transcripts.transcript_name.total.tsv")
f_meta <- paste0("results/blastocyst_counts.transcripts.transcript_name.meta.tsv")
#f_meta <- paste0("expression/mouse_blastocyst_transcript_meta.celltype.tsv")
counts <- read.table(f_expr, sep = '\t', header = T, row.names = 1)
meta = read.table(f_meta, sep = '\t', header = T)
rownames(meta) <- colnames(counts)
strains = meta$Strain
for (i in 1:length(meta$Strain)){
    s <- strains[i]
    if (endsWith(s, "X")){
        s <- substr(s, 1, nchar(s) - 1)
    }
    strains[i] <- s
}

# create Seruat object
obj <- CreateSeuratObject(counts = counts, project = "Blastocyst", min.cells = 5, meta.data = meta)
obj$strain <- strains
obj <- subset(obj, nCount_RNA >= 2000 & nFeature_RNA > 1000)
obj <- subset(obj, CellTypeFromGene != "Morula" & CellTypeFromGene != "Unknown" & CellTypeFromGene != "")
# obj <- subset(obj, strain == "Early-blastocyst" | strain == "Mid-blastocyst" | strain == "Late-blastocyst" | strain == "Morula" | strain == "Early-blastocystX" | strain == "Mid-blastocystX" | strain == "Late-blastocystX" | strain == "mESC")
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
obj <- FindClusters(obj, resolution = 0.7)
obj <- RunUMAP(obj, dims = 1:5)
obj@reductions$umap@cell.embeddings[,1] <- -1 * obj@reductions$umap@cell.embeddings[,1]

# UMAP

pdf(paste0(outdir, "/umap_groupby_clusters.pdf"), width=6, height=6)
DimPlot(obj, reduction="umap", pt.size=1, group.by="seurat_clusters", label=TRUE)
dev.off()

pdf(paste0(outdir, "/umap_groupby_celltype_from_gene.pdf"), width=6, height=6)
DimPlot(obj, reduction="umap", pt.size=1, group.by="CellTypeFromGene")
dev.off()

# annotate

new.cluster.ids <- c("ICM-2", "EPI", "TE-2", "PE", "ICM-1", "muralTE", "polarTE", "TE-1")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

pdf(paste0(outdir, "/umap_groupby_celltype.pdf"), width=4.5, height=4.5)
cols <- c("#F8766D", "#7CAE00", "#C77CFF", "#00BFC4", "#00BE67", "#CD969B", "#839600", "#00A9FF")
DimPlot(obj, reduction = "umap", pt.size=1, label=T, label.box=F, cols=cols) & NoAxes() & NoLegend()
dev.off()

# output

write.table(obj@meta.data, paste0(outdir, "/metadata.tsv"), sep = "\t")
write.table(obj@active.ident, paste0(outdir, "/active_ident.tsv"), sep = "\t")
saveRDS(obj, paste0(outdir, "/seurat.rds"))

# re-load

obj <- readRDS(paste0(outdir, "/seurat.rds"))

# marker isoforms

pdf(paste0(outdir, "/seurat_markers.pdf"), width=14, height=4.5)
FeaturePlot(obj, features = c("Sox2-201", "Tdgf1-201", "Fgf4-201", "Sox17-201", "Sox17-204", "Gata6-201", "Cdx2-201", "Gata3-201", "Eomes-201", "Eomes-202", "Eomes-203", "Elf5-201"), ncol=6, pt.size=0.2, label=F) & NoAxes()
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(obj, only.pos = TRUE)
write.csv(markers, paste0(outdir, "/markers.csv"))
markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 8) %>% ungroup() -> top10
pdf(paste0(outdir, "/heatmap_markers.pdf"), width=7, height=9)
DoHeatmap(obj, features = top10$gene, size=3, group.colors=cols) # + NoLegend()
dev.off()

# Moral TE
DotPlot(obj, features = c('Adam19-201', 'Ahnak-202', 'Akr1b8-201', 'Anxa2-201', 'Aqp3-201', 'BC051665-201', 'BC053393-201', 'Basp1-201', 'C430049B03Rik-201', 'C430049B03Rik-204', 'Capg-202', 'Cdc42se2-201', 'Chmp5-201', 'Cited2-201', 'Cldn4-201', 'Crxos-201', 'Crxos-204', 'Crxos-205', 'Ctsll3-201', 'Dag1-204', 'Dkkl1-201', 'Dkkl1-204', 'Dmkn-201', 'Dppa1-201', 'Dppa1-202', 'Efhd2-201', 'Emp2-201', 'Entpd1-201', 'Fabp3-201', 'Folr1-202', 'Fxyd4-203', 'Fxyd4-207', 'Glipr2-201', 'Gm2a-201', 'Gm4926-201', 'Gnai2-201', 'Gnai2-202', 'Irx3-202', 'Krt18-201', 'Krt19-201', 'Krt7-201', 'Krt8-201', 'Lgals1-201', 'Lgals1-202', 'Lgals3-201', 'Lgals3-204', 'Litaf-201', 'Mbnl3-201', 'Ndrg1-201', 'Nfkbia-201', 'Nppb-201', 'Plac1-201', 'Plac1-202', 'Plac8-201', 'Psap-202', 'Psap-205', 'Ptgr1-201', 'Rab18-205', 'Rhou-201', 'S100a11-201', 'S100a6-201', 'Sct-202', 'Sct-204', 'Serpinb6c-201', 'Serpinb6c-202', 'Sh3bgrl3-201', 'Sin3b-202', 'Slc1a1-201', 'Slc6a14-201', 'Smpdl3a-201', 'Tagln2-202', 'Tfap2a-201', 'Tfap2a-203', 'Tfap2c-201', 'Tinagl1-203', 'Tmsb4x-203', 'Tpm1-209', 'Tpm1-211', 'Tpm1-212', 'Tpm1-213', 'Tspan8-201', 'Tspan8-202', 'Tspan8-203', 'Zfp703-202')) + RotatedAxis()

# Polar TE
DotPlot(obj, features = c('Aprt-201', 'Ass1-201', 'Ddah1-201', 'Gsto1-201', 'Ly6a-201', 'Ly6a-203', 'Ly6a-204')) + RotatedAxis()

# Final
pdf(paste0(outdir, "/mte_pte_markers.pdf"), width=7.2, height=4)
DotPlot(obj, features = c('Basp1-201', 'BC051665-201', 'Cited2-201', 'Dmkn-201', 'Gm4926-201', 'Krt7-201', 'Krt19-201', 'Plac1-201', 'Aprt-201', 'Ass1-201', 'Ddah1-201', 'H19-202')) + RotatedAxis() + xlab(NULL) + ylab(NULL)
dev.off()

# 
tmp <- subset(obj, seurat_clusters == 5 | seurat_clusters == 6)
# VlnPlot(tmp, features = c('Gata3-201', 'Rxra-201')) + RotatedAxis()
pdf(paste0(outdir, "/gata3.pdf"), width=6.5, height=4)
VlnPlot(tmp, features = c('Gata3-201', 'Dmkn-201', 'Plac1-201', 'H19-202'), ncol=4) & RotatedAxis() & xlab(NULL) & NoLegend()
dev.off()


# DotPlot(obj, features = c("Snhg5-202", "Snhg5-225", "Snhg5-226", "Vma21-202", "Vma21.novel.116053", "Rpl41-201", "Rpl41-206", "Rpl41.novel.12445", "Srebf2-201", "Srebf2-205", "Prrc1-201", "Prrc1-205", "Prrc1.novel.52373", "Actn1-201", "Actn1-202", "Folr1-201", "Folr1-202", "Folr1.novel.98548", "H3f3b-203", "H3f3b.novel.22137", "H3f3b.novel.22144", "Ccdc12-201", "Ccdc12-210", "Pstk-201", "Pstk.novel.100777", "Pdgfa-201", "Pdgfa-202"  )) + RotatedAxis() + xlab(NULL) + ylab(NULL)
