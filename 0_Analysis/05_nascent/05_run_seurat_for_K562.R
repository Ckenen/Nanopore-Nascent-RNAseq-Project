#!/usr/bin/env Rscript

library(Seurat)

outdir <- "seurat_outputs"
f_expr <- "data/K562_counts.tsv"
f_meta <- "data/K562_meta.tsv"
counts <- read.table(f_expr, sep = '\t', header = T, row.names = 1)
meta = read.table(f_meta, sep = '\t', header = T)
rownames(meta) <- colnames(counts)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# create Seruat object

m <- CreateSeuratObject(counts = counts, project = "NanoNASCseq_K562", min.cells = 10, meta.data = meta)
m <- subset(m, nCount_RNA >= 4000 & nFeature_RNA > 2000)

# VlnPlot(m, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "Strain")
# FeatureScatter(m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Strain")

m <- NormalizeData(m, normalization.method = "LogNormalize", scale.factor = 10000)
m <- FindVariableFeatures(m, selection.method = "vst", nfeatures = 2000)

# top10 <- head(VariableFeatures(m), 10)
# plot1 <- VariableFeaturePlot(m)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, max.overlaps=100)
# plot2

m <- ScaleData(m, features = rownames(m))
m <- RunPCA(m, features = VariableFeatures(m))

# VizDimLoadings(m, dims = 1:6, reduction = "pca")
# DimPlot(m, reduction = "pca", group.by = "Strain")
# DimHeatmap(m, dims = 1:6, cells = 30, balanced = TRUE)

m <- CellCycleScoring(m, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(m, reduction = "pca", group.by = "Phase")

m <- JackStraw(m, num.replicate = 100, prop.freq=0.05)
m <- ScoreJackStraw(m, dims = 1:20)
JackStrawPlot(m, dims = 1:20)
ElbowPlot(m)

m <- FindNeighbors(m, dims = 1:4)
m <- FindClusters(m, resolution = 0.5)
m <- RunUMAP(m, dims = 1:4, )

DimPlot(m, reduction = "umap", pt.size = 1, group.by="seurat_clusters", label=TRUE)
DimPlot(m, reduction = "umap", pt.size = 1, group.by="Phase", label=TRUE)

# Cell-cycle

m <- RunPCA(m, features = c(s.genes, g2m.genes))
m <- JackStraw(m, num.replicate = 100, prop.freq=0.05)
m <- ScoreJackStraw(m, dims = 1:20)
JackStrawPlot(m, dims = 1:20)
ElbowPlot(m)

m <- FindNeighbors(m, dims = 1:4)
m <- FindClusters(m, resolution = 0.5)
m <- RunUMAP(m, dims = 1:4, )

DimPlot(m, reduction = "umap", pt.size = 1, group.by="seurat_clusters", label=TRUE)
DimPlot(m, reduction = "umap", pt.size = 1, group.by="Phase", label=TRUE)

FeaturePlot(m, features = c("E2F1", "E2F2", "AURKA"), ncol=3)
VlnPlot(m, features = c("AURKA", "BUB1", "CDC20", "CDCA2", "CDK1", "CENPE", "MKI67", "TOP2A", "BRCA1", "MYBL2"), group.by = "Phase")

VlnPlot(m, features = c(cc.genes$s.genes[1:6]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$s.genes[7:12]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$s.genes[13:18]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$s.genes[19:24]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$s.genes[25:30]), group.by = "Phase")

VlnPlot(m, features = c(cc.genes$g2m.genes[1:6]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$g2m.genes[7:12]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$g2m.genes[13:18]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$g2m.genes[19:24]), group.by = "Phase")
VlnPlot(m, features = c(cc.genes$g2m.genes[25:30]), group.by = "Phase")

VlnPlot(m, features = c("UNG", "DTL", "HELLS", "ATAD2", "KIF11", "CKAP2", "CENPE", "CDCA3"), group.by = "Phase", ncol = 4)

m1 <- subset(m, s4U == 50 & TCRatio > 0.008)
DimPlot(m1, reduction = "umap", pt.size = 1, group.by="Phase", label=TRUE)

markers <- FindAllMarkers(m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )

# $s.genes
#  [1] "MCM5"     "PCNA"     "TYMS"     "FEN1"     "MCM2"     "MCM4"    
#  [7] "RRM1"     "UNG"      "GINS2"    "MCM6"     "CDCA7"    "DTL"     
# [13] "PRIM1"    "UHRF1"    "MLF1IP"   "HELLS"    "RFC2"     "RPA2"    
# [19] "NASP"     "RAD51AP1" "GMNN"     "WDR76"    "SLBP"     "CCNE2"   
# [25] "UBR7"     "POLD3"    "MSH2"     "ATAD2"    "RAD51"    "RRM2"    
# [31] "CDC45"    "CDC6"     "EXO1"     "TIPIN"    "DSCC1"    "BLM"     
# [37] "CASP8AP2" "USP1"     "CLSPN"    "POLA1"    "CHAF1B"   "BRIP1"   
# [43] "E2F8"    

# $g2m.genes
#  [1] "HMGB2"   "CDK1"    "NUSAP1"  "UBE2C"   "BIRC5"   "TPX2"    "TOP2A"  
#  [8] "NDC80"   "CKS2"    "NUF2"    "CKS1B"   "MKI67"   "TMPO"    "CENPF"  
# [15] "TACC3"   "FAM64A"  "SMC4"    "CCNB2"   "CKAP2L"  "CKAP2"   "AURKB"  
# [22] "BUB1"    "KIF11"   "ANP32E"  "TUBB4B"  "GTSE1"   "KIF20B"  "HJURP"  
# [29] "CDCA3"   "HN1"     "CDC20"   "TTK"     "CDC25C"  "KIF2C"   "RANGAP1"
# [36] "NCAPD2"  "DLGAP5"  "CDCA2"   "CDCA8"   "ECT2"    "KIF23"   "HMMR"   
# [43] "AURKA"   "PSRC1"   "ANLN"    "LBR"     "CKAP5"   "CENPE"   "CTCF"   
# [50] "NEK2"    "G2E3"    "GAS2L3"  "CBX5"    "CENPA"  

write.csv(m@meta.data, paste0(outdir, "/meta.data.csv"))

# Below is draft

# UMAP

pdf(paste0(outdir, "/seurat_clusters.pdf"), width=6, height=6)
DimPlot(obj, reduction = "umap", pt.size = 1, group.by="seurat_clusters", label=TRUE)
dev.off()

pdf(paste0(outdir, "/seurat_strains.pdf"), width=6, height=6)
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

# annotate

new.cluster.ids <- c("ICM", "EPI", "TE-2", "TE-1", "TE-3", "PE-2", "PE-1", "Unknown", "Morula")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

pdf(paste0(outdir, "/seurat_annotated_clusters.pdf"), width=6, height=6)
DimPlot(obj, reduction = "umap", pt.size=1, label=TRUE)
dev.off()

# output

f_seurat_meta <- paste0(outdir, "/seurat_metadata.tsv")
f_seurat_ident <- paste0(outdir, "/seurat_active_ident.tsv")
write.table(obj@meta.data, f_seurat_meta, sep = "\t")
write.table(obj@active.ident, f_seurat_ident, sep = "\t")

saveRDS(obj, paste0(outdir, "/seurat_annotated_obj.rds"))

# re-load

obj <- readRDS(paste0(outdir, "/seurat_annotated_obj.rds"))

pdf(paste0(outdir, "/feature_s4u.pdf"), width=6, height=6)
FeaturePlot(obj, features = c("s4U"), pt.size = 1)
dev.off()

obj_t <- subset(obj, s4U == 400 & Time == 3)
pdf(paste0(outdir, "/feature_nascent_rna_ratio.pdf"), width=6, height=6)
FeaturePlot(obj_t, features = c("NascentRNAs..."), pt.size = 1)
dev.off()

# filter

obj1 <- subset(obj, seurat_clusters != 7 & seurat_clusters != 8)

pdf(paste0(outdir, "/seurat_filtered_strains.pdf"), width=4.5, height=4.5)
DimPlot(obj1, reduction = "umap", pt.size = 1, group.by="strain")
dev.off()

pdf(paste0(outdir, "/seurat_filtered_annotated_clusters.pdf"), width=4.5, height=4.5)
DimPlot(obj1, reduction = "umap", pt.size = 1, label=TRUE)
dev.off()

pdf(paste0(outdir, "/seurat_filtered_markers.pdf"), width=9, height=4.5)
FeaturePlot(obj1, features = c("Nanog", "Sox2", "Gata4", "Sox17", "Cdx2", "Eomes", "Elf5", "Gata3"), ncol=4, pt.size = 1)
dev.off()

pdf(paste0(outdir, "/filtered_feature_s4u.pdf"), width=4.5, height=4.5)
FeaturePlot(obj1, features = c("s4U"), pt.size = 1)
dev.off()

obj_t <- subset(obj1, s4U == 400 & Time == 3)
pdf(paste0(outdir, "/filtered_feature_nascent_rna_ratio.pdf"), width=4.5, height=4.5)
FeaturePlot(obj_t, features = c("NascentRNAs..."), pt.size = 1)
dev.off()