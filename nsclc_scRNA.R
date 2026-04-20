# ============================================================
# scRNA-seq Analysis: NSCLC Dataset
# Corrected for Seurat v5 + current R best practices
# ============================================================

# Load libraries
library(Seurat)       # v5+
library(tidyverse)
library(hdf5r)
library(ggplot2)

# ── Load the NSCLC dataset ───────────────────────────────────
nsclc.sparse.m <- Read10X_h5(
  filename = "C:/Zain/Codes/R/Projects/scRNA Analaysis/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"
)
str(nsclc.sparse.m)

cts <- nsclc.sparse.m$`Gene Expression`

# ── Initialize Seurat object (raw counts) ───────────────────
nsclc.seurat.obj <- CreateSeuratObject(
  counts     = cts,
  project    = "NSCLC",
  min.cells  = 3,
  min.features = 200
)
str(nsclc.seurat.obj)


# ── 1. QC ────────────────────────────────────────────────────
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(
  nsclc.seurat.obj,
  pattern = "^MT-"
)

VlnPlot(
  nsclc.seurat.obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(
  nsclc.seurat.obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) + geom_smooth(method = "lm")


# ── 2. Filtering ─────────────────────────────────────────────
nsclc.seurat.obj <- subset(
  nsclc.seurat.obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 2500 &
           percent.mt  < 5
)


# ── 3. Normalize data ────────────────────────────────────────
# FIX: explicitly specify normalization.method and scale.factor (v5 best practice)
nsclc.seurat.obj <- NormalizeData(
  nsclc.seurat.obj,
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)


# ── 4. Identify highly variable features ────────────────────
nsclc.seurat.obj <- FindVariableFeatures(
  nsclc.seurat.obj,
  selection.method = "vst",
  nfeatures        = 2000
)

top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# FIX: must assign combined plot to a variable and then print it,
#      otherwise LabelPoints output is silently discarded
plot1        <- VariableFeaturePlot(nsclc.seurat.obj)
labeled_plot <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(labeled_plot)


# ── 5. Scaling ───────────────────────────────────────────────
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)


# ── 6. Linear dimensionality reduction (PCA) ────────────────
# FIX: Removed the invalid manual chunking + CreateDimReducObject block.
#
#   WHY IT WAS WRONG:
#     (a) CreateDimReducObject() returns a DimReduc object, NOT a Seurat
#         object — assigning it back to nsclc.seurat.obj destroyed the object.
#     (b) Running PCA independently on each chunk produces incompatible
#         loading vectors; merging cell embeddings across chunks is
#         statistically meaningless and gives incorrect results.
#     (c) Seurat v5 handles large datasets efficiently on its own via
#         irlba (sparse truncated SVD) — manual chunking is unnecessary.
#
#   SOLUTION: run RunPCA directly on the full object.

nsclc.seurat.obj <- RunPCA(
  nsclc.seurat.obj,
  features = VariableFeatures(object = nsclc.seurat.obj)
)

# FIX: In Seurat v5 the correct accessor syntax uses [[ ]] not @reductions
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# Determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)


# ── 7. Clustering ────────────────────────────────────────────
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj,  resolution = 0.5)

# Visualize clusters
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# Setting identity of clusters
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.5"
print(Idents(nsclc.seurat.obj))


# ── 8. Non-linear dimensionality reduction (UMAP) ───────────
# If umap-learn is not installed:
#   reticulate::py_install(packages = "umap-learn")
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

DimPlot(nsclc.seurat.obj, reduction = "umap", label = TRUE)


# ── Optional: multiple clustering resolutions ───────────────
# nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj,
#                                   resolution = c(0.1, 0.3, 0.7, 1))
# View(nsclc.seurat.obj@meta.data)
