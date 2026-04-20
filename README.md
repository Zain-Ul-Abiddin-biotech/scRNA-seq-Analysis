# 🧬 Single-Cell RNA-seq Analysis of NSCLC — Step-by-Step Seurat v5 Pipeline

> **Who is this for?**
> This guide is written for **beginners in R** who want to understand, run, and customize a single-cell RNA sequencing (scRNA-seq) analysis pipeline. Every section explains *what* the code does, *why* it's done, and *what you can change* to suit your own data.

---

## 📋 Table of Contents

1. [What is this pipeline?](#1-what-is-this-pipeline)
2. [Requirements — what you need to install](#2-requirements--what-you-need-to-install)
3. [The Dataset](#3-the-dataset)
4. [Step-by-Step Code Walkthrough](#4-step-by-step-code-walkthrough)
   - [Step 0 — Load Libraries](#step-0--load-libraries)
   - [Step 1 — Load the Data](#step-1--load-the-data)
   - [Step 2 — Quality Control (QC)](#step-2--quality-control-qc)
   - [Step 3 — Filtering](#step-3--filtering)
   - [Step 4 — Normalize the Data](#step-4--normalize-the-data)
   - [Step 5 — Find Variable Features](#step-5--find-variable-features)
   - [Step 6 — Scale the Data](#step-6--scale-the-data)
   - [Step 7 — PCA (Linear Dimensionality Reduction)](#step-7--pca-linear-dimensionality-reduction)
   - [Step 8 — Clustering](#step-8--clustering)
   - [Step 9 — UMAP (Non-linear Dimensionality Reduction)](#step-9--umap-non-linear-dimensionality-reduction)
5. [Tweaking the Pipeline for Your Data](#5-tweaking-the-pipeline-for-your-data)
6. [Common Errors and Fixes](#6-common-errors-and-fixes)
7. [Project Structure](#7-project-structure)

---

## 1. What is this pipeline?

This R script performs a standard **single-cell RNA sequencing (scRNA-seq)** analysis on a **Non-Small Cell Lung Cancer (NSCLC)** dataset.

In simple terms: you have thousands of individual cells, and for each cell you know how active (expressed) each gene is. This pipeline helps you:

- **Clean** the data (remove dead or broken cells)
- **Normalize** the data (make cells comparable to each other)
- **Find patterns** (which genes vary the most between cells?)
- **Group cells** into clusters (cells that look similar go together)
- **Visualize** the clusters on a 2D plot (UMAP)

The goal is to discover **what types of cells** are present in the tumor sample.

---

## 2. Requirements — what you need to install

### R version
This pipeline requires **R 4.3 or higher** and **Seurat v5**.

### Installing packages
Run this once in R before running the script:

```r
# Install Seurat v5
install.packages("Seurat")

# Install tidyverse (data manipulation + ggplot2)
install.packages("tidyverse")

# Install hdf5r (needed to read .h5 files)
install.packages("hdf5r")

# Install ggplot2 (plotting)
install.packages("ggplot2")

# Optional: for UMAP support
reticulate::py_install(packages = "umap-learn")
```

> ⚠️ **Windows users:** Installing `hdf5r` may require HDF5 libraries.
> Download them from: https://www.hdfgroup.org/downloads/hdf5/

---

## 3. The Dataset

This pipeline uses the **20k NSCLC (Non-Small Cell Lung Cancer) Dissociated Tumor Cells** dataset from 10x Genomics.

| Property | Value |
|---|---|
| File format | `.h5` (HDF5 sparse matrix) |
| Cells (approx.) | ~20,000 |
| Sequencing | 3' gene expression (Next GEM Multiplex) |
| Source | [10x Genomics Dataset Page](https://www.10xgenomics.com/datasets) |

**To change the dataset**, update this line in the script with your own file path:

```r
# 👇 Change this path to point to your own .h5 file
nsclc.sparse.m <- Read10X_h5(filename = "C:/Your/Path/To/your_file.h5")
```

If your data is **not in .h5 format** (e.g. it's a folder with `matrix.mtx`, `barcodes.tsv`, `features.tsv`), use:

```r
# For 10x Genomics standard output folder format
cts <- Read10X(data.dir = "C:/Your/Path/To/filtered_feature_bc_matrix/")
```

---

## 4. Step-by-Step Code Walkthrough

---

### Step 0 — Load Libraries

```r
library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)
```

**What this does:** Loads the toolboxes (packages) the script needs.

| Package | Purpose |
|---|---|
| `Seurat` | The main scRNA-seq analysis toolkit |
| `tidyverse` | Data wrangling helpers (includes `dplyr`, `ggplot2`) |
| `hdf5r` | Lets R read `.h5` files |
| `ggplot2` | Beautiful plots |

---

### Step 1 — Load the Data

```r
nsclc.sparse.m <- Read10X_h5(filename = "...your_file.h5")
cts <- nsclc.sparse.m$`Gene Expression`

nsclc.seurat.obj <- CreateSeuratObject(
  counts       = cts,
  project      = "NSCLC",
  min.cells    = 3,
  min.features = 200
)
```

**What this does:**
1. Reads the raw count matrix from the `.h5` file — this is a table of cells × genes
2. Extracts just the Gene Expression data (the `.h5` file may also contain protein data)
3. Creates a **Seurat object** — the main container that holds your data and results throughout the pipeline

**Parameters you can tweak:**

| Parameter | Default | What it controls |
|---|---|---|
| `project` | `"NSCLC"` | Just a label — change to your project name |
| `min.cells` | `3` | A gene must appear in at least this many cells to be kept. Increase to filter rare genes. |
| `min.features` | `200` | A cell must express at least this many genes to be kept. Increase to remove near-empty cells. |

---

### Step 2 — Quality Control (QC)

```r
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")

VlnPlot(nsclc.seurat.obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
```

**What this does:**
- `percent.mt` — calculates what percentage of each cell's reads come from **mitochondrial genes** (genes starting with `MT-`). A very high percentage usually means the cell is dying or damaged.
- `VlnPlot` — shows violin plots of three QC metrics:
  - `nFeature_RNA` = how many unique genes a cell has
  - `nCount_RNA` = total number of RNA molecules detected
  - `percent.mt` = % mitochondrial reads
- `FeatureScatter` — a scatter plot to check that cells with more counts also have more genes (they should correlate)

**What to look for in the plots:**
- Cells with **very low** `nFeature_RNA` → likely empty droplets or dead cells
- Cells with **very high** `nFeature_RNA` → likely doublets (two cells captured together)
- Cells with **high** `percent.mt` → likely dying cells

> 💡 **Tip for mouse data:** Change `"^MT-"` to `"^mt-"` (lowercase) since mouse mitochondrial genes use lowercase.

---

### Step 3 — Filtering

```r
nsclc.seurat.obj <- subset(
  nsclc.seurat.obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 2500 &
           percent.mt < 5
)
```

**What this does:** Removes low-quality cells based on what you saw in the QC plots.

**Parameters you can tweak** — these are the most important numbers to adjust for your own data:

| Filter | Default | When to change |
|---|---|---|
| `nFeature_RNA > 200` | 200 | Increase if your data has many empty droplets |
| `nFeature_RNA < 2500` | 2500 | Increase (e.g. to 5000) for datasets with naturally complex cells like neurons |
| `percent.mt < 5` | 5% | Some tissues (e.g. heart muscle) have naturally high mitochondrial content — you may need to raise this to 10–25% |

> 💡 **How to decide your thresholds:** Look at the violin plots from Step 2. The cutoffs should be placed where the distribution separates "good" cells from outliers.

---

### Step 4 — Normalize the Data

```r
nsclc.seurat.obj <- NormalizeData(
  nsclc.seurat.obj,
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
```

**What this does:** Different cells capture different total amounts of RNA — some cells just happen to have more molecules sequenced than others. Normalization corrects for this so cells are comparable.

The formula used is:
```
normalized value = log( (raw count / total counts per cell) × 10,000 + 1 )
```

**Parameters you can tweak:**

| Parameter | Default | Notes |
|---|---|---|
| `normalization.method` | `"LogNormalize"` | Standard choice for most datasets. Alternatives: `"CLR"` for CITE-seq protein data |
| `scale.factor` | `10000` | Rarely needs changing |

---

### Step 5 — Find Variable Features

```r
nsclc.seurat.obj <- FindVariableFeatures(
  nsclc.seurat.obj,
  selection.method = "vst",
  nfeatures        = 2000
)

top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
labeled_plot <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(labeled_plot)
```

**What this does:** Identifies the **2000 most variable genes** — genes that differ the most between cells. These are the most informative genes for distinguishing cell types. Less variable genes are mostly just noise.

**Parameters you can tweak:**

| Parameter | Default | Notes |
|---|---|---|
| `nfeatures` | `2000` | Increase to `3000–5000` for more complex datasets with many cell types |
| `selection.method` | `"vst"` | Best default choice — uses variance stabilizing transformation |

The plot shows genes on a variability scale. The labeled top 10 are often well-known marker genes for the cell types in your sample.

---

### Step 6 — Scale the Data

```r
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
```

**What this does:** Scales each gene so it has a **mean of 0 and standard deviation of 1**. This prevents highly expressed genes from dominating the analysis just because they have large absolute values.

> ⚠️ **Note:** Scaling `all.genes` is thorough but uses more memory. For large datasets (50k+ cells), you can scale only variable features:
> ```r
> nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj)  # defaults to variable features only
> ```

**Optional — regress out unwanted variation** (e.g. cell cycle effects):
```r
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj,
                               vars.to.regress = "percent.mt")
```

---

### Step 7 — PCA (Linear Dimensionality Reduction)

```r
nsclc.seurat.obj <- RunPCA(
  nsclc.seurat.obj,
  features = VariableFeatures(object = nsclc.seurat.obj)
)

print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(nsclc.seurat.obj)
```

**What this does:** PCA (Principal Component Analysis) compresses the data from 2000 gene dimensions down to a smaller number of "principal components" (PCs). Each PC captures a pattern of variation in the data.

- `print(...)` — shows the top genes driving each PC
- `DimHeatmap(...)` — heatmap showing which cells and genes drive PC1
- `ElbowPlot(...)` — **the most important plot here** — shows you how many PCs are meaningful

**How to read the Elbow Plot:**
Look for the "elbow" — the point where the curve flattens out. The number of PCs at the elbow is what you should use in the next steps (typically 10–20).

```r
# 👇 Change dims = 1:15 in the next steps based on your ElbowPlot
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)
```

---

### Step 8 — Clustering

```r
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = 0.5)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.5"
```

**What this does:**
1. `FindNeighbors` — builds a graph connecting similar cells based on their PCA coordinates
2. `FindClusters` — groups connected cells into clusters using the Louvain algorithm
3. `DimPlot` — visualizes the clusters (though UMAP in Step 9 gives a better view)

**The most important parameter to tweak — `resolution`:**

| `resolution` | Effect |
|---|---|
| `0.1 – 0.3` | Fewer, broader clusters (good for coarse cell type detection) |
| `0.5` | Default — a balanced starting point |
| `0.7 – 1.0` | More clusters (good for finding rare subtypes) |
| `1.0+` | Many fine-grained clusters (risk of over-splitting) |

> 💡 **Tip:** Run multiple resolutions and compare:
> ```r
> nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj,
>                                   resolution = c(0.1, 0.3, 0.5, 0.7, 1.0))
> ```

---

### Step 9 — UMAP (Non-linear Dimensionality Reduction)

```r
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap", label = TRUE)
```

**What this does:** UMAP projects the high-dimensional PCA data into **2 dimensions** that you can plot. Cells that are similar to each other end up close together on the UMAP plot. Each dot is one cell, colored by its cluster.

**What to do after you get your UMAP:**
- Each cluster represents a group of similar cells
- You can identify what cell type each cluster is by looking at marker genes:
```r
# Find marker genes for every cluster
markers <- FindAllMarkers(nsclc.seurat.obj, only.pos = TRUE, min.pct = 0.25)

# Find markers for a specific cluster (e.g. cluster 3)
cluster3_markers <- FindMarkers(nsclc.seurat.obj, ident.1 = 3, min.pct = 0.25)
```

- Then visualize known marker genes:
```r
# e.g. EPCAM marks epithelial/cancer cells, CD3D marks T cells
FeaturePlot(nsclc.seurat.obj, features = c("EPCAM", "CD3D", "CD68", "PECAM1"))
```

---

## 5. Tweaking the Pipeline for Your Data

Here's a quick reference for the most common adjustments:

| What you want to change | Parameter to edit | Location in script |
|---|---|---|
| Input file path | `filename = "..."` | Step 1 |
| Minimum genes per cell | `min.features = 200` | Step 1 |
| Mitochondrial gene filter | `percent.mt < 5` | Step 3 |
| Max genes per cell | `nFeature_RNA < 2500` | Step 3 |
| Number of variable genes | `nfeatures = 2000` | Step 5 |
| Number of PCs to use | `dims = 1:15` | Steps 7, 8, 9 |
| Clustering resolution | `resolution = 0.5` | Step 8 |
| Mouse vs human mito genes | `"^MT-"` → `"^mt-"` | Step 2 |

---

## 6. Common Errors and Fixes

| Error message | Likely cause | Fix |
|---|---|---|
| `could not find function "Read10X_h5"` | Seurat not loaded | Run `library(Seurat)` |
| `unable to open file` | Wrong file path | Check your path — use `/` not `\` on Windows, or use `\\` |
| `cannot allocate vector of size X Gb` | Not enough RAM | Scale only variable features (see Step 6); reduce `nfeatures` |
| `Python / UMAP not found` | umap-learn not installed | Run `reticulate::py_install("umap-learn")` |
| `object of type 'S4' is not subsettable` | Object overwritten | Restart R, rerun from the top — likely caused by the old chunking bug |
| Plots appear blank or in wrong window | RStudio plot pane issue | Run `dev.new()` before plotting, or use `print(plot_object)` |

---

## 7. Project Structure

```
📁 your-project/
│
├── 📄 nsclc_scRNA_corrected.R       ← Main analysis script
├── 📄 README.md                     ← This file
│
└── 📁 data/
    └── 20k_NSCLC_...matrix.h5      ← Raw data file (download separately)
```

> ⚠️ **Do not upload your `.h5` data file to GitHub** — it is several gigabytes in size. Add it to your `.gitignore`:
> ```
> *.h5
> *.h5ad
> ```

---

## 📚 Further Reading

- [Seurat v5 official vignettes](https://satijalab.org/seurat/articles/get_started.html)
- [10x Genomics dataset portal](https://www.10xgenomics.com/datasets)
- [Orchestrating Single-Cell Analysis with Bioconductor (OSCA)](https://bioconductor.org/books/release/OSCA/) — free online textbook

---

## 🤝 Contributing

Feel free to fork this repository and adapt the pipeline for your own dataset. If you find a bug or improvement, open an Issue or Pull Request.

---

*Pipeline corrected and documented for Seurat v5 compatibility.*
