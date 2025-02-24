library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)
plan("multisession", workers = 25)

To Load te file for another downstreaming task
path <- "/home/fenosoa/scratch/Maya_Project/pipeline_Region/out3/UMAP_xenium3.obj.RData"
temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

xenium.obj <- temp_env$xenium_merge.obj

# SCTransform performs normalization and variance stabilization on the data.
# It replaces the log-normalization step and can help mitigate the effects of technical noise.
print("==================================================================================")
print("# Perform Principal Component Analysis (PCA) on the normalized data")
xenium.obj <- RunPCA(
  xenium.obj,             
  npcs = 30,              # Number of principal components to compute
  features = rownames(xenium.obj)  # Features (genes) to include in the PCA
)
# PCA reduces the dimensionality of the dataset while retaining the most important information.
# Let's compute UMAP using 30 principal components based on all features (genes).
print("==================================================================================")
print("# Compute a UMAP (Uniform Manifold Approximation and Projection) for visualization")
xenium.obj <- RunUMAP(
  xenium.obj,             # The Seurat object
  dims = 1:30             # Use the first 30 principal components for UMAP calculation
)
# It helps identify and visualize clusters or patterns in the data.
print("==================================================================================")
print("# Identify cell neighbors based on the PCA-reduced dimensions")
xenium.obj <- FindNeighbors(
  xenium.obj,             # The Seurat object
  reduction = "pca",      # Use PCA-reduced data for neighbor identification
  dims = 1:30             # Use the first 30 principal components
)
# FindNeighbors builds a shared nearest neighbor (SNN) graph, which is used for clustering.
# It identifies which cells are most similar based on their PCA representation.
print("==================================================================================")
print("# Cluster the cells using the SNN graph")
xenium.obj <- FindClusters(
  xenium.obj,             
  resolution = 0.3        # Clustering resolution (controls the number of clusters)
)
# FindClusters assigns cells to discrete clusters based on the SNN graph.
# The `resolution` parameter affects the granularity of the clustering:
# - Lower resolution (e.g., 0.3) results in fewer, larger clusters.
# - Higher resolution (e.g., 1.0) results in more, smaller clusters.

# After running this workflow, the Seurat object will have normalized data, dimensionality-reduced representations (PCA, UMAP),
UMAP_xenium_merge.obj <- xenium.obj
save(UMAP_xenium_merge.obj, file = "/home/fenosoa/scratch/Maya_Project/pipeline_Region/out3/UMAP_xenium_merge.obj.RData")
