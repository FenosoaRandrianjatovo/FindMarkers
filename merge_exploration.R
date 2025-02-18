library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)
plan("multisession", workers = 25)

setwd("/scratch/fenosoa/Maya_Project/merge_data")

path <- "output_merge/xenium_merge.obj.RData"
temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

merged_xenium <- temp_env$xenium_merge.obj

plot <- VlnPlot(merged_xenium, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ggsave("merge_VlnPlot.png", plot = plot, width = 20, height = 15, dpi = 300)


# merged_xenium@images
# $fov
# Spatial coordinates for 153315 cells and 10029 molecules
#  First 10 molecules: A2ML1, AAMP, AAR2, AARSD1, ABAT, ABCA1
# Default segmentation boundary: centroids 
# Associated assay: Xenium 
# Key: Xenium_ 

# $fov.2
# Spatial coordinates for 146573 cells and 10029 molecules
#  First 10 molecules: A2ML1, AAMP, AAR2, AARSD1, ABAT, ABCA1
# Default segmentation boundary: centroids 
# Associated assay: Xenium 
# Key: fov2_ 

# $fov.3
# Spatial coordinates for 204093 cells and 10029 molecules
#  First 10 molecules: A2ML1, AAMP, AAR2, AARSD1, ABAT, ABCA1
# Default segmentation boundary: centroids 
# Associated assay: Xenium 
# Key: fov3_ 


plot <- ImageDimPlot(merged_xenium, fov = "fov.2", molecules = c("AKT3", "KHK", "CD", "PDCD4"), nmols = 20000)
ggsave("merge_ImageDimPlot_fov_2.png", plot = plot, width = 20, height = 15, dpi = 300)


options(future.globals.maxSize = 1000 * 1024^3)




# SCTransform performs normalization and variance stabilization on the data.
# It replaces the log-normalization step and can help mitigate the effects of technical noise.
print("==================================================================================")

print("# SCTransform performs normalization and variance stabilization on the data.")

xenium.obj <- SCTransform(
  merged_xenium,             # This object containing spatial transcriptomics data
  assay = "Xenium",        # Specify the assay to normalize and scale (e.g., "Xenium")
)

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

print("==================================================================================")
plot4 <- DimPlot(xenium.obj)
ggsave("UMAP_Plot_merged.png", plot = plot4, width = 20, height = 15, dpi = 300)

print("==================================================================================")


print("Let's save the Object for downstreaming task")

umap_xenium_merged.obj <- xenium.obj
save(umap_xenium_merged.obj, file = "umap_xenium_merged.obj.RData")
