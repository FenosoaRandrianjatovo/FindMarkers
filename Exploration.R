# module load StdEnv/2023
# module spider arrow/18.1.0
# module load arrow/18.1.0
# module spider gdal
# module spider proj
# module spider geos
# module spider udunits
# module load  gcc/12.3  r/4.4.0


library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)
plan("multisession", workers = 25)

print("Set the working directory")
setwd("/home/fenosoa/scratch/Maya_Project/pipeline_Region/out3")

path <- "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_3__20241108__233230"


xenium.obj <- LoadXenium(path, fov = "fov")

print("# remove cells with 0 counts")
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

print("==================================================================================")
genes <- rownames(xenium.obj)

print("Save Genes as a text file")
# Save as a text file
write.table(genes, file = "genes_region3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


options(future.globals.maxSize = 480 * 1024^3)  # Set the limit to 400 GB

print("==================================================================================")
print("SCTransform begin")

# Normalize and scale the data using SCTransform
xenium.obj <- SCTransform(
  xenium.obj,             # This object containing spatial transcriptomics data
  assay = "Xenium",        # Specify the assay to normalize and scale (e.g., "Xenium")
)

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

print("==================================================================================")
plot4 <- DimPlot(xenium.obj)
ggsave("/home/fenosoa/scratch/Maya_Project/pipeline_Region/out3/UMAP_Plot.png", plot = plot4, width = 20, height = 15, dpi = 300)

print("==================================================================================")


print("Let's save the Object for downstreaming task")

xenium3.obj <- xenium.obj
save(xenium3.obj, file = "/home/fenosoa/scratch/Maya_Project/pipeline_Region/out3/UMAP_xenium3.obj.RData")

# To Load te file for another downstreaming task
# path <- "/home/fenosoa/scratch/Maya_Project/pipeline_Region/out3/UMAP_xenium3.obj.RData"
# temp_env <- new.env()
# print("xenium.obj :   Loading begin")
# load(path, envir = temp_env)

# xenium.obj <- temp_env$xenium3.obj

plot4 <- DimPlot(xenium.obj)
ggsave("UMAP_R3_Plot.png", plot = plot4, width = 20, height = 15, dpi = 300)


plot2 <- ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)
ggsave("ImageDimPlot_R3.png", plot = plot2, width = 20, height = 15, dpi = 300)

ggsave("ImageDimPlot_VS_UMAP_R3.png", plot = plot2+ plot4, width = 30, height = 20, dpi = 300)


groups <- levels(xenium.obj)

print("Nomber of Cluster")
print(length(groups))

n <- length(groups) 
# We have 15 Cluster for region 3
# Create the "ident" folder if it does not exist
if (!dir.exists("ident")) {
  dir.create("ident")
}

# Loop through ident values from 0 to 15 
print("Begin to save ImageDimPlot for each ident")
for (ident in 0:(n-1)) {
  # Generate the plot
  p <- ImageDimPlot(xenium.obj, fov = "fov", cols = "red", 
                    cells = WhichCells(xenium.obj, idents = ident))
  
  # Define the filename dynamically
  filename <- paste0("ident/ImageDimPlot_for_ident", ident, ".png")
  

  ggsave(filename, plot = p, width = 30, height = 20, dpi = 300)
  print(paste("Saved plot:", filename))
}

print("Begin to save ImageFeaturePlot for each ident")
p1 <- ImageFeaturePlot(xenium.obj, features = "ATP1A1")
p2 <- ImageDimPlot(xenium.obj, molecules = "ATP1A1", nmols = 10000, alpha = 0.3, mols.cols = "red")

ggsave("ImageDimPlot_for_melecule_ATP1A1.png", plot = p1 + p2, width = 30, height = 20, dpi = 300)

# To run FindMarkers we need presto
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# Let's find all the Markers of the 15 clusters.


# Create the "markers" folder if it does not exist
if (!dir.exists("markers")) {
  dir.create("markers")
}
print("Begin to save markers gene for each ident")
# Initialize an empty list to store marker genes
markers_list <- list()
marker.0 <- FindMarkers(xenium.obj, ident.1 = "0")
# Loop through ident values from 0 to 11
for (ident in 0:(n-1)) {
  # Find marker genes for the current ident
  markers <- FindMarkers(xenium.obj, ident.1 = as.character(ident))
  
  # Store the result in a dynamically named variable
  assign(paste0("markers.", ident), markers, envir = .GlobalEnv)
  
  # Also store it in a list for easy future reference
  markers_list[[as.character(ident)]] <- markers
  
  # Define the filename dynamically
  filename <- paste0("markers/markers.", ident, ".csv")
  
  # Save the markers data as a CSV file
  write.csv(markers, file = filename, row.names = TRUE)
  print(paste("Saved markers for ident", ident, "to", filename))
}
# Create the "markers" folder if it does not exist
if (!dir.exists("Feature")) {
  dir.create("Feature")
}


print("Begin to save FeaturePlot and ImageFeaturePlot for each ident")
# Define marker indices (0 to 11)
marker_indices <- 0:(n-1)

# Loop through marker objects from markers.0 to markers.11
for (i in 0:(n-1)) {
  # Construct the variable name dynamically
  marker_var <- get(paste0("markers.", i), envir = .GlobalEnv)
  
  # Ensure the object exists and has at least 6 rows
  if (!is.null(marker_var) && nrow(marker_var) >= 6) {
    
    # Loop through the first six selected genes
    for (j in 1:6) {
      selected_gene <- rownames(marker_var)[j]
      
      # Generate FeaturePlot and ImageFeaturePlot for the current gene
      feature_plot <- FeaturePlot(xenium.obj, features = selected_gene)
      image_feature_plot <- ImageFeaturePlot(xenium.obj, features = selected_gene)
      
      # Define the filename dynamically
      filename <- paste0("Feature/FeaturePlot_ImageFeaturePlot_", selected_gene, "_ident", i, ".png")
      
      # Save the ImageFeaturePlot
      ggsave(filename, plot = image_feature_plot + feature_plot, width = 50, height = 40, dpi = 300, limitsize = FALSE)
      
      print(paste("Saved plot:", filename))
    }
  } else {
    print(paste("Skipping markers.", i, " - Not enough genes or NULL"))
  }
}

# Plot a marker gene’s expression: VlnPlot
if (!dir.exists("VlnPLot")) {
  dir.create("VlnPLot")
}
# Loop through markers.0 to markers.11
print("Begin to save VlnPlot for each ident with the first 6 genes")
for (i in 0:(n-1)) {
  marker_var <- get(paste0("markers.", i))  # Get marker dataset dynamically
  
  if (!is.null(marker_var) && nrow(marker_var) >= 6) {
    for (j in 1:6) {  # Loop through the first 6 genes
      gene_name <- rownames(marker_var)[j]

      # Generate violin plot
      te <- VlnPlot(xenium.obj, features = gene_name, pt.size = 0.1)

      # Define filename dynamically
      filename <- paste0("VlnPLot/markers.", i, "_gene_", gene_name, "_VlnPlot.png")

      # Save the plot
      ggsave(filename, plot = te, width = 20, height = 15, dpi = 300)

      print(paste("Saved plot:", filename))
    }
  } else {
    print(paste("Skipping markers.", i, " - Not enough genes or NULL"))
  }
}
if (!dir.exists("VlnPlot_Feature_Count")) {
  dir.create("VlnPlot_Feature_Count")
}
print("Begin to save VlnPlot for nFeature_Xenium and nCount_Xenium")

test <- VlnPlot(xenium.obj, features = "nFeature_Xenium", ncol = 1, pt.size = 0.1)
ggsave("VlnPlot_Feature_Count/VlnPlot_nFeature_Xenium_after_Clustering.png", plot = test, width = 20, height = 15, dpi = 300)

test1 <- VlnPlot(xenium.obj, features = "nCount_Xenium", ncol = 1, pt.size = 0.1)
ggsave("VlnPlot_Feature_Count/VlnPlot_nCount_Xenium_after_Clustering.png", plot = test1, width = 20, height = 15, dpi = 300)


before <- VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), group.by = "orig.ident")
ggsave("VlnPlot_Feature_Count/VlnPlot_nFeature_Xenium_nCount_Xenium_before_Clustering.png", plot = before, width = 20, height = 15, dpi = 300)





# Problem with Cropping after doin the following command 
# install.packages("sf", dependencies = TRUE)
# DONE (spatstat)
# ERROR: dependencies ‘sf’, ‘units’ are not available for package ‘stars’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/stars’
# ERROR: dependencies ‘leaflet’, ‘raster’, ‘sf’ are not available for package ‘leafem’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leafem’
# ERROR: dependencies ‘leaflet’, ‘sf’ are not available for package ‘leafgl’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leafgl’
# ERROR: dependency ‘leaflet’ is not available for package ‘leaflegend’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leaflegend’
# ERROR: dependency ‘leaflet’ is not available for package ‘leafsync’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leafsync’
# ERROR: dependencies ‘sf’, ‘lwgeom’, ‘stars’, ‘units’ are not available for package ‘tmaptools’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/tmaptools’
# ERROR: dependencies ‘leafem’, ‘leaflet’, ‘leafpop’, ‘raster’, ‘satellite’, ‘sf’ are not available for package ‘mapview’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/mapview’
# ERROR: dependencies ‘leafem’, ‘leafgl’, ‘leaflegend’, ‘leaflet’, ‘leafsync’, ‘sf’, ‘stars’, ‘tmaptools’, ‘units’ are not available for package ‘tmap’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/tmap’

# The downloaded source packages are in
# 	‘/tmp/RtmpDLzjBp/downloaded_packages’
# There were 18 warnings (use warnings() to see them)


# To run crop we need The package "sf" is required to overlay spatial information on the image.
# create a Crop
# cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1750, 3000), y = c(3750, 5250), coords = "plot")
# # set a new field of view (fov)
# xenium.obj[["tumor"]] <- cropped.coords


plot2 <- SpatialFeaturePlot(xenium.obj, features = "nCount_Spatial") + theme(legend.position = "right")
