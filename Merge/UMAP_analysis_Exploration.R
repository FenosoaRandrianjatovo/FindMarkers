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
setwd("/home/fenosoa/scratch/Maya_Project/FindMarkers/images")

path <- normalizePath("/home/fenosoa/scratch/Maya_Project/merge_data/UMAP_xenium_merge.obj.RData")

temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

xenium.obj <- temp_env$UMAP_xenium_merge.obj

groups <- levels(xenium.obj)

print("Nomber of Cluster")
print(length(groups))

n <- length(groups) 

print("==================================================================================")
plot4 <- DimPlot(xenium.obj)
ggsave("UMAP_merge_Plot.png", plot = plot4, width = 20, height = 15, dpi = 300)

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









