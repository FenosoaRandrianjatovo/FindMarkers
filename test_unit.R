library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)
plan("multisession", workers = 25)

# To Load te file for another downstreaming task
path <- "/home/fenosoa/scratch/Maya_Project/merge_data/output_merge/xenium_merge.obj.RData"
temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

xenium.obj <- temp_env$xenium_merge.obj

counts1 <- as(merged_xenium@assays$Xenium@layers$counts.1, "dgCMatrix")
counts2 <- as(merged_xenium@assays$Xenium@layers$counts.2, "dgCMatrix")
counts3 <- as(merged_xenium@assays$Xenium@layers$counts.3, "dgCMatrix")

# Combine them column-wise
merged_counts <- cbind(counts1, counts2, counts3)

# Expected: 5001 x 503981
print(dim(merged_counts))

# Remove the old separate count layers to avoid confusion
xenium.obj@assays$Xenium@layers$counts.1 <- NULL
xenium.obj@assays$Xenium@layers$counts.2 <- NULL
xenium.obj@assays$Xenium@layers$counts.3 <- NULL

xenium.obj@assays$Xenium@layers <- list("counts.1" = merged_counts)

print(dim(xenium.obj@assays$Xenium@layers$counts.1)) # Expected: 5001 x 503981

# Normalize and scale the data using SCTransform
xenium.obj_test <- SCTransform(
  xenium.obj,             # This object containing spatial transcriptomics data
  assay = "Xenium",        # Specify the assay to normalize and scale (e.g., "Xenium")
)















