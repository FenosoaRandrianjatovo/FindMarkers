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

path1 <- "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_1__20241108__233229 "
path2 <-  "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_2__20241108__233230"
path3 <-"/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_3__20241108__233230"

if (!dir.exists("output_merge")) {
  dir.create("output_merge")
}

merge_3_XeniumObjects <- function(path1, path2, path3, 
                               output_file = "output_merge/xenium_merge.obj.RData", 
                               fov = "fov") {
  # Load Xenium object for Patient 1
  xenium.obj1 <- LoadXenium(path1, fov = fov)
  message("# Removing cells with 0 counts from xenium.obj1")
  xenium.obj1 <- subset(xenium.obj1, subset = nCount_Xenium > 0)
  
  # Load Xenium object for Patient 2
  xenium.obj2 <- LoadXenium(path2, fov = fov)
  message("# Removing cells with 0 counts from xenium.obj2")
  xenium.obj2 <- subset(xenium.obj2, subset = nCount_Xenium > 0)
  
  # Load Xenium object for Patient 3
  xenium.obj3 <- LoadXenium(path3, fov = fov)
  message("# Removing cells with 0 counts from xenium.obj3")
  xenium.obj3 <- subset(xenium.obj3, subset = nCount_Xenium > 0)
  
  # Merge the three objects
  xenium_merge.obj <- merge(xenium.obj1, 
                            y = list(xenium.obj2, xenium.obj3), 
                            add.cell.ids = c("Patient1", "Patient2", "Patient3"))
  
  message("# Save the merged object to file")
  save(xenium_merge.obj, file = output_file)
  
  # Return the merged object
  return(xenium_merge.obj)
}

