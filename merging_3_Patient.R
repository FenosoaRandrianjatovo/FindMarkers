# Source for this code: https://satijalab.org/seurat/archive/v4.3/merge
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

path1 <- "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_1__20241108__233229"
path2 <-  "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_2__20241108__233230"
path3 <-"/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_3__20241108__233230"

if (!dir.exists("output_merge")) {
  dir.create("output_merge")
}
options(future.globals.maxSize = 520 * 1024^3)

merge_3_XeniumObjects <- function(path1, path2, path3, 
                               output_file = "output_merge/xenium_merge.obj_merge.data_T.RData", 
                               fov = "fov") {
  # Record overall start time
  overall_start <- Sys.time()
  message("Process started at: ", overall_start)

  # Load Xenium object for Patient 1
  xenium.obj1 <- LoadXenium(path1, fov = fov)
  message("# Removing cells with 0 counts from xenium.obj1")
  xenium.obj1 <- subset(xenium.obj1, subset = nCount_Xenium > 0)
  message("# Normalize and scale the data using SCTransform")
  xenium.obj1 <- SCTransform(
    xenium.obj1,             # This object containing spatial transcriptomics data
    assay = "Xenium",        # Specify the assay to normalize and scale (e.g., "Xenium")
  )
  
  # Load Xenium object for Patient 2
  xenium.obj2 <- LoadXenium(path2, fov = fov)
  message("# Removing cells with 0 counts from xenium.obj2")
  xenium.obj2 <- subset(xenium.obj2, subset = nCount_Xenium > 0)
  message("# Normalize and scale the data using SCTransform")
  xenium.obj2 <- SCTransform(
    xenium.obj2,             # This object containing spatial transcriptomics data
    assay = "Xenium",        # Specify the assay to normalize and scale (e.g., "Xenium")
  )
  
  message("# Load Xenium object for Patient 3")
  xenium.obj3 <- LoadXenium(path3, fov = fov)
  message("# Removing cells with 0 counts from xenium.obj3")
  xenium.obj3 <- subset(xenium.obj3, subset = nCount_Xenium > 0)
  message("# Normalize and scale the data using SCTransform")
  xenium.obj3 <- SCTransform(
    xenium.obj3,             # This object containing spatial transcriptomics data
    assay = "Xenium",        # Specify the assay to normalize and scale (e.g., "Xenium")
  )

   # Merge all the objects into one large Seurat object
  message("Merging all Xenium objects...")
  merge_start <- Sys.time()
  
  # Merge the three objects
  xenium_merge.obj <- merge(xenium.obj1, 
                            y = c(xenium.obj2, xenium.obj3), 
                            add.cell.ids = c("Patient1", "Patient2", "Patient3"), 
                            project = "XeniumMergeObject",
                            merge.data = TRUE)
  
  

  merge_end <- Sys.time()
  merge_duration <- merge_end - merge_start
  message("Merging completed in: ", round(as.numeric(merge_duration, units = "secs"), 2), " seconds")

  message("# Save the merged object to file")
  
  save(xenium_merge.obj, file = output_file)

  overall_end <- Sys.time()

  total_duration <- overall_end - overall_start

  message("Total process completed in: ", round(as.numeric(total_duration, units = "secs"), 2), " seconds")
  
  # Return the merged object
  return(xenium_merge.obj)
}

# merged_xenium <- merge_3_XeniumObjects(path1, path2, path3, 
#                                     output_file = "output_merge/xenium_merge.obj.RData", 
#                                     fov = "fov")

merged_xenium <- merge_3_XeniumObjects(path1, path2, path3, 
                                    output_file = "output_merge/xenium_merge.obj_merge.data_T.RData", 
                                    fov = "fov")

mergeAllXeniumObjects <- function(paths, fov = "fov", output_file = "output_merge/xenium_merge_all.obj.RData") {

  # Record overall start time
  overall_start <- Sys.time()
  # Number of patients (number of file paths)
  num_patients <- length(paths)
  
  # Initialize a list to store each loaded Xenium object
  xenium_list <- vector("list", num_patients)
  
  # Create unique patient IDs to prefix cell names (e.g., "Patient1", "Patient2", ...)
  patient_ids <- paste0("Patient", seq_len(num_patients))
  
  message("# Loop through each path, load the Xenium object, and remove cells with 0 counts")
  
  for (i in seq_along(paths)) {
    message("Loading data for ", patient_ids[i], " from: ", paths[i])
    xenium_obj <- LoadXenium(paths[i], fov = fov)
    
    message("Removing cells with 0 counts for ", patient_ids[i])
    xenium_obj <- subset(xenium_obj, subset = nCount_Xenium > 0)
    
    # Store the processed object in the list
    xenium_list[[i]] <- xenium_obj
  }
  
  # Merge all the objects into one large Seurat object
  message("Merging all Xenium objects...")
  merged_obj <- merge(xenium_list[[1]], y = xenium_list[-1], add.cell.ids = patient_ids)
  
  # Save the merged object to the specified file path
  message("Saving merged object to: ", output_file)
  save(merged_obj, file = output_file)
  
  overall_end <- Sys.time()

  total_duration <- overall_end - overall_start

  message("Total process completed in: ", round(as.numeric(total_duration, units = "secs"), 2), " seconds")
  
  return(merged_obj)
}


# # Example: paths vector containing 300 file paths
# paths <- c("path/to/patient1/data",
#            "path/to/patient2/data",
#            # ...
#            "path/to/patient300/data")

# # Merge all Xenium objects and save the result
# merged_xenium <- mergeXeniumObjects(paths, fov = "fov", 
#                                     output_file = "mypath/xenium_merge.obj.RData")


# > dim(LayerData(merged_xenium))
# Warning: data layer is not found and counts layer is used
# [1]   5001 153315
# Warning message:
# In LayerData.Assay5(object = object[[assay]], layer = layer, ...) :
#   multiple layers are identified by counts.1 counts.2 counts.3
#  only the first layer is used
# > dim(merged_xenium@assays$Xenium["counts.1"])
# [1]   5001 153315
# > dim(merged_xenium@assays$Xenium["counts.2"])
# [1]   5001 146573
# > dim(merged_xenium@assays$Xenium["counts.3"])
# [1]   5001 204093





