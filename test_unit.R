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


