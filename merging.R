library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)
plan("multisession", workers = 25)

path1 <- "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_1__20241108__233229 "
path2 <-  "/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_2__20241108__233230"
path3 <-"/home/fenosoa/scratch/Maya_Project/data/output-XETG00325__0029986__Region_3__20241108__233230"
