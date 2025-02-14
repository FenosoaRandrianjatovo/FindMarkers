# Differential Expression Analysis using FindMarkers in Seurat

## Overview
This repository provides a workflow for identifying differentially expressed genes (DEGs) in a spatial transcriptomics dataset processed using the Seurat package in R. The analysis focuses on region of interest  1, which contains 12 clusters, and generates marker gene lists for each cluster using the `FindMarkers` function.

## Data Structure
The analysis is based on an object (`xenium.obj`) containing spatially resolved transcriptomic data. The script systematically applies `FindMarkers` to identify DEGs for each cluster, generating 12 marker files:

```
marker.0.csv
marker.1.csv
...
marker.11.csv
```
Each file contains a list of marker genes along with associated statistical metrics for a given cluster.

## Algorithm: Differential Expression Analysis using FindMarkers

The `FindMarkers` function in Seurat is used to identify genes that are differentially expressed between a given cluster (`ident.1`) and all other clusters in the dataset. The general algorithm follows these steps:

1. **Cluster Identification**
   - Each cell in the dataset is assigned a cluster identity (from 0 to 11 for ROI 1).
   - The `FindMarkers` function is applied to each cluster (`ident.1 = "X"`) to compare gene expression levels in that cluster versus all other clusters.

2. **Statistical Testing**
   - For  `FindMarkers` function we  use the default statistical tests:
     - Wilcoxon Rank Sum Test 

3. **Filtering and Ranking Genes**
   - The function computes several key metrics:
     - **p_val**: The raw p-value for differential expression.
     - **avg_log2FC**: The log2 fold change in gene expression between the cluster and the rest.
     - **pct.1**: The percentage of cells in the target cluster expressing the gene.
     - **pct.2**: The percentage of cells in all other clusters expressing the gene.
     - **p_val_adj**: Adjusted p-value for multiple testing correction (Benjamini-Hochberg method).
   - Genes are ranked based on adjusted p-value and log2 fold change.

4. **Exporting Marker Genes**
   - The identified DEGs are saved as CSV files (`marker.X.csv`), where each file corresponds to one cluster.
   - These files can be used for downstream analyses such as pathway enrichment, cell-type annotation, or spatial mapping of gene expression.

## Example Usage
```r
marker.0 <- FindMarkers(xenium.obj, ident.1 = "0")
write.csv(marker.0, "marker.0.csv")
```
This process is repeated for all clusters from `0` to `11` to generate the complete set of marker gene lists.

## Dependencies
- R (>= 4.0)
- Seurat (>= 4.3)
- arrow
- dplyr

## Conclusion
The resulting marker genes offer insights into cell-type-specific expression patterns within ROI 1, ROI 2 and ROI 3
