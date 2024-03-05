# scRNAseq-for-EVS-in-hOrganoids
Based on Seurat and Cell2location in R
### For reference
cortex_sc = readRDS(glue::glue("C:\\Users\\leh\\OneDrive - TUNI.fi\\Documents\\Data\\cap_3\\outs\\reference_gut\\Full_obj_log_counts_soupx_v2_DUOdata.RDS"))
library(dplyr)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Seurat)
logcounts(cortex_sc) = assay(cortex_sc, "X")
reference = as.Seurat(cortex_sc, "X")
reference_SCT <- Seurat::SCTransform(reference, assay = NULL) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

  ### For loading h5 file
  library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(hdf5r)

# Loading h5 file in R

sample1 = Read10X_h5(filename = "C://Users//leh//OneDrive - TUNI.fi//Documents//Data//scRNA_2024//sample1//filtered_feature_bc_matrix.h5") 
obj = CreateSeuratObject(counts = sample1)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
library(ggplot2)
DimPlot(obj, group.by = "ident", split.by = 'orig.ident', label = T, label.size = 8, sizes.highlight = 8,  label.color = "black", pt.size = 2) + theme(legend.position = "right", text = element_text(size = 30))

## Checking with the reference
anchors <- FindTransferAnchors(reference = reference_SCT, query = obj, normalization.method = "SCT") # 289 anchors
predictions.assay <- TransferData(anchorset = anchors, refdata = reference_SCT$category, prediction.assay = TRUE,
                                  weight.reduction = obj[["pca"]], dims = 1:25)
obj[["predictions"]] <- predictions.assay
Idents(obj) = predictions.assay
DimPlot(obj, group.by = "predictions", split.by = 'orig.ident', label = T, label.size = 8, sizes.highlight = 8,  label.color = "black", pt.size = 2) + theme(legend.position = "right", text = element_text(size = 30))


