
library(Seurat)

library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(1234)

pancreas.query <- readRDS(snakemake@input[[1]]) #our data set
hpap           <- readRDS(snakemake@input[[2]]) #HPAP nondiabetics dataset

#Note the original HPAP variable @ meta.data$Cell Type were relabelled celltype

pancreas.query[["RNA"]] <- JoinLayers(pancreas.query[["RNA"]])

DefaultAssay(hpap) <- "RNA"

pancreas.anchors <- FindTransferAnchors(reference = hpap, query = pancreas.query,
     reference.assay = "RNA", query.assay = "RNA", dims = 1:20, reduction = "pcaproject")

predictions <- TransferData(anchorset = pancreas.anchors, refdata = hpap$celltype,
    dims = 1:20)

pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)


pancreas.query@active.assay = "RNA"

saveRDS(pancreas.query, snakemake@output[[1]])