library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(presto)
set.seed(1234)


alldata <- readRDS(snakemake@input[[1]])
str(alldata)

Idents(object = alldata) <- "integrated_clusters.cca19202223"
table(Idents(alldata))

alldata <- SplitObject(alldata, split.by = "type")
str(alldata)
snRNA <- alldata[["SnRNA"]]
scRNA <- alldata[["ScRNA"]]
str(snRNA)
str(scRNA)
Idents(object = snRNA) <- "integrated_clusters.cca19202223"
table(Idents(snRNA))
Idents(object = scRNA) <- "integrated_clusters.cca19202223"
table(Idents(scRNA))

markers_genes <- FindAllMarkers(scRNA, logfc.threshold = 1, min.pct = 0.25, test.use = "wilcox", assay = "RNA", only.pos = TRUE)
write.table(markers_genes, snakemake@output[[1]])
#only marker genes with FDR < 0.05 were considered, even though they are present in the final list
markers_genes <- FindAllMarkers(snRNA, logfc.threshold = 1, min.pct = 0.25, test.use = "wilcox", assay = "RNA", only.pos = TRUE)
write.table(markers_genes, snakemake@output[[2]])
#only marker genes with FDR < 0.05 were considered, even though they are present in the final list

