library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(presto)
set.seed(1234)


alldata <- readRDS(snakemake@input[[1]])
str(alldata)

pancreas.query <- SplitObject(alldata, split.by = "type")
str(pancreas.query)
scRNA <- pancreas.query[["ScRNA"]]
snRNA <- pancreas.query[["SnRNA"]]


#############Analysis including stressed beta
Idents(object = scRNA) <- "annotation_stressed_beta"
table(Idents(scRNA))

Idents(object = snRNA) <- "annotation_stressed_beta"
table(Idents(snRNA))

    # Compute differentiall expression scRNA-seq, incl stressed beta
markers_genes <- FindAllMarkers(scRNA, logfc.threshold = 1, min.pct = 0.5, test.use = "wilcox", assay = "RNA", only.pos = TRUE)
write.table(markers_genes, snakemake@output[[1]])

    # Compute differentiall expression snRNA-seq, incl stressed beta
markers_genes <- FindAllMarkers(snRNA, logfc.threshold = 1, min.pct = 0.5, test.use = "wilcox", assay = "RNA", only.pos = TRUE)
write.table(markers_genes, snakemake@output[[2]])



Idents(object = scRNA) <- "annotation_no_stressed_beta"
table(Idents(scRNA))

Idents(object = snRNA) <- "annotation_no_stressed_beta"
table(Idents(snRNA))

    # Compute differentiall expression scRNA-seq, dataset without stressed beta
markers_genes <- FindAllMarkers(scRNA, logfc.threshold = 1, min.pct = 0.5, test.use = "wilcox", assay = "RNA", only.pos = TRUE)
write.table(markers_genes, snakemake@output[[3]])
    # Compute differentiall expression snRNA-seq, dataset without stressed beta
    # Compute differentiall expression
markers_genes <- FindAllMarkers(snRNA, logfc.threshold = 1, min.pct = 0.5, test.use = "wilcox", assay = "RNA", only.pos = TRUE)
write.table(markers_genes, snakemake@output[[4]])
