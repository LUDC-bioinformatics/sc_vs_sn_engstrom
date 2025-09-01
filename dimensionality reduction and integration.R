library(Seurat)
set.seed(1234)
sessionInfo()

alldata <- readRDS(snakemake@input[[1]])
str(alldata)
alldata[["RNA"]] <- split(alldata[["RNA"]], f = alldata$orig.ident)


# run standard anlaysis workflow
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata)
alldata <- ScaleData(alldata, vars.to.regress = c("nFeature_RNA", "percent_mito"))
alldata <- RunPCA(alldata)

alldata <- FindNeighbors(alldata, dims = 1:30, reduction = "pca")
alldata <- FindClusters(alldata, resolution = 1, cluster.name = "unintegrated_clusters")


alldata <- RunUMAP(alldata, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

png(snakemake@output[[1]], width=2000, height=1000)
DimPlot(alldata, reduction = "umap.unintegrated", group.by = c("orig.ident", "unintegrated_clusters", "type", "predicted.annotation.l1"))
dev.off()

alldata <- IntegrateLayers(object = alldata, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE)



# re-join layers after integration
alldata[["RNA"]] <- JoinLayers(alldata[["RNA"]])

alldata <- FindNeighbors(alldata, reduction = "integrated.cca", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 1, cluster.name = "integrated_clusters.cca")

alldata <- RunUMAP(alldata, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap_integrated.cca")

png(snakemake@output[[2]], width=2000, height=1000)
DimPlot(alldata, reduction = "umap_integrated.cca", group.by = c("orig.ident", "integrated_clusters.cca", "predicted.annotation.l1"))
dev.off()

saveRDS(alldata, snakemake@output[[3]])



