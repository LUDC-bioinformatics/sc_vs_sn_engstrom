library(Seurat)
set.seed(1234)
sessionInfo()

alldata <- readRDS(snakemake@input[[1]])
str(alldata)

#remove low quality clusters
alldata[["RNA"]] <- JoinLayers(alldata[["RNA"]])
Idents(object = alldata) <- "integrated_clusters.cca"
alldata <- subset(alldata, idents = c("19", "20", "22", "23"), invert = TRUE)

#split per sample (n=8)
alldata[["RNA"]] <- split(alldata[["RNA"]], f = alldata$orig.ident)


# run standard anlaysis workflow on nonintegrated samples
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata)
alldata <- ScaleData(alldata, vars.to.regress = c("nFeature_RNA", "percent_mito"))
alldata <- RunPCA(alldata, reduction.name = "pca19202223")
alldata <- FindNeighbors(alldata, reduction = "pca19202223", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 1, cluster.name = "unintegrated_clusters19202223")


alldata <- RunUMAP(alldata, dims = 1:30, reduction = "pca19202223", reduction.name = "umap.unintegrated19202223")

png(snakemake@output[[1]], width=2000, height=1000)
DimPlot(alldata, reduction = "umap.unintegrated19202223", group.by = c("orig.ident", "unintegrated_clusters19202223", "type", "predicted.annotation.l1"))
dev.off()

#integrate samples again
alldata <- IntegrateLayers(object = alldata, method = CCAIntegration, orig.reduction = "pca19202223", new.reduction = "integrated.cca19202223",
    verbose = FALSE)

# re-join layers after integration and make a new UMAP
alldata[["RNA"]] <- JoinLayers(alldata[["RNA"]])
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000)
alldata <- ScaleData(alldata, vars.to.regress = c("nFeature_RNA", "percent_mito"), verbose = FALSE)
alldata <- RunPCA(alldata, npcs = 30, verbose = FALSE)
alldata <- FindNeighbors(alldata, reduction = "integrated.cca19202223", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 1, cluster.name = "integrated_clusters.cca19202223")

alldata <- RunUMAP(alldata, dims = 1:30, reduction = "integrated.cca19202223", reduction.name = "umap_integrated.cca19202223")

png(snakemake@output[[2]], width=2000, height=1000)
DimPlot(alldata, reduction = "umap_integrated.cca19202223", group.by = c("orig.ident", "integrated_clusters.cca19202223", "predicted.annotation.l1"))
dev.off()

Layers(alldata)
DefaultAssay(alldata) <- "RNA"

png(snakemake@output[[3]], width=1000, height=1500)
FeaturePlot(alldata, reduction = "umap_integrated.cca19202223", features = c("INS", "GCG","nCount_RNA","nFeature_RNA","percent_mito","percent_ribo"))
dev.off()


saveRDS(alldata, snakemake@output[[4]])


a <- dim(alldata)
write.table(a, snakemake@output[[5]])


