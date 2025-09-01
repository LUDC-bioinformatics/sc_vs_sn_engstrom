library(SoupX)
library(Seurat)
library(Matrix)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
set.seed(1234)

#Per sample basis, each sample has a folder in Cellranger
filt.matrix <- Read10X_h5(snakemake@input[[1]],use.names = T)
raw.matrix  <- Read10X_h5(snakemake@input[[2]],use.names = T)

srat  <- CreateSeuratObject(counts = filt.matrix)
str(srat)

sample_name <- basename(dirname(snakemake@input[[1]])) 
srat$orig.ident = sample_name

#SoupX
sc  <- SoupChannel(raw.matrix, filt.matrix)
sc
srat    <- NormalizeData(srat, normalization.method = "LogNormalize")
srat    <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
srat    <- ScaleData(srat)
srat    <- RunPCA(srat, features = VariableFeatures(object = srat))
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, resolution = 1.0, verbose = T)

str(srat)
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings

str(meta)
meta$seurat_clusters

rownames(meta)

str(meta)
umap    <- srat@reductions$umap@cell.embeddings
sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)
head(meta)
png(snakemake@output[[1]], width=1400, height=1400)
sc  <- autoEstCont(sc)
dev.off()

png(snakemake@output[[2]], width=1400, height=1400)
plotMarkerDistribution(sc)
dev.off()

a <- head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20)
write.table(a, snakemake@output[[3]])

str(sc)
b <- sc$fit$rhoEst
write.table(b, snakemake@output[[4]])

adj.matrix  <- adjustCounts(sc, roundToInt = T)

alldata = CreateSeuratObject(adj.matrix, min.cells = 10, min.features = 200, project = sample_name)
alldata$orig.ident <- sample_name

b <- dim(alldata)
write.table(b, snakemake@output[[5]])

saveRDS(alldata, snakemake@output[[6]])
