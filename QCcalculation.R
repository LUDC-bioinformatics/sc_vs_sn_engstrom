library(Seurat)
library(Matrix)


alldata <- readRDS(snakemake@input[[1]])
str(alldata)
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")


alldata$log10GenesPerUMI <- log10(alldata$nFeature_RNA) / log10(alldata$nCount_RNA)
str(alldata)
saveRDS(alldata, snakemake@output[[1]])

