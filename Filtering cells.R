library(Seurat)
library(Matrix)


alldata <- readRDS(snakemake@input[[1]])#one sample is read at a time

##filtering cells
alldata <- subset(x = alldata, 
                         subset= (nCount_RNA >= 500) & 
                           (nFeature_RNA >= 500) & 
                           (log10GenesPerUMI > 0.80) & 
                           (percent_mito < 5) &
                           (percent_ribo < 35))


saveRDS(alldata, snakemake@output[[1]])

