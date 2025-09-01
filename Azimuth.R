library(Seurat)
library(Azimuth)
library(patchwork)
set.seed(1234)

sessionInfo()
alldata <- readRDS(snakemake@input[[1]])
print("1)")

print("3)")

alldata <- RunAzimuth(alldata, reference = "pancreasref")
print("4)")

saveRDS(alldata, snakemake@output[[1]])
