library(Seurat)
library(Matrix)


# Load and check dimensions
alldata <- readRDS(snakemake@input[[1]])#one sample is read at a time

# Filter MALAT1
alldata <- alldata[ !Features(alldata) %in% "MALAT1", ]
#
saveRDS(alldata, snakemake@output[[1]])



