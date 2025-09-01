

library(Seurat)
library(Matrix)

#sc = cells
#mosc = nuclei
s230315354sc <- readRDS(snakemake@input[[1]])

s231031363sc <- readRDS(snakemake@input[[2]])

s210409330sc <- readRDS(snakemake@input[[3]])

s220518345sc <- readRDS(snakemake@input[[4]])

s230403354mosc <- readRDS(snakemake@input[[5]])

s231116363mosc <- readRDS(snakemake@input[[6]])

s230412330mosc <- readRDS(snakemake@input[[7]])


s230215345mosc <- readRDS(snakemake@input[[8]])

#ändra
alldata <- merge(s230315354sc, y= c(s231031363sc, s210409330sc, s220518345sc, s230403354mosc,s231116363mosc,s230412330mosc ,s230215345mosc ), add.cell.ids = c( "s230315354sc", "s231031363sc", "s210409330sc", "s220518345sc", "s230403354mosc","s231116363mosc","s230412330mosc" ,"s230215345mosc"), project = "snsc", merge.data = TRUE)
saveRDS(alldata, snakemake@output[[1]])
