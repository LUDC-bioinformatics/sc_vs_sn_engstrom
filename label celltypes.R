library(Seurat)
library(dplyr)
library(ggplot2)
set.seed(1234)


alldata <- readRDS(snakemake@input[[1]])
str(alldata)

# Annotate clusters, including stressed beta
alldata$clu <- alldata$integrated_clusters.cca19202223
summary(alldata$clu)

Idents(object = alldata) <- "clu"
table(Idents(alldata))

new.cluster.ids <- c("beta", "beta","acinar", "beta","alpha","ductal","alpha","alpha","beta","alpha+beta", "stressed beta", "ductal", "delta", "alpha","stressed beta","mesenchymal cells","gamma","endothelial","stressed beta","immune")
stopifnot(length(levels(alldata)) == length(new.cluster.ids))
names(new.cluster.ids) <-levels(alldata)
alldata <- RenameIdents(alldata, new.cluster.ids)
new.cluster.ids

alldata$annotation_stressed_beta <- Idents(alldata)

png(snakemake@output[[1]], width=1400, height=1400)
dimplot <- DimPlot(alldata, reduction = "umap_integrated.cca19202223", group.by = "annotation_stressed_beta", split.by = "type", label = TRUE, pt.size = 0.5,  label.color = "black", cols = c(RColorBrewer::brewer.pal(8,'Set2'), RColorBrewer::brewer.pal(12,'Paired'), RColorBrewer::brewer.pal(12,'Set3') ), label.size = 10, repel = T) + NoAxes()
dimplot +
    ggtitle("Manually annotated celltype") +
    theme(plot.title = element_text(size = 28), 
    legend.text = element_text(size = 24),       # Font size for plot title
    strip.text = element_text(size = 18)     # Font size for facet labels
  )
dev.off()

# Annotate clusters, Stressed beta are now labelled beta


alldata$clu3 <- alldata$integrated_clusters.cca19202223
summary(alldata$clu3)

Idents(object = alldata) <- "clu3"
table(Idents(alldata))


new.cluster.ids <- c("beta", "beta","acinar", "beta","alpha","ductal","alpha","alpha","beta","alpha+beta", "beta", "ductal", "delta", "alpha","beta","mesenchymal cells","gamma","endothelial","beta","immune")


stopifnot(length(levels(alldata)) == length(new.cluster.ids))
names(new.cluster.ids) <- levels(alldata)
alldata <- RenameIdents(alldata, new.cluster.ids)
new.cluster.ids

alldata$annotation_no_stressed_beta <- Idents(alldata)


png(snakemake@output[[2]], width=1400, height=1400)
dimplot <- DimPlot(alldata, reduction = "umap_integrated.cca19202223", group.by = "annotation_no_stressed_beta", split.by = "type", label = TRUE, pt.size = 0.5,  label.color = "black", cols = c(RColorBrewer::brewer.pal(8,'Set2'), RColorBrewer::brewer.pal(12,'Paired'), RColorBrewer::brewer.pal(12,'Set3') ), label.size = 10, repel = T) + NoAxes()
dimplot +
    ggtitle("Manually annotated celltype 3") +
    theme(plot.title = element_text(size = 28), 
    legend.text = element_text(size = 24),       # Font size for plot title
    strip.text = element_text(size = 18)     # Font size for facet labels
  )
dev.off()


saveRDS(alldata, snakemake@output[[3]])
