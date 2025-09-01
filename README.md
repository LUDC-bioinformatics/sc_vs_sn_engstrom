# sc_vs_sn_engstrom
Seurat code for manuscript "Evaluating cell-specific gene expression using single-cell and single-nuclei RNA-sequencing data from human pancreatic islets of the same donors"

Seurat scripts:

#SoupX and load samples
SoupX_and_read_Seurat.R

############QC and filtering#######

*QC calculation 
QCcalculation.R

*Filtering cells
Filtering cells.R

*Removing MALAT
removing malat.R

*Merge samples  
merge samples.R

*Dimensionality reduction and integration  
dimensionality reduction and integration.R

*reclustering after removing low quality clusters  
reclustering_and_integration.R

*Manual annotation
marker_gene_analysis_clusters.R

*label celltypes 
label celltypes.R

*Azimuth annotation 
Azimuth.R

* HPAP annotation
celltype annotation HPAP

*Identification of cell type specific marker genes
Marker_gene_analysis_cell_types.R
