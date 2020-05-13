library(dplyr)
library(Seurat)
library(patchwork)
#
#
#### import scRNA-seq data of healthy controls ####
healthy.data <- Read10X(data.dir = "cellranger_output_for_healthy_controls/")
healthy <- CreateSeuratObject(counts = healthy.data, project = "healthy", min.cells = 3, min.features = 200)
healthy[["percent.mt"]] <- PercentageFeatureSet(healthy, pattern = "^MT-")
healthy_filtered <- subset(healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
healthy_doublet <- read.table('healthy_controls_doublets.txt', sep=',', check.names=F, row.names=1)
healthy_doublet <- healthy_doublet[healthy_doublet$predicted_doublets=='True',]
healthy_cells <- setdiff(colnames(healthy_filtered), rownames(healthy_doublet[1]))
healthy_filtered <- SubsetData(healthy_filtered, cells=healthy_cells)
healthy_filtered@meta.data['batch'] <- healthy_filtered@meta.data['orig.ident']
#
#
#### import scRNA-seq data of COVID-19 patients ####
ncov.data <- Read10X(data.dir = "cellranger_output_for_COVID19_patients/")
ncov <- CreateSeuratObject(counts = ncov.data, project = "covid19", min.cells = 3, min.features = 200)
ncov[["percent.mt"]] <- PercentageFeatureSet(ncov, pattern = "^MT-")
ncov_filtered <- subset(ncov, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
cellinfo <- read.table("COVID19_cells_info.csv", sep='\t', header=TRUE, quote="", check.names=FALSE, row.names=1)
ncov_filtered <- AddMetaData(ncov_filtered, metadata = cellinfo[colnames(ncov_filtered),], col.name=c('batch'))
ncov_doublet <- read.table('COVID19_patients_doublets.txt', sep=',', check.names=F, row.names=1)
ncov_doublet <- ncov_doublet[ncov_doublet$predicted_doublets=='True',]
ncov_cells <- setdiff(colnames(ncov_filtered), rownames(ncov_doublet[1]))
ncov_filtered <- SubsetData(ncov_filtered, cells=ncov_cells)
#
#
#### integrated cells from healthy controls and COVID-19 patients ####
overlap_genes <- intersect(rownames(ncov_filtered), rownames(pbmc_filtered))
pbmc.list <- c(pbmc_filtered, ncov_filtered)
names(pbmc.list) = c('healthy', 'covid19')
pbmc.anchors <- FindIntegrationAnchors(object.list=pbmc.list, dims = 1:40)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:40, features.to.integrate=overlap_genes)
#
#
#### normalization, clustering, UMAP, and plotting marker genes ####
pbmc.integrated <- ScaleData(pbmc.integrated, features = rownames(pbmc.integrated))
pbmc.integrated <- RunPCA(pbmc.integrated, npcs=50, verbose=F)
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:50)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.3)
pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:50)
DimPlot(pbmc.integrated, reduction = "umap", label=TRUE)
FeaturePlot(pbmc.integrated, features=c('PTPRC', 'CD14', 'FCGR3A', 'CD3D', 'CD4', 'IL7R', 'CCR7', 'CD8A', 'PRDM1', 'MKI67', 
                                        'TRGC1', 'NKG7', 'CD79A', 'CD38', 'CD1C', 'CLEC4C', 'PPBP', 'CD34'), 
            slot='data', min.cutoff=0, max.cutoff='q90')
#
#
#### save data ####
save(pbmc.integrated, file='integrated.allgenes.RData')
write.table(pbmc.integrated[['RNA']]@data, file='RNA.data.csv',sep=',', quote=F)
write.table(pbmc.integrated[['integrated']]@data, file='integrated.data.csv', sep=',', quote=F)
write.table(pbmc.integrated@reductions$umap@cell.embeddings, file='umap.csv', sep=',', quote=F)
write.table(pbmc.integrated@meta.data, file='meta_data.csv', sep=',')
#
#