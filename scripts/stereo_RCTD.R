sto <- readRDS("sto/sto_bin50.rds")
ref <- dior::read_h5ad(file='adata_ovary_combined_fullprocessed_annotated.h5ad', target.object = 'seurat') #https://doi.org/10.6084/m9.figshare.c.6930214

Idents(ref)<-"Celltypes"
counts <- ref@assays$RNA@counts
counts <-floor(counts)

cluster <- as.factor(adata$Celltypes)
names(cluster) <- colnames(adata)
nUMI <- adata$n_counts
names(nUMI) <- colnames(adata)
reference <- Reference(counts, cluster, nUMI)

sto <-subset(sto,nFeature_Spatial>1000&percent.mito<5)
coords <- data.frame(x=sto@meta.data$x,y=sto@meta.data$y)
coords[is.na(colnames(coords))] <- NULL
rownames(coords)<-rownames(sto@meta.data)
library(spacexr,lib = "../rpackage/")
counts <- sto@assays$Spatial@counts
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
sto <- AddMetaData(sto, metadata = RCTD@results$results_df)
DimPlot(sto,reduction = "spatial",group.by = "first_type")
