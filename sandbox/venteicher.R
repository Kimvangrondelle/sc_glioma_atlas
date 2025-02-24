#file to process extra sample for member of the department -- not included in the research

library(Matrix)

data <- read.table("../../mnt/neuro-genomic-1-ro/single_cell_data/GSE89567_Venteicher/GSE89567_IDH_A_processed_data.txt", header = TRUE)

# to get the sample ids. 

ids <- unique(sub("([^.][_]*[^._]+).*", "\\1", colnames(data)))
ids

sample_id <- list(42, 61, 43, 45, 44, 56, 64, 57, 103, 107)

seurat_list <- list()


for (id in sample_id) {
  sample_data <- data[, grep(id, colnames(data))]
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = paste0("venteicher_", id))
  seurat_obj$sample_id <- id
  seurat_list[[as.character(id)]] <- seurat_obj
}

combined <- merge(seurat_list[[1]], y= seurat_list[2:length(seurat_list)], add.cell.ids=names(seurat_list))

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

combined <- RunPCA(combined, npcs=30)

combined <- RunHarmony(combined, group.by="sample_id", assay.use="RNA")
vent <- combined
vent[["percent.mt"]] <- PercentageFeatureSet(vent, pattern = "^MT")
VlnPlot(vent, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


vent <- subset(x=vent, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA > 200 & nCount_RNA < 22000 & percent.mt < 1.5)
vent <- NormalizeData(vent)
vent <- FindVariableFeatures(vent, selection.method = "vst", nfeatures = 2000)
vent <- ScaleData(vent, features = rownames(vent))
vent <- RunPCA(vent, features = VariableFeatures(object = vent))

DimPlot(vent, reduction = "pca") 
ElbowPlot(vent, ndims = 45)

#12
d <- 12
vent <- FindNeighbors(vent, dims = 1:d)
vent <- FindClusters(vent, resolution = 0.8)
vent <- RunUMAP(vent, dims = 1:d)
DimPlot(vent, reduction = "umap")
saveRDS(vent, file = "output/venteicher/vent_processed.rds")

