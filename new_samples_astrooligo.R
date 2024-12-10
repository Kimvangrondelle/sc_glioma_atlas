astros <- readRDS("../../mnt/neuro-genomic-1-ro/Glimmunology/LGG_project/data_analysis/GLASS/3PR.RDS")

DimPlot(astros, reduction = "umap")

saveRDS(astros, "output/hijfte/3pr.rds")




#11136 ----
diaz11136 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119535_SF11136_LGG", gene.column = 1)
diaz11136 <- CreateSeuratObject(counts = diaz11136, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz11136[["percent.mt"]] <- PercentageFeatureSet(diaz11136, pattern = "^MT")
VlnPlot(diaz11136, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz11136 <- subset(x=diaz11136, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA > 200 & nCount_RNA < 30000 & percent.mt < 10)
diaz11136 <- NormalizeData(diaz11136)
diaz11136 <- FindVariableFeatures(diaz11136, selection.method = "vst", nfeatures = 2000)
diaz11136 <- ScaleData(diaz11136, features = rownames(diaz11136))
diaz11136 <- RunPCA(diaz11136, features = VariableFeatures(object = diaz11136))

DimPlot(diaz11136, reduction = "pca") 
ElbowPlot(diaz11136, ndims = 45)

d_11136 <- 15
diaz11136 <- FindNeighbors(diaz11136, dims = 1:d_11136)
diaz11136 <- FindClusters(diaz11136, resolution = 0.8)
diaz11136 <- RunUMAP(diaz11136, dims = 1:d_11136)
DimPlot(diaz11136, reduction = "umap")

saveRDS(diaz11136, "output/diaz_astrooligo/diaz11136_astro.rds")

markers_11136 <- FindAllMarkers(diaz11136, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_11136 %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)


diaz12017 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119536_SF12017_LGG", gene.column = 1)
diaz12017 <- CreateSeuratObject(counts = diaz12017, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz12017[["percent.mt"]] <- PercentageFeatureSet(diaz12017, pattern = "^MT")
VlnPlot(diaz12017, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz12017 <- subset(x=diaz12017, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & nCount_RNA > 200 & nCount_RNA < 37500 & percent.mt < 7.5)
diaz12017 <- NormalizeData(diaz12017)
diaz12017 <- FindVariableFeatures(diaz12017, selection.method = "vst", nfeatures = 2000)
diaz12017 <- ScaleData(diaz12017, features = rownames(diaz12017))
diaz12017 <- RunPCA(diaz12017, features = VariableFeatures(object = diaz12017))

DimPlot(diaz12017, reduction = "pca") 
ElbowPlot(diaz12017, ndims = 45)

d_12017 <- 12
diaz12017 <- FindNeighbors(diaz12017, dims = 1:d_12017)
diaz12017 <- FindClusters(diaz12017, resolution = 0.8)
diaz12017 <- RunUMAP(diaz12017, dims = 1:d_12017)
DimPlot(diaz12017, reduction = "umap")

saveRDS(diaz12017, "output/diaz_astrooligo/diaz12017_astro.rds")



diaz11612 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119539_SF11612_ODG", gene.column = 1)
diaz11612 <- CreateSeuratObject(counts = diaz11612, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz11612[["percent.mt"]] <- PercentageFeatureSet(diaz11612, pattern = "^MT")
VlnPlot(diaz11612, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz11612 <- subset(x=diaz11612, subset = nFeature_RNA > 200 & nFeature_RNA < 2700 & nCount_RNA > 200 & nCount_RNA < 7500 & percent.mt < 27)
diaz11612 <- NormalizeData(diaz11612)
diaz11612 <- FindVariableFeatures(diaz11612, selection.method = "vst", nfeatures = 2000)
diaz11612 <- ScaleData(diaz11612, features = rownames(diaz11612))
diaz11612 <- RunPCA(diaz11612, features = VariableFeatures(object = diaz11612))

DimPlot(diaz11612, reduction = "pca") 
ElbowPlot(diaz11612, ndims = 45)

d_11612 <- 11
diaz11612 <- FindNeighbors(diaz11612, dims = 1:d_11612)
diaz11612 <- FindClusters(diaz11612, resolution = 0.8)
diaz11612 <- RunUMAP(diaz11612, dims = 1:d_11612)
DimPlot(diaz11612, reduction = "umap")

saveRDS(diaz11612, "output/diaz_astrooligo/diaz11612_oligo.rds")



diaz11949 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119538_SF11949_ODG", gene.column = 1)
diaz11949 <- CreateSeuratObject(counts = diaz11949, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz11949[["percent.mt"]] <- PercentageFeatureSet(diaz11949, pattern = "^MT")
VlnPlot(diaz11949, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz11949 <- subset(x=diaz11949, subset = nFeature_RNA > 200 & nFeature_RNA < 2600 & nCount_RNA > 200 & nCount_RNA < 29000 & percent.mt < 29)
diaz11949 <- NormalizeData(diaz11949)
diaz11949 <- FindVariableFeatures(diaz11949, selection.method = "vst", nfeatures = 2000)
diaz11949 <- ScaleData(diaz11949, features = rownames(diaz11949))
diaz11949 <- RunPCA(diaz11949, features = VariableFeatures(object = diaz11949))

DimPlot(diaz11949, reduction = "pca") 
ElbowPlot(diaz11949, ndims = 45)

d_11949 <- 9
diaz11949 <- FindNeighbors(diaz11949, dims = 1:d_11949)
diaz11949 <- FindClusters(diaz11949, resolution = 0.8)
diaz11949 <- RunUMAP(diaz11949, dims = 1:d_11949)
DimPlot(diaz11949, reduction = "umap")

saveRDS(diaz11949, "output/diaz_astrooligo/diaz11949_oligo.rds")