# 10022 ----
diaz10022 <- Read10X(data.dir = "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119521_SF10022_GBM")
diaz10022 <- CreateSeuratObject(counts = diaz10022, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz10022[["percent.mt"]] <- PercentageFeatureSet(diaz10022, pattern = "^MT")
VlnPlot(diaz10022, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz10022 <- subset(x=diaz10022, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & nCount_RNA > 4000 & nCount_RNA < 25000 & percent.mt < 1.5)
diaz10022 <- NormalizeData(diaz10022)
diaz10022 <- FindVariableFeatures(diaz10022, selection.method = "vst", nfeatures = 2000)
diaz10022 <- ScaleData(diaz10022, features = rownames(diaz10022))
diaz10022 <- RunPCA(diaz10022, features = VariableFeatures(object = diaz10022))

DimPlot(diaz10022, reduction = "pca") 
ElbowPlot(diaz10022, ndims = 45)

d_10022 <- 15
diaz10022 <- FindNeighbors(diaz10022, dims = 1:d_10022)
diaz10022 <- FindClusters(diaz10022, resolution = 0.8)
diaz10022 <- RunUMAP(diaz10022, dims = 1:d_10022)
DimPlot(diaz10022, reduction = "umap")

markers_10022 <- FindAllMarkers(diaz10022, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_10022 %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)

# 12264 ----
diaz12264 <- Read10X(data.dir = "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119525_SF12264_GBM")
diaz12264 <- CreateSeuratObject(counts = diaz12264, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz12264[["percent.mt"]] <- PercentageFeatureSet(diaz12264, pattern = "^MT")
VlnPlot(diaz12264, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz12264 <- subset(x=diaz12264, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 20)
diaz12264 <- NormalizeData(diaz12264)
diaz12264 <- FindVariableFeatures(diaz12264, selection.method = "vst", nfeatures = 2000)
diaz12264 <- ScaleData(diaz12264, features = rownames(diaz12264))
diaz12264 <- RunPCA(diaz12264, features = VariableFeatures(object = diaz12264))

DimPlot(diaz12264, reduction = "pca") 
ElbowPlot(diaz12264, ndims = 45)

d_12264 <- 12
diaz12264 <- FindNeighbors(diaz12264, dims = 1:d_12264)
diaz12264 <- FindClusters(diaz12264, resolution = 0.8)
diaz12264 <- RunUMAP(diaz12264, dims = 1:d_12264)
DimPlot(diaz12264, reduction = "umap")

markers_12264 <- FindAllMarkers(diaz12264, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_12264 %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)

# 4297 ----
diaz4297 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119527_SF4297_GBM")
diaz4297 <- CreateSeuratObject(counts = diaz4297, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz4297[["percent.mt"]] <- PercentageFeatureSet(diaz4297, pattern = "^MT")
VlnPlot(diaz4297, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz4297 <- subset(x=diaz4297, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 1)
diaz4297 <- NormalizeData(diaz4297)
diaz4297 <- FindVariableFeatures(diaz4297, selection.method = "vst", nfeatures = 2000)
diaz4297 <- ScaleData(diaz4297, features = rownames(diaz4297))
diaz4297 <- RunPCA(diaz4297, features = VariableFeatures(object = diaz4297))

DimPlot(diaz4297, reduction = "pca") 
ElbowPlot(diaz4297, ndims = 45)

d_4297 <- 13
diaz4297 <- FindNeighbors(diaz4297, dims = 1:d_4297)
diaz4297 <- FindClusters(diaz4297, resolution = 0.8)
diaz4297 <- RunUMAP(diaz4297, dims = 1:d_4297)
DimPlot(diaz4297, reduction = "umap")

diaz4297 <- saveRDS(diaz4297, "output/diaz/diaz4297.rds")
d4297 <- readRDS("output/diaz/diaz4297.rds")

markers_4297 <- FindAllMarkers(diaz4297, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_4297 %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)

# 11644 ----
diaz11644 <- Read10X(data.dir = "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119534_SF11644_GBM", gene.column = 1)
diaz11644 <- CreateSeuratObject(counts = diaz11644, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz11644[["percent.mt"]] <- PercentageFeatureSet(diaz11644, pattern = "^MT")
VlnPlot(diaz11644, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz11644 <- subset(x=diaz11644, subset = nFeature_RNA > 1200 & nFeature_RNA < 8000 & nCount_RNA > 200 & nCount_RNA < 45000 & percent.mt < 5)
diaz11644 <- NormalizeData(diaz11644)
diaz11644 <- FindVariableFeatures(diaz11644, selection.method = "vst", nfeatures = 2000)
diaz11644 <- ScaleData(diaz11644, features = rownames(diaz11644))
diaz11644 <- RunPCA(diaz11644, features = VariableFeatures(object = diaz11644))

DimPlot(diaz11644, reduction = "pca") 
ElbowPlot(diaz11644, ndims = 45)

d_11644 <- 13
diaz11644 <- FindNeighbors(diaz11644, dims = 1:d_11644)
diaz11644 <- FindClusters(diaz11644, resolution = 0.8)
diaz11644 <- RunUMAP(diaz11644, dims = 1:d_11644)
DimPlot(diaz11644, reduction = "umap")

markers_11644 <- FindAllMarkers(diaz11644, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_11644 %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)

# 11681 ----
diaz11681 <- Read10X(data.dir = "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4658373_SF11681_GBM")
diaz11681 <- CreateSeuratObject(counts = diaz11681, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz11681[["percent.mt"]] <- PercentageFeatureSet(diaz11681, pattern = "^MT")
VlnPlot(diaz11681, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz11681 <- subset(x=diaz11681, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & nCount_RNA > 200 & nCount_RNA < 38000 & percent.mt < 10)
diaz11681 <- NormalizeData(diaz11681)
diaz11681 <- FindVariableFeatures(diaz11681, selection.method = "vst", nfeatures = 2000)
diaz11681 <- ScaleData(diaz11681, features = rownames(diaz11681))
diaz11681 <- RunPCA(diaz11681, features = VariableFeatures(object = diaz11681))

DimPlot(diaz11681, reduction = "pca") 
ElbowPlot(diaz11681, ndims = 45)

d_11681 <- 14
diaz11681 <- FindNeighbors(diaz11681, dims = 1:d_11681)
diaz11681 <- FindClusters(diaz11681, resolution = 0.8)
diaz11681 <- RunUMAP(diaz11681, dims = 1:d_11681)
DimPlot(diaz11681, reduction = "umap")

markers_11681 <- FindAllMarkers(diaz11681, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_11681 %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)