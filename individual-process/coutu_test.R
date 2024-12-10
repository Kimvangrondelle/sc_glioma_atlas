# BT338 combined ----
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT338_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou338_1", min.cells=3, min.features=200)

object_2 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT338_2of2.filtered_gene_matrices")
object_2 <- CreateSeuratObject(counts = object_2, project = "cou338_2", min.cells=3, min.features=200)

objects.combined <- merge(object_1, y = object_2, add.cell.ids = c("1/2", "2/2"), project = "coutu338")
objects.combined

cou338 <- objects.combined
cou338 <- JoinLayers(cou338)
cou338[["percent.mt"]] <- PercentageFeatureSet(cou338, pattern = "^MT")

VlnPlot(cou338, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou338 <- subset(x= cou338, subset = nFeature_RNA > 1000 & nFeature_RNA < 7800 & nCount_RNA > 200 & nCount_RNA < 50000 & percent.mt < 2)

cou338 <- NormalizeData(cou338)
cou338 <- FindVariableFeatures(cou338, selection.method = "vst", nfeatures = 2000)
top20_338 <- head(VariableFeatures(cou338), 20)

cou338 <- ScaleData(cou338, features = rownames(cou338))
cou338 <- RunPCA(cou338, features = VariableFeatures(object=cou338))
DimPlot(cou338, reduction = "pca") + NoLegend()
ElbowPlot(cou338, ndims=45)

d_338 <- 12 
cou338 <- FindNeighbors(cou338, dims = 1:d_338)
cou338 <- FindClusters(cou338, resolution = 0.8)
cou338 <- RunUMAP(cou338, dims=1:d_338)
DimPlot(cou338, reduction = "umap")

all_markers_338 <- FindAllMarkers(cou338, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_338 <- all_markers_338 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top5_markers_338)


# BT338 1/2 ----
sid_338 <- 'coutu_sample_338'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT338_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou338", min.cells=3, min.features=200)

cou338 <- object_1
cou338[["percent.mt"]] <- PercentageFeatureSet(cou338, pattern = "^MT")

VlnPlot(cou338, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou338 <- subset(x= cou338, subset = nFeature_RNA > 1000 & nFeature_RNA < 7800 & nCount_RNA > 200 & nCount_RNA < 50000 & percent.mt < 2)

cou338 <- NormalizeData(cou338)
cou338 <- FindVariableFeatures(cou338, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(cou338), 20)

cou338 <- ScaleData(cou338, features = rownames(cou338))
cou338 <- RunPCA(cou338, features = VariableFeatures(object=cou338))
DimPlot(cou338, reduction = "pca") + NoLegend()
ElbowPlot(cou338, ndims=45)

d_338 <- 12 
cou338 <- FindNeighbors(cou338, dims = 1:d_338)
cou338 <- FindClusters(cou338, resolution = 0.8)

cou338 <- RunUMAP(cou338, dims=1:d_338)
DimPlot(cou338, reduction = "umap")



# BT338 2/2 ----
sid_3382 <- 'coutu_sample_3382'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT338_2of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou3382", min.cells=3, min.features=200)

cou3382 <- object_1
cou3382[["percent.mt"]] <- PercentageFeatureSet(cou3382, pattern = "^MT")

VlnPlot(cou3382, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou3382 <- subset(x= cou3382, subset = nFeature_RNA > 1000 & nFeature_RNA < 7800 & nCount_RNA > 200 & nCount_RNA < 50000 & percent.mt < 2)

cou3382 <- NormalizeData(cou3382)
cou3382 <- FindVariableFeatures(cou3382, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(cou3382), 20)

cou3382 <- ScaleData(cou3382, features = rownames(cou3382))
cou3382 <- RunPCA(cou3382, features = VariableFeatures(object=cou3382))
DimPlot(cou3382, reduction = "pca") + NoLegend()
ElbowPlot(cou3382, ndims=45)

d_3382 <- 12 
cou3382 <- FindNeighbors(cou3382, dims = 1:d_3382)
cou3382 <- FindClusters(cou3382, resolution = 0.8)

cou3382 <- RunUMAP(cou3382, dims=1:d_3382)
DimPlot(cou3382, reduction = "umap")



# BT363 1/2 ----
sid_363 <- 'coutu_sample_363'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT363_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou363", min.cells=3, min.features=200)

cou363 <- object_1
cou363[["percent.mt"]] <- PercentageFeatureSet(cou363, pattern = "^MT")

VlnPlot(cou363, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou363 <- subset(x= cou363, subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 2)

cou363 <- NormalizeData(cou363)
cou363 <- FindVariableFeatures(cou363, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(cou363), 20)

cou363 <- ScaleData(cou363, features = rownames(cou363))
cou363 <- RunPCA(cou363, features = VariableFeatures(object=cou363))
DimPlot(cou363, reduction = "pca") + NoLegend()
ElbowPlot(cou363, ndims=45)

d_363 <- 12 
cou363 <- FindNeighbors(cou363, dims = 1:d_363)
cou363 <- FindClusters(cou363, resolution = 0.8)

cou363 <- RunUMAP(cou363, dims=1:d_363)
DimPlot(cou363, reduction = "umap")



# BT363 2/2 ----
sid_3632 <- 'coutu_sample_3632'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT363_2of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou3632", min.cells=3, min.features=200)

cou3632 <- object_1
cou3632[["percent.mt"]] <- PercentageFeatureSet(cou3632, pattern = "^MT")

VlnPlot(cou3632, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou3632 <- subset(x= cou3632, subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 2)

cou3632 <- NormalizeData(cou3632)
cou3632 <- FindVariableFeatures(cou3632, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(cou3632), 20)

cou3632 <- ScaleData(cou3632, features = rownames(cou3632))
cou3632 <- RunPCA(cou3632, features = VariableFeatures(object=cou3632))
DimPlot(cou3632, reduction = "pca") + NoLegend()
ElbowPlot(cou3632, ndims=45)

d_3632 <- 12 
cou3632 <- FindNeighbors(cou3632, dims = 1:d_3632)
cou3632 <- FindClusters(cou3632, resolution = 0.8)

cou3632 <- RunUMAP(cou3632, dims=1:d_3632)
DimPlot(cou3632, reduction = "umap")


# BT363 combined ----
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT363_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou363_1", min.cells=3, min.features=200)

object_2 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT363_2of2.filtered_gene_matrices")
object_2 <- CreateSeuratObject(counts = object_2, project = "cou363_2", min.cells=3, min.features=200)

objects.combined <- merge(object_1, y = object_2, add.cell.ids = c("1/2", "2/2"), project = "coutu363")

cou363 <- objects.combined
cou363 <- JoinLayers(cou363)
cou363[["percent.mt"]] <- PercentageFeatureSet(cou363, pattern = "^MT")

VlnPlot(cou363, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou363 <- subset(x= cou363, subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 2)

cou363 <- NormalizeData(cou363)
cou363 <- FindVariableFeatures(cou363, selection.method = "vst", nfeatures = 2000)
top20_363 <- head(VariableFeatures(cou363), 20)

cou363 <- ScaleData(cou363, features = rownames(cou363))
cou363 <- RunPCA(cou363, features = VariableFeatures(object=cou363))
DimPlot(cou363, reduction = "pca") + NoLegend()
ElbowPlot(cou363, ndims=45)

d_363 <- 12 
cou363 <- FindNeighbors(cou363, dims = 1:d_363)
cou363 <- FindClusters(cou363, resolution = 0.8)
cou363 <- RunUMAP(cou363, dims=1:d_363)
DimPlot(cou363, reduction = "umap")

all_markers_363 <- FindAllMarkers(cou363, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_363 <- all_markers_363 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top5_markers_363)



# BT364 1/2 ----
sid_364 <- 'coutu_sample_364'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou364", min.cells=3, min.features=200)

cou364 <- object_1
cou364[["percent.mt"]] <- PercentageFeatureSet(cou364, pattern = "^MT")

VlnPlot(cou364, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou364 <- subset(x= cou364, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 200 & nCount_RNA < 35000 & percent.mt < 2)

cou364 <- NormalizeData(cou364)
cou364 <- FindVariableFeatures(cou364, selection.method = "vst", nfeatures = 2000)

cou364 <- ScaleData(cou364, features = rownames(cou364))
cou364 <- RunPCA(cou364, features = VariableFeatures(object=cou364))
DimPlot(cou364, reduction = "pca") + NoLegend()
ElbowPlot(cou364, ndims=45)

d_364 <- 14 
cou364 <- FindNeighbors(cou364, dims = 1:d_364)
cou364 <- FindClusters(cou364, resolution = 0.8)

cou364 <- RunUMAP(cou364, dims=1:d_364)
DimPlot(cou364, reduction = "umap")



# BT364 2/2 ----
sid_3642 <- 'coutu_sample_3642'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_2of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou3642", min.cells=3, min.features=200)

cou3642 <- object_1
cou3642[["percent.mt"]] <- PercentageFeatureSet(cou3642, pattern = "^MT")

VlnPlot(cou3642, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou3642 <- subset(x= cou3642, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 200 & nCount_RNA < 35000 & percent.mt < 2)

cou3642 <- NormalizeData(cou3642)
cou3642 <- FindVariableFeatures(cou3642, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(cou3642), 20)

cou3642 <- ScaleData(cou3642, features = rownames(cou3642))
cou3642 <- RunPCA(cou3642, features = VariableFeatures(object=cou3642))
DimPlot(cou3642, reduction = "pca") + NoLegend()
ElbowPlot(cou3642, ndims=45)

d_3642 <- 14 
cou3642 <- FindNeighbors(cou3642, dims = 1:d_3642)
cou3642 <- FindClusters(cou3642, resolution = 0.8)

cou3642 <- RunUMAP(cou3642, dims=1:d_3642)
DimPlot(cou3642, reduction = "umap")



# BT364 combined ----
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou364_1", min.cells=3, min.features=200)

object_2 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_2of2.filtered_gene_matrices")
object_2 <- CreateSeuratObject(counts = object_2, project = "cou364_2", min.cells=3, min.features=200)

objects.combined <- merge(object_1, y = object_2, add.cell.ids = c("1/2", "2/2"), project = "coutu364")

cou364 <- objects.combined
cou364 <- JoinLayers(cou364)
cou364[["percent.mt"]] <- PercentageFeatureSet(cou364, pattern = "^MT")

VlnPlot(cou364, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou364 <- subset(x= cou364, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 200 & nCount_RNA < 35000 & percent.mt < 2)

cou364 <- NormalizeData(cou364)
cou364 <- FindVariableFeatures(cou364, selection.method = "vst", nfeatures = 2000)
top20_364 <- head(VariableFeatures(cou364), 20)

cou364 <- ScaleData(cou364, features = rownames(cou364))
cou364 <- RunPCA(cou364, features = VariableFeatures(object=cou364))
DimPlot(cou364, reduction = "pca") + NoLegend()
ElbowPlot(cou364, ndims=45)

d_364 <- 14 
cou364 <- FindNeighbors(cou364, dims = 1:d_364)
cou364 <- FindClusters(cou364, resolution = 0.8)
cou364 <- RunUMAP(cou364, dims=1:d_364)
DimPlot(cou364, reduction = "umap")

all_markers_364 <- FindAllMarkers(cou364, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_364 <- all_markers_364 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top5_markers_364)



# BT397 1/2 ----
sid_397 <- 'coutu_sample_397'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT397_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou397", min.cells=3, min.features=200)

cou397 <- object_1
cou397[["percent.mt"]] <- PercentageFeatureSet(cou397, pattern = "^MT")

VlnPlot(cou397, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou397 <- subset(x= cou397, subset = nFeature_RNA > 700 & nFeature_RNA < 4300 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 2)

cou397 <- NormalizeData(cou397)
cou397 <- FindVariableFeatures(cou397, selection.method = "vst", nfeatures = 2000)

cou397 <- ScaleData(cou397, features = rownames(cou397))
cou397 <- RunPCA(cou397, features = VariableFeatures(object=cou397))
DimPlot(cou397, reduction = "pca") + NoLegend()
ElbowPlot(cou397, ndims=45)

d_397 <- 12 
cou397 <- FindNeighbors(cou397, dims = 1:d_397)
cou397 <- FindClusters(cou397, resolution = 0.8)

cou397 <- RunUMAP(cou397, dims=1:d_397)
DimPlot(cou397, reduction = "umap")


# BT397 2/2 ----
sid_3972 <- 'coutu_sample_3972'
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT397_2of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou3972", min.cells=3, min.features=200)

cou3972 <- object_1
cou3972[["percent.mt"]] <- PercentageFeatureSet(cou3972, pattern = "^MT")

VlnPlot(cou3972, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou3972 <- subset(x= cou3972, subset = nFeature_RNA > 700 & nFeature_RNA < 4500 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 2)

cou3972 <- NormalizeData(cou3972)
cou3972 <- FindVariableFeatures(cou3972, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(cou3972), 20)

cou3972 <- ScaleData(cou3972, features = rownames(cou3972))
cou3972 <- RunPCA(cou3972, features = VariableFeatures(object=cou3972))
DimPlot(cou3972, reduction = "pca") + NoLegend()
ElbowPlot(cou3972, ndims=45)

d_3972 <- 13 
cou3972 <- FindNeighbors(cou3972, dims = 1:d_3972)
cou3972 <- FindClusters(cou3972, resolution = 0.8)

cou3972 <- RunUMAP(cou3972, dims=1:d_3972)
DimPlot(cou3972, reduction = "umap")


# BT397 combined ----
object_1 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT397_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = "cou397_1", min.cells=3, min.features=200)

object_2 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT397_2of2.filtered_gene_matrices")
object_2 <- CreateSeuratObject(counts = object_2, project = "cou397_2", min.cells=3, min.features=200)

objects.combined <- merge(object_1, y = object_2, add.cell.ids = c("1/2", "2/2"), project = "coutu397")

cou397 <- objects.combined
cou397 <- JoinLayers(cou397)
cou397[["percent.mt"]] <- PercentageFeatureSet(cou397, pattern = "^MT")

VlnPlot(cou397, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

cou397 <- subset(x= cou397, subset = nFeature_RNA > 700 & nFeature_RNA < 4300 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 2)

cou397 <- NormalizeData(cou397)
cou397 <- FindVariableFeatures(cou397, selection.method = "vst", nfeatures = 2000)
top20_397 <- head(VariableFeatures(cou397), 20)

cou397 <- ScaleData(cou397, features = rownames(cou397))
cou397 <- RunPCA(cou397, features = VariableFeatures(object=cou397))
DimPlot(cou397, reduction = "pca") + NoLegend()
ElbowPlot(cou397, ndims=45)

d_397 <- 13 
cou397 <- FindNeighbors(cou397, dims = 1:d_397)
cou397 <- FindClusters(cou397, resolution = 0.8)
cou397 <- RunUMAP(cou397, dims=1:d_397)
DimPlot(cou397, reduction = "umap")

all_markers_397 <- FindAllMarkers(cou397, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_397 <- all_markers_397 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top5_markers_397)