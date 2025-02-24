# it was tried to integrate the samples per dataset. However they were not used in the research. 
couturier <- list.files(path= "output/couturier", pattern = "*.rds", full.names = TRUE)

objs <- list()
#loop through samples
for (file in couturier) {
  obj <- readRDS(file = file)
  obj$sample_id <- file
  objs[[length(objs) + 1]] <- obj
}
#merging of the samples
if (length(objs) > 1) {
  combined <- merge(objs[[1]], y= objs[2:length(objs)], add.cell.ids=names(objs))
}  
 
#preprocess
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

combined <- RunPCA(combined, npcs=30)
#integration of the samples on sample_id
combined <- RunHarmony(combined, group.by="sample_id", assay.use="RNA")
coutu <- combined
#same process as individual samples
coutu[["percent.mt"]] <- PercentageFeatureSet(coutu, pattern = "^MT")
VlnPlot(coutu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#subsetting, normalizing, scaling and pca
coutu <- subset(x=coutu, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA > 200 & nCount_RNA < 45000 & percent.mt < 2)
coutu <- NormalizeData(coutu)
coutu <- FindVariableFeatures(coutu, selection.method = "vst", nfeatures = 2000)
coutu <- ScaleData(coutu, features = rownames(coutu))
coutu <- RunPCA(coutu, features = VariableFeatures(object = coutu))

DimPlot(coutu, reduction = "pca") 
ElbowPlot(coutu, ndims = 45)

#12
d <- 15 #optimal number of pcs
coutu <- FindNeighbors(coutu, dims = 1:d)
coutu <- FindClusters(coutu, resolution = 0.8)
coutu <- RunUMAP(coutu, dims = 1:d)
DimPlot(coutu, reduction = "umap")
saveRDS(coutu, file = "output/integrated/integrated_couturier.rds")

integrated_coutu <- readRDS(file = "output/integrated/integrated_couturier.rds")
DimPlot(integrated_coutu, reduction = "umap", group.by = "sample_id")

  
diaz <- list.files(path= "output/diaz", pattern = "*.rds", full.names = TRUE)


yuanie <- list.files(path= "output/yuan", pattern = "*.rds", full.names = TRUE)

objs <- list()
for (file in yuanie) {
  obj <- readRDS(file = file)
  obj$sample_id <- file
  objs[[length(objs) + 1]] <- obj
}
if (length(objs) > 1) {
  combined <- merge(objs[[1]], y= objs[2:length(objs)], add.cell.ids=names(objs))
}  
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

combined <- RunPCA(combined, npcs=30)

combined <- RunHarmony(combined, group.by="sample_id", assay.use="RNA")
yuan <- combined

yuan[["percent.mt"]] <- PercentageFeatureSet(yuan, pattern = "^MT")
VlnPlot(yuan, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

yuan <- subset(x=yuan, subset = nFeature_RNA > 200 & nFeature_RNA < 4100 & nCount_RNA > 200 & nCount_RNA < 10000 & percent.mt < 2)
yuan <- NormalizeData(yuan)
yuan <- FindVariableFeatures(yuan, selection.method = "vst", nfeatures = 2000)
yuan <- ScaleData(yuan, features = rownames(yuan))
yuan <- RunPCA(yuan, features = VariableFeatures(object = yuan))

DimPlot(yuan, reduction = "pca") 
ElbowPlot(yuan, ndims = 45)

#12
d <- 15
yuan <- FindNeighbors(yuan, dims = 1:d)
yuan <- FindClusters(yuan, resolution = 0.8)
yuan <- RunUMAP(yuan, dims = 1:d)
DimPlot(yuan, reduction = "umap")
saveRDS(yuan, file = "output/integrated/integrated_yuan.rds")


diazie <- list.files(path= "output/diaz", pattern = "*.rds", full.names = TRUE)

objs <- list()
for (file in diazie) {
  obj <- readRDS(file = file)
  obj$sample_id <- file
  objs[[length(objs) + 1]] <- obj
}
if (length(objs) > 1) {
  combined <- merge(objs[[1]], y= objs[2:length(objs)], add.cell.ids=names(objs))
}  
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

combined <- RunPCA(combined, npcs=30)

combined <- RunHarmony(combined, group.by="sample_id", assay.use="RNA")
diaz <- combined

diaz[["percent.mt"]] <- PercentageFeatureSet(diaz, pattern = "^MT")
VlnPlot(diaz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz <- subset(x=diaz, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA > 200 & nCount_RNA < 45000 & percent.mt < 20)
diaz <- NormalizeData(diaz)
diaz <- FindVariableFeatures(diaz, selection.method = "vst", nfeatures = 2000)
diaz <- ScaleData(diaz, features = rownames(diaz))
diaz <- RunPCA(diaz, features = VariableFeatures(object = diaz))

DimPlot(diaz, reduction = "pca") 
ElbowPlot(diaz, ndims = 45)

#12
d <- 14
diaz <- FindNeighbors(diaz, dims = 1:d)
diaz <- FindClusters(diaz, resolution = 0.8)
diaz <- RunUMAP(diaz, dims = 1:d)
DimPlot(diaz, reduction = "umap")
saveRDS(diaz, file = "output/integrated/integrated_diaz.rds")
