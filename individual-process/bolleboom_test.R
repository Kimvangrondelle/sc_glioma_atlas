# load data
bboom <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/Bolleboom-Gao/H243_GBM/H243_filtered_feature_bc_matrix")
# create seurat object and prefilter based on min.cells and min.features.
bboom <- CreateSeuratObject(counts = bboom, project = "sc_glioma_atlas", min.cells=3, min.features=200)
# calculate percent mitochondrial genes
bboom[["percent.mt"]] <- PercentageFeatureSet(bboom, pattern = "^MT")
# make violin plots to determine parameters for further analysis
VlnPlot(bboom, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# subset on found parameters
bboom <- subset(x=bboom, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 200 & nCount_RNA < 18000 & percent.mt < 1)
# process data: Normalizing, Finding variable features, Scaling and PCA
bboom <- NormalizeData(bboom)
bboom <- FindVariableFeatures(bboom, selection.method = "vst", nfeatures = 2000)
bboom <- ScaleData(bboom, features = rownames(bboom))
bboom <- RunPCA(bboom, features = VariableFeatures(object = bboom))

# Elbow plot to determine optimal number of principal components
DimPlot(bboom, reduction = "pca") 
ElbowPlot(bboom, ndims = 45)


d <- 16 # -- optimal principal components
# clustering
bboom <- FindNeighbors(bboom, dims = 1:d)
bboom <- FindClusters(bboom, resolution = 0.8)
bboom <- RunUMAP(bboom, dims = 1:d)
# umap to see how cells cluster together
DimPlot(bboom, reduction = "umap")
# save object
# saveRDS(bboom, file = "output/bolleboom/H243.rds")
# find DE genes per cluster using find all markers
markers <- FindAllMarkers(bboom, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top10)


