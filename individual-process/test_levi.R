
# load data
sid <- 'van_Hijfte_Sample_Y'
object_1 <- Read10X(data.dir = "../../../mnt/neuro-genomic-1-ro/Glimmunology/LGG_project/10x_singlenuc_RNAseq/Optimisation_10x_snRNAseq/Glioma_Y_and_O/Levi2_Glioma_Y/outs/raw_feature_bc_matrix")
# create seurat object and prefilter based on min.cells and min.features.
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "glioma_glim")
# calculate percent mitochondrial genes
mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
# make violin plots to determine parameters for further analysis
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 
# subset on found parameters
object_1 <- subset(x= object_1, subset = nFeature_RNA > 1400 & nFeature_RNA < 4500 & nCount_RNA > 2200 & nCount_RNA < 14000 & percent.mito < 0.025)
# process data: Normalizing, Finding variable features, Scaling and PCA
object_1 <- NormalizeData(object_1)
object_1 <- FindVariableFeatures(object_1, selection.method = "vst", nfeatures = 2000)
object_1 <- ScaleData(object_1, features = all.genes_levi)

object_1 <- RunPCA(object_1, features = VariableFeatures(object=object_1))

# Elbow plot to determine optimal number of principal components
DimPlot(object_1, reduction = "pca") + NoLegend()
ElbowPlot(object_1, ndims=45)

d_levi <- 16 # -- optimal principal components
# clustering 
object_1 <- FindNeighbors(object_1, dims = 1:d_levi)
object_1 <- FindClusters(object_1, resolution = 0.8)
# umap to see how cells cluster together
object_1 <- RunUMAP(object_1, dims=1:d_levi)
DimPlot(object_1, reduction = "umap")
# find DE genes per cluster using find all markers
all_markers_levi <- FindAllMarkers(object_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_levi <- all_markers_levi %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_levi)


