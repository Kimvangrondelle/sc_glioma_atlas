library(Matrix)
library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)

# for dot plot met goede genes -------
DotPlot(
  yuan32,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))



# code prepare 25 ----
sid_25 <- 'yuan_sample_PJ025'
yuan25 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758474_PJ025")

yuan25 <- CreateSeuratObject(counts=yuan25, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan25[["percent.mt"]] <- PercentageFeatureSet(yuan25, pattern = "^MT")

VlnPlot(yuan25, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan25 <- subset(x= yuan25, subset = nFeature_RNA > 700 & nFeature_RNA < 3500 & nCount_RNA > 200 & nCount_RNA < 9000 & percent.mt < 2)

yuan25 <- NormalizeData(yuan25)

yuan25 <- FindVariableFeatures(yuan25, selection.method = "vst", nfeatures = 2000)
top20_25 <- head(VariableFeatures(yuan25), 20)
# top20_25

plot1_25 <- VariableFeaturePlot(yuan25)
plot2_25 <- LabelPoints(plot = plot1_25, points = top20_25)
# plot1_25 + plot2_25

# all.genes_25 <- rownames(yuan25)
yuan25 <- ScaleData(yuan25, features = rownames(yuan25))

yuan25 <- RunPCA(yuan25, features = VariableFeatures(object=yuan25))

# print(yuan25[["pca"]], dims=1:5, nfeatures=5)
DimPlot(yuan25, reduction = "pca") + NoLegend()
ElbowPlot(yuan25, ndims=45)

d_25 <- 9 
yuan25 <- FindNeighbors(yuan25, dims = 1:d_25)
yuan25 <- FindClusters(yuan25, resolution = 0.8)
#have to find optimal number of resolution
# yuan_25_num_clus <- yuan25

yuan25 <- RunUMAP(yuan25, dims=1:d_25)
DimPlot(yuan25, reduction = "umap")

all_markers_25 <- FindAllMarkers(yuan25, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_25 <- all_markers_25 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(top5_markers_25)



# assign celltype to cluster nested ----

new.cluster.id.25 <- c("tumor", "dividing tumor", "macrophage", "astrocyte", "oligodendrocyte", "pericyte", "oligodendrocyte", "endothelial", "macrophage")
names(new.cluster.id.25) <- levels(yuan25)
yuan25 <- RenameIdents(yuan25, new.cluster.id.25)
DimPlot(yuan25, reduction = "umap", label = TRUE)

# code prepare 30 ----
library(Matrix)
library(Seurat)
library(patchwork)
sid_30 <- 'yuan_sample_PJ030'
yuan30 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758475_PJ030")

yuan30 <- CreateSeuratObject(counts=yuan30, project = "sc_glioma_atlas", min.cells=3, min.features=200)
yuan30[["percent.mt"]] <- PercentageFeatureSet(yuan30, pattern = "^MT")

VlnPlot(yuan30, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan30 <- subset(x= yuan30, subset = nFeature_RNA > 700 & nFeature_RNA < 3600 & nCount_RNA > 200 & nCount_RNA < 5000 & percent.mt < 2)
yuan30 <- NormalizeData(yuan30)

yuan30 <- FindVariableFeatures(yuan30, selection.method = "vst", nfeatures = 2000)
top20_30 <- head(VariableFeatures(yuan30), 20)

# plot1_30 <- VariableFeaturePlot(yuan30)
# plot2_30 <- LabelPoints(plot = plot1_30, points = top20_30)
# plot1_30 + plot2_30

# all.genes_30 <- rownames(yuan30)
yuan30 <- ScaleData(yuan30, features = rownames(yuan30))
yuan30 <- RunPCA(yuan30, features = VariableFeatures(object=yuan30))

DimPlot(yuan30, reduction = "pca") + NoLegend()
ElbowPlot(yuan30, ndims=45)

d_30 <- 11 
yuan30 <- FindNeighbors(yuan30, dims = 1:d_30)
yuan30 <- FindClusters(yuan30, resolution = 0.8)
#have to find optimal number of resolution

yuan30 <- RunUMAP(yuan30, dims=1:d_30)
DimPlot(yuan30, reduction = "umap")

all_markers_30 <- FindAllMarkers(yuan30, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_30 <- all_markers_30 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_30)



new.cluster.id.30 <- c("tumor", "tumor", "dividing tumor", "macrophage", "oligodendrocyte", "tumor", "endothelial", "pericyte", "neuron", "oligodendrocyte")
names(new.cluster.id.30) <- levels(yuan30)
yuan30 <- RenameIdents(yuan30, new.cluster.id.30)
DimPlot(yuan30, reduction = "umap", label = TRUE)




# code prepare 32 ----
sid_32 <- 'yuan_sample_PJ032'
yuan32 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758476_PJ032")

yuan32 <- CreateSeuratObject(counts=yuan32, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan32[["percent.mt"]] <- PercentageFeatureSet(yuan32, pattern = "^MT")

VlnPlot(yuan32, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan32 <- subset(x= yuan32, subset = nFeature_RNA > 700 & nFeature_RNA < 4000 & nCount_RNA > 100 & nCount_RNA < 9000 & percent.mt < 2)

yuan32 <- NormalizeData(yuan32)

yuan32 <- FindVariableFeatures(yuan32, selection.method = "vst", nfeatures = 2000)
top20_32 <- head(VariableFeatures(yuan32), 20)
# top20_32

plot1_32 <- VariableFeaturePlot(yuan32)
plot2_32 <- LabelPoints(plot = plot1_32, points = top20_32)
plot1_32 + plot2_32

# all.genes_32 <- rownames(yuan32)
yuan32 <- ScaleData(yuan32, features = rownames(yuan32))

yuan32 <- RunPCA(yuan32, features = VariableFeatures(object=yuan32))

DimPlot(yuan32, reduction = "pca") + NoLegend()
ElbowPlot(yuan32, ndims=45)

d_32 <- 7 
yuan32 <- FindNeighbors(yuan32, dims = 1:d_32)
yuan32 <- FindClusters(yuan32, resolution = 0.8)
#have to find optimal number of resolution

yuan32 <- RunUMAP(yuan32, dims=1:d_32)
DimPlot(yuan32, reduction = "umap")

all_markers_32 <- FindAllMarkers(yuan32, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_32 <- all_markers_32 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_32)


DotPlot(
  yuan32,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))

new.cluster.id.32 <- c("macrophage", "tumor", "tumor", "dividing tumor")
names(new.cluster.id.32) <- levels(yuan32)
yuan32 <- RenameIdents(yuan32, new.cluster.id.32)
DimPlot(yuan32, reduction = "umap", label = TRUE)




# code prepare 35 ----
sid_35 <- 'yuan_sample_PJ035'
yuan35 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758477_PJ035")

yuan35 <- CreateSeuratObject(counts=yuan35, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan35[["percent.mt"]] <- PercentageFeatureSet(yuan35, pattern = "^MT")

VlnPlot(yuan35, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan35 <- subset(x= yuan35, subset = nFeature_RNA > 700 & nFeature_RNA < 2200 & nCount_RNA > 200 & nCount_RNA < 3900 & percent.mt < 2)

yuan35 <- NormalizeData(yuan35)

yuan35 <- FindVariableFeatures(yuan35, selection.method = "vst", nfeatures = 2000)
top20_35 <- head(VariableFeatures(yuan35), 20)

plot1_35 <- VariableFeaturePlot(yuan35)
plot2_35 <- LabelPoints(plot = plot1_35, points = top20_35)
# plot1_35 + plot2_35

# all.genes_35 <- rownames(yuan35)
yuan35 <- ScaleData(yuan35, features = rownames(yuan35))

yuan35 <- RunPCA(yuan35, features = VariableFeatures(object=yuan35))

DimPlot(yuan35, reduction = "pca") + NoLegend()
ElbowPlot(yuan35, ndims=45)

d_35 <- 9 
yuan35 <- FindNeighbors(yuan35, dims = 1:d_35)
yuan35 <- FindClusters(yuan35, resolution = 0.8)
#have to find optimal number of resolution

yuan35 <- RunUMAP(yuan35, dims=1:d_35)
DimPlot(yuan35, reduction = "umap")

all_markers_35 <- FindAllMarkers(yuan35, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_35 <- all_markers_35 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_35)



DotPlot(
  yuan32,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))

new.cluster.id.35 <- c("tumor", "tumor", "endothelial", "tumor", "macrophage", "neuron", "astrocyte", "endothelial", "pericyte")
names(new.cluster.id.35) <- levels(yuan35)
yuan35 <- RenameIdents(yuan35, new.cluster.id.35)
DimPlot(yuan35, reduction = "umap", label = TRUE)


# code prepare 48 ----
sid_48 <- 'yuan_sample_PJ048'
yuan48 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2940098_PJ048")

yuan48 <- CreateSeuratObject(counts=yuan48, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan48[["percent.mt"]] <- PercentageFeatureSet(yuan48, pattern = "^MT")

VlnPlot(yuan48, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan48 <- subset(x= yuan48, subset = nFeature_RNA > 700 & nFeature_RNA < 3400 & nCount_RNA > 200 & nCount_RNA < 8000 & percent.mt < 2)

yuan48 <- NormalizeData(yuan48)

yuan48 <- FindVariableFeatures(yuan48, selection.method = "vst", nfeatures = 2000)
top20_48 <- head(VariableFeatures(yuan48), 20)

plot1_48 <- VariableFeaturePlot(yuan48)
plot2_48 <- LabelPoints(plot = plot1_48, points = top20_48)
# plot1_48 + plot2_48

# all.genes_48 <- rownames(yuan48)
yuan48 <- ScaleData(yuan48, features = rownames(yuan48))

yuan48 <- RunPCA(yuan48, features = VariableFeatures(object=yuan48))

DimPlot(yuan48, reduction = "pca") + NoLegend()
ElbowPlot(yuan48, ndims=45)

d_48 <- 15
yuan48 <- FindNeighbors(yuan48, dims = 1:d_48)
yuan48 <- FindClusters(yuan48, resolution = 0.8)
#have to find optimal number of resolution

yuan48 <- RunUMAP(yuan48, dims=1:d_48)
DimPlot(yuan48, reduction = "umap")

all_markers_48 <- FindAllMarkers(yuan48, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_48 <- all_markers_48 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_48)


DotPlot(
  yuan48,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))

new.cluster.id.48 <- c("dividing tumor", "neuron", "oligodendrocyte", "neuron", "dividing tumor", "oligodendrocyte", "astrocyte", "astrocyte", "pericyte", "pericyte")
names(new.cluster.id.48) <- levels(yuan48)
yuan48 <- RenameIdents(yuan48, new.cluster.id.48)
DimPlot(yuan48, reduction = "umap", label = TRUE)


# code prepare 16 ----
sid_16 <- 'yuan_sample_PJ016'
yuan16 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758471_PJ016")

yuan16 <- CreateSeuratObject(counts=yuan16, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan16[["percent.mt"]] <- PercentageFeatureSet(yuan16, pattern = "^MT")

VlnPlot(yuan16, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan16 <- subset(x= yuan16, subset = nFeature_RNA > 1000 & nFeature_RNA < 6200 & nCount_RNA > 200 & nCount_RNA < 22000 & percent.mt < 2)

yuan16 <- NormalizeData(yuan16)

yuan16 <- FindVariableFeatures(yuan16, selection.method = "vst", nfeatures = 2000)
top20_16 <- head(VariableFeatures(yuan16), 20)

plot1_16 <- VariableFeaturePlot(yuan16)
plot2_16 <- LabelPoints(plot = plot1_16, points = top20_16)
plot1_16 + plot2_16

# all.genes_16 <- rownames(yuan16)
yuan16 <- ScaleData(yuan16, features = rownames(yuan16))

yuan16 <- RunPCA(yuan16, features = VariableFeatures(object=yuan16))

DimPlot(yuan16, reduction = "pca") + NoLegend()
ElbowPlot(yuan16, ndims=45)

d_16 <- 11
yuan16 <- FindNeighbors(yuan16, dims = 1:d_16)
yuan16 <- FindClusters(yuan16, resolution = 0.8)
#have to find optimal number of resolution

yuan16 <- RunUMAP(yuan16, dims=1:d_16)
DimPlot(yuan16, reduction = "umap")

all_markers_16 <- FindAllMarkers(yuan16, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_16 <- all_markers_16 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_16)



DotPlot(
  yuan16,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))

new.cluster.id.16 <- c("macrophage", "dividing tumor", "mesenchymal", "neuron", "oligodendrocyte", "neuron", "astrocyte", "dividing tumor")
names(new.cluster.id.16) <- levels(yuan16)
yuan16 <- RenameIdents(yuan16, new.cluster.id.16)
DimPlot(yuan16, reduction = "umap", label = TRUE)


# code prepare 17 ----
sid_17 <- 'yuan_sample_PJ017'
yuan17 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758472_PJ017")

yuan17 <- CreateSeuratObject(counts=yuan17, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan17[["percent.mt"]] <- PercentageFeatureSet(yuan17, pattern = "^MT")

VlnPlot(yuan17, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan17 <- subset(x= yuan17, subset = nFeature_RNA > 700 & nFeature_RNA < 4200 & nCount_RNA > 200 & nCount_RNA < 12000 & percent.mt < 2)

yuan17 <- NormalizeData(yuan17)

yuan17 <- FindVariableFeatures(yuan17, selection.method = "vst", nfeatures = 2000)
top20_17 <- head(VariableFeatures(yuan17), 20)

plot1_17 <- VariableFeaturePlot(yuan17)
plot2_17 <- LabelPoints(plot = plot1_17, points = top20_17)
# plot1_17 + plot2_17

# all.genes_17 <- rownames(yuan17)
yuan17 <- ScaleData(yuan17, features = rownames(yuan17))

yuan17 <- RunPCA(yuan17, features = VariableFeatures(object=yuan17))

DimPlot(yuan17, reduction = "pca") + NoLegend()
ElbowPlot(yuan17, ndims=45)

d_17 <- 15
yuan17 <- FindNeighbors(yuan17, dims = 1:d_17)
yuan17 <- FindClusters(yuan17, resolution = 0.8)
#have to find optimal number of resolution

yuan17 <- RunUMAP(yuan17, dims=1:d_17)
DimPlot(yuan17, reduction = "umap")

all_markers_17 <- FindAllMarkers(yuan17, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_17 <- all_markers_17 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_17)



DotPlot(
  yuan17,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))

new.cluster.id.17 <- c("macrophage", "astrocyte", "macrophage", "macrophage", "dividing tumor", "astrocyte", "Tcell")
names(new.cluster.id.17) <- levels(yuan17)
yuan17 <- RenameIdents(yuan17, new.cluster.id.17)
DimPlot(yuan17, reduction = "umap", label = TRUE)


# code prepare 18 ----
sid_18 <- 'yuan_sample_PJ018'
yuan18 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758473_PJ018")

yuan18 <- CreateSeuratObject(counts=yuan18, project = "sc_glioma_atlas", min.cells=3, min.features=200)

yuan18[["percent.mt"]] <- PercentageFeatureSet(yuan18, pattern = "^MT")

VlnPlot(yuan18, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

yuan18 <- subset(x= yuan18, subset = nFeature_RNA > 1000 & nFeature_RNA < 4800 & nCount_RNA > 200 & nCount_RNA < 12000 & percent.mt < 2)

yuan18 <- NormalizeData(yuan18)

yuan18 <- FindVariableFeatures(yuan18, selection.method = "vst", nfeatures = 2000)
top20_18 <- head(VariableFeatures(yuan18), 20)

plot1_18 <- VariableFeaturePlot(yuan18)
plot2_18 <- LabelPoints(plot = plot1_18, points = top20_18)
# plot1_18 + plot2_18

yuan18 <- ScaleData(yuan18, features = rownames(yuan18))

yuan18 <- RunPCA(yuan18, features = VariableFeatures(object=yuan18))

DimPlot(yuan18, reduction = "pca") + NoLegend()
ElbowPlot(yuan18, ndims=45)

d_18 <- 14
yuan18 <- FindNeighbors(yuan18, dims = 1:d_18)
yuan18 <- FindClusters(yuan18, resolution = 0.8)
#have to find optimal number of resolution

yuan18 <- RunUMAP(yuan18, dims=1:d_18)
DimPlot(yuan18, reduction = "umap")

all_markers_18 <- FindAllMarkers(yuan18, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# For each cluster, select the top 5 markers based on avg_log2FC
top5_markers_18 <- all_markers_18 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)%>%
  summarise(genes = paste(gene, collapse = ", "))
print(top5_markers_18)


DotPlot(
  yuan18,
  features = c("OLIG1", "TMEM144", "FLI1", "GFAP", "PTPRZ1", "RBFOX3", "CD163","P2RY12", "EGFR", "TOP2A"))

new.cluster.id.18 <- c("endothelial", "tumor", "oligodendrocyte", "oligodendrocyte", "dividing tumor", "astrocyte", "neuron", "pericyte", "macrophage")
names(new.cluster.id.18) <- levels(yuan18)
yuan18 <- RenameIdents(yuan18, new.cluster.id.18)
DimPlot(yuan18, reduction = "umap", label = TRUE)




# combined ----
anchors_1 <- FindIntegrationAnchors(object.list = list(yuan25, yuan30, yuan32, yuan35, yuan48, yuan16, yuan17, yuan18), dims=1:30)
combined <- IntegrateData(anchorset = anchors_1, dims=1:30)
combined <- ScaleData(combined)

combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)

#DimPlot(combined, reduction = "umap")
DimPlot(combined, reduction = "umap", group.by = "orig.ident")
DimPlot(combined, reduction = "umap")




