#testing the merging of two samples

object_1 <- Read10X("../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_1of2.filtered_gene_matrices")
object_1 <- CreateSeuratObject(counts = object_1, project = paste0("int"), min.cells=3, min.features=200)

object_2 <- Read10X("../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_2of2.filtered_gene_matrices")
object_2 <- CreateSeuratObject(counts = object_2, project = paste0("int"), min.cells=3, min.features=200)

objects.combined <- merge(object_1, y = object_2, add.cell.ids = c("1/2", "2/2"), project = "coutu338")
objects.combined
seurat_obj <- objects.combined
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT")



seurat_obj <- subset(x = seurat_obj, 
                     subset = nFeature_RNA > 1000 & 
                       nFeature_RNA < 7000 & 
                       nCount_RNA > 200 & 
                       nCount_RNA < 35000 & 
                       percent.mt < 2)

#checking whether cells were present in both samples
print("1/2_AAACGGGTCGTCGTTC-1" %in% rownames(seurat_obj@assays$RNA@cells@.Data))
print("2/2_AAACGGGTCGTCGTTC-1" %in% rownames(seurat_obj@assays$RNA@cells@.Data))
print("1/2_AGGCCGTTCTAAGCCA-1" %in% rownames(seurat_obj@assays$RNA@cells@.Data))
print("2/2_AGGCCGTTCTAAGCCA-1" %in% rownames(seurat_obj@assays$RNA@cells@.Data))
print("1/2_CCCATACCACGCGAAA-1" %in% rownames(seurat_obj@assays$RNA@cells@.Data))
print("2/2_CCCATACCACGCGAAA-1" %in% rownames(seurat_obj@assays$RNA@cells@.Data))


