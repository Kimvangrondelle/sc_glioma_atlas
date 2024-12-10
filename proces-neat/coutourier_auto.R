j <- 1
process_datasets <- function(data_dirs, subset_params, dims_set, resolution = 0.8) {
  # Create an empty list to store Seurat objects
  seurat_objects <- list()
  
  for (i in list(1, 3, 5, 7)) {
    print(i)
    object_1 <- Read10X(data.dir = data_dirs[[i]])
    object_1 <- CreateSeuratObject(counts = object_1, project = paste0("int_", i), min.cells=3, min.features=200)
    
    object_2 <- Read10X(data.dir = data_dirs[[i+1]])
    object_2 <- CreateSeuratObject(counts = object_2, project = paste0("int_", i+1), min.cells=3, min.features=200)
    
    objects.combined <- merge(object_1, y = object_2, add.cell.ids = c("1/2", "2/2"), project = "coutu338")
    objects.combined
    seurat_obj <- objects.combined
    seurat_obj <- JoinLayers(seurat_obj)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT")
    
    print(j)
    subset_params_i <- subset_params[[j]]
    
    seurat_obj <- subset(x = seurat_obj, 
                         subset = nFeature_RNA > subset_params_i$nFeature_lower & 
                           nFeature_RNA < subset_params_i$nFeature_upper & 
                           nCount_RNA > subset_params_i$nCount_lower & 
                           nCount_RNA < subset_params_i$nCount_upper & 
                           percent.mt < subset_params_i$percent_mt)
    
    # Normalize data and find variable features
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features=all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object=seurat_obj))
    
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims_set[[j]])
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    
    seurat_obj <- RunUMAP(seurat_obj, dims=1:dims_set[[j]])
    
    # Add the processed Seurat object to the list
    seurat_objects[[j]] <- seurat_obj
    j <- j + 1
    
  }
  return(seurat_objects)
}

data_dirs <- list(
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT338_1of2.filtered_gene_matrices",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT338_2of2.filtered_gene_matrices",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT363_1of2.filtered_gene_matrices", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT363_2of2.filtered_gene_matrices", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_1of2.filtered_gene_matrices",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT364_2of2.filtered_gene_matrices",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT397_1of2.filtered_gene_matrices", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/EGAS00001004422_Couturier/BT397_2of2.filtered_gene_matrices"
)

subset_params <- list(
  # waardes kloppen nog niet
  list(nFeature_lower = 1000, nFeature_upper = 7800, nCount_lower = 200, nCount_upper = 50000, percent_mt = 2),
  list(nFeature_lower = 1000, nFeature_upper = 5500, nCount_lower = 200, nCount_upper = 20000, percent_mt = 2),
  list(nFeature_lower = 1000, nFeature_upper = 7000, nCount_lower = 200, nCount_upper = 35000, percent_mt = 2),
  list(nFeature_lower = 700, nFeature_upper = 4300, nCount_lower = 200, nCount_upper = 20000, percent_mt = 2)
  )

dims_set <- list(12, 12, 14, 13)


objects <- process_datasets(data_dirs, subset_params, dims_set)

