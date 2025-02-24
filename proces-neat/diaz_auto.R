j<- 1
process_datasets <- function(data_dirs, subset_params, dims_set, resolution = 0.8) {
  # Create an empty list to store Seurat objects
  seurat_objects <- list()
  
  # Loop through each dataset and apply the same pipeline
  for (i in data_dirs) {
    print(i)
    if (i == "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119534_SF11644_GBM") {
      data <- Read10X(data.dir = i, gene.column = 1)
    } else {
      data <- Read10X(data.dir = i)
    }
    seurat_obj <- CreateSeuratObject(counts = data, project = paste0("test_int_", i), min.cells=3, min.features=200)
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
    #clustering
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims_set[[j]])
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    seurat_obj <- RunUMAP(seurat_obj, dims=1:dims_set[[j]])
    
    seurat_objects[[j]] <- seurat_obj
    j <- j + 1
    
  }
  return(seurat_objects)
}  

data_dirs <- list(
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119521_SF10022_GBM",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119525_SF12264_GBM",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119527_SF4297_GBM", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4119534_SF11644_GBM", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/scRNA_GSM4658373_SF11681_GBM"
)

# Define subset parameters for each dataset
subset_params <- list(
  list(nFeature_lower = 2000, nFeature_upper = 6000, nCount_lower = 4000, nCount_upper = 25000, percent_mt = 1.5),
  list(nFeature_lower = 200, nFeature_upper = 10000, nCount_lower = 500, nCount_upper = 50000, percent_mt = 20),
  list(nFeature_lower = 200, nFeature_upper = 5500, nCount_lower = 500, nCount_upper = 15000, percent_mt = 1),
  list(nFeature_lower = 1200, nFeature_upper = 8000, nCount_lower = 200, nCount_upper = 45000, percent_mt = 5),
  list(nFeature_lower = 500, nFeature_upper = 7000, nCount_lower = 200, nCount_upper = 38000, percent_mt = 10)
)
#optimal number of principal components
dims_set <- list(15, 12, 13, 13, 14)

# Call the function with 5 datasets
objects <- process_datasets(data_dirs, subset_params, dims_set)

