j<- 1
process_datasets <- function(data_dirs, subset_params, dims_set, resolution = 0.8) {
  seurat_objects <- list()
  
  # Loop through each dataset and apply the same pipeline
  for (i in data_dirs) {
    print(i)
    data <- Read10X(data.dir = i)
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
    
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims_set[[j]])
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    seurat_obj <- RunUMAP(seurat_obj, dims=1:dims_set[[j]])
    
    seurat_objects[[j]] <- seurat_obj
    j <- j + 1
    
  }
  return(seurat_objects)
}  
# Example usage
data_dirs <- list(
  "../../mnt/neuro-genomic-1-ro/single_cell_data/Bolleboom-Gao/H243_GBM/H243_filtered_feature_bc_matrix"
)

# Define subset parameters for each dataset
subset_params <- list(
  list(nFeature_lower = 200, nFeature_upper = 6000, nCount_lower = 200, nCount_upper = 18000, percent_mt = 1)
)

dims_set <- list(16)

# # Call the function with 8 datasets
objects <- process_datasets(data_dirs, subset_params, dims_set)

