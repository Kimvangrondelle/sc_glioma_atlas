# Function to process and integrate multiple datasets, and rename clusters
# order is 25,30,32,35,48,16,17,18
library(Seurat)
library(dplyr)
j<- 1
process_datasets <- function(data_dirs, subset_params, dims_set, resolution = 0.8) {
  # Create an empty list to store Seurat objects
  seurat_objects <- list()
  
  # Loop through each dataset and apply the same pipeline
  for (i in data_dirs) {
    print(i)
    # Read 10X data
    data <- Read10X(data.dir = i)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = data, project = paste0("test_int_", i), min.cells=3, min.features=200)
    
    # Calculate mitochondrial gene percentage
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT")
    
    # Subset the data based on given parameters for each dataset
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
    seurat_obj <- RunTSNE(seurat_obj, dims = 1:dims_set[[j]])
    
    # Add the processed Seurat object to the list
    seurat_objects[[j]] <- seurat_obj
    j <- j + 1
    
  }
  return(seurat_objects)
  }  
# Example usage
data_dirs <- list(
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758474_PJ025",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758475_PJ030",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758476_PJ032", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758477_PJ035", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2940098_PJ048",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758471_PJ016",
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758472_PJ017", 
  "../../../mnt/neuro-genomic-1-ro/single_cell_data/GSE103224_Yuan/GSM2758473_PJ018"
)

# Define subset parameters for each dataset
subset_params <- list(
  list(nFeature_lower = 700, nFeature_upper = 3500, nCount_lower = 200, nCount_upper = 9000, percent_mt = 2),
  list(nFeature_lower = 700, nFeature_upper = 3600, nCount_lower = 200, nCount_upper = 5000, percent_mt = 2),
  list(nFeature_lower = 700, nFeature_upper = 4000, nCount_lower = 100, nCount_upper = 9000, percent_mt = 2),
  list(nFeature_lower = 700, nFeature_upper = 2200, nCount_lower = 200, nCount_upper = 3900, percent_mt = 2),
  list(nFeature_lower = 700, nFeature_upper = 3400, nCount_lower = 200, nCount_upper = 8000, percent_mt = 2),
  list(nFeature_lower = 1000, nFeature_upper = 6200, nCount_lower = 200, nCount_upper = 22000, percent_mt = 2),
  list(nFeature_lower = 700, nFeature_upper = 4200, nCount_lower = 200, nCount_upper = 12000, percent_mt = 2),
  list(nFeature_lower = 1000, nFeature_upper = 4800, nCount_lower = 200, nCount_upper = 12000, percent_mt = 2)# Add parameters for the other 6 datasets here
)

dims_set <- list(9, 11, 7, 9, 15, 11, 15, 14)

# # Call the function with 8 datasets
objects <- process_datasets(data_dirs, subset_params, dims_set)

