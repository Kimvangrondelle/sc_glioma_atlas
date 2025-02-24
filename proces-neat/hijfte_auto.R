# Function to process and integrate multiple datasets, and rename clusters

library(Seurat)
library(dplyr)
j <- 1
process_datasets <- function(data_dirs, subset_params, dims_set, resolution = 0.8) {
  # Create an empty list to store Seurat objects
  seurat_objects <- list()
  
  # Loop through each dataset and apply the same pipeline
  for (i in data_dirs) {
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
    #clustering and visualization
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
  "../../mnt/neuro-genomic-1-ro/Glimmunology/LGG_project/10x_singlenuc_RNAseq/Optimisation_10x_snRNAseq/Glioma_Y_and_O/Levi2_Glioma_Y/outs/raw_feature_bc_matrix"
)

# Define subset parameters for each dataset
subset_params <- list(
  list(nFeature_lower = 1400, nFeature_upper = 4500, nCount_lower = 2200, nCount_upper = 14000, percent_mt = 2.5)# Add parameters for the other 6 datasets here
)
# optimal number of principal components
dims_set <- list(16)

# Call the function 
objects <- process_datasets(data_dirs, subset_params, dims_set)


