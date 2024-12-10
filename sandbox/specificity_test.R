
gene_exp_matr_sample_clus <- function(object) { #get matrix with average gene expres per cluster
  seurat_obj <- readRDS(file = paste0("output/", object))
  #get cluster ids --- later this will be the cell types. 
  cluster_ids <- unique(seurat_obj@meta.data$seurat_clusters)
  average_expression_list <- list()
  for (j in cluster_ids) {
    print(j)
    #subset cluster
    cluster_subset <- subset(seurat_obj, idents = j)
    #get expression matrix for the cluster
    gene_expression_matrix <- GetAssayData(cluster_subset, slot = "data")
    #get average expression over all cells in the cluster
    average_expression <- rowMeans(gene_expression_matrix)
    average_expression_list[[as.character(j)]] <- average_expression
  }
    # bind the average expression of all clusters into one matrix
  average_expression_matrix <- do.call(cbind, average_expression_list)
  # set colnames with the corresponding celtype names 
  colnames(average_expression_matrix) <- paste0("Cluster_", cluster_ids)  
  return(average_expression_matrix)
}

mat16 <- gene_exp_matr_sample_clus("yuan16.rds")


tau <- function(average_exp) {  # 1 is cluster ex == max ex
  n_genes <- nrow(average_exp)
  n_clus <- ncol(average_exp)
  
  tau_mat <- matrix(0, nrow = n_genes, ncol = n_clus)
  
  rownames(tau_mat) <- rownames(average_exp)
  colnames(tau_mat) <- colnames(average_exp)
  
  for (gene in 1:n_genes) {
    gene_ex <- average_exp[gene, ]
    #get the max average expression for a gene
    max_ex <- max(gene_ex)
    for (cluster in 1:n_clus) {
      if (max_ex > 0) { # make sure no division zero
        #divide every expression by the max_ex
        norm_ex <- (gene_ex[cluster] / max_ex)
        tau_mat[gene, cluster] <- norm_ex
      } else {
        tau_mat[gene, cluster] <- 0
      }
    }
  }
  return(tau_mat)
}

tau_val_mat <- tau(mat16)
tau_val_df <- as.data.frame(tau_val_mat)


get_cluster_higest_tau <- function(row) {
  max_tau <- max(row)
  max_col <- names(row)[which.max(row)]
  return(c(Max_Tau = max_tau, Cluster = max_col))
}

res <- t(apply(tau_val_df, 1, get_cluster_higest_tau))
res_df <- as.data.frame(res)
#res_df <- res_df %>% arrange(Cluster)
print(res_df)