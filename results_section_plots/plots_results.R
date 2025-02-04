seurat_object <- readRDS(file = "output/hijfte/3pr.rds")
num_cells <- ncol(seurat_object)
num_genes <- nrow(seurat_object)
# 
# 
cat("Number of genes:", num_genes, "\n")
cat("Number of cells:", num_cells, "\n")
# 
# 
num_clusters <- length(unique(seurat_object$seurat_clusters))

cat("Number of clusters:", num_clusters, "\n")
# 
# # DimPlot(seurat_object, reduction = "umap") + ggtitle("UMAP visualization of sample Diaz-11612")
# DimPlot(seurat_object, reduction = "umap", group.by = "celltype") + labs(title = "Annotated UMAP visualization of sample Hijfte-y ")
# 


cell_types <- seurat_object$celltype
cell_type_counts <- table(cell_types)

print(cell_type_counts)


library(pheatmap)

data <- read.csv("output/celltypes/hijfte/y/spec_prop.zscore_tau.csv", row.names = 1)
spread_scores <- data$tau_score_vst

specificity_matrix <- data[, !colnames(data) %in% c("tau_score_vst","maximum","max_cluster","logmax")]

# sorted_indices <- order(spread_scores, decreasing = TRUE)
# sorted_specificity_matrix <- specificity_matrix[sorted_indices, ]
# sorted_spread_scores <- spread_scores[sorted_indices]
# 
# annotation <- data.frame(tau = sorted_spread_scores)
# rownames(annotation) <- rownames(sorted_specificity_matrix)
# 
# pheatmap(sorted_specificity_matrix, annotation_row = annotation, cluster_rows = FALSE, scale = "row", main = "Gene specificity Heatmap", color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = FALSE)


own <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv", row.names = 1)
spread_scores <- own$tau_score_vst

long_own <- own %>%
  pivot_longer(cols = c("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC", "sum_Tcell", "sum_Astro", "sum_Neuron", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure"), 
               names_to = "CellType", 
               values_to = "SpecificityScore")

ggplot(long_own, aes(x = tau_score_vst, y = SpecificityScore, color = CellType)) +
  geom_point(pch=19, alpha=0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Scatterplot of Tau vs. Own Specificity Score by Cell Type for all samples combined",
       x = "Tau score",
       y = "Specificity Score") +
  theme_minimal() 

bayes <- read.csv("output/celltypes/combined_all/spec_bayes.csv", row.names = 1)
bayes["tau_score_vst"] <- own$tau_score_vst

long_bayes <- bayes %>%
  pivot_longer(cols = c("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC", "sum_Tcell", "sum_Astro", "sum_Neuron", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure"), 
               names_to = "CellType", 
               values_to = "SpecificityScore")

ggplot(long_bayes, aes(x = tau_score_vst, y = SpecificityScore, color = CellType)) +
  geom_point(pch=19, alpha=0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Scatterplot of Tau vs. Bayes Specificity Score by Cell Type for all samples combined",
       x = "Tau score",
       y = "Specificity Score") +
  theme_minimal()
