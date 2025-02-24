#extract number of genes and cells per sample
seurat_object <- readRDS(file = "output/hijfte/3pr.rds")
num_cells <- ncol(seurat_object)
num_genes <- nrow(seurat_object)
# 
# 
cat("Number of genes:", num_genes, "\n")
cat("Number of cells:", num_cells, "\n")
# 
# extract number of clusters per cell type
num_clusters <- length(unique(seurat_object$seurat_clusters))

cat("Number of clusters:", num_clusters, "\n")
# 
# # DimPlot(seurat_object, reduction = "umap") + ggtitle("UMAP visualization of sample Diaz-11612")
# DimPlot(seurat_object, reduction = "umap", group.by = "celltype") + labs(title = "Annotated UMAP visualization of sample Hijfte-y ")
# 

#extract number of cell types and number of cells per cell type. 
cell_types <- seurat_object$celltype
cell_type_counts <- table(cell_types)

print(cell_type_counts)


#lists of celltypes per sample: Hijfte-y, diaz12017, diaz11612 and combined. 
c("GDiv.tumor", "Endo", "Per", "TAM", "GTumor", "Oligo", "Astro", "Neuron")
c("ADiv.tumor", "Oligo", "Per", "TAM", "ATumor")
c("O?", "Oligo", "OTumor")
c("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC", "sum_Tcell", "sum_Astro", "sum_Neuron", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure")

own <- read.csv("output/celltypes/hijfte/y/spec_prop.zscore_tau.csv", row.names = 1)
spread_scores <- own$tau_score_vst
#make long dataframe so that each gene appears as times as number of cell types present. -- to make scatterplot and color on cell type
long_own <- own %>%
  pivot_longer(cols = c("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC", "sum_Tcell", "sum_Astro", "sum_Neuron", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure"), 
               names_to = "CellType", 
               values_to = "SpecificityScore")
#scatterplot of tau vs specificity score ZPEX
ggplot(long_own, aes(x = tau_score_vst, y = SpecificityScore, color = CellType)) +
  geom_point(pch=19, alpha=0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  ggpubr::stat_cor(method= "spearman")+ #to get correlation per cell type
  # scale_y_continuous(limits = c(-250, 250)) +
  labs(title = "Scatterplot of Tau vs. Own Specificity Score by Cell Type for combined sample",
       x = "Tau score",
       y = "ZPEX score") +
  theme_minimal() 

#calculate overall spearman correlation between tau score and zpex score
own_scores <- long_own$SpecificityScore
tau_scores <- long_own$tau_score_vst

combined <- as.data.frame(cbind(own_scores, tau_scores))
combined <- combined %>% replace(is.na(.), 0)
cor(combined$own_scores, combined$tau_scores, method = "spearman")



bayes <- read.csv("output/celltypes/hijfte/y/spec_bayes.csv", row.names = 1)
bayes <- as.data.frame(-log10(as.matrix(bayes)))
bayes["tau_score_vst"] <- own$tau_score_vst
#make long dataframe so that each gene appears as times as number of cell types present. -- to make scatterplot and color on cell type
long_bayes <- bayes %>%
  pivot_longer(cols = c("GDiv.tumor", "Endo", "Per", "TAM", "GTumor", "Oligo", "Astro", "Neuron"), 
               names_to = "CellType", 
               values_to = "SpecificityScore")
#scatterplot of tau vs specificity score bayes
ggplot(long_bayes, aes(x = tau_score_vst, y = SpecificityScore, color = CellType)) +
  geom_point(pch=19, alpha=0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  ggpubr::stat_cor(method= "spearman")+ #to get correlation per cell type
  scale_x_continuous(limits = c(-0.15, 0.8))+
  labs(title = "Scatterplot of Tau vs. Bayes Specificity Score by Cell Type for sample hijfte-y",
       x = "Tau score",
       y = "-log10(Specificity Score Bayes)") +
  theme_minimal()

#calculate overall spearman correlation between tau score and bayes score
bayes_scores <- long_bayes$SpecificityScore
tau_scores <- long_bayes$tau_score_vst

combined <- as.data.frame(cbind(bayes_scores, tau_scores))
combined <- combined %>% replace(is.na(.), 0)

cor(combined$bayes_scores, combined$tau_scores, method = "spearman")

