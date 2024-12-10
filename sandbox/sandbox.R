diaz10022 <- Read10X(data.dir = "../../mnt/neuro-genomic-1-ro/single_cell_data/GSE138794_Diaz/snRNA_GSM4119521_SF10022_GBM")
diaz10022 <- CreateSeuratObject(counts = diaz10022, project = "sc_glioma_atlas", min.cells=3, min.features=200)
diaz10022[["percent.mt"]] <- PercentageFeatureSet(diaz10022, pattern = "^MT")
VlnPlot(diaz10022, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

diaz10022 <- subset(x=diaz10022, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & nCount_RNA > 4000 & nCount_RNA < 25000 & percent.mt < 1.5)
diaz10022 <- NormalizeData(diaz10022)
diaz10022 <- FindVariableFeatures(diaz10022, selection.method = "vst", nfeatures = 2000)
diaz10022 <- ScaleData(diaz10022, features = rownames(diaz10022))
diaz10022 <- RunPCA(diaz10022, features = VariableFeatures(object = diaz10022))

DimPlot(diaz10022, reduction = "pca") 
ElbowPlot(diaz10022, ndims = 45)

d_10022 <- 15
diaz10022 <- FindNeighbors(diaz10022, dims = 1:d_10022)
diaz10022 <- FindClusters(diaz10022, resolution = 0.8)
diaz10022 <- RunUMAP(diaz10022, dims = 1:d_10022)
DimPlot(diaz10022, reduction = "umap")

markers_10022 <- FindAllMarkers(diaz10022, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers_10022 %>%
  group_by(cluster) %>%
  top_n(n=-10, wt= p_val) %>%
  summarise(genes = paste(gene, collapse = ", "))

filtered_genes <- markers_10022 %>%
  filter(p_val < 0.05) %>%              # Filter for adjusted p-values < 0.05
  arrange(cluster, desc(avg_log2FC))         # Sort by cluster and log fold change (descending)

# Now combine the top 10 genes for each cluster
top10 <- filtered_genes %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p_val) %>%            # Select top 10 genes with the lowest p-values
  summarise(genes = paste(gene, collapse = ", "))

top10

cluster_counts <- table(diaz10022$seurat_clusters)

# Zet de resultaten om naar een data.frame voor beter overzicht
cluster_counts_df <- as.data.frame(cluster_counts)



object <- readRDS(file = "output/hijfte/hijfte-y.rds")
agg_object <- AggregateExpression(object, group.by = "seurat_clusters", return.seurat = TRUE, normalization.method = "LogNormalize")

cells_per_cluster <- table(Idents(object))
counts <- agg_object@assays$RNA@layers$counts
rownames(counts) <- Features(agg_object)
counts_vst <- counts %>%
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>%
  DESeq2::vst(blind=T)
vst_counts <- as.data.frame(counts_vst@assays@data@listData)
tau_score_vst <- apply(vst_counts, 1, calc_tau)
counts_tau_vst <- data.frame(vst_counts, tau_score_vst)

z_vst <- t(apply(vst_counts, 1, z_score))
prop_vst <- t(apply(vst_counts, 1, proportion))
spec <- round(prop_vst * z_vst, 3)
spec_tau <- data.frame(spec, tau_score_vst)

spec_tau$maximum <- apply(spec_tau[ , !(colnames(spec_tau) %in% "tau_score_vst")], 1, max)
spec_tau$logmax <- log(spec_tau$maximum)



# plotlog <- ggplot(sorted_spec, aes(x=tau_score_vst, y=logmax)) +
#   geom_point() + # Show dots
#   geom_text(
#     label=rownames(sorted_spec), 
#     nudge_x = 0, nudge_y = 0, 
#     check_overlap = T
#   )
library(plotly)

fig <- plot_ly(type = 'scatter', mode = 'markers') 
fig <- fig %>%
  add_trace(
    x = spec_tau$tau_score_vst, 
    y = spec_tau$logmax,
    text = rownames(spec_tau),
    hoverinfo = 'text',
    marker = list(color='green'),
    showlegend = F
  )

top10_genes <- unlist(strsplit(top10$genes, ", "))

# Modify the plot
plotlog <- ggplot(spec_tau, aes(x=tau_score_vst, y=logmax)) +
  # Color points based on whether the gene is in top10 or not
  geom_point(aes(color = ifelse(rownames(spec_tau) %in% top10_genes, "top10", "other"))) +
  
  # Add labels with nudging
  # geom_text(
  #   label = rownames(spec_tau),
  #   nudge_x = 0, nudge_y = 0,
  #   check_overlap = TRUE
  # ) +
  
  # Define custom colors for top10 and other genes
  scale_color_manual(values = c("top10" = "red", "other" = "grey")) +
  
  # Add labels and theme
  labs(color = "Gene Category") +
  theme_minimal()

# Display the plot
print(plotlog)


diaz10022 <- readRDS(file = "output/diaz/diaz10022.rds")
dim(diaz10022)

pseudo_diaz10022 <- AggregateExpression(diaz10022, group.by = "seurat_clusters", return.seurat = TRUE, normalization.method = "LogNormalize")

dim(pseudo_diaz10022)

diaz_10022_counts <- pseudo_diaz10022@assays$RNA@layers$counts
rownames(diaz_10022_counts) <- Features(pseudo_diaz10022)
cells_per_cluster <- table(Idents(diaz10022))
norm_counts <- sweep(diaz_10022_counts, 2, cells_per_cluster, FUN = "/")

vst_counts <- vst(as.matrix(diaz_10022_counts))
vst_counts
tau_score_vst <- apply(vst_counts, 1, calc_tau)
tau_score_vst
hist(tau_score_vst)
diaz_10022_vst <- data.frame(vst_counts, tau_score_vst)
diaz_10022_vst
sorted_diaz <- diaz_10022_vst[order(diaz_10022_vst$tau_score_vst, decreasing = TRUE),]



cv_score_vst <- t(apply(vst_counts, 1, cv_new))
dim(cv_score_vst)
prop_vst <- t(apply(vst_counts, 1, function(x) x / sum(x)))
dim(prop_vst)
prop_cv_vst <- round(prop_vst * cv_score_vst, 3)
prop_cv_vst
prop_10022_vst <- data.frame(prop_cv_vst, tau_score_vst)
sorted_prop_10022 <- prop_10022_vst[order(prop_10022_vst$tau_score_vst, decreasing = TRUE),]

sorted_prop_10022$maximum <- apply(sorted_prop_10022, 1, max)
sorted_prop_10022

ggplot(sorted_prop_10022, aes(x=tau_score_vst, y=maximum)) +
  geom_point() + # Show dots
  geom_text(
    label=rownames(sorted_prop_10022), 
    nudge_x = 0, nudge_y = 0, 
    check_overlap = T
  )

library(plotly)

fig <- plot_ly(type = 'scatter', mode = 'markers') 
fig <- fig %>%
  add_trace(
    x = sorted_prop_10022$tau_score_vst, 
    y = sorted_prop_10022$maximum,
    text = rownames(sorted_prop_10022),
    hoverinfo = 'text',
    marker = list(color='green'),
    showlegend = F
  )

fig
# vst_subset_1_7_10_13 <- as.data.frame(vst_counts[, c("1", "7", "10", "13")]) 
# cv_subset <- t(apply(vst_subset_1_7_10_13, 1, cv_new))
# prop_subset <- t(apply(vst_subset_1_7_10_13, 1, function(x) x/sum(x)))
# prop_cv_subset <- prop_subset * cv_subset
# vst_subset_1_7_10_13$tau <- apply(vst_subset_1_7_10_13, 1, calc_tau)
# vst_subset_1_7_10_13
# vst_subset_1_7_10_13 <- vst_subset_1_7_10_13[order(vst_subset_1_7_10_13$tau, decreasing = TRUE),]

# prop_cv_subset$tau <- t(apply(vst_subset_1_7_10_13, 1, calc_tau))
# sorted_prop_cv <- prop_cv_subset[order(prop_cv_subset$tau),]

pseudo_markers_diaz10022 <- FindAllMarkers(pseudo_diaz10022, test.use = "DESeq2", min.cells.group = 0, min.cells.feature = 0, 
                                           #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                           return.thresh = 1.1, logfc.threshold = 0, min.pct = 0)


pseudo_markers_diaz10022 %>%  dplyr::group_by(cluster) %>% dplyr::select(gene) %>% dplyr::tally()


head(pseudo_markers_diaz10022)
dim(pseudo_markers_diaz10022)


top10 <- pseudo_markers_diaz10022 %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = p_val) %>%
  summarise(genes = paste(gene, collapse = ", "))
top10

cv_score <- t(apply(diaz_10022_counts, 1, cv_new))
dim(cv_score)
prop <- t(apply(diaz_10022_counts, 1, function(x) x / sum(x)))
dim(prop)
prop_cv <- prop * cv_score
prop_cv

tau_score_diaz <- apply(diaz_10022_counts, 1, calc_tau)
diaz_10022_specificity <- data.frame(prop_cv, tau_score_diaz)
diaz_10022_specificity

diaz_10022_counts

print(rowSums(diaz_10022_counts) > 40000)



set.seed(123)  # For reproducibility

# Simulate counts for 100 genes across 6 samples (clusters)
vst_test <- data.frame(
  row.names = paste0("Gene", 1:100),  # 20 genes
  Cluster1 = rpois(100, lambda = 500),  # Random counts with mean = 500
  Cluster2 = rpois(100, lambda = 45),  # Random counts with mean = 45
  Cluster3 = rpois(100, lambda = 55),  # Random counts with mean = 55
  Cluster4 = rpois(100, lambda = 60),  # Random counts with mean = 60
  Cluster5 = rpois(100, lambda = 52),  # Random counts with mean = 52
  Cluster6 = rpois(100, lambda = 48)   # Random counts with mean = 48
)

# Display the first few rows of the simulated dataset
vst_test

# Create the condition (metadata) for the samples (clusters)
colData <- data.frame(
  cond = as.factor(c("c1", "c1", "c2", "c2", "c3", "c3"))
)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(vst_test),
  colData = colData,
  design = ~cond
)

# Perform the variance stabilizing transformation (VST)
vst_counts <- varianceStabilizingTransformation(dds, blind = TRUE)

# Check the transformed counts
vst_transformed <- assay(vst_counts)
vst_transformed

counts.vst <- vst_test %>%
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>%
  DESeq2::vst(blind=T, nsub = 50)

counts.vst@assays@data@listData
head(vst_test)
head(counts.vst@assays@data@listData)