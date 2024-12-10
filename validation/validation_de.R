counts <- read.csv("counts_combined_summed.csv")
rownames(counts) <- counts$X
counts$X <- NULL
object <- CreateSeuratObject(counts)
Idents(object) <- rownames(object@assays$RNA@cells)

markers <- FindAllMarkers(object, test.use = "DESeq2", min.cells.group = 0, min.cells.feature = 0,
                          #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                          return.thresh = 1.1, logfc.threshold = 0.25, min.pct = 0.1) 

markers_by_cluster <- split(markers, markers$cluster)
top_100_genes_per_cluster <- list()

for (cluster in names(markers_by_cluster)) {
  cluster_markers <- markers_by_cluster[[cluster]]
  cluster_markers_sort <- cluster_markers[order(-cluster_markers$avg_log2FC), ]
  top_100_genes_per_cluster[[cluster]] <- cluster_markers_sort[1:200, ]
}

print(top_100_genes_per_cluster)



sample <- "yuan/48"  # for example: bolleboom/H243
file <- paste0("output/celltypes/", sample ,"/spec_prop.zscore_tau.csv")

specificity_scores_sample <- read.csv(file, row.names = 1)
specificity_scores_sample <- as.data.frame(specificity_scores_sample)
specificity_scores_sample <- subset(specificity_scores_sample, select= -c(tau_score_vst, maximum, max_cluster, logmax))

celltypes_in_sample <- list(colnames(specificity_scores_sample))

# sample_astro <- specificity_scores_sample[, "Astro", drop=FALSE]
# sample_astro <- head(sample_astro[order(-sample_astro$Astro), , drop = FALSE], 200)
# sample_oligo <- specificity_scores_sample[, "Oligo", drop=FALSE]
# sample_oligo <- head(sample_oligo[order(-sample_oligo$Oligo), , drop = FALSE], 200)
# sample_tam <- specificity_scores_sample[, "TAM", drop=FALSE]
# sample_tam <- head(sample_tam[order(-sample_tam$TAM), , drop = FALSE], 200)
# sample_neuron <- specificity_scores_sample[, "Neuron", drop=FALSE]
# sample_neuron <- head(sample_neuron[order(-sample_neuron$Neuron), , drop = FALSE], 200)
sample_tumor <- specificity_scores_sample[, "GTumor", drop=FALSE]
sample_tumor <- head(sample_tumor[order(-sample_tumor$GTumor), , drop = FALSE], 200)
sample_div <- specificity_scores_sample[, "GDiv.tumor", drop=FALSE]
sample_div <- head(sample_div[order(-sample_div$GDiv.tumor), , drop = FALSE], 200)
# sample_endo <- specificity_scores_sample[, "Endo", drop=FALSE]
# sample_endo <- head(sample_endo[order(-sample_endo$Endo), , drop = FALSE], 200)
# sample_per <- specificity_scores_sample[, "Per", drop=FALSE]
# sample_per <- head(sample_per[order(-sample_per$Per), , drop = FALSE], 200)
sample_opc <- specificity_scores_sample[, "OPC", drop=FALSE]
sample_opc <- head(sample_opc[order(-sample_opc$OPC), ], 200)
# sample_tcell <- specificity_scores_sample[, "Tcell", drop=FALSE]
# sample_tcell <- head(sample_tcell[order(-sample_tcell$Tcell), ], 200)

# print("Astro ") 
# print(length(intersect(top_100_genes_per_cluster$sum_Astro$gene, rownames(sample_astro))))
# print("Oligo ") 
# print(length(intersect(top_100_genes_per_cluster$sum_Oligo$gene, rownames(sample_oligo))))
# print("TAM ") 
# print(length(intersect(top_100_genes_per_cluster$sum_TAM$gene, rownames(sample_tam))))
# print("Neuron ") 
# print(length(intersect(top_100_genes_per_cluster$sum_Neuron$gene, rownames(sample_neuron))))
print("Tumor ") 
print(length(intersect(top_100_genes_per_cluster$sum_Tumor$gene, rownames(sample_tumor))))
print("Div ")
print(length(intersect(top_100_genes_per_cluster$sum_Div_tumor$gene, rownames(sample_div))))
# print("Endo ")
# print(length(intersect(top_100_genes_per_cluster$sum_Endo$gene, rownames(sample_endo))))
# print("Per ")
# print(length(intersect(top_100_genes_per_cluster$sum_Per$gene, rownames(sample_per))))
print("OPC ")
print(length(intersect(top_100_genes_per_cluster$sum_OPC$gene, rownames(sample_opc))))
# print("Tcell ")
# print(length(intersect(top_100_genes_per_cluster$sum_Tcell$gene, rownames(sample_tcell))))



