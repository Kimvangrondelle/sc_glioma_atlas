counts <- read.csv("counts_combined_summed_all.csv")
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



sample <- "hijfte/pr3"  # for example: bolleboom/H243
file <- paste0("output/celltypes/", sample ,"/spec_prop.zscore_tau.csv")

specificity_scores_sample <- read.csv(file, row.names = 1)
specificity_scores_sample <- as.data.frame(specificity_scores_sample)
specificity_scores_sample <- subset(specificity_scores_sample, select= -c(tau_score_vst, maximum, max_cluster, logmax))

celltypes_in_sample <- list(colnames(specificity_scores_sample))

sample_astro <- specificity_scores_sample[, "Astro", drop=FALSE]
sample_astro <- head(sample_astro[order(-sample_astro$Astro), , drop = FALSE], 200)
sample_oligo <- specificity_scores_sample[, "Oligo", drop=FALSE]
sample_oligo <- head(sample_oligo[order(-sample_oligo$Oligo), , drop = FALSE], 200)
sample_tam <- specificity_scores_sample[, "TAM", drop=FALSE]
sample_tam <- head(sample_tam[order(-sample_tam$TAM), , drop = FALSE], 200)
sample_neuron <- specificity_scores_sample[, "Neuron", drop=FALSE]
sample_neuron <- head(sample_neuron[order(-sample_neuron$Neuron), , drop = FALSE], 200)
sample_gtumor <- specificity_scores_sample[, "GTumor", drop=FALSE]
sample_gtumor <- head(sample_gtumor[order(-sample_gtumor$GTumor), , drop = FALSE], 200)
sample_gdiv <- specificity_scores_sample[, "GDiv.tumor", drop=FALSE]
sample_gdiv <- head(sample_gdiv[order(-sample_gdiv$GDiv.tumor), , drop = FALSE], 200)
sample_endo <- specificity_scores_sample[, "Endo", drop=FALSE]
sample_endo <- head(sample_endo[order(-sample_endo$Endo), , drop = FALSE], 200)
sample_per <- specificity_scores_sample[, "Per", drop=FALSE]
sample_per <- head(sample_per[order(-sample_per$Per), , drop = FALSE], 200)
sample_opc <- specificity_scores_sample[, "OPC", drop=FALSE]
sample_opc <- head(sample_opc[order(-sample_opc$OPC), ], 200)
sample_tcell <- specificity_scores_sample[, "Tcell", drop=FALSE]
sample_tcell <- head(sample_tcell[order(-sample_tcell$Tcell), ], 200)
sample_atumor <- specificity_scores_sample[, "ATumor", drop=FALSE]
sample_atumor <- head(sample_atumor[order(-sample_atumor$ATumor), ], 200)
sample_otumor <- specificity_scores_sample[, "OTumor", drop=FALSE]
sample_otumor <- head(sample_otumor[order(-sample_otumor$OTumor), ], 200)
sample_adiv <- specificity_scores_sample[, "ADiv.tumor", drop=FALSE]
sample_adiv <- head(sample_adiv[order(-sample_adiv$ADiv.tumor), ], 200)
sample_onot <- specificity_scores_sample[, "O.", drop=FALSE]
sample_onot <- head(sample_onot[order(-sample_onot$O.), ], 200)

print("Astro ")
print(length(intersect(top_100_genes_per_cluster$sum_Astro$gene, rownames(sample_astro))))
print("Oligo ")
print(length(intersect(top_100_genes_per_cluster$sum_Oligo$gene, rownames(sample_oligo))))
print("TAM ")
print(length(intersect(top_100_genes_per_cluster$sum_TAM$gene, rownames(sample_tam))))
print("Neuron ")
print(length(intersect(top_100_genes_per_cluster$sum_Neuron$gene, rownames(sample_neuron))))
print("Tumor ") 
print(length(intersect(top_100_genes_per_cluster$sum_GTumor$gene, rownames(sample_gtumor))))
print("Div ")
print(length(intersect(top_100_genes_per_cluster$sum_GDiv_tumor$gene, rownames(sample_gdiv))))
print("Endo ")
print(length(intersect(top_100_genes_per_cluster$sum_Endo$gene, rownames(sample_endo))))
print("Per ")
print(length(intersect(top_100_genes_per_cluster$sum_Per$gene, rownames(sample_per))))
print("OPC ")
print(length(intersect(top_100_genes_per_cluster$sum_OPC$gene, rownames(sample_opc))))
print("Tcell ")
print(length(intersect(top_100_genes_per_cluster$sum_Tcell$gene, rownames(sample_tcell))))
print("ATumor")
print(length(intersect(top_100_genes_per_cluster$sum_ATumor$gene, rownames(sample_atumor))))
print("OTumor")
print(length(intersect(top_100_genes_per_cluster$sum_OTumor$gene, rownames(sample_otumor))))
print("ADiv")
print(length(intersect(top_100_genes_per_cluster$sum_ADiv_tumor$gene, rownames(sample_adiv))))
print("O.")
print(length(intersect(top_100_genes_per_cluster$sum_O_notsure$gene, rownames(sample_onot))))




