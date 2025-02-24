#partly used in final thesis
#checked DE genes per cell type from combined sample with marker genes found by McKenzie

counts <- read.csv("counts_combined_summed_all.csv")
rownames(counts) <- counts$X
counts$X <- NULL
counts <- as.data.frame(counts)
counts <- subset(counts, select= c(sum_Endo, sum_TAM, sum_Astro, sum_Oligo, sum_Neuron))
object <- CreateSeuratObject(counts)
Idents(object) <- rownames(object@assays$RNA@cells)




markers <- FindAllMarkers(object, test.use = "DESeq2", min.cells.group = 0, min.cells.feature = 0,
                          #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                          return.thresh = 1.1, logfc.threshold = 0.25, min.pct = 0.1) 
markers <- read.csv("output/celltypes/combined_all/DESeq_output.csv")

markers_by_cluster <- split(markers, markers$cluster)
top_100_genes_per_cluster <- list()

for (cluster in names(markers_by_cluster)) {
  cluster_markers <- markers_by_cluster[[cluster]]
  cluster_markers_sort <- cluster_markers[order(-cluster_markers$avg_log2FC), ]
  top_100_genes_per_cluster[[cluster]] <- cluster_markers_sort[1:500, ]
}

print(top_100_genes_per_cluster)

#load in data from McKenzie
kenzie <- read_excel("validation/41598_2018_27293_MOESM2_ESM.xlsx", sheet = "top_human_specificity")
kenzie <- as.data.frame(kenzie)
colnames(kenzie) <- kenzie[2, ]
kenzie <- kenzie[-2, ]

celltypes_kenzie <- split(kenzie, kenzie$Celltype)

#extract and sort data on grand_mean
kenzie_ast <- celltypes_kenzie$ast
rownames(kenzie_ast) <- kenzie_ast$gene 
kenzie_ast$gene <- NULL 
kenzie_ast$grand_mean <- as.numeric(kenzie_ast$grand_mean)
kenzie_ast <- kenzie_ast[order(-kenzie_ast$grand_mean), ]

kenzie_oli <- celltypes_kenzie$oli
rownames(kenzie_oli) <- kenzie_oli$gene 
kenzie_oli$gene <- NULL 
kenzie_oli$grand_mean <- as.numeric(kenzie_oli$grand_mean)
kenzie_oli <- kenzie_oli[order(-kenzie_oli$grand_mean), ]


kenzie_neu <- celltypes_kenzie$neu
rownames(kenzie_neu) <- kenzie_neu$gene 
kenzie_neu$gene <- NULL 
kenzie_neu$grand_mean <- as.numeric(kenzie_neu$grand_mean)
kenzie_neu <- kenzie_neu[order(-kenzie_neu$grand_mean), ]

kenzie_mic <- celltypes_kenzie$mic
rownames(kenzie_mic) <- kenzie_mic$gene 
kenzie_mic$gene <- NULL 
kenzie_mic$grand_mean <- as.numeric(kenzie_mic$grand_mean)
kenzie_mic <- kenzie_mic[order(-kenzie_mic$grand_mean), ]

kenzie_end <- celltypes_kenzie$end
rownames(kenzie_end) <- kenzie_end$gene 
kenzie_end$gene <- NULL 
kenzie_end$grand_mean <- as.numeric(kenzie_end$grand_mean)
kenzie_end <- kenzie_end[order(-kenzie_end$grand_mean), ]

#intersect top 100 DE genes with the top n markers from McKenzie
print(length(intersect(top_100_genes_per_cluster$sum_Astro$gene, head(rownames(kenzie_ast), 500))))
print(length(intersect(top_100_genes_per_cluster$sum_Oligo$gene, head(rownames(kenzie_oli), 500))))
print(length(intersect(top_100_genes_per_cluster$sum_Neuron$gene, head(rownames(kenzie_neu), 500))))
print(length(intersect(top_100_genes_per_cluster$sum_TAM$gene, head(rownames(kenzie_mic), 500))))
print(length(intersect(top_100_genes_per_cluster$sum_Endo$gene, head(rownames(kenzie_end), 500))))