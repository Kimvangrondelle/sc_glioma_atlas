# load in zpex data (combined) and only select the cell type columns so not the colmns "tau score vst etc)
df <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv", row.names = 1)
df <- df[, !colnames(df) %in% c("tau_score_vst","maximum","max_cluster","logmax")]
head(df)
#filter out non-coding genes
filtered_genes <- rownames(df)[!grepl("^(RP11|RP3|RP13|RP4|RP5|RP1-)", rownames(df))]
df <- df[filtered_genes, ]
#select top 20 per cell type
top_genes_per_celltype <- sapply(colnames(df), function(celltype) {
  top_genes <- rownames(df[order(-df[[celltype]]), ][1:20, , drop = FALSE])
  paste(celltype, ":", paste(top_genes, collapse = ", "))
})
#make txt file out of it. 
writeLines(top_genes_per_celltype, "top_20_genes_per_celltype.txt")

# load in DESeq2 output combined sample and filter out non-coding genes
markers <- read.csv("DESeq_output_combined_all.csv", row.names = 1)
filtered_genes <- rownames(markers)[!grepl("^(RP11|RP3|RP13|RP4|RP5|RP1-)", rownames(markers))]
markers <- markers[filtered_genes, ]

# split per cell type and select top 2 lowest p-val per cell type. 
markers_by_cluster <- split(markers, markers$cluster)

top_100_genes_per_cluster <- list()

for (cluster in names(markers_by_cluster)) {
  cluster_markers <- markers_by_cluster[[cluster]]
  cluster_markers_sort <- cluster_markers[order(cluster_markers$p_val), ]
  # top_100_genes_per_cluster[[cluster]] <- noquote((cluster_markers_sort$gene[1:20]))
  top_100_genes_per_cluster[[cluster]] <- cluster_markers_sort[cluster_markers_sort$p_val < 0.01, ]$gene
}

# writeLines(top_100_genes_per_cluster, "top_20_genes_lowest_pval_deseq.txt")

print(top_100_genes_per_cluster)
# significant_gene_counts <- sapply(top_100_genes_per_cluster, nrow)
# print(significant_gene_counts)