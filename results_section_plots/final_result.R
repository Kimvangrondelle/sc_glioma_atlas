df <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv", row.names = 1)
df <- df[, !colnames(df) %in% c("tau_score_vst","maximum","max_cluster","logmax")]
head(df)

filtered_genes <- rownames(df)[!grepl("^(RP11|RP3|RP13|RP4|RP5|RP1-)", rownames(df))]
df <- df[filtered_genes, ]

top_genes_per_celltype <- sapply(colnames(df), function(celltype) {
  top_genes <- rownames(df[order(-df[[celltype]]), ][1:20, , drop = FALSE])
  paste(celltype, ":", paste(top_genes, collapse = ", "))
})

writeLines(top_genes_per_celltype, "top_20_genes_per_celltype.txt")


markers <- read.csv("DESeq_output_combined_all.csv", row.names = 1)
filtered_genes <- rownames(markers)[!grepl("^(RP11|RP3|RP13|RP4|RP5|RP1-)", rownames(markers))]
markers <- markers[filtered_genes, ]

markers_by_cluster <- split(markers, markers$cluster)

top_100_genes_per_cluster <- list()

for (cluster in names(markers_by_cluster)) {
  cluster_markers <- markers_by_cluster[[cluster]]
  cluster_markers_sort <- cluster_markers[order(cluster_markers$p_val), ]
  top_100_genes_per_cluster[[cluster]] <- noquote(rownames(cluster_markers_sort[1:20, ]))
}

writeLines(top_100_genes_per_cluster, "top_20_genes_lowest_pval_deseq.txt")

print(top_100_genes_per_cluster)