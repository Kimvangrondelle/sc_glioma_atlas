#barplots van elk celltype met zelfde volgorde als in plot naast de plot

# corrplot: 
# https://github.com/yhoogstrate/recursiveCorPlot/blob/master/R/function_recursiveCorPlot.R

#ggcorplot
# https://github.com/yhoogstrate/glass-od/blob/main/scripts/func_ggcorrplot.R

partially <- readRDS("wies_glass_validatie/dge.partially.paired.clusters.Rds")

hclust <- readRDS("wies_glass_validatie/glass_nl__hclust.Rds")

h <- readRDS("wies_glass_validatie/h.Rds")



bulk <- hclust |>
  head(n=200) |>
  t() |>
  cor()



a = corrplot::corrplot(bulk, order="hclust",  tl.cex=0.25, tl.pos="l")
corrplot <- corrplot::corrplot(bulk, order="hclust",  tl.cex=0.25, tl.pos="l")
a$corrPos$y
a$corr |> rownames()
rownames(a$corr)

png("corrplot.png", width = 1000, height = 1000)
par(mar = c(0, 0, 0, 0))
corrplot::corrplot(bulk, order = "hclust", tl.cex = 0.25, tl.pos = "l")
dev.off()


genes_corplot <- sapply(strsplit(rownames(a$corr), "_"), '[', 2)

# specificity <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv")
# specificity <- read.csv("output/celltypes/hijfte/y/spec_prop.zscore_tau.csv")
specificity <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv")
specificity <- as.data.frame(specificity)
specificity <- subset(specificity, select= -c(tau_score_vst, maximum, max_cluster, logmax))


genes_with_scores <- data.frame(Gene = genes_corplot) %>% 
  left_join(specificity, by = c("Gene" = "X"))

genes_with_scores[is.na(genes_with_scores)] <- 0
genes_with_scores$Gene <- factor(genes_with_scores$Gene, levels = genes_corplot)


barplot_list <- list()
# celltypes <- colnames(specificity)[-1]
celltypes <- list("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC")
celltypes <- list("sum_Astro", "sum_Tcell", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure", "sum_Neuron")


for (celltype in celltypes) {
  p <- ggplot(genes_with_scores, aes_string(x="Gene", y=celltype)) +
    geom_col(fill = "skyblue") +
    ggtitle(paste(celltype)) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", angle = 90, size=10)) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Specificity \n score") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(genes_with_scores$Gene))) +   # Reverse the order of the y-axis
    theme(axis.text.y = element_blank())
  barplot_list[[celltype]] <- p
}



corrplot_image <- ggdraw() + draw_image("corrplot.png", scale=1.03)
barplots_combined <- wrap_plots(barplot_list, ncol = length(barplot_list))

final_plot <- plot_grid(
  corrplot_image,
  barplots_combined,
  ncol = 2,  # Place plots side by side
  align = "hv",  # Align horizontally and vertically
  rel_heights = c(1.25, 1)  # Ensure both columns have equal width
)
print(final_plot)







genes_corplot <- sapply(strsplit(rownames(a$corr), "_"), '[', 2)

deseq <- read.csv("DESeq2_output_all_genes")
deseq <- as.data.frame(deseq)
deseq <- subset(deseq, select= -c(avg_log2FC, pct.1, pct.2, p_val_adj, gene))


deseq <- pivot_wider(deseq, names_from = cluster, values_from = p_val)
deseq <- deseq %>%
  mutate(across(where(is.numeric), ~ -log10(.)))

genes_with_des <- data.frame(Gene = genes_corplot) %>% 
  left_join(deseq, by = c("Gene" = "X"))

genes_with_des[is.na(genes_with_des)] <- 0
genes_with_des$Gene <- factor(genes_with_des$Gene, levels = genes_corplot)


barplot_list <- list()
# celltypes <- colnames(specificity)[-1]
celltypes <- list("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC")
celltypes <- list("sum_Astro", "sum_Tcell", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure", "sum_Neuron")


for (celltype in celltypes) {
  p <- ggplot(genes_with_des, aes_string(x="Gene", y=celltype)) +
    geom_col(fill = "lightpink") +
    ggtitle(paste(celltype)) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", angle = 90, size=10)) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("-log10(p_val)") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(genes_with_des$Gene))) +   # Reverse the order of the y-axis
    theme(axis.text.y = element_blank())
  
  barplot_list[[celltype]] <- p
}



corrplot_image <- ggdraw() + draw_image("corrplot.png", scale=1.03)
barplots_combined <- wrap_plots(barplot_list, ncol = length(barplot_list))

final_plot <- plot_grid(
  corrplot_image,
  barplots_combined,
  ncol = 2,  # Place plots side by side
  align = "hv",  # Align horizontally and vertically
  rel_heights = c(1.2, 1)
)
print(final_plot)





genes_corplot <- sapply(strsplit(rownames(a$corr), "_"), '[', 2)

deseq <- read.csv("output/celltypes/combined_all/spec_bayes.csv", row.names = 1)
deseq <- as.data.frame(-log10(as.matrix(deseq)))
deseq <- deseq %>% tibble::rownames_to_column("X")

genes_with_des <- data.frame(Gene = genes_corplot) %>% 
  left_join(deseq, by = c("Gene" = "X"))

genes_with_des[is.na(genes_with_des)] <- 0
genes_with_des$Gene <- factor(genes_with_des$Gene, levels = genes_corplot)


barplot_list <- list()
# celltypes <- colnames(specificity)[-1]
celltypes <- list("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC")
celltypes <- list("sum_Astro", "sum_Tcell", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure", "sum_Neuron")


for (celltype in celltypes) {
  p <- ggplot(genes_with_des, aes_string(x="Gene", y=celltype)) +
    geom_col(fill = "lightgreen") +
    ggtitle(paste(celltype)) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", angle = 90, size=10)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("probabilistic \n score") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(genes_with_des$Gene))) +   # Reverse the order of the y-axis
    theme(axis.text.y = element_blank())
  
  barplot_list[[celltype]] <- p
}



corrplot_image <- ggdraw() + draw_image("corrplot.png", scale=1.0)
barplots_combined <- wrap_plots(barplot_list, ncol = length(barplot_list))

final_plot <- plot_grid(
  corrplot_image,
  barplots_combined,
  ncol = 2,  # Place plots side by side
  align = "hv",  # Align horizontally and vertically
  rel_heights = c(1.1, 1)
)
print(final_plot)








genes_corplot <- sapply(strsplit(rownames(a$corr), "_"), '[', 2)

deseq <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv")


genes_with_des <- data.frame(Gene = genes_corplot) %>% 
  left_join(deseq, by = c("Gene" = "X"))

genes_with_des[is.na(genes_with_des)] <- 0
genes_with_des$Gene <- factor(genes_with_des$Gene, levels = genes_corplot)


barplot_list <- list()
# celltypes <- colnames(specificity)[-1]
celltypes <- list("sum_GDiv_tumor", "sum_Endo", "sum_Per", "sum_TAM", "sum_GTumor", "sum_Oligo", "sum_OPC")
celltypes <- list("sum_Astro", "sum_Tcell", "sum_ATumor", "sum_OTumor", "sum_ADiv_tumor", "sum_O_notsure", "sum_Neuron")
celltypes <- list("tau_score_vst")

for (celltype in celltypes) {
  p <- ggplot(genes_with_des, aes_string(x="Gene", y=celltype)) +
    geom_col(fill = "black") +
    ggtitle(paste(celltype)) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", angle = 90, size=10)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("tau score") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(genes_with_des$Gene))) +   # Reverse the order of the y-axis
    theme(axis.text.y = element_blank())
  
  barplot_list[[celltype]] <- p
}



corrplot_image <- ggdraw() + draw_image("corrplot.png", scale=1.0)
barplots_combined <- wrap_plots(barplot_list, ncol = length(barplot_list))

final_plot <- plot_grid(
  corrplot_image,
  barplots_combined,
  ncol = 2,  # Place plots side by side
  align = "hv",  # Align horizontally and vertically
  rel_heights = c(1.1, 1)
)
print(final_plot)


specificity <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv")
specificity <- as.data.frame(specificity)
specificity <- subset(specificity, select= -c(tau_score_vst, maximum, max_cluster, logmax))


genes_with_scores <- data.frame(Gene = genes_corplot) %>%
  left_join(specificity, by = c("Gene" = "X"))

genes_with_scores[is.na(genes_with_scores)] <- 0
genes_with_scores$Gene <- factor(genes_with_scores$Gene, levels = genes_corplot)
rownames(genes_with_scores) <- genes_with_scores$Gene
genes_with_scores$Gene <- NULL
genes_with_scores <- as.data.frame(lapply(genes_with_scores, as.numeric))
rownames(genes_with_scores) <- genes_corplot

# deseq <- read.csv("output/celltypes/combined_all/spec_bayes.csv", row.names = 1)
# deseq <- as.data.frame(-log10(as.matrix(deseq)))
# deseq <- deseq %>% tibble::rownames_to_column("X")
# 
# genes_with_scores <- data.frame(Gene = genes_corplot) %>% 
#   left_join(deseq, by = c("Gene" = "X"))
# 
# genes_with_scores[is.na(genes_with_scores)] <- 0
# genes_with_scores$Gene <- factor(genes_with_scores$Gene, levels = genes_corplot)
# rownames(genes_with_scores) <- genes_with_scores$Gene
# genes_with_scores$Gene <- NULL
# genes_with_scores <- as.data.frame(lapply(genes_with_scores, as.numeric))
# rownames(genes_with_scores) <- genes_corplot



left_block_indices <- 1:75          
middle_block_indices <- 76:143 
right_block_indices <- 144:200


left_block_genes <- genes_corplot[left_block_indices]
middle_block_genes <- genes_corplot[middle_block_indices]
second_middle <- genes_corplot[second_middle_indices]
right_block_genes <- genes_corplot[right_block_indices]

test_specificity_block <- function(scores, gene_names, block_genes) {
  # Subset block scores
  block_scores <- scores[gene_names %in% block_genes]
  other_scores <- scores[!gene_names %in% block_genes]
  
  # Remove NAs
  block_scores <- na.omit(block_scores)
  other_scores <- na.omit(other_scores)
  
  # Check if both groups have at least two values
  if (length(block_scores) < 2 || length(other_scores) < 2) {
    return(NA)  # Return NA if the test cannot be performed
  }
  
  # Perform Wilcoxon rank-sum test
  test_result <- wilcox.test(block_scores, other_scores, alternative = "greater")
  
  return(test_result$p.value)
}

# Iterate over each cell type and block
results <- lapply(colnames(genes_with_scores), function(cell_type) {
  scores <- genes_with_scores[[cell_type]]
  gene_names <- rownames(genes_with_scores)  # Assuming rownames are genes
  
  data.frame(
    Cell_Type = cell_type,
    Block = c("Left", "Middle", "Right"),
    P_Value = c(
      test_specificity_block(scores, gene_names, left_block_genes),
      test_specificity_block(scores, gene_names, middle_block_genes),
      test_specificity_block(scores, gene_names, right_block_genes)
    )
  )
})

# Combine results into a single data frame
p_value_results <- do.call(rbind, results)

# Adjust p-values for multiple testing
p_value_results$Adjusted_P_Value <- p.adjust(p_value_results$P_Value, method = "BH")

# Show significant results (adjusted p-value < 0.05), ignoring NA values
significant_results <- p_value_results[!is.na(p_value_results$Adjusted_P_Value) & p_value_results$Adjusted_P_Value < 0.05, ]
print(significant_results)