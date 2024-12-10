vst_counts <- read.csv("output/celltypes/combined/vst_counts.csv", row.names = 1)
vst_counts <- as.matrix(vst_counts)  # Flatten into a matrix
mode(vst_counts) <- "numeric"

vst_norm <- vst_counts / rowSums(vst_counts)


gene_means <- rowMeans(vst_norm)
gene_vars <- apply(vst_norm, 1, var)


posterior_probs <- function(celltype_expression, mean, var) {
  likelihood <- dnorm(celltype_expression, mean = mean, sd = sqrt(var))
  prior <- mean / sum(mean)
  posterior <- likelihood * prior
  return(posterior)
}

specificity_scores <- matrix(0, nrow = nrow(vst_norm), ncol= ncol(vst_norm))
rownames(specificity_scores) <- rownames(vst_norm)
colnames(specificity_scores) <- colnames(vst_norm)

for (gene_idx in 1:nrow(vst_norm)) {
  gene_expression <- vst_norm[gene_idx, ]
  gene_mean <- gene_means[gene_idx]
  gene_var <- gene_vars[gene_idx]
  
  specificity_scores[gene_idx, ] <- sapply(1:ncol(vst_norm), function(celltype_idx)
    posterior_probs(celltype_expression = gene_expression[celltype_idx], 
                    mean = gene_mean, 
                    var = gene_var))
}

specificity_scores <- specificity_scores / rowSums(specificity_scores)