gene_expression <- data.frame(
  Gene = c("Gene2", "Gene5", "Gene10"), 
  CellType_A = c(5, 90, 800),
  CellType_B = c(200, 92, 50),
  CellType_C = c(10, 89, 20),
  CellType_D = c(220, 91, 10),
  CellType_E = c(15, 88, 100)
)

print(gene_expression)

rownames(gene_expression) <- gene_expression$Gene
gene_expression <- gene_expression[!colnames(gene_expression) %in% c("Gene") ]

tau <- apply(gene_expression, 1, calc_tau)
counts_tau <- data.frame(gene_expression, tau)

z <- t(apply(gene_expression, 1, z_score))
prop <- t(apply(gene_expression, 1, proportion))
spec <- round(prop * z, 3)
spec_tau <- data.frame(spec, tau)



gene_expression <- as.matrix(gene_expression)  # Flatten into a matrix
mode(gene_expression) <- "numeric"
vst_norm <- gene_expression / rowSums(gene_expression)
gene_means <- rowMeans(vst_norm)
gene_vars <- apply(vst_norm, 1, var)

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
specificity_scores <- -log10(specificity_scores)