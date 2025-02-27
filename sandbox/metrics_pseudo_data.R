#pseudo data and application of metrics on this data for sups in thesis. 

#pseudo data of three genes, one widespread expression, one specific expression and one with 2 high expression values. 
gene_expression <- data.frame(
  Gene = c("Gene2", "Gene5", "Gene10"), 
  CellType_A = c(5, 90, 800),
  CellType_B = c(200, 92, 50),
  CellType_C = c(10, 89, 20),
  CellType_D = c(220, 91, 10),
  CellType_E = c(15, 88, 100)
)

print(gene_expression)
#make first column the rownames
rownames(gene_expression) <- gene_expression$Gene
gene_expression <- gene_expression[!colnames(gene_expression) %in% c("Gene") ]
#apply tau on expression
tau <- apply(gene_expression, 1, calc_tau)
counts_tau <- data.frame(gene_expression, tau)
#apply zpex score on expression
z <- t(apply(gene_expression, 1, z_score))
prop <- t(apply(gene_expression, 1, proportion))
spec <- round(prop * z, 3)
spec_tau <- data.frame(spec, tau)


#apply bayes on expression -- make numeric, normalize 
gene_expression <- as.matrix(gene_expression)  # Flatten into a matrix
mode(gene_expression) <- "numeric"
vst_norm <- gene_expression / rowSums(gene_expression)
#calculate mean and variance per gene
gene_means <- rowMeans(vst_norm)
gene_vars <- apply(vst_norm, 1, var)

specificity_scores <- matrix(0, nrow = nrow(vst_norm), ncol= ncol(vst_norm))
rownames(specificity_scores) <- rownames(vst_norm)
colnames(specificity_scores) <- colnames(vst_norm)

for (gene_idx in 1:nrow(vst_norm)) {
  gene_expression <- vst_norm[gene_idx, ]
  gene_mean <- gene_means[gene_idx]
  gene_var <- gene_vars[gene_idx]
  #per row, a score is given for each cell type
  specificity_scores[gene_idx, ] <- sapply(1:ncol(vst_norm), function(celltype_idx)
    posterior_probs(celltype_expression = gene_expression[celltype_idx], 
                    mean = gene_mean, 
                    var = gene_var))
}
#scores are normalized to sum to 1
specificity_scores <- specificity_scores / rowSums(specificity_scores)
#scores are log10 transformed -- as the genes were most specific got the highest scores 
# -- this was logical as specific means no evidence for uniformity what was expected. 
specificity_scores <- -log10(specificity_scores)


#test tau on aggregated and transformed data
gene_expression <- data.frame(
  Gene = paste0("Gene", 1:10), 
  Astro = c(0, 0, 32, 6, 0, 37, 0, 13, 0, 1),
  Endo = c(0, 0, 37, 8, 5, 35, 0, 39, 0, 1),
  Gdiv = c(1, 11, 510, 34, 30, 391, 3, 341, 2, 1),
  Gtumor = c(2, 28, 2767, 400, 170, 2232, 13, 1690, 9, 5),
  Neuron = c(0, 1, 21, 12, 4, 30, 0, 22, 0, 0),
  Oligo = c(0, 18, 713, 121, 59, 704, 2, 376, 3, 4),
  Per = c(0, 0, 48, 6, 7, 66, 0, 25, 0, 0),
  TAM = c(0, 17, 648, 67, 32, 479, 3, 348, 1, 2)
)


gene_expression <- data.frame(
  Gene = paste0("Gene", 1:10), 
  Astro = c(1.08, 1.08, 7.48, 5.20, 1.08, 7.69, 1.08, 6.23, 1.08, 3.20),
  Endo = c(1.08, 1.08, 7.31, 5.23, 4.64, 7.23, 1.08, 7.39, 1.08, 2.97),
  Gdiv = c(1.67, 2.92, 7.55, 4.01, 3.88, 7.18, 2.09, 6.99, 1.91, 1.67),
  Gtumor = c(1.44, 2.39, 7.56, 4.94, 3.93, 7.26, 1.99, 6.87, 1.84, 1.65),
  Neuron = c(1.08, 3.42, 7.23, 6.45, 4.99, 7.73, 1.08, 7.29, 1.08, 1.08),
  Oligo = c(1.08, 3.02, 7.49, 5.08, 4.21, 7.47, 1.77, 6.59, 1.92, 2.04),
  Per = c(1.08, 1.08, 7.45, 4.66, 4.86, 7.90, 1.08, 6.54, 1.08, 1.08),
  TAM = c(1.08, 3.22, 7.77, 4.71, 3.85, 7.34, 2.05, 6.89, 1.64, 1.87)
)