gene_expression <- data.frame(
  Gene = paste0("Gene", 1:10), 
  CellType_A = c(100, 5, 30, 500, 90, 10, 70, 15, 55, 800),
  CellType_B = c(105, 200, 35, 20, 92, 400, 75, 300, 58, 50),
  CellType_C = c(102, 10, 32, 10, 89, 5, 72, 10, 57, 20),
  CellType_D = c(98, 220, 28, 5, 91, 450, 68, 330, 52, 10),
  CellType_E = c(103, 15, 33, 400, 88, 8, 74, 12, 56, 700)
)

print(gene_expression)

# whether the expression is widespread or specific to a cell type
#moest sws erbij
calc_tau <- function(x) {
  max_x <- max(x)
  tau <- sum(1-(x/max_x)) / (length(x)-1)
  return(tau)
}


# genes <- rownames(pseudo_diaz10022@assays$RNA@features@.Data)
# pseudo_diaz10022_counts <- pseudo_diaz10022@assays$RNA@layers$counts
# 
# rownames(pseudo_diaz10022_counts) <- genes

pseudo_diaz1002_with_tau <- data.frame(diaz_10022_counts, tau_score_diaz)
pseudo_diaz1002_with_tau #summed counts + tau score

# calc_tsi <- function(x){ # was moeilijk met threshold zetten 
#   maxx <- max(x)
#   tsi <- maxx / sum(x)
#   return(tsi)
# }
# tsi_score <- apply(gene_expression[, -1], 1, calc_tsi)#lapply
# tsi_score
# cv_score <- apply(gene_expression[, -1], 1, function(x) sd(x) / mean(x))

# gen_ex <- data.frame(gene_expression, tau_score)
# gen_ex <- data.frame(gen_ex, tsi_score)
# gen_ex <- data.frame(gen_ex, cv_score)
# print(gen_ex)

# z_norm -- positive when the expression is higher than the mean
# -- negative when the expression is lower than the mean
#-- magnitude tells how much higher or lower the expression is compared across the different cell types
# z_norm <- t(apply(gene_expression[, -1], 1, function(x) scale(x)))
# z_norm <- data.frame(Gene = gene_expression$Gene, z_norm)
# print(z_norm)

# z_score_new <- function(expression){
#   sapply(seq_along(expression), function(i) {
#     other_items <- expression[-i]
#     mean_other <- mean(other_items)
#     sd_other <- sd(other_items)
#     (expression[i] - mean_other) / sd_other
#   })
# }
normalized <- t(apply(gene_expression[, -1], 1, z_score))
prop <- t(apply(gene_expression[,-1], 1, function(x) x/sum(x)))
spec <- round(prop * normalized, 3)
log(spec)

# gene_expression
# normalized

# # proportion of expression per cell type
# prop <- t(apply(gene_expression[, -1], 1, function(x) x / sum(x)))
# prop <- data.frame(Gene = gene_expression$Gene, prop)
# print(prop)

#variance coefficient -- how much variation between the expression in the cell types
# proportion of expression
# multiplied proportion with the variance coefficient
# Higher values when high variance 
# low variance is less interesting. 

tau_score <- apply(gene_expression[, -1], 1, calc_tau)
tau_score

cv_score <- apply(gene_expression[, -1], 1, function(x) sd(x) / mean(x))
prop <- t(apply(gene_expression[, -1], 1, function(x) x / sum(x)))

prop_cv <- prop * cv_score
prop_tau <- prop * tau_score
print(prop_cv)
print(prop_tau) <- data.frame(prop_cv, tau_score)
rownames(result) <- rownames(gene_expression)
result


pem_log <- function(expression, epsilon = 1e-6) {
  pem_scores <- log2((expression / sum(expression)) + epsilon)
  return(pem_scores)
}
specificity_pem <- t(apply(gene_expression[, -1], 1, pem_log))
specificity_pem

specificity_with_cv <- function(expression) {
  cv <- sd(expression) / mean(expression)
  specificity_scores <- (expression / sum(expression)) * cv
  return(specificity_scores)
}

specificity_score_cv <- t(apply(gene_expression[, -1], 1, specificity_with_cv))
specificity_score_cv <- data.frame(Gene = gene_expression$Gene, specificity_score_cv)
print(specificity_score_cv)

