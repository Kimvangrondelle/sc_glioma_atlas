
files5 <- list.files(path= "output/diaz/", pattern = "*.rds", full.names = TRUE)[-6]
files1 <- list.files(path= "output/yuan/", pattern = "*.rds", full.names = TRUE)
files2 <- list.files(path= "output/couturier/", pattern = "*.rds", full.names = TRUE)
files3 <- list.files(path= "output/bolleboom/", pattern = "*.rds", full.names = TRUE)[1]
files4 <- list.files(path= "output/hijfte/", pattern = "*.rds", full.names = TRUE)[-2]
files6 <- list.files(path= "output/diaz_astrooligo/", pattern = "*.rds", full.names = TRUE)

rds_files <- c(files5, files1, files2, files3, files4, files6)


objects <- list()
files_all <- c(files5, files1, files2, files3, files4, files6)
for (file in files_all) {
  object <- readRDS(file)
  objects[[length(objects) + 1]] <- object}


generate_pseudobulk <- function(seurat_list, num_cells, num_samples, condition_celltype, condition_fraction, excluded_celltypes) {
  
  # Combineer metadata van alle Seurat-objecten
  all_cells <- do.call(rbind, lapply(seurat_list, function(obj) {
    data.frame(Cell = colnames(obj), CellType = obj$celltype, Sample = obj@project.name)
  }))
  
  if (!is.null(excluded_celltypes)) {
    all_cells <- all_cells %>% filter(!CellType %in% excluded_celltypes)
  }
  
  set.seed(42)  # Reproduceerbaarheid
  
  # Hulpfunctie om een pseudo-bulk sample te genereren
  create_sample <- function(selected_cells, seurat_list) {
    selected_matrices <- lapply(seurat_list, function(obj) {
      subset_cells <- intersect(colnames(obj), selected_cells$Cell)
      
      if (length(subset_cells) > 0) {
        counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
        return(Matrix::rowSums(counts_matrix[, subset_cells, drop = FALSE]))
      } else {
        return(NULL)
      }
    })
    
    # Verwijder NULL waarden en combineer expressie
    selected_matrices <- Filter(Negate(is.null), selected_matrices)
    if (length(selected_matrices) == 0) {
      stop("Geen overlappende cellen gevonden in de Seurat-objecten.")
    }
    
    return(Reduce("+", selected_matrices) / length(selected_matrices))  # Gemiddelde expressie per gen
  }
  
  # Maak 10 pseudo-bulk samples voor de controle groep (uniforme verdeling)
  control_samples <- replicate(num_samples, {
    sampled_cells <- all_cells %>%
      group_by(CellType) %>%
      sample_n(size = round(num_cells / n_distinct(all_cells$CellType)), replace = TRUE)
    create_sample(sampled_cells, seurat_list)
  }, simplify = FALSE)
  
  # Maak 10 pseudo-bulk samples voor de conditie groep (50% specifiek celtype)
  condition_samples <- replicate(num_samples, {
    target_count <- round(num_cells * condition_fraction)
    remaining_count <- num_cells - target_count
    
    target_cells <- all_cells %>%
      filter(CellType == condition_celltype) %>%
      sample_n(target_count, replace = TRUE)
    
    remaining_cells <- all_cells %>%
      filter(CellType != condition_celltype) %>%
      sample_n(remaining_count, replace = TRUE)
    
    condition_cells <- bind_rows(target_cells, remaining_cells)
    create_sample(condition_cells, seurat_list)
  }, simplify = FALSE)
  
  # Zet de samples om in een matrix voor DESeq2
  all_samples <- c(control_samples, condition_samples)
  pseudobulk_matrix <- do.call(cbind, all_samples)
  colnames(pseudobulk_matrix) <- c(paste0("control_", 1:num_samples), paste0("condition_", 1:num_samples))
  
  pseudobulk_matrix <- round(pseudobulk_matrix)
  pseudobulk_matrix[pseudobulk_matrix < 0] <- 0 
  
  # Metadata voor DESeq2
  colData <- data.frame(
    sample = colnames(pseudobulk_matrix),
    condition = rep(c("control", "condition"), each = num_samples),
    row.names = colnames(pseudobulk_matrix)
  )
  
  return(list(counts = pseudobulk_matrix, colData = colData))
}


result <- generate_pseudobulk(objects, num_cells = 1000, num_samples = 10, condition_celltype = "Astro", condition_fraction = 0.5, excluded_celltypes = NULL)
counts_matrix <- result$counts
colData <- result$colData

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

res_sorted <- res[order(res$pvalue), ]
filtered_genes <- rownames(res_sorted)[!grepl("^(RP11|RP3|RP13|RP4|RP5|RP1-)", rownames(res_sorted))]
res_sorted <- res_sorted[filtered_genes, ]
res_significant <- res_sorted[res_sorted$pvalue < 0.01, ]
pseudo_sig <- rownames(res_significant)

sig_des <- top_100_genes_per_cluster[["sum_Astro"]]
overlap <- intersect(pseudo_sig, sig_des)



t_scores_df <- data.frame(Gene = rownames(res), t_score = res@listData$stat)

df <- read.csv("output/celltypes/combined_all/spec_prop.zscore_tau.csv", row.names = 1)
subset_df <- data.frame(Gene = rownames(df), Sum_Astro = df$sum_Astro, row.names = NULL)

merged_df <- inner_join(t_scores_df, subset_df, by = "Gene")
rownames(merged_df) <- merged_df$Gene

overlap_genes <- merged_df[rownames(merged_df) %in% overlap, ]

ggplot(merged_df, aes(x = t_score, y = Sum_Astro)) +
  geom_point() +
  # geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Scatterplot of own score vs t-score deseq Astro",
       x = "t_score",
       y = "own score for astro") +
  theme_minimal() 
cor(merged_df$t_score, merged_df$Sum_Astro)

# intersect(pseudo_sig, score_o)
# intersect(pseudo_sig, des_o)

# 
# ranks1 <- match(overlap, pseudo_sig)
# ranks2 <- match(overlap, sig_des)
# 
# spearman_corr <- cor(ranks2, ranks1, method = "spearman")
# kendall_corr <- cor(ranks2, ranks1, method = "kendall")
# 
# cat("Spearman's rank correlation:", spearman_corr, "\n")
# cat("Kendall's tau correlation:", kendall_corr, "\n")

hijfte <- readRDS("output/hijfte/hijfte-y.rds")
cell_counts_by_type <- table(hijfte@meta.data$celltype)
