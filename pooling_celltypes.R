
agg_celltype <- function(files, dirs) {
  for (i in 1:length(files)) {
    object <- readRDS(file = paste0("output/hijfte/", files[i]))
    agg_object <- AggregateExpression(object, group.by = "celltype", return.seurat = TRUE)
    Idents(agg_object) <- rownames(agg_object@assays$RNA@cells)
    
    
    counts <- agg_object@assays$RNA@layers$counts
    colnames(counts) <- rownames(agg_object@assays$RNA@cells@.Data)
    rownames(counts) <- Features(agg_object)
    # 
    # agg_object@assays$RNA@layers$counts <- ceiling(counts)
    # counts <- ceiling(counts)
    counts_vst <- counts %>%
      DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>%
      DESeq2::vst(blind=T)
    vst_counts <- as.data.frame(counts_vst@assays@data@listData)
    tau_score_vst <- apply(vst_counts, 1, calc_tau)
    counts_tau_vst <- data.frame(vst_counts, tau_score_vst)
    
    z_vst <- t(apply(vst_counts, 1, z_score))
    prop_vst <- t(apply(vst_counts, 1, proportion))
    spec <- round(prop_vst * z_vst, 3)
    spec_tau <- data.frame(spec, tau_score_vst)
    
    spec_tau$maximum <- apply(spec, 1, max)
    spec_tau$max_cluster <- colnames(spec)[apply(spec, 1, which.max)]
    spec_tau$logmax <- log(spec_tau$maximum)
    
    specific_genes <- subset(x = spec_tau, subset = logmax >= 0 & tau_score_vst >= 0.6)
    specific_genes_celltype <- data.frame(celltype = specific_genes$max_cluster) 
    rownames(specific_genes_celltype) <- rownames(specific_genes)
    
    markers <- FindAllMarkers(agg_object, test.use = "DESeq2", min.cells.group = 0, min.cells.feature = 0,
                              #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                              return.thresh = 1.1, logfc.threshold = 0.25, min.pct = 0.1)
    top10 <- markers %>%
      group_by(cluster) %>%
      top_n(n = -10, wt = p_val) %>%
      summarise(genes = paste(gene, collapse = ", "))

    top10_genes <- unlist(strsplit(top10$genes, ", "))
    # common_genes <- intersect(specific_genes, top10_genes)

    vst_counts <- as.matrix(vst_counts)  # Flatten into a matrix
    mode(vst_counts) <- "numeric"
    vst_norm <- vst_counts / rowSums(vst_counts)
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
    
    
    plotlog <- ggplot(spec_tau, aes(x=tau_score_vst, y=logmax)) +
      # Color points based on whether the gene is in top10 or not
      geom_point()+
      geom_segment(y=0, x = 0.6, xend = 0.9, linetype = "dashed", color = "black") +  
      geom_segment(x=0.6, y=0, yend=2, linetype = "dashed", color = "black") +
      
     
      # Add labels and theme
      labs(
        title = "Relationship between Tau Score and Log Maximum Specificity",
        x = "Tau Score",
        y = "Log Maximum Expression",
        color = "Gene Category",
        caption = "Dashed lines represent thresholds used for analysis, genes within this corner are said to \n be cell type-specific. The value Log Maximum Specificity is the log transformed highest \n specificity score for a gene, using the z*prop metric. Tau represents the spread of \n expression. The red dots represent the genes that were DE according to DESeq analysis on \n the raw, per cluster aggregated counts."
      ) +
      theme_minimal()
    
    
    plotclus <- ggplot(spec_tau, aes(x=tau_score_vst, y=logmax)) +
      # Color points based on whether the gene is in top10 or not
      geom_point(aes(color = max_cluster)) +
      geom_segment(y=0, x = 0.6, xend = 0.9, linetype = "dashed", color = "black") +  
      geom_segment(x=0.6, y=0, yend=2, linetype = "dashed", color = "black") +
      
      # Define custom colors for top10 and other genes
      #scale_color_viridis_d(name = "Max Cluster") +
      
      # Add labels and theme
      labs(
        title = "Relationship between Tau Score and Log Maximum Specificity",
        x = "Tau Score",
        y = "Log Maximum Expression",
        color = "Gene Category",
        caption = "Colored based on which cluster has highest specificity"
      ) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set3")
    
    
    plot <- ggplot(spec_tau, aes(x=tau_score_vst, y=maximum)) +
      geom_point() + # Show dots
      geom_text(
        label=rownames(spec_tau), 
        nudge_x = 0, nudge_y = 0, 
        check_overlap = T
      )
    
    # markers$diffexpressed <- "NO"
    # markers$diffexpressed[markers$avg_log2FC > 4 & markers$p_val < 0.01] <- "UP"
    # markers$diffexpressed[markers$avg_log2FC < -2.5 & markers$p_val < 0.05] <- "DOWN"
    # 
    # markers$gene <- rownames(markers)
    # 
    # volcano <- ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val), color = diffexpressed)) +
    #   geom_point() +
    #   scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
    #                      labels = c("Downregulated", "Not significant", "Upregulated")) +
    #   ggtitle("Differentially Expressed Genes") +
    #   geom_text_repel(data = subset(markers, diffexpressed != "NO"), # Label only "UP" and "DOWN" genes
    #                   aes(label = gene), max.overlaps = Inf, box.padding = 0.3) +
    #   theme_minimal() +
    #   facet_wrap(~ cluster)
    
    if (!dir.exists(dirs[i])) {
      dir.create(dirs[i], recursive = TRUE)
    }
    
    write.csv(spec_tau, file = file.path(dirs[i], "spec_prop.zscore_tau.csv"))
    write.csv(specificity_scores, file = file.path(dirs[i], "spec_bayes.csv"))
    png(file = file.path(dirs[i], "plot_logmaxvstau.png"), width = 960, height = 960)
    print(plotlog)  # Ensure the plot is printed inside the png device
    dev.off()
    png(file = file.path(dirs[i], "plot_clus.png"), width = 960, height = 960)
    print(plotclus)  # Ensure the plot is printed inside the png device
    dev.off()
    png(file = file.path(dirs[i], "plot_maxvstau.png"), width = 960, height = 960)
    print(plot)  # Ensure the plot is printed inside the png device
    dev.off()
    # png(file = file.path(dirs[i], "volcano.png"), width = 960, height = 960)
    # print(volcano)  # Ensure the plot is printed inside the png device
    # dev.off()
    
    formatted_gene_list <- paste(top10$cluster, top10$genes, sep = "\t")
    writeLines(formatted_gene_list, con = file.path(dirs[i], "top10_genes_per_cluster.txt"))
    write.csv(specific_genes_celltype, file = file.path(dirs[i], "spec_gene_celltype.csv"))
    # write.csv(markers, file = file.path(dirs[i], "DESeq_output.csv"))
    
    cat("Files saved successfully in:", dirs[i], "\n")
    
  }
}  

files <- list.files(path= "output/diaz/", pattern = "*.rds")
out_dirs <- list.dirs(path = "output/celltypes/diaz")[- 1] 
agg_celltype(files, out_dirs)

files <- list.files(path= "output/yuan/", pattern = "*.rds")
out_dirs <- list.dirs(path = "output/celltypes/yuan")[- 1] 
agg_celltype(files, out_dirs)

files <- list.files(path= "output/couturier/", pattern = "*.rds")
out_dirs <- list.dirs(path = "output/celltypes/couturier")[- 1] 
agg_celltype(files, out_dirs)

files <- list.files(path= "output/bolleboom/", pattern = "*.rds")
out_dirs <- list.dirs(path = "output/celltypes/bolleboom")[- 1] 
agg_celltype(files, out_dirs)

files <- list.files(path= "output/hijfte/", pattern = "*r.rds")
out_dirs <- list.dirs(path = "output/celltypes/hijfte")[- 1] [1]

agg_celltype(files, out_dirs)

files <- list.files(path= "output/diaz_astrooligo/", pattern = "*.rds")
out_dirs <- list.dirs(path = "output/celltypes/diaz_astrooligo")[- 1] 
agg_celltype(files, out_dirs)





z_score <- function(expr){
  sapply(seq_along(expr), function(i) {
    sd_val <- sd(expr[-i])
    if (sd_val == 0) { #nog naar kijken, 
      return((expr[i] - mean(expr[-i])) / 1)  
    } else {
      return((expr[i] - mean(expr[-i])) / sd_val)
    }
  })
}

proportion <- function(expr) {
  sapply(seq_along(expr), function(i){
    total <- sum(expr)
    return(expr[i] / total)
  })
}

calc_tau <- function(x) {
  max_x <- max(x)
  tau <- sum(1-(x/max_x)) / (length(x)-1)
  return(tau)
}

seuratobj <- CreateSeuratObject(agg_expr)
seuratobj@assays$RNA$counts
seuratobj$orig.ident <- colnames(agg_expr$RNA)
head(seuratobj@meta.data$orig.ident)
head(seuratobj)
Idents(seuratobj) <- "orig.ident"
head(seuratobj)
markers <- FindAllMarkers(seuratobj, test.use = "DESeq2", min.cells.group = 0, min.cells.feature = 0,
                          #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                          return.thresh = 1.1, logfc.threshold = 0.25, min.pct = 0.1)
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = p_val) %>%
  summarise(genes = paste(gene, collapse = ", "))

top10_genes <- unlist(strsplit(top10$genes, ", "))
common_genes <- intersect(specific_genes, top10_genes)
