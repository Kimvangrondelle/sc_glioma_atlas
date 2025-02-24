save_objects <- function(samples, seurat_objects) {
  for (a in 1:length(samples)){
    saveRDS(objects[[a]], file = paste0("output/bolleboom/", samples[a], ".rds"))
  }
}

samples <- list("H243")
# samples <- list("yuan25", "yuan30") etc --- f.e. cou...
save <- save_objects(samples, objects)



markers <- function(list_of_files) {
  markers_per_clus_per_sam <- list()
  for (i in 1:length(list_of_files)) {
    object <- readRDS(file = paste0("../output/diaz/", files[i]))
    
    all_markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    print(colnames(all_markers))
    top20_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) %>%
      summarise(genes = paste(gene, collapse = ", "))
    markers_per_clus_per_sam[[files[i]]] <- top20_markers
  }
  return(markers_per_clus_per_sam)
}

files <- list.files(path= "../output/diaz/", pattern = "*.rds")
mar_per_clus <- markers(files)


# mar_per_clus[[--filename--]][["genes"]][--clusternummer--] om de genen op te roepen
# mar_per_clus[["yuan16.rds"]][["genes"]][1] voor de top20 genen van cluster 0 in yuan16
print(mar_per_clus[["hijfte-y.rds"]][["genes"]])