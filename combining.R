files5 <- list.files(path= "output/diaz/", pattern = "*.rds", full.names = TRUE)
files1 <- list.files(path= "output/yuan/", pattern = "*.rds", full.names = TRUE)
files2 <- list.files(path= "output/couturier/", pattern = "*.rds", full.names = TRUE)
files3 <- list.files(path= "output/bolleboom/", pattern = "*.rds", full.names = TRUE)
files4 <- list.files(path= "output/hijfte/", pattern = "*.rds", full.names = TRUE)
files6 <- list.files(path= "output/diaz_astrooligo/", pattern = "*.rds", full.names = TRUE)


dfs <- list()
files_all <- c(files5, files1, files2, files3, files4, files6)
for (file in files_all) {
  object <- readRDS(file)
  agg_object <- AggregateExpression(object, group.by = "celltype", return.seurat = TRUE)
  counts <- agg_object@assays$RNA@layers$counts
  rownames(counts) <- Features(agg_object)
  colnames(counts) <- rownames(agg_object@assays$RNA@cells@.Data) 
  counts <- as.data.frame(counts)
  counts <- rownames_to_column(counts, var = "gene")
  dfs[[length(dfs) + 1]] <- counts
}

merged_df <- Reduce(function(x, y) {
  
  merged <- merge(x, y, by = "gene", all = TRUE)
  
  
  return(merged)
}, dfs)

rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL
merged_df[is.na(merged_df)] <- 0

# colnames(merged_df) <- gsub("\\..*", "", colnames(merged_df))
colnames(merged_df) <- make.names(colnames(merged_df), unique = TRUE)
# Sum counts per row by celltype prefix
df_sum <- merged_df %>%
  rowwise() %>%
  mutate(
    sum_GDiv_tumor = sum(c_across(starts_with("GDiv")), na.rm = TRUE),
    sum_Endo = sum(c_across(starts_with("Endo")), na.rm = TRUE),
    sum_Per = sum(c_across(starts_with("Per")), na.rm = TRUE),
    sum_TAM = sum(c_across(starts_with("TAM")), na.rm = TRUE),
    sum_Neuron = sum(c_across(starts_with("Neuron")), na.rm = TRUE),
    sum_Oligo = sum(c_across(starts_with("Oligo")), na.rm = TRUE),
    sum_OPC = sum(c_across(starts_with("OPC")), na.rm = TRUE),
    sum_Tcell = sum(c_across(starts_with("Tcell")), na.rm = TRUE),
    sum_Astro = sum(c_across(starts_with("Astro")), na.rm = TRUE), 
    sum_GTumor = sum(c_across(starts_with("GTumor")), na.rm = TRUE), 
    sum_ATumor = sum(c_across(starts_with("ATumor")), na.rm = TRUE), 
    sum_OTumor = sum(c_across(starts_with("OTumor")), na.rm = TRUE), 
    sum_ADiv_tumor = sum(c_across(starts_with("ADiv")), na.rm = TRUE), 
    sum_O_notsure = sum(c_across(starts_with("O?")), na.rm = TRUE)
  ) %>%
  select(sum_GDiv_tumor, sum_Endo, sum_Per, sum_TAM, sum_GTumor, sum_Oligo, sum_OPC, sum_Tcell, sum_Astro, sum_Neuron, sum_ATumor, sum_OTumor, sum_ADiv_tumor, sum_O_notsure)

rownames(df_sum) <- rownames(merged_df)

write.csv(df_sum, "counts_combined_summed.csv")

