#not used in final thesis
# compare specificity scores of cell type in combined sample to same sample in the individual samples

list_celltypes <- list("sum_GTumor", "sum_ATumor", "sum_OTumor", "sum_GDiv tumor", "sum_ADiv tumor", "sum_Oligo", "sum_TAM",
                       "sum_Per", "sum_Endo", "sum_Astro", "sum_Neuron", "sum_OPC", "sum_Tcell", "sum_O_notsure")

# list_celltypes <- list("sum_Tumor", "sum_Div_tumor", "sum_Oligo", "sum_TAM", "sum_Per", 
#                        "sum_Endo", "sum_Astro", "sum_Neuron", "sum_OPC", "sum_Tcell")

#prepare combined sample
pseudobulk <- read.csv("counts_combined_summed_all.csv", row.names = 1)
pseudobulk <- as.data.frame(pseudobulk)

pseudobulk <- sweep(pseudobulk, 1, rowSums(pseudobulk), FUN = "/")

#prepare individual sample
sample <- "hijfte/pr3"  # for example: bolleboom/H243
file <- paste0("output/celltypes/", sample ,"/spec_prop.zscore_tau.csv")

specificity_scores_sample <- read.csv(file, row.names = 1)
specificity_scores_sample <- as.data.frame(specificity_scores_sample)
specificity_scores_sample <- subset(specificity_scores_sample, select= -c(tau_score_vst, maximum, max_cluster, logmax))

rows_with_zeros <- specificity_scores_sample %>%
  filter(if_all(everything(), ~ . == 0))

rows_remove <- rownames(rows_with_zeros)
celltypes_in_sample <- list(colnames(specificity_scores_sample))
specificity_scores_sample <- specificity_scores_sample[!rownames(specificity_scores_sample) %in% rows_remove, ]

#only use genes present in both combined as individual sample
common <- intersect(rownames(pseudobulk), rownames(specificity_scores_sample))

pseudo_filtered <- pseudobulk[common, , drop= FALSE]
specificity_sample_filtered <- specificity_scores_sample[common, , drop = FALSE]


#combined sample filtered on cell type
bulk_astro <- pseudo_filtered[, "sum_Astro", drop=FALSE]
bulk_oligo <- pseudo_filtered[, "sum_Oligo", drop=FALSE]
bulk_tam <- pseudo_filtered[, "sum_TAM", drop=FALSE]
bulk_gtumor <- pseudo_filtered[, "sum_GTumor", drop=FALSE]
bulk_gdiv <- pseudo_filtered[, "sum_GDiv_tumor", drop=FALSE]
bulk_neuron <- pseudo_filtered[, "sum_Neuron", drop=FALSE]
bulk_opc <- pseudo_filtered[, "sum_OPC", drop=FALSE]
bulk_endo <- pseudo_filtered[, "sum_Endo", drop=FALSE]
bulk_per <- pseudo_filtered[, "sum_Per", drop=FALSE]
bulk_tcell <- pseudo_filtered[, "sum_Tcell", drop=FALSE]
bulk_atumor <- pseudo_filtered[, "sum_ATumor", drop=FALSE]
bulk_otumor <- pseudo_filtered[, "sum_OTumor", drop=FALSE]
bulk_adiv <- pseudo_filtered[, "sum_ADiv_tumor", drop=FALSE]
bulk_onot <- pseudo_filtered[, "sum_O_notsure", drop=FALSE]

#individual sample filtered on cell type
sample_astro <- specificity_sample_filtered[, "Astro", drop=FALSE]
sample_oligo <- specificity_sample_filtered[, "Oligo", drop=FALSE]
sample_tam <- specificity_sample_filtered[, "TAM", drop=FALSE]
sample_neuron <- specificity_sample_filtered[, "Neuron", drop=FALSE]
sample_gtumor <- specificity_sample_filtered[, "GTumor", drop=FALSE]
sample_gdiv <- specificity_sample_filtered[, "GDiv.tumor", drop=FALSE]
sample_endo <- specificity_sample_filtered[, "Endo", drop=FALSE]
sample_per <- specificity_sample_filtered[, "Per", drop=FALSE]
sample_opc <- specificity_sample_filtered[, "OPC", drop=FALSE]
sample_tcell <- specificity_sample_filtered[, "Tcell", drop=FALSE]
sample_atumor <- specificity_sample_filtered[, "ATumor", drop=FALSE]
sample_otumor <- specificity_sample_filtered[, "OTumor", drop=FALSE]
sample_adiv <- specificity_sample_filtered[, "ADiv.tumor", drop=FALSE]
sample_onot <- specificity_sample_filtered[, "O.", drop=FALSE]

#spearman correlation of scores for cell type between combined and individual
cor(sample_astro, bulk_astro, method = "spearman")
cor(sample_oligo, bulk_oligo, method = "spearman")
cor(sample_tam, bulk_tam, method = "spearman")
cor(sample_neuron, bulk_neuron, method = "spearman")
cor(sample_gtumor, bulk_gtumor, method = "spearman")
cor(sample_gdiv, bulk_gdiv, method = "spearman")
cor(sample_endo, bulk_endo, method = "spearman")
cor(sample_per, bulk_per, method = "spearman")
cor(sample_opc, bulk_opc, method = "spearman")
cor(sample_tcell, bulk_tcell, method = "spearman")
cor(sample_atumor, bulk_atumor, method = "spearman")
cor(sample_otumor, bulk_otumor, method = "spearman")
cor(sample_adiv, bulk_adiv, method = "spearman")
cor(sample_onot, bulk_onot, method = "spearman")

# df_astro <- cbind(bulk_astro, sample_astro)
# df_oligo <- cbind(bulk_oligo, sample_oligo)
# df_tam <- cbind(bulk_tam, sample_tam)
# df_neuron <- cbind(bulk_neuron, sample_neuron)
# df_gtumor <- cbind(bulk_tumor, sample_tumor)
# df_gdiv <- cbind(bulk_div, sample_div)
# df_endo <- cbind(bulk_endo, sample_endo)
# df_per <- cbind(bulk_per, sample_per)
# df_opc <- cbind(bulk_opc, sample_opc)

#plot scores
ggplot(df_astro, aes(x = sum_Astro, y = Astro)) +
  geom_point() +
  labs(x= "bulk counts", y= "specificity score", title = "hijfte Astro") +
  xlim(-1, 1)+
  theme_minimal()