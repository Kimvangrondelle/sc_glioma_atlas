#compare genes from McKenzie with genes highest scores per cell type

# list_celltypes <- list("sum_GTumor", "sum_ATumor", "sum_OTumor", "sum_GDiv tumor", "sum_ADiv tumor", "sum_Oligo", "sum_TAM", 
#                        "sum_Per", "sum_Endo", "sum_Astro", "sum_Neuron", "sum_OPC", "sum_Tcell", "sum_O_notsure")

list_celltypes <- list("sum_Tumor", "sum_Div_tumor", "sum_Oligo", "sum_TAM", "sum_Per", 
                       "sum_Endo", "sum_Astro", "sum_Neuron", "sum_OPC", "sum_Tcell")


bulk <- read.csv("counts_combined_summed_all.csv", row.names = 1)
bulk <- as.data.frame(bulk)

bulk <- sweep(bulk, 1, rowSums(bulk), FUN = "/")


sample <- "combined_all"  # for example: bolleboom/H243
file <- paste0("output/celltypes/", sample ,"/spec_prop.zscore_tau.csv")
file <- paste0("output/celltypes/combined/spec_tau.csv")

specificity_scores_sample <- read.csv(file, row.names = 1)
specificity_scores_sample <- as.data.frame(specificity_scores_sample)
specificity_scores_sample <- subset(specificity_scores_sample, select= -c(tau_score_vst, maximum, max_cluster, logmax))

#normalized counts combined sample
bulk_astro <- bulk[, "sum_Astro", drop=FALSE] 
bulk_astro <- bulk_astro[order(bulk_astro$sum_Astro, decreasing = TRUE), , drop = FALSE]
bulk_oligo <- bulk[, "sum_Oligo", drop=FALSE]
bulk_oligo <- bulk_oligo[order(bulk_oligo$sum_Oligo, decreasing = TRUE), , drop = FALSE]
bulk_tam <- bulk[, "sum_TAM", drop=FALSE]
bulk_tam <- bulk_tam[order(bulk_tam$sum_TAM, decreasing = TRUE), , drop = FALSE]
bulk_tumor <- bulk[, "sum_Tumor", drop=FALSE]
bulk_tumor <- bulk_tumor[order(bulk_tumor$sum_Tumor, decreasing = TRUE), , drop = FALSE]
bulk_div <- bulk[, "sum_Div_tumor", drop=FALSE]
bulk_div <- bulk_div[order(bulk_div$sum_Div_tumor, decreasing = TRUE), , drop = FALSE]
bulk_neuron <- bulk[, "sum_Neuron", drop=FALSE]
bulk_neuron <- bulk_neuron[order(bulk_neuron$sum_Neuron, decreasing = TRUE), , drop = FALSE]
bulk_opc <- bulk[, "sum_OPC", drop=FALSE]
bulk_opc <- bulk_opc[order(bulk_opc$sum_OPC, decreasing = TRUE), , drop = FALSE]
bulk_endo <- bulk[, "sum_Endo", drop=FALSE]
bulk_endo <- bulk_endo[order(bulk_endo$sum_Endo, decreasing = TRUE), , drop = FALSE]
bulk_per <- bulk[, "sum_Per", drop=FALSE]
bulk_per <- bulk_per[order(bulk_per$sum_Per, decreasing = TRUE), , drop = FALSE]
bulk_tcell <- bulk[, "sum_Tcell", drop=FALSE]
bulk_tcell <- bulk_tcell[order(bulk_tcell$sum_Tcell, decreasing = TRUE), , drop = FALSE]

# for specificity scores on combined
combined_astro <- specificity_scores_sample[, "sum_Astro", drop=FALSE] 
combined_astro <- combined_astro[order(combined_astro$sum_Astro, decreasing = TRUE), , drop = FALSE]
combined_oligo <- specificity_scores_sample[, "sum_Oligo", drop=FALSE]
combined_oligo <- combined_oligo[order(combined_oligo$sum_Oligo, decreasing = TRUE), , drop = FALSE]
combined_tam <- specificity_scores_sample[, "sum_TAM", drop=FALSE]
combined_tam <- combined_tam[order(combined_tam$sum_TAM, decreasing = TRUE), , drop = FALSE]
combined_tumor <- specificity_scores_sample[, "sum_Tumor", drop=FALSE]
combined_tumor <- combined_tumor[order(combined_tumor$sum_Tumor, decreasing = TRUE), , drop = FALSE]
combined_div <- specificity_scores_sample[, "sum_Div_tumor", drop=FALSE]
combined_div <- combined_div[order(combined_div$sum_Div_tumor, decreasing = TRUE), , drop = FALSE]
combined_neuron <- specificity_scores_sample[, "sum_Neuron", drop=FALSE]
combined_neuron <- combined_neuron[order(combined_neuron$sum_Neuron, decreasing = TRUE), , drop = FALSE]
combined_opc <- specificity_scores_sample[, "sum_OPC", drop=FALSE]
combined_opc <- combined_opc[order(combined_opc$sum_OPC, decreasing = TRUE), , drop = FALSE]
combined_endo <- specificity_scores_sample[, "sum_Endo", drop=FALSE]
combined_endo <- combined_endo[order(combined_endo$sum_Endo, decreasing = TRUE), , drop = FALSE]
combined_per <- specificity_scores_sample[, "sum_Per", drop=FALSE]
combined_per <- combined_per[order(combined_per$sum_Per, decreasing = TRUE), , drop = FALSE]
combined_tcell <- specificity_scores_sample[, "sum_Tcell", drop=FALSE]
combined_tcell <- combined_tcell[order(combined_tcell$sum_Tcell, decreasing = TRUE), , drop = FALSE]


# for specificity scores on individual sample
sample_astro <- specificity_scores_sample[, "Astro", drop=FALSE]
sample_astro <- sample_astro[order(sample_astro$Astro, decreasing = TRUE), , drop = FALSE]
sample_oligo <- specificity_scores_sample[, "Oligo", drop=FALSE]
sample_oligo <- sample_oligo[order(sample_oligo$Oligo, decreasing = TRUE), , drop = FALSE]
sample_tam <- specificity_scores_sample[, "TAM", drop=FALSE]
sample_tam <- sample_tam[order(sample_tam$TAM, decreasing = TRUE), , drop = FALSE]
sample_neuron <- specificity_scores_sample[, "Neuron", drop=FALSE]
sample_neuron <- sample_neuron[order(sample_neuron$Neuron, decreasing = TRUE), , drop = FALSE]
sample_tumor <- specificity_scores_sample[, "GTumor", drop=FALSE]
sample_tumor <- sample_tumor[order(sample_tumor$GTumor, decreasing = TRUE), , drop = FALSE]
sample_div <- specificity_scores_sample[, "GDiv.tumor", drop=FALSE]
sample_div <- sample_div[order(sample_div$GDiv.tumor, decreasing = TRUE), , drop = FALSE]
sample_endo <- specificity_scores_sample[, "Endo", drop=FALSE]
sample_endo <- sample_endo[order(sample_endo$Endo, decreasing = TRUE), , drop = FALSE]
sample_per <- specificity_scores_sample[, "Per", drop=FALSE]
sample_per <- sample_per[order(sample_per$Per, decreasing = TRUE), , drop = FALSE]
sample_opc <- specificity_scores_sample[, "OPC", drop=FALSE]
sample_opc <- sample_opc[order(sample_opc$OPC, decreasing = TRUE), , drop = FALSE]
sample_tcell <- specificity_scores_sample[, "Tcell", drop=FALSE]
sample_tcell <- sample_tcell[order(sample_tcell$Tcell, decreasing = TRUE), , drop = FALSE]




kenzie <- read_excel("validation/41598_2018_27293_MOESM2_ESM.xlsx", sheet = "top_human_specificity")
kenzie <- as.data.frame(kenzie)
colnames(kenzie) <- kenzie[2, ]
kenzie <- kenzie[-2, ]

celltypes_kenzie <- split(kenzie, kenzie$Celltype)

#extract and sort data on grand_mean
kenzie_ast <- celltypes_kenzie$ast
rownames(kenzie_ast) <- kenzie_ast$gene 
kenzie_ast$gene <- NULL 
kenzie_ast$grand_mean <- as.numeric(kenzie_ast$grand_mean)
kenzie_ast <- kenzie_ast[order(-kenzie_ast$grand_mean), ]

kenzie_oli <- celltypes_kenzie$oli
rownames(kenzie_oli) <- kenzie_oli$gene 
kenzie_oli$gene <- NULL 
kenzie_oli$grand_mean <- as.numeric(kenzie_oli$grand_mean)
kenzie_oli <- kenzie_oli[order(-kenzie_oli$grand_mean), ]


kenzie_neu <- celltypes_kenzie$neu
rownames(kenzie_neu) <- kenzie_neu$gene 
kenzie_neu$gene <- NULL 
kenzie_neu$grand_mean <- as.numeric(kenzie_neu$grand_mean)
kenzie_neu <- kenzie_neu[order(-kenzie_neu$grand_mean), ]

kenzie_mic <- celltypes_kenzie$mic
rownames(kenzie_mic) <- kenzie_mic$gene 
kenzie_mic$gene <- NULL 
kenzie_mic$grand_mean <- as.numeric(kenzie_mic$grand_mean)
kenzie_mic <- kenzie_mic[order(-kenzie_mic$grand_mean), ]

kenzie_end <- celltypes_kenzie$end
rownames(kenzie_end) <- kenzie_end$gene 
kenzie_end$gene <- NULL 
kenzie_end$grand_mean <- as.numeric(kenzie_end$grand_mean)
kenzie_end <- kenzie_end[order(-kenzie_end$grand_mean), ]

#overlap between list McKenzie and highest normalized counts combined
print(length(intersect(head(rownames(kenzie_ast), 200), head(rownames(bulk_astro), 200))))
print(length(intersect(head(rownames(kenzie_oli), 200), head(rownames(bulk_oligo), 200))))
print(length(intersect(head(rownames(kenzie_neu), 200), head(rownames(bulk_neuron), 200))))
print(length(intersect(head(rownames(kenzie_mic), 200), head(rownames(bulk_tam), 200))))
print(length(intersect(head(rownames(kenzie_end), 200), head(rownames(bulk_endo), 200))))

#overlap between list McKenzie and highest spec score individual sample
print(length(intersect(head(rownames(kenzie_ast), 200), head(rownames(sample_astro), 200))))
print(length(intersect(head(rownames(kenzie_oli), 200), head(rownames(sample_oligo), 200))))
print(length(intersect(head(rownames(kenzie_neu), 200), head(rownames(sample_neuron), 200))))
print(length(intersect(head(rownames(kenzie_mic), 200), head(rownames(sample_tam), 200))))
print(length(intersect(head(rownames(kenzie_end), 200), head(rownames(sample_endo), 200))))

#overlap between list McKenzie and highest spec score combined
print(length(intersect(head(rownames(kenzie_ast), 200), head(rownames(combined_astro), 200))))
print(length(intersect(head(rownames(kenzie_oli), 200), head(rownames(combined_oligo), 200))))
print(length(intersect(head(rownames(kenzie_neu), 200), head(rownames(combined_neuron), 200))))
print(length(intersect(head(rownames(kenzie_mic), 200), head(rownames(combined_tam), 200))))
print(length(intersect(head(rownames(kenzie_end), 200), head(rownames(combined_endo), 200))))