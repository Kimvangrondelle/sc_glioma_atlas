# loading the data
hijfte <- readRDS(file = "output/hijfte/hijfte-y.rds")
# making extra metadata column with again cluster numbers
hijfte$celltype <- as.character((hijfte$seurat_clusters))
# replace cluster numbers with the cell type
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(0, 1, 2, 3, 8, 10, 12), "GTumor", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(4, 6, 13, 18), "Oligo", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(5, 9, 14), "TAM", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(7, 11), "GDiv tumor", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(15), "Per", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(16), "Endo", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(17), "Astro", hijfte$celltype)
hijfte$celltype <- ifelse(hijfte$seurat_clusters %in% c(19), "Neuron", hijfte$celltype)
# plots to see original clusters and annotated clusters
hijfteclus <- DimPlot(hijfte, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Hijfte")
hijftecell <- DimPlot(hijfte, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Hijfte")

#save UMAPs as jpg -- both as un as annotated
ggsave(filename = "annotation/output/hijfte/hijfte-y-clus.jpg", height = 5, width = 10, plot = hijfteclus, quality = 50)
ggsave(filename = "annotation/output/hijfte/hijfte-y-cell.jpg", height = 5, width = 10, plot = hijftecell, quality = 50)
# save object
saveRDS(hijfte, file = "output/hijfte/hijfte-y.rds")

bolleboom <- readRDS(file = "output/bolleboom/H243.rds")
bolleboom$celltype <- as.character(bolleboom$seurat_clusters)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(0, 1, 19), "Oligo", bolleboom$celltype)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(2, 3, 4, 16), "GTumor", bolleboom$celltype)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(5, 6, 7, 8, 10, 12, 13, 14, 17, 18, 20), "Neuron", bolleboom$celltype)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(15), "Astro", bolleboom$celltype)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(21), "GDiv tumor", bolleboom$celltype)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(9), "OPC", bolleboom$celltype)
bolleboom$celltype <- ifelse(bolleboom$seurat_clusters %in% c(11), "TAM", bolleboom$celltype)


bolclus <- DimPlot(bolleboom, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Bolleboom")
bolcell <- DimPlot(bolleboom, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Bolleboom")

ggsave(filename = "annotation/output/bolleboom-clus.jpg", height = 5, width = 10, plot = bolclus, quality = 50)
ggsave(filename = "annotation/output/bolleboom-cell.jpg", height = 5, width = 10, plot = bolcell, quality = 50)
saveRDS(bolleboom, file = "output/bolleboom/H243.rds")

diaz10022 <- readRDS(file = "output/diaz/diaz10022.rds")
diaz10022$celltype <- as.character(diaz10022$seurat_clusters)
diaz10022$celltype <- ifelse(diaz10022$seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 7, 12, 15), "GTumor", diaz10022$celltype)
diaz10022$celltype <- ifelse(diaz10022$seurat_clusters %in% c(6), "Oligo", diaz10022$celltype)
diaz10022$celltype <- ifelse(diaz10022$seurat_clusters %in% c(9, 13), "TAM", diaz10022$celltype)
diaz10022$celltype <- ifelse(diaz10022$seurat_clusters %in% c(10), "GDiv tumor", diaz10022$celltype)
diaz10022$celltype <- ifelse(diaz10022$seurat_clusters %in% c(14), "Neuron", diaz10022$celltype)
diaz10022$celltype <- ifelse(diaz10022$seurat_clusters %in% c(8, 11), "OPC", diaz10022$celltype)

clus10022 <- DimPlot(diaz10022, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Diaz10022")
cell10022 <- DimPlot(diaz10022, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Diaz10022")

ggsave(filename = "annotation/output/diaz/diaz10022-clus.jpg", height = 5, width = 10, plot = clus10022, quality = 50)
ggsave(filename = "annotation/output/diaz/diaz10022-cell.jpg", height = 5, width = 10, plot = cell10022, quality = 50)
saveRDS(diaz10022, file = "output/diaz/diaz10022.rds")

diaz11644 <- readRDS(file = "output/diaz/diaz11644.rds")
diaz11644$celltype <- as.character(diaz11644$seurat_clusters)
diaz11644$celltype <- ifelse(diaz11644$seurat_clusters %in% c(0, 2, 3, 5, 6, 7, 8, 12, 13), "GTumor", diaz11644$celltype)
diaz11644$celltype <- ifelse(diaz11644$seurat_clusters %in% c(14), "Oligo", diaz11644$celltype)
diaz11644$celltype <- ifelse(diaz11644$seurat_clusters %in% c(1, 4, 10, 11), "TAM", diaz11644$celltype)
diaz11644$celltype <- ifelse(diaz11644$seurat_clusters %in% c(9), "GDiv tumor", diaz11644$celltype)
diaz11644$celltype <- ifelse(diaz11644$seurat_clusters %in% c(15), "Endo", diaz11644$celltype)

clus11644 <- DimPlot(diaz11644, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Diaz11644")
cell11644 <- DimPlot(diaz11644, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Diaz11644")

ggsave(filename = "annotation/output/diaz/diaz11644-clus.jpg", height = 5, width = 10, plot = clus11644, quality = 50)
ggsave(filename = "annotation/output/diaz/diaz11644-cell.jpg", height = 5, width = 10, plot = cell11644, quality = 50)
saveRDS(diaz11644, file = "output/diaz/diaz11644.rds")

diaz11681 <- readRDS(file = "output/diaz/diaz11681.rds")
diaz11681$celltype <- as.character(diaz11681$seurat_clusters)
diaz11681$celltype <- ifelse(diaz11681$seurat_clusters %in% c(1, 2, 6, 9), "GTumor", diaz11681$celltype)
diaz11681$celltype <- ifelse(diaz11681$seurat_clusters %in% c(7), "Oligo", diaz11681$celltype)
diaz11681$celltype <- ifelse(diaz11681$seurat_clusters %in% c(0, 4, 5, 10), "TAM", diaz11681$celltype)
diaz11681$celltype <- ifelse(diaz11681$seurat_clusters %in% c(3), "GDiv tumor", diaz11681$celltype)
diaz11681$celltype <- ifelse(diaz11681$seurat_clusters %in% c(8), "Tcell", diaz11681$celltype)

clus11681 <- DimPlot(diaz11681, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Diaz11681")
cell11681 <- DimPlot(diaz11681, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Diaz11681")

ggsave(filename = "annotation/output/diaz/diaz11681-clus.jpg", height = 5, width = 10, plot = clus11681, quality = 50)
ggsave(filename = "annotation/output/diaz/diaz11681-cell.jpg", height = 5, width = 10, plot = cell11681, quality = 50)
saveRDS(diaz11681, file = "output/diaz/diaz11681.rds")

diaz12264 <- readRDS(file = "output/diaz/diaz12264.rds")
diaz12264$celltype <- as.character(diaz12264$seurat_clusters)
diaz12264$celltype <- ifelse(diaz12264$seurat_clusters %in% c(0, 4, 5, 6, 7, 8, 10), "GTumor", diaz12264$celltype)
diaz12264$celltype <- ifelse(diaz12264$seurat_clusters %in% c(2, 3), "Oligo", diaz12264$celltype)
diaz12264$celltype <- ifelse(diaz12264$seurat_clusters %in% c(1, 9, 11), "TAM", diaz12264$celltype)
diaz12264$celltype <- ifelse(diaz12264$seurat_clusters %in% c(12), "Tcell", diaz12264$celltype)

clus12264 <- DimPlot(diaz12264, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Diaz12264")
cell12264 <- DimPlot(diaz12264, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Diaz12264")

ggsave(filename = "annotation/output/diaz/diaz12264-clus.jpg", height = 5, width = 10, plot = clus12264, quality = 50)
ggsave(filename = "annotation/output/diaz/diaz12264-cell.jpg", height = 5, width = 10, plot = cell12264, quality = 50)
saveRDS(diaz12264, file = "output/diaz/diaz12264.rds")

diaz4297 <- readRDS(file = "output/diaz/diaz4297.rds")
diaz4297$celltype <- as.character(diaz4297$seurat_clusters)
diaz4297$celltype <- ifelse(diaz4297$seurat_clusters %in% c(0, 2, 3, 5, 6, 7, 9, 10, 12), "GTumor", diaz4297$celltype)
diaz4297$celltype <- ifelse(diaz4297$seurat_clusters %in% c(8), "Oligo", diaz4297$celltype)
diaz4297$celltype <- ifelse(diaz4297$seurat_clusters %in% c(1, 4, 11), "TAM", diaz4297$celltype)

clus4297 <- DimPlot(diaz4297, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters Diaz4297")
cell4297 <- DimPlot(diaz4297, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes Diaz4297")

ggsave(filename = "annotation/output/diaz/diaz4297-clus.jpg", height = 5, width = 10, plot = clus4297, quality = 50)
ggsave(filename = "annotation/output/diaz/diaz4297-cell.jpg", height = 5, width = 10, plot = cell4297, quality = 50)
saveRDS(diaz4297, file = "output/diaz/diaz4297.rds")

coutu338 <- readRDS(file = "output/couturier/couturier338.rds")
coutu338$celltype <- as.character(coutu338$seurat_clusters)
coutu338$celltype <- ifelse(coutu338$seurat_clusters %in% c(0, 3, 4, 5, 6, 9), "GTumor", coutu338$celltype)
coutu338$celltype <- ifelse(coutu338$seurat_clusters %in% c(8), "Oligo", coutu338$celltype)
coutu338$celltype <- ifelse(coutu338$seurat_clusters %in% c(10), "TAM", coutu338$celltype)
coutu338$celltype <- ifelse(coutu338$seurat_clusters %in% c(1, 2), "Per", coutu338$celltype)
coutu338$celltype <- ifelse(coutu338$seurat_clusters %in% c(7), "GDiv tumor", coutu338$celltype)

clus338 <- DimPlot(coutu338, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters coutu338")
cell338 <- DimPlot(coutu338, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes coutu338")

ggsave(filename = "annotation/output/couturier/coutu338-clus.jpg", height = 5, width = 10, plot = clus338, quality = 50)
ggsave(filename = "annotation/output/couturier/coutu338-cell.jpg", height = 5, width = 10, plot = cell338, quality = 50)
saveRDS(coutu338, file = "output/couturier/couturier338.rds")

coutu363 <- readRDS(file = "output/couturier/couturier363.rds")
coutu363$celltype <- as.character(coutu363$seurat_clusters)
coutu363$celltype <- ifelse(coutu363$seurat_clusters %in% c(0, 2, 3, 4, 6, 7), "GTumor", coutu363$celltype)
coutu363$celltype <- ifelse(coutu363$seurat_clusters %in% c(5, 8), "Oligo", coutu363$celltype)
coutu363$celltype <- ifelse(coutu363$seurat_clusters %in% c(1), "TAM", coutu363$celltype)
coutu363$celltype <- ifelse(coutu363$seurat_clusters %in% c(10), "Endo", coutu363$celltype)
coutu363$celltype <- ifelse(coutu363$seurat_clusters %in% c(9), "GDiv tumor", coutu363$celltype)

clus363 <- DimPlot(coutu363, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters coutu363")
cell363 <- DimPlot(coutu363, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes coutu363")

ggsave(filename = "annotation/output/couturier/coutu363-clus.jpg", height = 5, width = 10, plot = clus363, quality = 50)
ggsave(filename = "annotation/output/couturier/coutu363-cell.jpg", height = 5, width = 10, plot = cell363, quality = 50)
saveRDS(coutu363, file = "output/couturier/couturier363.rds")

coutu364 <- readRDS(file = "output/couturier/couturier364.rds")
coutu364$celltype <- as.character(coutu364$seurat_clusters)
coutu364$celltype <- ifelse(coutu364$seurat_clusters %in% c(0, 1, 4, 5, 7, 8, 10, 14), "GTumor", coutu364$celltype)
coutu364$celltype <- ifelse(coutu364$seurat_clusters %in% c(2, 12), "Oligo", coutu364$celltype)
coutu364$celltype <- ifelse(coutu364$seurat_clusters %in% c(3, 6, 9, 11), "TAM", coutu364$celltype)
coutu364$celltype <- ifelse(coutu364$seurat_clusters %in% c(13), "GDiv tumor", coutu364$celltype)

clus364 <- DimPlot(coutu364, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters coutu364")
cell364 <- DimPlot(coutu364, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes coutu364")

ggsave(filename = "annotation/output/couturier/coutu364-clus.jpg", height = 5, width = 10, plot = clus364, quality = 50)
ggsave(filename = "annotation/output/couturier/coutu364-cell.jpg", height = 5, width = 10, plot = cell364, quality = 50)
saveRDS(coutu364, file = "output/couturier/couturier364.rds")

coutu397 <- readRDS(file = "output/couturier/couturier397.rds")
coutu397$celltype <- as.character(coutu397$seurat_clusters)
coutu397$celltype <- ifelse(coutu397$seurat_clusters %in% c(5), "GTumor", coutu397$celltype)
coutu397$celltype <- ifelse(coutu397$seurat_clusters %in% c(10), "Oligo", coutu397$celltype)
coutu397$celltype <- ifelse(coutu397$seurat_clusters %in% c(0, 1, 2, 3, 4, 6, 7, 9), "TAM", coutu397$celltype)
coutu397$celltype <- ifelse(coutu397$seurat_clusters %in% c(8), "GDiv tumor", coutu397$celltype)
coutu397$celltype <- ifelse(coutu397$seurat_clusters %in% c(11), "Tcell", coutu397$celltype)

clus397 <- DimPlot(coutu397, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters coutu397")
cell397 <- DimPlot(coutu397, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes coutu397")

ggsave(filename = "annotation/output/couturier/coutu397-clus.jpg", height = 5, width = 10, plot = clus397, quality = 50)
ggsave(filename = "annotation/output/couturier/coutu397-cell.jpg", height = 5, width = 10, plot = cell397, quality = 50)
saveRDS(coutu397, file = "output/couturier/couturier397.rds")

yuan18 <- readRDS(file = "output/yuan/yuan18.rds")
yuan18$celltype <- as.character(yuan18$seurat_clusters)
yuan18$celltype <- ifelse(yuan18$seurat_clusters %in% c(0, 1, 2, 3, 5, 6), "GTumor", yuan18$celltype)
yuan18$celltype <- ifelse(yuan18$seurat_clusters %in% c(4), "GDiv tumor", yuan18$celltype)
yuan18$celltype <- ifelse(yuan18$seurat_clusters %in% c(7), "TAM", yuan18$celltype)
yuan18$celltype <- ifelse(yuan18$seurat_clusters %in% c(8), "Endo", yuan18$celltype)
yuan18$celltype <- ifelse(yuan18$seurat_clusters %in% c(9), "Oligo", yuan18$celltype)

clus18 <- DimPlot(yuan18, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters yuan18")
cell18 <- DimPlot(yuan18, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes yuan18")

ggsave(filename = "annotation/output/yuan/yuan18-clus.jpg", height = 5, width = 10, plot = clus18, quality = 50)
ggsave(filename = "annotation/output/yuan/yuan18-cell.jpg", height = 5, width = 10, plot = cell18, quality = 50)
saveRDS(yuan18, file = "output/yuan/yuan18.rds")

yuan25 <- readRDS(file = "output/yuan/yuan25.rds")
yuan25$celltype <- as.character(yuan25$seurat_clusters)
yuan25$celltype <- ifelse(yuan25$seurat_clusters %in% c(0, 2, 3, 4, 6, 8), "GTumor", yuan25$celltype)
yuan25$celltype <- ifelse(yuan25$seurat_clusters %in% c(1), "GDiv tumor", yuan25$celltype)
yuan25$celltype <- ifelse(yuan25$seurat_clusters %in% c(9), "TAM", yuan25$celltype)
yuan25$celltype <- ifelse(yuan25$seurat_clusters %in% c(5), "Per", yuan25$celltype)
yuan25$celltype <- ifelse(yuan25$seurat_clusters %in% c(7), "Endo", yuan25$celltype)

clus25 <- DimPlot(yuan25, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters yuan25")
cell25 <- DimPlot(yuan25, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes yuan25")

ggsave(filename = "annotation/output/yuan/yuan25-clus.jpg", height = 5, width = 10, plot = clus25, quality = 50)
ggsave(filename = "annotation/output/yuan/yuan25-cell.jpg", height = 5, width = 10, plot = cell25, quality = 50)
saveRDS(yuan25, file = "output/yuan/yuan25.rds")

yuan35 <- readRDS(file = "output/yuan/yuan35.rds")
yuan35$celltype <- as.character(yuan35$seurat_clusters)
yuan35$celltype <- ifelse(yuan35$seurat_clusters %in% c(0, 1, 2, 3, 5), "GTumor", yuan35$celltype)
yuan35$celltype <- ifelse(yuan35$seurat_clusters %in% c(6), "GDiv tumor", yuan35$celltype)
yuan35$celltype <- ifelse(yuan35$seurat_clusters %in% c(4), "TAM", yuan35$celltype)
yuan35$celltype <- ifelse(yuan35$seurat_clusters %in% c(7), "Endo", yuan35$celltype)
yuan35$celltype <- ifelse(yuan35$seurat_clusters %in% c(8), "Per", yuan35$celltype)


clus35 <- DimPlot(yuan35, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters yuan35")
cell35 <- DimPlot(yuan35, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes yuan35")

ggsave(filename = "annotation/output/yuan/yuan35-clus.jpg", height = 5, width = 10, plot = clus35, quality = 50)
ggsave(filename = "annotation/output/yuan/yuan35-cell.jpg", height = 5, width = 10, plot = cell35, quality = 50)
saveRDS(yuan35, file = "output/yuan/yuan35.rds")


yuan48 <- readRDS(file = "output/yuan/yuan48.rds")
yuan48$celltype <- as.character(yuan48$seurat_clusters)
yuan48$celltype <- ifelse(yuan48$seurat_clusters %in% c(0, 1, 2, 3, 6, 7, 8), "GTumor", yuan48$celltype)
yuan48$celltype <- ifelse(yuan48$seurat_clusters %in% c(4), "GDiv tumor", yuan48$celltype)
yuan48$celltype <- ifelse(yuan48$seurat_clusters %in% c(5, 9), "OPC", yuan48$celltype)

clus48 <- DimPlot(yuan48, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters yuan48")
cell48 <- DimPlot(yuan48, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes yuan48")

ggsave(filename = "annotation/output/yuan/yuan48-clus.jpg", height = 5, width = 10, plot = clus48, quality = 50)
ggsave(filename = "annotation/output/yuan/yuan48-cell.jpg", height = 5, width = 10, plot = cell48, quality = 50)
saveRDS(yuan48, file = "output/yuan/yuan48.rds")


pr3 <- readRDS(file = "output/hijfte/3pr.rds")
pr3$celltype <- as.character(pr3$seurat_clusters)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(0, 1, 7, 8, 10, 14, 15, 16, 17, 20, 21), "ATumor", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(2, 6), "OPC", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(3, 12, 18, 22), "TAM", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(4, 9, 11, 27), "Astro", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(5, 23), "ADiv tumor", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(13, 19, 24), "Oligo", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(25, 26), "Neuron", pr3$celltype)
pr3$celltype <- ifelse(pr3$seurat_clusters %in% c(28), "Endo", pr3$celltype)

cluspr <- DimPlot(pr3, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters pr3")
cellpr <- DimPlot(pr3, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes pr3")

ggsave(filename = "annotation/output/hijfte/pr3-clus.jpg", height = 5, width = 10, plot = cluspr, quality = 50)
ggsave(filename = "annotation/output/hijfte/pr3-cell.jpg", height = 5, width = 10, plot = cellpr, quality = 50)
saveRDS(pr3, file = "output/hijfte/3pr.rds")


diaz11136 <- readRDS(file = "output/diaz_astrooligo/diaz11136_astro.rds")
diaz11136$celltype <- as.character(diaz11136$seurat_clusters)
diaz11136$celltype <- ifelse(diaz11136$seurat_clusters %in% c(3, 4, 5), "ATumor", diaz11136$celltype)
diaz11136$celltype <- ifelse(diaz11136$seurat_clusters %in% c(8), "TAM", diaz11136$celltype)
diaz11136$celltype <- ifelse(diaz11136$seurat_clusters %in% c(0, 1, 2, 6, 7), "Oligo", diaz11136$celltype)
diaz11136$celltype <- ifelse(diaz11136$seurat_clusters %in% c(9), "Endo", diaz11136$celltype)

clus11136 <- DimPlot(diaz11136, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters diaz11136")
cell11136 <- DimPlot(diaz11136, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes diaz11136")

ggsave(filename = "annotation/output/diaz_astrooligo/diaz11136-clus.jpg", height = 5, width = 10, plot = clus11136, quality = 50)
ggsave(filename = "annotation/output/diaz_astrooligo/diaz11136-cell.jpg", height = 5, width = 10, plot = cell11136, quality = 50)
saveRDS(diaz11136, file = "output/diaz_astrooligo/diaz11136_astro.rds")


diaz12017 <- readRDS(file = "output/diaz_astrooligo/diaz12017_astro.rds")
diaz12017$celltype <- as.character(diaz12017$seurat_clusters)
diaz12017$celltype <- ifelse(diaz12017$seurat_clusters %in% c(1, 2, 4, 5, 6, 9), "ATumor", diaz12017$celltype)
diaz12017$celltype <- ifelse(diaz12017$seurat_clusters %in% c(0, 8), "TAM", diaz12017$celltype)
diaz12017$celltype <- ifelse(diaz12017$seurat_clusters %in% c(3, 7), "Oligo", diaz12017$celltype)
diaz12017$celltype <- ifelse(diaz12017$seurat_clusters %in% c(11), "Per", diaz12017$celltype)
diaz12017$celltype <- ifelse(diaz12017$seurat_clusters %in% c(10), "ADiv tumor", diaz12017$celltype)

clus12017 <- DimPlot(diaz12017, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters diaz12017")
cell12017 <- DimPlot(diaz12017, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes diaz12017")

ggsave(filename = "annotation/output/diaz_astrooligo/diaz12017-clus.jpg", height = 5, width = 10, plot = clus12017, quality = 50)
ggsave(filename = "annotation/output/diaz_astrooligo/diaz12017-cell.jpg", height = 5, width = 10, plot = cell12017, quality = 50)
saveRDS(diaz12017, file = "output/diaz_astrooligo/diaz12017_astro.rds")


diaz11612 <- readRDS(file = "output/diaz_astrooligo/diaz11612_oligo.rds")
diaz11612$celltype <- as.character(diaz11612$seurat_clusters)
diaz11612$celltype <- ifelse(diaz11612$seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6), "OTumor", diaz11612$celltype)
diaz11612$celltype <- ifelse(diaz11612$seurat_clusters %in% c(7), "Oligo", diaz11612$celltype)
diaz11612$celltype <- ifelse(diaz11612$seurat_clusters %in% c(8, 9), "O?", diaz11612$celltype)

clus11612 <- DimPlot(diaz11612, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters diaz11612")
cell11612 <- DimPlot(diaz11612, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes diaz11612")

ggsave(filename = "annotation/output/diaz_astrooligo/diaz11612-clus.jpg", height = 5, width = 10, plot = clus11612, quality = 50)
ggsave(filename = "annotation/output/diaz_astrooligo/diaz11612-cell.jpg", height = 5, width = 10, plot = cell11612, quality = 50)
saveRDS(diaz11612, file = "output/diaz_astrooligo/diaz11612_oligo.rds")

diaz11949 <- readRDS(file = "output/diaz_astrooligo/diaz11949_oligo.rds")
diaz11949$celltype <- as.character(diaz11949$seurat_clusters)
diaz11949$celltype <- ifelse(diaz11949$seurat_clusters %in% c(0, 1, 4), "OTumor", diaz11949$celltype)
diaz11949$celltype <- ifelse(diaz11949$seurat_clusters %in% c(5), "Oligo", diaz11949$celltype)
diaz11949$celltype <- ifelse(diaz11949$seurat_clusters %in% c(3), "OPC", diaz11949$celltype)
diaz11949$celltype <- ifelse(diaz11949$seurat_clusters %in% c(2), "O?", diaz11949$celltype)

clus11949 <- DimPlot(diaz11949, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters diaz11949")
cell11949 <- DimPlot(diaz11949, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes diaz11949")

ggsave(filename = "annotation/output/diaz_astrooligo/diaz11949-clus.jpg", height = 5, width = 10, plot = clus11949, quality = 50)
ggsave(filename = "annotation/output/diaz_astrooligo/diaz11949-cell.jpg", height = 5, width = 10, plot = cell11949, quality = 50)
saveRDS(diaz11949, file = "output/diaz_astrooligo/diaz11949_oligo.rds")



#integrated files -- not used
vent <- readRDS(file = "output/integrated/vent_processed.rds")
vent$celltype <- as.character(vent$seurat_clusters)
vent$celltype <- ifelse(vent$seurat_clusters %in% c(0, 2, 3, 4, 5, 6, 7, 8, 9, 12), "Tumor", vent$celltype)
vent$celltype <- ifelse(vent$seurat_clusters %in% c(1, 10), "TAM", vent$celltype)
vent$celltype <- ifelse(vent$seurat_clusters %in% c(11), "Div Tumor", vent$celltype)
vent$celltype <- ifelse(vent$seurat_clusters %in% c(13), "Olig", vent$celltype)

clusvent <- DimPlot(vent, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters vent")
cellvent <- DimPlot(vent, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes vent")

ggsave(filename = "annotation/output/venteicher-clus.jpg", height = 5, width = 10, plot = clusvent, quality = 50)
ggsave(filename = "annotation/output/venteicher-cell.jpg", height = 5, width = 10, plot = cellvent, quality = 50)
saveRDS(vent, file = "output/venteicher/vent_processed.rds")



diaz_in <- readRDS(file = "output/integrated/integrated_diaz.rds")
diaz_in$celltype <- as.character(diaz_in$seurat_clusters)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(0, 6, 10, 12, 13, 14, 15, 16, 17, 19, 22, 24, 25), "Tumor", diaz_in$celltype)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(1, 11, 18, 20, 27), "TAM", diaz_in$celltype)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(8, 23), "Div Tumor", diaz_in$celltype)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(3, 5, 9, 17), "Olig", diaz_in$celltype)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(2, 4, 21), "Neuron", diaz_in$celltype)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(7), "Astro", diaz_in$celltype)
diaz_in$celltype <- ifelse(diaz_in$seurat_clusters %in% c(26), "Endoper", diaz_in$celltype)

clusdiaz_in <- DimPlot(diaz_in, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters diaz_in")
celldiaz_in <- DimPlot(diaz_in, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes diaz_in")

ggsave(filename = "annotation/output/diaz-in-clus.jpg", height = 5, width = 10, plot = clusdiaz_in, quality = 50)
ggsave(filename = "annotation/output/diaz-in-cell.jpg", height = 5, width = 10, plot = celldiaz_in, quality = 50)


coutu_in <- readRDS(file = "output/integrated/integrated_couturier.rds")
coutu_in$celltype <- as.character(coutu_in$seurat_clusters)
coutu_in$celltype <- ifelse(coutu_in$seurat_clusters %in% c(1, 2, 3, 4, 6, 8, 12, 14, 16), "Tumor", coutu_in$celltype)
coutu_in$celltype <- ifelse(coutu_in$seurat_clusters %in% c(0, 5, 7, 9, 15, 20), "TAM", coutu_in$celltype)
coutu_in$celltype <- ifelse(coutu_in$seurat_clusters %in% c(17, 18), "Div Tumor", coutu_in$celltype)
coutu_in$celltype <- ifelse(coutu_in$seurat_clusters %in% c(11, 13, 19), "Olig", coutu_in$celltype)
coutu_in$celltype <- ifelse(coutu_in$seurat_clusters %in% c(10), "Per", coutu_in$celltype)
coutu_in$celltype <- ifelse(coutu_in$seurat_clusters %in% c(21), "Endoper", coutu_in$celltype)

cluscoutu_in <- DimPlot(coutu_in, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters coutu_in")
cellcoutu_in <- DimPlot(coutu_in, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes coutu_in")

ggsave(filename = "annotation/output/coutu-in-clus.jpg", height = 5, width = 10, plot = cluscoutu_in, quality = 50)
ggsave(filename = "annotation/output/coutu-in-cell.jpg", height = 5, width = 10, plot = cellcoutu_in, quality = 50)


yuan_in <- readRDS(file = "output/integrated/integrated_yuan.rds")
yuan_in$celltype <- as.character(yuan_in$seurat_clusters)
yuan_in$celltype <- ifelse(yuan_in$seurat_clusters %in% c(0, 1, 4, 5, 6, 7, 8, 11, 14), "Tumor", yuan_in$celltype)
yuan_in$celltype <- ifelse(yuan_in$seurat_clusters %in% c(9), "TAM", yuan_in$celltype)
yuan_in$celltype <- ifelse(yuan_in$seurat_clusters %in% c(2, 3, 10), "Div Tumor", yuan_in$celltype)
yuan_in$celltype <- ifelse(yuan_in$seurat_clusters %in% c(15), "Olig", yuan_in$celltype)
yuan_in$celltype <- ifelse(yuan_in$seurat_clusters %in% c(12), "Per", yuan_in$celltype)
yuan_in$celltype <- ifelse(yuan_in$seurat_clusters %in% c(13), "Endoper", yuan_in$celltype)

clusyuan_in <- DimPlot(yuan_in, reduction = "umap", label = TRUE, pt.size = .4, group.by = "seurat_clusters") + labs(title = "clusters yuan_in")
cellyuan_in <- DimPlot(yuan_in, reduction = "umap", label = TRUE, pt.size = .4, group.by = "celltype") + labs(title = "celltypes yuan_in")

ggsave(filename = "annotation/output/yuan-in-clus.jpg", height = 5, width = 10, plot = clusyuan_in, quality = 50)
ggsave(filename = "annotation/output/yuan-in-cell.jpg", height = 5, width = 10, plot = cellyuan_in, quality = 50)
