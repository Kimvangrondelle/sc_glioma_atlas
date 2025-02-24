#test aggregat expression and tau on this aggregated counts. 
object <- readRDS(file = "output/hijfte/hijfte-y.rds")
object <- readRDS(file = "output/hijfte/3pr.rds")

rownames(object@assays$RNA@layers$counts) <- rownames(object@assays$RNA@features@.Data)
head(object@assays$RNA@layers$counts, 10)
agg_object <- AggregateExpression(object, group.by = "celltype", return.seurat = TRUE)
Idents(agg_object) <- rownames(agg_object@assays$RNA@cells)


counts <- agg_object@assays$RNA@layers$counts
colnames(counts) <- rownames(agg_object@assays$RNA@cells@.Data)
rownames(counts) <- Features(agg_object)
counts <- as.data.frame(counts)
# 
# agg_object@assays$RNA@layers$counts <- ceiling(counts)
# counts <- ceiling(counts)
counts_vst <- counts %>%
  DESeq2::DESeqDataSetFromMatrix( data.frame(cond = as.factor(paste0('c',round(runif(ncol(.)))+1) )), ~cond) %>%
  DESeq2::vst(blind=T)
vst_counts <- as.data.frame(counts_vst@assays@data@listData)

tau <- apply(vst_counts, 1, calc_tau)


counts_tau <- data.frame(vst_counts, tau)

bayes <- read.csv("output/celltypes/combined/spec_bayes.csv", header = TRUE)
bayes1 <- read.csv("output/celltypes/combined/spec_tau.csv", header = TRUE)