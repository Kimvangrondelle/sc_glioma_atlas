library(Seurat)
library(tidyverse)
library(biomaRt)
library(infercnv)
#define sample ID and paths for output files
sid = "3pr"
path1 <- paste0("infercnv_out/hijfte-3pr/infercnv_",sid,"_out_pdf")
path2 <- paste0("infercnv_out/hijfte-3pr/infercnv_",sid,"_processed.Rds")
path3 <- paste0("infercnv_out/hijfte-3pr/infercnv_",sid,".Rds")
path4 <- paste("output/hijfte/3pr.rds")

#check if inferCNV output already exists, if not, generate it
if(!file.exists(path1)) {
  if(!file.exists(path2)) {
    if(!file.exists(path3)) {
      #load object
      object <- readRDS(file = path4)
      #create annotated cluster labels
      object$annotated_clusters <- as.factor(paste0(as.character(object$seurat_clusters), ". ", object$celltype))
      
      # rcm <- object@assays$RNA@layers$counts
      rcm <- object@assays$RNA@counts
      # str(rcm)
      # colnames(rcm) <- colnames(object)
      # rownames(rcm) <- rownames(object)
      colnames(rcm) <- object@assays$RNA@counts@Dimnames[[2]]
      rownames(rcm) <- object@assays$RNA@counts@Dimnames[[1]]
      # create annotation dataframe
      rcm <- as.data.frame(rcm)
      af <- data.frame(
        V1 = colnames(rcm)) |>
        dplyr::left_join(data.frame(V2 = object$annotated_clusters) |>
                           tibble::rownames_to_column("V1"), by=c('V1'='V1')) |> 
        tibble::column_to_rownames('V1') |> 
        dplyr::mutate(V2 = as.character(V2)) 
      #ensure column names and annotation rownames match
      stopifnot(colnames(rcm) == rownames(af))
      #keep only common cells
      common_cells <- intersect(colnames(rcm), rownames(af))
      rcm <- rcm[, common_cells, drop = FALSE]
      af <- af[common_cells, , drop = FALSE]
      
      identical(colnames(rcm), rownames(af))
      #define reference groups (non-malignant cell types)
      refgroups = data.frame(annotated_clusters = as.character(object$annotated_clusters),
                             celltype = as.character(object$celltype)) |> 
        dplyr::filter(celltype %in% c("TAM","Neuron","Tcell","Bcell", "Per","Endo","Oligo", "OPC", "Astro", "?", "O?")) |> 
        #dplyr::filter(celltype %in% c("Oligo", "TAM")) |> 
        dplyr::pull(annotated_clusters) |> 
        unique()
      refgroups <- as.character(refgroups)
      #retrieve gene information from ensemble
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      
      
      gene_info <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                         filters = "hgnc_symbol", values = rownames(rcm), mart = ensembl) %>% 
        dplyr::rename(chr = chromosome_name) %>% 
        dplyr::rename(start = start_position) %>% 
        dplyr::rename(end = end_position) %>% 
        dplyr::rename(gene = hgnc_symbol) %>% 
        dplyr::filter(grepl("^H(G|S)", chr) == F) %>% 
        dplyr::mutate(chr = paste0("chr",chr)) %>% 
        dplyr::filter(!duplicated(gene)) %>% 
        tibble::column_to_rownames('gene')
      #create infercnv object
      infercnv_obj = infercnv::CreateInfercnvObject(
        raw_counts_matrix= rcm,
        gene_order_file=gene_info,
        
        annotations_file = af,
        ref_group_names=refgroups # group names for only the non-malignants
      )
      #save infercnv object
      saveRDS(infercnv_obj, file=path3)
      # rm(rcm, af, infercnv_obj)
      # gc()
      
    } else {
      print(paste0("File: ", path3, " already present - skipping re-generation"))
    }
    
    #run infercnv analysis
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff=0.1,
                                  out_dir = path1,
                                  cluster_by_groups=TRUE,
                                  denoise=T,
                                  HMM=T,
                                  output_format="pdf")
    saveRDS(infercnv_obj, file=path2)
    # rm(infercnv_obj)
    # gc()
  }
}


#remove temporary files generated during infercnv processing
system(paste0("rm ",path1,"/*.txt"))
system(paste0("rm ",path1,"/*.dat"))
system(paste0("rm ",path1,"/*_obj"))
