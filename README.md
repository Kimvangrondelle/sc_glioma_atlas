# SC Glioma Atlas
Kim van Grondelle; september 2024 - february 2025

Master thesis SC Glioma Atlas for Bioinformatics and Systems Biology at VU Amsterdam. The research was conducted at the Neurology department of Erasmus MC Rotterdam. 

## Abstract
Gliomas are aggressive brain tumors characterized by a highly heterogeneous tumor microenvironment (TME). This complex tumor microenvironment plays critical roles in tumor progression and therapy resistance. Single-cell RNA sequencing is used to investigate the glioma microenvironment as this allows to analyze transcriptional profiles of individual cells.

The aim of the study is to determine and formulate the specificity of the expression for each gene for all glioma TME cell types. 'Which markers / genes are to which extent specific for which cell types present in the microenvironment of gliomas?' and 'Are there markers / genes specific for cell states, such as dividing or hypoxic?' are the research questions for this study. 

ScRNAseq datasets of three glioma types were analyzed and processed using Seurat. (1) To annotate the clusters found during processing, known marker genes were used. (2) Novel gene specificity scoring methods were developed to identify cell type-specific genes. (3) An external dataset was used to validate the different scoring methods used, including a newly developed method, a probabilistic scoring method, and DESeq2 to find differentially expressed genes. 

The results show per cell type the genes with the highest specificity score, together with the genes being the most differentially expressed per cell type. Genes found to be specific for certain cell types showed high specificity scores using the newly designed metric for these same cell types, confirming the robustness of the scoring approach. Furthermore, copy number variation analysis supported the annotation of tumor cells and highlighted intratumoral heterogeneity. 

The study provided a data resource as foundation for insights into the cellular composition and dynamics within gliomas, perhaps even guiding therapeutic strategies and enhancing our understanding of gliomas. 

## Structure of files
The structure consists of multiple directories and some single R files. 

### directories
- Annotation: Contains files used to annotate the clusters.
- Individual-process: Contains the files used to preprocess the samples individually to extract the parameters. 
- Process-neat: Contains the files to preprocess the samples per paper in a single function using the functions extracted from the individual process files. 
- Results_section_plots: Contains files with code to produce scatterplots used in the thesis and the lists of final result. 
- Sandbox: contains files used to test pieces of code
- Validation: Contains files used for different ways of validation
- Wies_glass_validatie: contains files used for validation using external data (correlation plot)

### files
- combining.R was used to combine all count matrices together into one big matrix. 
- infercnv.R was used to get a copy number profile of each sample. 
- pooling.R was used to pool the combined sample, to run DESeq2 and to extract the specificity scores. 
- pooling_celltypes.R was used to pool the individual samples, to run DESeq2 and to extract the specificity scores. 
- specificity_scoring_functions.R contains the functions used to calculate the specificity scores (these functions are used in the pooling files).

### output
- top_20_genes_DE_DESeq2.txt contains the 20 genes most differentially expressed per celltype. 
- top_20_genes_highest_specscore.txt contains the 20 genes with the highest specificity scores per cell type.
- Within the scores_combined_all directory, the files containing the scores for the zpex and bayes scores are stored. 

## Packages
- R -- 4.4.1 
- infercnv -- 1.21.0
- biomaRt -- 2.60.1  
- Matrix -- 1.7-1
- corrplot -- 0.95  
- Seurat -- 5.0.0 
- DESeq2 -- 1.44.0 
- tidyverse -- 2.0.0
- dplyr -- 1.1.4   
- cowplot -- 1.1.3
- ggplot2 -- 3.5.1  
- GridExtra -- 2.3
- ggpubr -- 0.6.0   
- patchwork -- 1.3.0
- ggrepel -- 0.9.6  


## Data used
The data that was used to conduct this research was coming from 5 different studies and was available in a 10X format: counts per gene per cell. 
- Couturier--EGAS00001004422 (4 samples): Couturier, C., Ayyadhury, S., & Le, e. a., PU. (2020). Single-cell rna-seq reveals that glioblstoma recapitulates a normal neurodevelopmental hierarchy. Nature communications. 
- Yuan--GSE103224 (4 samples): Yuan, J., Levitin, H., & Frattini, V. e. a. (2018). Single-cell transcriptome analysis of lineage diversity in high-grade glioma. Genome Med.
- Diaz--GSE138794 (9 samples): Wang, L., Babikir, H., Müller, S., Yagnik, G., Shamardani, K., Catalan, F., Kohanbash, G., Alvarado, B., Di Lullo, E., Kriegstein, A., Shah, S., Wadhwa, H., Chang, S., Phillips, J., Aghi, M., & Diaz, A. (2019). The phenotypes of proliferating glioblastoma cells reside on a single axis of variation. Cancer Discov..  
- 2 samples were used from; Hoogstrate, Y., Draaisma, K., & Ghisai, S. e. a. (2023). Transcriptome analysis reveals tumor microenvironment changes in glioblastoma. Cancer Cell.
- a single sample was used from: Ghisai, S., van Hijfte, L., & Vallentgoed, W. e. a. (2024). Epigenetic landscape reorganisation and reactivation of embryonic development genes are associated with malignancy in idh-mutant astrocytoma. Acta Neuropathologica.. This last sample, "Hijfte-3pr", consists of 6 samples integrated into one object.

## Workflow
![workflow](https://github.com/user-attachments/assets/fc8f8a05-35c9-41b9-9544-bdd202368317)
Figure 1: Workflow followed within the study. 

## How to use
In figure 1, the workflow followed within the study is shown. Below a description was added which files must be used to follow the workflow. 
1. "individual-process/" to preprocess the samples stored as 10X format to extract parameters.
1. "individual-process/" to actually preprocess the samples using the found parameters. Objects were stored in a convenient location. 
1. "annotation/marker_plot.R" to analyze the expression of marker genes in the found clusters. 
1. "annotation/annotations.R" to store the annotations of the clusters in the metadata of the objects. 
1. "infercnv.R" to extract copy number profiles to refine the initial annotations. -- annotations were updated again with "annotation/annotations.R"
1. "pooling/celltypes" to pool the individual samples on cell type. DESeq2 was run on the pooled object to extract the DE genes per cluster. On the pooled counts, the specificity metrics (from "specificity_scoring_functions.R") were applied to get the specificity scores per gene per cell type. -- Files with scores are stored, together with the 20 genes most differentially expressed per cell type. 
1. After all samples had been processed individually, "combining.R" was run to pool all samples together into one big matrix. 
1. "pooling.R" was run to extract the DE genes per cell type in the combined sample, and again the specificity scoring metrics were applied. Files with scores are stored, together with the 20 genes most differentially expressed per cell type.
1. "wies_glass_validation" was used to perform validation analysis on the extracted specificity scores.

Figure 2 shows the annotated UMAP visualization of sample Hijfte-y to give an idea of the analyses that were conducted within the research. 
![umap_annotated_hijfte-y](https://github.com/user-attachments/assets/644c25d7-fad0-4232-8462-cd2c080176c8)
Figure 2: Annotated UMAP visualization of sample Hijfte-y

## Output
The output of the sc-glioma-atlas are per gene the specificity scores for each cell type. These scores can be used to select genes specific for each cell type and can be used for gene set enrichment analysis to further investigate the genes and their functions. The two top-20 lists represent per cell type the 20 genes with the highest specificity score or the genes that were most differentially expressed between the cell types. 
