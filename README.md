# SC glioma atlas

## Abstract
Gliomas are aggressive brain tumors characterized by a highly heterogeneous tumor microenvironment (TME). This complex tumor microenvironment plays critical roles in tumor progression and therapy resistance. Single-cell RNA sequencing is used to investigate the glioma microenvironment as this allows to analyze transcriptional profiles of individual cells.

The aim of the study is to determine and formulate the specificity of the expression for each gene for all glioma TME cell types. 'Which markers / genes are to which extent specific for which cell types present in the microenvironment of gliomas?' and 'Are there markers / genes specific for cell states, such as dividing or hypoxic?' are the research questions for this study. 

ScRNAseq datasets of three glioma types were analyzed and processed using Seurat. (1) To annotate the clusters found during processing, known marker genes were used. (2) Novel gene specificity scoring methods were developed to identify cell type-specific genes. (3) An external dataset was used to validate the different scoring methods used, including a newly developed method, a probabilistic scoring method, and DESeq2 to find differentially expressed genes. 

The results show per cell type the genes with the highest specificity score, together with the genes being the most differentially expressed per cell type. Genes found to be specific for certain cell types showed high specificity scores using the newly designed metric for these same cell types, confirming the robustness of the scoring approach. Furthermore, copy number variation analysis supported the annotation of tumor cells and highlighted intratumoral heterogeneity. 

The study provided a data resource as foundation for insights into the cellular composition and dynamics within gliomas, perhaps even guiding therapeutic strategies and enhancing our understanding of gliomas. 

## Structure of files
The structure consists of multiple directories and some single R files. 
- Annotation: Contains files used to annotate the clusters.
- Individual-process: Contains the files used to preprocess the samples individually to extract the parameters. 
- Process-neat: Contains the files to preprocess the samples per paper in a single function using the functions extracted from the individual process files. 
- Results_section_plots: Contains files with code to produce scatterplots used in the thesis and the lists of final result. 
- Sandbox: contains files used to test pieces of code
- Validation: Contains files used for different ways of validation
- Wies_glass_validatie: contains files used for validation using external data (correlation plot)

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

## How to use
From raw data till scores. 