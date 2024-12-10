# need CellRanger installed on system
# will result in the matrix, features and barcodes files 
# nog ff verder uitzoeken hoe dit werkt
# schijnt lang te duren bij een file 
# reference is GRCh38 

cellranger count --id=run_test_y_levi \
--fastqs=/mnt/home/user.name/yard/run_cellranger_count/pbmc_1k_v3_fastqs \
--sample=pbmc_1k_v3 \
--transcriptome=/mnt/home/user.name/yard/run_cellranger_count/refdata-gex-GRCh38-2020-A


cellranger count --id=run_test_y_levi  --fastqs=/home/r105682/yard/fastq_files_y_levi --transcriptome=/home/r105682/yard/apps/refdata-gex-GRCh38-2024-A --create-bam=false

cd cellranger-9.0.0/
  export PATH=/home/r105682/yard/apps/cellranger-9.0.0:$PATH