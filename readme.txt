Raw reads available at NCBI SRA project PRJNA1248589
Raw and normalized read counts are available at NCBI GEO GSE296766

trimmomatic.sh was used to trim and filter raw reads

hisat2-build.sh was used to build index from reference genome (NCBI assembly GCF_016920845.1)

hisat2-align.sh was used to align reads to reference genome

htseq.sh was used to generate counts for all files

Count tables were downloaded to a local device and all remaining analysis was performed in R (v. 4.4.1) and RStudio (v. 2024.04.2+764)

Filterandcombine.R was used to merge all counts files into a single table and combine counts of paired and unpaired read files for each sample

DEseq2.R was used to generate normalized read counts table and identify differentially expressed genes

gene_ontology.R was used to create list of genes with GO terms for stickleback and pull p values for differential expression for each gene

figures_tables_stas.R was used to subset data into tables, generate figures for publication, and perform various statistical analyses