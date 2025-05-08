Raw reads available at NCBI SRA project PRJNA1248589

trimmomatic.sh was used to trim and filter raw reads

trimmed reads were aligned with hisat 2:
    hisat2-build.sh was used to build index
    hisat2-paired.sh and hisat2-unpaired.sh were used to align reads to reference genome

sam-to-bam.sh was used to convert sam alignment files to bam format

htseq.sh was used to generate counts for all files

Count tables were downloaded to a local machine and all remaining analysis was performed in R (v. 4.4.1) and RStudio (v. 2024.04.2+764)

Filterandcombine.R was used to merge all counts files into a single table and combine counts of paired and unpaired read files for each sample

DEseq2.R was used to generate normalized read counts table and identify differentially expressed genes

gene_ontology.R was used to create list of genes with GO terms for stickleback and pull p values for differential expression for each gene

figures_tables.R was used to generate figures for publication and subset data into tables