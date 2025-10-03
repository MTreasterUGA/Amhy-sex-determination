library(tidyverse)

################################################################################
##### Make stickleback GO table ################################################
################################################################################

### Read in GO terms from gene2go table from NCBI
go_terms = read.delim(gzfile("gene2go.gz"), header = FALSE, 
                      sep = "\t", comment.char = "#", 
                      col.names = c("Taxonomy ID", "Gene ID", "GO Term ID", 
                                    "Evidence Code", "Qualifier", "GO Term", 
                                    "PMID", "GO Category"))

### keep only go terms for threespine stickleback, ID 69293

go_terms = gene2go[gene2go$`Taxonomy ID` == 69293,]

### pull gene name, symbol, and chromosome from NCBI table of stickleback genes
genes_all = read_tsv("ncbi_genelist.tsv")

Name = character()
Symbol = character()
Chromosome = character()
for (i in go_terms$`Gene ID`){
  Name = c(Name, genes_all[match(i,genes_all$'Gene ID'),"Name"])
  Symbol = c(Symbol, genes_all[match(i,genes_all$'Gene ID'),"Symbol"])
  Chromosome = c(Chromosome, genes_all[match(i,genes_all$'Gene ID'),"Chromosome"])
}
Name = as.matrix(Name, ncol=1)
Chromosome = as.matrix(Chromosome, ncol=1)
Symbol = as.matrix(Symbol, ncol=1)
go_terms = cbind(go_terms, Name, Symbol, Chromosome)
rm(i, Name, Chromosome, Symbol)

go_terms = go_terms[, c(9, 10 , 11, 1, 2, 3, 4, 5, 6, 7, 8)]

### save table for future reference or 
write.table(as.matrix(go_terms), file = "GO_terms_gac.tsv", 
            quote = F, sep = "\t", row.names = F)

##### pull go terms from saved file ############################################
go_terms = read.table("GO_terms_gac.tsv", sep = "\t", header = T, 
                      fill = TRUE, quote = "\"")

################################################################################
##### complile go terms for degs ###############################################
################################################################################

### read in list of degs from file
degs = read.table("DEseq_degs_all_stages.tsv", header = T, 
                  quote = "\"", sep = "\t")

### make table of available GO terms for all degs
go_degs = go_terms[go_terms$Gene.ID %in% degs$Gene.ID,]

### add p values
go_degs = go_degs |>
  left_join(degs |> select(Gene.ID, s17, s20, s23, d0, d4, d8), by = "Gene.ID")

write.table(go_degs, file = "GO_degs.tsv", quote = F, sep = "\t", row.names = F)
