library(tidyverse)

#### Read in GO terms from gene2go table from NCBI

go_terms = read.delim(gzfile("gene2go.gz"), header = FALSE, sep = "\t", comment.char = "#",
                     col.names = c("Taxonomy ID", "Gene ID", "GO Term ID", "Evidence Code", "Qualifier", "GO Term", "PMID", "GO Category"))

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

#### save table for future reference or pull go terms from saved file
write.table(as.matrix(go_terms), file = "GO_terms_gac.tsv", quote = F, sep = "\t", row.names = F)

go_terms = read.table("GO_terms_gac.tsv", sep = "\t", header = T, fill = TRUE, quote = "\"")

#### read in list of degs from files
degs_17_auto = read.table("degs_17_auto.tsv", sep = "\t", row.names = 1, quote = "\"")
degs_20_auto = read.table("degs_20_auto.tsv", sep = "\t", row.names = 1, quote = "\"")
degs_23_auto = read.table("degs_23_auto.tsv", sep = "\t", row.names = 1, quote = "\"")
degs_d0_auto = read.table("degs_d0_auto.tsv", sep = "\t", row.names = 1, quote = "\"")
degs_d4_auto = read.table("degs_d4_auto.tsv", sep = "\t", row.names = 1, quote = "\"")
degs_d8_auto = read.table("degs_d8_auto.tsv", sep = "\t", row.names = 1, quote = "\"")

#### compile list of Gene IDs from all deg tables
degs = c(degs_17_auto$id, degs_20_auto$id, degs_23_auto$id, degs_d0_auto$id, degs_d4_auto$id, degs_d8_auto$id)
degs = unique(degs)

#### make table of GO terms for all available degs and record adjusted p value for significant stages
go_degs = go_terms[go_terms$Gene.ID %in% degs,]

s17 = s20 = s23 = d0 = d4 = d8 = matrix(nrow=0, ncol=1)
for (i in go_degs$Symbol){
  s17 = rbind(s17, pvals_17[i,]$padj)
  s20 = rbind(s20, pvals_20[i,]$padj)
  s23 = rbind(s23, pvals_23[i,]$padj)
  d0 = rbind(d0, pvals_d0[i,]$padj)
  d4 = rbind(d4, pvals_d4[i,]$padj)
  d8 = rbind(d8, pvals_d8[i,]$padj)
  }
go_degs = cbind(go_degs, s17, s20, s23, d0, d4, d8)
rm(i, s17, s20, s23, d0, d4, d8)

write.table(go_degs, file = "GO_degs.tsv", quote = F, sep = "\t", row.names = F)



######### Go enrichment with topGO

if (!requireNamespace("BiocManager", quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install()
BiocManager::install("topGO")

library(BiocManager)
library(topGO)

# make list of all GO terms in stickleback with p values for each stage

go_terms = read.table("GO_terms_gac.tsv", sep = "\t", header = T, fill = TRUE, quote = "\"")
go_terms$s17[1] = NA
go_terms$s20[1] = NA
go_terms$s23[1] = NA
go_terms$d0[1] = NA
go_terms$d4[1] = NA
go_terms$d8[1] = NA

get_padj = function(pvalues, geneid) {
  if (geneid %in% pvalues$id) {return(pvalues$padj[which(pvalues$id == geneid)])}
  else {return(NA)}
}

for (i in 1:nrow(go_terms)) {
  ID = go_terms$Gene.ID[i]
  go_terms$s17[i] = get_padj(pvals_17, ID)
  go_terms$s20[i] = get_padj(pvals_20, ID)
  go_terms$s23[i] = get_padj(pvals_23, ID)
  go_terms$d0[i] = get_padj(pvals_d0, ID)
  go_terms$d4[i] = get_padj(pvals_d4, ID)
  go_terms$d8[i] = get_padj(pvals_d8, ID)
}

write.table(go_terms, file = "GO_pvals.tsv", quote = F, sep = "\t", row.names = F)

### 