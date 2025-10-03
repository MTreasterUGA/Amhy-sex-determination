if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

{
  library(DESeq2)
  library(tidyverse)
}

################################################################################
##### count data and coldata for deseq2 ########################################
################################################################################

### file with count data
input = "counts_raw_summed.tsv"

### output of coldata file
output = "counts_raw_summed_coldata.txt"

counts_raw = read.table(input, header = T, row.names = 1)
samples = colnames(counts_raw)
coldata = data.frame(counts_raw = samples)

##### sex and stage as collapsed factor
### pull group from first three characters of sample name
group = vector()
for (i in 1:nrow(coldata)){
  group = c(group, substr(coldata[i,1],1,3))
}

coldata = cbind(coldata,group)
coldata = column_to_rownames(coldata,var="counts_raw")

write.table(coldata,file=output,quote=F,row.names=T)

# variable cleanup
rm(group, i, input, output, samples)

##### manually read in coldata and/or counts tables ############################

counts_raw = read.table("counts_raw_summed.tsv")
coldata = read.table("counts_raw_summed_coldata.txt",header=T,row.names=1)

################################################################################
##### DEG analysis #############################################################
################################################################################

### check that names and order from counts and coldata table match
if(identical(rownames(coldata), colnames(counts_raw))) {
  print("samples match")
} else(print("samples in counts and coldata tables DO NOT match"))

### make deseq data set object (dds)

coldata$group = factor(coldata$group)
dds = DESeqDataSetFromMatrix(countData = counts_raw, 
                             colData = coldata, 
                             design = ~ group) # sex and stage collapsed

### remove all genes with fewer than 10 reads across all samples
dds = dds[rowSums(counts(dds)) >= 10,]

### generate normalized count table
dds = estimateSizeFactors(dds)
counts_norm = counts(dds, normalized=TRUE)
write.table(counts_norm, "counts_normalized.tsv", sep = "\t")

### run degs analysis
dds = DESeq(dds)
res = results(dds)

### use contrast to make specific comparisons between groups of samples. 
### set p to threshold that will be used for deg significance
### DEseq2 will use this to perform p value correction

p = 0.05
genes_all = read_tsv("ncbi_genelist.tsv")

for (i in (c(17, 20, 23,"d0","d4","d8"))){
  degs = as.data.frame(results(dds, 
            contrast = c("group", paste0("m",i), paste0("f",i)), alpha = p))
  degs = degs[!is.na(degs$padj),]
    ### add chromosome and gene name from NCBI table
    symbols = as.list(row.names(degs))
    genename = character()
    chromosome = character()
    id = character()
    for (j in symbols){
        genename = c(genename, genes_all[match(j,genes_all$Symbol),"Name"])
        chromosome = c(chromosome, 
                       genes_all[match(j,genes_all$Symbol),"Chromosome"])
        id = c(id, genes_all[match(j,genes_all$Symbol),"Gene ID"])
      }
    genename = as.matrix(genename, ncol=1)
    chromosome = as.matrix(chromosome, ncol=1)
    id = as.matrix(id, ncol=1)
    degs = cbind(degs,genename,chromosome,id)
  assign(paste0("pvals_", i), degs)
  degs = degs[degs$padj <= p,]
  assign(paste0("degs_", i),degs)
  rm(degs, genename, chromosome, symbols, i, j, id)
}

# variable cleanup
rm(p)

### make tables of degs excluding chromosome 19 and Y
degs_17_auto = degs_17[!(degs_17$chromosome %in% c(19, "Y")),]
degs_20_auto = degs_20[!(degs_20$chromosome %in% c(19, "Y")),]
degs_23_auto = degs_23[!(degs_23$chromosome %in% c(19, "Y")),]
degs_d0_auto = degs_d0[!(degs_d0$chromosome %in% c(19, "Y")),]
degs_d4_auto = degs_d4[!(degs_d4$chromosome %in% c(19, "Y")),]
degs_d8_auto = degs_d8[!(degs_d8$chromosome %in% c(19, "Y")),]

### export all pvalues to tsv files
sets = ls(pattern = "pvals_")
for (i in sets){
  write.table(as.matrix(get(i)), 
              file = paste0("DEseq_",i,".tsv"), quote = F, sep ="\t")
}

### export list of degs to tsv files
sets = ls(pattern = "deg")
for (i in sets){
  write.table(as.matrix(get(i)), 
              file = paste0(i,".tsv"), quote = F, sep ="\t")
}

rm(i, sets)
  

### make list of autosomal degs

### make a list of all of the gene names from the deg tables
degs = c(row.names(degs_17_auto), rownames(degs_20_auto), 
         rownames(degs_23_auto), rownames(degs_d0_auto), 
         rownames(degs_d4_auto), rownames(degs_d8_auto))
degs = unique(degs)

### subset the full table of all genes to only DEGS
degs = genes_all[genes_all$Symbol %in% degs,]
degs = cbind("Symbol" = degs$Symbol, 
             "Name" = degs$Name, 
             "Chromosome" = degs$Chromosome, 
             "Gene ID" = degs$`Gene ID`)
degs = as.data.frame(degs)
degs = arrange(degs, tolower(Symbol))

### add p value column for each stage
degs = mutate(degs, s17 = pvals_17[degs$Symbol,"padj"])
degs = mutate(degs, s20 = pvals_20[degs$Symbol,"padj"])
degs = mutate(degs, s23 = pvals_23[degs$Symbol,"padj"])
degs = mutate(degs, d0 = pvals_d0[degs$Symbol,"padj"])
degs = mutate(degs, d4 = pvals_d4[degs$Symbol,"padj"])
degs = mutate(degs, d8 = pvals_d8[degs$Symbol,"padj"])

### count how many different stages a gene is a deg in
degs = mutate(degs, stages = rowSums((degs[,5:10] < 0.05)))

write_tsv(degs, file = "DEseq_degs_all_stages.tsv",)

rm(list=ls())

