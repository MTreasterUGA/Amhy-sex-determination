if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(tidyverse)

### generate coldata from counts file
{
  # read in raw counts file, set input to name of counts file in this director. Set output to desired name of coldata file for deseq2
  
  input = "counts_raw_summed.txt"
  output = "counts_raw_summed_coldata.txt"
  counts_raw = read.table(input)
  
  # remove rows of unaligned counts
  counts_raw = counts_raw[!(row.names(counts_raw) %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual")),]
  
  samples = colnames(counts_raw)
  coldata = data.frame(counts_raw = samples)
  
  
  ############################sex and stage as collapsed factor
  # pull group from first three characters of sample name
  group = vector()
  for (i in 1:nrow(coldata)){
    group = c(group, substr(coldata[i,1],1,3))
  }
  
  coldata = cbind(coldata,group)
  coldata = column_to_rownames(coldata,var="counts_raw")
  
  write.table(coldata,file=output,quote=F,row.names=T)
  
  # variable cleanup
  rm(group, i, input, output, samples)

  # check that names and order from counts and coldata table match
}

### to manually read in coldata and/or counts tables
counts_raw = read.table("counts_raw_summed.txt")
coldata = read.table("counts_raw_summed_coldata.txt",header=T,row.names=1)

### DEG analysis
{
  if(all(rownames(coldata) != colnames(counts_raw))) {
    print("samples in counts and coldata tables do not match")
    stop()
    print("this text should never appear")
  }
  
  # convert variable from characters to factors, make deseq data set object (dds)
  
  #coldata$sex = factor(coldata$sex)
  #coldata$stage = factor(coldata$stage)
  coldata$group = factor(coldata$group)
  
  dds = DESeqDataSetFromMatrix(countData = counts_raw, 
                               colData = coldata, 
                               #design = ~ sex + stage + sex:stage) # if sex and stage as separate factors
                               design = ~ group) # for sex and stage collapsed into group
  
  #remove all genes with fewer than 10 reads across all samples
  dds = dds[rowSums(counts(dds)) >= 10,]
  
  #generate normalized count table
  dds = estimateSizeFactors(dds)
  counts_norm = counts(dds, normalized=TRUE)
  
  write.table(counts_norm, "counts_normalized.txt")
  write.table(counts_norm, "counts_normalized.tsv", sep = "\t")
  
  #run degs analysis
  
  dds = DESeq(dds)
  res = results(dds)
  
  # use contrast to make specific comparisons between groups of samples. set p to threshold that will be used for deg significance. DEseq2 will use this to perform p value correction
  
  p = 0.05
  genes_all = read_tsv("ncbi_genelist.tsv")
  
  for (i in (c(17, 20, 23,"d0","d4","d8"))){
    degs = as.data.frame(results(dds, contrast = c("group", paste0("m",i), paste0("f",i)), alpha = p))
    degs = degs[!is.na(degs$padj),]
      ### add chromosome and gene name from NCBI table
      symbols = as.list(row.names(degs))
      genename = character()
      chromosome = character()
      id = character()
      for (j in symbols){
          genename = c(genename, genes_all[match(j,genes_all$Symbol),"Name"])
          chromosome = c(chromosome, genes_all[match(j,genes_all$Symbol),"Chromosome"])
          id = c(id, genes_all[match(j,genes_all$Symbol),"Gene ID"])
        }
      genename = as.matrix(genename, ncol=1)
      chromosome = as.matrix(chromosome, ncol=1)
      id = as.matrix(id, ncol=1)
      degs = cbind(degs,genename,chromosome,id)
    assign(paste0("pvals_", i), degs)
    degs = degs[degs$padj <= p,]
    assign(paste0("degs_", i),degs)
    rm(degs, genename, chromosome, symbols, i, j)
  }
  
  # variable cleanup
  rm(p, id)
  
  # make tables of degs excluding chromosome 19 and Y
  degs_17_auto = degs_17[!(degs_17$chromosome %in% c(19, "Y")),]
  degs_20_auto = degs_20[!(degs_20$chromosome %in% c(19, "Y")),]
  degs_23_auto = degs_23[!(degs_23$chromosome %in% c(19, "Y")),]
  degs_d0_auto = degs_d0[!(degs_d0$chromosome %in% c(19, "Y")),]
  degs_d4_auto = degs_d4[!(degs_d4$chromosome %in% c(19, "Y")),]
  degs_d8_auto = degs_d8[!(degs_d8$chromosome %in% c(19, "Y")),]
  
  ### export all pvalues to tsv files
  sets = ls(pattern = "pvals_")
  for (i in sets){
    write.table(as.matrix(get(i)), file = paste0("DEseq_",i,".tsv"), quote = F, sep ="\t")
  }
  
  ### export list of degs to tsv files
  sets = ls(pattern = "deg")
  for (i in sets){
    write.table(as.matrix(get(i)), file = paste0(i,".tsv"), quote = F, sep ="\t")
  }
  
  rm(i, sets)
}
  

### make list of autosomal degs
{
  
  # make a list of all of the gene names from the deg tables, discarding duplicates
  degs = c(row.names(degs_17_auto), rownames(degs_20_auto), rownames(degs_23_auto), rownames(degs_d0_auto), rownames(degs_d4_auto), rownames(degs_d8_auto))
  degs = unique(degs)
  
  # subset the full table of genes to the symbols, name, chromosome, and id of genes found in the list of autosomal degs
  degs = genes_all[genes_all$Symbol %in% degs,]
  degs = cbind(degs$Symbol, degs$Name, degs$Chromosome, degs$`Gene ID`)
  colnames(degs) = c("symbol", "name", "chromosome", "id")
  degs = column_to_rownames(as.data.frame(degs), var = "symbol")
  degs = degs[sort(rownames(degs)),]
  
  # make new columns in the table of degs with p values from the deg table for each stage, NA for not significant genes
  degs$s17 = pvals_17[rownames(degs),"padj"]
  degs$s20 = pvals_20[rownames(degs),"padj"]
  degs$s23 = pvals_23[rownames(degs),"padj"]
  degs$d0 = pvals_d0[rownames(degs),"padj"]
  degs$d4 = pvals_d4[rownames(degs),"padj"]
  degs$d8 = pvals_d8[rownames(degs),"padj"]
  
  # count how many different stages a gene is a deg in
  degs$stages = rowSums((degs[,4:9] < 0.05))
  
  write_tsv(rownames_to_column(degs, "symbol"), file = "DEseq_degs_all_stages.tsv",)
}

