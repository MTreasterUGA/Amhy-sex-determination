
library(tidyverse)

################################################################################
##### Merge counts tables from each sample into one table ######################
################################################################################

### read all files in raw_counts folder
files = list.files("raw_counts/")

### combine all files into single table if they have the same row names
for (i in 1:length(files)){
  # on first iteration, make merged table from first counts file
  if (i == 1){
    merged = read.table(paste0("raw_counts/",files[i]),row.names = 1)
    colnames(merged) = paste0(files[i])
  }
  # on each iteration, check if rownames (genes) are the same for the count file
  else{
    data = read.table(paste0("raw_counts/",files[i]),row.names = 1)
    colnames(data) = paste0(files[i])
    # if rownames are the same, add counts to merged table
    if (all(row.names(merged) == row.names(data))) {
      merged = cbind(merged,data)
    }
    # if rownames are different, print which file has mismatched rownames
    else{
      print(paste0("ERROR  ", files[i], " has mismatched row names"))
    }
  }
}

### sort table by rownames alphabetically
merged = rownames_to_column(merged)
merged = merged[order(merged$rowname),]
row.names(merged) = NULL
merged = column_to_rownames(merged)

### write merged table to new file in this directory
write.table(merged, file="counts_raw.tsv", quote=F, sep="\t" )

### variable cleanup
rm(list=ls())


################################################################################
##### Sum paired and unpaired read counts for each sample into one column ######
################################################################################

### import counts table and remove rows where all values are zero

data=read.table("counts_raw.tsv",header=T,row.names=1)
data=data[rowSums(data)>0,]

### combine paired and unpaired read counts

### check if adjacent 3 columns are the same sample

samplesin = colnames(data)
for (i in 1:(ncol(data)/3)){
  if
  (substr(samplesin[3*i-2],1,regexpr("_",samplesin[3*i-2])-1) != 
   substr(samplesin[3*i-1],1,regexpr("_",samplesin[3*i-1])-1) & 
   substr(samplesin[3*i-2],1,regexpr("_",samplesin[3*i-2])-1) != 
   substr(samplesin[3*i],1,regexpr("_",samplesin[3*i])-1)) {
  print("error, some adjacent columns not from same sample")  
  stop()
  print("this text should never appear")
  }
  else {print(paste0("sample ",i, " matched"))}
}

### generate sample name for each set from characters before first underscore 
samplesout = vector()
for (i in 1:(ncol(data)/3)){
  n = samplesin[3*i-1]
  samplesout = c(samplesout,c(substr(n,1,regexpr("_",n)-1)))
}

### sum every three adjacent columns, i.e. 1+2+3, 4+5+6, 7+8+9, etc
counts_raw = matrix(nrow = nrow(data))
for (i in 1:(ncol(data)/3)){
  counts_raw = cbind(counts_raw,data[3*i-2]+data[3*i-1]+data[3*i])
}
counts_raw = counts_raw[,-1]


### rename columns to new sample names write file to output file
colnames(counts_raw) = samplesout

### summarize aligned and unaligned reads per sample
unaligned = counts_raw[(1:5),]
aligned = colSums(counts_raw[-(1:5),])
total = colSums(counts_raw)

align_rate = rbind(unaligned, aligned = aligned, total = total)
write.table(align_rate, file = "alignment_rate.tsv", quote=F, sep = "\t")

### save counts for genes only
counts_raw = counts_raw[-(1:5),]
write.table(counts_raw,file="counts_raw_summed.tsv",quote=F,sep="\t")

### variable cleanup
rm(list=ls())

