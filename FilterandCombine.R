library(tidyverse)

###### Merge multiple count tables for different samples into one large table
{
  # read all files in raw_counts folder
  files = list.files("raw_counts/")
  
  # combine all files into single table if they have the same row names
  for (i in 1:length(files)){
    # on first iteration, make merged table from first counts file
    if (i == 1){
      merged = read.table(paste0("raw_counts/",files[i]),row.names = 1)
      colnames(merged) = paste0(files[i])
    }
    # on subsequent iterations, check if rownames (genes) are the same for the next count file
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
  # sort table by rownames alphabetically
  merged = rownames_to_column(merged)
  merged = merged[order(merged$rowname),]
  row.names(merged) = NULL
  merged = column_to_rownames(merged)
  
  # write merged table to new file in this directory
  write.table(merged, file="counts_raw.txt",quote =F)
  
  # variable cleanup
  rm(i, data, merged, files)
}

##### Combines counts from paired and unpaired read files for each sample into one column

{
  # change input to file in this directory, set output to desired name for exported file
  input="counts_raw.txt"
  output="counts_raw_summed.txt"
  
  # import counts table and remove rows where all values are zero
  
  data=read.table(input,header=T,row.names=1)
  data=data[rowSums(data)>0,]
  
  ### combine paired and unpaired read counts
  
  # check if adjacent columns are the same sample by comparing string before the fist underscore. if all pairs of columns do not match, stop running and give an error
  
  samplesin = colnames(data)
  for (i in 1:(ncol(data)/2)){
    if
    (substr(samplesin[2*i-1],1,regexpr("_",samplesin[2*i-1])-1) != substr(samplesin[2*i],1,regexpr("_",samplesin[2*i])-1)){
    print("error, some adjacent columns not from same sample")  
    stop()
    print("this text should never appear")
    }
    else {print(paste0("sample ",i, " matched"))}
  }
  
  # read in existing sample names for paired and unpaired sets. Generate one new name for each combined set from the characters before the first underscore. 
  samplesout = vector()
  for (i in 1:(ncol(data)/2)){
    n = samplesin[2*i-1]
    samplesout = c(samplesout,c(substr(n,1,regexpr("_",n)-1)))
  }
  
  # sum each pair of adjacent columns, i.e., sum columns 1 and 2, 3 and 4, 5 and 6, etc
  counts_raw = matrix(nrow = nrow(data))
  for (i in 1:(ncol(data)/2)){
    counts_raw = cbind(counts_raw,data[2*i-1]+data[2*i])
  }
  counts_raw = counts_raw[,-1]
  
  
  # rename columns to new sample names generated in samplesout and write file to output file
  colnames(counts_raw) = samplesout
  write.table(counts_raw,file=output,quote=F)
  
  # variable cleanup
  rm(i, data, n, counts_raw, input, output, samplesin, samplesout)
}

