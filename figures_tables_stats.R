# color pallet    hex         rgb
# black           '#000000'   0, 0, 0
# white           '#ffffff'   255, 255, 255
# light gray      '#aaaaaa'   170, 170, 170
# dark gray       '#555555'   85, 85, 85
# orange          '#e69f00'   230, 159, 0
# sky blue        '#56b4e9'   86, 180, 233
# bluish green    '#009e73'   0, 158, 115
# yellow          '#f0e442'   240, 228, 66
# blue            '#0072b2'   0, 114, 178
# vermillion      '#d55e00'   213, 94, 0
# reddish purple  '#cc79a7'   204, 121, 167

c('#000000', '#ffffff', '#aaaaaa', '#555555', '#e69f00', '#56b4e9', '#009e73', '#f0e442', '#0072b2', '#d55e00', '#cc79a7')



{
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(ggpubr)
}

########################### plots of normalized counts #########################

#### read in full normalized counts table
counts_norm = read.table("counts_normalized.txt",header = T, row.names = 1)

{
### pull sample names from counts table and sex from first character of sample name
  samples = colnames(counts_norm)
  sex = vector()
  for (i in 1:length(samples)){
    sex = c(sex, substr(samples[i],1,1))
  }
  
  ### pull stage from characters 2-3 of sample name
  stage = vector()
  for (i in 1:length(samples)){
    stage = c(stage, substr(samples[i],2,3))
  }
  
  ### combine samples sex and stage into one table. Make sex and stage readable with gsub
  samples = as.data.frame(cbind(samples, sex, stage))
  samples$sex = samples$sex %>% gsub("m", "male", .) %>% gsub("f","female", .)
  samples$stage = samples$stage %>% gsub("17", "stage 17", .) %>% gsub("20", "stage 20", .) %>% gsub("23", "stage 23", .) %>% gsub("d0", "0 dah", .) %>% gsub("d4", "4 dah", .) %>% gsub("d8", "8 dah", .)
  samples$stage = factor(samples$stage, levels = c("stage 17", "stage 20", "stage 23", "0 dah", "4 dah", "8 dah"))
  
  # variable cleanup
  rm(i, sex, stage)
  
  ### read in list of genes of interest from table in working directory
  gois = read.table("genes_of_interest.tsv", header = T, sep = "\t")
  
  ### make pdf of counts for each gene in genelist
  for (i in 1:nrow(gois)){
    counts = as.numeric(counts_norm[gois[i,2],])
    data = cbind(samples, counts)
    data = data[order(data$counts),]
    plot = ggplot(data) +
      aes(x = stage, y = counts, shape = sex, fill = sex) +
      geom_point(size = 2, position = position_dodge(width = 0.2), show.legend = F) +
      theme_light() +
      theme(text = element_text(size = 8, family = "sans"),
            axis.title.y = element_text(margin = margin(r=15)),
            plot.title = element_text(face = "italic", margin = margin(r=10)),
            ) +
      scale_shape_manual(values = c(21,22)) +
      scale_fill_manual(values = c('#e69f00', '#0072b2')) +
      scale_x_discrete(labels = c("s17", "s20", "s23", "d0", "d4", "d8")) +
      labs(x = NULL, y = NULL, title = gois[i,1]) +
      scale_y_continuous(breaks = pretty_breaks(n = 4))
    # save plot in figures directory                  
    ggsave(paste0("figures/counts_",gois[i,1],".pdf"), plot, width =1.75, height =1.5 , units = "in", dpi = 600)
  }
  
  ### pull p values for all genes from table
  
  pvalues = matrix(nrow = 0, ncol = 8)
  for (i in 1:nrow(gois)){
    pvalues = rbind(pvalues, c(gois[i,1], gois[i,2], pvals_17[gois[i,2],]$padj, pvals_20[gois[i,2],]$padj, pvals_23[gois[i,2],]$padj, pvals_d0[gois[i,2],]$padj, pvals_d4[gois[i,2],]$padj, pvals_d8[gois[i,2],]$padj))
  }
  
  colnames(pvalues) = c("gene", "id", "stage 17", "stage 20", "stage 23", "0 dah", "4 dah", "8 dah")
  
  write.table(pvalues, "figures/counts_pvalues.tsv", quote = F, sep = "\t", row.names = F)
  rm(i, pvalues, gois)
}

## plot counts for males only amha vs amhy

#pull counts for amha and amhy. name column accordingly
amha_counts = cbind(samples, as.numeric(counts_norm["amh",]))
colnames(amha_counts)[4] = "counts"
amha_counts = amha_counts[19:36,]
amha_counts$gene = "amha"

amhy_counts = cbind(samples, as.numeric(counts_norm["LOC120812167",]))
colnames(amhy_counts)[4] = "counts"
amhy_counts = amhy_counts[19:36,]
amhy_counts$gene = "amhy"

amh_counts = rbind(amh_counts, amhy_counts)

plot = 
  ggplot(amh_counts) +
  aes(x = stage, y = counts, shape = gene, fill = gene) +
  geom_point(size = 2, position = position_dodge(width = 0.2), show.legend = F) +
  theme_light() +
  theme(text = element_text(size = 8, family = "sans"),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(face = "italic", margin = margin(r=10)),
  ) +
  scale_shape_manual(values = c(24,25)) +
  scale_fill_manual(values = c('#009e73', '#56b4e9')) +
  scale_x_discrete(labels = c("s17", "s20", "s23", "d0", "d4", "d8")) +
  labs(x = NULL, y = NULL, title = "XY amha and amhy") +
  scale_y_continuous(breaks = pretty_breaks(n = 3))
# save plot in figures directory                  
ggsave(paste0("figures/counts_male_amh.pdf"), plot, width =1.75, height =1.5 , units = "in", dpi = 600)


######################### amh combined counts ##############################

## make "samples" table from above section ^^^

amha = as.numeric(counts_norm["amh",])
amhy = as.numeric(counts_norm["LOC120812167",])
amhy_length_adj = amhy*(3895/2496)
total.amh = as.numeric(amha) + as.numeric(amhy_length_adj)
amha = cbind(samples, amha)
amhy = cbind (samples, amhy)
amhy_length_adj = cbind(samples,amhy_length_adj)
total.amh = cbind (samples, total.amh)

### t tests for combined amh levels in males vs females. manually add p value to final table
t.test(total.amh ~ sex, total.amh[total.amh$stage == "stage 17",], alternative = "two.sided")
t.test(total.amh ~ sex, total.amh[total.amh$stage == "stage 20",], alternative = "two.sided")
t.test(total.amh ~ sex, total.amh[total.amh$stage == "stage 23",], alternative = "two.sided")
t.test(total.amh ~ sex, total.amh[total.amh$stage == "0 dah",], alternative = "two.sided")
t.test(total.amh ~ sex, total.amh[total.amh$stage == "4 dah",], alternative = "two.sided")
t.test(total.amh ~ sex, total.amh[total.amh$stage == "8 dah",], alternative = "two.sided")

#### reorganize data. Put sample stage and sex in top two rows. Under each row, have a column with three samples and the counts of the gene at that stage and sex
temptable = rbind(c("stage 17","stage 20","stage 23","0 dah","4 dah","8 dah","stage 17","stage 20","stage 23","0 dah","4 dah","8 dah"),c("female","female","female","female","female","female","male","male","male","male","male","male"))

amha = rbind(temptable, matrix(as.numeric(amha[,4]), ncol=12, nrow=3))
amha = as.data.frame(cbind(c("stage", "sex", "count1", "count2", "count3"),amha))

amhy = rbind(temptable, matrix(as.numeric(amhy[,4]), ncol=12, nrow=3))
amhy = as.data.frame(cbind(c("stage", "sex", "count1", "count2", "count3"),amhy))

amhy_length_adj = rbind(temptable, matrix(as.numeric(amhy_length_adj[,4]), ncol=12, nrow=3))
amhy_length_adj = as.data.frame(cbind(c("stage", "sex", "count1", "count2", "count3"),amhy_length_adj))

total.amh = rbind(temptable, matrix(as.numeric(total.amh[,4]), ncol=12, nrow=3))
total.amh = as.data.frame(cbind(c("stage", "sex", "count1", "count2", "count3"),total.amh))

### means and SD to each table
amha[6,1] = "mean"
amha[7,1] = "sd"
for (i in 2:ncol(amha)) {
  amha[6, i] = mean(as.numeric(amha[3:5,i]))
  amha[7, i] = sd(as.numeric(amha[3:5,i]))
}
write.table(amha, "amha.tsv", sep = "\t", quote = F, row.names = F)

amhy[6,1] = "mean"
amhy[7,1] = "sd"
for (i in 2:ncol(amhy)) {
  amhy[6, i] = mean(as.numeric(amhy[3:5,i]))
  amhy[7, i] = sd(as.numeric(amhy[3:5,i]))
}
write.table(amhy, "amhy.tsv", sep = "\t", quote = F, row.names = F)

amhy_length_adj[6,1] = "mean"
amhy_length_adj[7,1] = "sd"
for (i in 2:ncol(amhy_length_adj)) {
  amhy_length_adj[6, i] = mean(as.numeric(amhy_length_adj[3:5,i]))
  amhy_length_adj[7, i] = sd(as.numeric(amhy_length_adj[3:5,i]))
}
write.table(amhy, "amhy_length_adj.tsv", sep = "\t", quote = F, row.names = F)

total.amh[6,1] = "mean"
total.amh[7,1] = "sd"
for (i in 2:ncol(total.amh)) {
  total.amh[6, i] = mean(as.numeric(total.amh[3:5,i]))
  total.amh[7, i] = sd(as.numeric(total.amh[3:5,i]))
}
write.table(total.amh, "totalamh.tsv", sep = "\t", quote = F, row.names = F)

### anova for difference in amhy expression across stages

male.amhy = cbind(samples, counts_norm["LOC120812167",])
male.amhy = male.amhy[male.amhy$sex == "male",]
colnames(male.amhy)[4] = "counts"

res.aov = aov(counts ~ stage, data = male.amhy)
summary(res.aov)
TukeyHSD(res.aov)

################################# Volcano Plots ################################
 ### log fold change is female/male, positive value is higher in female, negative is higher in male
{
  sets = ls(pattern = "pvals_")
  
  for (i in 1:length(sets)) {
    data = get(sets[i])
    # remove X and Y chromosome genes
    data = data[!(data$chromosome %in% c(19, "Y")),]
    # make column indicating if adjusted p value is significant or not
    data$significant = ifelse(data$padj < 0.05, ifelse(data$stat < 0, "male", "female"), "ns")
    plot = ggplot(data) + 
      aes(x = log2FoldChange, y = -log10(pvalue), color = significant) +
      geom_point(show.legend = F, size = 1) + 
      theme_light() +
      theme(text = element_text(size = 8, family = "sans"),
            axis.title.y = element_text(margin = margin(r=5)),
            plot.title = element_text(size = 10, margin = margin(b=5))
            ) +
      scale_color_manual(values = c("male" = '#0072b2', "female" = '#e69f00', "ns" = '#dcdcdc')) +
      labs(title = sets[i])
    # save plot in figures directory  
    ggsave(paste0("figures/volcano_",sets[i],".pdf"), plot, width = 3.6, height = 2.8 , units = "in", dpi = 600)
  }
  
  rm(i, sets, data, plot)
}
########################### PCA Plots ##########################################
{
  genes_all = read_tsv("ncbi_genelist.tsv")
  genes_y = subset(genes_all, genes_all$Chromosome == "Y")
  genes_x = subset(genes_all, genes_all$Chromosome == 19)
  genes_auto = subset(genes_all, !(genes_all$Chromosome %in% c(19, "Y")))
  
  counts_norm_y = subset(counts_norm, row.names(counts_norm) %in% genes_y$Symbol)
  counts_norm_x = subset(counts_norm, row.names(counts_norm) %in% genes_x$Symbol)
  counts_norm_auto = subset(counts_norm, row.names(counts_norm) %in% genes_auto$Symbol)
  
  sets = ls(pattern = "counts_norm")
  
  for (i in 1:length(sets)){
    
    # transpose table so genes are columns and samples are rows
    counts_transposed = t(get(sets[i]))
    
    # pca analysis
    pca_res = prcomp(counts_transposed, scale = T)
    pca_summ = summary(pca_res)
    
    # make table with PCA results of the top two components
    pca_data = data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])
    pca_data = cbind(pca_data, samples)
    
    # make pca plot
    plot = ggplot(pca_data) +
      aes(x = PC1, y = PC2, shape = stage, color = sex) +
      geom_point(size = 2, show.legend = F) +
      theme_light() +
      theme(text = element_text(size = 8, family = "sans"),
            axis.title.y = element_text(margin = margin(r=5)),
            plot.title = element_text(size = 10, margin = margin(b=5))
            ) +
      #geom_text_repel(color = '#000000', point.padding = 1, box.padding = .5, force = 2, min.segment.length = 0.1) +
      scale_color_manual(values = c('#e69f00', '#0072b2')) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) +
      #scale_shape_manual(values = c("\u25a0", "\u25b2", "\u25bc", "\u25c6", "\u25cf", "\u2726")) +
      guides(fill=guide_legend(override.aes=list(shape=21))) +
      labs(title = sets[i],
           x = paste0("PC1  (",pca_summ$importance[2,1],")"),
           y = paste0("PC2  (",pca_summ$importance[2,2],")")
           )
     
    
    # save pdf of pca plot
    ggsave(paste0("figures/pca_",sets[i],".pdf"), plot, width = 3.6, height = 3.6 , units = "in", dpi = 600)
  }
  
  rm(sets, i, counts_transposed, pca_res, pca_summ, plot)
}

######################## hatch rates ###########################################

{
  hatch_rates = read.csv("Hatch Rates F2.csv")
  hatch_total = hatch_rates[hatch_rates$genotype %in% c("Total XX", "Total XY"),]
  hatch_rates = hatch_rates[hatch_rates$genotype %in% c("XX", "XY"),]
  
  plot = ggplot(hatch_rates) +
    aes(x = genotype, y = eggs, fill = genotype) +
    geom_point(size = 3, shape = 21, position = position_dodge2(width = .1), show.legend = F) +
    labs(title = "clutch size", y = NULL) +
    theme_light() +
    theme(text = element_text(size = 8, family = "sans"),
          axis.title.y = element_text(margin = margin(r=5)),
          plot.title = element_text(size = 10, margin = margin(b=5))
          ) +
    scale_color_manual(values = c('#e69f00', '#0072b2')) + 
    scale_y_continuous(expand = c(0, NA), limits = c(0, 150))
  ggsave("figures/rate_clutch.pdf", plot, width = 1.5, height = 2 , units = "in", dpi = 600)
  
  plot = ggplot(hatch_rates) +
    aes(x = genotype, y = fertilized.total, fill = genotype) +
    geom_point(size = 3, shape = 21, position = position_dodge2(width = .1), show.legend = F) +
    labs(title = "fertilizaton rate", y = NULL) +
    theme_light() +
    theme(text = element_text(size = 8, family = "sans"),
          axis.title.y = element_text(margin = margin(r=5)),
          plot.title = element_text(size = 10,margin = margin(b=5))
    ) +
    scale_color_manual(values = c('#e69f00', '#0072b2'))
  ggsave("figures/rate_fert.pdf", plot, width = 1.5, height = 2 , units = "in", dpi = 600)
  
  plot = ggplot(hatch_rates) +
    aes(x = genotype, y = hatched.total, fill = genotype) +
    geom_point(size = 3, shape = 21, position = position_dodge2(width = .1), show.legend = F) +
    labs(title = "hatched:total", y = NULL) +
    theme_light() +
    theme(text = element_text(size = 8, family = "sans"),
          axis.title.y = element_text(margin = margin(r=5)),
          plot.title = element_text(size = 10, margin = margin(b=5))
    ) +
    scale_color_manual(values = c('#e69f00', '#0072b2'))
  ggsave("figures/rate_hatch.pdf", plot, width = 1.5, height = 2 , units = "in", dpi = 600)
  
  plot = ggplot(hatch_rates) +
    aes(x = genotype, y = hatched.fertilized, fill = genotype) +
    geom_point(size = 3, shape = 21, position = position_dodge2(width = .1), show.legend = F) +
    labs(title = "hatched:fertilized", y = NULL) +
    theme_light() +
    theme(text = element_text(size = 8, family = "sans"),
          axis.title.y = element_text(margin = margin(r=5)),
          plot.title = element_text(size = 10, margin = margin(b=5))
    ) +
    scale_color_manual(values = c('#e69f00', '#0072b2'))
  ggsave("figures/rate_hatchfert.pdf", plot, width = 1.5, height = 2 , units = "in", dpi = 600)
  
  rm(plot, hatch_rates, hatch_total)
}

compare_means(formula = eggs ~ genotype, data = hatch_rates, method = "wilcox.test")

# test for difference in fertilization rate
prop.test(x = hatch_total$fertilized, n = hatch_total$eggs, alternative = "two.sided")

# test for difference in hatched eggs out of total
prop.test(x = hatch_total$hatched, n = hatch_total$eggs, alternative = "two.sided")

# test for difference in hatched eggs out of fertilized
prop.test(x = hatch_total$hatched, n = hatch_total$fertilized, alternative = "two.sided")

##################### XY x XY genotypes ##########################

##### 36 YY embryos out of 148 total
chisq.test(x = c(36, 112), p = c(1/4, 3/4))

#### 92 total, 36 XX, 27 XYp, 29 XYm
chisq.test(x = c(36, 27, 29), p = c(1/3, 1/3, 1/3))

##################### dnds calculation ########################

library(ape)

ape_dnds = as.data.frame(list.files())
colnames(ape_dnds) = "file"

for (i in 1:length(list.files())) {
  alignment = read.dna(list.files()[i], format = "fasta")
  alignment = del.colgapsonly(alignment, threshold = 0.1)
  ape_dnds$dn.ds[i] = as.numeric(dnds(alignment))
}

library(seqinr)

seqinr_dnds = as.data.frame(list.files())
colnames(seqinr_dnds) = "file"

for (i in 1:length(list.files())) {
  alignment = read.alignment(list.files()[i], format = "fasta")
  kaks = kaks(alignment)
  seqinr_dnds$ka[i] = as.numeric(kaks[1])
  seqinr_dnds$ks[i] = as.numeric(kaks[2])
}
seqinr_dnds$ka.ks = seqinr_dnds$ka/seqinr_dnds$ks

alignment = read.alignment("alignments/gac_amha_vs_gac_amhy_cds.fasta", format = "fasta")
kaks(test)
list.files()






