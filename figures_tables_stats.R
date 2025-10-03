# color pallet    hex         rgb
# black           '#000000'   0,    0,    0
# white           '#ffffff'   255,  255,  255
# light gray      '#aaaaaa'   170,  170,  170
# dark gray       '#555555'   85,   85,   85
# orange          '#e69f00'   230,  159,  0
# sky blue        '#56b4e9'   86,   180,  233
# bluish green    '#009e73'   0,    158,  115
# yellow          '#f0e442'   240,  228,  66
# blue            '#0072b2'   0,    114,  178
# vermillion      '#d55e00'   213,  94,   0
# reddish purple  '#cc79a7'   204,  121,  167

c('#000000', '#ffffff', '#aaaaaa', '#555555', '#e69f00', 
  '#56b4e9', '#009e73', '#f0e442', '#0072b2', '#d55e00', '#cc79a7')

{
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(emmeans)
  library(betareg)
  library(ape)
}

################################################################################
##### plots of normalized counts ###############################################
################################################################################

### read in full normalized counts table
counts_norm = read.table("counts_normalized.txt",header = T, row.names = 1)

### pull sample names from counts table
samples = colnames(counts_norm)

### sex from first character of sample name
sex = vector()
for (i in 1:length(samples)){
  sex = c(sex, substr(samples[i],1,1))
}
subs = c("m" = "male",
         "f" = "female")
for (i in names(subs)){
  sex = gsub(i, subs[i], sex)
}

### stage from characters 2-3 of sample name
stage = vector()
for (i in 1:length(samples)){
  stage = c(stage, substr(samples[i],2,3))
}
subs = c("17" = "stage 17",
         "20" = "stage 20",
         "23" = "stage 23",
         "d0" = "0 dph",
         "d4" = "4 dph",
         "d8" = "8 dph")
for (i in names(subs)){
  stage = gsub(i, subs[i], stage)
}

### combine samples, sex, stage table. Make sex and stage readable with gsub
samples = as.data.frame(cbind(samples, sex, stage))
samples$stage = factor(samples$stage, levels = (c("stage 17",
                                                  "stage 20",
                                                  "stage 23",
                                                  "0 dph",
                                                  "4 dph",
                                                  "8 dph")))

### cleanup variables
rm (sex, stage, subs, i)

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
  ggsave(paste0("figures/counts_",gois[i,1],".pdf"), plot,
         width =1.75, height =1.5 , units = "in", dpi = 600)
  rm(data, counts)
}

### read in p value tables
for (i in c("pvals_17",
            "pvals_20",
            "pvals_23",
            "pvals_d0",
            "pvals_d4",
            "pvals_d8")){
  assign(i, read.table(paste0("DEseq_", i, ".tsv"),header = T,
                       row.names = 1, sep = "\t", quote = "\""))
}

### pull p values for all genes from table
pvalues = matrix(nrow = 0, ncol = 8)
colnames(pvalues) = c("gene",
                      "id",
                      "stage 17",
                      "stage 20",
                      "stage 23",
                      "0 dph",
                      "4 dph",
                      "8 dph")
for (i in 1:nrow(gois)){
  pvalues = rbind(pvalues, c(gois[i,1],
                             gois[i,2],
                             pvals_17[gois[i,2],]$padj,
                             pvals_20[gois[i,2],]$padj,
                             pvals_23[gois[i,2],]$padj,
                             pvals_d0[gois[i,2],]$padj,
                             pvals_d4[gois[i,2],]$padj,
                             pvals_d8[gois[i,2],]$padj))
}

write.table(pvalues, "figures/counts_pvalues.tsv",
            quote = F, sep = "\t", row.names = F)
rm(i, pvalues, gois)

##### plot counts for males only amha vs amhy #####

### pull counts for amha and amhy. name column accordingly
counts_amha = cbind(samples, counts = as.numeric(counts_norm["amh",]))
counts_amha$gene = "amha"

counts_amhy = cbind(samples, counts = as.numeric(counts_norm["LOC120812167",]))
counts_amhy$gene = "amhy"

### scale amhy values for difference in transcript length
counts_amhy_la = mutate(counts_amhy, counts = counts_amhy$counts * (3895/2496))

counts_amh = rbind(counts_amha, counts_amhy_la)

plot = 
  ggplot(counts_amh[counts_amh$sex == "male",]) +
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
ggsave("figures/counts_male_amh.pdf", plot,
       width =1.75, height =1.5 , units = "in", dpi = 600)

rm(counts_amha, counts_amhy, counts_amhy_la, counts_amh)

##### compare amh combined counts #####

counts_amh = mutate(samples, 
                  amha = as.numeric(counts_norm["amh",]),
                  amhy = as.numeric(counts_norm["LOC120812167",])) %>%
            mutate(amhy_la = .$amhy * (3895/2496)) %>%
            mutate(amh_total = .$amha + .$amhy_la)


### reorganize data for human readibility

total_amh = cbind(rbind(stage = levels(samples$stage), sex = "female"),
                  rbind(stage = levels(samples$stage), sex = "male"))

amha = matrix(as.numeric(counts_amh[,"amha"]), ncol=12, nrow=3, byrow = FALSE)
rownames(amha) = c("amha1", "amha2", "amha3")
amha = rbind(amha, amha_mean = NA,amha_stdev = NA)
for (i in 1:ncol(amha)) {
  amha["amha_mean", i] = mean(amha[1:3, i])
  amha["amha_stdev", i] = sd(amha[1:3, i])
}

amhy = matrix(as.numeric(counts_amh[,"amhy"]), ncol=12, nrow=3, byrow = FALSE)
rownames(amhy) = c("amhy1", "amhy2", "amhy3")
amhy = rbind(amhy, amhy_mean = NA,amhy_stdev = NA)
for (i in 1:ncol(amhy)) {
  amhy["amhy_mean", i] = mean(amhy[1:3, i])
  amhy["amhy_stdev", i] = sd(amhy[1:3, i])
}

amhy_la = matrix(as.numeric(counts_amh[,"amhy_la"]), ncol=12, nrow=3, byrow = FALSE)
rownames(amhy_la) = c("amhy_la1", "amhy_la2", "amhy_la3")
amhy_la = rbind(amhy_la, amhy_la_mean = NA,amhy_la_stdev = NA)
for (i in 1:ncol(amhy_la)) {
  amhy_la["amhy_la_mean", i] = mean(amhy_la[1:3, i])
  amhy_la["amhy_la_stdev", i] = sd(amhy_la[1:3, i])
}

amh_total = matrix(as.numeric(counts_amh[,"amh_total"]), ncol=12, nrow=3, byrow = FALSE)
rownames(amh_total) = c("amh_total1", "amh_total2", "amh_total3")
amh_total = rbind(amh_total, amh_total_mean = NA,amh_total_stdev = NA)
for (i in 1:ncol(amh_total)) {
  amh_total["amh_total_mean", i] = mean(amh_total[1:3, i])
  amh_total["amh_total_stdev", i] = sd(amh_total[1:3, i])
}

total_amh = rbind(total_amh, amha, amhy, amhy_la, amh_total)


### t tests for combined amh levels in females < males 

total_amh = rbind(total_amh, p_amh_total = c( "-","-","-","-","-","-",
        t.test(amh_total ~ sex, counts_amh[counts_amh$stage == "stage 17",], 
               alternative = "less")$p.value,
        t.test(amh_total ~ sex, counts_amh[counts_amh$stage == "stage 20",], 
               alternative = "less")$p.value,
        t.test(amh_total ~ sex, counts_amh[counts_amh$stage == "stage 23",], 
               alternative = "less")$p.value,
        t.test(amh_total ~ sex, counts_amh[counts_amh$stage == "0 dph",], 
               alternative = "less")$p.value,
        t.test(amh_total ~ sex, counts_amh[counts_amh$stage == "4 dph",], 
               alternative = "less")$p.value,
        t.test(amh_total ~ sex, counts_amh[counts_amh$stage == "8 dph",], 
               alternative = "less")$p.value
))

### t test for amhy != 1/2 amha

male.amha = counts_amh[,c(1, 2, 3, 4)]
male.amha = male.amha[male.amha$sex == "male",]
male.amha$amha = male.amha$amha/2
male.amha = mutate(male.amha, gene = "amha")
colnames(male.amha)[4] = "counts"

male.amhy = counts_amh[,c(1, 2, 3, 6)]
male.amhy = male.amhy[male.amhy$sex == "male",]
male.amhy = mutate(male.amhy, gene = "amhy")
colnames(male.amhy)[4] = "counts"

male.amh = rbind(male.amha, male.amhy)

total_amh = rbind(total_amh, p_amhy_0.5amh = c( "-","-","-","-","-","-",
t.test(counts ~ gene, male.amh[male.amh$stage == "stage 17",],
       alternative = "two.sided")$p.value,
t.test(counts ~ gene, male.amh[male.amh$stage == "stage 20",],
       alternative = "two.sided")$p.value,
t.test(counts ~ gene, male.amh[male.amh$stage == "stage 23",],
       alternative = "two.sided")$p.value,
t.test(counts ~ gene, male.amh[male.amh$stage == "0 dph",],
       alternative = "two.sided")$p.value,
t.test(counts ~ gene, male.amh[male.amh$stage == "4 dph",],
       alternative = "two.sided")$p.value,
t.test(counts ~ gene, male.amh[male.amh$stage == "8 dph",],
       alternative = "two.sided")$p.value
))

write.table(total_amh, "amh_all_counts.tsv", sep = "\t", row.names = T, col.names = F)

### anova for difference in amhy expression across stages

res.aov = aov(counts ~ stage, data = male.amhy)
summary(res.aov)
TukeyHSD(res.aov)

rm(amh_total, amha, amhy, amhy_la, counts_amh, i, male.amh, male.amha,
   male.amhy, plot, res.aov, total_amh)

################################################################################
##### Volcano Plots ############################################################
################################################################################

 ### lfc is female/male, positive = higher in female, negative = higher in male

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
  ggsave(paste0("figures/volcano_",sets[i],".pdf"), plot,
         width = 3.6, height = 2.8 , units = "in", dpi = 600)
}

rm(i, sets, data, plot)

################################################################################
##### PCA Plots ################################################################
################################################################################

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

rm(list = ls())

################################################################################
##### stats for crosses ########################################################
################################################################################

###### hatch rates in XX females and XY female siblings

hatch_rates = read.csv("Hatch Rates F2.csv")
hatch_total = hatch_rates[hatch_rates$genotype %in% c("Total XX", "Total XY"),]
hatch_rates = hatch_rates[hatch_rates$genotype %in% c("XX", "XY"),]
hatch_rates$genotype = factor(hatch_rates$genotype, levels = c("XX","XY"))

hatch_rates = hatch_rates |>
  group_by(genotype) |>
  arrange(eggs, .by_group = TRUE) |>
  mutate(offset_x = as.numeric(genotype) + (row_number() - mean(row_number())) * 0.05) |>
  ungroup()
plot = ggplot(hatch_rates) +
  aes(x = offset_x, y = eggs, fill = genotype) +
  geom_point(size = 2, shape = 21, show.legend = F) +
  labs(title = "clutch size", y = NULL, x = "genotype") +
  theme_light() +
  theme(text = element_text(size = 8, family = "sans"),
        axis.title.y = element_text(margin = margin(r=5)),
        plot.title = element_text(size = 10, margin = margin(b=5)),
        panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = c('#e69f00', '#d55e00')) + 
  scale_y_continuous(expand = c(0, NA), limits = c(0, 150)) +
  scale_x_continuous(breaks = 1:2, labels = levels(hatch_rates$genotype), limits = c(0.60, 2.4))
ggsave("figures/rate_clutch.pdf", plot,
       width = 1.5, height = 1.6 , units = "in", dpi = 600)

hatch_rates = hatch_rates |>
  group_by(genotype) |>
  arrange(fertilized.total, .by_group = TRUE) |>
  mutate(offset_x = as.numeric(genotype) + (row_number() - mean(row_number())) * 0.05) |>
  ungroup()
plot = ggplot(hatch_rates) +
  aes(x = offset_x, y = fertilized.total, fill = genotype) +
  geom_point(size = 2, shape = 21, show.legend = F) +
  labs(title = "viabilty rate", y = NULL, x = "genotype") +
  theme_light() +
  theme(text = element_text(size = 8, family = "sans"),
        axis.title.y = element_text(margin = margin(r=5)),
        plot.title = element_text(size = 10, margin = margin(b=5)),
        panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = c('#e69f00', '#d55e00')) + 
  #scale_y_continuous(expand = c(0, 1)) +
  scale_x_continuous(breaks = 1:2, labels = levels(hatch_rates$genotype), limits = c(0.60, 2.4))
ggsave("figures/rate_fert.pdf", plot,
       width = 1.5, height = 1.6 , units = "in", dpi = 600)

### test for difference in clutch size
compare_means(formula = eggs ~ genotype, data = hatch_rates, method = "wilcox.test")

### test for difference in fertilization rate
model = betareg(fertilized.total ~ genotype, data = hatch_rates)
summary(model)

rm(plot, hatch_rates, hatch_total)

##### XY x XY cross genotypes

### 36 YY embryos out of 148 total
chisq.test(x = c(36, 112), p = c(1/4, 3/4))

### 92 total, 36 XX, 27 XYp, 29 XYm
chisq.test(x = c(36, 27, 29), p = c(1/3, 1/3, 1/3))

################################################################################
##### color quantification #####################################################
################################################################################

color = read.table("color_data.tsv", header = T, sep = "\t")
color$genotype = factor(color$genotype, levels = (c("XXWT", "XXTG", "XYKO", "XYWT")))

color = color |>
  group_by(genotype) |>
  arrange(prop.dark, .by_group = TRUE) |>
  mutate(offset_x = as.numeric(genotype) + (row_number() - mean(row_number())) * 0.05) |>
  ungroup()
plot = ggplot(color) +
  aes(x = offset_x, y = prop.dark, fill = genotype) +
  geom_point(size = 2, shape = 21, show.legend = F) +
  labs(title = "prop.dark", y = NULL, x = "genotype") +
  theme_light() +
  theme(text = element_text(size = 8, family = "sans"),
        axis.title.y = element_text(margin = margin(r=5)),
        plot.title = element_text(size = 10, margin = margin(b=5)),
        panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = c('#e69f00', '#56b4e9', '#d55e00', '#0072b2')) +
  scale_x_continuous(breaks = 1:4, labels = levels(color$genotype), limits = c(0.7, 4.3)) +
  ylim(0,1)
ggsave("figures/prop_dark.pdf", plot, width = 2, height = 2 , units = "in", dpi = 600)

color = color |>
  group_by(genotype) |>
  arrange(prop.blue, .by_group = TRUE) |>
  mutate(offset_x = as.numeric(genotype) + (row_number() - mean(row_number())) * 0.05) |>
  ungroup()
plot = ggplot(color) +
  aes(x = offset_x, y = prop.blue, fill = genotype, shape = genotype) +
  geom_point(size = 2, shape = 21, show.legend = F) +
  labs(title = "prop.blue", y = NULL, x = "genotype") +
  theme_light() +
  theme(text = element_text(size = 8, family = "sans"),
        axis.title.y = element_text(margin = margin(r=5)),
        plot.title = element_text(size = 10, margin = margin(b=5)),
        panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = c('#e69f00', '#56b4e9', '#d55e00', '#0072b2')) +
  scale_x_continuous(breaks = 1:4, labels = levels(color$genotype), limits = c(0.7, 4.3)) +
  ylim(0,1)
ggsave("figures/prop_blue.pdf", plot, width = 2, height = 2 , units = "in", dpi = 600)

##### beta regression for pairwise comparisons

model = betareg(prop.blue ~ genotype | genotype, data = color)
#summary(model)
emmeans(model, pairwise ~ genotype, type = "response")

model = betareg(prop.dark ~ genotype | genotype, data = color)
#summary(model)
emmeans(model, pairwise ~ genotype, type = "response")

################################################################################
##### dnds calculation #########################################################
################################################################################

ape_dnds = as.data.frame(list.files("alignments/"))
colnames(ape_dnds) = "file"

for (i in 1:nrow(ape_dnds)) {
  alignment = read.dna(paste0("alignments/", ape_dnds[i, 1]), format = "fasta")
  alignment = del.colgapsonly(alignment, threshold = 0.1)
  ape_dnds$dn.ds[i] = as.numeric(dnds(alignment))
}

write.table(ape_dnds, file = "dnds.tsv", sep = "\t", row.names = F)
