require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(dplyr)
require(reshape2)
require(preprocessCore)


samples <- read.table("sample-list.csv", header = TRUE) 

#recoding of the groups for simplicity:
samples$treatment %>% 
  recode(control_control_condition_control_additive = "intact",
         control_carrier_control_additive = "carrier",
         stem_cells_control_condition_control_additive = "SC_native",
         stem_cells_graphenea_control_additive = "SC_grafenea",
         stem_cells_go_au_control_additive = "SC_GO_AU"
         ) -> samples$treatment

samples$group <- paste(samples$damage, samples$treatment, sep = "_")

samples %>% select(sample_id, file, damage, treatment, group) -> samples
grafenea_sc_samples <- scan("samples-grafens-stem-cells.txt", what=character(), sep="\n")

samples %>% filter(
  samples$file%in% grafenea_sc_samples) -> samples_grafenea

#the groups included in this analysis are: 
table(samples_grafenea$group)

#### FIRST ANALYSIS - TWO-WAY ANOVA

#factor 1 - damage
#factor 2 - treatment

#the intact group cannot be considered in this model 


fpkm <- data.frame(read.table('genes.fpkm_table.tsv',header = TRUE, sep = "\t"))
rnor.gene.list <- data.frame(read.delim("rnor.genes.tsv"))
gene.name <- rnor.gene.list$Gene.name[match(fpkm$tracking_id, rnor.gene.list$Gene.stable.ID)]
fpkm <- data.frame(gene.name, fpkm)
rownames(fpkm) <- fpkm$tracking_id
colnames(fpkm)[2] <- c("transcript.ID")
rm(rnor.gene.list)

fpkms.normalised <- data.matrix(fpkm[,c(-1,-2)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised)
fpkms.log <- log2(fpkms.normalised + 1)
rm(fpkms.normalised)

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA
fpkms.log <- data.frame(fpkms.log)
colnames(fpkms.log) <- colnames(fpkm[,c(3:117)])
rownames(fpkms.log) <- rownames(fpkm)

#ANOVA analysis 2-way:
results <- data.frame(fpkm$gene.name, fpkm$transcript.ID, row.names = rownames(fpkm))

two.way <- function(counts, damage, treatment) {
  if (is.na(counts[1])) {
    c(rep(NA, 3))
  } else {
    aov(as.numeric(counts) ~ as.factor(damage)*as.factor(treatment)) %>%
      summary() %>% unlist() -> res
  
    res[c("Pr(>F)1","Pr(>F)2","Pr(>F)3")]
  }
}


apply(fpkms.log[,match(samples_grafenea$file, colnames(fpkms.log))],
      1,
      two.way,
      damage=samples_grafenea$damage,
      treatment=samples_grafenea$treatment) %>% 
  t() %>%
  apply(.,
        2,
        p.adjust,
        method = 'fdr'
  ) %>%
  cbind(results, .) -> results

### perform post-hoc tests vs mia ###

stat_paired_t <- function(x) {
  if (is.na(x[1])) { c(NA, NA, NA, NA, NA)
  } else { 
    pairwise.t.test(
      as.numeric(x),
      samples_grafenea$group,
      p.adjust.method = 'none') %>% 
      unlist() -> res
    res[(c("p.value1", "p.value13", "p.value12", "p.value10", "p.value11"))] %>%
      p.adjust(method = 'bonferroni')
  }}

apply(fpkms.log[,match(samples_grafenea$file, colnames(fpkms.log))],
      1,
      stat_paired_t) %>% t() %>%
  cbind(results, .) -> results

colnames(results) <- c("gene_symbol",
                       "transcript_id",
                       "FDR_damage",
                       "FDR_treatment",
                       "FDR_int",
                       "t_test_intact_vs_mia",
                       "t_test_nacl_vs_mia",
                       "t_test_SC_vs_mia",
                       "t_test_SCGOAU_vs_mia",
                       "t_test_SCgrafenea_vs_mia")

#order: intact vs mia (1), nacl vs mia (13), miaSC vs mia (12), miaSC_GOAU vs mia (10), MIASC_grafenea vs mia (11)

#add a max dispersion factor:

get.max.dispersion <- function(x) {
  x %>% cbind(., samples_grafenea$group[order(samples_grafenea$damage, samples_grafenea$group)]) %>% 
    data.frame %>% rownames_to_column() %>%  
    set_colnames(c("file", "counts", "group")) %>%
    group_by(group) %>% 
    summarise(mean.to.sd = mean(as.numeric(counts))/sd(as.numeric(counts))) %>%
    select(mean.to.sd) %>% max() -> y
  return(y)
}

#compute selected fold changes between important groups:
fpkms.log[match(
  results$transcript_id,
  rownames(fpkms.log)),
  samples_grafenea$file[order(samples_grafenea$damage, samples_grafenea$group)]] %>%
  apply(., 1, get.max.dispersion) -> results$max.mean.to.sd

### get direction and fold changes for important comparisons ### (note - the intact group did not go into statistics)

compute_folds <- function(df, sample_info, sample_names_re, factors) {
  sample_info_grouped <-
    sample_info %>%
    group_by(across(all_of(factors)))
  groups <-
    sample_info_grouped %>%
    group_split()
  group_keys <-
    sample_info_grouped %>%
    group_keys() %>% 
    unite('name', everything(), remove = FALSE)
  
  for (tb_idx in 1:length(groups)) {
    tb <- groups[[tb_idx]]
    tb_key <- group_keys[tb_idx, 'name'] %>% unlist()
    df <-
      df %>%
      mutate(
        "{tb_key}_mean" := apply(
          select(., matches(sample_names_re)),
          1,
          function(x, names) {mean(x[names])},
          names = tb$id
        )
      )
  }
  df
}

samples_grafenea$id <- samples_grafenea$file

compute_folds(
  fpkms.log,
  samples_grafenea,
  'Cai',
  factors = c('group')
) -> group_means

### add information about direction and fold


results$fold_intact_vs_mia <- abs(group_means$intact_intact_mean-group_means$MIA_carrier_mean)
results$dir_intact_vs_mia <- ((group_means$intact_intact_mean-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "UP", "TRUE" = "DOWN")

results$fold_mia_vs_sc <- abs(group_means$MIA_SC_native_mean-group_means$MIA_carrier_mean)
results$dir_mia_vs_sc <- ((group_means$MIA_SC_native_mean-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "DOWN", "TRUE" = "UP")

results$fold_mia_vs_sc_go_au <- abs(group_means$MIA_SC_GO_AU_mean-group_means$MIA_carrier_mean)
results$dir_mia_vs_sc_go_au <- ((group_means$MIA_SC_GO_AU_mean-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "DOWN", "TRUE" = "UP")
  
results$fold_mia_vs_graf <- abs(group_means$MIA_SC_grafenea_mean-group_means$MIA_carrier_mean)
results$dir_mia_vs_graf <- ((group_means$MIA_SC_grafenea-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "DOWN", "TRUE" = "UP")


# bind results with counts
results <- cbind(results,fpkms.log)

###select three genes verified in PCR for export:
results %>% filter(results$gene_symbol %in% c('Comp', 'Mmp3', 'Mmp13')) -> selected
write_csv(selected, 'pcr_genes_results.csv')


###save top regulated genes:
results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_SC_vs_mia, t_test_SCGOAU_vs_mia, t_test_SCgrafenea_vs_mia), ~ . < 0.05)) %>%
  write_csv(., 'SC_regulated_genes_results.csv')

results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_SC_vs_mia, t_test_SCGOAU_vs_mia, t_test_SCgrafenea_vs_mia), ~ . < 0.05)) %>%
  select(gene_symbol) %>%
  write.table(row.names = FALSE, quote = FALSE)

### plot top results ###
reordering <- order(samples_grafenea$damage, samples_grafenea$treatment)


to_plot <- results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_SC_vs_mia, t_test_SCGOAU_vs_mia, t_test_SCgrafenea_vs_mia), ~ . < 0.05)) %>%
  select(samples_grafenea$file[reordering])


rownames(to_plot) <- results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_SC_vs_mia, t_test_SCGOAU_vs_mia, t_test_SCgrafenea_vs_mia), ~ . < 0.05)) %>%
  select(gene_symbol) %>% t()


#### actual plotting ####

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(as.character(samples_grafenea$group[reordering]))

col.labels <- c(rep("", 3),
                group.names[1], rep(" ", 6), 
                group.names[2], rep(" ", 6),
                group.names[3], rep(" ", 5),
                group.names[4], rep(" ", 5),
                group.names[5], rep(" ", 4),
                group.names[6], rep(" ", 4),
                group.names[7], rep(" ", 2),
                group.names[8], rep(" ", 3),
                group.names[9], rep(" ", 2))


cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


to_plot %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to_plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    main="",
    key=TRUE,
    keysize = 0.3,
    Colv = FALSE,
    scale="row",
    colsep = c(6,15,20,25,31,37,41,44),
    sepwidth = c(0.3,0.3),
    labRow=rownames(to_plot),
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.5,
    cexCol = 1,
    offsetCol = 0,
    margins = c(10, 10)
  )



######################################################################################
#################### clean workspace and start EVS analysis ##########################

samples <- read.table("sample-list.csv", header = TRUE) 

#recoding of the groups for simplicity:
samples$treatment %>% 
  recode(control_control_condition_control_additive = "intact",
         control_carrier_control_additive = "carrier",
         stem_cells_control_condition_control_additive = "SC_native",
         stem_cells_graphenea_control_additive = "SC_grafenea",
         EVS_graphenea_control_additive = "EVS_grafenea",
         EVS_go_au_control_additive = "EVS_GO_AU",
         stem_cells_go_au_control_additive = "SC_GO_AU",
         EVS_control_condition_control_additive = "EVS_native"
  ) -> samples$treatment

samples$group <- paste(samples$damage, samples$treatment, sep = "_")

samples %>% select(sample_id, file, damage, treatment, group) -> samples

grafenea_evs_samples <- scan("samples-grafens-evs.txt", what=character(), sep="\n")

samples %>% filter(
  samples$file%in% grafenea_evs_samples) -> samples_evs_grafenea

#the groups included in this analysis are: 
table(samples_evs_grafenea$group)

#### FIRST ANALYSIS - TWO-WAY ANOVA
#samples_evs_grafenea %>% filter(group %in% c())

#filter the samples_evs_grafenea$group to have groups that can be involved in two-way 



#factor 1 - damage
#factor 2 - treatment

#the intact group cannot be considered in this model 


fpkm <- data.frame(read.table('genes.fpkm_table.tsv',header = TRUE, sep = "\t"))
rnor.gene.list <- data.frame(read.delim("rnor.genes.tsv"))
gene.name <- rnor.gene.list$Gene.name[match(fpkm$tracking_id, rnor.gene.list$Gene.stable.ID)]
fpkm <- data.frame(gene.name, fpkm)
rownames(fpkm) <- fpkm$tracking_id
colnames(fpkm)[2] <- c("transcript.ID")
rm(rnor.gene.list)

fpkms.normalised <- data.matrix(fpkm[,c(-1,-2)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised)
fpkms.log <- log2(fpkms.normalised + 1)
rm(fpkms.normalised)

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA
fpkms.log <- data.frame(fpkms.log)
colnames(fpkms.log) <- colnames(fpkm[,c(3:117)])
rownames(fpkms.log) <- rownames(fpkm)

#ANOVA analysis 2-way:
results <- data.frame(fpkm$gene.name, fpkm$transcript.ID, row.names = rownames(fpkm))

two.way <- function(counts, damage, treatment) {
  if (is.na(counts[1])) {
    c(rep(NA, 3))
  } else {
    aov(as.numeric(counts) ~ as.factor(damage)*as.factor(treatment)) %>%
      summary() %>% unlist() -> res
    
    res[c("Pr(>F)1","Pr(>F)2","Pr(>F)3")]
  }
}


apply(fpkms.log[,match(samples_evs_grafenea$file, colnames(fpkms.log))],
      1,
      two.way,
      damage=samples_evs_grafenea$damage,
      treatment=samples_evs_grafenea$treatment) %>% 
  t() %>%
  apply(.,
        2,
        p.adjust,
        method = 'fdr'
  ) %>%
  cbind(results, .) -> results

### perform post-hoc tests vs mia ###


pairwise.t.test(
  as.numeric(fpkms.log[10,samples_evs_grafenea$file]),
  samples_evs_grafenea$group,
  p.adjust.method = 'none') %>% unlist()


stat_paired_t <- function(x) {
  if (is.na(x[1])) { c(NA, NA, NA, NA, NA)
  } else { 
    pairwise.t.test(
      as.numeric(x),
      samples_evs_grafenea$group,
      p.adjust.method = 'none') %>% 
      unlist() -> res
    res[(c("p.value1", "p.value11", "p.value10", "p.value8", "p.value9"))] %>%
      p.adjust(method = 'bonferroni')
  }}

apply(fpkms.log[,match(samples_evs_grafenea$file, colnames(fpkms.log))],
      1,
      stat_paired_t) %>% t() %>%
  cbind(results, .) -> results

colnames(results) <- c("gene_symbol",
                       "transcript_id",
                       "FDR_damage",
                       "FDR_treatment",
                       "FDR_int",
                       "t_test_intact_vs_mia",
                       "t_test_nacl_vs_mia",
                       "t_test_EVS_vs_mia",
                       "t_test_EVS_GOAU_vs_mia",
                       "t_test_EVS_grafenea_vs_mia")

#add a max dispersion factor:

get.max.dispersion <- function(x) {
  x %>% cbind(., samples_evs_grafenea$group[order(samples_evs_grafenea$damage, samples_evs_grafenea$group)]) %>% 
    data.frame %>% rownames_to_column() %>%  
    set_colnames(c("file", "counts", "group")) %>%
    group_by(group) %>% 
    summarise(mean.to.sd = mean(as.numeric(counts))/sd(as.numeric(counts))) %>%
    select(mean.to.sd) %>% max() -> y
  return(y)
}

#compute selected fold changes between important groups:
fpkms.log[match(
  results$transcript_id,
  rownames(fpkms.log)),
  samples_evs_grafenea$file[order(samples_evs_grafenea$damage, samples_evs_grafenea$group)]] %>%
  apply(., 1, get.max.dispersion) -> results$max.mean.to.sd

### get direction and fold changes for important comparisons ### (note - the intact group did not go into statistics)

compute_folds <- function(df, sample_info, sample_names_re, factors) {
  sample_info_grouped <-
    sample_info %>%
    group_by(across(all_of(factors)))
  groups <-
    sample_info_grouped %>%
    group_split()
  group_keys <-
    sample_info_grouped %>%
    group_keys() %>% 
    unite('name', everything(), remove = FALSE)
  
  for (tb_idx in 1:length(groups)) {
    tb <- groups[[tb_idx]]
    tb_key <- group_keys[tb_idx, 'name'] %>% unlist()
    df <-
      df %>%
      mutate(
        "{tb_key}_mean" := apply(
          select(., matches(sample_names_re)),
          1,
          function(x, names) {mean(x[names])},
          names = tb$id
        )
      )
  }
  df
}

samples_evs_grafenea$id <- samples_evs_grafenea$file

compute_folds(
  fpkms.log,
  samples_evs_grafenea,
  'Cai',
  factors = c('group')
) -> group_means

colnames(group_means)
### add information about direction and fold


results$fold_intact_vs_mia <- abs(group_means$intact_intact_mean-group_means$MIA_carrier_mean)
results$dir_intact_vs_mia <- ((group_means$intact_intact_mean-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "UP", "TRUE" = "DOWN")

results$fold_mia_vs_sc <- abs(group_means$MIA_EVS_native_mean-group_means$MIA_carrier_mean)
results$dir_mia_vs_sc <- ((group_means$MIA_EVS_native_mean-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "DOWN", "TRUE" = "UP")

results$fold_mia_vs_sc_go_au <- abs(group_means$MIA_EVS_GO_AU_mean-group_means$MIA_carrier_mean)
results$dir_mia_vs_sc_go_au <- ((group_means$MIA_EVS_GO_AU_mean-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "DOWN", "TRUE" = "UP")

results$fold_mia_vs_graf <- abs(group_means$MIA_EVS_grafenea_mean-group_means$MIA_carrier_mean)
results$dir_mia_vs_graf <- ((group_means$MIA_EVS_grafenea-group_means$MIA_carrier_mean) > 0) %>%
  as.character() %>%
  recode("FALSE" = "DOWN", "TRUE" = "UP")


# bind results with counts
results <- cbind(results,fpkms.log)

###select three genes verified in PCR for export:
results %>% filter(results$gene_symbol %in% c('Comp', 'Mmp3', 'Mmp13')) -> selected
write_csv(selected, 'evs_pcr_genes_results.csv')

table(samples_evs_grafenea$group)

###save top regulated genes:

results %>% 
  filter(FDR_damage < 0.1) %>% nrow()

results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_EVS_vs_mia, t_test_EVS_GOAU_vs_mia, t_test_EVS_grafenea_vs_mia), ~ . < 0.05)) %>%
  write_csv(., 'EVS_regulated_genes_results.csv')

results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_EVS_vs_mia, t_test_EVS_GOAU_vs_mia, t_test_EVS_grafenea_vs_mia), ~ . < 0.05)) %>%
  select(gene_symbol) %>%
  write.table(row.names = FALSE, quote = FALSE)



### plot top results ###
reordering <- order(samples_evs_grafenea$damage, samples_evs_grafenea$treatment)

to_plot <- results %>%
  filter(results$gene_symbol %in% c('Comp', 'Mmp3', 'Mmp13')) %>%
  select(samples_evs_grafenea$file[reordering])

rownames(to_plot) <- results %>% 
  filter(results$gene_symbol %in% c('Comp', 'Mmp3', 'Mmp13')) %>% 
  select(gene_symbol) %>% t()


to_plot <- results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_EVS_vs_mia, t_test_EVS_GOAU_vs_mia, t_test_EVS_grafenea_vs_mia), ~ . < 0.05)) %>%
  select(samples_evs_grafenea$file[reordering])

rownames(to_plot) <- results %>% 
  filter(FDR_damage < 0.1) %>% 
  filter(if_any(c(t_test_EVS_vs_mia, t_test_EVS_GOAU_vs_mia, t_test_EVS_grafenea_vs_mia), ~ . < 0.05)) %>%
  select(gene_symbol) %>% t()


#### actual plotting ####

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(as.character(samples_evs_grafenea$group[reordering]))

col.labels <- c(rep("", 3),
                group.names[1], rep(" ", 6), 
                group.names[2], rep(" ", 6),
                group.names[3], rep(" ", 4),
                group.names[4], rep(" ", 4),
                group.names[5], rep(" ", 4),
                group.names[6], rep(" ", 4),
                group.names[7], rep(" ", 2),
                group.names[8], rep(" ", 3),
                group.names[9], rep(" ", 2))


cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


to_plot %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to_plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    main="",
    key=TRUE,
    keysize = 0.3,
    Colv = FALSE,
    scale="row",
    colsep = c(6,15,20,25,30,36,41,44),
    sepwidth = c(0.3,0.3),
    labRow=rownames(to_plot),
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 1,
    cexCol = 1,
    offsetCol = -0.1,
    margins = c(10, 10)
  )





