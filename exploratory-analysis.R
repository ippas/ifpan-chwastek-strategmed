require(edgeR)
require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(dplyr)
require(reshape2)



#create sample info
samples <- data.frame(read.table('samples.table.tsv',header = TRUE, sep = "\t"))[,c(1,2)]


samples$file %>% str_split("/") %>% data.frame() -> temp.sample.names
temp.sample.names[8,] %>% t() %>% str_replace_all(".cxb", "") -> samples$file
rm(temp.sample.names)

samples$animal <- samples$file %>% str_replace_all("Cai", "")

groups <- data.frame(read.table('samples-treatments.csv',header = TRUE, sep = "\t"))


samples$damage <- groups$damage[match(samples$animal, groups$animal)]

samples$main <- groups$main_ingredient[match(samples$animal, groups$animal)]
samples$condition <- groups$culture_condition[match(samples$animal, groups$animal)]
samples$additive <- groups$additive[match(samples$animal, groups$animal)]

samples$group <- paste(samples$damage, samples$main, samples$condition, samples$additive, sep="_")
samples <- samples[samples$main != "EVS_less",]

samples$treatment <- paste(samples$main, samples$condition, samples$additive, sep="_")

groups <- data.frame(read.table('samples-treatments.csv',header = TRUE, sep = "\t"))

vf <- data.frame(read.table('von-frey.csv', header = TRUE, sep = "\t"))
samples$von.frey.14 <- vf$D14[match(samples$animal, vf$rat)]
samples$von.frey.42 <- vf$D42[match(samples$animal, vf$rat)]


ggplot(samples, aes(x=group, y=von.frey.42, color=as.numeric(as.factor(group)))) + geom_boxplot() + ylim(c(0,50)) + geom_point(size = 3)


write.table(samples, "sample-list.csv", row.names = FALSE) 



fpkm <- data.frame(read.table('genes.fpkm_table.tsv',header = TRUE, sep = "\t"))
rnor.gene.list <- data.frame(read.delim("rnor.genes.tsv"))
gene.name <- rnor.gene.list$Gene.name[match(fpkm$tracking_id, rnor.gene.list$Gene.stable.ID)]
fpkm <- data.frame(gene.name, fpkm)
rownames(fpkm) <- fpkm$tracking_id
colnames(fpkm)[2] <- c("transcript.ID")

rm(rnor.gene.list)

#normalize the distribution 
fpkms.normalised <- data.matrix(fpkm[,c(-1,-2)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised,copy=FALSE)
fpkms.log <- log2(fpkms.normalised + 1)

rm(fpkms.normalised)

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA


#remove reads from samples of EVS smaller group:
fpkms.log <- fpkms.log[, !colnames(fpkms.log) %in% c(samples$file[samples$main == "EVS_less"])]


#ANOVA analysis:

results <- data.frame(fpkm$gene.name, fpkm$transcript.ID, row.names = rownames(fpkm))

#one-way:

stat.one.way <- function(counts, groups) {
  ifelse(
    is.na(counts[1]),
    NA,
    unlist(summary(aov(counts ~ as.factor(groups))))[9]) }


#perform one-way stat on all groups on damage

apply(fpkms.log[,match(samples$file, colnames(fpkms.log))],
      1,
      stat.one.way,
      groups=samples$damage) %>% 
  p.adjust(method="fdr") -> results$fdr.one.way.damage


#filter results 
results %>% filter(fdr.one.way.damage < 0.1) -> results.filtered


#perform one-way stat on all groups on groups

apply(fpkms.log[results.filtered$fpkm.transcript.ID,match(samples$file, colnames(fpkms.log))],
      1,
      stat.one.way,
      groups=samples$group) %>% 
  p.adjust(method="fdr") -> results.filtered$fdr.one.way.group

results.filtered %>% filter(fdr.one.way.group < 0.1) -> results.filtered

#perform fold analysis - where there is na-cl group - use NaCl group. Otherwise use the intact group


fold.change.special <- function(x, group, ctrl) {
  abs(mean(as.numeric(x[match(
            samples$file[samples$group == group], colnames(fpkms.log))]))
      -
        mean(as.numeric(x[match(
          samples$file[samples$group == ctrl], colnames(fpkms.log))])))
}

#where there is na-cl group - use NaCl group. Otherwise use the intact group

treatments <- sort(unique(samples$treatment))

samples %>% filter(damage == "NaCl") %>% select(treatment) %>% as.vector() %>% t %>% sort() %>%
  unique() -> treatments.with.nacl

treatments <- treatments[!(treatments %in% treatments.with.nacl)]
treatments <- treatments[-c(1)]

for (s in treatments.with.nacl) {
  a <- unique(samples$group[samples$treatment == s])[1]
  b <- unique(samples$group[samples$treatment == s])[2]
  counts <- fpkms.log[results.filtered$fpkm.transcript.ID,]
  results.filtered %>% mutate(!!sym(s) := apply(counts, 1, fold.change.special, group = a, ctrl = b)) -> results.filtered
}
  
for (s in treatments) {
  a <- unique(samples$group[samples$treatment == s])
  b <- unique(samples$group[samples$damage == "intact"])
  counts <- fpkms.log[results.filtered$fpkm.transcript.ID,]
  results.filtered %>% mutate(!!sym(s) := apply(counts, 1, fold.change.special, group = a, ctrl = b)) -> results.filtered
}

all.treatments <- c(treatments.with.nacl, treatments)


results.filtered %>% filter_at(., all.treatments, any_vars(. > 2)) -> results.filtered

#add gene-max dispersion:

get.max.dispersion <- function(x) {
  x %>% cbind(., samples$group[order(samples$damage, samples$group)]) %>% 
  data.frame %>% rownames_to_column() %>%  
  set_colnames(c("file", "counts", "group")) %>% group_by(group) %>% summarise(mean.to.sd = mean(as.numeric(counts))/sd(as.numeric(counts))) %>%
  select(mean.to.sd) %>% max() -> y
  return(y)
}


fpkms.log[match(
  results.filtered$fpkm.transcript.ID,
  rownames(fpkms.log)),
  samples$file[order(samples$damage, samples$group)]] %>%
  apply(., 1, get.max.dispersion) -> results.filtered$max.mean.to.sd


# results.filtered %>% filter(max.mean.to.sd > 10) -> results.filtered

# average the results per group:

fpkms.log[match(
  results.filtered$fpkm.transcript.ID,
  rownames(fpkms.log)),
  samples$file[order(samples$damage, samples$group)]] %>% t() %>%
  data.frame() %>% mutate(group = samples$group[match(rownames(.), samples$file)]) %>%
  group_by(group) %>% summarise_all(list(mean)) %>% t() -> to.plot

colnames(to.plot) <- to.plot[1,]

to.plot <- to.plot[-1,] %>% data.frame

gene.names.to.plot <- rownames(to.plot)
to.plot <- mutate_all(to.plot, as.numeric)



#plot results:

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(samples$group[order(samples$damage, samples$group)])

cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}



to.plot %>% data.matrix() %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 2.5) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    # scale="row",
    #colsep = as.numeric(column.separators),
    sepwidth = c(0.3,0.3),
    labRow=fpkm$gene.name[match(gene.names.to.plot, rownames(fpkm))],
    labCol=group.names,         
    srtCol = 45,
    cexRow = 0.5,
    cexCol = 0.8,
    offsetCol = -0.95,
    margins = c(18,5)
  )









#export gene list

write.table(fpkm$gene.name[match(gene.names.to.plot, rownames(fpkm))], 'top-gene-list', 
            quote=FALSE, 
            row.names = FALSE,
            col.names = FALSE)









# calculate correlations of counts

correlation <- function(x) {
  y <- samples$von.frey.14[order(samples$damage, samples$group)]
  cor.test(x, y, method="pearson") %>% unlist %>% as.numeric() -> temp
  return(temp[3])
}

fpkms.log[match(
  results.filtered$fpkm.transcript.ID,
  rownames(fpkms.log)),
  samples$file[order(samples$damage, samples$group)]] %>%
  apply(., 1, correlation) -> results.filtered$corr.14

correlation <- function(x) {
  y <- samples$von.frey.42[order(samples$damage, samples$group)]
  cor.test(x, y, method="pearson") %>% unlist %>% as.numeric() -> temp
  return(temp[3])
}

fpkms.log[match(
  results.filtered$fpkm.transcript.ID,
  rownames(fpkms.log)),
  samples$file[order(samples$damage, samples$group)]] %>%
  apply(., 1, correlation) -> results.filtered$corr.42

#add counts

cbind(results.filtered, to.plot) -> results.to.export
write.table(results.to.export, "ifpan-chwastek-strategmed-results.tsv", row.names = FALSE)


#plot correlation of counts


results.filtered %>% filter(corr.14 < 0.1) -> results.filtered

to.plot <- fpkms.log[match(
  results.filtered$fpkm.transcript.ID,
  rownames(fpkms.log)),
  samples$file[order(samples$damage, samples$group)]] %>% na.omit()


to.plot %>% apply(1, scale) %>% t %>% apply(1, cut.threshold, threshold = 2.5) %>% t %>% 
  `colnames<-`(colnames(to.plot)) %>% t %>% data.frame %>% 
  mutate(sampleid = rownames(.)) %>% 
  melt(id.vars=c("sampleid")) %>%
  mutate(von.frey.42 = samples$von.frey.42[match(sampleid, samples$file)],
  von.frey.14 = samples$von.frey.14[match(sampleid, samples$file)],
  gene.name = results.filtered$fpkm.gene.name[match(variable, results.filtered$fpkm.transcript.ID)]) -> to.plot.corr

ggplot(to.plot.corr, aes(x = value, y = von.frey.14, color = gene.name)) + 
  geom_point() + 
  stat_smooth(data = to.plot.corr, aes(x = value, y = von.frey.14), method = "lm", formula = y ~ x, size = 0.2, color = "black")
  
  
  
  


