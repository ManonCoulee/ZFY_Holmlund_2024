####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'reshape', 'TCseq', 'statmod',
  'GenomicFeatures','FactoMineR')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

wd = '.'
sample_dir = "~/Documents/Mouse/RNAseq/Samples"
sample_dir = "~/Documents/Human/RNAseq/Samples/SPZ_mRNA"
sample_dir = "/media/mcoulee/F29A8FEE9A8FADA3/Hayden/RNAseq/BAM"
dirs = list.dirs(sample_dir, full.names = T, recursive = F)[1]
# files = list.files(dirs,"RPG")
files = list.files(sample_dir,"RPG")

# files = list.files(dirs,"RPG")

rout = "../MSCI/CPM_0"
rout = "."

LFC = log2(1.5)

CT = list.dirs(sample_dir, full.names = F, recursive = F)
nk = 8

## Load SamplePlan data descriptor
SamplePlan = read.table(paste0(wd,"/SamplePlan.tsv"),sep = "\t", h = T, row.names = 1)
SamplePlan$SampleType = factor(SamplePlan$SampleType)
# SamplePlan$CellType = factor(SamplePlan$CellType,levels = CT)
SamplePlan$CellType = factor(SamplePlan$CellType)
# SamplePlan$CellType = factor(SamplePlan$CellType,levels = c("Kit_m_spike","Kit_p_spike"))
SamplePlan$SamplePool = factor(SamplePlan$SamplePool)
SamplePlan$CellTypeRed = factor(SamplePlan$CellTypeRed)
SamplePlan$BatchEffect = factor(SamplePlan$BatchEffect)

####################################################################################################
##                                              SUMMARY STATS
####################################################################################################

## Mapping summary from STAR log
starlog = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)),
  function(path){
    x = scan(paste0(path,".STAR.log"), what = 'character', sep = "\n")[c(6,11,12,5,8,23,25,32)]
    as.numeric(do.call(rbind, strsplit(x = x, split = "|\t", fixed = T))[, 2])
  }, simplify = T), row.names = c('average_read_length', '#mapped_splice_sites',
    '#annotated_splice_sites', '#sequenced_reads', '#uniquely_mapped_reads',
    '#multi_mapped_reads_kept', '#multi_mapped_reads_discarded', '#chimeric_reads'))

starlog['#unmapped_reads',] = starlog['#sequenced_reads',] - colSums(starlog[c(grep('mapped_reads',
  rownames(starlog), value = T), '#chimeric_reads'),])
names(starlog) = rownames(SamplePlan)
write.table(file = file.path(rout, 'STAR-stats.tsv'),
  x = data.frame(Feature = rownames(starlog), starlog, check.names = F), sep = "\t",
  row.names = F, quote = F)

## Plot STAR summary statistics
o = c('unmapped','uniquely mapped','multi-mapped (kept)','multi-mapped (discarded)')
o = melt(data.frame(Read = factor(x = o, levels = o),
  100*t(t(starlog[c(9,5:7),])/as.numeric(starlog[4,]))))
p = ggplot(
  data = o, aes(
    x = factor(x = variable, labels = with(SamplePlan, paste(CellType, SampleType, SamplePool))),
    y = value, fill = Read, order = Read
  )) +
  geom_bar(stat = 'identity', width = .95) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_fill_manual(name = NULL, values = brewer.pal(9, 'Set1')[c(1, 3:5)]) +
  labs(
    title = 'Distribution of sequenced reads (STAR)',
    y = 'Relative read counts (%)'
  ) +
  theme_minimal() + theme(
    plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
    legend.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5)),
    legend.position = 'bottom', legend.direction = 'horizontal',
    axis.title = element_text(size = rel(1.5)), axis.title.x = element_blank(),
    axis.text = element_text(size = rel(1.5)), axis.text.x = element_text(angle = 45, vjust = 1,
      hjust = 1)
)
# ggsave(file.path(rout, 'STAR-stats-1.pdf'), plot = p, width = 10, height = 8)
ggsave(file.path(rout, 'STAR-stats-1.png'), device = 'png', plot = p, width = 10, height = 8,
  dpi = 150)

## Plot distribution gene identify
o = c('known', 'unknown')
o = melt(data.frame(Splice = factor(x = o, levels = o),100*t(t(rbind(starlog[3,],
  apply(starlog[3:2,],2,diff)))/as.numeric(starlog[2,]))))
p = ggplot(
  data = o, aes(
    x = factor(x = variable, labels = with(SamplePlan, paste(CellType, SampleType, SamplePool))),
    y = value, fill = Splice, order = Splice)) +
  geom_bar(stat = 'identity', width = .95) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(name = 'Gencode v38', values = brewer.pal(9,'Set1')[c(3,1)]) +
  # scale_fill_manual(name = 'Gencode vM19', values = brewer.pal(9,'Set1')[c(3,1)]) +
  labs(title = 'Distribution of splice sites (STAR)',y = 'Relative #covered sites (%)') +
  theme_minimal() + theme(
    plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.position = 'bottom', legend.direction = 'horizontal',
    axis.title = element_text(size = rel(1.5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = rel(1.5)),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
)
# ggsave(file.path(rout, 'STAR-stats-2.pdf'), plot = p, width = 10, height = 8)
ggsave(file.path(rout, 'STAR-stats-2.png'), device = 'png', plot = p, width = 10, height = 8,
  dpi = 150)

## Mapping summary from idxstats
idxstats = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)),
  function(path){
    x = read.table(paste0(path,'.STAR.idxstats'), sep = "\t", header = F,row.names = 1)
    print(head(x))
    x[grep(pattern = 'chr', x = rownames(x)), 2]
  }, simplify = T), row.names = c(paste0('chr', c(1:22, 'X', 'Y')), 'MT'))
names(idxstats) = rownames(SamplePlan)

write.table(file.path(rout, 'IDX-stats.tsv'),
  x = data.frame(CHR = rownames(idxstats),idxstats, check.names = F),
  sep = "\t", row.names = F, quote = F)

## Plot SAMTools idxstats
p = ggplot(
  data = melt(data.frame(Map = factor(x = rownames(idxstats), levels = rownames(idxstats)),
    100*t(t(idxstats)/colSums(idxstats)))),
  aes(x = factor(x = variable, labels = with(SamplePlan, paste(CellType, SampleType, SamplePool))),
    y = value, fill = Map, order = Map)) +
  geom_bar(stat = 'identity', width = .95) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = c(rev(brewer.pal(12,'Paired')[-11]), brewer.pal(12,'Paired')[-1]),
    name = 'Color map') +
  labs(title = 'Distribution of mapped reads (SAMTools)',y = 'Relative read counts (%)') +
  theme_minimal() + theme(
    plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.position = 'right', legend.direction = 'vertical',
    axis.title = element_text(size = rel(1.5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = rel(1.5)),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
# ggsave(file.path(rout, 'IDX-stats.pdf'), plot = p, width = 12, height = 8)
ggsave(file.path(rout,'IDX-stats.png'), device = 'png', plot = p, width = 12, height = 8, dpi = 150)

####################################################################################################
##                                           ANNOTATION
####################################################################################################

## Gencode gene counts from STAR
x = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)), function(path){
  read.table(paste0(path, '.RPG.tsv'), sep = "\t", header = F, skip = 4)[,4]
}, simplify = T),
  row.names = read.table(paste0(SamplePlan[1,"SamplePath"],"/",rownames(SamplePlan)[1],".RPG.tsv"),
  sep = "\t", header = F, skip = 4)[,1])

write.table(x,paste0(rout,"/genecounts_raw.tsv"), quote  = F, sep = "\t", col.names = T, row.names = T)

## Gene counts data normalization with CPM
xx = cpm(y = x, normalized.lib.sizes = T, log = F)
write.table(xx,paste0(rout,"/genecounts_cpm.tsv"), quote  = F, sep = "\t", col.names = T, row.names = T)

min_sample = 2
min_expression = 1

## Filter-out low expressed genes (cpm > 1 in minimum 2 samples)
expressed_genes = rowSums(xx > min_expression) >= min_sample
rlog_sample = SamplePlan
rownames(rlog_sample) = paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan))

genecounts = list(
  ## Raw read counts (as estimated by STAR)
  raw = x[expressed_genes, ],
  ## For expression data visualization
  cpm = cpm(y = x[expressed_genes, ], lib.size = colSums(x), normalized.lib.sizes = T, log = F),
  ## For performing PCA analysis
  rlog = assay(rlog(DESeqDataSetFromMatrix(countData = x,
      colData = rlog_sample,
      design = ~1)))[expressed_genes, ]
)

## Distribution read alignment in each cell type
# read_distribution = genecounts$cpm
# table_variance = data.frame()
# distribution = lapply(CT, function(x) {
#   print(x)
#   read_distribution_sample = read_distribution[,grepl(x,colnames(read_distribution))]
#   list_condition = lapply(c("KO","Ctl"), function(condition) {
#     read_distribution_sample_condition = read_distribution_sample[,grepl(condition,
#       colnames(read_distribution_sample), ignore.case = T)]
#     read_distribution_all = unlist(as.list(read_distribution_sample_condition))
#     read_distribution_variance = var(read_distribution_all)
#     print(paste0(condition," : ",read_distribution_variance))
#
#     return(read_distribution_all)
#   })
#   return(list_condition)
# })

## Load gene annotations (from Genecode transcripts.fa)
gene_annotation = data.frame(do.call(rbind, strsplit(x = gsub(
  pattern = "^>", replacement = '', x = unlist(system(paste('zgrep -P "^>"',
    file.path('../gencode.vM19.transcripts.fa.gz')), intern = T)),
    # file.path('~/Documents/Annotations/gencode.v38.transcripts.fa.gz')), intern = T)),
  perl = T), split = '|', fixed = T)), stringsAsFactors = F)
names(gene_annotation) = c('tx_ENS', 'gene_ENS', 'gene_OTT', 'tx_OTT', 'tx_name', 'gene_name',
  'tx_length', 'gene_type')
gene_annotation[gene_annotation == '-'] = NA
gene_annotation$tx_length = as.numeric(gene_annotation$tx_length)
bed_annotation = read.table("../gencode.vM19.annotation.bed", sep = "\t",
  h = F)
# bed_annotation = read.table("~/Documents/Annotations/gencode.v38.annotation.bed", sep = "\t",
#   h = F)[,1:6]
colnames(bed_annotation) = c("chr","start","end","strand","score","gene_ENS") ## Mouse
# colnames(bed_annotation) = c("chr","start","end","gene_ENS","score","strand") ## Human

gene_annotation = merge(bed_annotation,gene_annotation, by = "gene_ENS")

## Estimating gene_length by using GenomicFeatures
txdb = makeTxDbFromGFF(file.path('../gencode.vM19.annotation.gtf'),
# txdb = makeTxDbFromGFF(file.path('~/Documents/Annotations/gencode.v38.annotation.gtf'),
  format = 'gtf')
exons.list.per.gene = exonsBy(x = txdb, by = 'gene')
## Same but much faster
exonic.gene.sizes = sum(width(reduce(exons.list.per.gene)))
genecounts$rpkm = rpkm(y = x, gene.length = exonic.gene.sizes[match(names(exonic.gene.sizes),
  rownames(x))], normalized.lib.sizes = T, log = F)[expressed_genes, ]

write.table(genecounts$rpkm,"genecounts_rpkm.tsv", quote  = F, sep = "\t", col.names = T,
  row.names = T)

####################################################################################################
##                                        FILES GENERATION
####################################################################################################

## Create file for next GSEA analysis
table = data.frame(genecounts$cpm)
colnames(table) = rownames(SamplePlan)
table$Name = rownames(genecounts$raw)
table$Description = NA

for(dir in unique(SamplePlan$CellType)) {
  name = rownames(SamplePlan[SamplePlan$CellType == dir,])
  print(name)

  x_CTL = table[,name[grepl("WT",name, ignore.case = T)]]
  x_CTL$gene_expression = rowSums(x_CTL)/length(x_CTL)
  x_CTL$gene_ENS = rownames(x_CTL)
  write.table(x_CTL,paste0(dir,"_WT_gene_expression.tsv"),row.names = F, col.names = T, quote = F,
    sep = "\t")

  x_KO = table[,name[grepl("KO",name)]]
  x_KO$gene_expression = rowSums(x_KO)/length(x_KO)
  x_KO$gene_ENS = rownames(x_KO)
  write.table(x_KO,paste0(dir,"_KO_gene_expression.tsv"),row.names = F, col.names = T, quote = F,
    sep = "\t")

  xx = table[,c("Name","Description",name)]
  # xx = table[,c("Name","Description",name[c(2,1,3,4,6,5)])] # SPIN1

  write.table(xx,paste0(dir,"_genecounts.gct"),row.names = F, col.names = T, quote = F,
    sep = "\t")
}

## Create a file for TC analysis
colnames(genecounts$raw) = rownames(SamplePlan)
write.table(genecounts$raw,paste0(rout,"/genecounts_filtered_raw.tsv"),row.names = T, col.names = T, quote = F,
sep = "\t")

# write.table(gene_annotation,"~/Documents/Annotations/human_gene_annotation.tsv",row.names = F,
#   col.names = T, quote = F,sep = "\t")
# write.table(exonic.gene.sizes,"~/Documents/Annotations/human_gene_length.tsv",row.names = T,
#   col.names= F, quote = F,sep = "\t")

####################################################################################################
##                                      SAMPLE DISTRIBUTION
####################################################################################################

## Perform a PCA to get an overview of overall sample distribution
x = PCA(X = t(genecounts$rlog), scale.unit = TRUE, graph = FALSE, ncp = 3)
SamplePCA = list(
  map = x$ind$coord,
  var = as.numeric(format(x$eig[1:3, 2], digits = 2, nsmall = 2, trim = TRUE)),
  cor = x$var$cor
)
colnames(SamplePCA$map) = colnames(SamplePCA$cor) = names(SamplePCA$var) = paste('Comp.', 1:3)
write.table(file.path(rout,'PCA-map.tsv'),
  x = data.frame(SampleID = rownames(SamplePCA$map), SamplePCA$map, check.names = F),
  sep = "\t", row.names = F, quote = F)

m_export = unique(subset(gene_annotation[, c(1:6)]))
m_export = unique(subset(gene_annotation[, c(2:6,grep('gene', x = names(gene_annotation)))],
  select = -gene_type))
m_export = m_export[match(rownames(SamplePCA$cor), m_export$gene_ENS), ]
write.table(file.path(rout, 'PCA-cor.tsv'),x = data.frame(m_export, SamplePCA$cor, check.names = F),
  sep = "\t", row.names = F, quote = F)

## Plot Individuals Factor map (first three pairwize component comparison)
invisible(lapply(1:length((L = list(c(1, 2), c(2, 3), c(1, 3)))), function(i){
  m = as.data.frame(SamplePCA$map[, -setdiff(1:3,L[[i]])]); names(m)[1:2] = c('x', 'y')
  p = ggplot(data = m, aes(x = x, y = y, shape = SamplePlan$SampleType,color=SamplePlan$CellType)) +
    geom_point(size = 10, alpha = .5) +
    geom_vline(xintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
    #geom_text(aes(x = x, y = y), label = SamplePlan$SamplePool, fontface = 'bold', size = 4,
    #  alpha = .5, hjust = .5, vjust = .5, show_guide = FALSE) +
    scale_color_manual(name = 'CellType', values = brewer.pal(9, 'Set1')) +
    scale_shape_manual(name = 'SampleType', values = 15:18) +
    labs(
      title = 'Individuals Factor Map (PCA)',
      subtitle = paste('Comp.', L[[i]][1], 'vs.', L[[i]][2],
        paste0('(rlog normalized read counts / cpm>',min_expression,' / n=',nrow(SamplePCA$cor),')')),
      x = paste('Comp.', L[[i]][1], paste0('(', SamplePCA$var[L[[i]][1]], '%)')),
      y = paste('Comp.', L[[i]][2], paste0('(', SamplePCA$var[L[[i]][2]], '%)'))
    ) +
    theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )
  # ggsave(filename = file.path(rout, paste0('PCA', L[[i]][1], '-', L[[i]][2], '.pdf')),
  #   plot = p, width = 12, height = 8)
  ggsave(file.path(rout, paste0('PCA', L[[i]][1], '-', L[[i]][2], '.png')), plot = p, width = 12,
    height = 8, device = 'png', dpi = 150)
}))

