####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'FactoMineR', 'reshape', 'statmod',
  'GenomicFeatures')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

wd = 'RNAseq_B1_B2_FC1.5'
sample_dir = "Samples"
dirs = list.dirs(sample_dir, full.names = T, recursive = F)[1]
files = list.files(sample_dir,"RPG")

rout = wd

LFC = log2(1.5)

# CT = list.dirs(sample_dir, full.names = F, recursive = F)[-c(1:2,5,10)]
CT = list.dirs(sample_dir, full.names = F, recursive = F)[1:2]
CT = list.dirs(sample_dir, full.names = F, recursive = F)[5]
nk = 8
min_expression = 1

## Load SamplePlan data descriptor
SamplePlan = read.table(paste0(wd,"/SamplePlan.tsv"),sep = "\t", h = T, row.names = 1)
SamplePlan$SampleType = factor(SamplePlan$SampleType)
SamplePlan$CellType = factor(SamplePlan$CellType)
# SamplePlan$CellType = factor(SamplePlan$CellType,levels = c("Kit_m_spike","Kit_p_spike"))
SamplePlan$SamplePool = factor(SamplePlan$SamplePool)
SamplePlan$CellTypeRed = factor(SamplePlan$CellTypeRed)
SamplePlan$SamplePath = "Samples"
SamplePlan$BatchEffect = factor(SamplePlan$BatchEffect)

## Load genecounts file
genecounts = read.table("RNAseq_B1_B2_FC1/genecounts_filtered_raw.tsv", sep = "\t", h = T, stringsAsFactors = F)
gene_annotation = read.table("RNAseq_B1_B2_FC1/ZFY1_KO_effect_within_RS_cells_volcano_table_FDR_PV.tsv", sep = "\t", h = T)[,1:8]
gene_annotation = read.table("../mouse_gene_annotation.tsv", sep = "\t", h = T,
  stringsAsFactors = F)
m_export = unique(subset(gene_annotation[, -grep('tx_', x = names(gene_annotation))],
  select = -gene_type))
rownames(m_export) = m_export$gene_ENS

## TEST NORMALIZATION
genecounts_CPM = read.table("genecounts_cpm.tsv", sep = "\t", h = T, stringsAsFactors = F)
genecounts_RPKM = read.table("genecounts_rpkm.tsv", sep = "\t", h = T, stringsAsFactors = F)

genecounts_CPM = genecounts_CPM[rownames(genecounts_RPKM),]

boxplot(log(genecounts))
boxplot(log(genecounts_CPM+0.01))
boxplot(log(genecounts_RPKM+0.01))
SamplePlan[c(1:18,31:42),]

genecounts_mel = read.table("../DE_analysis_3_factors_SC_SCII_RS_spike/genecounts_filtered_raw.tsv", sep = "\t", h = T, stringsAsFactors = F)
genecounts_CPM_mel = read.table("../DE_analysis_3_factors_SC_SCII_RS_spike/genecounts_cpm.tsv", sep = "\t", h = T, stringsAsFactors = F)
genecounts_CPM_mel = genecounts_CPM_mel[rownames(genecounts_mel),]
boxplot(log(genecounts_mel))
boxplot(log(genecounts_CPM_mel+0.01))


####################################################################################################
##                                              DE analysis
####################################################################################################

## Perform differential gene expression analysis
cell = SamplePlan$CellTypeRed
age = relevel(SamplePlan$CellType, ref = "YNG") ##Pierre
phenotype = relevel(SamplePlan$SampleType, ref = "CTL")
phenotype = relevel(SamplePlan$SampleType, ref = "Young")
phenotype = relevel(SamplePlan$SampleType, ref = "WT") ##Hayden
phenotype = relevel(SamplePlan$SampleType, ref = "CTL") ##Pierre

batch = relevel(SamplePlan$BatchEffect,ref = "B2") ##Hayden

## Create a model with all multifactor combinations
design = model.matrix(~phenotype+cell+phenotype:cell)
design = model.matrix(~phenotype)
design = model.matrix(~phenotype+cell+batch+phenotype:cell) ##Hayden
design = model.matrix(~phenotype + age + phenotype:age) ##Pierre

# colnames(design) = c("all","KOvsCTL.Kitm","CTL.Kitp","CTL.RS","CTL.SC","CTL.SCII",
#   "KOvsCTL.Kitp","KOvsCTL.RS","KOvsCTL.SC","KOvsCTL.SCII")
# colnames(design) = c("all","KOvsCTL.RS","CTL.SC","CTL.SCII","KOvsCTL.SC","KOvsCTL.SCII")
# colnames(design) = c("all","KOvsCTL.Kitm","CTL.Kitp","KOvsCTL.Kitp")
# colnames(design) = c("all","KOvsCTL.RS","CTL.SC","KOvsCTL.SC") # SPIN1
# colnames(design) = c("all","KDvsCTL") ##Laila data
# colnames(design) = c("all","KOvsCTL") ##SPINDOC
# colnames(design) = c("all","OldvsYoung") ##spermaging
colnames(design) = c("all","DKOvsWT.RS","Zfy1vsWT.RS","Zfy2vsWT.RS","CTL.SCI","CTL.SCII","batch",
 "DKOvsWT.SCI","Zfy1vsWT.SCI","Zfy2vsWT.SCI","DKOvsWT.SCII","Zfy1vsWT.SCII","Zfy2vsWT.SCII") ##Hayden
# colnames(design) = c("all","KOvsCTL.YNG","CTL.OLD","KOvsCTL.OLD") ##Pierre

# design = model.matrix(~phenotype+SamplePlan$CellType)
# colnames(design) = c("all","KOvsCTL","SPIN1")
# design = model.matrix(~phenotype)
# colnames(design) = c("all","KOvsCTL")

## Creating DGEList object
y = DGEList(counts = genecounts)
## Calculating TMM-based scaling factors
y = calcNormFactors(object = y, method = "RLE")
## Estimating data dispersion
y = estimateDisp(y = y, design = design)

## Plot BCV
o = order(y$AveLogCPM)
p = ggplot(data = data.frame(x = y$AveLogCPM, y = sqrt(y$tagwise.dispersion),
  z = sqrt(y$trended.dispersion), zz = sqrt(y$common.dispersion))) +
  geom_point(aes(x = x, y = y, color = 'Tagwise'), size = 4, alpha = .5) +
  geom_hline(aes(color = 'Common', yintercept = zz), size = 2, alpha = .8) +
  geom_line(aes(x = x[o], y = z[o], color = 'Trended'), size = 2, alpha = .8) +
  scale_color_manual(name = 'Data dispersion', values = c('Tagwise' = brewer.pal(9,'Set1')[2],
    'Common' = brewer.pal(9,'Set1')[1], 'Trended' = brewer.pal(9,'Set1')[3])) +
  labs(x = 'Average LogCPM estimates',y = 'Biological coefficient of variation',
    title = paste('BCV plot of DOT1L samples', paste0('(n=',nrow(SamplePlan),'*',nrow(y),')'))) +
  theme(
    plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.position = 'bottom', legend.direction = 'horizontal',
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5)),
    panel.border = element_rect(colour = 'grey50', fill = NA)
)
ggsave(file.path(rout, 'BCV.png'), plot = p, width = 10, height = 10, device = 'png', dpi = 150)

## Plot MDS
n_top = 500
o = plotMDS(x = y, plot = FALSE, top = n_top)
p = ggplot(data = data.frame(x = o$x, y = o$y), aes(x = x, y = y, shape = SamplePlan$SampleType,
    color = SamplePlan$CellType, alpha = SamplePlan$BatchEffect)) +
  geom_point(size = 10) +
  geom_vline(xintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
  geom_hline(yintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
  #geom_text(aes(x = x, y = y), label = SamplePlan$SamplePool, fontface = 'bold', size = 4,
  #  alpha = .5, hjust = .5, vjust = .5, show_guide = FALSE) +
  scale_color_manual(name = 'CellType', values = c(brewer.pal(9, 'Set1'))) +
  scale_shape_manual(name = 'SampleType', values = 15:18) +
  scale_alpha_manual(name = "Batch", values = c(0.5,1)) +
  labs(title = 'Multi-dimensional scaling plot (MDS)',
    subtitle = paste0('Pairwize gene selection / cpm>',min_expression,' / n=',n_top,' top genes'),
    x = 'Leading LogFC dim. 1', y = 'Leading LogFC dim. 2') +
  theme(
    plot.title = element_text(size = 35, lineheight = 2, vjust = 1, face = 'bold'),
    plot.subtitle = element_text(size = 25, face = 'bold'),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    panel.border = element_rect(colour = 'grey50', fill = NA)
)
ggsave(file.path(rout, 'MDS_batch.png'), plot = p, width = 12, height = 8, device = 'png', dpi = 300)

fit = glmQLFit(y = y, design = design, robust = T)
cmp = makeContrasts(
  # 'DOT1L KO effect within Kit- cells' = KOvsCTL.Kitm,
  # 'DOT1L KO effect within Kit+ cells' = KOvsCTL.Kitp - KOvsCTL.Kitm,
  # 'DOT1L KO effect within RS cells' = KOvsCTL.RS - KOvsCTL.Kitm,
  # 'DOT1L KO effect within SC cells' = KOvsCTL.SC - KOvsCTL.Kitm,
	# 'DOT1L KO effect within SCII cells' = KOvsCTL.SCII - KOvsCTL.Kitm,
  # 'DOT1L Kit+ vs Kit- within CTL' = CTL.Kitp,
  # 'DOT1L Kit+ vs Kit- within KO' = (KOvsCTL.Kitp+CTL.Kitp),
  # 'DOT1L SC vs Kit+ within CTL' = CTL.SC - CTL.Kitp,
  # 'DOT1L SC vs Kit+ within KO' = (KOvsCTL.SC+CTL.SC) - (KOvsCTL.Kitp+CTL.Kitp),
  # 'DOT1L SCII vs SC within CTL' = CTL.SCII - CTL.SC,
  # 'DOT1L SCII vs SC within KO' = (KOvsCTL.SCII+CTL.SCII) - (KOvsCTL.SC+CTL.SC),
  # 'DOT1L RS vs SCII within CTL' = CTL.RS - CTL.SCII,
  # 'DOT1L RS vs SCII within KO' = (KOvsCTL.RS+CTL.RS) - (KOvsCTL.SCII+CTL.SCII),
  ## Model Kit-/Kit+
  # 'DOT1L KO effect within Kit- cells' = KOvsCTL.Kitm,
  # 'DOT1L KO effect within Kit+ cells' = KOvsCTL.Kitm + KOvsCTL.Kitp,
  # 'DOT1L Kit+ vs Kit- within CTL' = CTL.Kitp,
  # 'DOT1L Kit+ vs Kit- within KO' = (KOvsCTL.Kitp+CTL.Kitp),
  ## Model SC/SCII/RS
  # 'DOT1L KO effect within RS cells' = KOvsCTL.RS,
  # 'DOT1L KO effect within SC cells' = KOvsCTL.RS + KOvsCTL.SC,
	# 'DOT1L KO effect within SCII cells' = KOvsCTL.RS + KOvsCTL.SCII,
  # 'DOT1L SCII vs SC within CTL' = CTL.SCII - CTL.SC,
  # 'DOT1L SCII vs SC within KO' = (KOvsCTL.SCII+CTL.SCII) - (KOvsCTL.SC+CTL.SC),
  # 'DOT1L RS vs SCII within CTL' = - CTL.SCII,
  # 'DOT1L RS vs SCII within KO' = - (KOvsCTL.SCII+CTL.SCII),
  ## Laila
  # 'KD effect' = KDvsCTL,
  ## SPIN1
  # 'SPIN1 KO effect within RS cells' = KOvsCTL.RS,
  # 'SPIN1 KO effect within SC cells' = KOvsCTL.RS + KOvsCTL.SC,
  ## SPINDOC
  # 'SPINDOC KO effect within testis cells' = KOvsCTL,
  ## spermaging
  # 'Age effect within testis cells' = OldvsYoung,
  ## Hayden
   'DKO effect within RS cells' = DKOvsWT.RS,
   'DKO effect within SCI cells' = DKOvsWT.RS - DKOvsWT.SCI,
   'DKO effect within SCII cells' = DKOvsWT.RS - DKOvsWT.SCII,
   'ZFY1 KO effect within RS cells' = Zfy1vsWT.RS,
   'ZFY1 KO effect within SCI cells' = Zfy1vsWT.RS - Zfy1vsWT.SCI,
   'ZFY1 KO effect within SCII cells' = Zfy1vsWT.RS - Zfy1vsWT.SCII,
   'ZFY2 KO effect within RS cells' = Zfy2vsWT.RS,
	 'ZFY2 KO effect within SCI cells' = Zfy2vsWT.RS - Zfy2vsWT.SCI,
	 'ZFY2 KO effect within SCII cells' = Zfy2vsWT.RS - Zfy2vsWT.SCII,
  ## Pierre
  #'KO effect within YNG cells' = KOvsCTL.YNG,
  #'KO effect within OLD cells' = KOvsCTL.OLD - KOvsCTL.YNG,
  #'OLD vs YNG within CTL cells' = CTL.OLD,
  #'OLD vs YNG within KO cells' = KOvsCTL.OLD + CTL.OLD,

  levels = design
)

thresholds = data.frame(
  FC = rep(1.5,length(colnames(cmp))),
  P = rep(.05,length(colnames(cmp))),
  row.names = colnames(cmp)
)

DEG = sapply(colnames(cmp), function(contrast){
  print(contrast)
  qlf = glmQLFTest(glmfit = fit, contrast = cmp[, contrast])
  FC = log2(thresholds[contrast,'FC'])
  # qlf = glmTreat(glmfit = fit, contrast = cmp[, contrast], lfc = FC)
  tt = with(topTags(object = qlf, n = NULL, sort.by = 'none', adjust.method = "BH"), table)
  tt$gl = 0
  tt$gl[tt$PValue <= thresholds[contrast, 'P']  & tt$logFC < -FC] = -1
	tt$gl[tt$PValue <= thresholds[contrast, 'P']  & tt$logFC > FC] = 1
  tt$gl[tt$FDR <= thresholds[contrast, 'P']  & tt$logFC < -FC] = -2
	tt$gl[tt$FDR <= thresholds[contrast, 'P']  & tt$logFC > FC] = 2

  tt$gl = factor(tt$gl, levels = c(-2,-1, 0, 1, 2),
  labels = c('Down-regulated (FDR < 5%)','Down-regulated (p-value < 5%)', 'Not regulated',
  'Up-regulated (p-value < 5%)', 'Up-regulated (FDR < 5%)'))
  # tt$gl = factor(tt$gl, levels = c(-2, 0, 2),
   #  labels = c('Down-regulated (FDR < 5%)','Not regulated', 'Up-regulated (FDR < 5%)'))
  # tt$gene_name = rownames(tt)
  # tt = merge(m_export,tt,by = "gene_name")
  tt
}, simplify = F)

invisible(lapply(names(DEG), function(x){
  write.table(file.path(rout,paste0(gsub(pattern = ' ', replacement = '_', x = x),
    '_volcano_table_FDR_PV.tsv')),
    x = data.frame(m_export[rownames(DEG[[x]]),],DEG[[x]], check.names = F), sep = "\t",
    row.names = F, quote = F)
}))

# tt = DEG[[colnames(cmp)[3]]]
# # RD2 = sample(1:nrow(tt),size = 1031)
# tt_DEG = tt[tt$gl != "Not regulated",]
# tt_DEG = tt[RD2,]
# tt_DEG$FDR = p.adjust(tt_DEG$PValue, "BH")
#
# write.table(file.path(rout,paste0(gsub(pattern = ' ', replacement = '_', x = names(DEG)[3]),
#     '_volcano_table_only_DEG.tsv')),
#     x = data.frame(m_export[rownames(tt_DEG),],tt_DEG, check.names = F), sep = "\t",
#     row.names = F, quote = F)

# bed = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed", sep = "\t", h = F)
# colnames(bed) = c("chr","start","end", "strand","score","gene_ENS")
#
# invisible(lapply(names(DEG), function(x){
#   table = data.frame(m_export[rownames(DEG[[x]]),],DEG[[x]], check.names = F)
#   DEG_merge = merge(table, bed, by.x = "gene_ENS",by.y = "gene_ENS")
#
#   write.table(file.path(rout,paste0(gsub(pattern = ' ', replacement = '_', x = x),
#     '_volcano_table.tsv')), x = DEG_merge, sep = "\t",
#     row.names = F, quote = F)
# }))

## Hayden
name = list.files(wd,pattern = "volcano_table_FDR_PV.tsv", full.names = T)
DEG = lapply(name, function(file) {
  tt = read.table(file, sep = "\t", h = T)
  l = unlist(str_split(unlist(str_split(file,"/"))[2],"_"))
  tt$cat = l[1]
  tt$cell = l[5]
  return(tt)
})
DEG[[1]]$cell = "RS"
names(DEG) = colnames(cmp)

df = do.call("rbind",DEG)
df = df[grepl("FDR",df$gl),]
p = ggplot(df[df$gl != "Not regulated" & df$cell == "DKO",],
           aes(x = factor(gl),
               alpha = factor(cell, levels = c("SCI","SCII","RS")),
               fill = gl)) +
  geom_bar(stat = "count", position = "dodge") +
  xlab("") + ylab("number of genes") + labs(fill = "", alpha = "") +
  scale_fill_manual(values = rep(c("blue","red"),5)) +
  theme_bw() + theme(strip.background  = element_blank(),
                     text = element_text(size=30, angle = 0),
                     panel.grid.major = element_line(colour = "grey80"),
                     panel.border = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.minor.x=element_blank(),
                     panel.grid.major.x=element_blank(),
                     legend.position = "bottom") +
  guides(fill = "none")
ggsave('distribution_barplot_DKO.png',
       plot = p, width = 10, height = 8, device = 'png', dpi = 300)

x = names(DEG)[6]
## DEG / MD plots
invisible(lapply(names(DEG), function(x){
  tt = DEG[[x]]
  print(table(tt$gl))
  
  ## For remove PV
  tt[tt$gl == "Down-regulated (p-value < 5%)" | tt$gl == "Up-regulated (p-value < 5%)","gl"] = "Not regulated"
  
  
  p = ggplot(data = tt) +
    geom_point(aes(x = logFC, y = -log10(PValue), color = gl), size = 5, alpha = .5) +
    geom_vline(xintercept = c(-log2(thresholds[x,'FC']), log2(thresholds[x,'FC'])), alpha = .5,
      color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = -log10(thresholds[x,'P']), alpha = .5, color = 'grey50',
      linetype = 'dashed') +
    #scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(1)],
     # brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5)])) +
    # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(2,1)],
    #   brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5,6)])) +
    # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
    scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(9)])) +
    labs(subtitle = paste(x, paste0('(*n=', sum(tt$gl != 'NR'),')')),
      title = paste(paste0('FC>',thresholds[x,'FC']),paste0('& P<',100*thresholds[x,'P'],'%'),
        'data counts'),
      x = '', y = '') +
    theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )
  #ggsave(file.path(rout, paste0(gsub(pattern = ' ',replacement = '_', x = x),'_volcano_plot_FDR_PV.png')),
  #  plot = p, width = 10, height = 8, device = 'png', dpi = 300)
   ggsave(file.path(rout, paste0(gsub(pattern = ' ',replacement = '_', x = x),'_volcano_plot_FDR.png')),
   plot = p, width = 10, height = 8, device = 'png', dpi = 300)

  df = data.frame(gl = levels(tt$gl), n = as.numeric(table(tt$gl)), stringsAsFactors = F)
  df$gl = factor(df$gl, levels = df$gl)
  df$X = factor(df$gl, levels = df$gl[order(df$n, decreasing = T)])
  h = df$n
  p = ggplot(data = data.frame(df, height = h + .001)) +
    geom_bar(aes(x = X, y = height, fill = gl), stat = 'identity', width = .98) +
    geom_text(aes(x = X, color = gl, y = height, label = paste0(n, "\n")), fontface = 'bold',
      size = 5, lineheight = .3, hjust = .5, vjust = 0, show_guide = F) +
   scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(2,1)],
       brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5,6)])) +
   scale_fill_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(2,1)],
     brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5,6)])) +
   # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
   # scale_fill_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
    labs(subtitle = paste(x, paste0('(*n=', sum(tt$gl != 'NR'),')')),
      title = paste(paste0('FC>',thresholds[x,'FC']),paste0('& P<',100*thresholds[x,'P'],'%'),
        'data counts'),
      x = '', y = '') +
    theme_void() + theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.5)),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.5))
    )
  ggsave(file.path(rout, paste0(gsub(pattern = ' ', replacement = '_', x = x),
    '_volcano_table_FDR_PV.png')), plot = p, width = 8, height = 8, device = 'png', dpi = 300)
  # ggsave(file.path(rout, paste0(gsub(pattern = ' ', replacement = '_', x = x),
  #   '_volcano_table_FDR.png')), plot = p, width = 8, height = 8, device = 'png', dpi = 300)

  p = ggplot(data = tt) +
    geom_point(aes(x = logCPM, y = logFC, color = gl), size = 5, alpha = .5) +
    geom_hline(yintercept = c(-log2(thresholds[x, 'FC']), log2(thresholds[x, 'FC'])), alpha = .5,
      color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = 0, alpha = .5, color = '#FFFFFF') +
    #scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(1)],
     # brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5)])) +
    # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(2,1)],
    #     brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5,6)])) +
    # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(1)],
    #     brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5)])) +
    # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
    scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(9)])) +
    labs(subtitle = paste(x, paste0('(*n=', sum(tt$gl != 'NR'), ')')),
      title = paste(paste0('FC>',thresholds[x,'FC']),paste0('& P<',100*thresholds[x,'P'],'%'),
        'MD plot'),
      x = 'Average LogCPM estimates',y = 'Log Fold-change') +
    theme(
      # plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )
  #ggsave(file.path(rout, paste0(gsub(pattern = ' ', replacement = '_', x = x), '_MD_plot_FDR_PV.png')),
   # plot = p, width = 10, height = 8, device = 'png', dpi = 300)
   ggsave(file.path(rout, paste0(gsub(pattern = ' ', replacement = '_', x = x), '_MD_plot_FDR.png')),
     plot = p, width = 10, height = 8, device = 'png', dpi = 300)
}))


matrix_DEG = DEG[,c("gene_ENS","gl")]
matrix_DEG = merge(matrix_DEG,DEG_melina[,c("gene_ENS","gl")], by.x = "gene_ENS", by.y = "gene_ENS")
DEG_RLE$gl = as.character(DEG_RLE$gl)
matrix_DEG = merge(matrix_DEG,DEG_RLE[,c("gene_ENS","gl")], by.x = "gene_ENS", by.y = "gene_ENS")
matrix_DEG[matrix_DEG == "Not regulated"] = 0
matrix_DEG[matrix_DEG == "Up-regulated"] = 1
matrix_DEG[matrix_DEG == "Down-regulated"] = -1
colnames(matrix_DEG) = c("gene","DEG_TMM","DEG_melina","DEG_RLE")
matrix_DEG$DEG_TMM = as.integer(matrix_DEG$DEG_TMM)
matrix_DEG$DEG_melina = as.integer(matrix_DEG$DEG_melina)
matrix_DEG$DEG_RLE = as.integer(matrix_DEG$DEG_RLE)

pheatmap(matrix_DEG[,c(2:4)])


###

RS = read.table("DOT1L_KO_effect_within_RS_cells_volcano_table_FDR.tsv",sep = "\t", h = T)
RS_anno = RS[RS$gl != "Not regulated",c(2:4,14)]
RS_anno$gl = factor(RS_anno$gl, labels = c("DR","UR"))
# RS_anno$type = "gene"
chrom = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", sep = ";", h = T)[,1:2]
chrom$Chromosome = paste0("chr",chrom$Chromosome)
# colnames(chr) = c("chr","end")
# chr$start = 0
# chr$type = "chr"

centro = read.table("~/Documents/Annotations/mouse_centromere_location.tsv",sep = "\t", h = F)[,c(2:4,8)]
colnames(centro) = c("chr","start","end","gl")
sat = read.table("~/Documents/Annotations/mouse_microsat_location.tsv",sep = "\t", h = F)[,c(2:4)]
sat$type = "microsatellite"
colnames(sat) = c("chr","start","end","gl")
sat$chr = factor(sat$chr, levels = chrom$Chromosome)
RS_anno = rbind(RS_anno,centro)
# RS_anno = rbind(RS_anno,sat)

RS_anno = merge(RS_anno, chrom, by.x = "chr", by.y = "Chromosome", all.x = T)

p = ggplot(RS_anno) +
  geom_bar(aes(x = Total_length, y = factor(chr, levels = rev(chrom$Chromosome))),
    stat = "identity", position = position_dodge(),
    fill = "white", color = "grey") +
  geom_point(aes(x = end, y = chr, color = gl), alpha = 0.5, size = 3) +
  scale_color_manual(values = c("blue","red","black","green")) +
  xlab("chromosome size (bp)") + ylab("") + labs(color = "") +
  theme_bw() + theme(
    text = element_text(size=35, angle = 0),
    legend.position = "bottom",
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank())

ggsave('RS_DEG_gene_position_chr.png',  plot = p, width = 10, height = 8, device = 'png', dpi = 600)
