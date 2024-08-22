####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('UpSetR','ggplot2', 'gtable','stringr','RColorBrewer', 'edgeR', 'DESeq2',
  'gridExtra','ggalluvial') #'clusterProfiler'
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                            INITIALIZATION
####################################################################################################

dir = c("../../DE_analysis_3_factors_Kit","../../DE_analysis_3_factors_SC_SCII_RS_spike")
dir = "RNAseq_B1_B2_FC1"

rout = file.path(".")
samples = c("Kit-","Kit+","RS","SCI","SCII")
# samples = c("RS","SCI","SCII")
samples = c("DKO","ZFY1","ZFY2")

####################################################################################################
##                                                 DATA
####################################################################################################


list_files = list.files(dir,
  pattern = "SCII_cells_volcano_table_FDR_PV.tsv", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x, sep = "\t", h = T, stringsAsFactors = F)
})

gene_length = read.table("~/Documents/Annotations/mouse_gene_length.tsv",stringsAsFactors = F,
  sep = "\t", h = F)
annotation = read.table("~/Documents/Annotations/mouse_gene_annotation.tsv",stringsAsFactors = F,
  h = T, sep = '\t')
bed = read.table("../gencode.vM19.annotation.bed",h=F, sep = '\t')
colnames(bed) = c("chr","start","end", "strand","score","gene_ENS")

####################################################################################################
##                                          CORRELATION PLOT
####################################################################################################

# create_correlation_RS_RSspike_plot = function(RS,RS_spike,out) {
#   ## Separate data in three categories to compare them
#   spike_up = RS_spike[which(grepl("Up-regulated",RS_spike$gl)),]
#   rs_up = RS[which(grepl("Up-regulated",RS$gl)),]
#   spike_down = RS_spike[which(grepl("Down-regulated",RS_spike$gl)),]
#   rs_down = RS[which(grepl("Down-regulated",RS$gl)),]
#
#   all_genes = merge(RS,RS_spike,by.x="gene_ENS",by.y="gene_ENS")
#   up_genes = merge(rs_up,spike_up,by.x="gene_ENS",by.y="gene_ENS")
#   down_genes = merge(rs_down,spike_down,by.x="gene_ENS",by.y="gene_ENS")
#
#   correlation = round(cor.test(all_genes$logCPM.x,all_genes$logCPM.y,
#     method = "pearson")[4]$estimate,4)
#   p_value = cor.test(all_genes$logCPM.x,all_genes$logCPM.y, method = "pearson")[3]
#   if (p_value == 0){p_value = "<2.2e-16"}
#   print(paste0("All genes - R² : ",correlation,"; p-value : ",p_value))
#   correlation_up = round(cor.test(up_genes$logCPM.x,up_genes$logCPM.y,
#     method = "pearson")[4]$estimate,4)
#   p_value_up = cor.test(up_genes$logCPM.x,up_genes$logCPM.y, method = "pearson")[3]
#   if (p_value_up == 0){p_value_up = "<2.2e-16"}
#   print(paste0("Up genes - R² : ",correlation_up,"; p-value : ",p_value_up))
#   correlation_down = round(cor.test(down_genes$logCPM.x,down_genes$logCPM.y,
#     method = "pearson")[4]$estimate,4)
#   p_value_down = cor.test(down_genes$logCPM.x,down_genes$logCPM.y, method = "pearson")[3]
#   if (p_value_down == 0){p_value_down = "<2.2e-16"}
#   print(paste0("Down genes - R² : ",correlation_down,"; p-value : ",p_value_down))
#
#
#   p = ggplot(all_genes ,aes(x = logCPM.x,y = logCPM.y,color = "Not-regulated"))+
#     geom_point()+
#     geom_point(data = up_genes, aes(x = logCPM.x,y = logCPM.y, color = "Up-regulated")) +
#     geom_point(data = down_genes, aes(x = logCPM.x,y = logCPM.y, color = "Down-regulated")) +
#     theme(text = element_text(size = 20)) +
#     labs(x = "log(CPM RNAseq)",y = "log(CPM RS Spike-in)",
#          title = "Correlation plot between deregulate common genes")+
#     scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9, 'Set1')[c(2, 9, 1)]),
#       drop = FALSE)
#   ggsave(paste0(out,"/Correlation_plot_deregulate_common_genes.png"), plot = p, width = 10,
#     height = 8, device = 'png', dpi = 150)
# }
#
# create_correlation_RS_RSspike_plot(list_data$RS, list_data$RS_spike,rout)

####################################################################################################
##                                      UPSET GENES DEREGULATED
####################################################################################################

gene = unique(c(list_data$SC$gene_ENS, list_data$RS$gene_ENS, list_data$SCII$gene_ENS,
  list_data$'Kit-'$gene_ENS, list_data$'Kit+'$gene_ENS))
# gene = unique(c(list_data$SC$gene_ENS, list_data$RS$gene_ENS, list_data$SCII$gene_ENS))

gene = unique(c(list_data$DKO$gene_ENS, list_data$ZFY1$gene_ENS, list_data$ZFY2$gene_ENS))

DKO = list_data$DKO[list_data$DKO$FDR <= 0.05,]
ZFY1 = list_data$ZFY1[list_data$ZFY1$FDR <= 0.05,]
ZFY2 = list_data$ZFY2[list_data$ZFY2$FDR <= 0.05,]

SSC_m = list_data$'Kit-'[list_data$'Kit-'$FDR <= 0.05,]
SSC_p = list_data$'Kit+'[list_data$'Kit+'$FDR <= 0.05,]

# df = rbind(data.frame(table(SSC_m$gl)),data.frame(table(SSC_p$gl)),data.frame(table(SCI$gl)),
#   data.frame(table(SCII$gl)),data.frame(table(RS$gl)))

df = rbind(data.frame(table(SSC_p$gl)),data.frame(table(SCI$gl)),
  data.frame(table(SCII$gl)),data.frame(table(RS$gl)))
df2 = df[2:4,]
df2$Freq = 0
df = rbind(df2,df)
df$name = rep(samples[c(1:2,4,5,3)],each = 2)
# df = rbind(data.frame(table(SC$gl)),data.frame(table(SCII$gl)),data.frame(table(RS$gl)))
# df$name = rep(samples[c(2,3,1)],each = 2)

## Distribution of DEG genes (barplot)
p = ggplot(df, aes(x = rev(Var1), y = Freq, alpha = factor(name, levels = samples[c(1:2,4,5,3)]))) +
  geom_bar(stat = "identity", position = "dodge", fill = rep(c("blue","red"),5)) +
  xlab("") + ylab("number of genes") + labs(fill = "", alpha = "") +
  scale_x_discrete(labels = c("Down-regulated","Up-regulated")) +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position = "bottom")
  # scale_fill_manual(values = c("blue","red"))
ggsave(filename = "DEG_distribution_barplot.png", plot = p,width = 10, height = 8,
  device = 'png', dpi = 150)

p2 = ggplot(df, aes(x = factor(name, levels = samples[c(1:2,4,5,3)]), y = Freq,fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") + ylab("number of genes") + labs(fill = "", alpha = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position = "bottom") +
  scale_fill_manual(values = c("blue","red"))
ggsave(filename = "DEG_distribution_barplot_V2.png", plot = p2,width = 10, height = 8,
  device = 'png', dpi = 150)

## Upset plot
RS_down = RS[RS$gl == "Down-regulated","gene_ENS"]
RS_up = RS[RS$gl == "Up-regulated","gene_ENS"]
SC_down = SC[SC$gl == "Down-regulated","gene_ENS"]
SC_up = SC[SC$gl == "Up-regulated","gene_ENS"]
SCII_down = SCII[SCII$gl == "Down-regulated","gene_ENS"]
SCII_up = SCII[SCII$gl == "Up-regulated","gene_ENS"]
SSC_p_down = SSC_p[SSC_p$gl == "Down-regulated","gene_ENS"]
SSC_p_up = SSC_p[SSC_p$gl == "Up-regulated","gene_ENS"]
SSC_m_down = SSC_m[SSC_m$gl == "Down-regulated","gene_ENS"]
SSC_m_up = SSC_m[SSC_m$gl == "Up-regulated","gene_ENS"]

up_down_table = data.frame("gene" = gene, "SSC_m_up" = 0, "SSC_m_down" = 0,"SSC_p_up" = 0,
  "SSC_p_down" = 0, "SC_up" = 0, "SC_down" = 0, "SCII_up" = 0, "SCII_down" = 0, "RS_up" = 0,
  "RS_down" = 0)
# up_down_table = data.frame("gene" = gene,"SC_up" = 0, "SC_down" = 0, "SCII_up" = 0, "SCII_down" = 0,
#   "RS_up" = 0,"RS_down" = 0)

up_down_list = list(SSC_m_up,SSC_m_down,SSC_p_up,SSC_p_down,SC_up,SC_down,SCII_up,SCII_down,
  RS_up,RS_down)
# up_down_list = list(SC_up,SC_down,SCII_up,SCII_down,RS_up,RS_down)

names(up_down_list) = colnames(up_down_table)[-1]

for(type in colnames(up_down_table)[-1]){
  print(type)
  pos = unlist(lapply(up_down_list[[type]], function(x) {
    which(up_down_table$gene == x)
  }))
  up_down_table[pos,type] = 1
}

# png("SC_SCII_RS_upset.png",width = 1500,height = 800)
# upset(up_down_table, sets = c("SC_up","SC_down", "SCII_up", "SCII_down", "RS_up", "RS_down"),
#   order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
#   sets.x.label = "Number gene\nper celltype",
#   sets.bar.color = rep(c("red","blue"), 3), text.scale = c(3,3,3,3,4,5),
#   number.angles = 0,
#   point.size = 6, line.size = 2)
# dev.off()

colnames(up_down_table) = c("gene","Kit- up","Kit- down","Kit+ up","Kit+ down","SCI up","SCI down",
  "SCII up","SCII down","RS up","RS down")

png("Kit_SC_SCII_RS_upset.png",width = 1500,height = 800)
upset(up_down_table, sets = c("Kit- up","Kit- down","Kit+ up","Kit+ down","SCI up","SCI down",
  "SCII up","SCII down","RS up","RS down"),
  order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
  sets.x.label = "Number gene\nper celltype",
  sets.bar.color = rep(c("red","blue"), 5),
  text.scale = c(3,3,3,3,4,2),
  number.angles = 15,
  point.size = 6, line.size = 2)
dev.off()

DKO_deg = DKO[DKO$gl != "Not regulated","gene_ENS"]
ZFY1_deg = ZFY1[ZFY1$gl != "Not regulated","gene_ENS"]
ZFY2_deg = ZFY2[ZFY2$gl != "Not regulated","gene_ENS"]
SSC_m_deg = SSC_m[SSC_m$gl != "Not regulated","gene_ENS"]
SSC_p_deg = SSC_p[SSC_p$gl != "Not regulated","gene_ENS"]

RS_FDR = RS[RS$FDR <= 0.05,"gene_ENS"]
SC_FDR = SC[SC$FDR <= 0.05,"gene_ENS"]
SCII_FDR = SCII[SCII$FDR <= 0.05,"gene_ENS"]
RS_PV = RS[RS$PValue <= 0.05,"gene_ENS"]
SC_PV = SC[SC$PValue <= 0.05,"gene_ENS"]
SCII_PV = SCII[SCII$PValue <= 0.05,"gene_ENS"]


deg_table = data.frame("gene" = gene,"SSC_m" = 0,"SSC_p" = 0,"SC" = 0, "SCII" = 0, "RS" = 0)
# deg_table = data.frame("gene" = gene,"SC" = 0, "SCII" = 0, "RS" = 0)
deg_table = data.frame("gene" = gene,"SC_PV" = 0, "SCII_PV" = 0, "RS_PV" = 0,
  "SC_FDR" = 0, "SCII_FDR" = 0, "RS_FDR" = 0)


deg_table = data.frame("gene" = gene,"DKO" = 0, "ZFY1" = 0, "ZFY2" = 0)

deg_list = list(SSC_m_deg,SSC_p_deg,SC_deg,SCII_deg,RS_deg)
# deg_list = list(SC_deg,SCII_deg,RS_deg)
deg_list = list(SC_PV,SCII_PV,RS_PV,SC_FDR,SCII_FDR,RS_FDR)



deg_list = list(DKO_deg,ZFY1_deg,ZFY2_deg)

names(deg_list) = colnames(deg_table)[-1]

for(type in colnames(deg_table)[-1]){
  print(type)
  pos = unlist(lapply(deg_list[[type]], function(x) {
    which(deg_table$gene == x)
  }))
  deg_table[pos,type] = 1
}

png("RS_upset_FDR.png",width = 1200,height = 800)
upset(deg_table, sets = c("DKO","ZFY1","ZFY2"),
  order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
  sets.x.label = "Number gene\nper celltype",
  text.scale = c(3,3,3,3,4,5),
  number.angles = 0,
  set_size.scale_max = 300,
  point.size = 6, line.size = 2)
dev.off()

colnames(deg_table) = c("gene","Kit-","Kit+","SCI","SCII","RS")

png("Kit_SC_SCII_RS_upset_reduce.png",width = 1600,height = 1080)
upset(deg_table, sets = c("Kit-","Kit+","SCI","SCII","RS"),
  order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
  sets.x.label = "Number gene\nper celltype",
  text.scale = c(3,3,3,3,4,5),
  number.angles = 0,
  point.size = 6, line.size = 2)
dev.off()

####################################################################################################
##                                    COMPARISON KIT-/+ WITH okada
####################################################################################################

SSC_okada_deg = list_data$SSC[list_data$SSC$gl != "Not regulated","gene_ENS"]
Kit_p_deg = list_data[["Kit+"]][list_data[["Kit+"]]$gl != "Not regulated","gene_ENS"]
Kit_m_deg = list_data[["Kit-"]][list_data[["Kit-"]]$gl != "Not regulated","gene_ENS"]

venn.diagram(x = list(Kit_m_deg,Kit_p_deg,SSC_okada_deg),
  category.names = c("Kit-","Kit+","SSC_okada"),
  filename = "SSC_DEG_comparison.png",
  output = T,
  imagetype = "png",
  lwd = 2,
  col = c("#3498DB","#884EA0","#2ECC71"),
  fill = c(alpha("#3498DB",0.3),alpha("#884EA0",0.3),alpha("#2ECC71",0.3)),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.pos = c(0,0,180),
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.055),
  cat.col = c("#3498DB","#884EA0","#2ECC71")
)

####################################################################################################
##                                    CLUSTERPROFILER ON DEG
####################################################################################################

library('org.Mm.eg.db')

create_enrichGO_plot = function(df,title) {
  df$ENTREZID = mapIds(org.Mm.eg.db, df$gene_name, 'ENTREZID','SYMBOL')
  df_without_na = na.omit(df$ENTREZID)
  go = enrichGO(df$ENTREZID,org.Mm.eg.db, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff  = 0.2,
    minGSSize = 1,universe = data$ENTREZID, pAdjustMethod = 'none')

  p = dotplot(go)
  ggsave(paste0(title,".png"), plot = p, device = "png", dpi = 300, height = 10, width = 20)

  write.table(go@result,paste0(title,".tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
}

lapply(names(list_data), function(name) {

  data = list_data[[name]]
  data$ENTREZID = mapIds(org.Mm.eg.db, data$gene_name, 'ENTREZID','SYMBOL')
  DEG = data[data$gl != "Not regulated",]
  up = data[data$gl == "Up-regulated",]
  down = data[data$gl == "Down-regulated",]

  create_enrichGO_plot(DEG,paste0(name,'_DEG_enrichGO'))
  create_enrichGO_plot(up,paste0(name,'_up_enrichGO'))
  create_enrichGO_plot(down,paste0(name,'_down_enrichGO'))
})

####################################################################################################
##                                    KIT PAPER FIGURES
####################################################################################################


Kit = read.table("DOT1L_KO_effect_within_RS_cells_volcano_table_FDR.tsv", h = T, sep = "\t")
Kit$name = ""
gene_interet = c(Kit[order(Kit$PValue, decreasing = F)[1:10],"gene_name"],"Casp8","Aldh1a1",
  "Aldh3b3","Aldh1a7")
l = unlist(lapply(gene_interet,function(gene) {
  print(gene)
  print(which(Kit$gene_name == gene))
}))
Kit[l,"name"] = gene_interet

p = ggplot(data = Kit,aes(x = logFC, y = -log10(PValue), color = gl,label = name)) +
    geom_point(size = 5, alpha = .5) +
    geom_text_repel(color = "black",nudge_x = 0.25,nudge_y = 0.1, max.overlaps = 20) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), alpha = .5,
      color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), alpha = .5, color = 'grey50',
      linetype = 'dashed') +
    scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
    labs(subtitle = paste("DOT1L KO effect within Kit- cells", paste0('(*n=', sum(Kit$gl != 'NR'),')')),
      title = paste(paste0('FC>1.5'),paste0('& P < 5%'),
        'V plot'),
      x = 'Log Fold-change',y = '-Log10(P)') +
    theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )
ggsave("DOT1L_KO_effect_within_Kit-_cells_volcano_plot_specific_genes.png",plot = p, device = "png",
  dpi = 300, height = 10, width = 10)

genecounts = read.table("RS_genecounts.gct", sep = "\t", h = T)
x = merge(genecounts,Kit, by.x = "Name", by.y = "gene_ENS")[c(1,3:12,25)]
x_DEG = x[x$gl != "Not regulated",]
x_DEG = x_DEG[order(x_DEG$gl, decreasing = T),]
rownames(x_DEG) = x_DEG$Name

anno_row = data.frame("Gene expression" = x_DEG$gl, row.names = x_DEG$Name)
anno_col = data.frame(Phenotype = rep(c("CTL","KO"),5), row.names = colnames(x_DEG[,2:11]))
anno_color = list("Gene.expression" = c("Down-regulated (FDR < 5%)" = "blue",
  "Up-regulated (FDR < 5%)" = "red")
  # Phenotype = c(CTL = "grey", KO = "black"))

pheatmap(x_DEG[,c(2,4,6,8,10,3,5,7,9,11)], scale = "row", cluster_rows = F, cluster_cols = F,
  annotation_row = anno_row, annotation_colors = anno_color,
   # annotation_col = anno_col,
  show_rownames = F, labels_col = c("CTL1","CTL2","CTL3","CTL4","CTL5","KO1","KO2","KO3","KO4","KO5"),
  filename = "RS_DEG_expression_heatmap.png", width = 11, height = 11, fontsize = 20)


####################################################################################################
##                                      ALLUVIAL
####################################################################################################

data_gene = unique(rbind(list_data$SCI[,1:5],list_data$SCII[,1:5],list_data$RS[,1:5])),
    # list_data$'Kit-'[,1:5],list_data$'Kit+'[,1:5]))
rownames(data_gene) = data_gene$gene_ENS

aluvial_plot = function(matrix, fill_color = "") {
  plot = ggplot(matrix,
    aes(y = freq, axis1 = factor(SCI), axis2 = factor(SCII), axis3 = factor(RS),
      fill = factor(motif), color = "1")) +
    # aes(y = freq, axis1 = factor(`Kit-`), axis2 = factor(`Kit+`), axis3 = factor(SCI),
    #   axis4 = factor(SCII), axis5 = factor(RS), fill = factor(motif), color = "1")) +
    geom_alluvium() +
    geom_stratum(width = 1/12, fill = fill_color, color = fill_color) +
    # scale_x_discrete(limits = c("Kit-","Kit+","SCI","SCII","RS"), expand = c(.1, .1)) +
    scale_x_discrete(limits = c("SCI","SCII","RS"), expand = c(.1, .1)) +
    scale_fill_grey(start = 1, end = 0) +
    scale_color_manual(values = "black") +
    ylab("") +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=25, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none")
    return(plot)
}

comparison_matrix = expand.grid('Kit-' = c("Up-regulated","Down-regulated","Not regulated"),
  'Kit+' = c("Up-regulated","Down-regulated","Not regulated"),
  SCI = c("Up-regulated","Down-regulated","Not regulated"),
  SCII = c("Up-regulated","Down-regulated","Not regulated"),
  RS = c("Up-regulated","Down-regulated","Not regulated")) #243 condition
comparison_matrix = expand.grid(
  SCI = c("Up-regulated","Down-regulated","Not regulated"),
  SCII = c("Up-regulated","Down-regulated","Not regulated"),
  RS = c("Up-regulated","Down-regulated","Not regulated"))

comparison_matrix$motif = paste0(
  # comparison_matrix$`Kit-`,
  # comparison_matrix$`Kit+`,
  comparison_matrix$SCI,
  comparison_matrix$SCII,
  comparison_matrix$RS)
rownames(comparison_matrix) = comparison_matrix$motif
comparison_matrix$freq = 0

gene_table = do.call("cbind",lapply(names(list_data)[3:5], function(x) {
  data = list_data[[x]]
  # nrow(data)
  rownames(data) = data$gene_ENS
  data_gene$gl = "Not regulated"
  data_gene[data$gene_ENS,"gl"] = data$gl
  data_gene$name = x
  return(data_gene)
}))

gene_table = gene_table[,colnames(gene_table) == "gl"]
colnames(gene_table) = names(list_data)[3:5]
gene_table = cbind(data_gene,gene_table)

gene_table$motif = paste0(
  # gene_table$`Kit-`,
  # gene_table$`Kit+`,
  gene_table$SCI,
  gene_table$SCII,
  gene_table$RS)

write.table(gene_table, "SC_SCII_RS_DEG_distribution_FDR_spike.tsv", col.names = T, row.names = F,
  sep  ="\t", quote = F)

frequence_motif = table(gene_table$motif)
comparison_matrix[names(frequence_motif),"freq"] = frequence_motif

## Remove empty motif
comparison_matrix = comparison_matrix[comparison_matrix$freq > 10 ,]
## Remove les not
comparison_matrix = comparison_matrix[-nrow(comparison_matrix),]

## Alluvial with all motif
# allu = aluvial_plot(comparison_matrix, fill_color = c(rep(c("grey","red"), 1),
#   rep(c("blue","grey","red"),4)))
allu = aluvial_plot(comparison_matrix, fill_color = rep(c("grey","blue","red"),3))

ggsave(filename = "SC_SCII_RS_alluvial_plot_FDR.png", plot = allu, width = 15, height = 5,
  device = 'png', dpi = 450)

## Alluvial of DEG dynamic change
# matrix_change_dynamic = comparison_matrix[c(1,3,5,9,11),]
# gene_table_change = do.call("rbind",lapply(matrix_change_dynamic$motif, function(motif) {
#   gene_table[gene_table$motif == motif,]
# }))
# gene_table_change$gene = rownames(gene_table_change)
# data = list_data$RS
# df = merge(data,gene_table_change, by.x = "gene_ENS", by.y = "gene")[,c(1:3,10:15)]
# bed = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed", sep = "\t", h = F)
# colnames(bed) = c("chr","start","end","strand","score","gene_ENS")
# df = merge(bed, df, by.x = "gene_ENS", by.y = "gene_ENS")
#
# write.table(df, "Kit_SC_SCII_RS_CHANGE_table.tsv", sep = "\t", row.names = F, col.names = F,
#   quote = F)
# change = aluvial_plot(matrix_change_dynamic, fill_color = c("blue","red","grey","blue","red","grey",
#   "blue","red","blue","red","blue","red"))
# ggsave(filename = "Kit_SC_SCII_RS_CHANGE_alluvial_plot.png", plot = change, width = 15, height = 5,
#   device = 'png', dpi = 450)

## Remove DEG change dynamic
# comparison_matrix = comparison_matrix[-c(1,3,5,9,11),]
## Alluvial of DOWN
comparison_down = comparison_matrix[grepl("Down-regulated",rownames(comparison_matrix)),]

down = aluvial_plot(comparison_down, fill_color = c(rep(c("grey"),2),
  rep(c("blue","grey"),3)))
ggsave(filename = "Kit_SC_SCII_RS_DOWN_alluvial_plot.png", plot = down, width = 15, height = 5,
  device = 'png', dpi = 450)

## Alluvial of UP
comparison_up = comparison_matrix[grepl("Up-regulated",rownames(comparison_matrix)),]

up = aluvial_plot(comparison_up, fill_color = rep(c("grey","red"),5))
ggsave(filename = "Kit_SC_SCII_RS_UP_alluvial_plot.png", plot = up, width = 15, height = 5,
  device = 'png', dpi = 450)


####################################################################################################
##                                      COMPARE CLUSTER
####################################################################################################

dir = c("SC","SCII","RS")
lapply(dir, function(d) {
  DEG = read.table(paste0(d,"_spike_volcano_table_FDR.tsv"), sep = "\t", h = T)
  list_table = lapply(unique(DEG$gl), function(x) {
    return(DEG[DEG$gl == x,"gene_name"])
  })
  names(list_table) = c("Not","Up","Down")
  compareGO = compareCluster(geneCluster = list_table, fun = "enrichGO",
    Org = "org.Mm.eg.db",ont = "BP", pAdjustMethod = "BH", keyType = "SYMBOL", qvalueCutoff = 0.05,
    pvalueCutoff = 1)
  print(d)
  png(paste0(d,"_GO_dotplot.png"), height = 400, width = 800)
  dotplot(compareGO, showCategory = 10)
  dev.off()
})
