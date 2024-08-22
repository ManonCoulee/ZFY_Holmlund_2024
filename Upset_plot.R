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
                 
dir = "RNAseq_B1_B2_FC1"

rout = file.path(".")
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
##                                      UPSET GENES DEREGULATED
####################################################################################################

gene = unique(c(list_data$DKO$gene_ENS, list_data$ZFY1$gene_ENS, list_data$ZFY2$gene_ENS))

DKO = list_data$DKO[list_data$DKO$FDR <= 0.05,]
ZFY1 = list_data$ZFY1[list_data$ZFY1$FDR <= 0.05,]
ZFY2 = list_data$ZFY2[list_data$ZFY2$FDR <= 0.05,]

DKO_deg = DKO[DKO$gl != "Not regulated","gene_ENS"]
ZFY1_deg = ZFY1[ZFY1$gl != "Not regulated","gene_ENS"]
ZFY2_deg = ZFY2[ZFY2$gl != "Not regulated","gene_ENS"]

deg_table = data.frame("gene" = gene,"DKO" = 0, "ZFY1" = 0, "ZFY2" = 0)
deg_list = list(DKO_deg,ZFY1_deg,ZFY2_deg)
names(deg_list) = colnames(deg_table)[-1]

for(type in colnames(deg_table)[-1]){
  print(type)
  pos = unlist(lapply(deg_list[[type]], function(x) {
    which(deg_table$gene == x)
  }))
  deg_table[pos,type] = 1
}

png("SCII_upset_FDR.png",width = 1200,height = 800)
upset(deg_table, sets = c("DKO","ZFY1","ZFY2"),
  order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
  sets.x.label = "Number gene\nper celltype",
  text.scale = c(3,3,3,3,4,5),
  number.angles = 0,
  set_size.scale_max = 300,
  point.size = 6, line.size = 2)
dev.off()
