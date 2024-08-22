####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2','stringr','RColorBrewer','gridExtra','ggalluvial', 'ggpubr')
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
names(list_data) = samples

bed = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed",h=F, sep = '\t')
colnames(bed) = c("chr","start","end", "strand","score","gene_ENS")
                 
####################################################################################################
##                              Gene expression distribution
####################################################################################################

create_genes_deregulate_FC_plot = function(name,out, listData) {
  data = listData[[name]]
  data$chromosome = "autosome"
  data[data$chr == "chrX","chromosome"] = "X"
  data[data$chr == "chrY","chromosome"] = "Y"
  data$cellType = name
  return(data)
}

listFC = lapply(samples, create_genes_deregulate_FC_plot, listData = list_data)
tableFC = do.call('rbind',listFC)
tableFC$cellType = factor(tableFC$cellType, levels = c("DKO","ZFY1","ZFY2"))
tableFC$logFC = as.integer(tableFC$logFC)
my_comparisons <- list(c("autosome", "X"),c("autosome", "Y"),c("X", "Y"))

p = ggplot(tableFC, aes(x = chromosome, y = logFC, fill = chromosome)) +
  geom_boxplot() +
  facet_grid(~cellType, switch = "both") +
  labs(fill = "") + xlab("") + ylab("log2FC KOvsCTL") + 
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = -1.5, linetype = "dashed") +
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = c("DodgerBlue3","orange","salmon")) +
  theme_bw() + theme(
    text = element_text(size=30, angle = 0),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(colour="white",fill="white",size=1.5, linetype="solid")) +
    stat_compare_means(aes(label=..p.adj..), label = "p.signif",method = "wilcox.test",
      comparisons = my_comparisons)

ggsave(filename = "SCII_FC.png", plot = p, width = 6, height = 5, device = 'png',
  dpi = 600)
