##Monocle 2
library(monocle)
library(RColorBrewer)
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
mytheme <- theme_bw() +
  theme(#legend.position = "right",
    axis.text.x=element_text(angle=0,hjust=0.5),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    strip.text = element_text(size=12,face="bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'))

#Set your Path
Path='/home/chengww/data/project/database'
setwd(Path)

tle4_sub <- readRDS('./scrna-atac/tle4_rna_atac.rds')
DefaultAssay(tle4_sub)<-'RNA'
set.seed(1995)
subobj <- tle4_sub
##data for monocle
RA_matrix <- as(as.matrix(subobj@assays$RNA@counts),'sparseMatrix')
p_data <-subobj@meta.data
f_data<-data.frame(row.names=rownames(subobj),gene_short_name=rownames(subobj))

##CDS object
pd<-new("AnnotatedDataFrame", data =p_data)
fd<-new("AnnotatedDataFrame", data =f_data)
cds <- newCellDataSet(RA_matrix,
                      phenoData =pd,
                      featureData =fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily=negbinomial.size())

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)

expressed_genes <- VariableFeatures(subobj) ##gene set1 
cds <- setOrderingFilter(cds, expressed_genes)  
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

p1<- plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.3,show_backbone=TRUE)+
  gg.theme 
p2 <- plot_cell_trajectory(cds,color_by="State", cell_size =0.3,show_backbone=TRUE)+
  scale_color_manual(values = type.col[c(1,5,10)]) + 
  theme(legend.position = "top")
p3 <- plot_cell_trajectory(cds, color_by = "sub",
                            cell_size =0.3,show_backbone=TRUE)+
  scale_color_manual(values = c('#F87A71','#00BEC3')) + 
  theme(legend.position = "top")+
  #facet_wrap(~celltype_new, nrow = 3)+
  gg.theme
p1/p3
p2
# plot_cell_trajectory(cds, color_by = "sub",cell_size =0.3,show_backbone=TRUE)+
#   facet_wrap(~sub, nrow = 1)+gg.theme
saveRDS(cds,'./cds_ex.Rds')





