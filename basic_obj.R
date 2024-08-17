library(ggplot2)
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
library(RColorBrewer)
type.col <- colorRampPalette(c(brewer.pal(n=3,name = 'Set3'),
                               brewer.pal(n=8,name='Dark2')))(19)
colorRampPalette()(15)


big.marker <- c('Cd34','Avp','Prss57',#造血干细胞 HSC
                'Klf1','Gypa','Slc4a1',#'CD36','TSPO2',#EPC
                'Hba1', 'Hbb', #Erythrocytes
                'Cxcl12','Col6a1','Col1a1','Col3a1',#Stromal cell
                'Cd3d','Cd3e','Cd4','Cd40lg','Sell',#naive
                'Ccr7','Lef1',#naive
                'Cd8a','Cd8b',#CD8+ T cells
                'Gzmk',#'GZMH',#CD8 mem
                'Gnly','Gzmb',#CD8 E
                'Klrf1','Klrc1','Ncam1',#NK cell#gnly'KLRB1',
                'Znf683',#NKT
                'Vpreb1','Tcl1a','Cd19','Cd79a','Ms4a1','Tcl1a',#B cells
                'Lyz', #Myeloid cells
                'Cd14', #classic Monocytes, 
                #'SEMA6B',
                'Fcgr3a', #non classic Monocytes
                'Retn','Mki67',#prMono
                #'ELANE',#gmp
                'Clec9a',##cDC1 
                "Clec10a",'Fcer1a','Cd1c',##cDC2
                "Il3ra","Lilra4",#pDC
                'Sdc1','Ighg1','Mzb1'#Plasma cells
)
big.marker <- c('CD34','AVP','PRSS57',#造血干细胞 HSC
                'KLF1','GYPA','SLC4A1',#'CD36','TSPO2',#EPC
                'HBA1', 'HBB', #Erythrocytes
                'CXCL12','COL6A1','COL1A1','COL3A1',#Stromal cell
                'CD3D','CD3E','CD4','CD40LG','SELL',#naive
                'CCR7','LEF1',#naive
                'IL7R','TIMP1',#'S100A4','GPR183',#c MEMORY #CCR7low
                'CXCR3','GZMK','GZMA',#CD4 TEM
                'CD8A','CD8B',#CD8+ T cells
                'GZMK',#'GZMH',#CD8 mem
                'GNLY','GZMB',#CD8 E
                'KLRF1','KLRC1','NCAM1',#NK cell#gnly'KLRB1',
                'ZNF683',#NKT
                'RORC','ZBTB16','SLC4A10','KLRB1',#'TRAV1',#',#MAIT 
                'IL2RA','FOXP3','CTLA4',#Treg
                'VPREB1','TCL1A','CD19','CD79A','MS4A1',#B cells
                'TNFRSF13B',#B memory
                #'SEC11C','XBP1',#MEMORY B
                'LYZ', #Myeloid cells
                'CD14', #classic Monocytes, 
                #'SEMA6B',
                'FCGR3A', #non classic Monocytes
                'RETN','MKI67',#prMono
                #'ELANE',#gmp
                'CLEC9A',##cDC1 
                "CLEC10A",'FCER1A','CD1C',##cDC2
                "IL3RA","LILRA4",#pDC
                'SDC1','IGHG1','MZB1'#Plasma cells
)
