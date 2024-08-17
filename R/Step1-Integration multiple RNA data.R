##Load Packages
library(RColorBrewer)
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(EnsDb.Mmusculus.v79)

#Set your Path
Path='/home/chengww/data/project/database'
setwd(Path)
type.col <- colorRampPalette(c( brewer.pal(n=12,name='Paired'),
                                brewer.pal(n=12,name = 'Set3')))(39)


##The same level of cellular annotation was performed prior to integration--------------------------------------------
#Load data
#1
ISSAAC_mCortex <-readRDS('./scrna-atac/E-MTAB-11264/ISSAAC_mCortex_allobj.rds')
DefaultAssay(ISSAAC_mCortex)<-'RNA'
ISSAAC_mCortex$orig.ident <-'ISSAAC_mCortex'
big.marker <- c(
  'Rbfox3',
  'Slc32a1',
  'Slc17a7',
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act
  'Rorb',#Ex-L4/5 IT
  'Slc17a8',#Ex-L5 NP
  'Cxcl14',#Ex-L5 NP Cxcl14
  'Npr3',#Ex-L5-PT
  'Foxp2','Tle4',#Ex-L6 CT
  'Bmp3',#Ex-L6 IT Bmp3
  'Gfra1',#Ex-L6 IT Gfra1
  'Ccn2',#Ex-L6b
  'Ndst4',#Ex-PIR Ndst4
  'Tac1',#Misc
  'Drd2',#In-Drd2
  'Hap1','Scn5a',#In-Hap1
  'Cemip',#In-Pvalb
  'Sst','Chodl','Elfn1',#In-Sst
  'Tac1',#In-Tac1
  'Vip',#In-Vip/Lamp5
  'Nptx2',
  'Aqp4',#Astro
  'Pdgfra','Olig1',#OPC
  'Mal',#Oligo
  'Col1a2','Vtn',#VLMC
  'Cldn5','Mcam',#ENdo
  'Mrc1','Cd163',#macrophage
  'Kcnj8',#peri, pericyte
  'Acta2'#SMC, smooth muscle cell
  
)
DefaultAssay(ISSAAC_mCortex) <-'RNA'
DotPlot(ISSAAC_mCortex,features = big.marker %>% unique(),
        cluster.idents = F,group.by = 'celltype')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
ISSAAC_mCortex$celltype <- as.character(ISSAAC_mCortex$celltype)
ISSAAC_mCortex$celltype <- ifelse(ISSAAC_mCortex$celltype %in% c("In-Tac1"),"In-Drd1 Tac1",ISSAAC_mCortex$celltype)
ISSAAC_mCortex$celltype <- ifelse(ISSAAC_mCortex$celltype %in% c("Ex-L6 IT Oprk1"),"Ex-L6 IT Gfra1",ISSAAC_mCortex$celltype)

level = c("Ex-L2/3 IT","Ex-L2/3 IT Act","Ex-L4/5 IT",
          "Ex-L5 NP", "Ex-L5 NP Cxcl14","Ex-L5-PT", "Ex-L6 CT",
          "Ex-L6 IT Bmp3","Ex-L6 IT Gfra1","Ex-L6b","Ex-PIR Ndst4",
          "Misc","In-Drd1 Tac1","In-Drd2","In-Hap1" ,"In-Pvalb", "In-Sst",
          "In-Vip/Lamp5","Astro","OPC","Oligo","VLMC")
ISSAAC_mCortex$celltype <- factor(ISSAAC_mCortex$celltype,levels = level)
DimPlot(ISSAAC_mCortex,reduction = 'umap.rna',group.by  = 'celltype',label = T,cols = type.col)
DimPlot(ISSAAC_mCortex,reduction = 'wnn.umap',group.by  = 'celltype',label = T,repel = T,cols = type.col)+
  labs(title = 'ATAC celltype')+NoLegend()

#2
GSE126074_AdBrainCortex <- readRDS('./scrna-atac/GSE126074_AdBrainCortex/GSE126074_AdBrainCortex_allobj.rds')
DefaultAssay(GSE126074_AdBrainCortex)<-'ATAC'
peaks <- CallPeaks(GSE126074_AdBrainCortex,macs2.path = '/home/chengww/anaconda3/envs/MACS/bin/macs2')
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(GSE126074_AdBrainCortex),
  features = peaks,
  cells = colnames(GSE126074_AdBrainCortex)
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# create a new assay using the MACS2 peak set and add it to the Seurat object
GSE126074_AdBrainCortex[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments ="./scrna-atac/GSE126074_AdBrainCortex/fragments.sort.bed.gz",
  annotation = annotations
)
saveRDS(GSE126074_AdBrainCortex,file = './scrna-atac/GSE126074_AdBrainCortex/GSE126074_AdBrainCortex_peaks_allobj.rds')

DefaultAssay(GSE126074_AdBrainCortex)<-'RNA'
GSE126074_AdBrainCortex$orig.ident<-'GSE126074_AdBrainCortex'
big.marker <- c(
  'Rbfox3',
  'Slc32a1',
  'Slc17a7',
  'Rasgrf2',
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act
  'Rorb',#Ex-L4/5 IT
  'Rmst',
  'Thsd7a',
  'Il1rapl2',
  'Slc17a8',#Ex-L5 NP
  'Cxcl14',#Ex-L5 NP Cxcl14
  'Npr3',#Ex-L5-PT
  'Foxp2',#Ex-L6 CT
  'Bmp3',#Ex-L6 IT Bmp3
  'Gfra1',#Ex-L6 IT Gfra1
  'Ccn2',#Ex-L6b
  'Ndst4',#Ex-PIR Ndst4
  #Misc
  'Drd2',#In-Drd2
  'Hap1','Scn5a',#In-Hap1
  'Cemip',#In-Pvalb
  'Sst','Chodl','Elfn1',#In-Sst
  'Tac1',#In-Tac1
  'Vip',#In-Vip/Lamp5
  #'Lamp5',
  'Aqp4',#Astro
  'Pdgfra','Olig1',#OPC
  'Mal',#Oligo
  #'Col1a2','Vtn'#VLMC
  'Cldn5','Mcam'
)
DotPlot(GSE126074_AdBrainCortex,features = big.marker %>% unique(),
        cluster.idents = F,group.by = 'celltype')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
DimPlot(GSE126074_AdBrainCortex,group.by  = 'seurat_clusters',label = T,cols = type.col)
Idents(GSE126074_AdBrainCortex) <-'seurat_clusters'
GSE126074_AdBrainCortex <- RenameIdents(GSE126074_AdBrainCortex, 
                                        '0' = 'Ex-L2/3 IT',
                                        '1' = 'Ex-L4/5 IT','2' = 'Ex-L4/5 Il1rapl2','3' = 'Ex-L6 CT',
                                        '4' = 'Ex-L4/5 Thsd7a','5' = 'Ex-L5-PT','6' = 'In-Pvalb',
                                        '7' = 'In-Sst','8' = 'Astro','9' = 'In-Vip',
                                        '10' = 'Oligo','11' = 'Ex-L6 IT Gfra1','12' = 'Ex-L6 IT Bmp3',
                                        '13' = 'Ex-L5 NP'
)

GSE126074_AdBrainCortex$celltype <- Idents(GSE126074_AdBrainCortex)
DimPlot(GSE126074_AdBrainCortex,group.by  = 'celltype',label = T,cols = type.col)
level = c("Ex-L2/3 IT","Ex-L4/5 IT", 'Ex-L4/5 Il1rapl2','Ex-L4/5 Thsd7a','Ex-L5 NP','Ex-L5-PT',
          'Ex-L6 CT','Ex-L6 IT Gfra1','Ex-L6 IT Bmp3',
          'In-Pvalb', "In-Sst", "In-Vip","Astro","Oligo")
GSE126074_AdBrainCortex$celltype <- factor(GSE126074_AdBrainCortex$celltype,levels = level)
DimPlot(GSE126074_AdBrainCortex,reduction = 'umap.rna',group.by  = 'celltype',label = T,cols = type.col)
DimPlot(GSE126074_AdBrainCortex,reduction = 'wnn.umap',group.by  = 'celltype',label = T,repel = T,cols = type.col)+
  labs(title = 'ATAC celltype')+NoLegend()


#3
mouse_brain<- readRDS('./scrna-atac/mouse_brain_10x/mouse_brain_10x_allobj.rds')
DefaultAssay(mouse_brain)<-'ATAC'
peaks <- CallPeaks(mouse_brain,macs2.path = '/home/chengww/anaconda3/envs/MACS/bin/macs2')
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(mouse_brain),
  features = peaks,
  cells = colnames(mouse_brain)
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# create a new assay using the MACS2 peak set and add it to the Seurat object
mouse_brain[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments ="./scrna-atac/mouse_brain_10x/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
  annotation = annotations
)
saveRDS(mouse_brain,file = './scrna-atac/mouse_brain_10x/mouse_brain_10x_peaks_allobj.rds')

DefaultAssay(mouse_brain)<-'RNA'
big.marker <- c(
  'Rbfox3',
  'Slc32a1',
  'Slc17a7',
  'Slc30a3',
  'Cux2',
  'Rorb',
  'Deptor',#'Ex-L5 IT',
  'Scnn1a',
  'Rspo1',#L4
  'Hsd11b1',
  'Batf3',
  'Oprk1',#L6 IT
  'Osr1',
  'Car3',
  'Fam84b',
  'Chrna6',
  'Pvalb',
  'Pappa2',
  'SIc17a8',
  'Trhr',
  'Tshz2',
  'Rapgef3',
  'Trh',
  'Gpr139',
  'Nxph4',
  'Rprm',
  'Crym',
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act
  'Rorb',#Ex-L4/5 IT
  'Slc17a8',#Ex-L5 NP
  'Cxcl14',#Ex-L5 NP Cxcl14
  'Npr3',#Ex-L5-PT
  'Foxp2',#Ex-L6 CT
  'Bmp3',#Ex-L6 IT Bmp3
  'Gfra1',#Ex-L6 IT Gfra1
  'Ccn2',#Ex-L6b
  'Ndst4',#Ex-PIR Ndst4
  'Tac1',#Misc
  'Drd2',#In-Drd2
  'Hap1','Scn5a',#In-Hap1
  'Cemip',#In-Pvalb
  'Sst','Chodl','Elfn1',#In-Sst
  'Tac1',#In-Tac1
  'Vip',#In-Vip/
  'Pax6',#In-Sncg
  'Aqp4',#Astro
  'Pdgfra','Olig1',#OPC
  'Mal',#Oligo
  'Col1a2','Vtn',#VLMC
  'Cldn5','Mcam',#ENdo
  'Mrc1','Cd163',#macrophage
  'Kcnj8',#peri, pericyte
  'Acta2'#SMC, smooth muscle cell
)
DotPlot(mouse_brain,features = big.marker %>% unique(),
        cluster.idents = T,group.by = 'seurat_clusters')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
DimPlot(mouse_brain,group.by  = 'celltype',label = T,cols = type.col)

Idents(mouse_brain) <-'seurat_clusters'
mouse_brain <- RenameIdents(mouse_brain,
                            '2' = 'Ex-L2/3 IT', 
                            '31' = 'unknown',
                            '9' = 'Ex-L4 IT',
                            '4' = 'Ex-L5 IT',
                            '21' = 'Ex-L5 NP',
                            '32' ='Ex-L5 NP Cxcl14',
                            '16' ='Ex-L5-PT',
                            '3' ='Ex-L5/6 IT',
                            '5' = 'Ex-L6 CT',
                            '24' = 'Ex-L6 IT',
                            '11' = 'Ex-L6 IT Bmp3',
                            '20' = 'Ex-L6b',
                            '25' ='Ex-PIR Ndst4','28' ='Ex-PIR Ndst4',
                            '8' ='In-Drd2','26' ='In-Drd2',
                            '7'='In-Hap1','19'='In-Hap1',
                            '12'='In-Pvalb',
                            '15' = 'In-Sst',
                            '6'='In-Tac1','14'='In-Tac1',
                            '13' = 'In-Vip',
                            '27' = 'Astro',
                            '1' = 'In-Sncg','30' = 'In-Sncg', '33' = 'In-Sncg',
                            '10' = 'OPC',
                            '0' = 'Oligo','22' = 'Oligo',
                            '23' = 'VLMC','29' = 'VLMC',           
                            '18'='Endo',
                            '17'='Macrophage'
)
mouse_brain$celltype <- Idents(mouse_brain)
DimPlot(mouse_brain,reduction = 'umap.rna',group.by  = 'celltype',label = T,cols = type.col,repel = T)
DimPlot(mouse_brain,reduction = 'wnn.umap',group.by  = 'celltype',label = T,repel = T,cols = type.col)+
  labs(title = 'ATAC celltype')+NoLegend()


#4
mouse_brain_AD<- readRDS('./scrna-atac/mouse_brain_10x_AD/mouse_brain_10x_AD_allobj.rds')
DefaultAssay(mouse_brain_AD)<-'ATAC'
peaks <- CallPeaks(mouse_brain_AD,macs2.path = '/home/chengww/anaconda3/envs/MACS/bin/macs2')
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(mouse_brain_AD),
  features = peaks,
  cells = colnames(mouse_brain_AD)
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# create a new assay using the MACS2 peak set and add it to the Seurat object
mouse_brain_AD[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments ="./scrna-atac/mouse_brain_10x_AD/Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_atac_fragments.tsv.gz",
  annotation = annotations
)
saveRDS(mouse_brain_AD,file = './scrna-atac/mouse_brain_10x_AD/mouse_brain_10x_AD_peaks_allobj.rds')

mouse_brain_AD$orig.ident <- paste('mouse_brain_10x_AD_',
                          lapply(rownames(mouse_brain_AD@meta.data),
                                 function(x) unlist(strsplit(x,'-')[[1]][[2]])) %>%
                            unlist(),
                          sep = '')
table(mouse_brain_AD$orig.ident)
tmp <- data.frame(sample=c('AD_17p9_rep4','AD_17p9_rep5','AD_2p5_rep2','AD_2p5_rep3',
                           'AD_5p7_rep2','AD_5p7_rep6','WT_13p4_rep2','WT_13p4_rep5',
                           'WT_2p5_rep2','WT_2p5_rep7','WT_5p7_rep2','WT_5p7_rep3'),
                  orig.ident=c('mouse_brain_10x_AD_1','mouse_brain_10x_AD_2','mouse_brain_10x_AD_3',
                               'mouse_brain_10x_AD_4','mouse_brain_10x_AD_5','mouse_brain_10x_AD_6',
                               'mouse_brain_10x_AD_7','mouse_brain_10x_AD_8','mouse_brain_10x_AD_9',
                               'mouse_brain_10x_AD_10','mouse_brain_10x_AD_11','mouse_brain_10x_AD_12'))

tmp_meta <- mouse_brain_AD@meta.data %>% tibble::rownames_to_column('cellid')
mouse_brain_AD@meta.data <- left_join(tmp_meta,tmp) %>% tibble::column_to_rownames('cellid')


DefaultAssay(mouse_brain_AD)<-'RNA'
DotPlot(mouse_brain_AD,features = big.marker %>% unique(),
        cluster.idents = T,group.by = 'celltype')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
DimPlot(mouse_brain_AD,group.by  = 'seurat_clusters',label = T)

Idents(mouse_brain_AD) <-'seurat_clusters'
big.marker <- c(
  'Rbfox3',
  'Slc32a1',
  'Slc17a7',
  'Slc30a3',
  'Cux2',
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act
  'Rorb',
  'Deptor',#'Ex-L5 IT',
  'Scnn1a',
  'Rspo1',#L4
  'Hsd11b1',
  'Oprk1',#L6 IT
  'Fam84b',
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act
  'Rorb',#Ex-L4/5 IT
  'Slc17a8',#Ex-L5 NP
  'Cxcl14',#Ex-L5 NP Cxcl14
  'Npr3',#Ex-L5-PT
  'Foxp2',#Ex-L6 CT
  'Bmp3',#Ex-L6 IT Bmp3
  'Gfra1',#Ex-L6 IT Gfra1
  'Ccn2',#Ex-L6b
  'Ndst4',#Ex-PIR Ndst4
  'Tac1',#Misc
  'Drd2',#In-Drd2
  'Hap1','Scn5a',#In-Hap1
  'Cemip',#In-Pvalb
  'Sst','Chodl','Elfn1',#In-Sst
  'Tac1',#In-Tac1
  'Vip',#In-Vip/
  'Pax6',#In-Sncg
  'Aqp4',#Astro
  'Lamp5',
  'Pdgfra','Olig1',#OPC
  'Mal',#Oligo
  'Col1a2','Vtn',#VLMC
  'Cldn5','Mcam',#ENdo
  'Mrc1','Cd163',#macrophage
  'Kcnj8',#peri, pericyte
  'Acta2'#SMC, smooth muscle cell
)
Idents(mouse_brain_AD) <-'seurat_clusters'
mouse_brain_AD <- RenameIdents(mouse_brain_AD,
                               '8' = 'Ex-L2/3 IT',
                               '1' = 'Ex-L2/3 IT Act', '22' = 'Ex-L2/3 IT Act', 
                               '5' = 'Ex-L4 IT',
                               '26' = 'Ex-L5 NP',
                               '19' ='Ex-L5-PT','21' ='Ex-L5-PT',
                               '17' ='unknown',
                               '6' = 'Ex-L6 CT',
                               '7' ='Ex-L5/6 IT',
                               '27' ='Ex-L6 IT',
                               '28' = 'Ex-L6 IT Gfra1',
                               '14' ='Ex-PIR Ndst4',
                               '9' ='In-Drd2',
                               '0'='In-Hap1',
                               '15'='In-Pvalb',
                               '20' = 'In-Sst',
                               '29'='In-Tac1','11' = 'In-Tac1', '12' = 'In-Tac1', '23' = 'In-Tac1', 
                               '16' = 'In-Vip',
                               '32' = 'Astro/Sncg','24' = 'Astro/Sncg','2' = 'Astro/Sncg', '10' = 'Astro/Sncg',
                               '13' = 'OPC','30' = 'OPC',
                               '3' = 'Oligo','4' = 'Oligo',
                               '33' = 'VLMC','34' = 'VLMC',           
                               '25'='Endo',
                               '18'='Macrophage',
                               '35'='Peri',
                               '31'='SMC'
)
mouse_brain_AD$celltype <- Idents(mouse_brain_AD)
DimPlot(mouse_brain_AD,reduction = 'umap.rna',group.by  = 'celltype',label = T,cols = type.col,repel = T)
DimPlot(mouse_brain_AD,reduction = 'wnn.umap',group.by  = 'celltype',label = T,repel = T,cols = type.col)+
  labs(title = 'ATAC celltype')+NoLegend()



#5
mouse_brain_e18<- readRDS('./scrna-atac/mouse_brain_10x_e18/mouse_brain_10x_e18_allobj.rds')
DefaultAssay(mouse_brain_e18)<-'ATAC'
peaks <- CallPeaks(mouse_brain_e18,macs2.path = '/home/chengww/anaconda3/envs/MACS/bin/macs2')
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(mouse_brain_e18),
  features = peaks,
  cells = colnames(mouse_brain_e18)
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# create a new assay using the MACS2 peak set and add it to the Seurat object
mouse_brain_e18[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments ="./scrna-atac/mouse_brain_10x_e18/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz",
  annotation = annotations
)
saveRDS(mouse_brain_e18,file = './scrna-atac/mouse_brain_10x_e18/mouse_brain_10x_e18_peaks_allobj.rds')


DefaultAssay(mouse_brain_e18)<-'RNA'
#Idents(mouse_brain_e18) <- 'seurat_clusters'
# markers <- FindAllMarkers(mouse_brain_e18,
#                           only.pos = T,
#                           min.pct = 0.25,
#                           logfc.threshold = 0.25)
# top2 <- markers %>% group_by(cluster) %>% top_n(3,avg_log2FC)
##remove cluster 17, high expression MT genes
mouse_brain_e18 <- subset(mouse_brain_e18,seurat_clusters!=17)
mouse_brain_e18 <- FindVariableFeatures(mouse_brain_e18, nfeatures = 3000)
mouse_brain_e18 <- NormalizeData(mouse_brain_e18)
mouse_brain_e18 <- ScaleData(mouse_brain_e18)
mouse_brain_e18 <- RunPCA(mouse_brain_e18, npcs = 30)
mouse_brain_e18 <- RunUMAP(mouse_brain_e18, dims = 1:30)
mouse_brain_e18 <- FindNeighbors(mouse_brain_e18, dims = 1:30)
mouse_brain_e18 <- FindClusters(mouse_brain_e18, resolution = 1.5, algorithm = 3)
DimPlot(mouse_brain_e18,reduction = 'umap.rna',group.by  = 'seurat_clusters',
        label = T,cols = type.col,repel = T)

big.marker <- c(
  'Rbfox3',
  'Slc32a1',
  'Slc17a7',
  'Slc30a3',
  'Cux2',
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act
  'Rorb',
  'Deptor',#'Ex-L5 IT',
  'Rorb',#Ex-L4/5 IT           20
  'Foxp2',#Ex-L6 CT           15,0 
  'Gfra1',#Ex-L6 IT Gfra1        21
  'Fabp7','Slc1a3','Apoe',#Astro/RG          4
  'Top2a','Mki67','Hmgn2',#IP-Hmgn2       8
  'Gadd45g',#IP-Gadd45g                    22
  'Eomes','Unc5d',#Ip-Eomes                 7
  'Reln','Car10','Ndnf','Trp73',#CR(L1)     24
  'Cntn2','Plxna2','Itpr1',#Ex-L2/3-Cntn2  3
  'Cux1','Satb2','Dab1',#Ex-L2/3-Cux1        6
  'Tenm3','Kcnq5',#Ex-L3/4-Tenm3              2,5
  'Foxp1','Unc5c','Rorb',#Ex-L3/4/5-Foxp1    16               
  'Galntl6','Cntn6','Fezf2','Bcl11b',#Ex-L4/5-Galntl6  9
  'Grik2','Hs6st3','Epha6',#Ex-L5/6-Epha6     10,12
  'Sox5','Hs3st4','Foxp2','Tle4',#Ex-L6-Tle4   1,19
  'Erbb4','Nxph1',#In-Nxph1                   11,18
  'Adarb2',#In-Adarb2                           17
  'Pdgfra','Olig1',#OPC                         25
  'Kdr','Cldn5','Mcam',#Endo                    23
  'Vtn','Kcnj8'#Peri                            23
  #'Apbb1ip'#Mic                                
)
DotPlot(mouse_brain_e18,features = big.marker %>% unique(),
        cluster.idents = T,group.by = 'seurat_clusters')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
Idents(mouse_brain_e18) <-'seurat_clusters'
mouse_brain_e18 <- RenameIdents(mouse_brain_e18,
                                '20' = 'Ex-L4/5 IT',
                                '15'='Ex-L6 CT','0'='Ex-L6 CT',
                                '21' = 'Ex-L6 IT Gfra1 ',
                                '4' = 'Astro/RG',
                                '8' = 'IP-Hmgn2',
                                '22' = 'IP-Gadd45g',
                                '7' = 'Ip-Eomes',
                                '24' = 'CR(L1)',
                                '3' ='Ex-L2/3-Cntn2',
                                '6' ='Ex-L2/3-Cux1',
                                '2' = 'Ex-L3/4-Tenm3','5' = 'Ex-L3/4-Tenm3',
                                '16' ='Ex-L3/4/5-Foxp1',
                                '9' = 'Ex-L4/5-Galntl6',
                                '10' = 'Ex-L5/6-Epha6','12' = 'Ex-L5/6-Epha6',
                                '1' = 'Ex-L6-Tle4','19' = 'Ex-L6-Tle4',
                                '11' = 'In-Nxph1','18' = 'In-Nxph1',
                                '17' = 'In-Adarb2',    
                                '25'='OPC',
                                '23'='Endo/Peri'
)
mouse_brain_e18$celltype <- Idents(mouse_brain_e18)
type.col <- colorRampPalette(c( brewer.pal(n=8,name='Paired'),
                                brewer.pal(n=5,name = 'Set3')))(57)
DimPlot(mouse_brain_e18,reduction = 'umap.rna',group.by  = 'celltype',label = T,cols = type.col,repel = T)
DimPlot(mouse_brain_e18,reduction = 'wnn.umap',group.by  = 'celltype',label = T,repel = T,cols = type.col)+
  labs(title = 'ATAC celltype')+NoLegend()



##integration---RNA---------------------------------------------------------------
dalist <- list(DietSeurat(ISSAAC_mCortex,assays = 'RNA'),
               # DietSeurat(GSE126074_P0_BrainCortex,assays = 'RNA'),
               DietSeurat(GSE126074_AdBrainCortex,assays = 'RNA'),
               DietSeurat(mouse_brain,assays = 'RNA'),
               DietSeurat(mouse_brain_e18,assays = 'RNA'),
               DietSeurat(mouse_brain_AD,assays = 'RNA'))

dalist <- lapply(X = dalist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = dalist)
dalist <- lapply(X = dalist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = dalist,reduction = "rpca",
                                  dims = 1:30)# reference = c(1, 2), 
brain.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

brain.integrated <- ScaleData(brain.integrated, verbose = FALSE)
brain.integrated <- RunPCA(brain.integrated, verbose = FALSE)
brain.integrated <- FindNeighbors(brain.integrated, dims = 1:30)
brain.integrated <- FindClusters(brain.integrated, resolution = 1)
brain.integrated <- RunUMAP(brain.integrated, dims = 1:30)
brain.integrated$celltype <- ifelse(brain.integrated$celltype=='In-Drd1 Tac1','In-Tac1',
                                    brain.integrated$celltype)


p1 <- DimPlot(brain.integrated, group.by = "orig.ident",cols = type.col[c(1,6,10,17,29,33)])+
  labs(title = 'Integration_RNA_samples')
p2<- DimPlot(brain.integrated, group.by = "celltype",cols = type.col,label = T,repel = T)+
  labs(title = 'Integration_RNA_celltypes')
p1|p2
saveRDS(brain.integrated,'./scrna-atac/brain.integrated_rna.rds')
DimPlot(brain.integrated, group.by = "orig.ident",split.by = 'orig.ident',
        cols = type.col[c(1,6,10,17,29,33)],ncol = 3)+
  labs(title = '')+NoLegend()
##重新进行注释-------------------------------------------------------------------
DefaultAssay(brain.integrated) <-'RNA'
DotPlot(brain.integrated,features = big.marker %>% unique(),
        cluster.idents = T,group.by = 'seurat_clusters')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
DimPlot(brain.integrated,group.by  = 'seurat_clusters',label = T)
DefaultAssay(brain.integrated)<-'RNA'
Idents(brain.integrated) <- 'seurat_clusters'
markers <- FindAllMarkers(brain.integrated,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
top2 <- markers %>% group_by(cluster) %>% top_n(3,avg_log2FC)



big.marker <- c(
  'Rbfox3',
  'Slc32a1',
  'Slc17a7',#'Slc30a3',#兴奋性神经元
  'Ptgs2',#"Ex-L2/3 IT"#Ex-L2/3 IT Act----- 0
  'Cpne6',##Ex-L2/3-Cpne6------16
  'Cntn2','Plxna2',#'Itpr1',#Ex-L2/3-Cntn2------6
  'Dab1',#'Cux1','Satb2',#Ex-L2/3-Cux1--------25            
  'Rorb',#Ex-L4/5 IT                     3
  'Cntn5', #Ex-L4/5-Cntn5------------------7
  'Cntn6',#'Fezf2',#'Bcl11b',#'Galntl6',#Ex-L4/5-Cntn6---28  
  'Deptor',#'Ex-L5 IT'--------------------8
  'Slc17a8',#Ex-L5 NP                     29
  'Cxcl14',#Ex-L5 NP Cxcl14-----------------30
  'Npr3',#Ex-L5-PT                        13
  'Epha6',#'Hs6st3','Grik2',#Ex-L5/6-Epha6 ---24                    
  'Hs3st4','Foxp2','Tle4',#'Sox5',#Ex-L6-Tle4------2
  'Oprk1',#L6 IT                         31
  #'Gfra1',#Ex-L6 IT Gfra1                42 
  'Bmp3',#Ex-L6 IT Bmp3                  11
  'Ccn2',#Ex-L6b--------------------------23
  #'Ndst4',#Ex-PIR Ndst4-------------------5
  #'Erbb4','Nxph1',#In-Nxph1 ---------14
  'Drd2',#In-Drd2-----------------------12
  'Hap1',#In-Hap1--------4
  'Cemip',#In-Pvalb-------14
  'Sst','Elfn1',#'Chodl',#In-Sst-------20
  'Tac1',#In-Tac1-------------9
 'Scn5a',
 'Pbx3',#NEW In-Pbx3------18
 'Adarb2',#In-Adarb2-------------
  'Vip',#In-Vip/Adarb2-------------17
  'Top2a','Mki67','Hmgn2',#IP-Hmgn2---37
  'Reln','Car10','Ndnf',#'Trp73',#CR(L1)------42
 'Fabp7','Slc1a3','Apoe',#Astro-Fabp7---------27
 'Aqp4',#Astro---------
 'Lars2',#'Map1a',#Astro-Lars2---------------15
 'Grm3',#'Frmd4a',Astro-Grm3-------------10
 'Agt',#Astro-Agt--------------22
  'Pdgfra','Olig1',#OPC-------19
  'Mal',#Oligo----------1,21
  'Col1a2',#VLMC------------36
  'Vtn','Kcnj8',#peri, pericyte----------------38
  'Cldn5','Mcam',#ENdo--------34
  'Mrc1','Cd163',#macrophage---26
  'Acta2'#SMC, smooth muscle cell---39
)
DotPlot(brain.integrated_all,features = big.marker %>% unique(),
        cluster.idents =F ,group.by = 'celltype_new')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  #mytheme +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))

#markers
tmp <- subset(brain.integrated1,seurat_clusters %in% c(10,22,15,27))
Idents(tmp) <-'seurat_clusters'
markers_new <- FindAllMarkers(tmp,
                              only.pos = T,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

Idents(brain.integrated_all) <-'seurat_clusters'
brain.integrated_all <- RenameIdents(brain.integrated_all,
                                 '0' = 'Ex-L2/3 IT',
                                 '16'='Ex-L2/3-Cpne6',
                                 '6' = 'Ex-L2/3-Cntn2',
                                '25' = 'Ex-L2/3-Dab1',
                                  '3' = 'Ex-L4/5 IT',
                                 '7' = 'Ex-L4/5-Cntn5',
                                 '5' = 'Ex-L4/5-Cntn5/6',
                                  '28'='Ex-L4/5-Cntn6',
                                 '8' = 'Ex-L5 IT',
                                 '29' = 'Ex-L5 NP',
                                 '30' ='Ex-L5 NP Cxcl14',
                                 '13' = 'Ex-L5 PT',
                                 '24' = 'Ex-L5/6-Epha6',
                                 '2'='Ex-L6-Tle4',
                                 '31'='Ex-L6 IT',
                                 '11'= 'Ex-L6 IT Bmp3',
                                 '23'='Ex-L6b',
                                 '5' ='Ex-PIR Ndst4',
                                 '12' ='In-Drd2',
                                 '4'='In-Hap1',
                                 '14'='In-Pvalb',
                                 '20' = 'In-Sst',
                                 '9'='In-Tac1',
                                '41'='In-Tac1-Scn5a',
                                 '18'='In-Pbx3',
                                 '17'='In-Vip/Adarb2',
                                 '37' = 'IP-Hmgn2',
                                 '42' = 'CR(L1)',
                                 '27' = 'Astro-Fabp7',
                                 '15' = 'Astro-Lars2',
                                '10' = 'Astro-Grm3',
                                '22' = 'Astro-Agt',
                                '19' = 'OPC',
                                 '1' = 'Oligo','21' = 'Oligo',
                                 '40' = 'Astro/Oligo',
                                 '36' = 'VLMC', 
                                 '38'='Peri',
                                 '34'='Endo',
                                 '26'='Macrophage',
                                 '39'='SMC'
)
brain.integrated_all$celltype_new <- Idents(brain.integrated_all)
DimPlot(brain.integrated_all,reduction = 'umap',group.by  = 'celltype_new',label = T,cols = type.col,repel = T)+
  labs(title = 'New celltypes after integration (RNA)')
saveRDS(brain.integrated,'./scrna-atac/brain.integrated_rna.rds')

brain.integrated <-readRDS('./scrna-atac/brain.integrated_rna.rds')
##remove mix clusters
brain.integrated1 <- subset(brain.integrated,seurat_clusters !=32 &seurat_clusters !=33 &seurat_clusters !=35)
DefaultAssay(brain.integrated1) <-'integrated'
brain.integrated1 <- ScaleData(brain.integrated1, verbose = FALSE)
brain.integrated1 <- RunPCA(brain.integrated1, verbose = FALSE)
brain.integrated1 <- FindNeighbors(brain.integrated1, dims = 1:30)
brain.integrated1 <- RunUMAP(brain.integrated1, dims = 1:30)
DefaultAssay(brain.integrated1) <-'RNA'
DotPlot(brain.integrated1,features = big.marker %>% unique(),
        cluster.idents = T,group.by = 'celltype_new')+
  RotatedAxis() +
  xlab('') +
  ylab('') +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "#C13533") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))

p1 <- DimPlot(brain.integrated1, reduction = 'umap',group.by = "orig.ident",cols = type.col[c(1,6,10,17,29,33)])+
  labs(title = 'Integration_RNA_samples')
p2<- DimPlot(brain.integrated1,reduction = 'umap', group.by = "celltype_new",cols = type.col,label = T,repel = T)+
  labs(title = 'Integration_RNA_celltypes')
p1|p2
DimPlot(brain.integrated1, group.by = "orig.ident", reduction = 'umap',split.by = "orig.ident",cols = type.col[c(1,6,10,17,29,33)])+
  labs(title = 'Integration_RNA_samples')

##sample informations--------------------
brain.integrated1$raw_cellid <- ifelse(brain.integrated1$orig.ident=='GSE126074_AdBrainCortex',
                                       gsub("\\_2",'',rownames(brain.integrated1@meta.data)[brain.integrated1$orig.ident=='GSE126074_AdBrainCortex']),
                                       lapply(rownames(brain.integrated1@meta.data),function(x) strsplit(x,'_')[[1]][1]) %>% unlist())

sample <- data.frame(raw_cellid=rownames(mouse_brain_AD@meta.data),
                     sample=mouse_brain_AD$sample) 
sample<- sample[sample$raw_cellid %in% brain.integrated1$raw_cellid,]
brain.integrated1$sample=NULL
ad.meta <- brain.integrated1@meta.data[brain.integrated1$orig.ident=='mouse_brain_10x_AD',] %>% 
    left_join(.,sample,by='raw_cellid')
ad.meta$int_cellid <- paste(ad.meta$raw_cellid,'_5',sep = '')
ad.meta <- dplyr::select(ad.meta,c(int_cellid,sample))

brain.integrated1$int_cellid<- rownames(brain.integrated1@meta.data)
brain.integrated1@meta.data <- brain.integrated1@meta.data %>% 
  tibble::rownames_to_column('tmpid') %>% 
  left_join(ad.meta,by='int_cellid') %>% 
  tibble::column_to_rownames('tmpid')

brain.integrated1$sample <- ifelse(is.na(brain.integrated1$sample),
                                   brain.integrated1$orig.ident, brain.integrated1$sample )
ad.meta <- dplyr::filter(ad.meta,sample %in% c("WT_13p4_rep5","WT_2p5_rep7","WT_2p5_rep2",
                                      "WT_5p7_rep2","WT_13p4_rep2","WT_5p7_rep3"))
brain.integrated1$Patient <- ifelse(brain.integrated1$int_cellid %in% ad.meta$int_cellid,
                                   'mouse_brain_10x_WT',brain.integrated1$orig.ident )
brain.integrated1$AD <- ifelse(brain.integrated1$Patient=='mouse_brain_10x_AD','AD','others')
brain.integrated1$Patient <- factor(brain.integrated1$Patient,levels =
                                      c("GSE126074_AdBrainCortex","ISSAAC_mCortex",
                                        "mouse_brain_10x","mouse_brain_10x_WT",'mouse_brain_10x_AD',"mouse_brain_10x_e18") )
brain.integrated1$E18 <- ifelse(brain.integrated1$Patient=='mouse_brain_10x_e18','E18','others')
brain.integrated1$AD_E18 <- ifelse(brain.integrated1$Patient=='mouse_brain_10x_AD','AD',ifelse(brain.integrated1$Patient=='mouse_brain_10x_e18','E18','Others'))


##ISSAAC_mCortex have 2 reps---------
rep1 <-brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='ISSAAC_mCortex'][
  stringr::str_detect(brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='ISSAAC_mCortex'],'-1')]
rep2 <-brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='ISSAAC_mCortex'][
  stringr::str_detect(brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='ISSAAC_mCortex'],'-2')]
brain.integrated1$sample <-ifelse(brain.integrated1$raw_cellid %in% rep1
  ,'ISSAAC_mCortex_rep1',
  ifelse(brain.integrated1$raw_cellid %in% rep2,'ISSAAC_mCortex_rep2',brain.integrated1$sample ))
##GSE126074_AdBrainCortex have 12 reps-------
tmp <- data.frame(tmpid=c('09A','09B','09C','09D','09E','09F','09G','09H','09I',
                          '09J',  '09K',  '09L'),
                  tmpsample=c('GSE126074_AdBrainCortex_rep1','GSE126074_AdBrainCortex_rep2',
                              'GSE126074_AdBrainCortex_rep3','GSE126074_AdBrainCortex_rep4',
                              'GSE126074_AdBrainCortex_rep5','GSE126074_AdBrainCortex_rep6',
                              'GSE126074_AdBrainCortex_rep7','GSE126074_AdBrainCortex_rep8',
                              'GSE126074_AdBrainCortex_rep9','GSE126074_AdBrainCortex_rep10',
                              'GSE126074_AdBrainCortex_rep11','GSE126074_AdBrainCortex_rep12'
                              ))
brain.integrated1$tmpid <- substr(brain.integrated1$raw_cellid,1,3) 
brain.integrated1@meta.data <- tibble::rownames_to_column(brain.integrated1@meta.data,'cell') %>% 
  left_join(tmp,by='tmpid') %>% tibble::column_to_rownames('cell')
brain.integrated1$sample <- ifelse(brain.integrated1$sample=='GSE126074_AdBrainCortex',brain.integrated1$tmpsample,brain.integrated1$sample)

brain.integrated1@meta.data <- brain.integrated1@meta.data %>% 
  dplyr::select("orig.ident","nCount_RNA","nFeature_RNA",               
                "percent.mt","nCount_ATAC","nFeature_ATAC" ,             
                "TSS.enrichment","TSS.percentile","nucleosome_signal",          
                "nucleosome_percentile","blacklist_fraction","RNA_snn_res.1",              
                "ATAC_snn_res.1","wsnn_res.1","RNA_UMAP1",
                "RNA_UMAP2", "ATAC_UMAP1","ATAC_UMAP2",                 
                "batch",  "ATAC_Cluster", "RNA_Cluster_Annotation",     
                "RNA_Cluster_colour" ,"celltype",  "RNA.weight",                 
                "peaks.weight", "seurat_clusters" ,           
                "ATAC.weight",    "integrated_snn_res.1", "celltype_new",          
                "raw_cellid","Patient",  
                "AD", "E18", "int_cellid",  "sample",    "AD_E18"                       
                                 )

##save integrated RNA object
saveRDS(brain.integrated1,'./scrna-atac/brain.integrated_rna_sub.rds')

#DEGs for each cluster
DefaultAssay(brain.integrated1)<-'RNA'
Idents(brain.integrated1) <- 'celltype_new'
markers_new <- FindAllMarkers(brain.integrated1,
                              only.pos = T,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
top2 <- markers_new %>% group_by(cluster) %>% top_n(2,avg_log2FC)

fwrite(markers_new,
       './scrna-atac/all_celltype_new_markers.txt',sep='\t')
