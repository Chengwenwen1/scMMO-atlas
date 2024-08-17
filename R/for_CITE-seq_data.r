##Load Packages
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(RColorBrewer)
type.col <- colorRampPalette(c( brewer.pal(n=12,name='Paired'),
                                brewer.pal(n=12,name = 'Set3')))(39)


#Set your Path
Path='/home/chengww/data/project/database'
setwd(Path)
##load raw data
load("./GSE245311_multimodal_harmony_integrated_controls.RData")
load('./GSE245311_multimodal_harmony_integrated_gb.RData')
control <- refquery_merged
gb <- refquery_hm
dim(control@assays$ADT@counts)
head(control@assays$ADT@counts[1:5,1:5])
dim(control@assays$RNA@counts)
head(control@assays$RNA@counts[1:5,1:5])

p3 <- DimPlot(control, reduction = 'umap', group.by = 'predicted.celltype', label = TRUE,
              cols = type.col,
              repel = TRUE) + NoLegend()
p4 <- DimPlot(control, reduction = 'wnn.umap', group.by = 'predicted.celltype', label = TRUE,
              cols = type.col,
              repel = TRUE) 
p3 + p4

VlnPlot(control, features = c('nCount_RNA', 'nFeature_RNA',
                         'nCount_ADT', 'nFeature_ADT'), 
        group.by = 'orig.ident', sort = TRUE,
        ncol = 2,
        pt.size = 0) +
  NoLegend()
DefaultAssay(control) <-'ADT'
FeaturePlot(control,features = c("Hu.CD8"),
            reduction = 'wnn.umap', max.cutoff = 5, 
            cols = c("lightgrey","darkgreen"), ncol = 1)

CITE_seq_Brain_control_H@meta.data <- dplyr::select(CITE_seq_Brain_control_H@meta.data,
                                                    c('orig.ident',"nCount_RNA","nFeature_RNA",
                                                      "nCount_ADT", "nFeature_ADT",'celltype'))
CITE_seq_Brain_control_H$celltype <- factor(CITE_seq_Brain_control_H$celltype,
                                            levels = c(
                                              "HSPC","Eryth" , "Platelet",
                                              "CD14 Mono", "CD16 Mono","Plasmablast",                 
                                              "cDC1", "cDC2","pDC",          
                                              "B naive","B intermediate","B memory",                   
                                              "CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL",  
                                              "CD8 Naive","CD8 TCM","CD8 TEM",
                                              "CD4 Proliferating","CD8 Proliferating",
                                                "Treg", "MAIT", "dnT","gdT",                           
                                              "ILC", "NK" , "NK Proliferating" ,"NK_CD56bright"    
                                            ))

Idents(CITE_seq_Brain_control_H) <-'celltype'
DimPlot(CITE_seq_Brain_control_H, reduction = 'wnn.umap',  label = TRUE,
        cols = type.col,
        repel = TRUE) 
DefaultAssay(CITE_seq_Brain_control_H) <-'ADT'
FeaturePlot(CITE_seq_Brain_control_H,features = c("Hu.CD8"),
            reduction = 'wnn.umap', max.cutoff = 5, 
            cols = c("lightgrey","darkgreen"), ncol = 1)
DefaultAssay(CITE_seq_Brain_control_H) <-'RNA'
FeaturePlot(CITE_seq_Brain_control_H,features = c("CD8B"),
            reduction = 'wnn.umap', max.cutoff = 5, 
            cols = c("lightgrey","darkgreen"), ncol = 1)
saveRDS(CITE_seq_Brain_control_H,
        file = '/data/chengww/project/database/rna_protein/CITE_seq_Brain_control_H.rds')


#DEGs for eaxh cluster
control$predicted.celltype <- factor(control$predicted.celltype,
                                     levels = c(
                                       "HSPC","Eryth" , "Platelet",
                                       "CD14 Mono", "CD16 Mono","Plasmablast",                 
                                       "cDC1", "cDC2","pDC",          
                                       "B naive","B intermediate","B memory",                   
                                       "CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL",  
                                       "CD8 Naive","CD8 TCM","CD8 TEM",
                                       "CD4 Proliferating","CD8 Proliferating",
                                       "Treg", "MAIT", "dnT","gdT",                           
                                       "ILC", "NK" , "NK Proliferating" ,"NK_CD56bright"    
                                     ))
Idents(control) <-'predicted.celltype'
DefaultAssay(control)<-'RNA'
markers <- FindAllMarkers(control,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
fwrite(markers,'/home/chengww/data/project/database/database/data/scrna_protein/CITE_seq_Brain_control_H_markers.txt',sep='\t')


#Proteins for eaxh cluster
DefaultAssay(control) <-'ADT'
markers.peak <- FindAllMarkers(control,
                               only.pos = T,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
fwrite(markers.peak,
       '/home/chengww/data/project/database/database/data/scrna_protein/CITE_seq_Brain_control_H_markers.protein.txt',sep='\t')


####glioblastoma的数据处理---tumor-------------------------------------
gb <- refquery_hm
dim(gb@assays$ADT@counts)
head(gb@assays$ADT@counts[1:5,1:5])
dim(gb@assays$RNA@counts)
head(gb@assays$RNA@counts[1:5,1:5])

p3 <- DimPlot(gb, reduction = 'umap', group.by = 'predicted.celltype', label = TRUE,
              cols = type.col,
              repel = TRUE) + NoLegend()
p4 <- DimPlot(gb, reduction = 'wnn.umap', group.by = 'predicted.celltype', label = TRUE,
              cols = type.col,
              repel = TRUE) 
p3 + p4

VlnPlot(gb, features = c('nCount_RNA', 'nFeature_RNA',
                              'nCount_ADT', 'nFeature_ADT'), 
        group.by = 'predicted.celltype', sort = TRUE,
        ncol = 2,
        pt.size = 0) +
  NoLegend()
DefaultAssay(gb) <-'ADT'
FeaturePlot(gb,features = c("Hu.CD8"),
            reduction = 'wnn.umap', max.cutoff = 5, 
            cols = c("lightgrey","darkgreen"), ncol = 1)
CITE_seq_Brain_tumor_H <- gb
CITE_seq_Brain_tumor_H
CITE_seq_Brain_tumor_H@assays$SCT<-NULL
CITE_seq_Brain_tumor_H@reductions$umap <- NULL
CITE_seq_Brain_tumor_H@graphs$wsnn <-NULL
CITE_seq_Brain_tumor_H@graphs$wknn <- NULL
CITE_seq_Brain_tumor_H@graphs$RNA_nn <- NULL
CITE_seq_Brain_tumor_H@graphs$RNA_snn <- NULL
CITE_seq_Brain_tumor_H@neighbors$weighted.nn <- NULL
CITE_seq_Brain_tumor_H@reductions$pca<-NULL
CITE_seq_Brain_tumor_H@reductions$lsi<- NULL
CITE_seq_Brain_tumor_H@reductions$harmony<-NULL
CITE_seq_Brain_tumor_H@reductions$apca<- NULL
CITE_seq_Brain_tumor_H@reductions$spca<- NULL
CITE_seq_Brain_tumor_H@assays$prediction.score.celltype<-NULL
CITE_seq_Brain_tumor_H@assays$prediction.score.celltype2<-NULL
CITE_seq_Brain_tumor_H@assays$prediction.score.celltype3<-NULL
CITE_seq_Brain_tumor_H@assays$RNA@scale.data<-matrix()
CITE_seq_Brain_tumor_H@assays$RNA@counts<-matrix()
CITE_seq_Brain_tumor_H@assays$ADT@counts<-matrix()
CITE_seq_Brain_tumor_H@assays$ADT@scale.data<-matrix()
CITE_seq_Brain_tumor_H@assays$predicted_ADT<-NULL
CITE_seq_Brain_tumor_H@assays$ADT_biolegend<-NULL
CITE_seq_Brain_tumor_H@assays$predicted_predicted_ADT<-NULL
CITE_seq_Brain_tumor_H$celltype <- gb$predicted.celltype
CITE_seq_Brain_tumor_H$orig.ident <- gb$sample
CITE_seq_Brain_tumor_H@meta.data <- dplyr::select(CITE_seq_Brain_tumor_H@meta.data,
                                                    c('orig.ident',"nCount_RNA","nFeature_RNA",
                                                      "nCount_ADT", "nFeature_ADT",'celltype'))
CITE_seq_Brain_tumor_H$celltype <- factor(CITE_seq_Brain_tumor_H$celltype,
                                            levels = c(
                                              "HSPC","Eryth" , "Platelet",
                                              "CD14 Mono", "CD16 Mono","Plasmablast",                 
                                              "cDC1", "cDC2","pDC",          
                                              "B naive","B intermediate","B memory",                   
                                              "CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL",  
                                              "CD8 Naive","CD8 TCM","CD8 TEM",
                                              "CD4 Proliferating","CD8 Proliferating",
                                              "Treg", "MAIT", "dnT","gdT",                           
                                              "ILC", "NK" , "NK Proliferating" ,"NK_CD56bright"    
                                            ))

Idents(CITE_seq_Brain_tumor_H) <-'celltype'
DimPlot(CITE_seq_Brain_tumor_H, reduction = 'wnn.umap',  label = TRUE,
        cols = type.col,
        repel = TRUE) 
DefaultAssay(CITE_seq_Brain_tumor_H) <-'ADT'
FeaturePlot(CITE_seq_Brain_tumor_H,features = c("Hu.CD8"),
            reduction = 'wnn.umap', max.cutoff = 5, 
            cols = c("lightgrey","darkgreen"), ncol = 1)
DefaultAssay(CITE_seq_Brain_tumor_H) <-'RNA'
FeaturePlot(CITE_seq_Brain_tumor_H,features = c("CD8B"),
            reduction = 'wnn.umap', max.cutoff = 5, 
            cols = c("lightgrey","darkgreen"), ncol = 1)
saveRDS(CITE_seq_Brain_tumor_H,
        file = '/data/chengww/project/database/rna_protein/CITE_seq_Brain_tumor_H.rds')


#DEGs for eaxh cluster-------------------------------------------
gb$predicted.celltype <- factor(gb$predicted.celltype,
                                     levels = c(
                                       "HSPC","Eryth" , "Platelet",
                                       "CD14 Mono", "CD16 Mono","Plasmablast",                 
                                       "cDC1", "cDC2","pDC",          
                                       "B naive","B intermediate","B memory",                   
                                       "CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL",  
                                       "CD8 Naive","CD8 TCM","CD8 TEM",
                                       "CD4 Proliferating","CD8 Proliferating",
                                       "Treg", "MAIT", "dnT","gdT",                           
                                       "ILC", "NK" , "NK Proliferating" ,"NK_CD56bright"    
                                     ))
Idents(gb) <-'predicted.celltype'
DefaultAssay(gb)<-'RNA'
markers <- FindAllMarkers(gb,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
fwrite(markers,'/home/chengww/data/project/database/database/data/scrna_protein/CITE_seq_Brain_tumor_H_markers.txt',sep='\t')

#DProteins for eaxh cluster-------------------------------------------
DefaultAssay(gb) <-'ADT'
markers.peak <- FindAllMarkers(gb,
                               only.pos = T,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
fwrite(markers.peak,
       '/home/chengww/data/project/database/database/data/scrna_protein/CITE_seq_Brain_tumor_H_markers.protein.txt',sep='\t')

#Count the number of cells in each cell cluster-----------------------
cellnumber <- as.data.frame(table(gb$predicted.celltype))
colnames(cellnumber) <-c('ClusterID','Number')
fwrite(cellnumber,
       '/home/chengww/data/project/database/database/data/scrna_protein/CITE_seq_Brain_tumor_H_cellnumber.txt',sep='\t')

