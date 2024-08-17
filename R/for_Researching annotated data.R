##Load Packages
library(Signac)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(EnsDb.Mmusculus.v79)

#Set your Path
Path='/home/chengww/data/project/database'
setwd(Path)
# load the RNA and ATAC data
counts_gene <- read.csv("./scrna-atac/E-MTAB-11264/mCortex_all_gene_expression_matrix.csv",header  = T,row.names = 1)
fcounts_atac <- read.csv("./scrna-atac/E-MTAB-11264/mCortex_all_ATAC_peak_matrix.csv",header = T,row.names = 1)
colnames(counts_gene) <- gsub('\\.','-',colnames(counts_gene))
colnames(fcounts_atac) <- gsub('\\.','-',colnames(fcounts_atac))
fragpath <- "./scrna-atac/E-MTAB-11264/fragments_all.tsv.gz"

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"


# create a Seurat object containing the RNA adata
brain <- CreateSeuratObject(
  counts =counts_gene,
  assay = "RNA"
)

# create ATAC assay and add it to the object
brain[["ATAC"]] <- CreateChromatinAssay(
  counts = fcounts_atac,
  sep = c(":", "-"),
  fragments = fragpath
)

# add the gene information to the object
Annotation(brain[["ATAC"]]) <- annotations
head(brain@assays$ATAC@counts[,1:5])



##load raw metadata for annotation
metadata <- read.csv('./scrna-atac/E-MTAB-11264/Sample_Info.csv',header = T,row.names = 1)
meta2<- brain@meta.data %>% tibble::rownames_to_column('cellID')
metadata$cellID <- rownames(metadata)
brain@meta.data <- left_join(meta2 %>% dplyr::select(cellID),metadata,by='cellID') %>% tibble::column_to_rownames('cellID')
brain$celltype <- brain$RNA_Cluster_Annotation

VlnPlot(
  object = brain,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  group.by ='orig.ident', 
  pt.size = 0.1) &
  mytheme &
  NoLegend() & 
  xlab(label = '')


##RNA data analysis
DefaultAssay(brain) <- "RNA"
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
brain <- ScaleData(brain)
brain <- RunPCA(brain)
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 1)
brain <- RunUMAP(brain, dims = 1:30,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DimPlot(brain, group.by = "celltype", label = T,raster =F,
        cols = type.col,repel = T) + 
  ggtitle("RNA")


##ATAC data analysis
### We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(brain) <- "ATAC"##using peak matrix
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = "q0")
brain <- RunSVD(brain)
brain<- RunUMAP(brain, reduction = "lsi", dims = 2:30,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(brain, group.by = "celltype", label = T,cols = type.col,repel = T) + 
  ggtitle("ATAC")


##base WNN plot UMAP figure
brain <- FindMultiModalNeighbors(brain, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
brain <- RunUMAP(brain, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#brain <- FindClusters(brain, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

pp1 <- DimPlot(brain, reduction = "umap.rna", 
               group.by = "celltype",
               label = TRUE, 
               label.size = 4, 
               repel = TRUE,
               cols = type.col) + ggtitle("RNA")
pp2 <- DimPlot(brain, reduction = "umap.atac", 
               group.by = "celltype", 
               label = TRUE, 
               label.size = 4, 
               repel = TRUE,
               cols = type.col) + ggtitle("ATAC")
pp3 <- DimPlot(brain, reduction = "wnn.umap", 
               group.by = "celltype", 
               label = TRUE, 
               label.size = 4, 
               repel = TRUE,
               cols = type.col) + ggtitle("WNN")
pp1 + pp2 + pp3 & theme(plot.title = element_text(hjust = 0.5))

saveRDS(brain,file = '/home/chengww/data/project/database/E-MTAB-11264/ISSAAC_mCortex_allobj.rds')


#DEGs for eaxh cluster-----------------------
DefaultAssay(brain) <- "RNA"
Idents(brain) <- 'celltype'
markers <- FindAllMarkers(brain,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
fwrite(markers,'/home/chengww/data/project/database/database/data/ISSAAC_mCortex_markers.txt',sep='\t')


#Peaks for eaxh cluster-----------------------
DefaultAssay(brain) <-'peaks'
a<-nrow(brain)#169180
markers.peak <- FindAllMarkers(brain,
                               only.pos = T,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
fwrite(markers.peak,
       '/home/chengww/data/project/database/database/data/ISSAAC_mCortex_markers.peak.txt',sep='\t')


#Count the number of cells in each cell cluster-----------------------
cellnumber <- as.data.frame(table(brain$celltype))
colnames(cellnumber) <-c('ClusterID','Number')
fwrite(cellnumber,
       '/home/chengww/data/project/database/database/data/ISSAAC_mCortex_cellnumber.txt',sep='\t')


#Process fragment files to generate bw files as input to peaks browser-----------------------
#Save cluster names
celltype = brain@meta.data %>% dplyr::select('celltype')
head(celltype)
celltype$celltype = gsub(' ','',celltype$celltype)
celltype$celltype = gsub('-','',celltype$celltype)
celltype$celltype = gsub('/','',celltype$celltype)
fwrite(celltype,
       '/home/chengww/data/project/database/E-MTAB-11264/celltype.txt',
       sep = '\t',row.names = T,col.names = F)   

#Split the *fragments.tsv file by clusters
fragment <- read.table("./fragments_all.tsv", header=F, sep="\t")
celltype <- read.table('celltype.txt',sep='\t')
for (i in unique(celltype$V2)){
  cellid = celltype[which(celltype[2]==i),1]
  fragment.tmp = fragment[which(fragment[,4] %in% cellid),]
  fwrite(fragment.tmp,paste('./bwfile/',i,'.bed',sep = ''),sep='\t',col.name=F)
}
system('nohup sh /data3/chengww/project/database/human_brain_10x/bwfile/fragment_bw.sh &')


####Generating a cell activity matrix for ATAC data---------------------------------
DefaultAssay(pbm) <-'peaks'
rownames(pbm)
gene.activities <- GeneActivity(pbm)
head(gene.activities[1:5,1:5])
dim(gene.activities)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbm[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
pbm <- NormalizeData(
  object = pbm,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbm$nCount_RNA)
)

DefaultAssay(pbm) <- 'ACTIVITY'

FeaturePlot(
  object = pbm,
  reduction = 'wnn.umap',
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'CCR7'),
  pt.size = 0.1,
  cols =c("lightgrey",'#BE766E','#C13533'),
  #assay = 'atac_gene',
  max.cutoff = 'q95',
  ncol = 3
)

saveRDS(brain,file = '/home/chengww/data/project/database/E-MTAB-11264/ISSAAC_mCortex_allobj_activity.rds')
