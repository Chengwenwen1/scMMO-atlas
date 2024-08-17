##Load Packages
library(Signac)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
type.col <- colorRampPalette(c(brewer.pal(n=10,name = 'Paired'),
                               brewer.pal(n=8,name='Dark2')))(30)


#Set your Path
Path='/home/chengww/data/project/database'
setwd(Path)

# load the RNA and ATAC data
counts <- Read10X_h5("./scrna-atac/mouse_brain_10x_e18/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
fragpath <- "./scrna-atac/mouse_brain_10x_e18/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz"


# create a Seurat object containing the RNA adata
brain <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
brain[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragpath
)
brain
table(brain$orig.ident)
brain$orig.ident <- 'mouse_brain_10x_e18'


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
Annotation(brain[["ATAC"]] ) <- annotations

##QC
DefaultAssay(brain) <- "ATAC"
brain <- NucleosomeSignal(brain)
brain <- TSSEnrichment(brain)
brain$blacklist_fraction <- FractionCountsInRegion(
  object = brain,
  assay = 'ATAC',
  regions = blacklist_mm10
)
VlnPlot(
  brain,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)
brain <- subset(
  x = brain,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)

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

#Integration with scRNA-seq data
# label transfer from Allen brain 
#allen_brain.rds download from Seurat web
allen <-  readRDS("./scrna-atac/allen_brain.rds")
DefaultAssay(allen)<-'RNA'
# use the RNA assay in the brain-seq data for integration with scRNA-seq
DefaultAssay(brain) <- 'RNA'
transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = brain,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = brain[['pca']],
  dims = 1:30
)

brain <- AddMetaData(object = brain, metadata = predicted.labels)
DimPlot(brain,  label = T,cols = type.col,group.by = 'predicted.id') + 
  ggtitle("RNA")
brain$celltype <- brain$predicted.id
brain$celltype <- factor(brain$celltype,levels=unique(brain$celltype))
DimPlot(brain,  label = T,cols = type.col,group.by = 'celltype') + 
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

saveRDS(brain,file = '/home/chengww/data/project/database/scrna-atac/mouse_brain_10x_e18/mouse_brain_10x_e18_allobj.rds')

#DEGs for eaxh cluster-----------------------
DefaultAssay(brain) <- 'RNA'
Idents(brain) <- 'celltype'
markers <- FindAllMarkers(brain,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
fwrite(markers,'/home/chengww/data/project/database/database/data/mouse_brain_10x_e18_markers.txt',sep='\t')


#Peaks for eaxh cluster-----------------------
DefaultAssay(brain) <-'ATAC'
nrow(brain)#172193
markers.peak <- FindAllMarkers(brain,
                               only.pos = T,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
fwrite(markers.peak,
       '/home/chengww/data/project/database/database/data/mouse_brain_10x_e18_markers.peak.txt',sep='\t')

#Count the number of cells in each cell cluster-----------------------
cellnumber <- as.data.frame(table(brain$celltype))
colnames(cellnumber) <-c('ClusterID','Number')
fwrite(cellnumber,
       '/home/chengww/data/project/database/database/data/mouse_brain_10x_e18_cellnumber.txt',sep='\t')



#Process fragment files to generate bw files as input to peaks browser------------------------
#Save cluster names
celltype = brain@meta.data %>% dplyr::select('celltype')
head(celltype)
celltype$celltype = gsub(' ','',celltype$celltype)
celltype$celltype <-gsub('/','',celltype$celltype)
fwrite(celltype,
       '/home/chengww/data/project/database/mouse_brain_10x_e18/celltype.txt',
       sep = '\t',row.names = T,col.names = F)   

#Split the *fragments.tsv file by clusters
fragment <- read.table("./e18_mouse_brain_fresh_5k_atac_fragments.tsv", header=F, sep="\t")
celltype <- read.table('celltype.txt',sep='\t')
for (i in unique(celltype$V2)){
  cellid = celltype[which(celltype[2]==i),1]
  fragment.tmp = fragment[which(fragment[,4] %in% cellid),]
  fwrite(fragment.tmp,paste('./bwfile/',i,'.bed',sep = ''),sep='\t',col.name=F)
}
system('nohup sh /data/chengww/project/database/mouse_brain_10x_e18/bwfile/fragment_bw.sh &')



####Generating a cell activity matrix for ATAC data---------------------------------
DefaultAssay(brain) <-'ATAC'
rownames(brain)
gene.activities <- GeneActivity(brain)
head(gene.activities[1:5,1:5])
dim(gene.activities)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
brain[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)

saveRDS(brain,file = '/home/chengww/data/project/database/scrna-atac/mouse_brain_10x_e18/mouse_brain_10x_e18_allobj_activity.rds')

