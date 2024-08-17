##Load Packages
library(Signac)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#Set your Path
Path='/home/chengww/data/project/database'
setwd(Path)

# load the RNA and ATAC data
counts <- Read10X_h5("./scrna-atac//human_brain_10x/human_brain_3k_filtered_feature_bc_matrix.h5")
fragpath <- "./scrna-atac/human_brain_10x/human_brain_3k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

brain <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
brain[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
brain$orig.ident <- 'human_brain_3k'

##QC
DefaultAssay(brain) <- "ATAC"
brain <- NucleosomeSignal(brain)
brain <- TSSEnrichment(brain)

VlnPlot(
  object = brain,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  group.by ='orig.ident', 
  pt.size = 0.1) &
  mytheme &
  NoLegend() & 
  xlab(label = '')

# filter out low quality cells
brain <- subset(
  x = brain,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)


##peak calling
# call peaks using MACS2
peaks <- CallPeaks(brain,macs2.path = '/home/chengww/anaconda3/envs/MACS/bin/macs2')

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(brain),
  features = peaks,
  cells = colnames(brain)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
brain[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

##RNA data analysis
DefaultAssay(brain) <- "RNA"
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
brain <- ScaleData(brain)
brain <- RunPCA(brain)
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 1)
brain <- RunUMAP(brain, dims = 1:30,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

##markers for celltypes annotation
marker<-c(
  'GFAP','SLC1A2','SLC4A4','RYR3',#'SLC24A4',#astrocytes (), 
  'GRM4','GAD1','GRIK2','RBFOX1',#cerebellar granule  (Gran)
  'CEMIP','DCN','COL1A2','COL1A1','COL3A1',#'ACTA2',#fibroblast (), 
  'CSF1R','APBB1IP','C3',#microglia () 
  'CAMK2A','SLC17A7','SLC17A6',#neuronal (),
  'PLP1','MOBP','MOG',#oligodendrocytes (),
  #'SOX2',#progenitors
  'PCDH15','PDGFRA'#OPC
)

DotPlot(brain,features =unique(marker),cluster.idents = F,
        group.by = 'seurat_clusters' )+coord_flip()

#Remove the mix subgroups and repeat the clustering process above
brain <- subset(brain, seurat_clusters !=11)

##Annotation cell clusters
n=length(unique(brain$seurat_clusters))
celltype=data.frame(ClusterID=0:(n-1),celltype=0:(n-1),stringsAsFactors = F)
celltype[celltype$ClusterID %in% c(6),2]='fibroblast'
celltype[celltype$ClusterID %in% c(3),2]='astrocytes'
celltype[celltype$ClusterID %in% c(7),2]='neuronal'
celltype[celltype$ClusterID %in% c(0,1,2,10),2]='oligodendrocytes'
celltype[celltype$ClusterID %in% c(9),2]='microglia'
celltype[celltype$ClusterID %in% c(5),2]='cerebellar granule1'
celltype[celltype$ClusterID %in% c(8),2]='cerebellar granule2'
celltype[celltype$ClusterID %in% c(4),2]='precursor cells'
brain@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  brain@meta.data[which(brain@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(brain$celltype)

##UMAP shows cell subpopulations
type.col <- colorRampPalette(c( brewer.pal(n=8,name='Paired')))(8)
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

##save data object
saveRDS(brain,file = '/home/chengww/data3/project/database/scrna-atac/human_brain_10x/human_brain_10x_allobj.rds')


#DEGs for eaxh cluster-----------------------
DefaultAssay(brain) <- 'RNA'
Idents(brain) <- 'celltype'
markers <- FindAllMarkers(brain,
                          only.pos = T,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
fwrite(markers,'/home/chengww/data3/project/database/database/data/human_brain_10x_markers.txt',sep='\t')

#Peaks for eaxh cluster-----------------------
DefaultAssay(brain) <-'ATAC'
nrow(brain) # Num.peaks = 134030
markers.peak <- FindAllMarkers(brain,
                               only.pos = T,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
fwrite(markers.peak,
       '/home/chengww/data3/project/database/database/data/human_brain_10x_markers.peak.txt',sep='\t')

#Count the number of cells in each cell cluster-----------------------
cellnumber <- as.data.frame(table(brain$celltype))
colnames(cellnumber) <-c('ClusterID','Number')
fwrite(cellnumber,
       '/home/chengww/data3/project/database/database/data/human_brain_10x_cellnumber.txt',sep='\t')


#Process fragment files to generate bw files as input to peaks browser------------------------
#Save cluster names
celltype = brain@meta.data %>% dplyr::select('celltype')
head(celltype)
celltype$celltype = gsub(' ','',celltype$celltype)
#celltype$celltype <-gsub('/','',celltype$celltype)
fwrite(celltype,
       '/home/chengww/data3/project/database/scrna-atac/human_brain_10x/celltype.txt',
       sep = '\t',row.names = T,col.names = F)   
#Split the *fragments.tsv file by clusters
fragment <- read.table("./human_brain_3k_atac_fragments.tsv", header=F, sep="\t")
celltype <- read.table('celltype.txt',sep='\t')
for (i in unique(celltype$V2)){
  cellid = celltype[which(celltype[2]==i),1]
  fragment.tmp = fragment[which(fragment[,4] %in% cellid),]
  fwrite(fragment.tmp,paste('./bwfile/',i,'.bed',sep = ''),sep='\t',col.name=F)
}

system('nohup sh /data3/chengww/project/database/human_brain_10x/bwfile/fragment_bw.sh &')




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
  scale.factor = median(pbm$nCount_RNA)
)
##Check gene activity
DefaultAssay(brain) <- 'ACTIVITY'
FeaturePlot(
  object = brain,
  reduction = 'wnn.umap',
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'CCR7'),
  pt.size = 0.1,
  cols =c("lightgrey",'#BE766E','#C13533'),
  #assay = 'atac_gene',
  max.cutoff = 'q95',
  ncol = 3
)

saveRDS(brain,file = '/home/chengww/data/project/database/scrna-atac/human_brain_10x/human_brain_10x_allobj_activity.rds')
