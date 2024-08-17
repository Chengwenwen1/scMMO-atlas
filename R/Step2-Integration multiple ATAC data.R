##ATAC integration--------------------------------------------------------------
DefaultAssay(ISSAAC_mCortex)<-'ATAC'
DefaultAssay(mouse_brain)<-'ATAC'
DefaultAssay(GSE126074_AdBrainCortex)<-'ATAC'
DefaultAssay(mouse_brain_AD)<-'ATAC'
#DefaultAssay(GSE126074_P0_BrainCortex)<-'ATAC'
DefaultAssay(mouse_brain_e18)<-'ATAC'

##Selection of cell IDs after RNA data filtering-------------
ISSAAC_mCortex$raw_cellid <- rownames(ISSAAC_mCortex@meta.data)
ISSAAC_mCortex <- subset(ISSAAC_mCortex,raw_cellid %in%
                        brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='ISSAAC_mCortex'])

GSE126074_AdBrainCortex$raw_cellid <- rownames(GSE126074_AdBrainCortex@meta.data)
GSE126074_AdBrainCortex <- subset(GSE126074_AdBrainCortex,raw_cellid %in%
                                         brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='GSE126074_AdBrainCortex'])

mouse_brain$raw_cellid <- rownames(mouse_brain@meta.data)
mouse_brain <- subset(mouse_brain,raw_cellid %in%
                             brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='mouse_brain_10x'])

mouse_brain_e18$raw_cellid <- rownames(mouse_brain_e18@meta.data)
mouse_brain_e18 <- subset(mouse_brain_e18,raw_cellid%in%
                                 brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='mouse_brain_10x_e18'])

mouse_brain_AD$raw_cellid <- rownames(mouse_brain_AD@meta.data)
mouse_brain_AD <- subset(mouse_brain_AD,raw_cellid%in%
                                brain.integrated1$raw_cellid[brain.integrated1$orig.ident=='mouse_brain_10x_AD'])

###Selection of cell IDs after RNA data filtering-------------
ISSAAC_mCortex.peak <-data.frame( V1=rownames(ISSAAC_mCortex@assays$ATAC@counts) )
ISSAAC_mCortex.peak <- tidyr::separate(ISSAAC_mCortex.peak,V1, into= c("chr","start",'end'),sep= "\\-") 

GSE126074_AdBrainCortex.peak <- data.frame( V1=rownames(GSE126074_AdBrainCortex@assays$ATAC@counts))
GSE126074_AdBrainCortex.peak <- tidyr::separate(GSE126074_AdBrainCortex.peak,V1, into= c("chr","start",'end'),sep= "\\-") 

mouse_brain.peak <-data.frame( V1=rownames(mouse_brain@assays$ATAC@counts) )
mouse_brain.peak <- tidyr::separate(mouse_brain.peak,V1, into= c("chr","start",'end'),sep= "\\-") 

mouse_brain_e18.peak <-data.frame( V1=rownames(mouse_brain_e18@assays$ATAC@counts) )
mouse_brain_e18.peak <- tidyr::separate(mouse_brain_e18.peak,V1, into= c("chr","start",'end'),sep= "\\-") 

mouse_brain_AD.peak <-data.frame( V1=rownames(mouse_brain_AD@assays$ATAC@counts) )
mouse_brain_AD.peak$V1 <-gsub('GRCh38-chr21','GRCh38:chr21',mouse_brain_AD.peak$V1)
mouse_brain_AD.peak <- tidyr::separate(mouse_brain_AD.peak,V1, into= c("chr","start",'end'),sep= "\\-") 

# convert to genomic ranges
ISSAAC_mCortex.peak <- makeGRangesFromDataFrame(ISSAAC_mCortex.peak)
GSE126074_AdBrainCortex.peak <- makeGRangesFromDataFrame(GSE126074_AdBrainCortex.peak)
mouse_brain.peak <- makeGRangesFromDataFrame(mouse_brain.peak)
mouse_brain_e18.peak <- makeGRangesFromDataFrame(mouse_brain_e18.peak)
mouse_brain_AD.peak <- makeGRangesFromDataFrame(mouse_brain_AD.peak)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(ISSAAC_mCortex.peak, GSE126074_AdBrainCortex.peak, 
                               mouse_brain.peak, mouse_brain_e18.peak,mouse_brain_AD.peak))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

##fragement
atac.frag <- list()
atac.frag[['ISSAAC_mCortex']] <- CreateFragmentObject('./scrna-atac/E-MTAB-11264/fragments_all.tsv.gz',
                                                               cells=colnames(ISSAAC_mCortex))
atac.frag[['GSE126074_AdBrainCortex']] <- CreateFragmentObject('./scrna-atac/GSE126074_AdBrainCortex/fragments.sort.bed.gz',
                                                               cells=colnames(GSE126074_AdBrainCortex))
atac.frag[['mouse_brain']] <- CreateFragmentObject('./scrna-atac/mouse_brain_10x/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz',
                                                   cells=colnames(mouse_brain))
atac.frag[['mouse_brain_AD']] <- CreateFragmentObject('./scrna-atac/mouse_brain_10x_AD/Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_atac_fragments.tsv.gz',
                                                      cells=colnames(mouse_brain_AD))
atac.frag[['mouse_brain_e18']] <- CreateFragmentObject('./scrna-atac/mouse_brain_10x_e18/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz',
                                                       cells=colnames(mouse_brain_e18))
counts <- list();
atac.assay <- list()
for (i in c('ISSAAC_mCortex','GSE126074_AdBrainCortex','mouse_brain','mouse_brain_e18','mouse_brain_AD')){#
  # quantify multiome peaks in the scATAC-seq dataset
  counts[[i]] <- FeatureMatrix(
    fragments = atac.frag[[i]],
    features = combined.peaks,
    cells = get(i) %>% colnames() %>% unique()
  )
  
  # create object
  atac.assay[[i]] <- CreateChromatinAssay(
    counts = counts[[i]],
    #min.features = 100,
    fragments = atac.frag[[i]]
  )
  
  z <- paste(i,'.atac',sep ='')
  
  assign(z, CreateSeuratObject(counts = atac.assay[[i]], assay = "ATAC"))
  #ad.atac <- subset(ad.atac, nCount_ATAC > 2000 & nCount_ATAC < 30000)
  
  # compute LSI
  assign(z, FindTopFeatures(get(z), min.cutoff = 10))
  assign(z, RunTFIDF(get(z)))
  assign(z, RunSVD(get(z)))
}
saveRDS(atac.assay,'./scrna-atac/atac.assay.Rds')

# merge
brain.combined <- merge(ISSAAC_mCortex.atac,
                        c(GSE126074_AdBrainCortex.atac,
                          mouse_brain.atac,
                          mouse_brain_e18.atac,
                          mouse_brain_AD.atac
                          #GSE126074_P0_BrainCortex.atac,
                        )#,
                        #add.cell.ids = c("ISSAAC_mCortex", "GSE126074_AdBrainCortex",'mouse_brain_10x' ,"mouse_brain_10x_e18", "mouse_brain_10x_AD")
)
brain.integrated1@meta.data[!duplicated(brain.integrated1@meta.data$orig.ident),]
brain.combined@meta.data[!duplicated(brain.combined@meta.data$orig.ident),]
# process the combined dataset
brain.combined <- FindTopFeatures(brain.combined, min.cutoff = 10)
brain.combined <- RunTFIDF(brain.combined)
brain.combined <- RunSVD(brain.combined)
brain.combined <- RunUMAP(brain.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(brain.combined, group.by = "orig.ident")
p1


integration.anchors <- FindIntegrationAnchors(
  object.list = list(ISSAAC_mCortex.atac, 
                     GSE126074_AdBrainCortex.atac,
                     mouse_brain.atac,
                     mouse_brain_e18.atac,
                     mouse_brain_AD.atac
                     #GSE126074_P0_BrainCortex.atac,
  ),
  anchor.features = rownames(ISSAAC_mCortex.atac),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = brain.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "celltype",cols = type.col,label = T)
p2
p1+p2

saveRDS(integrated,'./scrna-atac/brain.integrated_atac.rds')
