
brain.integrated1 <- readRDS('./scrna-atac/brain.integrated_rna.rds')
integrated<-readRDS('./scrna-atac/brain.integrated_atac.rds')

##merger RNA and ATAC data
brain.integrated_all <-brain.integrated1
brain.integrated_all[['ATAC']]<-integrated[['ATAC']]
brain.integrated_all@reductions$integrated_lsi <- integrated@reductions$integrated_lsi
brain.integrated_all@reductions$umap.atac <- integrated@reductions$umap

##WNN integrated--------------------------------------------
brain.integrated_all <- FindMultiModalNeighbors(brain.integrated_all, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:30, 2:35))
brain.integrated_all <- RunUMAP(brain.integrated_all, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#plot

pp1 <- DimPlot(brain.integrated_all, reduction = "umap", 
               group.by = "celltype_new",
               label = F, 
               label.size = 4, 
               repel = TRUE,
               cols = type.col) + ggtitle("RNA")+NoLegend()
pp2 <- DimPlot(brain.integrated_all, reduction = "umap.atac", 
               group.by = "celltype_new", 
               label = F, 
               label.size = 4, 
               repel = TRUE,
               cols = type.col) + ggtitle("ATAC")+NoLegend()
pp1/pp2
pp3 <- DimPlot(brain.integrated_all, reduction = "wnn.umap", 
               group.by = "celltype_new", 
               label = TRUE, 
               label.size = 4, 
               repel = TRUE,
               cols = type.col) + ggtitle("WNN")
pp1 + pp2 + pp3 & theme(plot.title = element_text(hjust = 0.5))

saveRDS(brain.integrated_all,'./scrna-atac/brain.integrated_rna_atac.rds')