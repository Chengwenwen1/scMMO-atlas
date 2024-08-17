##BULK RNA data analysis
setwd('/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_mm')
##library
library(data.table)
library(dplyr)
library(tibble)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
gg.theme <-
  theme_bw() +
  theme(#legend.position = "right",
    axis.text.x = element_text(angle=0,hjust=0.5),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(size=12,face="bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'))

gg.theme2 <-
  theme_bw() +
  theme(#legend.position = c(0.1,0.88),
    panel.border =element_rect(color = 'black',fill = NA,size=1),
    axis.text.x=element_text(),
    axis.text.y=element_text(),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(),
    strip.background = element_blank(),
    axis.line = element_blank())
###bulk RNA seq-----------1.DESEQ2----------------------------------------------
###使用的数据格式为原始counts
#exp matirx
bulk.counts <- fread('../../raw_data/bulk_data/ww_raw_bam/bulk_all_sample.count.txt',data.table = F)
gene.info <- bulk.counts[,1:6] 
gene.info$symbol <- mapIds(org.Hs.eg.db, keys = gene.info$Geneid, keytype = 'ENSEMBL',column = 'SYMBOL')
rownames(bulk.counts) <- bulk.counts$Geneid
bulk.counts <- bulk.counts[,-c(1:6)]
colnames(bulk.counts) <- colnames(bulk.counts)%>% substr(14,20)
bulk.counts$symbol <- gene.info$symbol
bulk.counts <- bulk.counts[!is.na(bulk.counts$symbol),]
bulk.counts <- bulk.counts %>% dplyr::group_by(symbol) %>% summarise_all(max) %>% as.data.frame()
rownames(bulk.counts) <- bulk.counts$symbol
bulk.counts$symbol<-NULL
##group information
bulk.meta <- fread('../../raw_data/mm_metadata.csv')
bulk.meta <- bulk.meta[1:21,]
#bulk.meta$Patient <- gsub('-','',bulk.meta$Patient) %>% paste('A',.,sep = '')
bulk.meta <- bulk.meta[which(bulk.meta$Patient %in% colnames(bulk.counts)),]
bulk.meta <- merge(data.frame(Patient = colnames(bulk.counts)), bulk.meta)
bulk.meta$Patientgroup <-paste(bulk.meta$Response,bulk.meta$Patient,bulk.meta$Gender,bulk.meta$Age,bulk.meta$SCISS,sep ='-')
colnames(bulk.counts) <- bulk.meta$Patientgroup
bulk.counts <-bulk.counts[,order(colnames(bulk.counts))]##分组顺序先CR后NR
group<- data.frame(row.names = colnames(bulk.counts),Response=substr(colnames(bulk.counts),1,2))
head(bulk.counts[,1:5])

##bulk.tpm matrix
gene.info <- left_join(data.frame(symbol=rownames(bulk.counts)),gene.info,by='symbol')
bulk.tpm <- t(t(bulk.counts/(gene.info$Length / 1000))/colSums(bulk.counts/(gene.info$Length / 1000)) * 1000000) %>% 
  as.data.frame()
head(bulk.tpm)
###矩阵文件和对应的样本信息已经给到
colnames(bulk.counts)
##其他PCA方法
library(edgeR)
bulk.counts<-bulk.counts[,-c(6,1,5,11,13,15)]
keep <- rowSums(cpm(bulk.counts)>1) >= 2
bulk.counts.keep <- bulk.counts[keep,]
bulk.cpm <- log2(cpm(bulk.counts.keep)+1)
head(bulk.cpm)
bulk.cpm <- as.data.frame(t(bulk.cpm)) %>% na.omit() 
bulk.pca <- prcomp(bulk.cpm,center = T,scale. = T)
head(bulk.pca$rotation)#特征向量，回归系数
head(bulk.pca$sdev)#特征值的开方
head(bulk.pca$x)#样本得分
summary(bulk.pca)##方差贡献率
bulk.summ <- summary(bulk.pca)
bulk.pca.fig <- as.data.frame(bulk.pca$x)
bulk.pca.fig$patient <- substr(rownames(bulk.pca.fig),1,10)
bulk.pca.fig$group <- substr(rownames(bulk.pca.fig),1,2)
xlab<-paste0("PC1(",round(bulk.summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(bulk.summ$importance[2,2]*100,2),"%)")
ggplot(data = bulk.pca.fig,aes(x=PC1,y=PC2,color=group))+
  #stat_ellipse(aes(fill=group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point(size=2.5)+
  scale_color_manual(values = resp.col)+
  geom_text_repel(aes(x=PC1,y=PC2,label=patient),
                  size=4.5,
                  segment.alpha=0.5,
                  point.padding=NA,
                  max.overlaps = 30)+
  labs(x=xlab,y=ylab,color="",title = 'PCA of cancer bulk data')+
  gg.theme2+
  theme(panel.border = element_rect(fill=NA,size = 1),
        legend.position ='none')+
  geom_vline(xintercept=c(0),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept =c(0),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05

##CR
bulk_dim2_down <- rownames(bulk.pca$rotation)[order(bulk.pca$rotation[,1],decreasing = F)][1:100]
##查看pc1的loading基因--------------------NR------------------------------------
bulk_dim2_up <- rownames(bulk.pca$rotation)[order(bulk.pca$rotation[,1],decreasing = T)][1:100]
##dim2的基因热图及富集分析
library(ComplexHeatmap)
bulk.tpm.log<- log2(bulk.tpm[,-c(6,1,5,11,13,15)]+1)#
fwrite(bulk.tpm.log,file = 'bulk.tpm.log.csv',sep = ',',row.names = T)
bulk.tpm.log <- dplyr::select(bulk.tpm.log,rownames(bulk.pca.fig)[order(bulk.pca.fig[,1],decreasing = F)])##按照PC1的样本顺序排列
bulk.htmap.data<- apply(bulk.tpm.log, 1, scale)  %>% as.data.frame() %>% 
  dplyr::select(c(bulk_dim2_up,bulk_dim2_down))%>%t()
colnames(bulk.htmap.data) <- substr(colnames(bulk.tpm.log),4,10)
##列注释
ha = HeatmapAnnotation(
  Response=substr(colnames(bulk.tpm.log),1,2),
  Gender=substr(colnames(bulk.tpm.log),12,12),
  SCISS=substr(colnames(bulk.tpm.log),17,25),
  Age=substr(colnames(bulk.tpm.log),14,15) %>% as.numeric(),
  col = list(Gender = c("F" = "#F1784B", "M" = "#889D5A"),
             Response = c('CR'="#35978F",'NR'="#BF812D"),
             Age=circlize::colorRamp2(c(70, 90), c("#FFFAF0", "#BE6065")),
             SCISS=c("Stage I"='#00DCCD',"Stage II"='#FF9289',"Stage III"='#EDAF00')
  ),
  na_col = "black"
)
gene <- c(postive.nfkb) %>%unique()##选择特定的基因注释
gene_pos <- which(rownames(bulk.htmap.data) %in% gene) 
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,
                                         labels = rownames(bulk.htmap.data)[gene_pos],
                                         labels_gp = gpar(fontsize = 8)))
Heatmap(bulk.htmap.data,
        col =colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(18:28)],
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = ha,
        show_row_names = F,
        row_split = c(rep('100 NR genes',100),rep('CR',100)),
        column_split =  c(rep('CR', 5),rep('NR',6)),
        right_annotation = row_anno,
        heatmap_legend_param = list(
          title= "Zscore", title_position = "topcenter", 
          legend_height=unit(4,"cm"), legend_direction="vertical"),
        column_title = "Loading genes in PC1", 
        column_title_gp = gpar(fontsize = 15, fontface = "bold")
)


##GO和KEGG富集分析--CR----------------------------------------------------------
bulk_dim2_down <- rownames(bulk.pca$rotation)[order(bulk.pca$rotation[,1],decreasing = F)][1:100]
bulk.go.enrich_CR <- clusterProfiler::enrichGO(gene = mapIds(org.Hs.eg.db, 
                                                             keys = bulk_dim2_down,
                                                             keytype = 'SYMBOL',
                                                             column = 'ENTREZID'),
                                               OrgDb = 'org.Hs.eg.db',
                                               ont = 'ALL',
                                               pvalueCutoff = 0.05,
                                               qvalueCutoff = 0.05,
                                               readable = T) 


#my_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(20:35)]
enrichplot::dotplot(bulk.go.enrich_CR,title='bulk cancer data of CR')+
  scale_color_gradient(low = "#C13533", high = "#2C1C40")+
  gg.theme2##CR
View(bulk.go.enrich_CR@result)

##GO和KEGG富集分析--NR----------------------------------------------------------
bulk_dim2_up <- rownames(bulk.pca$rotation)[order(bulk.pca$rotation[,1],decreasing = T)][1:100]
bulk.go.enrich <- clusterProfiler::enrichGO(gene = mapIds(org.Hs.eg.db, 
                                                          keys = bulk_dim2_up,
                                                          keytype = 'SYMBOL',
                                                          column = 'ENTREZID'),
                                            OrgDb = 'org.Hs.eg.db',
                                            ont = 'ALL',
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05,
                                            readable = T) 
# bulk.go.enrich <- clusterProfiler::simplify(
#   bulk.go.enrich,
#   cutoff = 0.7,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang",
#   semData = NULL
# )#GO去冗余
##挑选goterm进行展示
bulk.go.enrich1 <- bulk.go.enrich[
  bulk.go.enrich@result$Description %in% 
    c('positive regulation of cytokine production',
      'regulation of adaptive immune response',
      'interleukin-1 beta production',
      'myeloid cell differentiation',
      'I-kappaB kinase/NF-kappaB signaling',
      'positive regulation of I-kappaB kinase/NF-kappaB signaling',
      'positive regulation of NF-kappaB transcription factor activity',
      'regulation of phagocytosis',
      'interferon-gamma production',
      'positive regulation of interferon-gamma production',
      'regulation of inflammatory response',
      'regulation of tumor necrosis factor production',
      'defense response to virus'),]


bulk.go.enrich1$GeneRatio1 <- lapply(
  bulk.go.enrich1$GeneRatio,function(x){
    as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit(x,'/')[[1]][2])}
) %>% unlist() %>% round(digits = 3)

bulk.go.enrich1 <- bulk.go.enrich1[order(bulk.go.enrich1$GeneRatio1,decreasing=F),]
bulk.go.enrich1$Description <-factor(bulk.go.enrich1$Description,levels =bulk.go.enrich1$Description )
#my_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(15:30)]
ggplot(bulk.go.enrich1,aes(x=GeneRatio1,y=Description))+
  geom_point(aes(size=Count,color=p.adjust))+
  #scale_color_gradientn( colors=rev(my_palette))+
  scale_color_gradient(low = "#C13533", high = "#2C1C40")+
  gg.theme2+
  theme(text = element_text(size = 15),
        panel.border = element_rect(fill=NA,size = 1))+
  labs(x='GeneRatio',y='',title = 'GO terms of NR genes')

postive.nfkb <-c('CLEC7A','CCR7','TRIM8','CARD9','SLC44A2','MAP3K3','FLOT1','TLR2')

#FLOT1,PSEN1,SORL1,JAK2,CLEC7A,CYRIB,FGR,CEBPB,IL4R,CCR7,TLR2,RIOK3,CD14,CARD9,HK1,MNDA,SLC11A1


##差异基因
bulk.counts <- bulk.counts[, ]#!(colnames(bulk.counts) %in% c("CR-A026008",'CR-A015003','NR-A064005'))
group <- data.frame(row.names = colnames(bulk.counts),Response=substr(colnames(bulk.counts),1,2))
dds <- DESeqDataSetFromMatrix(countData=bulk.counts, 
                              colData=group, 
                              design=~Response)
dds <- dds[rowSums(counts(dds))>1,]
dds <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res<- results(dds)
summary(res)
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

### 1.3 volcano plot of mRNA
res2 <- res1 %>% na.omit() %>% rownames_to_column('gene_id') %>%
  mutate(change=factor(ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                              ifelse(log2FoldChange > 1 ,'NR high','CR high'),'NS'),
                       levels=c('NR high','CR high','NS')))

table(res2$change)
ggplot(res2,aes(x=log2FoldChange, y=-log10(padj), color=change))+
  geom_point(size=0.8)+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  theme_classic()+
  ylab('-log10 (P-adjusted)')+xlab('log2 Fold Change')+
  ggtitle('Response: NR vs CR')+
  geom_text_repel(
    data = res2[res2$padj<0.05 & abs(res2$log2FoldChange)>0.25,] ,
    aes(label = gene_id),
    size = 3,max.overlaps = 40,
    segment.color = "black", show.legend = FALSE)+ #添加关注的点的基因名
  geom_vline(xintercept=c(-0.5,0.5),lty='dashed',col="black",lwd=1) +
  geom_hline(yintercept = -log10(0.05),lty='dashed',col="black",lwd=1)+
  theme(
    legend.title = element_blank(),
    #legend.position = c(0.1, 0.7),
    panel.grid = element_blank(),
    legend.text=element_text(size = 12,family ="sans"),
    axis.title= element_text(size = 15,family ="sans"),
    plot.title = element_text(color="black", size=15,face="plain",hjust = 0.5,family ="sans"),
    axis.text= element_text(size=15,family ="sans"),
    axis.ticks = element_blank()
  ) 
##差异基因和pca基因的交集
#gene_pca_deg <-res2$gene_id[res2$change!=c('Not Significant')] [res2$gene_id[res2$change!=c('Not Significant')] %in% c(dim2_top)]
##GO和KEGG富集分析
bulk_deg_nr <- res2$gene_id[res2$change=='NR high']
bulk_deg_cr <- res2$gene_id[res2$change=='CR high']
bulk.go.enrich <- clusterProfiler::enrichGO(gene = mapIds(org.Hs.eg.db, 
                                                          keys = bulk_deg_cr,
                                                          keytype = 'SYMBOL',
                                                          column = 'ENTREZID'),
                                            OrgDb = 'org.Hs.eg.db',
                                            ont = 'ALL',
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05,
                                            readable = T) 

bulk.go.enrich <- clusterProfiler::simplify(
  bulk.go.enrich,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)#GO去冗余
#my_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(20:35)]
enrichplot::dotplot(bulk.go.enrich,showCategory = 10)+
  gg.theme+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(fill=NA,size = 1))+
  scale_color_gradientn( colors=rev(my_palette))+
  labs(title = 'DEG of CR')

View(bulk.go.enrich@result)




##单细胞数据每个病人的plasma cell的pca结果
library(Seurat)
library(edgeR)
mmpc <- subset(mm_pc,celltype %in% c('MMPC'))
mmpc$Patient <- paste(mmpc$Response,mmpc$Patient,mmpc$Gender,mmpc$Age,mmpc$SCISS,sep = '-')
name <- as.data.frame(table(mmpc$Patient))
mmpc <-subset(mmpc,Patient %in% name$Var1[name$Freq>10])
table(mmpc$Patient)
table(mmpc$celltype)
# sc.counts <- AverageExpression(mmpc,group.by = 'Patient',assays = 'RNA',slot = 'counts')
# sc.counts <- sc.counts$RNA
##按照加和处理矩阵数据
sc.counts <- data.frame(gene=rownames(mmpc))
for (i in unique(mmpc$Patient)){
  cellid <- rownames(mmpc@meta.data[mmpc@meta.data$Patient==i,])
  if (length(cellid)>1){
    counts <- apply(mm_pc@assays$RNA@counts[,cellid],1,sum) %>% as.data.frame() %>% rownames_to_column('gene')
  }else{
    counts <- data.frame(gene=rownames(mm_pc),z=mm_pc@assays$RNA@counts[,cellid])
  }
  
  colnames(counts)[2]<-i
  sc.counts <- left_join(sc.counts,counts,"gene")
}
sc.counts <- sc.counts %>% column_to_rownames('gene')
keep <- rowSums(cpm(sc.counts)>1) >= 2
sc.counts.keep <- sc.counts[keep,]
sc.counts.keep <- log2(cpm(sc.counts.keep)+1)
##select samples
#sc.counts.keep <- sc.counts.keep[,-c(3,4,11)]#"CR-A015003-M-80","CR-A171002-M-87"#-c(3,11)
sc.pca <- prcomp(t(sc.counts.keep),center = T,scale. = T)
head(sc.pca$rotation)#特征向量，回归系数
head(sc.pca$sdev)#特征值的开方
head(sc.pca$x)#样本得分
summary(sc.pca)##方差贡献率
sc.summ <- summary(sc.pca)
sc.pca.fig <- as.data.frame(sc.pca$x)
sc.pca.fig$group <- substr(rownames(sc.pca.fig),1,2) 
sc.pca.fig$patient <- substr(rownames(sc.pca.fig),1,10)
sc.xlab<-paste0("PC1(",round(sc.summ$importance[2,1]*100,2),"%)")
sc.ylab<-paste0("PC2(",round(sc.summ$importance[2,2]*100,2),"%)")
ggplot(data = sc.pca.fig,aes(x=PC1,y=PC2,color=group))+
  #stat_ellipse(aes(fill=group),type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point(size=2.5)+
  scale_color_manual(values = resp.col)+
  geom_text_repel(aes(x=PC1,y=PC2,label=patient),
                  size=4.5,
                  segment.alpha=0.5,
                  point.padding=NA,
                  max.overlaps = 30)+
  labs(x=sc.xlab,y=sc.ylab,color=" ",title = 'PCA of cancer cells pseudobulk')+
  gg.theme2+NoLegend()+
  theme(panel.border = element_rect(fill=NA,size = 1))+
  geom_vline(xintercept=c(0),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept =c(0),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05


##查看pc2的loading基因
sc_genes <- sc.pca$rotation %>% as.data.frame()
sc_genes[order(sc_genes[2],decreasing = F),] %>%select(PC2) %>% head()
sc_dim2_up <- rownames(sc.pca$rotation)[order(sc.pca$rotation[,2],decreasing = T)][1:100]
sc_dim2_down <- rownames(sc.pca$rotation)[order(sc.pca$rotation[,2],decreasing = F)][1:100]

##dim2的基因热图及富集分析
library(ComplexHeatmap)
#mmpc <- subset(mmpc,Patient != "CR-A015003-M-80" & Patient !="CR-A171002-M-87" & Patient !='CR-A006003-F')
sc.tpm <- AverageExpression(mmpc,group.by = 'Patient',assays = 'RNA',slot = 'data')
sc.tpm <- as.data.frame(sc.tpm$RNA)
fwrite(sc.tpm,file = './sc.tpm.csv',sep = ',',row.names = T)

sc.tpm.log <- log2(sc.tpm+1)
sc.tpm.log <- dplyr::select(sc.tpm.log,rownames(sc.pca.fig)[order(sc.pca.fig[,2],decreasing = T)])##按照PC2的样本顺序排列
sc.htmap.data<- apply(sc.tpm.log, 1, scale)  %>% as.data.frame() %>% 
  dplyr::select(c(sc_dim2_up,sc_dim2_down))%>%t()
colnames(sc.htmap.data) <- substr(colnames(sc.tpm.log),1,10)
##列注释
ha = HeatmapAnnotation(
  Age=substr(colnames(sc.tpm.log),14,15) %>% as.numeric(),
  SCISS=substr(colnames(sc.tpm.log),17,25),
  Gender=substr(colnames(sc.tpm.log),12,12),
  Response=substr(colnames(sc.tpm.log),1,2),
  col = list(Gender = c("F" = "#F1784B", "M" = "#889D5A"),
             Response = c('CR'="#35978F",'NR'="#BF812D"),
             Age=circlize::colorRamp2(c(70, 90), c("#FFFAF0", "#BE6065")),
             SCISS=c("Stage I"='#00DCCD',"Stage II"='#FF9289',"Stage III"='#EDAF00')
  ),
  na_col = "black"
)
##选择特定的基因注释
gene <- c(sc.niknfgenes) %>% unique()#sc.ikbnfgenes,sc.niknfgenes
gene_pos <- which(rownames(sc.htmap.data)%in%gene) 
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,
                                         labels = rownames(sc.htmap.data)[gene_pos],
                                         labels_gp = gpar(fontsize = 8)))
Heatmap(sc.htmap.data,
        col =colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(1:3,15:30,45:50)],
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = ha,
        show_row_names = F,
        row_split = c(rep(' ',100),rep(' ',100)),
        right_annotation = row_anno,
        heatmap_legend_param = list(
          title= "Zscore", title_position = "topcenter", 
          legend_height=unit(4,"cm"), legend_direction="vertical"),
        column_title = "Loading genes of PC2", 
        column_title_gp = gpar(fontsize = 15, fontface = "bold")
)


##GO和KEGG富集分析
sc_dim2_down <- rownames(sc.pca$rotation)[order(sc.pca$rotation[,2],decreasing = F)][1:100]#100
sc_dim2_up<- rownames(sc.pca$rotation)[order(sc.pca$rotation[,2],decreasing = T)][1:100]
sc.go.enrich <- clusterProfiler::enrichGO(gene = mapIds(org.Hs.eg.db, 
                                                        keys = sc_dim2_up,
                                                        keytype = 'SYMBOL',
                                                        column = 'ENTREZID'),
                                          OrgDb = 'org.Hs.eg.db',
                                          ont = 'ALL',
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff = 0.05,
                                          readable = T) 

View(sc.go.enrich@result)
# sc.go.enrich <- clusterProfiler::simplify(
#   sc.go.enrich,
#   cutoff = 0.7,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang",
#   semData = NULL
# )#GO去冗余
#my_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(25:35)]
enrichplot::dotplot(sc.go.enrich,showCategory = 10)+
  gg.theme+
  scale_color_gradient(low = "#C13533", high = "#2C1C40")+
  gg.theme2+
  labs(title = ' CR terms of pseudobulk data')##CR

enrichplot::dotplot(sc.go.enrich,showCategory = 10)+
  gg.theme+
  scale_color_gradient(low = "#C13533", high = "#2C1C40")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(fill=NA,size = 1))+
  #scale_color_gradientn( colors=rev(my_palette))+
  labs(title = 'GO terms of NR genes')


View(sc.go.enrich@result)
sc.go.enrich@result$geneID[sc.go.enrich@result$Description=='I-kappaB kinase/NF-kappaB signaling']
sc.niknfgenes <-c('BCL10','NOD2','ALPK1','RIPK1','MAP3K14','CFLAR','RELA')




##挑选goterm进行展示
sc.go.enrich.select <- sc.go.enrich[
  sc.go.enrich@result$Description %in% 
    c('I-kappaB kinase/NF-kappaB signaling',
      'regulation of extrinsic apoptotic signaling pathway',
      'positive regulation of interleukin-8 production',
      'NIK/NF-kappaB signaling',
      'CD4-positive, alpha-beta T cell activation',
      'T-helper 1 type immune response',
      'regulation of inflammatory response',
      'positive regulation of I-kappaB kinase/NF-kappaB signaling',
      'positive regulation of NF-kappaB transcription factor activity',
      'positive regulation of T cell receptor signaling pathway',
      'positive regulation of cytokine production involved in immune response',
      'leukocyte mediated immunity'),]


sc.go.enrich.select$GeneRatio1 <- lapply(
  sc.go.enrich.select$GeneRatio,function(x){
    as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit(x,'/')[[1]][2])}
) %>% unlist() %>% round(digits = 3)

sc.go.enrich.select <- sc.go.enrich.select[order(sc.go.enrich.select$GeneRatio1,decreasing=F),]
sc.go.enrich.select$Description <-factor(sc.go.enrich.select$Description,levels =sc.go.enrich.select$Description )
#my_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(15:30)]
ggplot(sc.go.enrich.select,aes(x=GeneRatio1,y=Description))+
  geom_point(aes(size=Count,color=p.adjust))+
  labs(x='GeneRatio',y='')+
  #scale_color_gradientn( colors=rev(my_palette))+
  scale_color_gradient(low = "#C13533", high = "#2C1C40")+
  gg.theme2+
  theme(text = element_text(size = 15),
        panel.border = element_rect(fill=NA,size = 1))+
  labs(x='',title = 'GO terms of NR genes')







###bulk RNA seq-----------2.edgeR----------------------------------------------
#rm(list = ls())
library(edgeR)
library(clusterProfiler)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(org.Hs.eg.db)
library(stringr)
library(stringi)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(data.table)
setwd('/home/chengww/data/project/multiple_myeloma/')
##counts data
counts <- fread('./raw_data/bulk_data/ww_raw_bam/bulk_all_sample.count.txt')
gene.inf <- counts[,c(1:6)]
gene.inf$symbol <- mapIds(org.Hs.eg.db, keys = gene.inf$Geneid, keytype = 'ENSEMBL',column = 'SYMBOL')
head(gene.inf)
counts <- merge(counts[,-c(2:6)],gene.inf[,c(1,7)],by='Geneid')
counts <- counts[,-1]
counts <- counts %>% group_by(symbol) %>% summarise_all(max) %>% as.data.frame() %>%na.omit()
counts <- column_to_rownames(counts,'symbol')
colnames(counts) <- colnames(counts) %>% substr(14,20)
head(counts) ##正式生成raw_counts matrix
##filter counts
keep <- rowSums(cpm(counts)>1) >= 2
counts2 <- counts[keep,]
cpm <- log2(cpm(counts2)+1)
head(cpm)
##获取样本分组信息
bulk.meta <- fread('./raw_data/mm_bulk.metadata.csv')
head(bulk.meta)
bulk.meta$Patient <- paste('A',gsub('-','',bulk.meta$Patient),sep = '')
bulk.meta <-left_join(data.frame(Patient=colnames(cpm)),bulk.meta)
bulk.meta$newgroup <- paste(bulk.meta$Response,bulk.meta$Patient,sep = '-')
colnames(cpm) <- bulk.meta$newgroup
cpm <- cpm[,order(colnames(cpm),decreasing = T)]#按照contorl组在前进行排序
colnames(counts2) <- bulk.meta$newgroup
counts2 <- counts2[,order(colnames(counts2),decreasing = T)]

##组间比较
group <- factor(substr(colnames(cpm),1,2),levels = c('NR','CR'))
data <- data.frame(expression=c(cpm),
                   sample=rep(colnames(cpm),each=nrow(cpm)))
ggplot(data,aes(x=sample,y=expression,fill=sample))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  xlab(NULL)+ylab('Log2(CPM+1)')
##PCA结果展示
dat1 <- cpm
dat1 <- as.data.frame(t(dat1)) %>% na.omit() 
dat1$group <- substr(rownames(dat1),1,2) %>% as.factor()
dat1_pca <- PCA(dat1[,-ncol(dat1)],graph = F)
ind <- get_pca_ind(dat1_pca)
fviz_pca_ind(dat1_pca, col.ind = 'cos2',
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T,
             legend.title='Groups')+theme_bw()
fviz_pca_ind(dat1_pca,
             # 只显示点而不显示文本，默认都显示
             #geom.ind = "point",
             # 设定分类种类
             col.ind = dat1$group,
             # 设定颜色
             palette = c("#FC4E07","#00AFBB"),
             # 添加椭圆 Concentration ellipses
             # addEllipses = TRUE,
             legend.title = "Groups",
             repel = T
)+theme_bw()
##查看变量贡献程度
dim2 <- dat1_pca[["var"]][["contrib"]] [,'Dim.2']
pcagenes <- dim2[order(dim2,decreasing = T)][1:2000]

###删除一个样本
dat1_filter <- dat1[!(rownames(dat1) %in% c("CR-A029001")),colnames(dat1) %in% c(sig.genes,'group')]
dat1_filter_pca <- PCA(dat1_filter[,-ncol(dat1_filter)],graph = F)
fviz_pca_ind(dat1_filter_pca,
             # 只显示点而不显示文本，默认都显示
             #geom.ind = "point",
             # 设定分类种类
             col.ind = dat1_filter$group,
             # 设定颜色
             palette = c("#FC4E07","#00AFBB"),
             # 添加椭圆 Concentration ellipses
             addEllipses = T,
             legend.title = "Groups",
             repel = T
)+theme_bw()


##edgeR的差异分析
group;colnames(counts2)
#表达值预处理和标准化
deg <-DGEList(counts = counts2,
              group=group)
deg$samples$lib.size
deg <-calcNormFactors(deg)#标准化

##计算差异倍数
design <- model.matrix(~group)#分组信息中一定要是对照组在前，处理组在后
#计算线性模型参数
deg <- estimateGLMCommonDisp(deg,design)
deg  <- estimateGLMTrendedDisp(deg,design)
deg <- estimateGLMTagwiseDisp(deg,design)
##拟合线性模型
fit <-glmFit(deg,design )
lrt <- glmLRT(fit,contrast = c(1,-1))
degs <- as.data.frame(topTags(lrt,n=nrow(deg)))
head(degs)
dim(degs)
degs$regulates <- ifelse(degs$logFC>1 & degs$PValue<0.05,
                         'up',ifelse(degs$logFC< c(-1) & degs$PValue<0.05,'down','normal'))
table(degs$regulates)
##绘制火山图
ggplot(degs,aes(x=logFC,y=-log10(PValue),color=regulates))+
  geom_point(alpha=0.5,size=1.8)+
  theme_set(theme_set(theme_bw(base_size=20))) + 
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c('blue','black','red'))+theme_bw()
##绘制热图
sig_genes <- degs[degs$regulates!='normal',] %>% rownames()
dat2 <- cpm[match(sig_genes,rownames(cpm)),]
head(dat2[,1:4])
dim(dat2)
group1 <- data.frame(row.names = colnames(cpm),group=group) 
library(pheatmap)
pheatmap(dat2,scale = 'row',show_colnames = T,show_rownames = F,
         cluster_cols = T,annotation_col = group1,
         main = "edgeR's DEGS")

##删除其中一个样本后差异分析#"CR-A029001"-----------
counts2_filter <- select(counts2,-c("CR-A029001"))
##edgeR的差异分析
colnames(counts2_filter)
new_group<- substr(colnames(counts2_filter),1,2) %>% factor(levels = c('NR','CR'))
#表达值预处理和标准化
deg <-DGEList(counts = counts2_filter,
              group=new_group)
deg$samples$lib.size
deg <-calcNormFactors(deg)#标准化

##计算差异倍数
design <- model.matrix(~new_group)#分组信息中一定要是对照组在前，处理组在后
#计算线性模型参数
deg <- estimateGLMCommonDisp(deg,design)
deg  <- estimateGLMTrendedDisp(deg,design)
deg <- estimateGLMTagwiseDisp(deg,design)
##拟合线性模型
fit <-glmFit(deg,design )
lrt <- glmLRT(fit,contrast = c(1,-1))
degs <- as.data.frame(topTags(lrt,n=nrow(deg)))
head(degs)
dim(degs)
degs$regulates <- ifelse(degs$logFC>1 & degs$PValue<0.05,
                         'up',ifelse(degs$logFC< c(-1) & degs$PValue<0.05,'down','normal'))
table(degs$regulates)
##绘制火山图
ggplot(degs,aes(x=logFC,y=-log10(PValue),color=regulates))+
  geom_point(alpha=0.5,size=1.8)+
  theme_set(theme_set(theme_bw(base_size=20))) + 
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c('blue','black','red'))+theme_bw()
##绘制热图
sig_genes <- degs[degs$regulates!='normal',] %>% rownames()
cpm2 <- cpm %>% as.data.frame() %>% select(-c("CR-A029001"))
dat2 <- cpm2[match(sig_genes,rownames(cpm2)),]
head(dat2[,1:4])
dim(dat2)
group1 <- data.frame(row.names = colnames(cpm2),group=new_group) 
library(pheatmap)
pheatmap(dat2,scale = 'row',show_colnames = T,show_rownames = F,
         cluster_cols = T,annotation_col = group1,
         main = "edgeR's DEGS")



###bulk RNA seq-----------3.limma ----------------------------------------------
### 使用的数据格式为log2（bulk.tpm+1）
#exp matirx
exp <- fread('../raw_data/bulk_data/ww_raw_bam/bulk_all_sample.count.txt',data.table = F)
gene.info <- exp[,1:6] 
gene.info$symbol <- mapIds(org.Hs.eg.db, keys = gene.info$Geneid, keytype = 'ENSEMBL',column = 'SYMBOL')

rownames(exp) <- exp$Geneid
bulk.counts <- exp[,-c(1:6)]
colnames(bulk.counts) <- colnames(bulk.counts)%>% substr(14,20)
exp.bulk.tpm<- t(t(bulk.counts/(gene.info$Length / 1000))/colSums(bulk.counts/(gene.info$Length / 1000)) * 1000000) %>% 
  as.data.frame()
identical(rownames(exp.bulk.tpm), gene.info$Geneid)

exp.bulk.tpm$symbol <- gene.info$symbol
exp.bulk.tpm <- exp.bulk.tpm[!is.na(exp.bulk.tpm$symbol),]
exp.bulk.tpm <- exp.bulk.tpm %>% dplyr::group_by(symbol) %>% summarise_all(max) %>% as.data.frame()
rownames(exp.bulk.tpm) <- exp.bulk.tpm$symbol
exp.bulk.tpm$symbol<-NULL
exp.bulk.tpm.log<-log2(exp.bulk.tpm + 1)
##group information
bulk.meta <- fread('../raw_data/mm_bulk.metadata.csv')
bulk.meta$Patient <- gsub('-','',bulk.meta$Patient) %>% paste('A',.,sep = '')
bulk.meta <- bulk.meta[which(bulk.meta$Patient %in% colnames(exp.bulk.tpm.log)),]
bulk.meta<-merge(data.frame(Patient = colnames(exp.bulk.tpm.log)), bulk.meta)


### DEG analysis for gene expression
grp.response<-if_else(bulk.meta$Response == 'CR', 1, 0)
mRNA.limma<-lmFit(object = exp.bulk.tpm.log, design = model.matrix(~grp.response)) %>%
  eBayes() %>% topTable(coef = 2, number = Inf, adjust.method = 'BH')
sig.mRNA<-mRNA.limma %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0.05) %>% rownames()