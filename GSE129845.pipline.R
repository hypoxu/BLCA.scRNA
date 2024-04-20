###
### GSE129845
### 
library(Seurat)
setwd('./Bladder.scRNAseq.GEO//GSE129845.Homo.Mice')
Sample1 <- Read10X('GSM3723357_Sample1')
Sample2 <- Read10X('GSM3723358_Sample2')
Sample3 <- Read10X('GSM3723359_Sample3')
Sample4 <- Read10X('GSM3723360_Sample4')  
Sample5 <- Read10X('GSM3723361_Sample5')  
#
Seurat1 <- CreateSeuratObject(counts = Sample1, project = "H1", min.cells = 1, min.features = 200)
Seurat2 <- CreateSeuratObject(counts = Sample2, project = "H2", min.cells = 1, min.features = 200)
Seurat3 <- CreateSeuratObject(counts = Sample3, project = "H3", min.cells = 1, min.features = 200)
Seurat4 <- CreateSeuratObject(counts = Sample4, project = "M1", min.cells = 1, min.features = 200)
Seurat5 <- CreateSeuratObject(counts = Sample5, project = "M2", min.cells = 1, min.features = 200)

###
Seurat1[["percent.mt"]] <- PercentageFeatureSet(Seurat1, pattern = "^MT-")
Seurat2[["percent.mt"]] <- PercentageFeatureSet(Seurat2, pattern = "^MT-")
Seurat3[["percent.mt"]] <- PercentageFeatureSet(Seurat3, pattern = "^MT-")
Seurat4[["percent.mt"]] <- PercentageFeatureSet(Seurat4, pattern = "^mt-")
Seurat5[["percent.mt"]] <- PercentageFeatureSet(Seurat5, pattern = "^mt-")
#
VlnPlot(Seurat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(Seurat5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

#
plot1 <- FeatureScatter(Seurat1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
FeatureScatter(Seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
FeatureScatter(Seurat2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
FeatureScatter(Seurat3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  


##  
SimpleRunSeuratPipline <- function(xx, nfeatures=3000, ndims=1:50){
   require(Seurat); require(dplyr)
   xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000) 
   xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = nfeatures)
   xx <- ScaleData(xx, features = rownames(xx))
   xx <- RunPCA(xx, features = VariableFeatures(object = xx), npcs=100)
   xx <- FindNeighbors(xx, dims = ndims) %>% FindClusters(resolution = c(1, 0.8, 0.6)) %>% RunUMAP(dims = ndims) %>% RunTSNE(dims = ndims)
   ## Cell Cycle analysis
   xx <- CellCycleScoring(xx, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
   xx
}
Homo01 <- SimpleRunSeuratPipline(Seurat1, ndims=1:25)
Homo02 <- SimpleRunSeuratPipline(Seurat2, ndims=1:25)
Homo03 <- SimpleRunSeuratPipline(Seurat3, ndims=1:25)
Mice01 <- SimpleRunSeuratPipline(Seurat4, ndims=1:25)
Mice02 <- SimpleRunSeuratPipline(Seurat5, ndims=1:25)


###
DoubletFinderWrapper <- function(seu_kidney, pcs=1:50, clusters='seurat_clusters', gt=FALSE){
   require(DoubletFinder)
   print('pK Identification')
   if (gt==FALSE){
   ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
     sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = pcs, sct = FALSE)
     sweep.stats_kidney    <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
     bcmvn_kidney          <- find.pK(sweep.stats_kidney)
   }else{
   ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
     sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = pcs, sct = FALSE)
     gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]    ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
     sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
     bcmvn_kidney <- find.pK(sweep.stats_kidney)
   }
   print ('Homotypic Doublet Proportion Estimate')
   homotypic.prop <- modelHomotypic(seu_kidney@meta.data[ ,clusters])  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
   nExp_poi       <- round(0.075*nrow(seu_kidney@meta.data))           ## Assuming 7.5% doublet formation rate - tailor for your dataset
   nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
   
   print ('Run DoubletFinder with varying classification stringencies')
   seu_kidney <- doubletFinder_v3(seu_kidney, PCs = pcs, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
   tmp        <- grep('pANN', colnames(seu_kidney@meta.data), value=T)
   seu_kidney <- doubletFinder_v3(seu_kidney, PCs = pcs, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = tmp, sct = FALSE)
   colnames(seu_kidney@meta.data)[(ncol(seu_kidney@meta.data)-2):ncol(seu_kidney@meta.data)] <- c('pANN', 'DF1', 'DF2')
   seu_kidney
}
#
Homo01 <- DoubletFinderWrapper(Homo01, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Homo02 <- DoubletFinderWrapper(Homo02, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Homo03 <- DoubletFinderWrapper(Homo03, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Mice01 <- DoubletFinderWrapper(Mice01, pcs=1:25, clusters='seurat_clusters', gt=FALSE)
Mice02 <- DoubletFinderWrapper(Mice02, pcs=1:25, clusters='seurat_clusters', gt=FALSE)

save(Homo01,Homo02,Homo03, Mice01,Mice02, file=file.path('Data.Analysis.Step01.Raw.20020719.RData'))
 
Homo01QC <- subset(Homo01, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Homo02QC <- subset(Homo02, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Homo03QC <- subset(Homo03, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Mice01QC <- subset(Mice01, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')
Mice02QC <- subset(Mice02, nFeature_RNA>=500 & nFeature_RNA<=4000 & percent.mt<=15 & DF2=='Singlet')

Homo01QC <- SimpleRunSeuratPipline(Homo01QC, ndims=1:20)
Homo02QC <- SimpleRunSeuratPipline(Homo02QC, ndims=1:20)
Homo03QC <- SimpleRunSeuratPipline(Homo03QC, ndims=1:20)
Mice01QC <- SimpleRunSeuratPipline(Mice01QC, ndims=1:20)
Mice02QC <- SimpleRunSeuratPipline(Mice02QC, ndims=1:20)

#
DimPlot(Homo01QC, label=T) + VlnPlot(Homo01QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Homo02QC, label=T) + VlnPlot(Homo02QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Homo03QC, label=T) + VlnPlot(Homo03QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Mice01QC, label=T) + VlnPlot(Mice01QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(Mice02QC, label=T) + VlnPlot(Mice02QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

 
###
MarkerJASN <- data.frame(Genes=c("KRT18","KRT19","KRT5","KRT17","KRT13","KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2",
                                 "VIM","S100A4","COL3A1","COL1A1","COL1A2","ACTA2","TAGLN","DES","CNN1","ACTG2","TPM2",
								 "SELE","PECAM1","VCAM1","CDH5","LYZ","MS4A7","CD14","CD209","CD3D","CD3E","MZB1","CD79A","GPM6A"),
                       Markers=c('Epi','Epi','Bas','Bas','BasInter','Umbr','UmbrInter','UmbrInter','UmbrInter','Umbr','Umbr',
					             'Interstitial','Fibr','Fibr','Fibr','Fibr','myoFibr','Muscle','Muscle','Muscle','Muscle','Muscle',
								 'Endo','Endo','Endo','Endo','Mono','Mono','Mono','Dendritic','Tcell','Tcell','Bcell','Bcell','Neuron'))
DoHeatmap(Homo01QC,features=MarkerJASN$Genes)
DoHeatmap(Homo02QC,features=MarkerJASN$Genes)
DoHeatmap(Homo03QC,features=MarkerJASN$Genes)
DoHeatmap(Mice01QC,features=stringr::str_to_title(MarkerJASN$Genes))
DoHeatmap(Mice02QC,features=stringr::str_to_title(MarkerJASN$Genes))


Homo01QC@meta.data$CC1 <- dplyr::recode(Homo01QC@meta.data$seurat_clusters, '0'='Fibr', '1'='TCells','2'='myoFibr', '3'='Endo', '4'='Mono','5'='Muscle','6'='Epi')
Homo02QC@meta.data$CC1 <- dplyr::recode(Homo02QC@meta.data$seurat_clusters, '0'='Fibr', '1'='Fibr','2'='myoFibr', '3'='Fibr', '4'='Bas','5'='Umbr','6'='Fibr','7'='Umbr', '8'='Fibr','9'='Mono', '10'='Mono','11'='TCells','12'='Endo')
Homo03QC@meta.data$CC1 <- dplyr::recode(Homo03QC@meta.data$seurat_clusters, '0'='Bas', '1'='Umbr','2'='UmbrInter', '3'='UmbrInter', '4'='UmbrInter',
                                            '5'='Muscle','6'='Fibr','7'='Fibr','8'='Endo','9'='TCells', '10'='Mono','11'='myoFibr','12'='Umbr','13'='BCell','14'='Fibr')
Mice01QC@meta.data$CC1 <- dplyr::recode(Mice01QC@meta.data$seurat_clusters, '0'='Bas', '1'='Umbr','2'='UmbrInter', '3'='Fibr', '4'='Bas','5'='Fibr','6'='UmbrInter','7'='myoFibr','8'='Umbr',
                                            '9'='Fibr', '10'='Mono','11'='unknown','12'='Neuron','13'='unknown','14'='Muscle','15'='Endo')
Mice02QC@meta.data$CC1 <- dplyr::recode(Mice02QC@meta.data$seurat_clusters, '0'='UmbrInter', '1'='UmbrInter','2'='Fibr', '3'='Bas', '4'='UmbrInter','5'='Bas',
                                            '6'='Fibr','7'='myoFibr','8'='Fibr','9'='Mono', '10'='Bas','11'='Fibr','12'='Umbr','13'='TCells','14'='Muscle','15'='Endo')

DimPlot(Homo01QC, label=T, group.by='CC1') + VlnPlot(Homo01QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='CC1') 
DimPlot(Homo02QC, label=T, group.by='CC1') + VlnPlot(Homo02QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='CC1') 
DimPlot(Homo03QC, label=T, group.by='CC1') + VlnPlot(Homo03QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='CC1') 
DimPlot(Mice01QC, label=T, group.by='CC1') + VlnPlot(Mice01QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='CC1') 
DimPlot(Mice02QC, label=T, group.by='CC1') + VlnPlot(Mice02QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='CC1') 


### 
HomoAll <- merge(Homo01QC, list(Homo02QC, Homo03QC))
MiceAll <- merge(Mice01QC, Mice02QC)
HomoAll <- SimpleRunSeuratPipline(HomoAll, ndims=1:20)
MiceAll <- SimpleRunSeuratPipline(MiceAll, ndims=1:20)

DimPlot(HomoAll, label=T) + VlnPlot(HomoAll, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(MiceAll, label=T) + VlnPlot(MiceAll, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
# 
DimPlot(HomoAll, label=T, group.by='orig.ident') 
DimPlot(HomoAll, label=T, group.by='CC1') 
DimPlot(HomoAll, label=T, group.by='seurat_clusters') 
VlnPlot(HomoAll, group.by='CC1', features=c('PIN1'))
# 
DimPlot(MiceAll, label=T, group.by='orig.ident') 
DimPlot(MiceAll, label=T, group.by='CC1') 
DimPlot(MiceAll, label=T, group.by='seurat_clusters') 
FeaturePlot(MiceAll, reduction='umap', order=T, features=c('Pin1'))
VlnPlot(MiceAll, group.by='CC1', features=c('Pin1'))


###

{## HMY整合数据
   library(harmony); library(SeuratWrappers); library(patchwork)
   HomoAll_HMY <- RunHarmony(HomoAll, group.by.vars = c("orig.ident"), assays='RNA')
   #ElbowPlot(HomoAll_HMY, ndims = 100, reduction='harmony')
   HomoAll_HMY <- RunUMAP(HomoAll_HMY, reduction = "harmony", dims = 1:25)
   HomoAll_HMY <- RunTSNE(HomoAll_HMY, reduction = "harmony", dims = 1:25)
   HomoAll_HMY <- FindNeighbors(HomoAll_HMY, reduction = "harmony", dims = 1:25) %>% FindClusters(reduction = "harmony", resolution = c(1, 0.8, 0.6, 0.4))
   
   DimPlot(HomoAll_HMY,  label=T)
   DimPlot(HomoAll_HMY,  group.by='orig.ident', label=T)
   DimPlot(HomoAll_HMY,  group.by='CC1', label=T)
   DimPlot(HomoAll_HMY,  group.by='ToSciTum2a', label=T)
   DimPlot(HomoAll_HMY,  group.by='ToSciNor2a', label=T)
   DimPlot(HomoAll_HMY, split.by='orig.ident', group.by='ToSciTum2a', label=T)
   
   FeaturePlot(HomoAll_HMY, reduction='umap', order=T, features=c('CA9','EPCAM','PECAM1','PTPRC'))
   FeaturePlot(HomoAll_HMY, reduction='umap', order=T, features=c('PIN1'))
   VlnPlot(HomoAll_HMY, group.by='ToSciTum2a', features=c('CA9','EPCAM','PECAM1','PTPRC'))
   
   FeaturePlot(HomoAll_HMY, reduction='umap', order=T, features=c('CA9','EPCAM','PECAM1','PTPRC')) +  scale_color_viridis_c()
   FeaturePlot(HomoAll_HMY, reduction='umap', order=T, features=c('nCount_RNA','nFeature_RNA','percent.mt')) +  scale_color_viridis_c()
   VlnPlot(HomoAll_HMY, group.by='ToSciTum2a', features=c('CA9','EPCAM','PECAM1','PTPRC'))
#   save(HomoAll_HMY, file='./Data.Analysis.2022/Data.scRNA.AllIn.01.HMY.RData')
}
  


### 
MarkerJASN <- data.frame(Genes=c("KRT18","KRT19","KRT5","KRT17","KRT13","KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2",
                                 "VIM","S100A4","COL3A1","COL1A1","COL1A2","ACTA2","TAGLN","DES","CNN1","ACTG2","TPM2",
								 "SELE","PECAM1","VCAM1","CDH5","LYZ","MS4A7","CD14","CD209","CD3D","CD3E","MZB1","CD79A","GPM6A"),
                       Markers=c('Epi','Epi','Bas','Bas','BasInter','Umbr','UmbrInter','UmbrInter','UmbrInter','Umbr','Umbr',
					             'Interstitial','Fibr','Fibr','Fibr','Fibr','myoFibr','Muscle','Muscle','Muscle','Muscle','Muscle',
								 'Endo','Endo','Endo','Endo','Mono','Mono','Mono','Dendritic','Tcell','Tcell','Bcell','Bcell','Neuron'))
DoHeatmap(MiceAll,features=stringr::str_to_title(MarkerJASN$Genes))
##
tmpMarkerMice <- c('Krt5','Krt18','Krt19','Krt20','Upk1a','Upk1b','Upk2','Upk3a','Upk3b','Col3a1','Col1a1','Col1a2','Vim',
                   'Acta2','Lyz2','Ms4a7','Cd14','Gpm6a','Tagln','Des','Cnn1','Actg2','Tpm2','Cd209a','Cd3d','Cd3e','Pecam1','Cdh5')
DoHeatmap(MiceAll,features=tmpMarkerMice)
#
MiceAll@meta.data$CC2 <- dplyr::recode(MiceAll@meta.data$seurat_clusters, '0'='Luminal', '1'='Basal','2'='Basal', '3'='Basal', '4'='Fibroblast','5'='Fibroblast',
                                            '6'='Basal','7'='Basal','8'='Bas.Lum.Mixed','9'='myoFibroblast', '10'='Bas.Lum.Mixed','11'='Fibroblast','12'='myoFibroblast','13'='Monocyte',
											'14'='DendriticCells','15'='Neurone','16'='TCells','17'='Muscle','18'='Endo')
##
VlnPlot(MiceAll, group.by='CC2', features='Pin1')


##
DimPlot(HomoAll_HMY, label=T )
DimPlot(HomoAll_HMY, label=T, group.by='CC1')
DoHeatmap(HomoAll_HMY,features=MarkerJASN$Genes)

tmpMarkerHomo <- c('KRT5','KRT17','KRT13','KRT18','KRT19','COL3A1','COL1A1','COL1A2','S100A4','ACTA2','VIM','ACTG2','DES','TPM2','TAGLN','CNN1','LYZ','MS4A7',
                    'CD14','SELE','VCAM1','PECAM1','CDH5','TNNT1','CD3D','CD3E','UPK3A','UPK1A','UPK1B','UPK2','PIGR','ADRA2A','HRH2','AVPR1A','MZB1','CD79A')
DoHeatmap(HomoAll_HMY,features=tmpMarkerHomo)

HomoAll_HMY@meta.data$CC2 <- dplyr::recode(HomoAll_HMY@meta.data$seurat_clusters, '0'='Intermediate', '1'='Basal','2'='Fibroblast', '3'='Luminal', '4'='Fibroblast','5'='Muscle',
                                            '6'='Monocyte','7'='Endo','8'='TCells','9'='myoFibroblast', '10'='TNNT1_Epis','11'='BCells','12'='Basal','13'='Fibroblast' )

##
MiceAll@meta.data$CC2     <- factor(MiceAll@meta.data$CC2, levels=c('Basal','Luminal','Bas.Lum.Mixed','Fibroblast','myoFibroblast','Muscle','Monocyte','DendriticCells','TCells','Endo','Neurone'))
HomoAll_HMY@meta.data$CC2 <- factor(HomoAll_HMY@meta.data$CC2, levels=c('Basal','TNNT1_Epis','Intermediate','Luminal','Fibroblast','myoFibroblast','Muscle','BCells','TCells','Monocyte','Endo'))
VlnPlot(MiceAll,     group.by='CC2', features='Pin1')
VlnPlot(HomoAll_HMY, group.by='CC2', features='PIN1')

## 
MiceAllAlra <- RunALRA(MiceAll)
HomoAllAlra <- RunALRA(HomoAll_HMY)

## 
library(ggsci); library(ggplot2)
CreatGGplotFromSeurat <- function(xx, assay=NA, features=NA, meta=NA){
   require(dplyr); DefaultAssay(xx) <- assay
   features <- features[features %in% rownames(xx)]
   if (length(features)<1){ stop('None of the features include in the assay') }
   print(features)
   tmpAssay <- GetAssay(xx[features, ], assay=assay)@data
   tmpAssay <- tmpAssay %>% as.data.frame() %>% t()
   if (is.na(meta)){meta = colnames(xx@meta.data)}
   x <- cbind(xx@meta.data[,meta], tmpAssay)
   colnames(x) <- gsub('-', '_', colnames(x))
   x
}

WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(MiceAllAlra, assay='RNA', features=c('Pin1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','Pin1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig03
ggsave(Fig03, file='./Mice.PIN1.RawExpr.pdf')

WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(MiceAllAlra, assay='alra', features=c('Pin1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','Pin1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig04
ggsave(Fig04, file='./Mice.PIN1.Alra.pdf')


WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(HomoAllAlra, assay='RNA', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','PIN1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig05
ggsave(Fig05, file='./Homo.PIN1.RawExpr.pdf')

WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(HomoAllAlra, assay='alra', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC2','PIN1')], id=c('CC2'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC2, value, fill=CC2)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig06
ggsave(Fig06, file='./Homo.PIN1.Alra.pdf')


##
WGCNA::sizeGrWindow(6,6)
DimPlot(MiceAllAlra, reduction='tsne', label=T, group.by='CC2') + NoLegend()  -> Fig11
DimPlot(HomoAllAlra, reduction='tsne', label=T, group.by='CC2') + NoLegend()  -> Fig21
#
DefaultAssay(MiceAllAlra) <- 'RNA'  ; FeaturePlot(MiceAllAlra, reduction='tsne', order=T, pt.size=1, features='Pin1') + scale_color_viridis_c() -> Fig12
DefaultAssay(MiceAllAlra) <- 'alra' ; FeaturePlot(MiceAllAlra, reduction='tsne', order=T, pt.size=1, features='Pin1') + scale_color_viridis_c() -> Fig13

DefaultAssay(HomoAllAlra) <- 'RNA'  ; FeaturePlot(HomoAllAlra, reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig22
DefaultAssay(HomoAllAlra) <- 'alra' ; FeaturePlot(HomoAllAlra, reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig23

#FeaturePlot(MiceAllAlra, features='Pin1') + scale_color_viridis_c()  -> Fig12
#FeaturePlot(HomoAllAlra, features='PIN1') + scale_color_viridis_c()  -> Fig22
#
#Nebulosa::plot_density(MiceAllAlra, reduction='tsne', "Pin1")  -> Fig13
#Nebulosa::plot_density(HomoAllAlra, reduction='tsne', "PIN1")  -> Fig23

ggsave(Fig11, file='./Mice.PIN1.TSNE.Clusters.pdf')
ggsave(Fig12, file='./Mice.PIN1.TSNE.RNA.pdf')
ggsave(Fig13, file='./Mice.PIN1.TSNE.Alra.pdf')
ggsave(Fig21, file='./Homo.PIN1.TSNE.Clusters.pdf')
ggsave(Fig22, file='./Homo.PIN1.TSNE.RNA.pdf')
ggsave(Fig23, file='./Homo.PIN1.TSNE.Alra.pdf')

