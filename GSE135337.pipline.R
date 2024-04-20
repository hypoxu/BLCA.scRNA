### 
### GSE135337
### 
library(Seurat); library(dplyr)
setwd('./Data.Homo.GSE135337.IJC')
# 
BC1 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4006644_BC1_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC2 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4006645_BC2_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC3 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4006646_BC3_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC4 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4006647_BC4_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC5 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4006648_BC5_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC6 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4751267_BC6_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BC7 <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM4751268_BC7_gene_cell_exprs_table.txt.gz") %>% as.data.frame()
BCN <- data.table::fread("./Data.Homo.GSE135337.IJC/Data.GEO/GSM5329919_BCN_gene_cell_exprs_table.xls.gz") %>% as.data.frame()
#
BC1 <- aggregate(BC1[ ,3:ncol(BC1)], by=list(BC1$Symbol), sum)
BC2 <- aggregate(BC2[ ,3:ncol(BC2)], by=list(BC2$Symbol), sum)
BC3 <- aggregate(BC3[ ,3:ncol(BC3)], by=list(BC3$Symbol), sum)
BC4 <- aggregate(BC4[ ,3:ncol(BC4)], by=list(BC4$Symbol), sum)
BC5 <- aggregate(BC5[ ,3:ncol(BC5)], by=list(BC5$Symbol), sum)
BC6 <- aggregate(BC6[ ,3:ncol(BC6)], by=list(BC6$Symbol), sum)
BC7 <- aggregate(BC7[ ,3:ncol(BC7)], by=list(BC7$Symbol), sum)
BCN <- aggregate(BCN[ ,3:ncol(BCN)], by=list(BCN$Symbol), sum)
#
rownames(BC1) <- BC1$Group.1; BC1 <- BC1[ ,-1]
rownames(BC2) <- BC2$Group.1; BC2 <- BC2[ ,-1]
rownames(BC3) <- BC3$Group.1; BC3 <- BC3[ ,-1]
rownames(BC4) <- BC4$Group.1; BC4 <- BC4[ ,-1]
rownames(BC5) <- BC5$Group.1; BC5 <- BC5[ ,-1]
rownames(BC6) <- BC6$Group.1; BC6 <- BC6[ ,-1]
rownames(BC7) <- BC7$Group.1; BC7 <- BC7[ ,-1]
rownames(BCN) <- BCN$Group.1; BCN <- BCN[ ,-1]

###
BC1st <- CreateSeuratObject(counts = BC1, project = "BC1", min.cells = 1, min.features = 200)
BC2st <- CreateSeuratObject(counts = BC2, project = "BC2", min.cells = 1, min.features = 200)
BC3st <- CreateSeuratObject(counts = BC3, project = "BC3", min.cells = 1, min.features = 200)
BC4st <- CreateSeuratObject(counts = BC4, project = "BC4", min.cells = 1, min.features = 200)
BC5st <- CreateSeuratObject(counts = BC5, project = "BC5", min.cells = 1, min.features = 200)
BC6st <- CreateSeuratObject(counts = BC6, project = "BC6", min.cells = 1, min.features = 200)
BC7st <- CreateSeuratObject(counts = BC7, project = "BC7", min.cells = 1, min.features = 200)
BCNst <- CreateSeuratObject(counts = BCN, project = "BCN", min.cells = 1, min.features = 200)

#
BC1st[["percent.mt"]] <- PercentageFeatureSet(BC1st, pattern = "^MT-")
BC2st[["percent.mt"]] <- PercentageFeatureSet(BC2st, pattern = "^MT-")
BC3st[["percent.mt"]] <- PercentageFeatureSet(BC3st, pattern = "^MT-")
BC4st[["percent.mt"]] <- PercentageFeatureSet(BC4st, pattern = "^MT-")
BC5st[["percent.mt"]] <- PercentageFeatureSet(BC5st, pattern = "^MT-")
BC6st[["percent.mt"]] <- PercentageFeatureSet(BC6st, pattern = "^MT-")
BC7st[["percent.mt"]] <- PercentageFeatureSet(BC7st, pattern = "^MT-")
BCNst[["percent.mt"]] <- PercentageFeatureSet(BCNst, pattern = "^MT-")
#
AddMetaData <- function(xx, Patient=NA, Gender=NA, Age=NA, Stage=NA, Types=NA, Source='GuangXi', addID){
  xx[['Source' ]] <- Source  #  
  xx[['Patient']] <- Patient #  
  xx[['Gender' ]] <- Gender  #  
  xx[['Age'    ]] <- Age     #  
  xx[['Stage'  ]] <- Stage   #  
  xx[['Types'  ]] <- Types   #  
  xx <- RenameCells(xx, add.cell.id = addID) # 
  xx
}
BC1st <- AddMetaData(BC1st, Patient='BC1', addID='BC1', Source='SH', Gender='M',  Age=NA, Stage='pTa', Types='T')
BC2st <- AddMetaData(BC2st, Patient='BC2', addID='BC2', Source='SH', Gender='M',  Age=NA, Stage='pT1', Types='T')
BC3st <- AddMetaData(BC3st, Patient='BC3', addID='BC3', Source='SH', Gender='M',  Age=NA, Stage='pT1', Types='T')
BC4st <- AddMetaData(BC4st, Patient='BC4', addID='BC4', Source='SH', Gender='M',  Age=NA, Stage='pT2', Types='T')
BC5st <- AddMetaData(BC5st, Patient='BC5', addID='BC5', Source='SH', Gender='M',  Age=NA, Stage='pT3', Types='T')
BC6st <- AddMetaData(BC6st, Patient='BC6', addID='BC6', Source='SH', Gender='M',  Age=NA, Stage='pTa', Types='T')
BC7st <- AddMetaData(BC7st, Patient='BC7', addID='BC7', Source='SH', Gender='M',  Age=NA, Stage='pTa', Types='T')
BCNst <- AddMetaData(BCNst, Patient='BCN', addID='BCN', Source='SH', Gender='M',  Age=NA, Stage='pT2', Types='N')

#
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
BC1st <- SimpleRunSeuratPipline(BC1st, ndims=1:25)
BC2st <- SimpleRunSeuratPipline(BC2st, ndims=1:25)
BC3st <- SimpleRunSeuratPipline(BC3st, ndims=1:25)
BC4st <- SimpleRunSeuratPipline(BC4st, ndims=1:25)
BC5st <- SimpleRunSeuratPipline(BC5st, ndims=1:25)
BC6st <- SimpleRunSeuratPipline(BC6st, ndims=1:25)
BC7st <- SimpleRunSeuratPipline(BC7st, ndims=1:25)
BCNst <- SimpleRunSeuratPipline(BCNst, ndims=1:25)
# 
DimPlot(BC1st, label=T); VlnPlot(BC1st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC2st, label=T); VlnPlot(BC2st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC3st, label=T); VlnPlot(BC3st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC4st, label=T); VlnPlot(BC4st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC5st, label=T); VlnPlot(BC5st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC6st, label=T); VlnPlot(BC6st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BC7st, label=T); VlnPlot(BC7st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
DimPlot(BCNst, label=T); VlnPlot(BCNst, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

# save(BC1st,BC2st,BC3st,BC4st,BC5st,BC6st,BC7st,BCNst, file='IJC.S001.ProccSeurat.RData')
MarkerJASN <- data.frame(Genes=c('EPCAM',"KRT18","KRT19","KRT5","KRT17","KRT13","KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2",
                                 "VIM","S100A4","COL3A1","COL1A1","COL1A2","ACTA2","TAGLN","DES","CNN1","ACTG2","TPM2",
								 "SELE","PECAM1","VCAM1","CDH5","LYZ","MS4A7","CD14","CD209","CD3D","CD3E","MZB1","CD79A","GPM6A"),
                       Markers=c('Epi','Epi','Epi','Bas','Bas','BasInter','Umbr','UmbrInter','UmbrInter','UmbrInter','Umbr','Umbr',
					             'Interstitial','Fibr','Fibr','Fibr','Fibr','myoFibr','Muscle','Muscle','Muscle','Muscle','Muscle',
								 'Endo','Endo','Endo','Endo','Mono','Mono','Mono','Dendritic','Tcell','Tcell','Bcell','Bcell','Neuron'))
#Reference component analysis of single-cell transcriptomes elucidates cellular heterogeneity in human colorectal tumors
NG_LHP <- data.frame(Genes=c('VIL1', 'KRT20', 'CLDN7', 'CDH1','SPARC', 'COL14A', 'COL3A1', 'DCN','CD38', 'MZB1', 'DERL3','TRBC2', 'CD3D', 'CD3E', 'CD3G',
                             'ITGAX', 'CD68', 'CD14', 'CCL3','KIT', 'TPSB2','PTPRC','PTPRB','PECAM1','EPCAM'), 
                    Labels=c(rep('Epi',4),rep('Fibr',4),rep('BCells',3),rep('TCells',4),rep('myeloid',4),rep('Mast',2),'Immune','Endo','Endo','Epi'))
DoHeatmap(BC1st, features=NG_LHP$Genes); dev.new(); DimPlot(BC1st, label=T)
DoHeatmap(BC2st, features=NG_LHP$Genes); dev.new(); DimPlot(BC2st, label=T)
DoHeatmap(BC3st, features=NG_LHP$Genes); dev.new(); DimPlot(BC3st, label=T)
DoHeatmap(BC4st, features=NG_LHP$Genes); dev.new(); DimPlot(BC4st, label=T)
DoHeatmap(BC5st, features=NG_LHP$Genes); dev.new(); DimPlot(BC5st, label=T)
DoHeatmap(BC6st, features=NG_LHP$Genes); dev.new(); DimPlot(BC6st, label=T)
DoHeatmap(BC7st, features=NG_LHP$Genes); dev.new(); DimPlot(BC7st, label=T)
DoHeatmap(BCNst, features=NG_LHP$Genes); dev.new(); DimPlot(BCNst, label=T)

RCCAll_HMY@meta.data$CC3 <- dplyr::recode(RCCAll_HMY@meta.data$seurat_clusters, '0'='TCells', '1'='TCells','2'='Epi', '3'='myoCAF', '4'='BCells', '5'='Mono',
                                                          '6'='Endo','7'='Epi','8'='Mono')
RCCAll_HMY@meta.data$CC3 <- dplyr::recode(RCCAll_HMY@meta.data$seurat_clusters, '0'='TCells', '1'='TCells','2'='Epi', '3'='myoCAF', '4'='BCells', '5'='Mono',
                                                          '6'='Endo','7'='Epi','8'='Immune','9'='myoCAF', '10'='Mono.Low', '11'='Mast','12'='Epi','13'='Epi',
														  '14'='TCells.Low','15'='Epi','16'='Fibr.Low' )

##
library(harmony); library(SeuratWrappers); library(patchwork)
#
BLCA_IJC <- merge(BC1st, list(BC2st,BC3st,BC4st,BC5st,BC6st,BC7st,BCNst))
BLCA_IJC <- SimpleRunSeuratPipline(BLCA_IJC, ndims=1:20)
#
BLCAIJC_HMY <- RunHarmony(BLCA_IJC, group.by.vars = c("orig.ident"), assays='RNA')
BLCAIJC_HMY <- RunUMAP(BLCAIJC_HMY, reduction = "harmony", dims = 1:25)
BLCAIJC_HMY <- RunTSNE(BLCAIJC_HMY, reduction = "harmony", dims = 1:25)
BLCAIJC_HMY <- FindNeighbors(BLCAIJC_HMY, reduction = "harmony", dims = 1:25) %>% FindClusters(reduction = "harmony", resolution = c(1, 0.8, 0.6, 0.4))
#
DimPlot(BLCAIJC_HMY, label=T, group.by='seurat_clusters') 
DoHeatmap(RCCAll_HMY,features=MarkerJASN$Genes); dev.new(); DimPlot(RCCAll_HMY, label=T)
#
BLCAIJC_HMY@meta.data$CC1 <- dplyr::recode(BLCAIJC_HMY@meta.data$seurat_clusters, '0'='Epi', '1'='Epi','2'='Epi', '3'='Epi', '4'='Epi', '5'='CAFs',
                                                          '6'='Epi','7'='Myeloid','8'='TCells','9'='Epi', '10'='BCells',
														  '11'='CAFs','12'='CAFs','13'='Epi','14'='Endo','15'='Epi' )
BLCAIJC_HMY@meta.data$CC2 <- dplyr::recode(BLCAIJC_HMY@meta.data$seurat_clusters, '0'='Epi', '1'='Epi','2'='Epi', '3'='Epi', '4'='Epi', '5'='CAFs',
                                                          '6'='Epi','7'='Myeloid','8'='TCells','9'='Epi', '10'='BCells',
														  '11'='CAFs','12'='CAFs','13'='Epi','14'='Endo','15'='Epi' )
BLCAIJC_HMY@meta.data$CC2 <- as.character(BLCAIJC_HMY@meta.data$CC2)
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC1' ] <- 'BC1.Epi'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC2' ] <- 'BC2.Hetero'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC3' ] <- 'BC3.Hetero'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC4' ] <- 'BC4.Basal'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC5' ] <- 'BC5.Epi'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC6' ] <- 'BC6.Epi'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BC7' ] <- 'BC7.Epi'
BLCAIJC_HMY@meta.data$CC2[ BLCAIJC_HMY@meta.data$CC2=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BCN' ] <- 'Nor.Epi'
#
BLCAIJC_HMY@meta.data$CC2 <- factor(BLCAIJC_HMY@meta.data$CC2, levels=c("Nor.Epi","BC1.Epi","BC2.Hetero","BC3.Hetero","BC4.Basal","BC5.Epi","BC6.Epi","BC7.Epi","Myeloid","BCells","TCells","CAFs","Endo"))
BLCAIJC_HMY@meta.data$CC3 <- as.character(BLCAIJC_HMY@meta.data$CC1)
BLCAIJC_HMY@meta.data$CC3[ BLCAIJC_HMY@meta.data$CC3=='Epi' & BLCAIJC_HMY@meta.data$orig.ident=='BCN' ] <- 'Normal.Epi'
BLCAIJC_HMY@meta.data$CC3[ BLCAIJC_HMY@meta.data$CC3=='Epi' & BLCAIJC_HMY@meta.data$orig.ident!='BCN' ] <- 'Tumor.Epi'
BLCAIJC_HMY@meta.data$CC3 <- factor(BLCAIJC_HMY@meta.data$CC3, levels=c("Normal.Epi","Tumor.Epi","Myeloid","BCells","TCells","CAFs","Endo"))

###
BLCA_BC1Epi <- subset(BLCAIJC_HMY, orig.ident=='BC1' & CC1=='Epi')
BLCA_BC2Epi <- subset(BLCAIJC_HMY, orig.ident=='BC2' & CC1=='Epi')
BLCA_BC3Epi <- subset(BLCAIJC_HMY, orig.ident=='BC3' & CC1=='Epi')
BLCA_BC4Epi <- subset(BLCAIJC_HMY, orig.ident=='BC4' & CC1=='Epi')
BLCA_BC5Epi <- subset(BLCAIJC_HMY, orig.ident=='BC5' & CC1=='Epi')
BLCA_BC6Epi <- subset(BLCAIJC_HMY, orig.ident=='BC6' & CC1=='Epi')
BLCA_BC7Epi <- subset(BLCAIJC_HMY, orig.ident=='BC7' & CC1=='Epi')
BLCA_BCNEpi <- subset(BLCAIJC_HMY, orig.ident=='BCN' & CC1=='Epi')

SimpleRunSeuratPipline2 <- function(xx, nfeatures=3000, ndims=1:50){
   require(Seurat); require(dplyr)
   xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = nfeatures)
   xx <- RunPCA(xx, features = VariableFeatures(object = xx), npcs=100)
   xx <- FindNeighbors(xx, dims = ndims) %>% FindClusters(resolution = c(1, 0.8, 0.6)) %>% RunUMAP(dims = ndims) %>% RunTSNE(dims = ndims)
   xx
}
BLCA_BC1Epi <- SimpleRunSeuratPipline2(BLCA_BC1Epi, ndims=1:25)
BLCA_BC2Epi <- SimpleRunSeuratPipline2(BLCA_BC2Epi, ndims=1:25)
BLCA_BC3Epi <- SimpleRunSeuratPipline2(BLCA_BC3Epi, ndims=1:25)
BLCA_BC4Epi <- SimpleRunSeuratPipline2(BLCA_BC4Epi, ndims=1:25)
BLCA_BC5Epi <- SimpleRunSeuratPipline2(BLCA_BC5Epi, ndims=1:25)
BLCA_BC6Epi <- SimpleRunSeuratPipline2(BLCA_BC6Epi, ndims=1:25)
BLCA_BC7Epi <- SimpleRunSeuratPipline2(BLCA_BC7Epi, ndims=1:25)
BLCA_BCNEpi <- SimpleRunSeuratPipline2(BLCA_BCNEpi, ndims=1:25)

#
IJC_Marker <- data.frame(Genes=c('KRT6A','KRT5','KRT14','KRT16','CDH3','CD44',
                                 'XBP1','UPK3A','UPK3B','UPK2','UPK1A','UPK1B','PPARG','KRT20','GATA3','FOXA1','FGFR3','ERBB3','ERBB2','CYP2J2','SNX31',
                                 'SMC4','MKI67','CDK1','TYMS','CDCA7'),
                        Labels=c(rep('Basal', 6), rep('Luminal', 15), rep('Cycle', 5)))
DoHeatmap(tmps, features=IJC_Marker$Genes);  
DoHeatmap(BLCAIJC_HMY, features=IJC_Marker$Genes);  

DimPlot(BLCA_BC1Epi, label=T) + DoHeatmap(BLCA_BC1Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BC2Epi, label=T) + DoHeatmap(BLCA_BC2Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BC3Epi, label=T) + DoHeatmap(BLCA_BC3Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BC4Epi, label=T) + DoHeatmap(BLCA_BC4Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BC5Epi, label=T) + DoHeatmap(BLCA_BC5Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BC6Epi, label=T) + DoHeatmap(BLCA_BC6Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BC7Epi, label=T) + DoHeatmap(BLCA_BC7Epi, features=IJC_Marker$Genes);  
DimPlot(BLCA_BCNEpi, label=T) + DoHeatmap(BLCA_BCNEpi, features=IJC_Marker$Genes);  
#

MarkerJASN <- data.frame(Genes=c('EPCAM',"KRT18","KRT19","KRT5","KRT17","KRT13","KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2",
                                 "VIM","S100A4","COL3A1","COL1A1","COL1A2","ACTA2","TAGLN","DES","CNN1","ACTG2","TPM2",
								 "SELE","PECAM1","VCAM1","CDH5","LYZ","MS4A7","CD14","CD209","CD3D","CD3E","MZB1","CD79A","GPM6A"),
                       Markers=c('Epi','Epi','Epi','Bas','Bas','BasInter','Umbr','UmbrInter','UmbrInter','UmbrInter','Umbr','Umbr',
					             'Interstitial','Fibr','Fibr','Fibr','Fibr','myoFibr','Muscle','Muscle','Muscle','Muscle','Muscle',
								 'Endo','Endo','Endo','Endo','Mono','Mono','Mono','Dendritic','Tcell','Tcell','Bcell','Bcell','Neuron'))
#
DimPlot(BLCA_BC1Epi, label=T) + DoHeatmap(BLCA_BC1Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BC2Epi, label=T) + DoHeatmap(BLCA_BC2Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BC3Epi, label=T) + DoHeatmap(BLCA_BC3Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BC4Epi, label=T) + DoHeatmap(BLCA_BC4Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BC5Epi, label=T) + DoHeatmap(BLCA_BC5Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BC6Epi, label=T) + DoHeatmap(BLCA_BC6Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BC7Epi, label=T) + DoHeatmap(BLCA_BC7Epi, features=MarkerJASN$Genes);  
DimPlot(BLCA_BCNEpi, label=T) + DoHeatmap(BLCA_BCNEpi, features=MarkerJASN$Genes);  
#
BLCA_BC1Epi@meta.data$CC4 <- dplyr::recode(BLCA_BC1Epi@meta.data$seurat_clusters, '0'='Luminal', '1'='Luminal','2'='Luminal', '3'='Luminal', '4'='Luminal', '5'='Luminal.Cycle', '6'='Luminal')
BLCA_BC2Epi@meta.data$CC4 <- dplyr::recode(BLCA_BC2Epi@meta.data$seurat_clusters, '0'='Luminal', '1'='Luminal','2'='Luminal.Cycle', '3'='Luminal.Cycle', '4'='Luminal', '5'='Basal', '6'='Luminal','7'='Luminal','8'='Luminal' )
BLCA_BC3Epi@meta.data$CC4 <- dplyr::recode(BLCA_BC3Epi@meta.data$seurat_clusters, '0'='Luminal', '1'='Luminal','2'='Luminal', '3'='Luminal', '4'='Luminal', '5'='Luminal.Cycle', '6'='Luminal','7'='Luminal','8'='Basal','9'='Luminal' )
BLCA_BC4Epi@meta.data$CC4 <- 'Basal'
BLCA_BC5Epi@meta.data$CC4 <- dplyr::recode(BLCA_BC5Epi@meta.data$seurat_clusters, '0'='Luminal', '1'='Luminal','2'='Luminal.Cycle', '3'='Luminal', '4'='Luminal', '5'='Luminal', '6'='Luminal','7'='Luminal','8'='Luminal' )
BLCA_BC6Epi@meta.data$CC4 <- dplyr::recode(BLCA_BC6Epi@meta.data$seurat_clusters, '0'='Basal', '1'='Basal','2'='Basal', '3'='Basal', '4'='Basal', '5'='Basal', '6'='Basal','7'='Basal', '8'='Basal','9'='Basal',
                                                                                  '10'='Basal', '11'='Basal','12'='Basal', '13'='Basal', '14'='Basal', '15'='Basal.Cycle')
BLCA_BC7Epi@meta.data$CC4 <- 'Luminal'
BLCA_BCNEpi@meta.data$CC4 <- dplyr::recode(BLCA_BCNEpi@meta.data$seurat_clusters, '0'='Luminal', '1'='Luminal','2'='Luminal', '3'='Luminal', '4'='Luminal', '5'='Luminal', '6'='Luminal','7'='Luminal',
                                           '8'='Basal', '9'='Cycle' )
tmp <- rbind(BLCA_BC1Epi@meta.data,BLCA_BC2Epi@meta.data,BLCA_BC3Epi@meta.data,BLCA_BC4Epi@meta.data,
             BLCA_BC5Epi@meta.data,BLCA_BC6Epi@meta.data,BLCA_BC7Epi@meta.data,BLCA_BCNEpi@meta.data)

BLCAIJC_HMY@meta.data$CC4 <- as.character( tmp[rownames(BLCAIJC_HMY@meta.data), 'CC4'] )
BLCAIJC_HMY@meta.data$CC4[is.na(BLCAIJC_HMY@meta.data$CC4)] <- as.character(BLCAIJC_HMY@meta.data$CC3)[is.na(BLCAIJC_HMY@meta.data$CC4)]

DimPlot(BLCAIJC_HMY, label=T) + DoHeatmap(BLCAIJC_HMY, features=MarkerJASN$Genes);  
DimPlot(BLCAIJC_HMY, label=T) + DoHeatmap(BLCAIJC_HMY, features=MarkerJASN$Genes);  

DimPlot(BLCAIJC_HMY, label=T, group.by='CC3') 
VlnPlot(BLCAIJC_HMY, features = c("PIN1"), group.by='CC4') 
VlnPlot(BLCAIJC_HMY, features = c("PIN1"), group.by='CC4',split.by='orig.ident') 
table(BLCAIJC_HMY@meta.data$orig.ident, BLCAIJC_HMY@meta.data$CC4)

#
library(DoMultiBarHeatmap)
DoMultiBarHeatmap(BLCAIJC_HMY, features=MarkerJASN$Genes, angle = 45, group.by='CC4', additional.group.by = 'orig.ident', slot='scale.data', assay='RNA') 

BLCAIJC_HMY@meta.data$CC4[which(BLCAIJC_HMY@meta.data$orig.ident=='BCN' & BLCAIJC_HMY@meta.data$CC4 == 'Luminal')] <- 'Normal.Lum'
BLCAIJC_HMY@meta.data$CC4[which(BLCAIJC_HMY@meta.data$orig.ident=='BCN' & BLCAIJC_HMY@meta.data$CC4 == 'Basal')]   <- 'Normal.Basal'
BLCAIJC_HMY@meta.data$CC4[which(BLCAIJC_HMY@meta.data$orig.ident=='BCN' & BLCAIJC_HMY@meta.data$CC4 == 'Cycle')]   <- 'Normal.Cycle'

DoMultiBarHeatmap(BLCAIJC_HMY, features=MarkerJASN$Genes, angle = 45, group.by='CC4', additional.group.by = 'orig.ident', slot='scale.data', assay='RNA') 

BLCAIJC_HMY@meta.data$CC4 <- factor(BLCAIJC_HMY@meta.data$CC4, levels=c("Normal.Basal","Normal.Lum","Normal.Cycle", "Basal","Basal.Cycle", "Luminal","Luminal.Cycle","BCells","Myeloid","TCells","CAFs","Endo"))
VlnPlot(BLCAIJC_HMY, features = c("PIN1"), group.by='CC4', pt.size=0.3) 



### imputation   
library(Seurat); library(dplyr); library(SeuratWrappers)
setwd('E:/Data.Analysis.SSPH/Data.BLCA.Task/Data.Homo.GSE135337.IJC')
load('Data.Analysis.Step04.HMY.20020719.RData')
BLCAIJC_Alra <- RunALRA(BLCAIJC_HMY)
#
## Plotting ##
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
##
WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(BLCAIJC_HMY, assay='RNA', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC4','PIN1')], id=c('CC4'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC4, value, fill=CC4)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + #facet_wrap(~variable, scales='free') + 
   geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  -> Fig01
ggsave(Fig01, file='./PIN1.RawExpr.pdf')
#
WGCNA::sizeGrWindow(7,5)
tmpMarker1 <- CreatGGplotFromSeurat(BLCAIJC_Alra, assay='alra', features=c('PIN1'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC4','PIN1')], id=c('CC4'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC4, value, fill=CC4)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2) -> Fig02  
ggsave(Fig02, file='./PIN1.Alra.pdf')
#
tmpMarker1 <- CreatGGplotFromSeurat(BLCAIJC_Alra, assay='alra', features=c('EPCAM'))
tmpMarker2 <- reshape::melt(tmpMarker1[,c('CC4','EPCAM')], id=c('CC4'))
LymCol2    <- c( pal_npg( "nrc")(9), '#DFC297', '#FEB24C', '#ACD851')  
ggplot(tmpMarker2, aes(CC4, value, fill=CC4)) + 
   geom_violin(position=position_dodge(0.8), scale='width') + geom_jitter(width=0.2, size=0.2) + 
   geom_boxplot(outlier.alpha=0, width=0.2, position=position_dodge(0.8), color='white') +
   theme_classic() + theme(legend.position = "none", validate = TRUE) + 
   theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color='black')) +
   xlab('') + ylab('Expression Level (logUMI)') + scale_color_manual(values=LymCol2) + scale_fill_manual(values=LymCol2)  


## Plotting ##
WGCNA::sizeGrWindow(6,4)
#DimPlot(BLCAIJC_Alra, reduction='umap', label=T, group.by='CC4') #+ NoLegend() # -> Fig11
DimPlot(subset(BLCAIJC_Alra, Types=='N'), reduction='tsne', label=T, group.by='CC4' ) -> Fig11
DimPlot(subset(BLCAIJC_Alra, Types=='T'), reduction='tsne', label=T, group.by='CC4' ) -> Fig12
#
WGCNA::sizeGrWindow(6,6)
DefaultAssay(BLCAIJC_Alra) <- 'RNA'  ; FeaturePlot(subset(BLCAIJC_Alra, Types=='N'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig13
DefaultAssay(BLCAIJC_Alra) <- 'RNA'  ; FeaturePlot(subset(BLCAIJC_Alra, Types=='T'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig14
DefaultAssay(BLCAIJC_Alra) <- 'alra' ; FeaturePlot(subset(BLCAIJC_Alra, Types=='N'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig15
DefaultAssay(BLCAIJC_Alra) <- 'alra' ; FeaturePlot(subset(BLCAIJC_Alra, Types=='T'), reduction='tsne', order=T, pt.size=1, features='PIN1') + scale_color_viridis_c() -> Fig16
#
WGCNA::sizeGrWindow(6,4); ggsave(Fig11, file='./Homo.PIN1.TSNE.Clusters.Nor.pdf')
WGCNA::sizeGrWindow(6,4); ggsave(Fig12, file='./Homo.PIN1.TSNE.Clusters.Tum.pdf')
WGCNA::sizeGrWindow(6,6); ggsave(Fig13, file='./Homo.PIN1.TSNE.RNA.Nor.pdf')
WGCNA::sizeGrWindow(6,6); ggsave(Fig14, file='./Homo.PIN1.TSNE.RNA.Tum.pdf')
WGCNA::sizeGrWindow(6,6); ggsave(Fig15, file='./Homo.PIN1.TSNE.Alra.Nor.pdf')
WGCNA::sizeGrWindow(6,6); ggsave(Fig16, file='./Homo.PIN1.TSNE.Alra.Tum.pdf')

















#Reference component analysis of single-cell transcriptomes elucidates cellular heterogeneity in human colorectal tumors
NG_LHP <- data.frame(Genes=c('VIL1', 'KRT20', 'CLDN7', 'CDH1','SPARC', 'COL14A', 'COL3A1', 'DCN','CD38', 'MZB1', 'DERL3','TRBC2', 'CD3D', 'CD3E', 'CD3G','ITGAX', 'CD68', 'CD14', 'CCL3','KIT', 'TPSB2','PTPRC','PTPRB','PECAM1','EPCAM'), 
                    Labels=c(rep('Epi',4),rep('Fibr',4),rep('BCells',3),rep('TCells',4),rep('myeloid',4),rep('Mast',2),'Immune','Endo','Endo','Epi'))
DoHeatmap(BLCAIJC_HMY, features=NG_LHP$Genes); dev.new(); DimPlot(BLCAIJC_HMY, label=T)

#save(RCC_All, RCCAll_HMY, RCCAll_ITG, RCCAll_MNN, file=file.path('Data.Analysis.Step02.Merged.20020719.RData'))


RandomSampleByCluster <- function(Cells, Source, MaxCells=1000){
   cellList <- lapply(split(Cells, Source), function(x){if(length(x)<MaxCells){ x }else{ sample(x, MaxCells) }})
   cells <- unlist(cellList);  cells
}
tmps <- BLCAIJC_HMY[ ,RandomSampleByCluster(rownames(BLCAIJC_HMY@meta.data), BLCAIJC_HMY$seurat_clusters, MaxCells=100)]
tmps <- ScaleData(tmps, features=unlist(MyeloidMarker_Kevin), assay='RNA')
#DoHeatmap(tmps, group.by='ClusterS2', features=unlist(MyeloidMarker_Kevin), assay='alra') 
DoHeatmap(tmps, features=NG_LHP$Genes); dev.new(); DimPlot(tmps, label=T)

   

##
#GeneAlias     <- data.frame(row.names=rownames(BLADDER), Symbol=rownames(BLADDER), alias=rownames(BLADDER))
#GeneAlias$New <- with(GeneAlias, limma::alias2SymbolTable(Symbol))
#
BLCA_IJC <- merge(BC1st, list(BC2st,BC3st,BC4st,BC5st,BC6st,BC7st,BCNst))


VlnPlot(BLCA_IJC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
#
plot1 <- FeatureScatter(Primary, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Primary, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
BLCA_IJC <- NormalizeData(BLCA_IJC, normalization.method = "LogNormalize", scale.factor = 10000)
BLCA_IJC <- FindVariableFeatures(BLCA_IJC, selection.method = "vst", nfeatures = 2500)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BLCA_IJC), 30)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BLCA_IJC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

##
BLCA_IJC <- ScaleData(BLCA_IJC, features = rownames(BLCA_IJC))
BLCA_IJC <- RunPCA(BLCA_IJC, features = VariableFeatures(object = BLCA_IJC))
DimPlot(BLCA_IJC, reduction = "pca", group.by='Treat')
DimPlot(BLCA_IJC, reduction = "pca", group.by='cell_type')

DimHeatmap(BLCA_IJC, dims = 1:24, cells = 500, balanced = TRUE)
#BLCA_IJC <- JackStraw(BLCA_IJC, dims=50, num.replicate = 100) %>% ScoreJackStraw(dims = 1:50)
#JackStrawPlot(BLCA_IJC, dims = 1:50)
ElbowPlot(BLCA_IJC, ndims = 50) 

BLCA_IJC <- FindNeighbors(BLCA_IJC, dims = 1:20) %>% FindClusters(resolution = 0.5)
BLCA_IJC <- RunUMAP(BLCA_IJC, dims = 1:20) %>% RunTSNE(dims = 1:20)
##
DimPlot(BLCA_IJC, reduction = "tsne", label=T) + FeaturePlot(BLCA_IJC, features = c("COL14A1", "ACTA2", "CD44", "ALDH1A3","EPCAM",'KRT18','KRT14'))

DimPlot(BLCA_IJC, reduction = "umap", group.by='orig.ident', label=T)   
DimPlot(BLCA_IJC, reduction = "umap", group.by='Types', label=T) 



save(BLCA_IJC, file='./Data.scRNA.AllIn.01.S01.RData')



library(Seurat); library(SeuratWrappers); library(SeuratDisk)
library(dplyr); library(ggplot2)
#setwd('E:/Data.Analysis.SSPH/Data.SSPH.Kidney.RAS/Data.scRNAseq/')
#load('./Data.Analysis.2022/Data.scRNA.AllIn.01.RData')
#
#{## 常规Seurat分析
#   SimpleRunSeuratPipline <- function(xx, nfeatures=3000, ndims=1:50){
#      xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000) 
#      xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = nfeatures)
#   #   xx <- ScaleData(xx, features = rownames(xx))
#      xx <- ScaleData(xx, features = VariableFeatures(xx))
#      xx <- RunPCA(xx, features = VariableFeatures(object = xx), npcs=100)
#      xx <- FindNeighbors(xx, dims = ndims) %>% FindClusters(resolution = c(1, 0.8, 0.6)) %>% RunUMAP(dims = ndims) %>% RunTSNE(dims = ndims)
#      gc(); 
#      xx
#   }
#   ccRCC <- SimpleRunSeuratPipline(ccRCC_CN, ndims=1:25, nfeatures=3500)
#   DimPlot(ccRCC, label=T)
#   DimPlot(ccRCC, group.by='orig.ident', label=T)
#   DimPlot(ccRCC, group.by='Types', label=T)
#   DimPlot(ccRCC, group.by='ToSciTum2a', label=T)
#   DimPlot(ccRCC, group.by='ToSciNor2a', label=T)
#   save(ccRCC, file='./Data.Analysis.2022/Data.scRNA.AllIn.01.S01.RData')
#   FeaturePlot(ccRCC, order=T, features=c('CA9'))
#   FeaturePlot(ccRCC, order=T, features=c('CA9','PTPRC','KDR','ACTA2'))
#}
{## HMY整合数据 BLCAITJ
   library(harmony); library(SeuratWrappers); library(patchwork)
   BLCA_HMY <- RunHarmony(BLCA_IJC, group.by.vars = c("orig.ident"), assays='RNA')
   #ElbowPlot(BLCA_HMY, ndims = 100, reduction='harmony')
   BLCA_HMY <- RunUMAP(BLCA_HMY, reduction = "harmony", dims = 1:25)
   BLCA_HMY <- RunTSNE(BLCA_HMY, reduction = "harmony", dims = 1:25)
   BLCA_HMY <- FindNeighbors(BLCA_HMY, reduction = "harmony", dims = 1:25) %>% FindClusters(reduction = "harmony", resolution = c(1, 0.8, 0.6, 0.4))
   
   DimPlot(BLCA_HMY,  label=T)
   DimPlot(BLCA_HMY,  group.by='orig.ident', label=T)
   DimPlot(BLCA_HMY,  group.by='Types', label=T)
   DimPlot(BLCA_HMY,  group.by='ToSciTum2a', label=T)
   DimPlot(BLCA_HMY,  group.by='ToSciNor2a', label=T)
   DimPlot(BLCA_HMY, split.by='orig.ident', group.by='ToSciTum2a', label=T)
   
   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('CA9','EPCAM','PECAM1','PTPRC'))
   VlnPlot(BLCA_HMY, group.by='ToSciTum2a', features=c('CA9','EPCAM','PECAM1','PTPRC'))
   
   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('CA9','EPCAM','PECAM1','PTPRC')) +  scale_color_viridis_c()
   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('nCount_RNA','nFeature_RNA','percent.mt')) +  scale_color_viridis_c()
   VlnPlot(BLCA_HMY, group.by='ToSciTum2a', features=c('CA9','EPCAM','PECAM1','PTPRC'))

   BLCA_HMY$MainCC <- dplyr::recode(BLCA_HMY$seurat_clusters, '9'='Lymp', '8'='Myel', '5'='Fibr', '12'='Fibr', '13'='Endo', 
                      '0'='Epis', '1'='Epis', '2'='Epis', '3'='Epis', '4'='Epis', '6'='Epis', '7'='Epis', '10'='Epis', '11'='Epis')
   DimPlot(BLCA_HMY, group.by='MainCC', label=T)
   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('EPCAM','KRT5','UPK2','ZFAS1'))  +  scale_color_viridis_c() 

   save(BLCA_HMY, file='./Data.scRNA.AllIn.01.HMY.RData')
}

{## 内置IntegrationAnchors整合数据
   BLCA_IJC.list <- SplitObject(BLCA_IJC, split.by = "orig.ident")
   BLCA_IJC.list <- lapply(X = BLCA_IJC.list, FUN = function(x) {
       x <- NormalizeData(x, verbose = FALSE)
       x <- FindVariableFeatures(x, verbose = FALSE)
   })
   features <- SelectIntegrationFeatures(object.list = BLCA_IJC.list)
   BLCA_IJC.list <- lapply(X = BLCA_IJC.list, FUN = function(x) {
       x <- ScaleData(x, features = features, verbose = FALSE)
       x <- RunPCA(x, features = features, verbose = FALSE)
   })
   anchors   <- FindIntegrationAnchors(object.list = BLCA_IJC.list, reduction = "rpca", dims = 1:50) #reference = c(7, 8), 
   BLCA_ITG <- IntegrateData(anchorset = anchors, dims = 1:25)
   BLCA_ITG <- ScaleData(BLCA_ITG, verbose = FALSE)
   BLCA_ITG <- RunPCA(BLCA_ITG, verbose = FALSE)
   BLCA_ITG <- RunUMAP(BLCA_ITG, dims = 1:25)
   BLCA_ITG <- RunTSNE(BLCA_ITG, dims = 1:25)
   rm(BLCA_IJC.list, features, anchors)
   DimPlot(BLCA_ITG, group.by = "orig.ident", label=T)
   DimPlot(BLCA_ITG, group.by = "Phase", label=T)
   DimPlot(BLCA_ITG, group.by = "Types", label=T)
   DimPlot(BLCA_ITG, group.by = "ToSciTum2a", label=T)
   DimPlot(BLCA_ITG, group.by = "ToSciNor2a", label=T)
   DimPlot(BLCA_ITG, reduction='umap', group.by='Drugs', label=T, split.by='Patients', ncol=2)
   FeaturePlot(BLCA_ITG, reduction='umap', order=T, features=c('CA9','EPCAM','PECAM1','PTPRC'))
   save(BLCA_ITG, file='./Data.scRNA.AllIn.01.ITG.RData')

}

{## MNN整合数据
   library(SeuratWrappers)
   BLCA_MNN <- RunFastMNN(object.list = SplitObject(BLCA_IJC, split.by = "orig.ident")) 
   BLCA_MNN <- RunUMAP(BLCA_MNN, reduction = "mnn", dims = 1:20)
   BLCA_MNN <- RunTSNE(BLCA_MNN, reduction = "mnn", dims = 1:20)
   BLCA_MNN <- FindNeighbors(BLCA_MNN, reduction = "mnn", dims = 1:20)
   BLCA_MNN <- FindClusters(BLCA_MNN, resolution = c(1, 0.8, 0.6, 0.4))
   #BLCAQC_MNN@meta.data$CellType1 <- BLCAQC_MNN@meta.data$cluster_name
   
   DimPlot(BLCA_MNN, label=T)
   DimPlot(BLCA_MNN, group.by='RNA_snn_res.0.6', label=T)
   DimPlot(BLCA_MNN, group.by='orig.ident', label=T)
   DimPlot(BLCA_MNN, group.by='Types', label=T)
   DimPlot(BLCA_MNN, group.by='ToSciTum2a', label=T)
   DimPlot(BLCA_MNN, group.by='ToSciNor2a', label=T)
   DimPlot(BLCA_MNN, group.by='ToSciTum2a', label=F, split.by='orig.ident', ncol=6)
   
   FeaturePlot(BLCA_MNN, order=T, features=c('CA9','EPCAM','PECAM1','PTPRC','ACTA2','KDR','MRC1'))
   FeaturePlot(BLCA_MNN, order=T, features=c('MED7','EPCAM','KDR','PTPRC'))
   VlnPlot(BLCA_MNN, group.by='ToSciTum2a', features=c('CA9','EPCAM','PECAM1','PTPRC'))
   DimPlot(BLCA_MNN, split.by='orig.ident', group.by='ToSciTum2a', label=T)
   save(BLCA_MNN, file='./Data.scRNA.AllIn.01.MNN.RData')
}



   FeaturePlot(BLCA_ITG, reduction='umap', order=T, features=c('COL1A1','EPCAM','PECAM1','PTPRC'))
   FeaturePlot(BLCA_MNN, reduction='umap', order=T, features=c('COL1A1','EPCAM','KDR','PTPRC'))
   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('COL1A1','EPCAM','KDR','PTPRC'))

Nebulosa::plot_density(BLCA_ITG, c('KRT15','EPCAM','PECAM1','PTPRC'), reduction='umap')
Nebulosa::plot_density(BLCA_MNN, c('KRT15','EPCAM','PECAM1','PTPRC'), reduction='umap')
Nebulosa::plot_density(BLCA_HMY, c('KRT15','EPCAM','PECAM1','PTPRC'), reduction='umap')






load("E:/Data.Analysis.SSPH/Data.BLCA.Task/Data.Homo.GSE135337.IJC/Data.scRNA.AllIn.01.HMY.RData")
library(Seurat); library(dplyr)
library(harmony); library(SeuratWrappers); library(patchwork)

DEGs   <- presto::wilcoxauc(BLCA_HMY, group_by='seurat_clusters', seurat_assay = 'RNA');  


load('E:/Data.Analysis.SSPH/Data.BLCA.Task/Homo_sapiens.GRCh38.100.chr.RData') # save(refGTF.Ens100, file=


subset(x, feature %in% subset(refGTF.Ens100, gene_biotype=='lncRNA')$gene_name)
          feature group   avgExpr     logFC statistic       auc          pval          padj    pct_in   pct_out
52749 MIR4435-2HG    10 0.6817538 0.5316223  17660792 0.7371209 3.726348e-199 9.816457e-197  61.19929 15.832071
50484    HOTAIRM1    10 0.7583012 0.3602393  15980916 0.6670067  8.118941e-53  3.460441e-51  69.31217 41.179951
47223       CRNDE    10 0.5689222 0.3059994  15758850 0.6577382  3.464190e-61  1.761500e-59  57.67196 26.737031
52323      MALAT1    10 5.8853552 0.2688116  14780707 0.6169128  9.748027e-22  1.581162e-20 100.00000 99.988167
59374      SNHG25    10 0.7369431 0.2683896  15114040 0.6308253  4.686455e-33  1.171078e-31  66.84303 41.281711
60822     TP53TG1    10 1.1418386 0.2391080  14430009 0.6022754  2.503711e-17  3.280026e-16  85.89065 72.429951
51881   LINC01133    10 0.2836318 0.2095847  14967596 0.6247131 8.251953e-107 8.591855e-105  32.80423  7.482961


   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('CA9','EPCAM','PECAM1','PTPRC')) +  scale_color_viridis_c()
   FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('nCount_RNA','nFeature_RNA','percent.mt')) +  scale_color_viridis_c()
   VlnPlot(BLCA_HMY, group.by='ToSciTum2a', features=c('CA9','EPCAM','PECAM1','PTPRC'))

   BLCA_HMY$MainCC <- dplyr::recode(BLCA_HMY$seurat_clusters, '9'='Lymp', '8'='Myel', '5'='Fibr', '12'='Fibr', '13'='Endo', 
                      '0'='Epis', '1'='Epis', '2'='Epis', '3'='Epis', '4'='Epis', '6'='Epis', '7'='Epis', '10'='Epis', '11'='Epis')
   DimPlot(BLCA_HMY, group.by='MainCC', label=T)
   DimPlot(BLCA_HMY, group.by='Types', label=T)
   DimPlot(BLCA_HMY, label=T)


FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('TP53TG1'))  +  scale_color_viridis_c() 


subset(DEGs, feature %in% subset(refGTF.Ens100, gene_biotype=='lncRNA')$gene_name & group==10 & logFC < -0.2 & pval<=0.01)


subset(DEGs, feature %in% subset(refGTF.Ens100, gene_biotype=='lncRNA')$gene_name & auc>0.75 & logFC>0.3 & pval<=0.01)
subset(DEGs, feature %in% subset(refGTF.Ens100, gene_biotype=='lncRNA')$gene_name & logFC < -0.2 & pval<=0.01)

FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('MIR4435-2HG'))  +  scale_color_viridis_c() 
FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('PCAT19'))  +  scale_color_viridis_c() 
FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('MEG3'))  +  scale_color_viridis_c() 
FeaturePlot(BLCA_HMY, reduction='umap', order=T, features=c('EPCAM'))  +  scale_color_viridis_c() 

VlnPlot(BLCA_HMY, group.by='MainCC', features=c('MED7'), split.by='Types')














#
DimPlot(BLCA_IJC, reduction = "tsne", label=T, group.by='Treat') + 
FeaturePlot(BLCA_IJC, reduction = "tsne", features = c("COL14A1", "ACTA2", "CD44", "ALDH1A1","EPCAM",'KRT18','KRT13','KRT5','KRT72','KRT14'))
#markers <- FindMarkers(BLCA_IJC, ident.1 = 12,  min.pct = 0.25)

DimPlot(BLCA_IJC, reduction = "umap", label=T) + 
FeaturePlot(BLCA_IJC, reduction = "umap", order=T, label=T, features = c("PTPRC", "PECAM1", "EPCAM", "HRAS", "KRT13","KRT5","COL14A1", "ACTA2", "CD44", "ALDH1A3",'KRT18','KRT14'))
#
#BLCA_IJC@meta.data$seurat_clusters
BLCA_IJC@meta.data$MainGroups <- 'unlabeled'
BLCA_IJC@meta.data$MainGroups[which(BLCA_IJC@meta.data$seurat_clusters %in% c(1,6,10))]       <- 'EpiTum'
BLCA_IJC@meta.data$MainGroups[which(BLCA_IJC@meta.data$seurat_clusters %in% c(7))]            <- 'EpiNor'
BLCA_IJC@meta.data$MainGroups[which(BLCA_IJC@meta.data$seurat_clusters %in% c(2,8,3,11,12))]  <- 'Immu'        ## 11->Bcell   12-> granulocytes,  TPSB2+
BLCA_IJC@meta.data$MainGroups[which(BLCA_IJC@meta.data$seurat_clusters %in% c(0))]            <- 'Endothelial' ## 'ENG','S1PR1','EMCN'
BLCA_IJC@meta.data$MainGroups[which(BLCA_IJC@meta.data$seurat_clusters %in% c(4,5,9,13))]     <- 'CAFs'        ##
DimPlot(BLCA_IJC, reduction = "tsne", label=T, group.by='cell_type') + 
DimPlot(BLCA_IJC, reduction = "tsne", label=T, group.by='MainGroups')
FeaturePlot(BLCA_IJC, reduction = "tsne", features = c("CD74", "CD72", "ICAM1", "ALDH1A3","EPCAM",'KRT18','KRT14'))
##
CellMarkers <- list()
CellMarkers[[ 'M1' ]]           <- c('CCR7','IRF5','CD80','CXCL9','CXCL10')
CellMarkers[[ 'M2' ]]           <- c('CD163','CXCL2','TGFBI','MS4A6A','F13A1','FOLR2','MRC1','CSF1R')
CellMarkers[['Inflammatory']]   <- c('IL1B','CEBPB','EGR1','PHLDA1','ZFP36','NR4A1','CD14','CCL3','CCL4','TLR2','NLRP3','NFKBIZ','NFKBID','IER3','TNF','CCL2','NFE2L2')
CellMarkers[['LC.LIKE']]        <- c('S100B','CD1A','CD1E','FCER1A','GSN')
#
CellMarkers[['CAFs']]           <- c('FAP','THY1','PDPN','MMP1','MMP2','MMP3','MMP7','MMP9','SFRP1','WNT2','WNT5A','CXCL1','CXCL2','CXCL3','CXCL6','CXCL12','CXCL14','CTGF','PDGFRA','PDGFRL','TWIIST2','TGFB3','COL1A1')
CellMarkers[['MatrixCAFS']]     <- c('COL1A1','COL13A1','ITGA8','TCF21','COL14A1')
CellMarkers[['MyoCAFs']]        <- c('ACTA2','ACTG2','TAGLN','MYLK','MYL9','MYL12A','MYL12B','MYL6','MYL6B','MYH10','MYH9','MYH11')
#
CellMarkers[['CD4.8']]          <- c('CD8A','CD8B','CD4')
CellMarkers[['Regulatory']]     <- c('FOXP3','IL2RA')
CellMarkers[['Cytotoxicity']]   <- c('PRF1','NKG7','GZMA','GZMB','GZMH','GZMK','IFNG')
CellMarkers[['Exhaustion']]     <- c('LAG3','PDCD1','CBLB','CXCL13','CTLA','TIGIT','PRDM1','SIRPG','ICOS','BTLA','HAVCR2')
CellMarkers[['Naive']]          <- c('CCR7','SELL','TCF7','LEF1')
CellMarkers[['MHC.class.I.II']] <- c('HLA-A','HLA-B','HLC-C','HLA-DRA','HLA-DRB1','HLA-DRB5')
##
fgsea::writeGmtPathways(CellMarkers, 'Pathways.CellMarkers.gmt')
##
library(VISION); DefaultAssay(BLCA_IJC) <- 'RNA'
ccRCC_Vis      <- Vision(BLCA_IJC, assay='RNA', dimRed = "tsne", signatures = 'Pathways.CellMarkers.gmt')
ccRCC_Vis      <- calcSignatureScores(ccRCC_Vis)
ccRCC_VisScore <- getSignatureScores(ccRCC_Vis)
BLCA_IJC[['VISION']] <- CreateAssayObject(t(ccRCC_VisScore))
rm(ccRCC_Vis); rm(ccRCC_VisScore); gc()
FeaturePlot(BLCA_IJC, order=T, features=names(CellMarkers))
#ssGSEA.BLCA.Markers <- gsva(as.matrix(BLCA_IJC@assays$RNA@data), CellMarkers, method="ssgsea", mx.diff=F, kcdf="Gaussian") 
#BLCA_IJC@meta.data   <- cbind(BLCA_IJC@meta.data, t(ssGSEA.BLCA.Markers))
#
DimPlot(BLCA_IJC, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC, reduction = "tsne", features = c("M1","M2","LC.LIKE","CAFs","MatrixCAFS","MyoCAFs","CD4.8","Regulatory","Cytotoxicity","Exhaustion","Naive","MHC.class.I.II","Inflammatory" ))

DimPlot(BLCA_IJC, reduction = "tsne", label=T) + FeaturePlot(BLCA_IJC, reduction = "tsne", features = CellMarkers[[ 'Inflammatory' ]])


DoHeatmap(BLCA_IJC, size = 4, features = c(
'KRT14','KRT5','KRT8','KRT18','ERSR1','PRLR',
'CCR7','IRF5','CD80','CXCL9','CXCL10',
'CD163','CXCL2','TGFBI','MS4A6A','F13A1','FOLR2','MRC1','CSF1R',
'IL1B','CEBPB','EGR1','PHLDA1','ZFP36','NR4A1','CD14','CCL3','CCL4','TLR2','NLRP3','NFKBIZ','NFKBID','IER3','TNF','CCL2','NFE2L2',
'S100B','CD1A','CD1E','FCER1A','GSN',
'FAP','THY1','PDPN','MMP1','MMP2','MMP3','MMP7','MMP9','SFRP1','WNT2','WNT5A','CXCL1','CXCL2','CXCL3','CXCL6','CXCL12','CXCL14','CTGF','PDGFRA','PDGFRL','TWIIST2','TGFB3','COL1A1',
'COL1A1','COL13A1','ITGA8','TCF21','COL14A1',
'ACTA2','ACTG2','TAGLN','MYLK','MYL9','MYL12A','MYL12B','MYL6','MYL6B','MYH10','MYH9','MYH11',
'CD8A','CD8B','CD4',
'FOXP3','IL2RA',
'PRF1','NKG7','GZMA','GZMB','GZMH','GZMK','IFNG',
'LAG3','PDCD1','CBLB','CXCL13','CTLA','TIGIT','PRDM1','SIRPG','ICOS','BTLA','HAVCR2',
'CCR7','SELL','TCF7','LEF1',
'HLA-A','HLA-B','HLC-C','HLA-DRA','HLA-DRB1','HLA-DRB5'
))
save(BLCA_IJC, CellMarkers, GeneAlias, file='GenomeMed.S002.Seurat.RData')

library(SeuratWrappers)
BLCA_IJCalra <- RunALRA(BLCA_IJC)
save(BLCA_IJC, BLCA_IJCalra, CellMarkers, GeneAlias, file='GenomeMed.S002.Seurat.RData')

#save.image("H:/Bio.Analysis.DataBase/BLCA_IJC.scRNA.GenomeMed/scRNAseq.RData")

###
### 针对每一个类群进行细致区分
BLCA_IJC_Epi <- subset(BLCA_IJC, seurat_clusters %in% c(2,6,8,7  ))  # EPCAM
BLCA_IJC_End <- subset(BLCA_IJC, seurat_clusters %in% c(0        ))  # PECAM1
BLCA_IJC_CAF <- subset(BLCA_IJC, seurat_clusters %in% c(4,5,9,13 ))  # COL1A1
BLCA_IJC_ImT <- subset(BLCA_IJC, seurat_clusters %in% c(1,10     ))  # PTPRC  c('CD2','CD3D','CD3E','CD3G','CD8A','CD8B','CD4','TNFRSF4','IL7R')
BLCA_IJC_ImM <- subset(BLCA_IJC, seurat_clusters %in% c(3        ))  # PTPRC  c('CD68','CD14','CD86','CD163','CSF1R','LYZ','FCGR3A','MS4A7')
BLCA_IJC_ImB <- subset(BLCA_IJC, seurat_clusters %in% c(11,12    ))  # PTPRC  c('CD24','CD79A') + c('CD9', 'CD63', 'CD81', 'CD82')

Throug.FVF.SD.RunPCA <- function(objs, nfeatures_n=2500){
  objs <- FindVariableFeatures(objs, selection.method = "vst", nfeatures = nfeatures_n)
  objs <- ScaleData(objs, features = rownames(objs))
  objs <- RunPCA(objs, features = VariableFeatures(object = objs))
  objs
}
BLCA_IJC_Epi <- Throug.FVF.SD.RunPCA(BLCA_IJC_Epi) # 上皮
BLCA_IJC_End <- Throug.FVF.SD.RunPCA(BLCA_IJC_End) # 内皮
BLCA_IJC_CAF <- Throug.FVF.SD.RunPCA(BLCA_IJC_CAF) # CAF + 肌细胞
BLCA_IJC_ImT <- Throug.FVF.SD.RunPCA(BLCA_IJC_ImT) # T细胞
BLCA_IJC_ImM <- Throug.FVF.SD.RunPCA(BLCA_IJC_ImM) # 巨噬细胞
BLCA_IJC_ImB <- Throug.FVF.SD.RunPCA(BLCA_IJC_ImB) # B细胞以及另一种

DimPlot(BLCA_IJC_Epi, reduction = "pca", group.by='Treat') + 
DimPlot(BLCA_IJC_End, reduction = "pca", group.by='Treat') + 
DimPlot(BLCA_IJC_CAF, reduction = "pca", group.by='Treat') + 
DimPlot(BLCA_IJC_ImT, reduction = "pca", group.by='Treat') + 
DimPlot(BLCA_IJC_ImM, reduction = "pca", group.by='Treat') + 
DimPlot(BLCA_IJC_ImB, reduction = "pca", group.by='Treat') 


DimHeatmap(BLCA_IJC_Epi, dims = 1:24, cells = 500, balanced = TRUE)

BLCA_IJC_Epi <- JackStraw(BLCA_IJC_Epi, dims=50, num.replicate = 100)
BLCA_IJC_Imm <- JackStraw(BLCA_IJC_Imm, dims=50, num.replicate = 100)
BLCA_IJC_End <- JackStraw(BLCA_IJC_End, dims=50, num.replicate = 100)
BLCA_IJC_CAF <- JackStraw(BLCA_IJC_CAF, dims=50, num.replicate = 100)
BLCA_IJC_Epi <- ScoreJackStraw(BLCA_IJC_Epi, dims = 1:50)
BLCA_IJC_Imm <- ScoreJackStraw(BLCA_IJC_Imm, dims = 1:50)
BLCA_IJC_End <- ScoreJackStraw(BLCA_IJC_End, dims = 1:50)
BLCA_IJC_CAF <- ScoreJackStraw(BLCA_IJC_CAF, dims = 1:50)

JackStrawPlot(BLCA_IJC_Epi, dims = 1:50) # 8
JackStrawPlot(BLCA_IJC_Imm, dims = 1:50) #16
JackStrawPlot(BLCA_IJC_End, dims = 1:50) # 2
JackStrawPlot(BLCA_IJC_CAF, dims = 1:50) # 5
ElbowPlot(BLCA_IJC_Epi, ndims = 50) + 
ElbowPlot(BLCA_IJC_Imm, ndims = 50) + 
ElbowPlot(BLCA_IJC_End, ndims = 50) + 
ElbowPlot(BLCA_IJC_CAF, ndims = 50) 

Through.FN.FC.umap.tsne <- function(objs, dims_n=1:16, dims_m=1:16, resolutions=0.6){
 objs <- FindNeighbors(objs, dims = dims_n)
 objs <- FindClusters(objs, resolution = resolutions)
 objs <- RunUMAP(objs, dims = dims_m)
 objs <- RunTSNE(objs, dims = dims_m) 
 objs
}
BLCA_IJC_Epi <- Through.FN.FC.umap.tsne(BLCA_IJC_Epi) # 上皮
BLCA_IJC_End <- Through.FN.FC.umap.tsne(BLCA_IJC_End) # 内皮
BLCA_IJC_CAF <- Through.FN.FC.umap.tsne(BLCA_IJC_CAF) # CAF + 肌细胞
BLCA_IJC_ImT <- Through.FN.FC.umap.tsne(BLCA_IJC_ImT) # T细胞
BLCA_IJC_ImM <- Through.FN.FC.umap.tsne(BLCA_IJC_ImM, dims_n=1:6,dims_m=1:6, resolutions=1) # 巨噬细胞
BLCA_IJC_ImB <- Through.FN.FC.umap.tsne(BLCA_IJC_ImB, dims_m=2) # B细胞以及另一种
#
DimPlot(BLCA_IJC_Epi, reduction = "tsne", label=T) + 
DimPlot(BLCA_IJC_End, reduction = "tsne", label=T) + 
DimPlot(BLCA_IJC_CAF, reduction = "tsne", label=T) + 
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
DimPlot(BLCA_IJC_ImM, reduction = "tsne", label=T) + 
DimPlot(BLCA_IJC_ImB, reduction = "tsne", label=T)  


FeaturePlot(BLCA_IJC_Epi, reduction = "tsne", features = c("PTPRC", "PECAM1", "EPCAM", "HRAS", "KRT13","KRT5","COL14A1", "ACTA2", "CD44", "ALDH1A3",'KRT18','KRT14'))
markers2 <- FindMarkers(BLCA_IJC_Epi, ident.1 = c(2,0,1,5), ident.2=4, min.pct = 0.25)

## 原文献CAF分三群
DimPlot(BLCA_IJC_CAF, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_CAF, reduction = "tsne", features = c('DCN','FAP','COL3A1','COL6A1','PDPN','COL1A2','MYH11','NF2F2','CRYAB','LMOD1','TPPP3','COL13A1','COL14A1'))
DimPlot(BLCA_IJC_CAF, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_CAF, reduction = "tsne", features = CellMarkers[['CAFs']])
DimPlot(BLCA_IJC_CAF, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_CAF, reduction = "tsne", features = CellMarkers[['MatrixCAFS']])
DimPlot(BLCA_IJC_CAF, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_CAF, reduction = "tsne", features = CellMarkers[['MyoCAFs']])
BLCA_IJC_CAF@meta.data$SubGroups <- 'CAF.NA'
BLCA_IJC_CAF@meta.data$SubGroups[which(BLCA_IJC_CAF@meta.data$RNA_snn_res.0.6 %in% c(2))]   <- 'CAF.COL13A1' 
BLCA_IJC_CAF@meta.data$SubGroups[which(BLCA_IJC_CAF@meta.data$RNA_snn_res.0.6 %in% c(4))]   <- 'CAF.COL14A1' 
BLCA_IJC_CAF@meta.data$SubGroups[which(BLCA_IJC_CAF@meta.data$RNA_snn_res.0.6 %in% c(1,3))] <- 'CAF.myoCAF' 
BLCA_IJC_CAF@meta.data$SubGroups[which(BLCA_IJC_CAF@meta.data$RNA_snn_res.0.6 %in% c(0))]   <- 'MuscleCell' 
## 原文献End仅一群，未细分，Endothelial 
DimPlot(BLCA_IJC_End, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_End, reduction = "tsne", features = c('PECAM1','VWF','CDH5','SELE','ENG','CD34'))
BLCA_IJC_End@meta.data$SubGroups <- 'Endothelial'
## T CELL  
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImT, reduction = "tsne", features = c('CD2','CD3D','CD3E','CD3G','CD8A','CD8B','CD4','TNFRSF4','IL7R'))
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImT, reduction = "tsne", features = CellMarkers[['Cytotoxicity']], order=T, max.cutoff=30)
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImT, reduction = "tsne", features = CellMarkers[['Regulatory']]) # Treg
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImT, reduction = "tsne", features = CellMarkers[['Exhaustion']]) # Treg
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImT, reduction = "tsne", features = CellMarkers[['Naive']])
DimPlot(BLCA_IJC_ImT, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImT, reduction = "tsne", features = c('CD2','CD3D','CD3E','CD3G','CD8A','CD8B','CD4','TNFRSF4','IL7R'))
BLCA_IJC_ImT@meta.data$SubGroups <- 'Tcells'
BLCA_IJC_ImT@meta.data$SubGroups[which(BLCA_IJC_ImT@meta.data$RNA_snn_res.0.6 %in% c(0))] <- 'naiveT' 
BLCA_IJC_ImT@meta.data$SubGroups[which(BLCA_IJC_ImT@meta.data$RNA_snn_res.0.6 %in% c(1))] <- 'CytoT' 
BLCA_IJC_ImT@meta.data$SubGroups[which(BLCA_IJC_ImT@meta.data$RNA_snn_res.0.6 %in% c(2))] <- 'Treg' 

## M1/M2
DimPlot(BLCA_IJC_ImM, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImM, reduction = "tsne", features = c('CD68','CD14','CD86','CD163','CSF1R','LYZ','FCGR3A','MS4A7'), order=T)
DimPlot(BLCA_IJC_ImM, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImM, reduction = "tsne", features = CellMarkers[[ 'M1' ]], order=T, max.cutoff=5)
DimPlot(BLCA_IJC_ImM, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImM, reduction = "tsne", features = CellMarkers[[ 'M2' ]], order=T)
DimPlot(BLCA_IJC_ImM, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImM, reduction = "tsne", features = CellMarkers[['Inflammatory']], order=T)
DimPlot(BLCA_IJC_ImM, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImM, reduction = "tsne", features = CellMarkers[['LC.LIKE']], order=T)
BLCA_IJC_ImM@meta.data$SubGroups <- 'Mac.NA'
BLCA_IJC_ImM@meta.data$SubGroups[which(BLCA_IJC_ImM@meta.data$RNA_snn_res.1 %in% c(3))]   <- 'M0' 
BLCA_IJC_ImM@meta.data$SubGroups[which(BLCA_IJC_ImM@meta.data$RNA_snn_res.1 %in% c(1))]   <- 'LC.like' 
BLCA_IJC_ImM@meta.data$SubGroups[which(BLCA_IJC_ImM@meta.data$RNA_snn_res.1 %in% c(0,2))] <- 'M2' 
# BLCA_IJC_ImB
# B cell   # 'CD9', 'CD63', 'CD81', 'CD82'
DimPlot(BLCA_IJC_ImB, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC, reduction = "tsne", features = c('CD24','CD79A'), order=T)
DimPlot(BLCA_IJC_ImB, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_ImB, reduction = "tsne", features = c('CD9', 'CD63', 'CD81', 'CD82'), order=T)
BLCA_IJC_ImB@meta.data$SubGroups <- 'B.NA'
BLCA_IJC_ImB@meta.data$SubGroups[which(BLCA_IJC_ImB@meta.data$RNA_snn_res.0.5 %in% c(12))]   <- 'ExosomeLikeT' 
BLCA_IJC_ImB@meta.data$SubGroups[which(BLCA_IJC_ImB@meta.data$RNA_snn_res.0.5 %in% c(11))]   <- 'Bcell' 

# Epi
DimPlot(BLCA_IJC_Epi, reduction = "tsne", label=T) + 
FeaturePlot(BLCA_IJC_Epi, reduction = "tsne", features = c("PTPRC", "PECAM1", "EPCAM", "HRAS", "KRT13","KRT5","COL14A1", "ACTA2", "CD44", "ALDH1A3",'KRT18','KRT14'))
BLCA_IJC_Epi@meta.data$SubGroups <- 'Epi.NA'
BLCA_IJC_Epi@meta.data$SubGroups[which(BLCA_IJC_Epi@meta.data$RNA_snn_res.0.6 %in% c(3))]       <- 'Epi.Lum' 
BLCA_IJC_Epi@meta.data$SubGroups[which(BLCA_IJC_Epi@meta.data$RNA_snn_res.0.6 %in% c(4))]       <- 'Epi.drug' 
BLCA_IJC_Epi@meta.data$SubGroups[which(BLCA_IJC_Epi@meta.data$RNA_snn_res.0.6 %in% c(0,1,2,5))] <- 'Epi.Bas' 


markers2 <- FindMarkers(BLCA_IJC_Epi, ident.1 = c(2,0,1,5), ident.2=4, min.pct = 0.25)


FinalCellGroups <- rbind(BLCA_IJC_Epi@meta.data[,c('MainGroups','SubGroups')],BLCA_IJC_CAF@meta.data[,c('MainGroups','SubGroups')],BLCA_IJC_End@meta.data[,c('MainGroups','SubGroups')],
                         BLCA_IJC_ImB@meta.data[,c('MainGroups','SubGroups')],BLCA_IJC_ImM@meta.data[,c('MainGroups','SubGroups')],BLCA_IJC_ImT@meta.data[,c('MainGroups','SubGroups')])
FinalCellGroups <- FinalCellGroups[colnames(BLCA_IJC),]
BLCA_IJC@meta.data$SubGroups <- FinalCellGroups$SubGroups

DimPlot(BLCA_IJC, reduction = "tsne", label=T)
save(BLCA_IJC,BLCA_IJC_CAF,BLCA_IJC_End,BLCA_IJC_Epi,BLCA_IJC_ImB,BLCA_IJC_ImM,BLCA_IJC_ImT,FinalCellGroups,
     CellMarkers,Expr1,Expr2,Expr3,Meta1,Meta2,Meta3,PDXHOMO,PDXMICE,Primary,refGTF.Ens100,ssGSEA.BLCA.Markers,
     Throug.FVF.SD.RunPCA,Through.FN.FC.umap.tsne    , file='GSE130001.Seurat.Clustering.RData')

### 以上，对该样本的分群分析结束

VlnPlot(BLCA_IJC, features = c("MED7"), slot = "counts", log = TRUE)
FeaturePlot(BLCA_IJC, reduction = "tsne", features = c('MED7'), order=T, max.cutoff=2)


cds <- as.cell_data_set(BLCA_IJC_Epi)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)








