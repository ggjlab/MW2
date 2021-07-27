setwd("/media/ggj/tom/Hider/jiace/")
##read dge file from the folder
temp<- list.files(pattern="*.dge.txt.gz")
name<- character()
for (i in 1:length(temp)) {
  message("loading DGE")
  name[i]<-unlist(strsplit(temp[i],"_dge"))[1]
  tempvalue<-read.table(temp[i],header = T,row.names = 1)
  assign(name[i],tempvalue)
  message(paste(name[i],"is loaded"))
  
}

##add annotation to the cell barcode
for(i in 1:length(temp)){
  dge<-get(name[i])
  colnames(dge)<-paste0(as.character(name[i]),".",colnames(get(name[i])))
  assign(name[i],dge)
}

###build dge more than 500 UMI
name_300more<-name
name_300more<-paste0(name_300more,sep="_","300more")
for (i in 1:length(name)) {
  dge<-get(name[i])
  temp<-dge[,colSums(dge)>=300]
  assign(name_300more[i],temp)
  message(paste(name[i],"is done"))
}

###build dge more than 300 UMI
name_300more<-name
name_300more<-paste0(name_300more,sep="_","300more")
for (i in 1:length(name)) {
  dge<-get(name[i])
  temp<-dge[,colSums(dge)>=300]
  assign(name_300more[i],temp)
  message(paste(name[i],"is done"))
}


#merge  (if you want merge more than one sample data)

combined<-merge(COL01_300more,COL43_300more,by = 'row.names',all = T)
rownames(combined)=combined$Row.names
combined<-combined[,-1]
combined<-merge(combined,COL03_500more,by = 'row.names',all = T)
rownames(combined)=combined$Row.names
combined<-combined[,-1]
combined<-merge(combined,COL04_500more,by = 'row.names',all = T)
rownames(combined)=combined$Row.names
combined<-combined[,-1]
combined[is.na(combined)] = 0

######################seurat
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

pbmc_01<-CreateSeuratObject(counts = COL01_300more, min.cells = 3,  
                            project = "COL01_300more")
pbmc_03<-CreateSeuratObject(counts = COL03_300more, min.cells = 3,  
                            project = "COL03_300more")
pbmc_41<-CreateSeuratObject(counts = COL41_300more, min.cells = 3,  
                            project = "COL41_300more")
pbmc_43<-CreateSeuratObject(counts = COL43_300more, min.cells = 3,  
                            project = "COL43_300more")
pbmc <- merge(pbmc_01, y = c(pbmc_03, pbmc_41,pbmc_43), project = "PBMC")

Idents(pbmc) = pbmc$orig.ident
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@assays$RNA), value = TRUE)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
par(mfrow = c(1, 2))
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, )
hv.genes <- head(VariableFeatures(pbmc),2000)
#pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
#                          x.low.cutoff = 0.015, y.cutoff = 0.4)
length(VariableFeatures(pbmc))

##### scale data
pbmc <- ScaleData(object =pbmc,features = hv.genes, vars.to.regress = c("nCount_RNA"))
#pbmc <- ScaleData(object =pbmc,features = hv.genes, vars.to.regress = c("percent.mt","nCount_RNA"))

##### runPCA
pbmc <- RunPCA(object =pbmc, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
               pcs.print = 1:5, genes.print = 5,npcs = 50)
ElbowPlot(object =pbmc, ndims = 50)#1:25

# run TSNE and FindCluster
pbmc <- FindNeighbors(pbmc,dims=1:20)
pbmc <- FindClusters(pbmc, resolution = 0.15)
pbmc <- RunTSNE(object =pbmc, reduction = "pca", dims= 1:20,
                reduction.name = "tSNE") 
pbmc <- RunUMAP(object =pbmc, reduction = "pca", dims= 1:20,
                reduction.name = "umap") 
pbmc$seurat_clusters = pbmc$RNA_snn_res.0.15
pbmc$seurat_clusters = as.integer(pbmc$seurat_clusters)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(length(levels(as.factor(pbmc$seurat_clusters))))
DimPlot(object = pbmc, reduction= "tSNE",label = T, pt.size = 1,label.size = 7,cols = col_flg) 
DimPlot(object = pbmc, reduction= "umap",group.by = 'seurat_clusters',label = F, pt.size = 0.1,label.size = 7,cols = col_flg) 

##### match barcode4 and treatments
cellname = colnames(pbmc)
library(stringi)
barcode4 = substr(cellname,25,34)
barcode4 = as.data.frame(barcode4)
rownames(barcode4) = cellname
pbmc@meta.data$barcode4 = barcode4$barcode4
metadata=pbmc@meta.data
a = table(barcode4)
a = rownames(a)
saveRDS(pbmc,file = './pbmc.RDS')
cluster_barcode = matrix(data = 0,nrow = 96,ncol = 5)
rownames(cluster_barcode)=a
for (i in 1:96){
  b = metadata[metadata$barcode4==a[i],]
  for (j in 1:5){
    cluster_barcode[i,j]=sum(b$seurat_clusters==j)
  }
}
cluster_barcode = as.data.frame(cluster_barcode)
cluster_barcode$barcode=rownames(cluster_barcode)
index=read.csv('../barcode4_index.csv',header = 0)
colnames(index)=c('index','barcode')
rownames(index)=index$barcode
library(plyr)
cluster_barcode = join(index,cluster_barcode,by='barcode')
write.csv(cluster_barcode,file = './cluster_barcode.csv',row.names = F)

##### annotate treatments with clustering result
a = read.csv('../merge/cluster_barcode.csv',row.names = 1)
anno = as.data.frame(pbmc$barcode4)
a$name = c('CH+SB','CH+SB','CH+RA+SB','CH+RA+SB','CH+LD+SB','CH+LD+SB','SB','SB',
           'CH+LD','CH+LD','CH+RA+LD','CH+RA+LD','CH+LD+P1+RA','CH+LD+P1+RA','LD','LD',
           'CH+IW','CH+IW','CH+RA+IW','CH+RA+IW','CH+LD+IW','CH+LD+IW','IW','IW',
           'CH+PU','CH+PU','CH+RA+PU','CH+RA+PU','CH+LD+PU','CH+LD+PU','PU','PU',
           'CH+CY','CH+CY','CH+RA+CY','CH+RA+CY','CH+LD+CY','CH+LD+CY','CY','CY',
           'CH+P1','CH+P1','CH+RA+P1','CH+RA+P1','CH+LD+P1','CH+LD+P1','P1','P1',
           'CH+P0','CH+P0','CH+RA+P0','CH+RA+P0','CH+LD+P0','CH+LD+P0','P0','P0',
           'CH+DA','CH+DA','CH+RA+DA','CH+RA+DA','CH+LD+DA','CH+LD+DA','DA','DA',
           'CH+FO','CH+FO','CH+RA+FO','CH+RA+FO','CH+LD+FO','CH+LD+FO','FO','FO',
           'CH+HA','CH+HA','CH+RA+HA','CH+RA+HA','CH+LD+HA','CH+LD+HA','HA','HA',
           'CH+RA','CH+RA','SA','SA','XA','XA','CH','CH',
           'empty','empty','LY','LY','RE','RE','RA','RA')
a$ident =   c(1,1,1,1,1,1,2,2,
              1,1,1,1,4,4,2,2,
              1,1,1,1,1,1,2,2,
              1,1,1,1,1,1,2,2,
              1,1,1,1,1,1,2,2,
              3,3,4,4,3,3,5,5,
              3,3,4,4,3,3,5,5,
              1,1,1,1,1,1,2,2,
              1,1,1,1,1,1,2,2,
              1,1,1,1,1,1,2,2,
              1,1,2,2,2,2,1,1,
              2,2,2,2,2,2,2,2
)
colnames(anno) = 'barcode'
anno$cellname = rownames(anno)
anno = join(anno,a,by = 'barcode')
write.csv(anno,file='../jiace/anno.csv')
pbmc@meta.data$ident = anno$ident
Idents(pbmc) = pbmc@meta.data$name
pbmc@meta.data$ident = as.character(pbmc@meta.data$ident)
pbmc$name = anno$name

##### Find and plot markers by different clusters
pbmc_3 = subset(x=pbmc,idents = c(2,3))
aa = FindAllMarkers(pbmc_3,only.pos = TRUE, min.pct = 0.25, 
                    logfc.threshold  = 0.25)
write.csv(aa,file = './merge/cluster3_marker.csv')

##### Generate trajectory analysis dataset for several specific treatments
library(SeuratData)
library(SeuratDisk)
pbmc_1 = subset(x=pbmc,idents = c('P0','P1','CH+P0','CH+P1','CH+RA+P0','CH+RA+P1'))
SaveH5Seurat(pbmc_1, filename = "pbmc_1.h5Seurat")
Convert("pbmc_1.h5Seurat", dest = "h5ad")

##### Plot trajectory analysis marker genes
gene_names = c('FAM20A','KRT8','SERPINB9','KRT18','STOM','L1TD1','FRAT2','TFAP2A','SEMA6A','TPM1','ID1','HAPLN1',
               'TP53I11','FGFR3','EZR','DNMT3B','TRIM24','PLS3','SLC16A1','AUTS2','SESN3','MALAT1','EEF2','GLUL',
               'PRRC2C','RAB3B','TRIM71','FUS','PTPRF','FGFR1',
               'USP9X','AASS','UGP2','MFGE8','ACTG1','PTMA','LITAF','ZFAND5','ESRG','FLNA',
               'FST', 'NKD1', 'VAT1L','VCAN','VIM','CCND1','CD24','SOX2',
               'TUBB2A','TUBB2B','DDIT4','PGK1',
               'PRTG', 'PCAT14', 'MSX1','ZNF703','SKAP2','HOXA1','DSP','ARL4C','NR6A1','MT-CO3',
               'SLC2A3','PDLIM1','IGFBP2')
top10 <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
Idents(pbmc_1)=pbmc_1$name
DoHeatmap(pbmc_1, features = top10$gene)
pbmc_1$name <- factor(x =pbmc_1$name, levels = c('P1','P0','CH+P1','CH+P0','CH+RA+P1','CH+RA+P0'))