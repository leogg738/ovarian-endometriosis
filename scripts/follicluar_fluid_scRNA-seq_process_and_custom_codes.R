rm(list = ls())
library(Seurat)
library(harmony)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(patchwork) 
library(ggpubr)
library(data.table)
library(RColorBrewer) 
library(decontX)
library(DoubletFinder)

dir='T7_LBSW-20231031-ScRNA-DJ-SH_6sample/'
samples=list.files(dir)
sceList = lapply(samples,function(pro){
  #pro=samples[1]
  print(pro)
  sce=CreateSeuratObject(counts =  Read10X(file.path(dir,pro) ) ,
                         project =  gsub('.matrix','',pro),
                         min.cells = 3,
                         min.features = 200)
  return(sce)
})
names(sceList)
samples = c("EMS-ZL","EMS-JJT","EMS-XJW","CON-LL","CON-XLY","CON-ZJ")
names(sceList) = samples
#estimate RNA contamination  
for(i in 1:length(sceList)){
  counts <- sceList[[i]]@assays$RNA@counts
  decontX_res <- decontX(counts)
  sceList[[i]]@meta.data$contamination =decontX_res$contamination
}
#estimate doublet
for (i in seq_along(sceList)) {
  sceList[[i]] <-sceList[[i]]%>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst",nfeatures = 2000)%>% 
    ScaleData(verbose = FALSE)%>%RunPCA(verbose = FALSE)%>%RunUMAP(dims = 1:20)%>% 
    FindNeighbors(reduction = 'pca', dims = 1:20)%>%FindClusters(resolution = 0.4) 
}
ks_detectDoublet <- function(obj,
                             dims,
                             estDubRate, 
                             ncores=1,
                             Homotypic=F,
                             annotation){
  #use DoubletFinder packages
  require(DoubletFinder)#2.0.4
  
  #select pK
  sweep.res.list <- paramSweep_v3(obj, PCs=dims,num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))
  
  
  #Doublet Proportion Estimate
  if(Homotypic==F){
    
    nExp_poi <- round(estDubRate * length(Cells(obj)))
    
  }else{
    
    homotypic.prop <- modelHomotypic(obj@meta.data[,annotation])
    nExp_poi <- round(estDubRate * length(Cells(obj)))
    nExp_poi  <- round(nExp_poi*(1-homotypic.prop))
  }
  
  # DoubletFinder:
  obj <- doubletFinder_v3(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE)
  
  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL
  
  return(obj)
}
for(i in seq_along(sceList)){
  sceList[[i]] <-ks_detectDoublet(sceList[[i]],dims = 1:20,estDubRate=0.05,
                                   ncores = 2,
                                   Homotypic=F, annotation="seurat_clusters")
}
#quality control
sce.all <- merge(sceList[[1]],y = sceList[-1],add.cell.ids = samples)

sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all,pattern = "^MT-")
sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all, pattern = "^HB[^(P)]")
sce.all <- CellCycleScoring(sce.all, s.features = cc.genes.updated.2019$s.genes, 
                            g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)

sce.all <- subset(sce.all, subset = 
                 nFeature_RNA > 200 & 
                 percent.mt < 35 &
                 DF.classify =="Singlet"&
                 contamination < 0.2)
sce.all <-sce.all%>%
  RunHarmony(reduction = "pca",dims.use = 1:25,group.by.vars = "orig.ident")%>%
  RunUMAP(reduction = "harmony",dims = 1:25)%>%
  FindNeighbors(reduction = "harmony",dims = 1:25)%>%
  FindClusters(resolution = seq(0.1,1,0.1))

clustree(sce.all)
DimPlot(sce.all,label = T)

sce.all$seurat_clusters <-sce.all$RNA_snn_res.0.3
Idents(sce.all)<-"seurat_clusters"
#annotation
sce.all@meta.data$celltype<-NA
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(0,2,10))]<-"Neutrophils"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(7))]<-"Granulosa1"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(1,5))]<-"Granulosa2"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(3))]<-"Granulosa3"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(4))]<-"Macro/Mono"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(8))]<-"DC"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(6))]<-"T/NKT"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(9))]<-"B"
sce.all@meta.data$celltype[which(sce.all@meta.data$seurat_clusters%in%c(11))]<-"Stromal"
sce.all@meta.data$celltype <- factor(sce.all@meta.data$celltype,levels = 
                                           c("Granulosa1","Granulosa2","Granulosa3","Stromal","Neutrophils","Macro/Mono","DC","T/NKT","B"))

DotPlot(sce.all,features = c("SERPINE2","GJA1","FOXL2","AMH","FSHR",
                             "STAR","CYP11A1","CDH2","FOXO1",
                             "DCN","LUM","COL1A1",
                             "PTPRC",
                             "CXCL8","FCGR3B","CSF3R",
                             "LYZ","CD163","CD86",
                             "CD1C","FCER1A","CLEC10A",
                             "CD3E","KLRD1","CD79A"
                             ),cols = c("RdYlBu"))+
  theme_bw()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size = 12))+
  theme(
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank())+
  theme(text = element_text(size = 14))+
  labs(x="",y="")

celltype_col <- c("#F98962","#ffa726","#ffcc80","#D62728", "#C5B0D5", "#1976d2","#bbdefb","#66bb6a","#98DF8A")
Idents(sce.all)<-"celltype"
DimPlot(sce.all,group.by = "celltype",cols = celltype_col)+NoAxes()+coord_fixed()
#celltype proportion
Ratio <- sce.all@meta.data %>%group_by(sample,celltype) %>%
  count() %>%
  group_by(sample) %>%
  mutate(Freq = n/sum(n)*100)
ggplot(Ratio, aes(x = sample, y = Freq, fill = celltype))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = celltype_col)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size = 12)
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="",y="Percentage(%)")+coord_flip()
#R o/e
meta<-sce.all@meta.data[,c("sample",'celltype')];
data<-table(meta[,c("sample",'celltype')])
data<-t(data)
res<-chisq.test(data)
ratio<-res$observed/res$expected;
pheatmap(ratio,
              color = c(colorRampPalette(c('#313695','white'))(30),
                        colorRampPalette(c('white','#ff7f00'))(70)),
              cluster_rows = F,
              cluster_cols = F,fontsize = 12)

#DEG
group_deg <-list()
celltype <- c("Granulosa1","Granulosa2","Granulosa3","Stromal","Neutrophils",
              "Macro/Mono","DC","T/NKT","B")
for(i in 1:length(celltype)){
  group_deg[[i]] <- FindMarkers(sce.all,ident.1 = "EMS",ident.2 = "CON",
                                group.by = "group",subset.ident = celltype[i])
  group_deg[[i]]$cluster <- celltype[i]
  group_deg[[i]]$gene <- rownames(group_deg[[i]])
  names(group_deg)[i]<-celltype[i]
}
library(scRNAtoolVis)
deg_table <- do.call(rbind,group_deg)
deg_table$cluster <- factor(deg_table$cluster,levels = 
                              c("Granulosa1","Granulosa2","Granulosa3","Stromal",
                                "Neutrophils","Macro/Mono","DC","T/NKT","B"))
a =  deg_table[deg_table$p_val_adj<0.05,]
markerVocalno(markers = a,
              topn = 3,
              labelCol = celltype_col[-4])
#Enrichment analysis
for(i in 1:length(group_deg)){
  group_deg[[i]] <- subset(group_deg[[i]],p_val_adj<0.05)
}
group_keggp <-list()
for(i in 1:length(group_deg)){
  group_keggp[[i]] <- bitr(group_deg[[i]][group_deg[[i]]$avg_log2FC>0,]$gene,"SYMBOL","ENTREZID",OrgDb = 'org.Hs.eg.db')
  group_keggp[[i]] <- enrichKEGG(gene = group_keggp[[i]]$ENTREZID,
                               keyType = "kegg",organism = "hsa",
                               pvalueCutoff  = 1
                               )
  names(group_keggp)[i] <- names(group_deg)[i]
  group_keggp[[i]]<-setReadable(group_keggp[[i]],OrgDb = 'org.Hs.eg.db',keyType = "ENTREZID")
}
group_keggn <-list()
for(i in 1:length(group_deg)){
  group_keggn[[i]] <- bitr(group_deg[[i]][group_deg[[i]]$avg_log2FC<0,]$gene,"SYMBOL","ENTREZID",OrgDb = 'org.Hs.eg.db')
  group_keggn[[i]] <- enrichKEGG(gene = group_keggn[[i]]$ENTREZID,
                                 keyType = "kegg",organism = "hsa",
                                 pvalueCutoff  = 1
  )
  names(group_keggn)[i] <- names(group_deg)[i]
  group_keggn[[i]]<-setReadable(group_keggn[[i]],OrgDb = 'org.Hs.eg.db',keyType = "ENTREZID")
}
##iron metabolism related genes expression
sce.all@meta.data$group_celltype <- paste(sce.all@meta.data$celltype,sce.all@meta.data$group,sep = "_")
DotPlot(sce.all,features = c("FTH1","FTL","TFRC",
                             "SLC11A2","SLC39A8","SLC39A14","CD44",
                             "STEAP1","STEAP2","STEAP3","STEAP4",
                             "LCN2",
                             'SLC40A1','HAMP',
                             'HMOX1','GLRX3','SLC7A11','SOD1','SOD2'),
        group.by = "group_celltype",assay = "RNA")+RotatedAxis()+coord_flip()+
  scale_color_gradient2(low = "navy",high = "firebrick",mid = "white")+
  theme(axis.text = element_text(size= 12,color = "black"),
        axis.text.x = element_text(angle = 90),
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)
  )+labs(x="",y="")
#cellchat
library(CellChat)
EMS <- sce.all%>%subset(group=="EMS")%>%
  subset(celltype%in%c("Granulosa1",  "Granulosa2",  "Granulosa3" , "Neutrophils", "Macro/Mono" , "DC","T/NKT","B"))
CON <- sce.all%>%subset(group=="CON")%>%
  subset(celltype%in%c("Granulosa1",  "Granulosa2",  "Granulosa3" , "Neutrophils", "Macro/Mono" , "DC","T/NKT","B"))
CON@meta.data$celltype <- as.character(CON@meta.data$celltype)
CON@meta.data$celltype <- factor(CON@meta.data$celltype,levels = c("Granulosa1",  "Granulosa2",  "Granulosa3" , "Neutrophils", "Macro/Mono" , "DC","T/NKT","B"))
CON_cellchat <- createCellChat(CON@assays$RNA@data,meta = CON@meta.data,group.by = "celltype")
levels(cellchat@idents)

CellChatDB <- CellChatDB.human
CON_cellchat@DB <- CellChatDB

CON_cellchat <- subsetData(CON_cellchat) 
CON_cellchat <- identifyOverExpressedGenes(CON_cellchat)
CON_cellchat <- identifyOverExpressedInteractions(CON_cellchat)

CON_cellchat <- computeCommunProb(CON_cellchat, type = "triMean")
CON_cellchat <- computeCommunProbPathway(CON_cellchat)
CON_cellchat <- aggregateNet(CON_cellchat)

CON.net <- subsetCommunication(CON_cellchat)
EMS.net <- subsetCommunication(EMS_cellchat)

cellchat_list <- list(CON = CON_cellchat,EMS=EMS_cellchat)
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat))
