mye <- subset(sce.all,celltype%in%c("Macro/Mono","DC"))
mye_col <- c("#1f77b4ff",'#aec7e8ff','#98df8aff',"#2ca02cff","#ffbb78ff","#ff7f0eff","#c5b0d5ff","#9467bdff","#AFAFAD")
mye <- subset(mye, percent.mt <20&nFeature_RNA<4000&nCount_RNA<20000&nFeature_RNA>500)
mye <- mye %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData()%>%
  RunPCA(npcs = 50)
mye <- mye%>%
  RunHarmony(reduction = "pca",dims.use = 1:20,group.by.vars = "orig.ident")%>%
  RunUMAP(reduction = "harmony",dims = 1:20)%>%
  FindNeighbors(reduction = "harmony",dims = 1:20)%>%
  FindClusters(resolution = seq(0.1,1,0.1))
mye$seurat_clusters <- mye$RNA_snn_res.0.4
#
Idents(mye)<-"seurat_clusters"
mye_marker <- FindAllMarkers(mye,only.pos = T)
myetop10 <- mye_marker%>%subset(p_val_adj<0.05)%>%group_by(cluster)%>%top_n(10,wt = avg_log2FC)
FeaturePlot(mye,features = c("S100A8","LYVE1","SPP1","CD1C","MARCO","FCGR3B","CD14","FCGR3A"),reduction = "tsne")
DotPlot(mye,features = c("S100A8","S100A9","FCN1","VCAN",
                         "FCGR3A","LILRA1","CDKN1C",
                         "THBS1","CLEC5A","CCL2",
                         "SPP1","LGMN","GPNMB",
                         "NR5A2","PTPRG","MAGI1",
                         "C1QC","MARCO","FN1",
                         "ATG7","LRMDA","DOCK4",
                         "CD1C","FCER1A","HLA-DQA1",
                         "TOP2A","MKI67",'PCNA'),cols = c("RdYlBu"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
  theme_bw()+labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text = element_text(size = 12,color = "black"))+
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  theme(text = element_text(size = 14))+coord_fixed()

mye@meta.data$mye_celltype <- NA
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(0))] <- "CD14_Mono"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(9))] <- "CD16_Mono"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(3))] <- "cDC"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(13))] <- "Cycling"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(1,4))] <- "SPP1_Mac"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(2,5))] <- "THBS1_Mac"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(11))] <- "NR5A2_Mac"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(6))] <- "C1QC_Mac"
mye@meta.data$mye_celltype[which(mye@meta.data$seurat_clusters%in%c(12))] <- "ATG7_Mac"
mye@meta.data$mye_celltype <- factor(mye@meta.data$mye_celltype,levels = c("CD14_Mono","CD16_Mono","THBS1_Mac","SPP1_Mac","NR5A2_Mac","C1QC_Mac","ATG7_Mac","cDC","Cycling"))
Idents(mye) <- "mye_celltype"
DimPlot(mye,cols = group_col)+theme_void()+coord_fixed()

mye_geneset <- read.csv("geneset/mye Signature Genes.csv")
mye <- AddModuleScore(mye,features = list(mye_geneset$M2[1:34]),name = "M2")
mye <- AddModuleScore(mye,features = list(mye_geneset$M1[1:16]),name = "M1")
mye <- AddModuleScore(mye,features = list(mye_geneset$Angiogenesis[1:25]),name = "Angiog")
mye <- AddModuleScore(mye,features = list(mye_geneset$Phagocytosis[1:4]),name = "Phago")
DotPlot(mye,features = list(M1=c('CD86','IL1B','KYNU','IRF1','IRF5','CD40','TNF'),
                            M2=c('CLEC7A',"CSTA",'CSTB','CSTC','CSTD','MSR1','CCL4','IL4R','TGFB1','FN1','LYVE1','VEGFB','MMP19','TNFSF8','TNFSF12'),
                            Angiogenesis= c("CD44","CXCR4","VCAN","FYN","VEGFA","EZH2","CCND2","E2F3","ITGAV","SPP1"),
                            Phagocytosis=mye_geneset$Phagocytosis[1:4]))+labs(x='',y='')+
  scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
#monocle
plot_cell_trajectory(cds_DGT, cell_size = 2.2, color_by = "Pseudotime",show_branch_points = T)+
  scale_color_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                    '#bef0b0', '#fdf4af', '#f9b64b',
                                    '#ec840e', '#ca443d', '#a51a49'), 
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = NA))

        
                                   
