library(monocle)
gc1 <- subset(sce.all,celltype=="Granulosa1")
seurat_to_monocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- GetAssayData(otherCDS, assay = assay, slot = slot)
    data <- data[rowSums(as.matrix(data)) != 0,]
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      } else {
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    } 
  } 
  return(monocle_cds)
}

monocle_fun <- seurat_to_monocle(gc1, assay = "RNA", slot = "counts")

#Estimate size factors and dispersions
monocle <- estimateSizeFactors(monocle_fun)
monocle <- estimateDispersions(monocle)

monocle <- detectGenes(monocle, min_expr = 0.1)
print(head(fDatamonocle)))
print(head(pData(monocle)))
expressed_genes <- row.names(subset(fData(monocle), num_cells_expressed >= 10))
pData(monocle)$Total_mRNAs <- Matrix::colSums(exprs(monocle))
monocle <- monocle[,pData(monocle)$Total_mRNAs < 1e6]

cds_DGT <- monocle
diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~group")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]
cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
plot_ordering_genes(cds_DGT)
save(cds_DGT, file = "cds_DGT.RData")
cds_DGT<- reduceDimension(cds_DGT, max_components = 2,reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT, root_state = 5, num_paths = NULL, reverse = F)
plot_cell_trajectory(cds_DGT, cell_size = 2.2, color_by = "group",show_branch_points = F)+
  facet_wrap(~group, nrow = 2)+
  scale_color_manual(values = group_col)+NoAxes()+coord_fixed()
df <- pData(cds_DGT)
ggpie(df[df$State%in%c(5,6),],group)+scale_fill_manual(values = group_col)
ggpie(df[df$State%in%c(1,2,3),],group)+scale_fill_manual(values = group_col)
ggpie(df[df$State%in%c(7),],group)+scale_fill_manual(values = group_col)
BEAM_res <- BEAM(cds_DGT, branch_point = 2, cores = 1,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res1 <-subset(BEAM_res,qval <1e-3)
library(ClusterGVis)
beam_plot<-plot_genes_branched_heatmap2(cds_DGT[row.names(BEAM_res1),],
                                        branch_point = 2,
                                        num_clusters = 3,
                                        cores = 1,
                                        use_gene_short_name = T,
                                        show_rownames = F)
visCluster(object = beam_plot,plot.type = "both",ctAnno.col = c("#F64B35FF","#4DBBD5FF","#00A087FF"))
library(org.Hs.eg.db)
enrich <- enrichCluster(object = beam_plot,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 4,
                        seed = 5201314,
                        add.gene = T)
enrich$logpvalue <- -log10(enrich$pvalue)
ggbarplot(enrich, y = "logpvalue", x = "Description",
          fill = "group", color = "group",rotate = T,
          xlab = "Cluster", ylab = "-log10(P-value)",
          palette = c("#F64B35FF","#4DBBD5FF","#00A087FF"))
#scFEA
human_moduleInfo <- read.csv("../scFEA/data/scFEA.human.moduleinfo.csv", header = T, row.names = 1)
adj_flux <- read.csv("../scFEA/output/gc_fluxres.csv",row.names = 1)
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))

df_averages <- as.data.frame(t(adj_flux)) %>%
  group_by(group = human_moduleInfo$SM_anno) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  column_to_rownames("group")
df_averages <- df_averages[,rownames(pData(cds_DGT))]
df_averages <- as.data.frame(t(df_averages))
pData(cds_DGT)$Glycolysis_TCA_cycle <- df_averages$Glycolysis_TCA_cycle
pData(cds_DGT)$Steroid_hormone_synthesis<-df_averages$Steroid_hormone_synthesis
plot_cell_trajectory(cds_DGT,color_by = "Glycolysis_TCA_cycle",show_branch_points = F)+scale_color_gsea()+NoAxes()+coord_fixed()
df_sub<-df[df$State%in%c(5,6,1,2,3),]
df_sub$tra_pathway <- "Trajectory2"
df_sub1<-df[df$State%in%c(5,6,7),]
df_sub1$tra_pathway <- "Trajectory1"
df_sub<-rbind(df_sub,df_sub1)
ggplot(df_sub,aes(x=Pseudotime,y=Steroid_hormone_synthesis,group = tra_pathway,color = tra_pathway))+
  geom_smooth(se=T)+
  scale_color_manual(values = c("#5cb3ba","#ff5e00"))+scale_fill_manual(values = c("#5cb3ba","#ff5e00"))+
  theme_bw() +
  theme(axis.title = element_text(size = 12,color ="black"), 
        axis.text = element_text(size= 12,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)
  )
