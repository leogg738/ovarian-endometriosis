library(data.table)
library(Seurat)

data_dir <- "samples"

files <- list.files(data_dir, full.names = TRUE)

file_prefixes <- unique(gsub("-.*", "", basename(files)))

for (prefix in file_prefixes) {
  dest_dir <- file.path(data_dir, prefix)
  dir.create(dest_dir, showWarnings = FALSE)
  prefix_files <- grep(prefix, files, value = TRUE)
  for (file in prefix_files) {
    new_name <- gsub(paste0("^", prefix, "-"), "", basename(file))
    file.rename(file, file.path(dest_dir, new_name))
  }
}

folders <- c("GSM8077816_young", "GSM8077817_young", "GSM8077818_young", 
             "GSM8077819_middle", "GSM8077820_middle", "GSM8077821_middle", 
             "GSM8077822_old", "GSM8077823_old", "GSM8077824_old")

for (folder in folders) {
  files <- list.files(folder, full.names = TRUE)

  for (file in files) {
    basename <- basename(file)
    if (grepl("-", basename)) {
      new_name <- sub("^[^-]*-", "", basename)

      new_file <- file.path(folder, new_name)
      file.rename(file, new_file)
    }
  }
}

lapply(folders, list.files)

rm(list = ls())

dir <- "./samples"
samples <- list.files(dir)
library(data.table)
sceList <- lapply(samples, function(pro) {
  print(pro)
  sce <- CreateSeuratObject(
    counts = Read10X(file.path(dir, pro)),
    project = gsub('./samples', '', pro)
  )
  return(sce)
})
library(decontX)
for(i in 1:length(sceList)){
  counts <- sceList[[i]]@assays$RNA@counts
  decontX_res <- decontX(counts)
  sceList[[i]]@meta.data$contamination =decontX_res$contamination
}
names(sceList)
samples 
names(sceList) = samples
sce.all <- merge(sceList[[1]],y = sceList[-1],add.cell.ids = samples)

sce.all[["percent.mt"]]<-PercentageFeatureSet(sce.all,pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all,pattern = "^RP[SL]")
VlnPlot(sce.all,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size = 0,group.by = "orig.ident")
VlnPlot(sce.all,features = c("percent.rp"),pt.size = 0,group.by = "orig.ident")
sce.all <- subset(sce.all, subset = 
                    nFeature_RNA < quantile(sce.all@meta.data$nFeature_RNA, 0.99)&
                    nFeature_RNA > 500&
                    nCount_RNA <quantile(sce.all@meta.data$nCount_RNA, 0.99)&
                    percent.mt < 20&
                    percent.rp <60)

sce_list <- SplitObject(sce.all, split.by = "orig.ident")
for (i in seq_along(sce_list)) {
  sce_list[[i]] <-sce_list[[i]]%>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst",nfeatures = 2000)%>% 
    ScaleData(verbose = FALSE)%>%RunPCA(verbose = FALSE)%>%RunUMAP(dims = 1:20)%>% 
    FindNeighbors(reduction = 'pca', dims = 1:20)%>%FindClusters(resolution = 0.4) 
}
library(DoubletFinder)
for(i in seq_along(sce_list)){
  sce_list[[i]] <-ks_detectDoublet(sce_list[[i]],dims = 1:20,estDubRate=0.05,
                                   ncores = 2,
                                   Homotypic=F, annotation="seurat_clusters")
}
all_exp <-  as.data.frame(AverageExpression(sce.all,group.by = "celltype_group"))
celltype_exp <- list()
for (i in 1:7) {
  start_col <- (i - 1) * 3 + 1
  end_col <- i * 3
  celltype_exp[[i]] <- all_exp[, start_col:end_col]
  names(celltype_exp)[i] <- str_split(colnames(celltype_exp[[i]])[1],pattern = "_",simplify = T)[,1]
  colnames(celltype_exp[[i]]) <- str_split(colnames(celltype_exp[[i]]),pattern = "_",simplify = T)[,2]
  celltype_exp[[i]] <- celltype_exp[[i]][,c("young","middle","old")]
  celltype_exp[[i]] <- as.matrix(celltype_exp[[i]])
}
library(vegan)
vegan <- cascadeKM(celltype_exp[[4]], inf.gr = 2, sup.gr = 10, iter = 1000)
plot(vegan, sortg=TRUE,grpmts.plot = TRUE)
VlnPlot(sce.all,"SLC40A1",split.by = "group",pt.size = 0)

mfuzz_class <- new('ExpressionSet',exprs = celltype_exp[[3]])
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
mfuzz_class <- standardise(mfuzz_class)
set.seed(1234)
mfuzz_cluster <- mfuzz(mfuzz_class, c = 8, m = mestimate(mfuzz_class))
mfuzz.plot2(mfuzz_class, cl=mfuzz_cluster,mfrow=c(2,4),centre=TRUE,x11=F,centre.lwd=0.2,
            time.labels = c("Y","M","O"))
protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- as.data.frame(cbind(protein_standard[names(protein_cluster), ], protein_cluster))
protein_standard_cluster$gene <- rownames(protein_standard_cluster)
protein_cluster <- split(protein_standard_cluster[,4:5],protein_cluster)
protein_cluster_iron <- list()
for(i in seq_along(protein_cluster)){
  protein_cluster_iron[[i]] <-protein_cluster[[i]][protein_cluster[[i]]$gene%in%ironhomeo_gene,]
}
library(org.Hs.eg.db)
protein_go <- list()
for(i in 1:8){
  protein_go[[i]] <- enrichGO(protein_cluster[[i]]$gene,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP",qvalueCutoff = 0.05)
}
library(viridis)
plotdata <- protein_go[[5]]@result%>%subset(ID%in%c("GO:0001819","GO:0048146","GO:0050727","GO:2001233","GO:0046058"))
plotdata <- protein_go[[1]]@result%>%subset(ID%in%c("GO:0043161","GO:0006119","GO:0010257","GO:0140053","GO:0016226"))
plotdata$Description <-factor(plotdata$Description,levels = rev(plotdata$Description))
p2 <- ggplot(data=plotdata, aes(x=Description,y=Count, fill=-log10(qvalue))) + geom_bar(stat="identity", width=0.9) + 
  coord_flip()  + ylab("")+theme(axis.text.y=element_text(color="black", size=12)) +
  scale_fill_gradient2(low = "#43a3ef",mid ="#f9e7a7",high ="#ef767b",midpoint = 6, #连续形legend颜色设置
                        guide = guide_colorbar(frame.colour = "black", #legend边框
                                               ticks.colour = NA))+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text = element_text(size = 12,colour = "black"))
p1/p2
pdf("c51_go.pdf",width = 10,height = 10)
acore <- acore(mfuzz_class,mfuzz_cluster,min.acore = 0.0)
tc <- rep(1:length(mfuzz_cluster$size),times=unlist(lapply(acore,nrow)))
acore <- do.call(rbind,acore);acore$cluster <- tc
acore <- acore %>% group_by(cluster)%>%arrange(desc(MEM.SHIP))%>%mutate(number = row_number())
library(ggrepel)
ggplot(acore[acore$cluster==1,],aes(number,MEM.SHIP))+geom_line()+
  geom_point(data = acore[acore$NAME%in%c("SLC40A1","ATG7","NCOA4","VDAC3","VDAC2","NOX4","SLC11A2","TP53","ACSL4","GSS"),], fill = "orange", size = 3,shape=21)+
  geom_text_repel(aes(label = NAME), 
            data = acore[acore$NAME%in%c("SLC40A1","ATG7","NCOA4","VDAC3","VDAC2","NOX4","SLC11A2","TP53","ACSL4","GSS"),],
             hjust =-0.2, angle = 45) +
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x="Genes ranked by Cluster1 contribution",y = "Gene weight")

ggplot(acore[acore$cluster==5,],aes(number,MEM.SHIP))+geom_line()+
  geom_point(data = acore[acore$NAME%in%c("SLC39A8","SLC39A14","FTH1","CP","PRNP","GCLM","GCLC","ACSL3"),], fill = "orange", size = 3,shape=21)+
  geom_text_repel(aes(label = NAME), 
                  data = acore[acore$NAME%in%c("SLC39A8","SLC39A14","FTH1","CP","PRNP","GCLM","GCLC","ACSL3"),],
                  hjust =-0.2, angle = 45) +
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x="Genes ranked by Cluster5 contribution",y = "Gene weight")
