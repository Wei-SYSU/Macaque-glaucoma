##Enviroment settings
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
#Set path
path = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path)
output_dir=paste(path,"/","All/",sep="")
dir.create(output_dir)
setwd(output_dir)
##load library
library(cowplot)
library(dplyr)
library(reshape2)
library(patchwork)
library(tidyverse)
library(Seurat)
library(harmony)
library(RColorBrewer);library(viridis)
library(Nebulosa)

#Merge data
dir_10X<-paste0(path,"/10X/")
mergedAno<-data.frame(filename10X=list.files(dir_10X))

##Batch analysis of each sample for quality control
mergedlist=list()
for(obname in mergedAno$filename10X){
  merged.data<-Read10X(data.dir = paste0(dir_10X,obname))
  head(merged.data)
  colnames(merged.data)=paste0(colnames(merged.data),"_",obname)
  dim(merged.data)
  merged.data[1:5,1:6]
  merged <- CreateSeuratObject(counts = merged.data, project =obname, min.cells = 3, min.features = 200)
  dim(merged)
  ##Calculate mito gene percentage
  mito.genes<-c("ND1","ND2","COX1","COX2","COX3","ATP6","ATP8","ND3","ND4","ND5","ND6","ND4L","CYTB")
  merged[["percent.mt"]] <- PercentageFeatureSet(merged, features = mito.genes)
  VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 60)
  dim(merged)
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
  all.genes <- rownames(merged)
  merged <- ScaleData(merged, features = all.genes)
  merged <- RunPCA(merged,features = VariableFeatures(object = merged))
  merged <- FindNeighbors(merged, reduction = "pca", dims = 1:20, verbose=FALSE, k.param = 10)
  merged <- FindClusters(merged, verbose = FALSE, resolution = 1)
  merged <- RunUMAP(merged, dims = 1:20)
  merged$orig_seuratcluster<-merged$seurat_clusters
  DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 2)
  VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1,ncol = 1)
  mergedlist[[obname]]=merged
}
##Generate merged file
merged=mergedlist[[1]]
for(i in 2:length(mergedAno$filename10X)){
  merged=merge(merged,mergedlist[[i]])
}

##runharmony and UMAP
if(T){
  merged <- NormalizeData(merged, verbose = FALSE)
  all.genes <- rownames(merged)
  merged <- ScaleData(merged, features = all.genes)
  #find cluster parameter
  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
  merged <- RunPCA(merged, npcs = 80,verbose = FALSE)
  # merged <- JackStraw(merged, verbose = TRUE)
  #merged <- ScoreJackStraw(merged, dims = 1:20)
  #JackStrawPlot(merged)
  #ElbowPlot(merged)
  #using Harmony to remove batch effect
  merged  <- RunHarmony(merged, "orig.ident", do.pcs=FALSE)
  merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20, verbose=FALSE, k.param = 10)
  merged <- FindClusters(merged, verbose = FALSE, resolution = 1)
  merged <- RunUMAP(merged,reduction = "harmony", dims = 1:20) 
  merged$orig_seuratcluster_harmony<-merged$seurat_clusters
  
  #存储
  save_plot(file = paste0(output_dir,"umap.png"),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 2))
  save_plot(file = paste0(output_dir,"orig.ident","_","umap.png"),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,split.by="orig.ident",label.size = 2))
  save_plot(file = paste0(output_dir,"orig.ident","_","umap2.png"),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,group.by="orig.ident",label.size = 2))
  save_plot(file = paste0(output_dir,"QC.png"),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1,ncol = 1),nrow=4,base_width=15)
  saveRDS(merged,file = paste0(path,"/all.Rdata",sep=""))
}

#Predict and manually name the clusters
if(T){
  merged@active.ident=merged$seurat_clusters
  #定义markers
  Rod = c(c("RHO", "NRL", "SAG","GNGT1", "NR2E3"),c("CNGA1"))
  Conelist = c(c("ARR3", "GNAT2", "GNGT2", "OPN1LW", "PDE6H"), c("GUCA1C"))
  HClist = c(c("LHX1", "ONECUT1", "ONECUT2", "CALB1"),c("SLC12A7","LPL"))
  BClist = c(c("CABP5","VSX2", "TRPM1", "VSX1"),c("LRTM1"))
  AClist = c(c("GAD1", "SLC6A9","C1QL2","TFAP2A"),c("SLC32A1" ))
  RGClist = c(c("RBPMS",  "POU4F2",  "THY1", "SLC17A6", "NEFL", "NEFM", "SNCG"),c("IRX2","SYT2"))
  Mullerlist = c(c("RLBP1"),c("FRZB","CRYAB","DAPL1","CP"))
  Microglialist = c(c("C1QA", "C1QB", "C1QC"),c("CXCL8","CFD","HLA-DRA","FCER1G","AIF1"))
  Astrocytelist = c(c("RGCC","ATP1A4","GYPC","GFAP","ANGPTL1"))
  EClist = c(c("RGS5", "MGP", "MYL9","IGFBP7"),c("APOLD1","HIGD1B","GNG11","NOTCH3"))
  Majorclusters<-list(predictname<-"Majorclusters",
                      Majorclusters<-c("PR","HC","BC","AC","RGC","Muller","Astrocyte","IC","EC"),
                      Majorclusterlist<-list(PR = c(c("RHO", "NRL", "SAG","GNGT1", "NR2E3"),c("CNGA1"),c("Arr3", "Pde6h","Opn1sw", "Opn1mw")),
                                             HC = c(c("LHX1", "ONECUT1", "ONECUT2", "CALB1"),c("SLC12A7","LPL")),
                                             BC = c(c("CABP5","VSX2", "TRPM1", "VSX1"),c("LRTM1")),
                                             AC = c(c("GAD1", "SLC6A9","C1QL2","TFAP2A"),c("SLC32A1" )),
                                             RGC =  c(c("RBPMS",  "POU4F2",  "THY1", "SLC17A6", "NEFL", "NEFM", "SNCG"),c("IRX2","SYT2")),
                                             Muller = c(c("RLBP1"),c("FRZB","CRYAB","DAPL1","CP")),
                                             Astrocyte = c(c("RGCC","ATP1A4","GYPC","GFAP","ANGPTL1")),
                                             IC=c(c("FCERAG","PTPRC"),c("C1QA", "C1QB", "C1QC"),c("CXCL8","CFD","HLA-DRA","FCER1G","AIF1")),
                                             EC = c(c("RGS5", "MGP", "MYL9","IGFBP7"),c("APOLD1","HIGD1B","GNG11","NOTCH3"))))
  save_plot(file = paste(output_dir,"markers.png",sep=""),
            DotPlot(merged, features = Majorclusters[[3]], group.by="seurat_clusters",dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  
  temp<-merged
  predicted=Predictcluster(temp,topredict=Majorclusters);rm(temp)
  count=table(merged@meta.data[,c("seurat_clusters","orig.ident")])
  write.csv(count,paste(output_dir,"count.csv",sep=""))
  count<-read.csv(file=paste(output_dir,"count.csv",sep=""))
  predicted=merge(predicted,count,by.x="seurat_clusters",by.y="X")
  predicted$seurat_clusters=as.numeric(predicted$seurat_clusters)
  predicted=predicted[order(predicted$seurat_clusters),]
  write.csv(predicted,paste(output_dir,"predicted.csv",sep=""))
  #用predicted1命名
  plotby=paste("predict1","seurat_clusters",sep="_")
  newidents<-predicted[["predict1"]]
  merged@active.ident<-merged$seurat_clusters
  names(newidents) <- levels(merged)
  merged <- RenameIdents(merged, newidents)
  #merged[[plotby]]<-merged@active.ident
  merged[[plotby]]<-paste(merged@active.ident,merged$seurat_clusters,sep="_")
  #Plotting and manual examination
  markertoplot = Majorclusters[[3]]
  save_plot(file = paste(output_dir,"umap_",plotby,"1.png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by=plotby)+ggtitle(obname),base_asp =1.5,base_height = 6)
  save_plot(file = paste(output_dir,"umap_",plotby,".png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by=plotby,split.by="orig.ident")+ggtitle(obname),base_asp =2,base_height = 6)
  save_plot(file = paste(output_dir,"umap_",plotby,"2.png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by="orig.ident")+ggtitle(obname),base_asp =1.5,base_height = 6)
  save_plot(file = paste(output_dir,"dot_",plotby,".png",sep=""),
            DotPlot(merged, features = markertoplot, group.by=plotby,dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  save_plot(file = paste(output_dir,"QC_",plotby,".png",sep=""),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by=plotby,pt.size=0,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,"QC_",plotby,"_2.png",sep=""),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    group.by=plotby,split.by="orig.ident",
                    pt.size=0,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,"dot_",plotby,"_IC.png",sep=""),
            DotPlot(merged, features = IC,  group.by=plotby,dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  
  ##name by manually adjusted clusters
  predicted=read.csv(paste(output_dir,"predicted.csv",sep=""))
  plotby=paste("predict3","seurat_clusters",sep="_")
  newidents<-predicted[["predict3"]]
  merged@active.ident<-merged$seurat_clusters
  names(newidents) <- levels(merged)
  merged <- RenameIdents(merged, newidents)
  #merged[[plotby]]<-merged@active.ident
  merged[[plotby]]<-paste(merged@active.ident,merged$seurat_clusters,sep="_")
  
  #画图
  markertoplot = Majorclusters[[3]]
  save_plot(file = paste(output_dir,"umap_",plotby,"1.png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by=plotby)+ggtitle(obname),base_asp =1.5,base_height = 6)
  save_plot(file = paste(output_dir,"umap_",plotby,".png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by=plotby,split.by="orig.ident")+ggtitle(obname),base_asp =2,base_height = 6)
  save_plot(file = paste(output_dir,"umap_",plotby,"2.png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by="orig.ident")+ggtitle(obname),base_asp =1.5,base_height = 6)
  save_plot(file = paste(output_dir,"dot_",plotby,".png",sep=""),
            DotPlot(merged, features = markertoplot, group.by=plotby,dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  save_plot(file = paste(output_dir,"QC_",plotby,".png",sep=""),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by=plotby,pt.size=0,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,"QC_",plotby,"_2.png",sep=""),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    group.by=plotby,split.by="orig.ident",
                    pt.size=0,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,"dot_",plotby,"_IC.png",sep=""),
            DotPlot(merged, features = IC,  group.by=plotby,dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  
  ##Add metadata
  plotby="Majorcluster"
  predicted=read.csv(paste(output_dir,"predicted.csv",sep=""))
  newidents<-predicted[["predict3"]]
  merged@active.ident<-merged$seurat_clusters
  names(newidents) <- levels(merged)
  merged <- RenameIdents(merged, newidents)
  #merged[[plotby]]<-merged@active.ident
  merged[[plotby]]<-merged@active.ident
  #画图
  markertoplot = Majorclusters[[3]]
  save_plot(file = paste(output_dir,"umap_",plotby,"1.png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by=plotby)+ggtitle(obname),base_asp =1.5,base_height = 6)
  save_plot(file = paste(output_dir,"umap_",plotby,".png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by=plotby,split.by="orig.ident")+ggtitle(obname),base_asp =2,base_height = 6)
  save_plot(file = paste(output_dir,"umap_",plotby,"2.png",sep=""),
            DimPlot(merged, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 4,cols=colors,repel=T,group.by="orig.ident")+ggtitle(obname),base_asp =1.5,base_height = 6)
  save_plot(file = paste(output_dir,"dot_",plotby,".png",sep=""),
            DotPlot(merged, features = markertoplot, group.by=plotby,dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  save_plot(file = paste(output_dir,"QC_",plotby,".png",sep=""),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by=plotby,pt.size=0,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,"QC_",plotby,"_2.png",sep=""),
            VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    group.by=plotby,split.by="orig.ident",
                    pt.size=0,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,"dot_",plotby,"_IC.png",sep=""),
            DotPlot(merged, features = IC, group.by=plotby,dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
}

saveRDS(merged,file = paste0(path,"/all.rds",sep=""))

##Second round of QC
sub=c("PR","HC","BC","AC","RGC","Muller","Astrocyte","EC","IC")

for(i in 1:length(sub)){
  merged@active.ident=merged$Majorcluster
  pbmc=subset(merged,idents=sub[[i]])
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  #find cluster parameter
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- RunPCA(pbmc, npcs = 80,verbose = FALSE)
  pbmc  <- RunHarmony(pbmc, "orig.ident", do.pcs=FALSE)
  pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:80, verbose=FALSE, k.param = 10)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 1)
  pbmc <- RunUMAP(pbmc, reduction = "harmony",dims = 1:80)
  #saveplot
  save_plot(file = paste0(output_dir,i,"umap.png"),
            DimPlot(pbmc, reduction = "umap",label = TRUE, pt.size = 0.2,label.size = 2))
  save_plot(file = paste0(output_dir,i,"orig.ident","_","umap.png"),
            DimPlot(pbmc, reduction = "umap",label = TRUE, pt.size = 0.2,split.by="orig.ident",label.size = 2))
  save_plot(file = paste0(output_dir,i,"orig.ident","_","umap2.png"),
            DimPlot(pbmc, reduction = "umap",label = TRUE, pt.size = 0.2,group.by="orig.ident",label.size = 2))
  save_plot(file = paste0(output_dir,i,"QC.png"),
            VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1,ncol = 1),nrow=4,base_width=15)
  save_plot(file = paste(output_dir,i,"_dot_",".png",sep=""),
            DotPlot(merged, features = markertoplot, cols = c("grey","blue"), dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  pbmc$tempcluster=pbmc$Majorcluster
  metalist[[i]]=pbmc@meta.data
}
merged$tempcluster=merged$Majorcluster
for(i in 1:length(sub)){
  ##remove the clusters expressing 
  clustertoremove=c()#according to the marker expression
  metalist[[i]][metalist[[i]]$seurat_clusters%in%clustertoremove, "tempcluster"]="removed"
  merged@meta.data[row.names(metalist[[i]]),"Majorcluster"]= metalist[[i]][row.names(metalist[[i]]),"tempcluster"]
}
##Annotation of IC
sub=c("IC")
for(i in 1:length(sub)){
  merged@active.ident=merged$Majorcluster
  pbmc=subset(merged,idents=sub[[i]])
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  #find cluster parameter
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- RunPCA(pbmc, npcs = 20,verbose = FALSE)
  pbmc  <- RunHarmony(pbmc, "orig.ident", do.pcs=FALSE)
  pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20, verbose=FALSE, k.param = 10)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.2)
  pbmc <- RunUMAP(pbmc, reduction = "harmony",dims = 1:20)
  save_plot(file = paste(output_dir,i,"dot_","2.png",sep=""),
            DotPlot(merged, features = c(c("FCERAG","PTPRC"),c("C1QA", "C1QB", "C1QC"),c("CXCL8","CFD","HLA-DRA","FCER1G","AIF1")),
                    dot.scale = 12) + RotatedAxis(), base_height = 15,base_asp = 1.5,bg="white")
  ##根据single R注释免疫细胞类群
  library(SingleR) 
  sce=pmbc
  sce_for_SingleR <- GetAssayData(sce, slot="data")
  clusters=sce@meta.data$seurat_clusters
  hsHemato <- HumanPrimaryCellAtlasData()
  pred.Immu <- SingleR(test = sce_for_SingleR, ref = hsHemato, labels = hsHemato$label.main,
                            method = "cluster", clusters = clusters, 
                            assay.type.test = "logcounts", assay.type.ref = "logcounts")
  cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                      hsHemato=pred.Immu$labels )
  rm(sce,sce_for_SingleR)
  ##Check the cellType manually and annotate 
  cellType$final=cellType$hsHemato; edit(cellType)
  ##add meta
  pbmc@active.ident=pbmc$seurat_clusters
  table(pbmc$seurat_clusters)  
  for(i in 1:nrow(celltype)){
    pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == celltype[i,"tempcluster"]),'final'] <- celltype[i,"final"]}
  merged@meta.data[row.names(pbmc@meta.data),"Majorcluster"]= pbmc@meta.data[row.names(pbmc@meta.data),"tempcluster"]
}

saveRDS(merged,file = paste0(path,"/all.Rdata",sep=""))
