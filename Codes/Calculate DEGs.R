##Enviroment settings
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
#Set path
path = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path)
output_dir=paste(path,"/","DEG/",sep="")
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

load(file=paste0(path,"/all.Rdata"))

DEGlistnew=list()
#---------------------all_Group-------------------
plotby="major";merged$major=factor(merged$major,levels=c("RGC","Muller","Microglia","BC","AC","PR","HC","EC","Astrocyte","Immune"));merged@active.ident=as.factor(merged$major)
table(merged$major)
table(merged)
Temp_markers=NULL
steps=NULL
for(k in names(table(merged$major))){
  sublist=subset(merged,idents=k)
  all.genes=rownames(sublist)
  sublist=ScaleData(sublist, features =all.genes)
  if(T){
    sublist$Group_major=paste(sublist$Group,sublist$major,sep="_")
    sublist@active.ident=as.factor(sublist$Group_major)
    type=paste("glau_",k," VS ","ctl_", k,sep="")
    print(paste("calculate DE", type))
    temp2=FindMarkers(sublist, 
                      ident.1 = paste("glau", k,sep="_"),
                      ident.2 = paste("ctl", k,sep="_"),
                      group.by="Group_major", 
                      assay="RNA",
                      test.use="MAST",
                      slot="scale.data",
                      random.seed=0,
                      logfc.threshold = 0, 
                      min.pct = 0.1)
    
    temp2$type=type
    temp2$gene=rownames(temp2)
    temp2$majortype=k
    temp2$minortype=k
    temp2$time="glau"
    temp2$lable="Group"
    #write.csv(temp2, file = paste(output_dir,"DE_",type,".csv",sep="")) 
    Temp_markers=rbind(temp2,Temp_markers)
  }
}
DEGlistnew[["major_Group"]]=Temp_markers

#---------------------all_Stage-------------------
plotby="major";merged$major=factor(merged$major,levels=c("RGC","Muller","Microglia","BC","AC","PR","HC","EC","Astrocyte","Immune"));merged@active.ident=as.factor(merged$major)
table(merged$major)
table(merged@active.ident)
Temp_markers=NULL
steps=NULL
for(k in names(table(merged$major))){
  sublist=subset(merged,idents=k)
  all.genes=rownames(sublist)
  sublist=ScaleData(sublist, features =all.genes)
  if(T){
    sublist$major_Stage=paste(sublist$major,sublist$Stage)
    sublist@active.ident=as.factor(sublist$major_Stage)
    for(j in names(table(sublist$Stage))[2:3]){ 
      type=paste(k, j,"VS","ctl")
      print(paste("calculate DE", type))
      temp2=FindMarkers(sublist, 
                        ident.1 = paste(k, j),
                        ident.2 = paste(k, "ctl"),
                        group.by="major_Stage", 
                        assay="RNA",
                        test.use="MAST",
                        slot="scale.data",
                        random.seed=0,
                        logfc.threshold = 0, 
                        min.pct = 0.1)
      
      temp2$type=type
      temp2$gene=rownames(temp2)
      temp2$majortype=k
      temp2$minortype=k
      temp2$time=j
      temp2$lable="Stage"
      #write.csv(temp2, file = paste(output_dir,"DE_",type,".csv",sep="")) 
      Temp_markers=rbind(temp2,Temp_markers)
    }
  }
}

DEGlistnew[["major_Stage"]]=Temp_markers
save(DEGlistnew,file=paste0(path,"/DEGlistnew.Rdata"))

###Count significant DEGs
DE=DEGlistnew[["major_Group"]]
DE[DE$avg_diff>0,"UorD"]="Up";DE[DE$avg_diff<0,"UorD"]="Down"
FC=0.25;genesig=DE[DE$p_val_adj<0.05&abs(DE$avg_diff)>FC,]
table(genesig[,c("type","UorD")])

DE=DEGlistnew[["major_Stage"]]
DE[DE$avg_diff>0,"UorD"]="Up";DE[DE$avg_diff<0,"UorD"]="Down"
FC=0.5;genesig=DE[DE$p_val_adj<0.05&abs(DE$avg_diff)>FC,]
table(genesig[,c("type","UorD")])

