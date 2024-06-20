##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

#设置路径
path = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path)
output_dir=paste(path,"/","All230816/",sep="")
dir.create(output_dir)
setwd(output_dir)

library(Seurat)
library(data.table)

load(file=paste0(path,"/all.Rdata"))

###----DEATH------------
library(msigdbr)
tmt=msigdbr_collections()
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") 
m_t2g2 <- msigdbr(species = "Homo sapiens", category = "H") 
m_t2g3 <- msigdbr(species = "Homo sapiens", category = "C5")
m_t2g=rbind(m_t2g,m_t2g2,m_t2g3)
m_t2g$gs_name=tolower(m_t2g$gs_name);m_t2g$gs_name=gsub("_"," ",m_t2g$gs_name);m_t2g$gs_name=capitalize(m_t2g$gs_name)
death=list()
for(i in c("death","apoptosis","ferroptosis","pyroptosis","necro")){
  death[[i]]=grep(i,m_t2g$gs_name,value=T)%>%unique
}
death[["PANoptosis"]]=death[2:5]%>%unlist()%>%unique

deathterm=list()
for(i in c("death","apoptosis","ferroptosis","pyroptosis","necro","PANoptosis")){
deathterm[[i]]=m_t2g[m_t2g$gs_name%in%death[[i]],]
deathterm[[i]]$gs_name=i
geneselect1=names(table(deathterm[[i]]$gene_symbol))[table(deathterm[[i]]$gene_symbol)>=floor(length(death[[i]])/2)]
geneselect2=names(table(deathterm[[i]]$gene_symbol))[table(deathterm[[i]]$gene_symbol)>3]
geneselect=c(geneselect1,geneselect2)
deathterm[[i]]=deathterm[[i]][deathterm[[i]]$gene_symbol%in%geneselect,]
deathterm[[i]]=distinct(deathterm[[i]],gene_symbol,.keep_all=T)
}

for (i in unlist(death)%>%unique()){
  deathterm[[i]]=m_t2g[m_t2g$gs_name%in%i,]
}
k.list=list()
for(i in names(deathterm)){
  k.list[[i]]=deathterm[[i]][,'gene_symbol']%>%unlist()%>%unique()
}

deathlist=data.frame('gs_name'=array(),'gene_symbol'=array());deathlist=deathlist[-1,]
for(i in names(deathterm)){
  tmt=deathterm[[i]][,c('gs_name','gene_symbol')]
  deathlist=rbind(deathlist,tmt)
}

if(T){
  library(UCell)
  library(irGSEA)
  pbmc3k.final=merged
  pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA",
                               slot = "data", seeds = 123, ncores = 1,
                               min.cells = 3, min.feature = 0,
                               custom = T, geneset = k.list, 
                               geneid = "symbol",
                               method = c("AUCell"),
                               aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                               kcdf = 'Gaussian') 
  Seurat::Assays(pbmc3k.final)
}

irGSEA_death_signature=pbmc3k.final@assays[c("AUCell")]
save(irGSEA_death_signature,deathterm,deathlist,file="E:/Weiyun/Rresults/final/irGSEA_death_signature.rdata")


#load(file="irGSEA_Disease2_PNAS2023.rdata")
##---------Density map----------------------
method="AUCell"
for(geneset in names(k.list)){
  #geneset= names(k.list)[1]
  scatterplot <- irGSEA.density.scatterplot(object = pbmc3k.final,
                                            method =method,
                                            adjust = 3,
                                            size =0.3,
                                            show.geneset = geneset,
                                            pal = "plasma",
                                            reduction = "tsne")+  theme(panel.border = element_rect(color = "black",fill = NA,size = 0.5)) +
    labs(x="tSNE_1", y ="tSNE_2")+ ggtitle(paste0(ob,": ",geneset))+theme(panel.border = element_blank(),
                                                                          #panel.border = element_blank(),
                                                                          #axis.title = element_blank(),
                                                                          axis.line = element_blank(),
                                                                          axis.ticks = element_blank(),
                                                                          axis.text.x  = element_blank(),
                                                                          axis.text.y = element_blank(),
                                                                          axis.title = element_blank(),
                                                                          plot.title =element_text(size=10,face="bold") )
  
  scatterplot  
  #ggsave(paste(output_dir,"Density_",markername,method,ob,geneset,".pdf",sep=""),width = 3,height =2,bg = "white",limitsize = FALSE)
  ggsave(paste(output_dir,"Density_",markername,method,ob,geneset,".png",sep=""),width = 3,height =2,bg = "white",limitsize = FALSE)
}

