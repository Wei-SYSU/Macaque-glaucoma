##Enviroment settings
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
#Set path
path = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path)
output_dir=paste(path,"/","GO/",sep="")
dir.create(output_dir)
setwd(output_dir)

load(file=paste0(path,"/DEGlistnew.Rdata"))
names(DEGlistnew)
##prepare enviroment
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(stringr)
library(DOSE)
library(RColorBrewer);library(viridis)
#Prepare DEGs for enrichment analaysis
order1=c("RGC","Microglia","Muller","PR","HC","BC","AC","EC","Astrocyte","Immune")
setname="major_Group"
DE=DEGlistnew[[1]];DE[DE$avg_diff>0,"UorD"]="up";DE[DE$avg_diff<0,"UorD"]="down";DE$UorD=factor(DE$UorD,levels = c("up","down"))
DE=DE[DE$majortype%in%order1,]
DE$majortype=factor(DE$majortype,levels=order1)
DE=DE[order(DE$majortype,DE$UorD,-abs(DE$avg_diff)),]


FC=0.25;symbols_list=DE[DE$p_val_adj<0.05&abs(DE$avg_diff)>FC,]; setname=paste("major_Group",FC);pro=paste(setname, "all")
gene_names <- symbols_list$gene
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#ENSMFAG=unique(grep("ENSMFAG",symbols_list$gene,value=T))
symbols_list=symbols_list[symbols_list$gene%in%entrez_genes$SYMBOL,]
gcSample=symbols_list[,c("gene","majortype")]
gcSample = split(gcSample$gene,gcSample$majortype)
#gcSample=symbols_list[,c("gene","majortype","time")]

#--------WKH------------------
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:WIKIPATHWAYS") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g2 <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g=rbind(m_t2g,m_t2g2)
m_t2g$gs_name=tolower(m_t2g$gs_name);m_t2g$gs_name=gsub("_"," ",m_t2g$gs_name);m_t2g$gs_name=capitalize(m_t2g$gs_name)
#write.csv(m_t2g[m_t2g$gs_name=="KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",],"KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.csv")
WKH <- compareCluster(gcSample, 
                      fun="enricher",
                      TERM2GENE=m_t2g,
                      pvalueCutoff = 0.05)
dotplot(WKH,showCategory =50,includeAll = T,size ="count",by="count", 
        font.size = 7,label_format =30)  + 
  scale_radius(rang=c(1,3))+
  scale_color_viridis(  begin = 0.9,
                        end = 0.1,)+
  #scale_color_gradientn(colours = rev(c('#336699','yellow','#D62728FF')))+ #颜色+
  scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))

#-------- Run full GO enrichment test for BP ---------------------------
formula_res_BP <- compareCluster(
  gcSample,
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont     = "BP",
  keyType = "SYMBOL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  pvalueCutoff  = 0.05)
# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego_BP <- simplify(
  formula_res_BP, 
  cutoff=0.8, 
  by="p.adjust", 
  select_fun=min
) 
library(enrichplot)
library(GOSemSim)
library(DOSE)
library(grid)
formula_res_BP_paired <- pairwise_termsim(lineage1_ego_BP)
#stage1=grep("early",names(formula_res_BP_paired@geneClusters),value=T)
#early=subset(formula_res_BP_paired,geneClusters%in%stage1)#[formula_res_BP_paired@geneClusters%in%stage1,]
if(T){
  cn=50
  
  pdf(paste(pro,'GO_BP_cluster_simplified_emapplot_group',cn,'2.pdf') ,width = 15,height = 15)
  print(cowplot::plot_grid(p1,ncol=1))
  dev.off() 
  emapplot(formula_res_BP_paired, showCategory = cn,
           #node_label="category",
           node_label="group",
           hilight.params = list(category =c("RGC","Muller","Microglia"), alpha_hilight = 1, alpha_no_hilight = 0.7),
           label_format = 40,repel=T,
           group_category=T,group_legend=T,cex_label_group=2,nWords=5,nCluster=10,
           clusterFunction=cluster::fanny,
           with_edge = T,
           shadowtext = F,
           legend_n = 5,
           min_edge=0.2,
           cex_line = 0.5, cex_pie2axis=6,
           #pie="count", 
           layout = "nicely"
           #layout =   'gem'
  ) +
    scale_fill_manual(breaks=order1,
                      values=celltype[order1,"color1"])
  pdf(paste(pro,'GO_BP_cluster_simplified_emapplot_category',cn,'.pdf') ,width = 15,height = 15)
  print(cowplot::plot_grid(p2,ncol=1))
  dev.off() 
  #    print(cowplot::plot_grid(p1, p2,p3,p4, ncol=2, labels=LETTERS[1:4]))
  #cowplot::plot_grid(p1, p2,p3,p4, ncol=2, labels=LETTERS[1:4])  
  cnetplot(formula_res_BP_paired,showCategory = cn, node_label="gene",
           foldChange = 0.5,
           shadowtext="none",
           colorEdge=T,
           legend_n = 5,
           min_edge=0.2,
           cex_line = 0.5, 
           cex_label_gene = 0.8) +
    scale_fill_manual(breaks=order1,
                      values=celltype[order1,"color1"])
}
