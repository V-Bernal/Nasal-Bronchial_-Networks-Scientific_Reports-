#|**********************************************************************;
#  This script reproduces the results in 
#  * Project: Gene network approach reveals co-expression patterns in nasal and bronchial epithelium (Imspkamp and Bernal et al).
#  * Author            : Kai Imkamp*, Victor Bernal*, 
#                        Marco Grzegorzcyk, Peter Horvatovich, Cornelis J Vermeulen, Irene H Heijink,
#                        Victor Guryev, Huib AM Kerstjens,Maarten van den Berge???, Alen Faiz???
#  * Date created      : 06 JUN 2019
#***********************************
# Description
# The aim is to explore to what extent nasal brushings can be used as a 
# proxy for bronchial brushings in the context of gene expression profiling.
# We performed genome wide gene expression profiling on matched nasal and bronchial
# epithelial brushes from 77 respiratory healthy individuals. To investigate differences and
# similarities among regulatory modules, network analysis was performed on correlated,
# differentially expressed and smoking-related genes using Gaussian Graphical Models.
#***************************************
# Prerequisite R library (CRAN)
# requires R libraries "GeneNet" , "stats4", "ggplot2", and "Hmisc".
#************************************************************************
# Parameters       
# * p = Number of variables (e.g. genes)   
# * n = Number of samples  
# * alpha = significance level
# 
#************************************************************************
#  * Revision History  : 26 JUN 2019
#  **********************************************************************
#  The shrunk probability density is presented here in [] 
#  For computation of empirical p-value via Monte Carlo see  [Martinez, W. L., & Martinez, A. R. (2007). Computational statistics handbook with MATLAB. Chapman and Hall/CRC.] 
#  The standard probability density is presented in [Fisher,R.A. (1924) The distribution of the partial correlation coefficient. Metron, 3,329-332.] 
#  The simulation of data, as well as estimation with the optimal shrinkage (lambda) is done with
#  [Sch?fer,J. and Strimmer,K. (2005a),  GeneNet 1.2.13. CRAN.]
#
#************************************************************************;


#..............................................
# 1) Take full data from bronch and nasal
# takes the genes from CORR file
# i) orders by adjusted pval CORR, 
# ii) keep the relevant genes only (CORR) and Standarize

rm(list=ls())
library(readxl)
#library( readr)

# /media/victor/VICTOR
setwd("H:/Semester7(15032017)/RUG/Kai")


expression_bronch <- read_excel("expression_bronch.xlsx")
expression_nasal <- read_excel("expression_nasal.xlsx")
ncol(expression_bronch) # first col are probesets

## separate the first column (gene names)
genes.list <- expression_bronch[,1]
data.frame( head(genes.list),head(expression_nasal[,1])  )
data.frame( tail(genes.list)  ,  tail(expression_nasal[,1])  )


expression_bronch <- expression_bronch[,2:ncol(expression_bronch)]
expression_nasal <- expression_nasal[,2:ncol(expression_nasal)]


### Focus on CORRELATED and DE genes
correlated_genes <- read_excel("correlated genes.xlsx")
CORR_genes<-correlated_genes[,1]

####### indexes of genes with expression data 
genes.list<-as.matrix(genes.list)

# bronqui
CORR_idix2<-match(as.matrix(CORR_genes),  as.matrix(genes.list))# gene list comes from expression
data.frame( c(unlist(CORR_genes))  , genes.list[CORR_idix2] )

# nasal
CORR_idix3<-match(as.matrix(CORR_genes),  as.matrix(genes.list))
data.frame( c(unlist(CORR_genes)),   genes.list[CORR_idix3])

#log2
#library("preprocessCore")
#data2<-normalize.quantiles( data2,copy=TRUE) #row is a probe.
#cite Bolstad et al, Bioinformatics (2003)


#.............. Standarize
# bronqui. here we remove the first column of gene names
bronch.z.corr<-scale (t( expression_bronch[CORR_idix2,] ))
#bronch.z<-t(bronch.z)

# nasal
nasal.z.corr<-scale (t( expression_nasal[CORR_idix3, ] ))
#nasal.zN<-t(nasals.zN)
#..............................................



##............... COLUMNS are genes!!!!!
# save bronchial
#write.table(bronch.z.corr,"shrinkage_correlatedgenes/Bronch/bronqui.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#write.table(nasal.z.corr,"shrinkage_correlatedgenes/Nasal/nasal.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#write.table(c(unlist(genes)),"shrinkage_correlatedgenes/genes.txt",sep="\t",row.names=FALSE,col.names=FALSE)
# save nasal



#sanity check
mean(bronch.z.corr[,1])
sd(nasal.z.corr[,1])
mean(nasal.z.corr[,1])
sd(nasal.z.corr[,1])



#################################################################
### 2)

library(huge)
library("GeneNet")
library(stats4)
#source("http://www.bioconductor.org/biocLite.R")
#class(biocLite)
#biocLite("limma")
library(limma)
#/media/victor/VICTOR/
source("H:/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")

set.seed(1)
# read data
#expression_bronch <- read.table("bronqui.txt")
genes <- CORR_genes #read.table("G:/Semester7(15032017)/RUG/Kai/shrinkage_correlatedgenes/genes.txt")
#correlated_genes <- read_excel("H:/Semester7(15032017)/RUG/Kai/correlated genes.xlsx")

n<-nrow(bronch.z.corr)
p=ncol(bronch.z.corr)

alpha=0.01

##.............. GGM Networks...........
# bronchial
GGM.corr.B<-ggm.estimate.pcor( bronch.z.corr, method = c("static"))
lambda.corr.B<-attr(GGM.corr.B,'lambda')
r.corr.B<-sm2vec(GGM.corr.B)

          out.corr.B = huge(as.matrix(bronch.z.corr),method="glasso")
          out.select.B = huge.select(out.corr.B, criterion = "stars", rep.num=10)
          #plot(out.select.B)
          rglasso.corr.B=sm2vec( cov2cor(as.matrix(out.select.B$opt.icov)) )
          sum(rglasso.corr.B!=0) 
          
# nasal
GGM.corr.N<-ggm.estimate.pcor( nasal.z.corr, method = c("static"))
lambda.corr.N<-attr(GGM.corr.N,'lambda')
r.corr.N<-sm2vec(GGM.corr.N)

          out.corr.N = huge(as.matrix(nasal.z.corr),method="glasso")
          out.select.N = huge.select(out.corr.N, criterion = "stars", rep.num=10 )
          #plot(out.select.B)
          rglasso.corr.N=sm2vec( cov2cor(as.matrix(out.select.N$opt.icov)) )
          sum(rglasso.corr.N!=0)  
     

#..........................................

##.............. # p values...........
# bronchial
genenet.corr.B <- network.test.edges(GGM.corr.B)
pval.corr.B<-p.shrunk(r.corr.B,ncol(bronch.z.corr),n,lambda.corr.B)
# nasal
genenet.corr.N <- network.test.edges(GGM.corr.N)
pval.corr.N<-p.shrunk(r.corr.N ,ncol(nasal.z.corr),n,lambda.corr.N)
#................................................




##.............. significat BH adjusted p values
# bronchial
adj.genenet.corr.B<-p.adjust(   genenet.corr.B$pval , method = "BH", n = length(genenet.corr.B$pval)  )
adj.shrunk.corr.B<-p.adjust(  pval.corr.B , method = "BH", n = length(pval.corr.B)  )

barplot(  c( sum(adj.genenet.corr.B <=alpha) ,
             sum(adj.shrunk.corr.B <=alpha), sum(rglasso.corr.B!=0)) ,col=c("red","blue","grey") ,names.arg=c("ENF", " Shrunk MLE","Glasso") )

# nasal
adj.genenet.corr.N<-p.adjust(   genenet.corr.N$pval , method = "BH", n = length(genenet.corr.N$pval)  )
adj.shrunk.corr.N<-p.adjust(  pval.corr.N , method = "BH", n = length(pval.corr.N)  )

barplot(  c( sum(adj.genenet.corr.N <=alpha) ,
             sum(adj.shrunk.corr.N <=alpha), sum(rglasso.corr.N!=0)) ,col=c("red","blue","grey") ,names.arg=c("ENF", " Shrunk MLE","Glasso") )

#write.table(vec2sm(adj.shrunk.corr.B),'BHShrunkCORRBronchi.txt',col.names = FALSE,row.names = FALSE)
#write.table(vec2sm(adj.shrunk.corr.N),'BHShrunkCORRNasal.txt',col.names = FALSE,row.names = FALSE)
#write.table(vec2sm(adj.genenet.corr),'BHgenenet.txt',col.names = FALSE,row.names = FALSE)

                  # Compare networks in each tissue
                  # Differences
                  length(setdiff(which(rglasso.corr.B!=0),which(rglasso.corr.N!=0)))
                  length(setdiff(which(adj.shrunk.corr.B<=alpha),which(adj.shrunk.corr.N<=alpha)))
                  # Similarities      
                  sum(which(rglasso.corr.B!=0)%in% which(rglasso.corr.N!=0))
                  sum(c(which(adj.shrunk.corr.B<=alpha))%in% c(which(adj.shrunk.corr.N<=alpha)))
                  
                  ab<- cbind(adj.genenet.corr.B <=alpha, adj.shrunk.corr.B<=alpha, rglasso.corr.B!=0) 
                  vennDiagram(ab, include = "both", circle.col = c("Blue"),
                              names = c("ENF Br", "Shrunk MLE Br", "Glasso Br"),cex = 1, counts.col = "black")# 
                  genes[sm.index(diag(p), diag = FALSE)[which(rowSums(ab)==3),]]
                  
                  # n
                  an<- cbind(adj.genenet.corr.N <=alpha, adj.shrunk.corr.N<=alpha, rglasso.corr.N!=0) 
                  vennDiagram(an, include = "both", circle.col = c("red"),
                              names = c("ENF Na", "Shrunk MLE Na", "Glasso Na"),cex = 1, counts.col = "black")# 
                  genes[sm.index(diag(p), diag = FALSE)[which(rowSums(an)==3),]  ,1]
                  
                  # Compare methods...............
                  # top edges
                  #top.glasso.sim<-c(which(rglasso.SM.B!=0)%in% which(rglasso.SM.N!=0))
                  
                  a<- cbind(adj.shrunk.corr.B <=alpha, 
                            rglasso.corr.B!=0,
                            rglasso.corr.N!=0,
                            adj.shrunk.corr.N<=alpha) 
                  
                  vennDiagram(a, include = "both", circle.col = c("Blue","Blue","red","red"),
                              names = c("GGM Br", 
                                        "Glasso Br",
                                        "Glasso Na", 
                                        "GGM Na"),cex = 0.75,lwd=2, counts.col = c("black"))  
                  
                  # connected in glasso
                  unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.B!=0),]] ) 
                  unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.N!=0),]] )
                  unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.B!=0) %in%  
                   which(rglasso.corr.N!=0),]] )
                  # connected in all
                  unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rowSums(a)==4),]] )


#....................................................


################################    
#3) PLOT GGMs

library("igraph")


## .....CORR
nodes.corr<-data.frame(CORR_genes)
adj.shrunk.corr.B<-vec2sm(adj.shrunk.corr.B); diag(adj.shrunk.corr.B) <-1
adj.shrunk.corr.N<-vec2sm(adj.shrunk.corr.N); diag(adj.shrunk.corr.N) <-1

temp.CORR.B<-adj.shrunk.corr.B ; temp.CORR.N<-adj.shrunk.corr.N
temp.CORR.B[temp.CORR.B<=alpha]<-1 ;temp.CORR.N[temp.CORR.N<=alpha]<-1
temp.CORR.B[temp.CORR.B<1]<-0 ;temp.CORR.N[temp.CORR.N<1]<-0
diag(temp.CORR.B)<-0 ;diag(temp.CORR.N)<-0


## ..................Merge.............
# 1 is bronchi, 
# 2 nasal, 3 both
temp.CORR.merged<-2*temp.CORR.N+temp.CORR.B

#......................................


#................CORR which genes overlap 
overlap<-which(temp.CORR.merged == 3, arr.ind=TRUE, useNames = FALSE)
temp.CORR.merged[ overlap]
unique(c( overlap[,1], overlap[,2]))
nodes.corr$...1[unique(c( overlap[,1], overlap[,2]))]
#write.table(x = nodes.corr$V1[unique(c( overlap[,1], overlap[,2]))],file =paste(alpha,"BH_overlap_Bronch&Nasal.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE )



##........... CORR Connected genes in bronch and nasal
sum(temp.CORR.B)/2
sum(temp.CORR.N)/2



#..........................................
# CORR Bronchi genes involved
conn.bronch<-which(temp.CORR.B == 1, arr.ind=TRUE, useNames = FALSE)
temp.CORR.B[ conn.bronch]
unique(c( conn.bronch[,1], conn.bronch[,2]))
length(nodes.corr$...1[unique(c( conn.bronch[,1], conn.bronch[,2]))])
# CORR Nasal  genes involved
conn.Nasal<-which(temp.CORR.N == 1, arr.ind=TRUE, useNames = FALSE)
temp.CORR.N[ conn.Nasal]
unique(c( conn.Nasal[,1], conn.Nasal[,2]))
length(nodes.corr$...1[unique(c( conn.Nasal[,1], conn.Nasal[,2]))])
#..........................................

############################################################################
## remove the unconnected genes
id.connected<-which( rowSums(temp.CORR.merged)!=0)
temp.CORR.merged<-temp.CORR.merged[id.connected,] 
temp.CORR.merged<-temp.CORR.merged[,id.connected] 
rowSums(temp.CORR.merged)!=0
names.corr<-nodes.corr[id.connected,]

#..... connected genes (go for GOenrichment)
#write.table(x = names.corr,file =paste(alpha,"BH_connected_Bronch&Nasal_corr.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE )


## igraph
g1.CORR <-graph_from_adjacency_matrix(as.matrix(temp.CORR.merged), mode = c( "undirected"), weighted = TRUE)

# Nodes
V(g1.CORR)$label <- as.matrix(names.corr)

colrs <- c("red", "blue")
colEdge<- c("blue", "red", "magenta")
V(g1.CORR)$label.color <- "darkgrey"
V(g1.CORR)$label.color[which( rowSums(temp.CORR.merged==3)>0)] <- "black"


# Edges
line.thickness <- c(0.2, 0.2, 3)
E(g1.CORR)$width <- line.thickness[E(g1.CORR)$weight]
E(g1.CORR)$arrow.mode <- 0


#layout
l.corr <-  layout_with_fr(g1.CORR) # x, y coordinates (N x 2) for the N nodes in the graph.
#l.corr <-layout.circle(g1.CORR)#layout_nicely(g1.CORR)

V(g1.CORR)$label.cex <- 0.5
V(g1.CORR)$label.cex[which( rowSums(temp.CORR.merged==3)>0)] <- 0.65

win.metafile(filename = paste(alpha,"BH_NetworkofBrnch&Nasal4.wmf"),width = 30,height = 30, pointsize = 40)
#pdf(paste(alpha,"BH_NetworkofBrnch&Nasal2.pdf"),width=8,height=8,paper='special')
#vertex.label.cex=0.5, 
plot(g1.CORR,vertex.shape="circle", vertex.frame.color="white",vertex.size=2,
     vertex.label.dist=-0.2,edge.color= colEdge[E(g1.CORR)$weight],layout=l.corr)
legend(x=-0.8, y=-0.8, c("Bronchial","Nasal","Overlap"), pch=21,
       col="#777777", pt.bg=colEdge, pt.cex=1, cex=1, bty="n", ncol=1)

dev.off()
#tiff(filename = paste(alpha,"BH_NetworkofBrnch&Nasal.tiff"),width = 7.5,height = 7.5,res = 300,units = "in")
#win.metafile(filename = paste(alpha,"BH_NetworkofBrnch&Nasal.wmf"),width = 30,height = 30, pointsize = 40)


#########################################################################
# Scatter
library(ggplot2)
## .....CORR
nodes.corr<-data.frame(CORR_genes)[,1]
#adj.shrunk.corr.B<-vec2sm(adj.shrunk.corr.B); 
diag(adj.shrunk.corr.B) <-1
#adj.shrunk.corr.N<-vec2sm(adj.shrunk.corr.N); 
diag(adj.shrunk.corr.N) <-1
# fix BUG name "FGF14+A2:A24"  <-"FGF14"
nodes.corr[match("FGF14+A2:A24",nodes.corr)]<-"FGF14"



id<-sm.index(matrix(0,length(nodes.corr),length(nodes.corr)), diag = FALSE)
df<-data.frame(nodes.corr[id[,1]],nodes.corr[id[,2] ])
edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )

# Scatter plot 
# create a data frame ["Edges","-log10 pval Bronchi","-log10 pval Nasal","Type"]
# Type cathegorical: 
# 0= edge is non sig both tissues,1= significant in one tissue, 2= sig in both
scat.pval<-data.frame(edge,-log10(sm2vec(as.matrix(adj.shrunk.corr.B))),-log10(sm2vec(as.matrix(adj.shrunk.corr.N))))
type<-data.frame(  c(scat.pval[,2]>=-log10(alpha)) + c(scat.pval[,3]>=-log10(alpha)) )
scat.pval[,4]<-type
colnames(scat.pval)<-c("Edges","Bronchi","Nasal","Type")

# keep labels only 2= sig in both AND most extreme -log10(alpha) >10  
# this is alpha 10E-11
scat.pval$Edges[ scat.pval[,4]<2]<-NA
scat.pval$Edges[ scat.pval[,2]<=12 & scat.pval[,3]<=12]<-NA

id<-match("TAS2R50-TAS2R31",scat.pval$Edges)
scat.pval$Edges[id]<-NA
sp <- ggplot(scat.pval, aes(x=Bronchi, y=Nasal))  +  geom_point(color = "red", size = 1.5)+
  theme_classic(base_size = 14, base_family ="")+ geom_vline(xintercept = -log10(alpha))+ geom_hline( yintercept = -log10(alpha))+
  xlim(min(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE), max(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE))+
  ylim(min(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE), max(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE))+
  geom_text(aes(label=scat.pval$Edges),size=4,nudge_y = 0.25,nudge_x = -2)+labs(x = "-log10 pval Bronchi",y = "-log 10 pval Nasal")

#tiff(paste(alpha,"scatterCO.tiff"),width = 7.5,height = 7.5,res = 300,units = "in")
sp 
#dev.off()


win.metafile(paste(alpha,"scatterCO.wmf"), width =7.5,height = 7.5, pointsize = 12)
sp
# ggplot(scat.pval, aes(x=Bronchi, y=Nasal))  +  geom_point(color = "red", size = 0.8)+
#   theme_classic(base_size = 12, base_family ="")+ geom_vline(xintercept = -log10(alpha),color = "gray")+ geom_hline( yintercept = -log10(alpha),color = "gray")+
#   xlim(min(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE), max(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE))+
#   ylim(min(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE), max(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE))+
#   geom_text(aes(label=scat.pval$Edges),size=3,nudge_y = 0.25,nudge_x = -2)+labs(x = "-log10 pval Bronchi",y = "-log 10 pval Nasal")

dev.off()

