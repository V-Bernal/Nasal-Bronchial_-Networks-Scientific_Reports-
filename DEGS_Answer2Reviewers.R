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
# similarities among regulatory modules, network analysis was performed on DEGSelated,
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
#  The standard probability density is presented in [Fisher,R.A. (1924) The distribution of the partial DEGSelation coefficient. Metron, 3,329-332.] 
#  The simulation of data, as well as estimation with the optimal shrinkage (lambda) is done with
#  [Sch?fer,J. and Strimmer,K. (2005a),  GeneNet 1.2.13. CRAN.]
#
#************************************************************************;


#..............................................
# 1) Take full data from bronch and nasal
# takes the genes from DEGS file
# i) orders by adjusted pval DEGS, 
# ii) keep the relevant genes only (DEGS) and Standarize

rm(list=ls())
library(readxl)
#library( readr)

setwd("/media/victor/VICTOR/Semester7(15032017)/RUG/Kai")


expression_bronch <- read_excel("expression_bronch.xlsx")
expression_nasal <- read_excel("expression_nasal.xlsx")
ncol(expression_bronch) # first col are probesets

## separate the first column (gene names)
genes.list <- expression_bronch[,1]
data.frame( head(genes.list),head(expression_nasal[,1])  )
data.frame( tail(genes.list)  ,  tail(expression_nasal[,1])  )


expression_bronch <- expression_bronch[,2:ncol(expression_bronch)]
expression_nasal <- expression_nasal[,2:ncol(expression_nasal)]


### Focus on DEGSELATED and DE genes
DEGS_genes <- read.delim2("difexpres.txt")
DEGS_genes<-as.character(DEGS_genes[,1])

####### indexes of genes with expression data 
genes.list<-as.matrix(genes.list)

# fix BUG name "FGF14+A2:A24"  <-"FGF14"
DEGS_genes<-as.character(DEGS_genes)
DEGS_genes[match("FGF14+A2:A24",DEGS_genes)]<-c("FGF14")

# bronqui
DEGS_idix2<-match(as.matrix(DEGS_genes),  as.matrix(genes.list))# gene list comes from expression
data.frame( DEGS_genes  , genes.list[DEGS_idix2] )

# nasal
DEGS_idix3<-match(as.matrix(DEGS_genes),  as.matrix(genes.list))
data.frame( DEGS_genes,   genes.list[DEGS_idix3])

#log2
#library("preprocessCore")
#data2<-normalize.quantiles( data2,copy=TRUE) #row is a probe.
#cite Bolstad et al, Bioinformatics (2003)


#.............. Standarize
# bronqui. here we remove the first column of gene names
bronch.z.DEGS<-scale (t( expression_bronch[DEGS_idix2,] ))
#bronch.z<-t(bronch.z)

# nasal
nasal.z.DEGS<-scale (t( expression_nasal[DEGS_idix3, ] ))
#nasal.zN<-t(nasals.zN)
#..............................................



##............... COLUMNS are genes!!!!!
# save bronchial
#write.table(bronch.z.DEGS,"shrinkage_DEGSelatedgenes/Bronch/bronqui.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#write.table(nasal.z.DEGS,"shrinkage_DEGSelatedgenes/Nasal/nasal.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#write.table(c(unlist(genes)),"shrinkage_DEGSelatedgenes/genes.txt",sep="\t",row.names=FALSE,col.names=FALSE)
# save nasal



#sanity check
mean(bronch.z.DEGS[,1])
sd(nasal.z.DEGS[,1])
mean(nasal.z.DEGS[,1])
sd(nasal.z.DEGS[,1])



#################################################################
### 2)

library(huge)
library("GeneNet")
library(stats4)
#source("http://www.bioconductor.org/biocLite.R")
#class(biocLite)
#biocLite("limma")
library(limma)
source("/media/victor/VICTOR/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")
setwd("/media/victor/VICTOR/Semester7(15032017)/RUG/Kai/SCientificReports")

set.seed(1)

# read data
#expression_bronch <- read.table("bronqui.txt")
genes <- DEGS_genes#read.table("H:/Semester7(15032017)/RUG/Kai/shrinkage_DEGSelatedgenes/genes.txt")
#DEGSelated_genes <- read_excel("H:/Semester7(15032017)/RUG/Kai/DEGSelated genes.xlsx")

n<-nrow(bronch.z.DEGS)
p=ncol(bronch.z.DEGS)


alpha=0.01

##.............. GGM Networks...........
# bronchial
GGM.DEGS.B<-ggm.estimate.pcor( bronch.z.DEGS, method = c("static"))
lambda.DEGS.B<-attr(GGM.DEGS.B,'lambda')
r.DEGS.B<-sm2vec(GGM.DEGS.B)

              out.DEGS.B = huge(as.matrix(bronch.z.DEGS),method="glasso")
              out.select.B = huge.select(out.DEGS.B, criterion = "stars", rep.num=10)
              #plot(out.select.B)
              rglasso.DEGS.B=sm2vec( cov2cor(as.matrix(out.select.B$opt.icov)) )
              sum(rglasso.DEGS.B!=0) 

# nasal
GGM.DEGS.N<-ggm.estimate.pcor( nasal.z.DEGS, method = c("static"))
lambda.DEGS.N<-attr(GGM.DEGS.N,'lambda')
r.DEGS.N<-sm2vec(GGM.DEGS.N)

              out.DEGS.N = huge(as.matrix(nasal.z.DEGS),method="glasso")
              out.select.N = huge.select(out.DEGS.N, criterion = "stars", rep.num=10 )
              #plot(out.select.B)
              rglasso.DEGS.N=sm2vec( cov2cor(as.matrix(out.select.N$opt.icov)) )
              sum(rglasso.DEGS.N!=0)  

#..........................................

##.............. # p values...........
# bronchial
genenet.DEGS.B <- network.test.edges(GGM.DEGS.B)
pval.DEGS.B<-p.shrunk(r.DEGS.B,ncol(bronch.z.DEGS),n,lambda.DEGS.B)
# nasal
genenet.DEGS.N <- network.test.edges(GGM.DEGS.N)
pval.DEGS.N<-p.shrunk(r.DEGS.N ,ncol(nasal.z.DEGS),n,lambda.DEGS.N)
#................................................

#save.image("degs.R")


##.............. significat BH adjusted p values
# bronchial
adj.genenet.DEGS.B<-p.adjust(   genenet.DEGS.B$pval , method = "BH", n = length(genenet.DEGS.B$pval)  )
adj.shrunk.DEGS.B<-p.adjust(  pval.DEGS.B , method = "BH", n = length(pval.DEGS.B)  )

barplot(  c( sum(adj.genenet.DEGS.B <=alpha) ,
             sum(adj.shrunk.DEGS.B <=alpha), sum(rglasso.DEGS.B!=0)) ,col=c("red","blue","grey") ,names.arg=c("ENF", " Shrunk MLE","Glasso") )

# nasal
adj.genenet.DEGS.N<-p.adjust(   genenet.DEGS.N$pval , method = "BH", n = length(genenet.DEGS.N$pval)  )
adj.shrunk.DEGS.N<-p.adjust(  pval.DEGS.N , method = "BH", n = length(pval.DEGS.N)  )

barplot(  c( sum(adj.genenet.DEGS.N <=alpha) ,
             sum(adj.shrunk.DEGS.N <=alpha), sum(rglasso.DEGS.N!=0)) ,col=c("red","blue","grey") ,names.arg=c("ENF", " Shrunk MLE","Glasso") )

#write.table(vec2sm(adj.shrunk.DEGS.B),'BHShrunkDEGSBronchi.txt',col.names = FALSE,row.names = FALSE)
#write.table(vec2sm(adj.shrunk.DEGS.N),'BHShrunkDEGSNasal.txt',col.names = FALSE,row.names = FALSE)
#write.table(vec2sm(adj.genenet.DEGS),'BHgenenet.txt',col.names = FALSE,row.names = FALSE)

              # Compare networks in each tissue
              # Differences
              length(setdiff(which(rglasso.DEGS.B!=0),which(rglasso.DEGS.N!=0)))
              length(setdiff(which(adj.shrunk.DEGS.B<=alpha),which(adj.shrunk.DEGS.N<=alpha)))
              # Similarities      
              sum(which(rglasso.DEGS.B!=0)%in% which(rglasso.DEGS.N!=0))
              sum(c(which(adj.shrunk.DEGS.B<=alpha))%in% c(which(adj.shrunk.DEGS.N<=alpha)))
              
              ab<- cbind(adj.genenet.DEGS.B <=alpha, adj.shrunk.DEGS.B<=alpha, rglasso.DEGS.B!=0) 
              vennDiagram(ab, include = "both", circle.col = c("Blue"),
                          names = c("ENF Br", "Shrunk MLE Br", "Glasso Br"),cex = 1, counts.col = "black")# 
              genes[sm.index(diag(p), diag = FALSE)[which(rowSums(ab)==3),]]
              
              # n
              an<- cbind(adj.genenet.DEGS.N <=alpha, adj.shrunk.DEGS.N<=alpha, rglasso.DEGS.N!=0) 
              vennDiagram(an, include = "both", circle.col = c("red"),
                          names = c("ENF Na", "Shrunk MLE Na", "Glasso Na"),cex = 1, counts.col = "black")# 
              genes[sm.index(diag(p), diag = FALSE)[which(rowSums(an)==3),]]
              
              # Compare methods...............
              # top edges
              #top.glasso.sim<-c(which(rglasso.SM.B!=0)%in% which(rglasso.SM.N!=0))
              
              a<- cbind(adj.shrunk.DEGS.B <=alpha, 
                        rglasso.DEGS.B!=0,
                        rglasso.DEGS.N!=0,
                        adj.shrunk.DEGS.N<=alpha) 
              vennDiagram(a, include = "both", circle.col = c("Blue","Blue","red","red"),
                          names = c("GGM Br", 
                                    "Glasso Br",
                                    "Glasso Na", 
                                    "GGM Na"),cex = .75,lwd=2, counts.col = c("black")) 
              
              # connected in glasso
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.B!=0),]] ) 
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.N!=0),]] )
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.B!=0) %in%  
                                                                  which(rglasso.corr.N!=0),]] )
              # connected in all
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rowSums(a)==4),]] )
              genes[unique(sm.index(diag(p), diag = FALSE)[which(rowSums(a)==0),])]
              cbind( genes[ sm.index(diag(p), diag = FALSE)[which(rowSums(a)==0),1]],
                     genes[ sm.index(diag(p), diag = FALSE)[which(rowSums(a)==0),2]])


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# We got empty networks for optimal lambda. We now check for another lambda
              
             #plot(out.select.B2)
              rglasso.DEGS.B2 =sm2vec( cov2cor(as.matrix(out.DEGS.B$icov[[2]])) )
              sum(rglasso.DEGS.B2!=0) 
              
              #plot(out.select.B2)
              rglasso.DEGS.N2 =sm2vec( cov2cor(as.matrix(out.DEGS.N$icov[[2]])) )
              sum(rglasso.DEGS.N2!=0)  
              
              a22<- cbind(adj.shrunk.DEGS.B <=alpha, 
                        rglasso.DEGS.B2!=0,
                        rglasso.DEGS.N2!=0,
                        adj.shrunk.DEGS.N<=alpha) 
              vennDiagram(a22, include = "both", circle.col = c("Blue","Blue","red","red"),
                          names = c("GGM Br", 
                                    "Glasso Br",
                                    "Glasso Na", 
                                    "GGM Na"),cex = .75,lwd=2, counts.col = c("black")) 
              
              # connected in glasso
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.B2!=0),]] ) 
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.N2!=0),]] )
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rglasso.corr.B2!=0) %in%  
                                                                  which(rglasso.corr.N2!=0),]] )
              # connected in all
              unique(genes$...1[sm.index(diag(p), diag = FALSE)[which(rowSums(a22)==4),]] )
              genes[unique(sm.index(diag(p), diag = FALSE)[which(rowSums(a22)==0),])]
              cbind( genes[ sm.index(diag(p), diag = FALSE)[which(rowSums(a22)==0),1]],
                     genes[ sm.index(diag(p), diag = FALSE)[which(rowSums(a22)==0),2]])
              
              
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#....................................................


################################    
#3) PLOT GGMs

library("igraph")


## .....DEGS
nodes.DEGS<-data.frame(DEGS_genes)
adj.shrunk.DEGS.B<-vec2sm(adj.shrunk.DEGS.B); diag(adj.shrunk.DEGS.B) <-1
adj.shrunk.DEGS.N<-vec2sm(adj.shrunk.DEGS.N); diag(adj.shrunk.DEGS.N) <-1

temp.DEGS.B<-adj.shrunk.DEGS.B ; temp.DEGS.N<-adj.shrunk.DEGS.N
temp.DEGS.B[temp.DEGS.B<=alpha]<-1 ;temp.DEGS.N[temp.DEGS.N<=alpha]<-1
temp.DEGS.B[temp.DEGS.B<1]<-0 ;temp.DEGS.N[temp.DEGS.N<1]<-0
diag(temp.DEGS.B)<-0 ;diag(temp.DEGS.N)<-0


## ..................Merge.............
# 1 is bronchi, 
# 2 nasal, 3 both
temp.DEGS.merged<-2*temp.DEGS.N+temp.DEGS.B

#......................................


#................DEGS which genes overlap 
overlap<-which(temp.DEGS.merged == 3, arr.ind=TRUE, useNames = FALSE)
temp.DEGS.merged[ overlap]
unique(c( overlap[,1], overlap[,2]))
nodes.DEGS$...1[unique(c( overlap[,1], overlap[,2]))]
#write.table(x = nodes.DEGS$V1[unique(c( overlap[,1], overlap[,2]))],file =paste(alpha,"BH_overlap_Bronch&Nasal.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE )



##........... DEGS Connected genes in bronch and nasal
sum(temp.DEGS.B)/2
sum(temp.DEGS.N)/2



#..........................................
# DEGS Bronchi genes involved
conn.bronch<-which(temp.DEGS.B == 1, arr.ind=TRUE, useNames = FALSE)
temp.DEGS.B[ conn.bronch]
unique(c( conn.bronch[,1], conn.bronch[,2]))
length(nodes.DEGS$...1[unique(c( conn.bronch[,1], conn.bronch[,2]))])
# DEGS Nasal  genes involved
conn.Nasal<-which(temp.DEGS.N == 1, arr.ind=TRUE, useNames = FALSE)
temp.DEGS.N[ conn.Nasal]
unique(c( conn.Nasal[,1], conn.Nasal[,2]))
length(nodes.DEGS$...1[unique(c( conn.Nasal[,1], conn.Nasal[,2]))])
#..........................................

############################################################################
## remove the unconnected genes
id.connected<-which( rowSums(temp.DEGS.merged)!=0)
temp.DEGS.merged<-temp.DEGS.merged[id.connected,] 
temp.DEGS.merged<-temp.DEGS.merged[,id.connected] 
rowSums(temp.DEGS.merged)!=0
names.DEGS<-nodes.DEGS[id.connected,]

#..... connected genes (go for GOenrichment)
#write.table(x = names.DEGS,file =paste(alpha,"BH_connected_Bronch&Nasal_DEGS.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE )


## igraph
g1.DEGS <-graph_from_adjacency_matrix(as.matrix(temp.DEGS.merged), mode = c( "undirected"), weighted = TRUE)

# Nodes
V(g1.DEGS)$label <- as.matrix(names.DEGS)

colrs <- c("red", "blue")
colEdge<- c("blue", "red", "magenta")
V(g1.DEGS)$label.color <- "black"


# Edges
line.thickness <- c(0.2, 0.2,1.5)
E(g1.DEGS)$width <- line.thickness[E(g1.DEGS)$weight]
E(g1.DEGS)$arrow.mode <- 0


#layout
l.DEGS <-  layout_with_fr(g1.DEGS) # x, y coordinates (N x 2) for the N nodes in the graph.
#l <-layout_nicely(g1)


#pdf(paste(alpha,"BH_NetworkofBrnch&Nasal.pdf"),width=8,height=8,paper='special')

plot(g1.DEGS,vertex.shape="circle", vertex.frame.color="white",vertex.size=2,vertex.label.cex=0.5 , 
     vertex.label.dist=-0.25,edge.color= colEdge[E(g1.DEGS)$weight],layout=l.DEGS)
legend(x=-1, y=-1, c("Bronchial","Nasal","Overlap"), pch=21,
       col="#777777", pt.bg=colEdge, pt.cex=0.8, cex=.8, bty="n", ncol=1)

#dev.off()
#tiff(filename = paste(alpha,"BH_NetworkofBrnch&Nasal.tiff"),width = 7.5,height = 7.5,res = 300,units = "in")
#win.metafile(filename = paste(alpha,"BH_NetworkofBrnch&Nasal.wmf"),width = 30,height = 30, pointsize = 40)


#########################################################################
# Scatter
library(ggplot2)
## .....DEGS
nodes.DEGS<-data.frame(DEGS_genes)[,1]
#adj.shrunk.corr.B<-vec2sm(adj.shrunk.corr.B); 
diag(adj.shrunk.DEGS.B) <-1
#adj.shrunk.corr.N<-vec2sm(adj.shrunk.corr.N); 
diag(adj.shrunk.DEGS.N) <-1
# fix BUG name "FGF14+A2:A24"  <-"FGF14"
nodes.DEGS[match("FGF14+A2:A24",nodes.DEGS)]<-"FGF14"



id<-sm.index(matrix(0,length(nodes.DEGS),length(nodes.DEGS)), diag = FALSE)
df<-data.frame(nodes.DEGS[id[,1]],nodes.DEGS[id[,2] ])
edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )

# Scatter plot 
# create a data frame ["Edges","-log10 pval Bronchi","-log10 pval Nasal","Type"]
# Type cathegorical: 
# 0= edge is non sig both tissues,1= significant in one tissue, 2= sig in both
scat.pval<-data.frame(edge,-log10(sm2vec(as.matrix(adj.shrunk.DEGS.B))),-log10(sm2vec(as.matrix(adj.shrunk.DEGS.N))))
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
  theme_classic(base_size = 14, base_family ="")+ geom_vline(xintercept = -log10(alpha), alpha=0.75)+ geom_hline( yintercept = -log10(alpha), alpha=0.75)+
  xlim(min(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE), max(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE))+
  ylim(min(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE), max(-log10(adj.shrunk.corr.B), -log10(adj.shrunk.corr.N),na.rm = TRUE))+
  geom_text(aes(label=scat.pval$Edges),size=4,nudge_y = 0.25,nudge_x = -2)+labs(x = "-log10 pval Bronchi",y = "-log 10 pval Nasal")

#tiff(paste(alpha,"scatterCO.tiff"),width = 7.5,height = 7.5,res = 300,units = "in")
sp 
#dev.off()


#win.metafile(paste(alpha,"scatterCO.wmf"), width =7.5,height = 7.5, pointsize = 12)

# ggplot(scat.pval, aes(x=Bronchi, y=Nasal))  +  geom_point(color = "red", size = 0.8)+
#   theme_classic(base_size = 12, base_family ="")+ geom_vline(xintercept = -log10(alpha),color = "gray")+ geom_hline( yintercept = -log10(alpha),color = "gray")+
#   xlim(min(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE), max(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE))+
#   ylim(min(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE), max(-log10(BHshrunk), -log10(BHshrunkNASAL),na.rm = TRUE))+
#   geom_text(aes(label=scat.pval$Edges),size=3,nudge_y = 0.25,nudge_x = -2)+labs(x = "-log10 pval Bronchi",y = "-log 10 pval Nasal")

#dev.off()
