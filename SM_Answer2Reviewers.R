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
# similarities among regulatory modules, network analysis was performed on SMelated,
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

#/media/victor/
setwd("H:/Semester7(15032017)/RUG/Kai")


goodorder <- read_excel("goodorder.xlsx")
expression_bronch <- read_excel("expression_bronch.xlsx")
expression_nasal <- read_excel("expression_nasal.xlsx")

genes.list <- expression_bronch[,1]
colnames(genes.list)<-c("genes expression file")

#same order than nasal: yes!
head(genes.list) # from bronch
head(expression_nasal[,1])
tail(genes.list)
tail(expression_nasal[,1])

# remove first col (are gene names)
expression_bronch[,1]
expression_nasal[,1]

expression_bronch <- expression_bronch[,2:ncol(expression_bronch)]
expression_nasal <- expression_nasal[,2:ncol(expression_nasal)]
### relevant genes... UGT1A6 with U!

genes<- list(
  "CYP1B1",
  "LOC344887",
  "CYP1A1",
  "SLC7A11",
  "UGT1A6",
  "STATH",
  "BPIFA2",
  "AKR1B10",
  "CRISP3",
  "TIMP3",
  "H19",
  "DNASE1L3",
  "SCEL",
  "TFF1",
  "UGT1A10",
  "FAM177B",
  "ALDH1A3",
  "C1QTNF9B-AS1",
  "TEX26",
  "SEC14L3",
  "MROH9",
  "CTNNAL1",
  "CYP2A13",
  "SLC28A2",
  "SAA1",
  "C3",
  "SAA2"
)
length(genes)
genes.list<-as.matrix(genes.list)
# match indexes of genes (express file) with expression data 
# bronqui
idix2<-match(genes,  as.matrix(genes.list))
c(unlist(genes))
genes.list[idix2]
# nasal
idix3<-match(genes,  as.matrix(genes.list))
c(unlist(genes))
genes.list[idix3]

## bronch and nasal were in same order!!

idix3==idix2

## non smoker choose the 27 genes
expression_bronch2<- expression_bronch[ idix2 , ]
expression_nasal2<- expression_nasal[ idix3 , ]


###################################################
## Standarize ## cols are genes!
###################################################
all.bronch<- scale(  t( expression_bronch2 ))
all.nasal<- scale( t(expression_nasal2))


#sanity check
mean(all.nasal[,2])
sd(all.nasal[,10])


## save
#write.table(all.bronch,"all_bronch.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#write.table(all.nasal,"all_nasal.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#write.table(genes.list[idix3],"genes.txt",sep="\t",row.names=FALSE,col.names=FALSE)
# it can also be genes.list[idix2]

##########################

library("GeneNet")
library("stats4")
library(huge)
# source("http://www.bioconductor.org/biocLite.R")
# class(biocLite)
# biocLite("limma")
library(limma)
source("H:/Semester7(15032017)/RUG/Shrinkage/shrinkagefunctions.R")

setwd("H:/Semester7(15032017)/RUG/Kai/SCientificReports")


set.seed(1)
alpha=0.05

# read data
#expression_bronch <- read.table("H:/Semester7(15032017)/RUG/Kai/shrinkage_S&NS/filter/all_bronch.txt")
genes <- read.table("H:/Semester7(15032017)/RUG/Kai/shrinkage_S&NS/filter/genes.txt")


n=nrow(all.bronch)
p=ncol(all.bronch) #  ==nrow(genes)



## GeneNet 

GGM.SM.B<-ggm.estimate.pcor(all.bronch, method = c("static"))
lambda.SM.B<-attr(GGM.SM.B,'lambda')
r.SM.B<-sm2vec(GGM.SM.B)

        out.SM.B = huge(as.matrix(all.bronch),method="glasso")
        out.select.B = huge.select(out.SM.B, criterion = "stars", rep.num=10)
        rglasso.SM.B=sm2vec( cov2cor(as.matrix(out.select.B$opt.icov)) )
        sum(rglasso.SM.B!=0) 
        
        out2.SM.B = huge(as.matrix(all.bronch),method="tiger")
        out2.select.B = huge.select(out2.SM.B, criterion = "stars", rep.num=10)
        rtiger.SM.B=sm2vec( cov2cor(as.matrix(out2.select.B$opt.icov)) )
        sum(rtiger.SM.B!=0) 

# nasal
GGM.SM.N<-ggm.estimate.pcor( all.nasal, method = c("static"))
lambda.SM.N<-attr(GGM.SM.N,'lambda')
r.SM.N<-sm2vec(GGM.SM.N)
        
        out.SM.N = huge(as.matrix(all.nasal),method="glasso")
        out.select.N = huge.select(out.SM.N, criterion = "stars", rep.num=10)
        rglasso.SM.N=sm2vec( cov2cor(as.matrix(out.select.N$opt.icov)) )
        sum(rglasso.SM.N!=0) 

        
        ##.............. # p values...........
        # bronchial
        genenet.SM.B <- network.test.edges(GGM.SM.B)
        pval.SM.B<-p.shrunk(r.SM.B,ncol(all.bronch),n,lambda.SM.B)
        # nasal
        genenet.SM.N <- network.test.edges(GGM.SM.N)
        pval.SM.N<-p.shrunk(r.SM.N ,ncol(all.nasal),n,lambda.SM.N)
        #................................................
        ##.............. significat BH adjusted p values
        # bronchial
        adj.genenet.SM.B<-p.adjust(   genenet.SM.B$pval , method = "BH", n = length(genenet.SM.B$pval)  )
        adj.shrunk.SM.B<-p.adjust(  pval.SM.B , method = "BH", n = length(pval.SM.B)  )
        
        barplot(  c( sum(adj.genenet.SM.B <=alpha) ,
                     sum(adj.shrunk.SM.B <=alpha), sum(rglasso.SM.B!=0)) ,col=c("red","blue","grey") ,names.arg=c("ENF", " Shrunk MLE","Glasso") )
        
        # nasal
        adj.genenet.SM.N<-p.adjust(   genenet.SM.N$pval , method = "BH", n = length(genenet.SM.N$pval)  )
        adj.shrunk.SM.N<-p.adjust(  pval.SM.N , method = "BH", n = length(pval.SM.N)  )
        
        barplot(  c( sum(adj.genenet.SM.N <=alpha) ,
                     sum(adj.shrunk.SM.N <=alpha), sum(rglasso.SM.N!=0)) ,col=c("red","blue","grey") ,names.arg=c("ENF", " Shrunk MLE","Glasso") )
        
        #write.table(vec2sm(adj.shrunk.SM.B),'BHShrunkSMBronchi.txt',col.names = FALSE,row.names = FALSE)
        #write.table(vec2sm(adj.shrunk.SM.N),'BHShrunkSMNasal.txt',col.names = FALSE,row.names = FALSE)
        #write.table(vec2sm(adj.genenet.SM),'BHgenenet.txt',col.names = FALSE,row.names = FALSE)
        
# Compare networks in each tissue
        # Differences
        length(setdiff(which(rglasso.SM.B!=0),which(rglasso.SM.N!=0)))
        length(setdiff(which(adj.shrunk.SM.B<=alpha),which(adj.shrunk.SM.N<=alpha)))
        # Similarities      
        sum(which(rglasso.SM.B!=0)%in% which(rglasso.SM.N!=0))
        sum(c(which(adj.shrunk.SM.B<=alpha))%in% c(which(adj.shrunk.SM.N<=alpha)))
        
        ab<- cbind(adj.genenet.SM.B <=alpha, adj.shrunk.SM.B<=alpha, rglasso.SM.B!=0) 
        vennDiagram(ab, include = "both", circle.col = c("Blue"),
                    names = c("ENF Br", "Shrunk MLE Br", "Glasso Br"),cex = 1, counts.col = "black")# 
        genes$V1[sm.index(diag(p), diag = FALSE)[which(rowSums(ab)==3),]]
        
        # n
        an<- cbind(adj.genenet.SM.N <=alpha, adj.shrunk.SM.N<=alpha, rglasso.SM.N!=0) 
        vennDiagram(an, include = "both", circle.col = c("red"),
                    names = c("ENF Na", "Shrunk MLE Na", "Glasso Na"),cex = 1, counts.col = "black")# 
        genes$V1[sm.index(diag(p), diag = FALSE)[which(rowSums(an)==3),]]
        
        # Compare methods...............
        # top edges
        #top.glasso.sim<-c(which(rglasso.SM.B!=0)%in% which(rglasso.SM.N!=0))
        
        a<- cbind(adj.shrunk.SM.B <=alpha, 
                  rglasso.SM.B!=0,
                  rglasso.SM.N!=0,
                  adj.shrunk.SM.N<=alpha) 
        vennDiagram(a, include = "both", circle.col = c("Blue","Blue","red","red"),
                    names = c("GGM Br", 
                              "Glasso Br",
                              "Glasso Na", 
                              "GGM Na"),cex = 0.75,lwd=2,counts.col = c("black") ) 
        
        
        genes$V1[unique(sm.index(diag(p), diag = FALSE)[which(rowSums(a)==4),])]
        cbind( genes$V1[ sm.index(diag(p), diag = FALSE)[which(rowSums(a)==4),1]],
                    genes$V1[ sm.index(diag(p), diag = FALSE)[which(rowSums(a)==4),2]])
        
        

        #....................................................
        ################################    
        #3) PLOT GGMs
        
        library("igraph")
        
        
        ## .....SM
        nodes.SM<-data.frame(as.matrix(genes))
        adj.shrunk.SM.B<-vec2sm(adj.shrunk.SM.B); diag(adj.shrunk.SM.B) <-1
        adj.shrunk.SM.N<-vec2sm(adj.shrunk.SM.N); diag(adj.shrunk.SM.N) <-1
        
        temp.SM.B<-adj.shrunk.SM.B ; temp.SM.N<-adj.shrunk.SM.N
        temp.SM.B[temp.SM.B<=alpha]<-1 ;temp.SM.N[temp.SM.N<=alpha]<-1
        temp.SM.B[temp.SM.B<1]<-0 ;temp.SM.N[temp.SM.N<1]<-0
        diag(temp.SM.B)<-0 ;diag(temp.SM.N)<-0
        
        
        ## ..................Merge.............
        # 1 is bronchi, 
        # 2 nasal, 3 both
        temp.SM.merged<-2*temp.SM.N+temp.SM.B
        
        #......................................
        
        
        #................SM which genes overlap 
        overlap<-which(temp.SM.merged == 3, arr.ind=TRUE, useNames = FALSE)
        temp.SM.merged[ overlap]
        unique(c( overlap[,1], overlap[,2]))
        nodes.SM[unique(c( overlap[,1], overlap[,2]))]
        #write.table(x = nodes.SM$V1[unique(c( overlap[,1], overlap[,2]))],file =paste(alpha,"BH_overlap_Bronch&Nasal.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE )
        
        
        
        ##........... SM Connected genes in bronch and nasal
        sum(temp.SM.B)/2
        sum(temp.SM.N)/2
        
        
        
        #..........................................
        # SM Bronchi genes involved
        conn.bronch<-which(temp.SM.B == 1, arr.ind=TRUE, useNames = FALSE)
        temp.SM.B[ conn.bronch]
        unique(c( conn.bronch[,1], conn.bronch[,2]))
        length(nodes.SM[unique(c( conn.bronch[,1], conn.bronch[,2])),1])
        # SM Nasal  genes involved
        conn.Nasal<-which(temp.SM.N == 1, arr.ind=TRUE, useNames = FALSE)
        temp.SM.N[ conn.Nasal]
        unique(c( conn.Nasal[,1], conn.Nasal[,2]))
        length(nodes.SM[unique(c( conn.Nasal[,1], conn.Nasal[,2])),1])
        #..........................................
        
        ############################################################################
        ## remove the unconnected genes
        id.connected<-which( rowSums(temp.SM.merged)!=0)
        temp.SM.merged<-temp.SM.merged[id.connected,] 
        temp.SM.merged<-temp.SM.merged[,id.connected] 
        rowSums(temp.SM.merged)!=0
        names.SM<-nodes.SM[id.connected,1]
        
        #..... connected genes (go for GOenrichment)
        #write.table(x = names.SM,file =paste(alpha,"BH_connected_Bronch&Nasal_SM.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE )
        
        
        ## igraph
        g1.SM <-graph_from_adjacency_matrix(as.matrix(temp.SM.merged), mode = c( "undirected"), weighted = TRUE)
        
        # Nodes
        V(g1.SM)$label <- as.matrix(names.SM)
        
        colrs <- c("red", "blue")
        colEdge<- c("blue", "red", "magenta")
        V(g1.SM)$label.color <- "black"
        
        
        # Edges
        line.thickness <- c(0.5, 0.5,3)
        E(g1.SM)$width <- line.thickness[E(g1.SM)$weight]
        E(g1.SM)$arrow.mode <- 0
        
        
        #layout
        l.SM <-  layout_with_fr(g1.SM) # x, y coordinates (N x 2) for the N nodes in the graph.
        #l <-layout_nicely(g1)
        
        
        #pdf(paste(alpha,"BH_NetworkofBrnch&Nasal.pdf"),width=8,height=8,paper='special')
        
        #plot(g1.SM,vertex.shape="circle", vertex.frame.color="white",vertex.size=3,vertex.label.cex=1 , 
        #     vertex.label.dist=-0.75,edge.color= colEdge[E(g1.SM)$weight],layout=l.SM)
        #legend(x=-1, y=-1, c("Bronchial","Nasal","Overlap"), pch=21,
        #       col="#777777", pt.bg=colEdge, pt.cex=1, cex=1, bty="n", ncol=1)
        
        #dev.off()
        #tiff(filename = paste(alpha,"BH_NetworkofBrnch&Nasal.tiff"),width = 7.5,height = 7.5,res = 300,units = "in")
        win.metafile(filename = paste(alpha,"BH_NetworkofBrnch&Nasal.wmf"),width = 30,height = 30, pointsize = 40)
        V(g1.SM)$label.cex <- 1
        V(g1.SM)$label.cex[which( rowSums(temp.SM.merged==3)>0)] <- 1
        plot(g1.SM,vertex.shape="circle", vertex.frame.color="white",vertex.size=2,
             vertex.label.dist=-0.2,edge.color= colEdge[E(g1.SM)$weight],layout=l.SM)
        legend(x=-0.9, y=-0.9, c("Bronchial","Nasal","Overlap"), pch=21,
               col="#777777", pt.bg=colEdge, pt.cex=1, cex=1, bty="n", ncol=1)
        dev.off()
        
        #########################################################################
        # Scatter
        library(ggplot2)
        ## .....SM
        nodes.SM<-data.frame(nodes.SM)
        colnames(nodes.SM)<-"names"
        #adj.shrunk.SM.B<-vec2sm(adj.shrunk.SM.B); 
        diag(adj.shrunk.SM.B) <-1
        #adj.shrunk.SM.N<-vec2sm(adj.shrunk.SM.N); 
        diag(adj.shrunk.SM.N) <-1
        # fix BUG name "FGF14+A2:A24"  <-"FGF14"
        #nodes.SM[match("FGF14+A2:A24",nodes.SM)]<-"FGF14"
        
        
        
        id<-sm.index(matrix(0,nrow(nodes.SM),nrow(nodes.SM)), diag = FALSE)
        df<-data.frame(unlist(nodes.SM$names[id[,1]]), unlist(nodes.SM$names[id[,2]]))
        edge<-data.frame( paste (df[,1],df[,2], sep = "-", collapse = NULL) )
        
        # Scatter plot 
        # create a data frame ["Edges","-log10 pval Bronchi","-log10 pval Nasal","Type"]
        # Type cathegorical: 
        # 0= edge is non sig both tissues,1= significant in one tissue, 2= sig in both
        scat.pval<-data.frame(edge,-log10(sm2vec(as.matrix(adj.shrunk.SM.B))),-log10(sm2vec(as.matrix(adj.shrunk.SM.N))))
        type<-data.frame(  c(scat.pval[,2]>=-log10(alpha)) + c(scat.pval[,3]>=-log10(alpha)) )
        scat.pval[,4]<-type
        colnames(scat.pval)<-c("Edges","Bronchi","Nasal","Type")
        
        # keep labels only 2= sig in both AND most extreme -log10(alpha) >1.30103  
        # this is alpha 0.05
        scat.pval$Edges[ scat.pval[,4]<2]<-NA
        scat.pval$Edges[ scat.pval[,2]<= -log10(alpha) & scat.pval[,3]<= -log10(alpha)]<-NA
        
        id<-match("TAS2R50-TAS2R31",scat.pval$Edges)
        scat.pval$Edges[id]<-NA
        sp <- ggplot(scat.pval, aes(x=Bronchi, y=Nasal))  +  geom_point(color = "red", size = 1.5)+
          theme_classic(base_size = 14, base_family ="")+ geom_vline(xintercept = -log10(alpha))+ geom_hline( yintercept = -log10(alpha))+
          xlim(min(-log10(adj.shrunk.SM.B), -log10(adj.shrunk.SM.N),na.rm = TRUE), max(-log10(adj.shrunk.SM.B), -log10(adj.shrunk.SM.N),na.rm = TRUE))+
          ylim(min(-log10(adj.shrunk.SM.B), -log10(adj.shrunk.SM.N),na.rm = TRUE), max(-log10(adj.shrunk.SM.B), -log10(adj.shrunk.SM.N),na.rm = TRUE))+
          geom_text(aes(label=scat.pval$Edges),size=4,nudge_y = -0.05,nudge_x = -0.5)+labs(x = "-log10 pval Bronchi",y = "-log 10 pval Nasal")
        
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
        