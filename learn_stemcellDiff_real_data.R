############################# Analysis for selected  genes in the stem Cell Differentiation geneset
library(countreg)
library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
source("scripts/KS_test.R")

###### load count matrix
data <- readRDS("data/counts.rds")
data.all <- t(data)

#### load groups
groups <- readRDS("data/groups.rds")
name.group <- as.character(unique(groups))

### load stem Cell Differentiation genes
name.genes <- read.delim( "data/stemCellDifferentiation.txt",header = FALSE)
name.genes <- unique( name.genes[-1,2])


##### consider only cells in HBC* group
  ind.group <- which(groups=="HBC*")
  Y <- t(data[,ind.group])
  
  ##### filter selected genes!
  ind.genes <- which(colnames(Y) %in% name.genes)
  Y.genes <- Y[,ind.genes]
  
  #### filter some genes with low mean
  meanY.genes  <- apply(Y.genes,2,mean)
  skip.genes <- which(meanY.genes>= 0.005)
  X.genes <- Y.genes [,skip.genes]
  
  
  ######## Preprocessing data
  ####### Step 1: normalizing with 95%quantile matching
  quanNorm <- .95
  p <- ncol(X.genes)
  dis <- apply(X.genes,1,quantile,quanNorm)
  toskip <- which(dis==0)
  if(length(toskip)>0){
    qnum <- mean(dis[-toskip])
    todoX <- X.genes[-toskip,]
    qX <- todoX/(dis[-toskip]%o%rep(1,p)/qnum)  
  }else{
    qnum <- mean(dis)
    todoX <- X.genes
    qX <- todoX/(dis%o%rep(1,p)/qnum) 
  }
  
  qX <- X.genes
  if (dim(qX)[1]>1000){
    var.cells <- apply(qX,1,var)
    or <- order(var.cells,decreasing=TRUE)
    fX <- qX[or[1:1000],]
  }else fX <- qX
  
  
  ######## Step 2: transform X to X^\alpha with alpha obtained from KS test
  ####finding optimal alpha by KS test
  alpha.ks <- seq(0.01,1,length=100)
  KSstat_est <- foreach(i = 1:dim(fX)[2], .combine = "cbind",.packages=c("iZID")) %dopar%{
    ks.stat <- rep(NA,length(alpha.ks)) 
    for (j in 1:length(alpha.ks)) {
      x <- fX[,i]^alpha.ks[j]
      p.nb <- 1-mean(x)/var(x)
      r <- mean(x)*(1-p.nb)/p.nb
      ks.stat[j] <- KS.statistic(x,r,p.nb)
    }
    ks.stat
  }
  sum.ks <- apply(KSstat_est,1, sum)
  ind.alpha <- which(sum.ks==min(sum.ks))
  
  #### transform X to X^\alpha with alpha obtained from KS test
  datatest <- floor(fX^alpha.ks[ind.alpha])
  gene.skip <- which(apply(datatest,2,mean)==0)
  if (length(gene.skip)>0) datatest <- datatest [,-gene.skip]
  

###################  Estimating gene networks by zinb1
  library(learn2count)
  
  n <- dim(datatest)[1]
  p <- dim(datatest)[2]
  datatest <- as.matrix(datatest,n,p)
  
  maxcard <- 3
  alpha <-  2*pnorm(n^.15,lower.tail=F)
  #####zinb1
  zinbPC1.time <- system.time(adj.zinb1 <- zinb1.noT(datatest,maxcard,alpha, extend=TRUE))
  save(adj.zinb1, zinbPC1.time,datatest, file = "3-stemCellDifferentiationtransform0595315result.RData" )
  
########### plot the resulting graph

library(igraph)

load("data/ALL_TF.rda")
TF.gene <- ALL_TF
#load("data/3-stemCellDifferentiationtransform0595315result.RData")

colnames(adj.zinb1) <- rownames(adj.zinb1) <- colnames(datatest)
node.levels <- apply(adj.zinb1, 2, sum)
ind.hub <- which(node.levels>=9)
print(colnames(adj.zinb1)[ind.hub])
ind.hubTF <- which(colnames(adj.zinb1)[ind.hub] %in% TF.gene)
print(colnames(adj.zinb1)[ind.hub][ind.hubTF])


skip <- which(apply(adj.zinb1,2,sum)==0)
graph <- graph.adjacency(adj.zinb1[-skip,-skip], mode=("undirected"))
graph$color[1:length(colnames(adj.zinb1))]<- "light blue"
graph$color[ind.hub]<-"orange"
graph$color[ind.hub][ind.hubTF] <- "red"
graph$color <- graph$color[-skip]
pdf("HBCactstemtrans.pdf",width=8,height=7,paper='special') 
plot(graph,layout=layout_with_fr,edge.arrow.size=0.05,
     vertex.size=5,vertex.size2 = 4,
     vertex.label.cex = 0.22,
     asp = 0.5,
     margin = 0,
     vertex.label.dist=0.05,
     vertex.frame.color=graph$color,
     vertex.color=graph$color,
     vertex.label.color ="black")
dev.off()


######## plot subgraph of Trp63

ind_Trp63 <- which(colnames(datatest)=="Trp63")
ind_subgraph <- which(adj.zinb1[,ind_Trp63]==1)
adj_sub <- adj.zinb1[c(ind_Trp63,ind_subgraph),c(ind_Trp63,ind_subgraph)]
graph <- graph.adjacency(adj_sub, mode=("undirected"))
graph$color[1:length(colnames(adj_sub))]<- "light blue"
graph$color[1]<-"red"

pdf("HBCactTrp63trans.pdf",width=13,height=9,paper='special') 
plot(graph,layout=layout_with_fr,edge.arrow.size=0.05,
     vertex.size=8,vertex.size2 = 4,
     vertex.label.cex = 0.8,
     asp = 0.5,
     margin = 0.01,
     vertex.label.dist=0.02,
     vertex.frame.color=graph$color,
     vertex.color=graph$color,
     vertex.label.color ="black")
dev.off()
