library(countreg)
library(foreach)
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
source("scripts/KS_test.R")

###### load count matrix
data <- readRDS("data/counts.rds")

#### load groups
groups <- readRDS("data/groups.rds")
name.group <- as.character(unique(groups))

### load TF genes
load(file="data/ALL_TF.rda")
TF.gene <- ALL_TF
###################
##### consider neuronal lineage
name.groupcond <- name.group[-which(name.group=="Sus")]

for (k  in 1: length(name.groupcond)) {
    ind.group <- which(groups==name.groupcond[k])
    length(ind.group)
    Y <- t(data[,ind.group])

    ##### filter TF genes!
    ind.TF <- which(colnames(Y) %in% TF.gene)
    Y.TF <- Y[,ind.TF]

    #### filter some genes with low mean
    meanY.TF  <- apply(Y.TF,2,mean)
    skip.TF <- which(meanY.TF>= 0.005)
    X.genes <- Y.TF [,skip.TF]

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

    if (dim(qX)[1]>1000){
        mean.cells <- apply(qX,1,mean)
        or <- order(mean.cells,decreasing=TRUE)
        fX <- qX[or[1:1000],]
    }else fX <- qX

    ######## Step 2: transform X to X^\alpha with alpha obtained from KS test
    ####finding optimal alpha by KS test
    alpha.ks <- seq(0.01,.5,length=100)
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
    name.file <- paste(k,"dataTFtranscelsalpha0595.RData", sep = "-")
    save(datatest,file=name.file)
}

###################  Estimating gene networks by zinb1
library(learn2count)

for (k in 1:5) {
    name.file <- paste(k,"dataTFtranscelsalpha0595.RData", sep = "-")
    load(name.file)
    n <- dim(datatest)[1]
    p <- dim(datatest)[2]
    datatest <- as.matrix(datatest,n,p)

    maxcard <- 3
    alpha <-  2*pnorm(n^.15,lower.tail=F)
    #####zinb1
    zinbPC1.time <- system.time(adj.zinb1 <- zinb1.noT(datatest,maxcard,alpha, extend=TRUE))
    name.fileresult <- paste(k,"dataTFtranscelsalpha059515result.RData", sep = "-")
    save(adj.zinb1,datatest, file = name.fileresult )
}
