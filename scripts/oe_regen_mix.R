library(countreg)
library(foreach)
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

######## KS tests
KS.statistic <- function(x,r,p){
    mle_ori <- nb.zihmle(x, r, p, type = "zi", lowerbound=1e-04,
                         upperbound=10^6)
    probs_ori <- mle_ori[3] + (1 - mle_ori[3]) * stats::pnbinom(0:max(x),
                                                                size = ceiling(mle_ori[1]), prob = mle_ori[2])
    step_ori <- stats::stepfun(0:max(x), c(0, probs_ori))
    z <- stats::knots(step_ori)
    dev <- c(0, (stats::ecdf(x))(z) - step_ori(z))
    max(abs(dev))
}

###### load count matrix
data <- readRDS("data/counts.rds")

### load TF genes
load(file  ="data/ALL_TF.rda")
TF.gene <- ALL_TF

#### load groups
groups <- readRDS("data/groups.rds")
name.group <- as.character(unique(groups))

##### We only consider neuronal lineage, so remove Sus group
name.groupcond <- name.group[-which(name.group=="Sus")]

########## filter and normalize data for each group
for (k  in 1: length(name.groupcond)) {
    ind.group <- which(groups==name.groupcond[k])
    Y <- t(data[,ind.group])

    ##### filter 500 TF genes with highest mean!
    ind.TF <- which(colnames(Y) %in% TF.gene)
    Y.TF <- Y[,ind.TF]
    mean.TF <- apply(Y.TF, 2, mean)
    or.mean <- order(mean.TF,decreasing = TRUE)
    YY.TF <- Y.TF[,or.mean[1:500]]
    ################# filter 1000 genes with highest mean (not TF genes) from the groups
    ##############    of 50% genes with most variance
    Y.nTF <- Y[,-ind.TF]
    PercentGenes <- .5
    p <- ncol(Y.nTF)
    tG <- floor(p*PercentGenes)
    vars <- apply(Y.nTF,2,var)
    or <- order(vars,decreasing=TRUE)
    YY.nTF <- Y.nTF[,or[1:tG]]
    meanYY.nTF  <- apply(YY.nTF,2,mean)
    ord <- order(meanYY.nTF ,decreasing=TRUE)
    YYY.nTF  <- YY.nTF [,ord[1:1000]]

    ########### normalizing data with 1500 genes
    X.genes <- cbind(YYY.nTF,YY.TF)
    #### Step1: normalizing with 95%quantile matching
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
    ##### filter 1000 cells with highest mean
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

    ####  transform X to X^\alpha with alpha obtained from KS test
    datatest <- floor(fX^alpha.ks[ind.alpha])
    gene.skip <- which(apply(datatest,2,mean)==0)
    if (length(gene.skip)>0) datatest <- datatest [,-gene.skip]
    name.file <- paste(k,"data1500transcels595rep.RData", sep = "-")
    save(datatest,file=name.file)
}

################### Estimating gene networks by zinb1

source("optimzinb.R")
source("PCzinb1noT.R")
for (k in 1:length(name.groupcond)) {
    name.file <- paste(k,"data1500transcels595rep.RData", sep = "-")
    load(name.file)
    n <- dim(datatest)[1]
    p <- dim(datatest)[2]
    datatest <- as.matrix(datatest,n,p)
    #### upper bound for conditional independence sets
    maxcard <- 3
    ## level of significance tests
    alpha <-  2*pnorm(n^.2,lower.tail=F)
    #####zinb1
    zinbPC1.time <- system.time(adj.zinb1 <- zinb1.noT(datatest,maxcard,alpha, extend=TRUE))
    name.fileresult <- paste(k,"data1500transcels595rep2result.RData", sep = "-")
    save(adj.zinb1,datatest, file = name.fileresult )

}
