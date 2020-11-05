rm (list=ls())
library("countreg")
library(HurdleNormal)
library(zinbgraph)
library(BiocParallel)
library(foreach)
# library(doParallel)
# cl <- makeCluster(7)
# registerDoParallel(cl)
# MC_CORES  <- 7
# options(mc.cores=MC_CORES)

# #setwd("/Users/KimHue/study/zinb/R/zinb_code")
# source("datasimPois.R")
# source("datasimNB.R")
# source("datasimzinb.R")
# source("optimzinb.R")
# source("optimnb.R")
# source("PCnb_optim.R")
# source("PCnbinom.R")
# source("PCPoisson.R")
# source("PCzinb1noT.R")
# source("PCzinb0noT.R")

source("analysis/result_scores.R")
source("analysis/run_algorithm.R")

#######load Adj matrix
load("analysis/randomsample10-03.RData")

##########
n <- 1000
p <- ncol(SAdj)
nsim <- 10
#### upper bound for cardinaities of conditional sets
maxcard <- 8
## level of the test
alpha <- 2*pnorm(n^.15,lower.tail=F)

result.al <- result_algorithms(n, SAdj, nsim, maxcard, alpha, model="poisson",mu=5,mu.nois=0.5,theta=0.5,pi=0.7)
round(apply(result.al,2,mean),2)[c(4,5,6,9,15,16,17,20,26,27,28,31,37,38,39,42,48,49,50,53,59,60,61,64)]

result.al <- result_algorithms(n, SAdj, nsim, maxcard, alpha, model="nb",mu=5,mu.nois=0.5,theta=0.5,pi=0.7)
round(apply(result.al,2,mean),2)[c(4,5,6,9,15,16,17,20,26,27,28,31,37,38,39,42,48,49,50,53,59,60,61,64)]

result.al <- result_algorithms(n, SAdj, nsim, maxcard, alpha, model="zinb",mu=5,mu.nois=0.5,theta=0.5,pi=0.7)
round(apply(result.al,2,mean),2)[c(4,5,6,9,15,16,17,20,26,27,28,31,37,38,39,42,48,49,50,53,59,60,61,64)]

