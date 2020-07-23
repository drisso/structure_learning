result_algorithms <- function(n, SAdj, nsim, maxcard, alpha, model,mu,mu.nois,theta,pi){
  W <- matrix(NA,nsim,66)
  set.seed(123)
  for (i in 1:nsim){
    #########simulated  data under differnt models
    if (model == "poisson"){
      X <- pois.simdata(n, p,SAdj,mu,mu.nois)
    }else if (model== "nb"){
      X <-  nbinom.Simdata (n, p,SAdj,mu,mu.nois,theta)
    }else
      X <- zinb.simdata(n,p,SAdj, mu,mu.nois,theta,pi)
    
    
    ##### estimate with PC-zinb model version 1 no T
    p <- ncol(X)
    
    zinbPC1.time = system.time(adj.zinb1 <- try(zinb1.noT(X,maxcard,alpha, extend=TRUE), silent = TRUE))
    if (class (adj.zinb1) !="try-error"){
      zinbv1 <- result(SAdj,adj.zinb1)
    }else{
      adj.zinb1 <- matrix(NA,p,p)
      zinbv1 <- result(SAdj,adj.zinb1)
    }
    ##### estimate with PC-zinb model version 0 no T
    zinbPC0.time = system.time(adj.zinb0 <- try(zinb0.noT(X,maxcard,alpha, extend=TRUE), silent = TRUE))
    if (class (adj.zinb0) !="try-error"){
      zinbv0 <- result(SAdj,adj.zinb0)
    }else{
      adj.zinb0 <- matrix(NA,p,p)
      zinbv0 <- result(SAdj,adj.zinb0)
    }
    
    ##########estimate with Poisson model
    pois.time = system.time(adj.pois <- try(pois.wald (X,maxcard,alpha,extend = TRUE),silent = TRUE))
    if (class (adj.pois)!="try-error"){
      Pc.pois = result(SAdj,adj.pois)
    }else{
      adj.pois = matrix(NA,p,p)
      Pc.pois = result(SAdj,adj.pois)
    }
    
    ########### estimate with negative binomial
    
    nbPC.time = system.time(adj.nb <- try(nb.wald (X,maxcard,alpha, extend=TRUE), silent = TRUE))
    if (class (adj.nb) !="try-error"){
      nb.PC <- result(SAdj,adj.nb)
    }else{
      adj.nb <- matrix(NA,p,p)
      nb.PC <- result(SAdj,adj.nb)
    }
    ########### estimate with negative binomial optim
    
    nbPCop.time = system.time(adjop.nb <- try(nbscale.noT(X,maxcard,alpha,extend=TRUE), silent = TRUE))
    if (class (adjop.nb) !="try-error"){
      nbop.PC <- result(SAdj,adjop.nb)
    }else{
      adjop.nb <- matrix(NA,p,p)
      nbop.PC <- result(SAdj,adjop.nb)
    }
    
    
    #estimation with zinb-Gauss
    Y <- log(X+1)
    hurdle_fit.time <- system.time(hurdle_fit<- try(fitHurdle(Y,nlambda =20, parallel=TRUE, 
                                                              control=list(debug=0)), silent = TRUE))
    if (class (hurdle_fit) !="try-error"){  
      ind <- max(which(hurdle_fit$BIC == min(hurdle_fit$BIC)))
      adj.Gaus <-as.matrix( hurdle_fit$adjMat[ind][[1]])
      adj.Gaus [adj.Gaus != 0] <- 1
      adj.Gaus <- adj.Gaus+t(adj.Gaus)
      adj.Gaus [adj.Gaus != 0] <- 1
      zinb.Gauss <- result(SAdj, adj.Gaus)
    }else{
      zinb.Gauss <- matrix(NA,p,p)
      zinb.Gauss <- result(SAdj, adj.Gaus)
    }
    
    W[i,] <- c(zinbv1,zinbPC1.time,zinbv0,zinbPC0.time,nbop.PC,
                nbPCop.time,nb.PC,nbPC.time, Pc.pois,pois.time,zinb.Gauss,hurdle_fit.time)
  }
  return(W)
}

