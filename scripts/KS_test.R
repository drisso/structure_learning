######## KS tests 
##### finding optimal alpha by applying Kolmogorov-Smirnov test
#### this alpha is use in power transformation: X to X^\alpha

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