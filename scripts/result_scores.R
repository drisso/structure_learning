### function to compute TP,FP,FN,PPV,Se, and F1 score
result = function(trueG ,estimatedG){
  TP <-  sum(trueG *estimatedG)
  FP <-  sum((estimatedG-trueG)==1)
  FN <- sum((estimatedG-trueG)==-1)
  PPV  <- TP/(TP+FP)
  Se <- TP/(TP+FN)
  F1 <- 2*(PPV*Se)/(PPV+Se)
  return (c( TP,FP,FN,PPV,Se,F1))
}