
###################################################
## Raw 2 x 2 contingency table analysis function ##
## Name: rbayes
##
## Syntax: x is a table specified as:
##         table(test(B), outcome(A))
##
## Table structure assuming binary data:
## 
##                         OUTCOME (A)
##             ---------------------------
##             |       | - (0) | + (1)   |
##             ---------------------------
##             | - (0) |  TN   |   FN    |
##   TEST (B)  ---------------------------  
##             | + (1) |  FP   |   TP    |
##             ---------------------------
##
##  NB. the table function in R sorts the binary
##  responses so the table looks odd compared to
##  traditional contingency tables
##
##



## Full 2 x 2 contingency table analysis function - primarily bayesian

rbayes <- function(x, ...){
  UseMethod("rbayes")
}

rbayes.default <- function(x,priors=.5){

  # convert table to data.frame for easier data extraction
  dat <- as.data.frame(x)
  
  tn <- dat[1,3] # TN
  fp <- dat[2,3] # FP
  fn <- dat[3,3] # FN
  tp <- dat[4,3] # TP
  sens <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  pvp <- tp/(fp+tp)
  pvn <- tn/(tn+fn)
  tot <- sum(dat[,3])
  eff <- (tp+tn)/tot # added on 7/12/06
  
  # Signal detection results

  sdt <- list(sensitivity=round(sens,2),specificity=round(spec,2),PPV=round(pvp,2),NPV=round(pvn,2),efficiency=round(eff,2))
  
  # Bayesian calculations
  # Let A = model goodness and B = Fit test

  # Good Models
  p.A.obs <- (fn+tp)/tot
  p.A.set <- priors # model priors
  p.notA.obs <- 1 - p.A.obs
  p.notA.set <- 1 - priors

  p.BgA <- tp/(tp+fn)
  p.BgnotA <- fp/(tn+fp)
  p.notBgA <- 1 - p.BgA # always based upon observed data
  p.notBgnotA <- 1 - p.BgnotA # likewise based upon observed data
  
  # Good Fit
  p.B.obs <- (fp+tp)/tot
  p.notB.obs <- 1 - p.B.obs
  p.B.set <- (p.BgA * p.A.set) + (p.BgnotA * p.notA.set)
  p.notB.set <- 1 - p.B.set

  ## the four posterior probabilities and their associated versions

  p.AgB.obs <- (p.BgA * p.A.obs)/ p.B.obs
  p.AgB.set <- (p.BgA * p.A.set) / p.B.set

  p.notAgB.obs <- (p.BgnotA * p.notA.obs) / p.B.obs
  p.notAgB.set <-  (p.BgnotA * p.notA.set) / p.B.set

  p.notAgnotB.obs <- (p.notBgnotA * p.notA.obs) / p.notB.obs
  p.notAgnotB.set <- (p.notBgnotA * p.notA.set) / p.notB.set

  p.AgnotB.obs <- (p.notBgA * p.A.obs) / p.notB.obs
  p.AgnotB.set <- (p.notBgA * p.A.set) / p.notB.set


  # results from bayesian analysis

  bayes <- list("comppriors"=round(p.A.obs,2),"specpriors"=round(p.A.set,2),"P(+Outcome|+Test) w/obs priors"=round(p.AgB.obs,2),"P(+Outcome|+Test) w/spec priors"=round(p.AgB.set,2),"P(-Outcome|+Test) w/obs priors"=round(p.notAgB.obs,2),"P(-Outcome|+Test) w/spec priors"=round(p.notAgB.set,2),"P(-Outcome|-Test) w/obs priors"=round(p.notAgnotB.obs,2),"P(-Outcome|-Test) w/spec priors"=round(p.notAgnotB.set,2),"P(+Outcome|-Test) w/obs priors"=round(p.AgnotB.obs,2),"P(+Outcome|-Test) w/spec priors"=round(p.AgnotB.set,2))

  # Now some additional stats to return, OR, RR, Kappa, Overall
  # fraction correct, mis-classification rate, NNT, Absolute Risk
  # reduction, Relative risk reduction, Positive LR, Negative LR,
  # Diagnostic OR, Error OR, Youden's J, NND, Forbes NMI, contingency
  # coef, adjusted contingency coef, phi coef, Yule's Q.

  # Added kappa and McNemar's test on 1/13/06 

  # SEE:  http://members.aol.com/johnp71/ctab2x2.html
  
  oddsratio <- (tp/fp)/(fn/tn)
  relrisk <- (tp/(tp+fp))/(fn/(fn + tn))
  relriskSE <- sqrt((fp/(tp*(tp+fp))) + (tn/(fn*(tn+fn))))
  relriskCI <- c(relrisk*exp(-1.96*relriskSE),relrisk*exp(1.96*relriskSE))
  kappa <- ((tp/tot + tn/tot) - ((tp/tot + fp/tot)*(tp/tot + fn/tot) + (fn/tot + tn/tot)*(tn/tot + fp/tot)))/(1 - ((tp/tot + fp/tot)*(tp/tot + fn/tot) + (fn/tot + tn/tot)*(tn/tot + fp/tot)))
  oafraccor <- (tp + tn)/tot
  missrate <- 1 - oafraccor
  NNT <- 1 / ((tp/(tp+fp)) - (fn/(fn + tn)))
  ARR <- (fn/(fn + tn)) - (tp/(tp+fp))
  RRR <- ARR/(fn/(fn + tn))
  PosLR <- sens/(1 - spec)
  NegLR <- (1 - sens) / spec
  DiagOR <- (sens/(1 - sens))/((1 - spec)/spec)
  YoudenJ <- sens + spec - 1
  NND <- 1/YoudenJ
  YulesQ <- (oddsratio - 1)/(oddsratio + 1)
  McNemars <- ((fn - fp)^2)/(fn + fp)
  McNemars.p <- round(pchisq(McNemars,df=1,lower.tail=F),digits=2)
  
  # Collect statistics for output and combine all three sets to display
  DxInds <- list(OR=oddsratio,RR=relrisk,RRSE=relriskSE,RRCI=relriskCI,Kappa=kappa,"Overall Fraction Correct"=oafraccor,"Miss Rate"=missrate,NNT=NNT,ARR=ARR,RRR=RRR, PosLR=PosLR, NegLR=NegLR, DiagOR=DiagOR, "Youden's J"=YoudenJ, NND=NND, "Yule's Q"=YulesQ,McNemars=McNemars,"McNemar's p"=McNemars.p)
  
  out <- list("Signal Detection Theory Results"=sdt,"Bayesian Results"=bayes,"Misc. Diagnostic Indicators"=DxInds)
  class(out) <- "rbayes"
  return(out)
}


summary.rbayes <- function(x){
  sdt <- x[[1]]
  bayes <- x[[2]]
  
  cat("\n Signal Detection Theory Results: \n\n")
  
}


#summary.ctt <- function(x){
#  itab <- data.frame(Mean=round(x$itemX,2),SD=round(x$itemSD,2),Rit=round(x$alpha$rit,2),Alpha.Del=round(x$alpha$alpha.rem,2),Discrim=round(x$discrim,2))
#  stab <- data.frame("Cronbach Alpha"=x$alpha$alpha,"Standardized Alpha"=x$alpha$zalpha,"Mean Inter-Item Corr"=x$alpha$Xrii)
#  cat("\n Classical Test Theory Item Statistics: \n\n")
#  print(itab,digits=2)
#  cat("\n Classical Test Theory Scale Statistics: \n\n")
#  print(stab,digits=2)
##   cat("\n\n Response Category Usage: \n\n")
##  print(x$catuse)
#  invisible(x)
#}

rbayes.plot <- function(x){

  # this function will plot the results of the rbayes function in terms of an ROC perhaps
}



## END rbayes
