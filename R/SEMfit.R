
###########################################################################
## SEM fit statistic computation for relative fit statistics
## 
##
##
##

SEMfit <- function(nullmodel,model){
  TLI <- ((summary(nullmodel)[[1]]/summary(nullmodel)[[2]]) - (summary(model)[[1]]/summary(model)[[2]])) / ((summary(nullmodel)[[1]]/summary(nullmodel)[[2]]) - 1)
  NFI <- (summary(nullmodel)[[1]] - summary(model)[[1]]) / summary(nullmodel)[[1]]
  RFI <- ((summary(nullmodel)[[1]]/summary(nullmodel)[[2]]) - (summary(model)[[1]]/summary(model)[[2]])) / (summary(nullmodel)[[1]]/summary(nullmodel)[[2]])
  IFI <- ((summary(nullmodel)[[1]]/summary(nullmodel)[[2]]) - (summary(model)[[1]]/summary(model)[[2]]))/(summary(nullmodel)[[2]]/summary(model)[[2]])
  PNFI <- (summary(model)[[2]]/summary(nullmodel)[[2]]) * NFI
  out <- list("Tucker-Lewis (>.90)"=TLI,"Bentler-Bonnet Normed Fit Index (>.90)"=NFI,"Bollen 1986 Relative Fit Index (>.90)"=RFI,"Bollen 1989a Incremental Fit Index (>.90)"=IFI,"Parsimonious Normed Fit Index (>.90)"=PNFI)
  out
}

