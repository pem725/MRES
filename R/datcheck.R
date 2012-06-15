
###############################
## Look at patterns in data  ##
## dat.check routine         ##
## added 9/22/05             ##
## latest revision:  10/5/06 ##
###############################


## Syntax
## x is matrix or data.frame of independent variables
## dv is a quoted character string that represents a single relevant dependent variable
##
## example:
## dat.check(mydat,dv="Outcome")


dat.check <- function(x,dv=NULL,mdp="./MissPatterns.pdf",hplot="./HistPlots.pdf",biplot="./BivariatePlots.pdf"){

  # reduce dataframe to only numeric
  indices<-1:dim(x)[2]
  indices<-as.numeric(na.omit(ifelse(indices*sapply(x,is.numeric),indices,NA)))
  x.num <- x[,indices]

  # separate ivs from specified dv
  dv.num <- match(dv,names(x.num))
  x.num.iv <- x.num[,-dv.num]
  x.num.dv <- x.num[,dv.num]
  
  # MD patterns
  pdf(file=mdp,bg="gray")
  Dmat <- as.matrix(1 * is.na(x.num))
  mdp <- Dmat %*% (2^(1:ncol(Dmat)))
  couple <- data.frame(Dmat,mdp) # combine Dmat and mdp
  nc <- couple[rev(order(mdp)),]
  dmat2 <- as.matrix(nc)
  killlast <- dmat2[,(-1 * ncol(dmat2))]
  nd <- t(killlast)
  xdim <- 1:nrow(Dmat)
  ydim <- 1:ncol(Dmat)
  image(z=nd,col=c("white","black"),axes=F,xlab="Variables",ylab="Observations")
##   axis(3,at=1:ncol(Dmat),labels=rev(names(Dmat)),line=0,las=2)
##   axis(1,at=1:ncol(Dmat),labels=rev(1:ncol(Dmat)))
##   axis(2,at=1:nrow(Dmat),labels=rev(Dmat[,1]),las=1)
  dev.off()

  # histograms of individual variables
  pdf(file=hplot)
  for (i in 1:ncol(x.num.iv)){
    hist(x.num.iv[,i],xlab=paste(names(x.num.iv)[i]),main=paste(names(x.num.iv)[i]),right=F)
  }
  if (length(dv) > 0){
    hist(x.num.dv,xlab="Dependent Variable Level",main=paste("Specified Dependent Variable:",names(x.num)[dv.num]),right=F)
  }
  dev.off()

  # bivariate relationships on complete cases - still a work in progress
  # pccc.out <- princomp(x.num[complete.cases(x.num),])

  pdf(file=biplot)
  for (i in 1:ncol(x.num.iv)){
    plot(x.num.dv~as.factor(x.num.iv[,i]),xlab=paste(names(x.num.iv)[i]),ylab=paste(names(x.num[dv.num])))
  }
  dev.off()
}
