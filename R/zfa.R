# Unit weighted factor score computation

zfa <- function(x,...){
  UseMethod("zfa")
}

zfa.default <- function(x,use="complete.obs"){
  outnames <- names(x)
  outdat <- data.frame(V1=rnorm(nrow(x))) # bogus container for output 
  for (i in 1:ncol(x)){
    outdat[,i] <- scale(x[,i]) # standardize variables
  }
  names(outdat) <- outnames
  scores <- rowMeans(outdat,na.rm=T) # compute means of z scores

  itcor <- rep(0,ncol(x))
  for (i in 1:ncol(x)){
    itcor[i] <- cor(x[,i], rowMeans(outdat[,-i],na.rm=T),use=use)
  }

  resid <- data.frame(V1=rep(0,nrow(x))) # bogus container for residual output
  for (i in 1:ncol(x)){
    resid[,i] <- resid(lm(scores~outdat[,i],na.action=na.exclude))
  }

  # calculate internal consistency of the zfa results
  c.alpha <- function(x){
    xmat <- as.matrix(na.omit(x))
    xmat.z <- scale(xmat)
    N <- dim(xmat)[2]
    rii <- cor(xmat,use="pairwise.complete.obs")
    Xrii <- (sum(rii) - nrow(rii))/(nrow(rii)^2 - nrow(rii))
    rit <- rep(NA,N)
    for (i in 1:N){
      rit[i] <- (cor(apply(xmat[,- i],1,sum),xmat[,i],use="pairwise.complete.obs"))
    }

    alpha.rem <- rep(NA,N)
    for (i in 1:N){
      alpha.rem[i] <- ((N-1)/(N-2)) * (1 - ((sum(diag(var(xmat[,-i])))) / (sum(var(xmat[,-i])))))
    }

    alpha <- (N/(N-1)) * (1 -((sum(diag(var(xmat,na.rm=T)))) / (sum(var(xmat,na.rm=T)))))
    zalpha <- (N/(N-1)) * (1 -((sum(diag(var(xmat.z,na.rm=T)))) / (sum(var(xmat.z,na.rm=T)))))

    res <- list(rii=rii,Xrii=Xrii,rit=rit,alpha.rem=alpha.rem,alpha=alpha,zalpha=zalpha)
    invisible(res)
  }

  alpha <- c.alpha(outdat)
  
  res <- list(zdat=outdat,scores=scores,itcor=itcor,resid=resid,alpha=alpha,outnames=outnames) # create list for output
  class(res) <- "zfa"
  invisible(res)
}

summary.zfa <- function(x){
  itcor <- data.frame(ItemTotal = x$itcor)
  row.names(itcor) <- x$outnames
  cat("\n Unit Weighted Factor Score Computation: \n\n")
  print(round(itcor,2))
  cat("\n ------------------------------------------ \n")
  cat(paste("\n Internal Consistency:",round(x$alpha$alpha,2),"\n"))
  cat("\n ------------------------------------------ \n")
  invisible(x)
}

residual.zfa <- function(x){
  return(x$resid)
}

plot.zfa <- function(x){
  plot(x$scores~x$zdat[,1],col=1,main="",ylab="",xlab="",xlim=c(-5,5),ylim=c(-5,5))
  abline(lm(x$scores~x$zdat[,1]),col=1,lwd=2)
  par(new=T)
  for (i in 2:ncol(x$zdat)){
    plot(x$scores~x$zdat[,i],col=i,main="",ylab="",xlab="",xlim=c(-5,5),ylim=c(-5,5),axes=F)
    abline(lm(x$scores~x$zdat[,i]),col=i,lwd=2)
    par(new=T)
  }
  abline(0,1)
}
