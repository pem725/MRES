## usage: x is the data.frame, items are the columns for the items,
## crit is the columns for the criterion variable.

carvedscales <- function(x,...){
  UseMethod("carvedscales")
}

carvedscales.default <- function(x,itms,crt,boot=F,Nobs=100,iter=100){
  library(gtools)
   
  Xtot <- rowSums(x[,itms],na.rm=T)
  crit <- x[,crt]
  items <- x[,itms]
  orig.cor <- cor(Xtot,crit)

  ## need to check the number of items
  Nitms <- length(itms)

  ## ESSENTIAL FUNCTIONS BELOW
  
  carvefcn <- function(items,crit,item.fin,total){
    crit.res <- rep(NA,nrow(item.fin))
    ts.res <- rep(NA,nrow(item.fin))
    for (i in 1:nrow(item.fin)){
      tmp <- items[,item.fin[i,]]
      if (is.null(ncol(tmp))){
        crit.res[i] <- 1-cor(tmp,crit)^2
        ts.res[i] <- 1-cor(tmp,total)^2
      } else {
        crit.res[i] <- 1-cor(rowSums(tmp),crit)^2
        ts.res[i] <- 1-cor(rowSums(tmp),total)^2
      }
    }
    out <- list(crit.res=crit.res,ts.res=ts.res)
    return(out)
  }

  carveNOboot <- function(items,crit){
    ## items <- x[,itms]
    ## crit <- x[,crt]
    total <- rowSums(items)
    item.sel <- permutations(2,ncol(items),c(T,F),repeats=T)
    item.fin <- item.sel[-1,] # only the permutations that matter
    out <- carvefcn(items,crit,item.fin,total) # call core fcn
    crit.res <- out[[1]]
    ts.res <- out[[2]]
    crit.best <- match(min(crit.res),crit.res)
    total.best <- match(min(ts.res[-length(ts.res)]),ts.res[-length(ts.res)])
    crit.best.items <- names(items)[item.fin[crit.best,]]
    total.best.items <- names(items)[item.fin[total.best,]]
    res <- list(crit.res=crit.res,total.res=ts.res,items=item.fin,item.names=names(items),crit=names(x)[crt],crit.best=crit.best,total.best=total.best,crit.best.items=crit.best.items,total.best.items=total.best.items,boot=F,orig.cor=orig.cor)
    return(res)
  }

  carveboot <- function(items,crit,Nobs,iter){
    out.dat <- data.frame(crit.best=rep(NA,iter),pattern=rep(NA,iter))
    for (i in 1:iter){
      sam <- sample(1:nrow(items),Nobs,T)
      itmp <- items[sam,]
      crtmp <- crit[sam]
      total <- rowSums(itmp)
      item.sel <- permutations(2,ncol(itmp),c(T,F),repeats=T)
      item.fin <- item.sel[-1,]
      out <- carvefcn(itmp,crtmp,item.fin,total) # call core fcn
      crit.res <- out[[1]]
#      cc.dat <- data.frame(cc=crit.cor,cc.abs=abs(crit.cor))
#      min.res <- cc.dat[order(cc.dat$cc.abs,decreasing=T),1][1]
      min.res <- min(crit.res)
      out.dat[i,] <- c(min.res,match(min.res,crit.res))
    }
    res <- list(out.dat=out.dat,items=item.fin,item.names=names(items),crit=names(x)[crt],boot=T,orig.cor=orig.cor)
    return(res)
  }

  ItmChop <- function(x){
    ## determine optimal cut points for item number by selecting the
    ## group with the smallest remainder.
    ocp.rem <- x
    iters <- 8:5
    for (i in 8:5){
      ocp.rem <- c(ocp.rem,x %% i)
    }
    Optgs <- iters[match(min(ocp.rem),ocp.rem[-1])]
    repeat {
      grp <- round(runif(x,1,round(x/Optgs)))
      if (max(table(grp)) < 10 & min(table(grp)) > 3) break
    }
    res <- data.frame(item=1:x,group=grp)
    return(res)
  }
  
  carveGA <- function(items,crit){ # begin carveGA - the genetic algorithm
    coreGA <- function(items,crit){
      ItmAssign <- ItmChop(ncol(items)) # see essential function above for details
      Nparcels <- max(ItmAssign$group)
      res <- NA # setup a storage container object for winning items
      for (i in 1:Nparcels){
        initItmSel <- ItmAssign$item[ItmAssign$group==i]
          assign(paste("tres",i,sep=""),carveNOboot(items[,initItmSel],crit))
          tmpout <- get(paste("tres",i,sep=""))
#          res <- c(res,initItmSel[tmpout$items[order(abs(tmpout$crit.cor),decreasing=T)[1],]])
           res <- c(res,initItmSel[tmpout$items[order(tmpout$crit.res)[1],]])
      }
      res <- res[-1]
      return(res)
    }
    res <- 1:ncol(items)
    repeat {
      res <- sort(coreGA(items[,res],crit))
      if (length(res) < 11) {
        break
      }
    }
    return(res) # return the item numbers that pass the GA
  } # end of carveGA - the genetic algorithm

  ## test carveGA
  test.cGA <- carveGA(items,crit)
  
  ### END ESSENTIAL FUNCTIONS
  
  ### Primary carvedscales function

  if (Nitms <= 10){
    if (boot==F){
      res <- carveNOboot(items,crit)
    } else {
      res <- carveboot(items,crit,Nobs,iter)
    }
  } else {
    if (boot==F){
      res <- carveNOboot(items[,carveGA(items,crit)],crit)
    } else {
      res <- carveboot(items[,carveGA(items,crit)],crit,Nobs,iter)
    }
  }
  
  class(res) <- "carvedscales"
  invisible(res)
  return(res)
}

summary.carvedscales <- function(x){
  if(x$boot==T){
    Xbar <- aggregate(x$out.dat$crit.best,by=list(x$out.dat$pattern),mean)
    N <- aggregate(x$out.dat$crit.best,by=list(x$out.dat$pattern),length)
    Imin <- aggregate(x$out.dat$crit.best,by=list(x$out.dat$pattern),min)
    Imax <- aggregate(x$out.dat$crit.best,by=list(x$out.dat$pattern),max)
    item.patts <- x$items[unique(x$out.dat$pattern),]
    item.out <- data.frame(Group.1=unique(x$out.dat$pattern),items=rep(NA,nrow(item.patts)))
    for (i in 1:nrow(item.patts)){
      item.out[i,2] <- paste(x$item.names[x$items[item.out[i,1],]],sep=" ",collapse=" ")
    }
    ## for (i in 1:nrow(item.out)){ ## NEED TO GET THE Nitems from each pattern MUST DO
    ## }  
    out <- merge(N,Xbar,by="Group.1")
    out <- merge(out,Imin,by="Group.1")
    out <- merge(out,Imax,by="Group.1")
    out <- merge(out,item.out,by="Group.1")
    names(out) <- c("Patt.Number","N","Mean","Min","Max","items")
    out[,3] <- round(out[,3],2)
    out[,4] <- round(out[,4],2)
    out[,5] <- round(out[,5],2)
    return(out)
    print(out)
  }
  if(x$boot==F){
    maxcor.crit <- max(x$crit.cor)
    maxcor.ts <- max(x$total.cor[-length(x$total.cor)])
    regcor <- x$total.cor[length(x$total.cor)]
    ## ZZZZ fix this 
  }
}

plot.carvedscales <- function(x){
  if (x$boot==T){
    dat <- x$out.dat[order(x$out.dat$pattern),]
    plot(crit.best~as.factor(pattern),data=dat,xlab="Pattern",ylab="Correlation",ylim=c(-1.2,1))
    abline(h=x$orig.cor,col="blue",lwd=2)
    legend("topright",legend=c("Full Scale Correlation"),lty=1,lwd=2,col=c("blue"))
  }
  
  if (x$boot==F){
    par(mfrow=c(1,2))
    hist(x$crit.cor,xlab="Correlation",main="Criterion")
    allcor <- x$crit.cor[length(x$crit.cor)]
    abline(v=allcor,col="red",lwd=2)
    abline(v=x$crit.cor[x$crit.best],col="blue",lwd=2)
    legend("topright",legend=c("Best","Total"),lty=1,lwd=2,col=c("blue","red"))

    hist(x$total.cor,xlab="Correlation",main="Total Score")
    abline(v=x$total.cor[x$total.best],col="blue",lwd=2)

    par(mfrow=c(1,1))
  }
}

## Here is my problem: I need to take all the permutations of a vector
## and apply a function to each permutation.  That is not difficult
## except for the fact that the number of permutations exceeds the
## memory capacity for R when the vector size exceeds 10.  Thus, I
## either need to find a way to store the permutations in an external
## file and then go through the external file - one permutation at a
## time - or create small "sets" of permutations that enable me to
## randomly select sets to compare.  The advantage of the first brute
## force method is that the solution would include all possible
## permuted samples and thus we would have greater confidence in the
## results; the downside is that it might take days to complete the
## analysis when the number of permutations is very high AND it might
## not be possible to run bootstrap estimates using that scenario
## simply due to time.  The advantage of the subsetting procedure is
## memory and time become non-issues but the permutations will not all
## be equally compared.  The best way to sample would be a random walk
## and then run the random walk procedure several times - perhaps
## several 1000 times - just to make sure that the results are not
## based solely on the sample of permuted combinations of items
## selected for analysis.

## Decision: I chose the selection method (option 2) because it allows
## me to stick with R for my coding environment.  I might consider
## FORTRAN if the problem gets too sticky with R.


