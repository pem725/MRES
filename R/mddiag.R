##########################################################################
## Missing Data Diagnostic Routines                                     ##
##                                                                      ##
## Purpose:  Analyze missing data patterns, plot patterns, test for     ##
##           for most stringent MCAR assumption                         ##
## Version 0.4                                                          ##
## Date Updated: 03/29/05                                               ##
## Author:  Patrick E. McKnight, Ph.D.                                  ##
##          Measurement Research Evaluation and Statistics Group (MRES) ##
##          Department of Psychology                                    ##
##          George Mason University                                     ##
##          Fairfax, VA                                                 ##
##          pmcknigh@gmu.edu                                            ##
##                                                                      ##
## Distributed under the GPL                                            ##
## Copyright (C) 2011 Patrick E. McKnight                               ##
##                                                                      ##
##    This program is free software; you can redistribute it and/or     ##
##    modify it under the terms of the GNU General Public License as    ##
##    published by the Free Software Foundation; either version 2 of    ##
##    the License, or (at your option) any later version.               ##
##                                                                      ##
##    This program is distributed in the hope that it will be useful,   ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of    ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
##    GNU General Public License for more details.                      ##
##                                                                      ##
##    You should have received a copy of the GNU General Public         ##
##    License along with this program; if not, write to the Free        ##
##    Software Foundation, Inc., 59 Temple Place, Suite 330,            ##
##    Boston, MA  02111-1307  USA                                       ##
##########################################################################

# Create bogus data
id <- seq(1:10)
age <- c(25,29,20,21,20,24,26,29,NA,30)
sex <- c(0,1,0,0,1,NA,1,1,0,0)
iv1 <- c(3,2,4,5,1,3,7,8,9,NA)
iv2 <- c(NA,3,8,9,1,2,3,5,7,9)
dv <- c(25,22,23,NA,15,16,22,NA,30,26)
miss1 <- data.frame(id,age,sex,iv1,iv2,dv)
miss2 <- as.matrix(miss1)

# Generic function

mddiag <- function(x, ...){
  if(is.null(class(x))) class(x) <- data.class(x)
  UseMethod("mddiag")
}

mddiag.default <- function(x){
  # Required libraries
  ##library(mvnmle) ## NOT AVAILABLE as of 9/16/13
##  library(norm)  - cannot compile norm on Gutsy

  # Input checking
  if (missing(x) || length(x) == 0) stop("'x' must be a either a data frame or vector")

  ##################################
  # Dummy code missing values as 1 #
  ##################################

  # x is the vector, data.frame, or data matrix
  # Dmat is the dummy matrix that results from the coding
  
  dat <- as.matrix(x)
  Dmat <- as.matrix(1 * is.na(dat))

  #############################################
  # Compute Different Amounts of Missing Data #
  #############################################
  
  # Compute complete case percent
  cc <- complete.cases(x)
  mode(cc) <- "single"
  ccsum <- sum(cc)
  ccperc <- mean(cc)*100

  # Available case method
  dt1 <- t(Dmat)
  varm <- as.integer((dt1 %*% rep(1,ncol(dt1))))
  name.list <- names(x)
  acperc <- (varm/(ncol(dt1)))*100
  
  # Sparse matrix method
  cols <- ncol(Dmat)
  rows <- nrow(Dmat)
  full <- cols * rows
  missing <- sum(Dmat)
  spmat <- (missing/full)*100

  # Ratio of sparse matrix percent to complete case percent
  cc2spmat <- spmat/ccperc

  #####################################
  #### Pattern detection and summary ##
  #####################################

  # Identify patterns
  # J = number of distinct patterns
  # mj = the number of cases for each pattern
  # Sj = a list of vectors identifying the cases (i.e., rownames) within each pattern
  # pj = the number of variables for each pattern
  # Dj = a list of matrices corresponding to the observed varaibles for each pattern
  # Yobsj = a list of vectors containing the arithmetic means for the observed variables for each pattern

  # Little's test only if specified, otherwise keep NA's from above (note I am using norm to estimate the ML parms)
##  ml.out <- mlest(dat)
    # getparam.norm(prelim.norm(dat),em.norm(prelim.norm(dat),showits=F))      # mlest(x)

  ## three lines below are broken by mlest and norm not working
  muhat <- 0 # ml.out$mu       # ml.out$muhat
  sigmahat <- 0 # ml.out$sigma    # ml.out$sigmahat
  sigmatilda <- 0 # (sigmahat*rows)/(rows - 1)

  mdp <- Dmat %*% (2^(1:ncol(Dmat)))
  J <- length(unique(mdp))
  mj <- as.vector(table(mdp))
  mj.mdp <- as.vector(sort(unique(mdp)))
  mdp.names <- data.frame(row.names(x),mdp)
  pj <- 0
  pj.1 <- cbind(Dsum=(ncol(Dmat) - rowSums(Dmat)),mdp)
  pj <- aggregate(pj.1[,1],list(pj.1[,2]),mean)[,2]

  # initialize lists and vectors
  Sj <- list(0)
  Dj <- list(0)
  Yobsj <- list(0)
  muobsj <- list(0)
  sigmaobsj <- list(0)
  patvars <- list(0)
  delvars <- list(0)
  dvector <- 0
  
  for(i in 1:J){
    obs <- c(1:length(mdp))
    Sj[[i]] <- obs[mdp==mj.mdp[i]]
    Dj.1 <- diag(1,cols,cols)
    patvars[[i]] <- as.data.frame(Dmat)[Sj[[i]][1],]

    varnums <- seq(1:cols)
    delv <- varnums*patvars[[i]]
    delvars[[i]] <- delv[delv > 0]

    if (length(delvars[[i]]==0)){
      Dj[[i]] <- Dj.1
    }
    if (length(delvars[[i]]) > 0){
      Dj[[i]] <- as.matrix(as.data.frame(Dj.1)[,-delvars[[i]]])
    }
    
    Xbar <- mean(x[Sj[[i]],])
    Yobsj[[i]] <- as.vector(Xbar[!is.na(Xbar)])
    muobsj[[i]] <- muhat %*% Dj[[i]]
    sigmaobsj[[i]] <- t(Dj[[i]]) %*% sigmatilda %*% Dj[[i]]
    dvector[i] <- mj[i] %*% (Yobsj[[i]] - muobsj[[i]]) %*% ((sigmaobsj[[i]])^(-1)) %*% t(Yobsj[[i]] - muobsj[[i]])
  }

  # Test of MCAR

  dsquared <- sum(dvector)
  df <- sum(pj) - cols
  chisq.prob <- 1 - pchisq(abs(dsquared),df=df)
    
  # Messy Index computation - ratio of patterns to missing cases
  # use patterns info from previous function to calculate the messy index and use original data
    
  mindex <- J/(rows - ccsum)

  res <- list(Dmat=Dmat,
              "Sample Size"=rows,
              "Complete Cases"=ccsum,
              "Percent of Complete Cases"=ccperc,
              "Available Case Percent"=acperc,
              "Sparse Matrix Amount"=spmat,
              "Complete cases to Sparse Matrix"=cc2spmat,
              "Number of Unique Patterns"=J,
              "Messy Missing Data Index"=mindex,
              "MCAR dvector"=dvector,
              "MCAR dsquared"=dsquared,
              "MCAR df"=df,
              "MCAR chisq.prob"=chisq.prob,
              mdp=mdp,
              mj=mj,
              mj.mdp=mj.mdp,
              mdp.names=mdp.names,
              Sj=Sj,
              pj=pj,
              Dj=Dj,
              Yobsj=Yobsj,
              ml.out=ml.out,
              muhat=muhat,
              sigmahat=sigmahat,
              sigmatilda=sigmatilda,
              muobsj=muobsj,
              sigmaobsj=sigmaobsj)

  class(res) <- "mddiag"
  invisible(res)
}

print.mddiag <- function(x){
  print(unlist(x[c(2:4,6:9,11:13)]),scipen=20)
  invisible(x)
}

summary.mddiag <- function(x){
  ## Need to put together a good summary for the diagnostics
  print(unlist(x))
  invisible(x)
}


# Graphical Diagnostic Procedures

# Display data matrix and color-coded missing data patterns

plot.mddiag <- function(x){
#  if (is.class(x))
    dat <- as.matrix(x)
    Dmat <- as.matrix(1 * is.na(dat))
    dmat1 <- x$Dmat
    mdp <- x$mdp
    couple <- data.frame(dmat1,mdp) # combine Dmat and mdp
    nc <- couple[rev(order(mdp)),]
    dmat2 <- as.matrix(nc)
    killlast <- dmat2[,(-1 * ncol(dmat2))]
    nd <- t(killlast)
    xdim <- 1:nrow(dmat1)
    ydim <- 1:ncol(dmat1)
    image(z=nd,col=c("black","white"),axes=F,xlab="Variables",ylab="Observations")
##  messed up axes for now, will improve shortly
#   axis(3,at=1:ncol(dmat1),labels=rev(names(dmat1)),line=0,las=2)
#   axis(1,at=1:ncol(dmat1),labels=rev(1:ncol(dmat1)))
#   axis(2,at=1:nrow(dmat1),labels=rev(dmat1[,1]),las=1)
#    output <- list(dmat1,mdp,couple,nc,dmat2,killlast,nd,xdim,ydim)
#    class(output) <- "mddiag"
#    output
  }

### THE FUNCTION BELOW WORKS PERFECTLY

plotMD <- function(x,col=c("black","white"),x.cex=.6,y.cex=.6,main="",lasx=2,xaxisat=1){
    dat <- as.matrix(x)
    vnames <- names(x)
    n <- nrow(x)
    p <- ncol(x)
    Dmat <- as.matrix(1 * is.na(dat))
    xdim <- 1:nrow(Dmat)
    ydim <- 1:ncol(Dmat)
    image(x=1:(p),y=1:n,z=t(Dmat)[,nrow(Dmat):1],col=col,axes=F,xlab="Variables",ylab="Observations",main=main)
    # xaxisat=1 is bottom, =3 is top - default is bottom
    axis(xaxisat, labels = vnames, at = 1:p, tck=0, lwd=0, las=lasx, cex.axis=x.cex)
    axis(2, lwd=0, labels = rev(row.names(x)[seq(0,n,by=round(n/10))]), las=2, at=seq(1,n,by=round(n/10)), pos=.7, hadj=1, cex.axis=y.cex)
}


## plotMD <- function(x,col=c("black","white"),x.cex=0.8,y.cex=0.8,main=""){
##     dat <- as.matrix(x)
##     vnames <- names(x)
##     n <- nrow(x)
##     p <- ncol(x)
##     Dmat <- as.matrix(1 * is.na(dat))
##     xdim <- 1:nrow(Dmat)
##     ydim <- 1:ncol(Dmat)
##     image(x=1:(p),y=1:n,z=t(Dmat)[,nrow(Dmat):1],col=col,axes=F,xlab="Variables",ylab="Observations",main=main)
##     axis(1, labels = vnames, at = 1:p, tck=0, lwd=0, las=1)
##     axis(2, lwd=0, labels = rev(row.names(x)[seq(0,n,by=round(n/10))]), las=2, at=seq(1,n,by=round(n/10)), pos=.7, hadj=1, cex.axis=y.cex)
## }



# Display patterns by a grouping variable in a dotchart (Cleveland, 1993)

plot.dotmiss <- function(x,gvar){
    misst <- table(x$mdp,gvar)
    dotchart(misst,main="Plot of missing data patterns for grouping variable")
  }

#####################################
# Perturb data with these functions #
#####################################

## cmcar <- function(x,prop){
##   outdat <- x
##   for(i in 1:as.integer(prop*prop*nrow(x)*ncol(x))){
##     obs.num <- sample(1:(nrow(x)),1)
##     var.num <- sample(1:(ncol(x)),1)
##     outdat[obs.num,var.num] <- NA
##   }
##   return(outdat)
## }


## new function - disregard the previous
createMCAR <- function(x,percent=.1,idvar=1){
  out <- x
  Nr <- nrow(x)
  Nc <- ncol(x)-1
  Nrm <- round(percent*Nr*Nc)
  sdat <- data.frame(r=rep(1:Nr,Nc),c=as.numeric(gl(Nc,Nr)))
  sdat <- sdat[sdat$c != idvar,]
  sdat <- sdat[sample(1:nrow(sdat),Nrm),]
  for(i in 1:nrow(sdat)){
    out[sdat[i,1],sdat[i,2]] <- NA
  }
  return(out)
}




cmar <- function(x,cvar,prop,rmvars){
  outdat <- x
  vtmp <- seq(1:ncol(x))
  vars <- vtmp[-cvar]
  cscores <- sort(outdat[,cvar])
  cut1 <- cscores[(prop*length(cscores))]
  cut2 <- cscores[(1 - prop)*length(cscores)]
  c1.diff <- ((length(cscores[cscores < cut1]))/(length(cscores))) - prop
  c2.diff <- ((length(cscores[cscores > cut2]))/(length(cscores))) - prop
  obs.dat <- data.frame(obs=seq(1:nrow(x)),selvar=x[,cvar])

  if (abs(c1.diff) <= abs(c2.diff)){
    rmvector <- subset(obs.dat,selvar < cut1)[,1]
  }
  else {
    rmvector <- subset(obs.dat,selvar > cut2)[,1]
  }

  for (i in 1:length(rmvector)){
    for(j in 1:as.integer(prop*length(rmvars))){
      var.num <- sample(rmvars,1)
      outdat[rmvector[i],var.num] <- NA
    }
  }
  return(outdat)
}

cmnar <- function(x,prop,rmvars){
  outdat <- x
  cutscores <- list(0)
  for(i in 1:length(rmvars)){
    scores <- sort(outdat[,rmvars[i]])
    cut1 <- scores[(prop*length(scores))]
    cut2 <- scores[(1 - prop)*length(scores)]
    c1.diff <- ((length(scores[scores < cut1]))/(length(scores))) - prop
    c2.diff <- ((length(scores[scores > cut2]))/(length(scores))) - prop
    obs.dat <- data.frame(obs=seq(1:nrow(x)),selvar=outdat[,rmvars[i]])

    if (abs(c1.diff) <= abs(c2.diff)){
      rmvector <- subset(obs.dat,selvar < cut1)[,1]
      outdat[rmvector,rmvars[i]] <- NA
    }
    else {
      rmvector <- subset(obs.dat,selvar > cut2)[,1]
      outdat[rmvector,rmvars[i]] <- NA
    }
  }
  return(outdat)
}

## a quick function to be used with psychometric instruments where multiple imputation via mice is performed

psyc.mice <- function(x){
  x.1 <- mean(rowSums(complete(x,1)))
  x.2 <- mean(rowSums(complete(x,2)))
  x.3 <- mean(rowSums(complete(x,3)))
  x.4 <- mean(rowSums(complete(x,4)))
  x.5 <- mean(rowSums(complete(x,5)))

  x.out <- c(x.1,x.2,x.3,x.4,x.5)
  mean.x <- mean(x.out)
  sd.x <- sd(x.out)
  return(list(mean.x,sd.x))
}

  

## END MDDIAG
