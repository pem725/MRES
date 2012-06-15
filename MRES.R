###########################################################################
##                                                                       ##
## Name:  MRES.R                                                         ##
##                                                                       ##
## Miscelaneous collection of functions used primarily in                ##
##  psychometric instrument development and evaluation.                  ##
##                                                                       ##
## Maintainer:  Patrick E. McKnight, Ph.D.                               ##
##              pmcknigh@gmu.edu                                         ##
##              George Mason University                                  ##
##              Fairfax, VA                                              ##
##                                                                       ##
## Authors:  Patrick E. McKnight                                         ##
##           Julius Najab                                                ##
##           Mei-Kuang Chen                                              ##
##           Pedro Wolf                                                  ##
##           Michael Menke                                               ##
##           Carrie Wiley                                                ##
##                                                                       ##
## Version Number:  0.9.7                                                ##
## Version Date:  5/23/2008                                              ##
##                                                                       ##
## Functions included:                                                   ##
##    zfa: create unit-weighted factor scores - edited 9/26 to be more   ##
##    thorough and produce better output and summary stats. Fixed itcor  ##
##    on 4/27/08 to exclude item from the total score (i.e., corrected   ##
##    item total correlation seems more appropriate for this type of     ##
##    analysis.                                                          ##
##                                                                       ##
##    ctt: generate classical test theory statistics                     ##
##                                                                       ##
##    ipu: individual and collective plot utility for rasch (bigsteps)   ##
##    and 2/3-PL (bilog/multilog) models                                 ##
##                                                                       ##
##    ioc: an alternative rasch model parameterization program based     ##
##    upon classical psychophysics models                                ##
##                                                                       ##
##    icc: intra-class correlation function for random and               ##
##    fixed-effects models                                               ##
##                                                                       ##
##    mddiag: missing data diagnosis program with both numerical and     ##
##    graphical routines                                                 ##
##                                                                       ##
##    rbayes: a function that generates a nearly complete set of         ##
##    statistics for 2x2 contingency tables.  Added McNemar's test,      ##
##    Kapp, weighted Kappa, Bowker's test of symmetry, and Cochran's     ##
##    Q test on 1/13/06.  ***RENAMED function*** and added summary       ##
##    function and cleaned up the general code to be more consistent     ##
##    with other R functions.                                            ##
##                                                                       ##
##    SEMfit: a function that computes the remainder of the relative     ##
##    fit statistics that the sem package does not provide               ##
##                                                                       ##
##    iia: intensive item analysis - a set of routines that produces     ##
##    different slices of the raw (observed) item scores                 ##
##                                                                       ##
##    dat.check:  routine of intensive data analysis where plots and     ##
##    other numerical methods help give me a better sense out of the     ##
##    dataset. added 9/22/05 edited to reduce graphic size on 10/6/06    ##
##                                                                       ##
##    plottwo:  a function that plots two Bigsteps/Winsteps item         ##
##    files so that an item by item comparison can be made for either    ##
##    two different administrations or for two different versions.       ##
##    added 10/17/05                                                     ##
##                                                                       ##
##    GTcoef:  a function that takes an aov output for a fully-crossed   ##
##    full factorial ANOVA and computes variance components and          ##
##    g-coefficients.  added 6/6/06.  Added the d-coefficient material   ##
##    necessary for the NSMD study.  Now we have a fuller g-study        ##
##    function. added 5/9/2008.                                          ##
##                                                                       ##
##    igc:  a function to compute individual growth curves.  added       ##
##    10/6/06.  edited on 2/10/07 and fixed bugs and generic functions.  ##
##    edited on 6/12/07 and fixed more bugs in the summary function.     ##
##                                                                       ##
##    corSum: a function to summarize correlation matrices in both       ##
##    graphical and numerical terms.  added 11/1/2006.  THERE IS AN      ##
##    ERROR IN THE SUMMARY FUNCTION noted on 5/9/2008.                   ##
###########################################################################

# Unit weighted factor score computation

zfa <- function(x,...){
  UseMethod("zfa")
}

zfa.default <- function(x,use="complete.obs"){
  outdat <- data.frame(V1=rnorm(nrow(x))) # bogus container for output 
  for (i in 1:ncol(x)){
    outdat[,i] <- scale(x[,i]) # standardize variables
  }
  scores <- rowMeans(outdat,na.rm=T) # compute means of z scores

  itcor <- rep(0,ncol(x))
  for (i in 1:ncol(x)){
    itcor[i] <- cor(x[,i], rowMeans(outdat[,-i],na.rm=T),use=use)
  }

  resid <- data.frame(V1=rep(0,nrow(x))) # bogus container for residual output
  for (i in 1:ncol(x)){
    resid[,i] <- resid(lm(scores~outdat[,i],na.action=na.exclude))
  }
  res <- list(zdat=outdat,scores=scores,itcor=itcor,resid=resid) # create list for output
  class(res) <- "zfa"
  invisible(res)
}

summary.zfa <- function(x){
  itcor <- data.frame(ItemTotal = x$itcor)
  cat("\n Unit Weighted Factor Score Computation: \n\n")
  print(round(itcor,2))
  invisible(x)
}

residual.zfa <- function(x){
  return(x$resid)
}

plot.zfa <- function(x,...){
  plot(x$scores~x$zdat[,1],col=1,main="",xlab="Standardized Item Scores",ylab="Unit-Weighted Factor Scores",xlim=c(-5,5),ylim=c(-5,5))
  par(new=T)
  for (i in 2:ncol(x$zdat)){
        plot(x$scores~x$zdat[,i],col=i,axes=F,main="",ylab="",xlab="",xlim=c(-5,5),ylim=c(-5,5))
    par(new=T)
  }
  abline(0,1)
}


# Classical Test Theory (ctt) function

ctt <- function(x, ...){
  UseMethod("ctt")
}

score.test <- function(x,key){
  key <- as.numeric(key)
  omat <- matrix(0,nrow(x),ncol(x))
  for (i in 1:ncol(x)){
    omat[,i] <- match(x[,i],key[i],0)
  }
  omat <- as.data.frame(omat)
  names(omat) <- names(x)
  return(omat)
}

ctt.discrim <- function(x,scorevar){
  cut <- mean(scorevar)
  cutbin <- rep(0,nrow(x))
  cutbin[scorevar > cut] <- 1
  D <- rep(0,ncol(x))
  for (i in 1:ncol(x)){
    out <- aggregate(x[,i],list(cutbin),mean)
    D[i] <- abs(out[2,2] - out[1,2])
  }
  return(D)
}

mean.response.cat <- function(x,scores){
  out <- aggregate(scores,list(x[,1]),mean,na.rm=T)
  for (i in 2:ncol(x)){
    out <- merge(out,aggregate(scores,list(x[,i]),mean,na.rm=T),by="Group.1",all=T)
  }
  out <- t(out)[-1,]
  row.names(out) <- names(x)
  out <- as.data.frame(out)
  return(out)
}

resp.cat.usage <- function(x){
  out <- table(x[,1])
  for (i in 2:ncol(x)){
    out <- rbind(out,table(x[,i]))
  }
  row.names(out) <- names(x)
  return(out)
}

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

##    crit <- rep(NA,N)
##    for (i in 1:N){
##      crit[i] <- cor(

   
   alpha.rem <- rep(NA,N)
   for (i in 1:N){
     alpha.rem[i] <- ((N-1)/(N-2)) * (1 - ((sum(diag(var(xmat[,-i])))) / (sum(var(xmat[,-i])))))
   }

   alpha <- (N/(N-1)) * (1 -((sum(diag(var(xmat,na.rm=T)))) / (sum(var(xmat,na.rm=T)))))
   zalpha <- (N/(N-1)) * (1 -((sum(diag(var(xmat.z,na.rm=T)))) / (sum(var(xmat.z,na.rm=T)))))
   

   res <- list(rii=rii,Xrii=Xrii,rit=rit,alpha.rem=alpha.rem,alpha=alpha,zalpha=zalpha)
   invisible(res)
 }


ctt.default <- function(x,scorevar=NULL,recode=F,key=NULL){
  # parse out items
 
  # set ID variable 
  # ifelse (id==0, id <- id, id <- seq(1:nrow(x)))

  # set scored object
  ifelse (recode==T, scored <- score.test(x,key), scored <- x)
  
  # set scores variable
  ifelse (recode==T,scores <- rowSums(scored,na.rm=T), ifelse (is.null(scorevar),scores <- rowSums(scored,na.rm=T),scores <- scorevar))

  # start ctt analyses
  alpha <- c.alpha(scored)
  ifelse ((max(as.matrix(scored),na.rm=T) - min(as.matrix(scored),na.rm=T) == 1), discrim <- ctt.discrim(scored,scores), discrim <- rep(NA,ncol(x)))

  # get other summary stats computed
  itemX <- mean(scored,na.rm=T)
  itemSD <- sd(scored,na.rm=T)
  itemVar <- var(scored,na.rm=T)

  # response category usage information
  mrc <- mean.response.cat(x,scores)
  catuse <- resp.cat.usage(x)
  
  res <- list(raw=x,scored=scored,scores=scores,alpha=alpha,discrim=discrim,itemX=itemX,itemSD=itemSD,itemVar=itemVar,mrc=mrc,catuse=catuse)
  class(res) <- "ctt"
  invisible(res)
}

summary.ctt <- function(x){
  itab <- data.frame(Mean=round(x$itemX,2),SD=round(x$itemSD,2),Rit=round(x$alpha$rit,2),Alpha.Del=round(x$alpha$alpha.rem,2),Discrim=round(x$discrim,2))
  stab <- data.frame("Cronbach Alpha"=x$alpha$alpha,"Standardized Alpha"=x$alpha$zalpha,"Mean Inter-Item Corr"=x$alpha$Xrii)
  cat("\n Classical Test Theory Item Statistics: \n\n")
  print(itab,digits=2)
  cat("\n Classical Test Theory Scale Statistics: \n\n")
  print(stab,digits=2)
#   cat("\n\n Response Category Usage: \n\n")
#   print(x$catuse)
  invisible(x)
}

print.ctt <- function(x){
  print(unlist(x),digits=2)
  invisible(x)
}

plot.ctt <- function(x){
  par(bg="grey")
  diff <- seq(0,1,by=.001)
  lcolors <- seq(20,657,by=5)
  prob <- diff*x$discrim[1] + (.5 - x$itemX[1]*x$discrim[1])
  plot(prob~diff,xlim=c(0,1),ylim=c(0,1),xlab="Respondent Ability",ylab="Probability of Endorsement",main="",type="l",col=lcolors[1])
  par(new=T)  
  for (i in 2:length(x$d)){
    prob <- diff*x$discrim[i] + (.5 - x$itemX[i]*x$discrim[i])
    plot(prob~diff,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",type="l",col=lcolors[i])
    par(new=T)
  }
}


##########################################################################
# This section contains two different functions, one to plot all ICC's
# for a test and one to plot individual item ICC with observed values
#----------------------------------------------------------------
# FUNCTION 1:  paicc or "Plot all item characteristic curves"
#----------------------------------------------------------------
#
# Syntax:
#
# data = an R dataframe
# nitmes = the number of items to read
# parms = the number of parameters to read
# legend = do you want a legend on the plot - very messy with a large number of items


paicc <- function(data=NULL, nitems=NULL, parms=NULL, legend=F){
	par(bg="grey")
	theta <- seq(-4, 4, by=.001)
	lcolors <- seq(20,657, by=5)

# Rasch Model Plots
        
if(parms==1){
  
	for(i in 1:nitems){
          	prob <- (1/(1+exp(-1.7*(theta - (data$V2[i])))))
              
		plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1),
		col=lcolors[i], main="ICC's for all items")
           
		par(new=T)
	}
        
        if(legend==T){
          legend(-4,1,legend=c(data$V15[1:nitems]),lty=c(rep(1,nitems)),col=c(lcolors[1:nitems]))
        }
}
        

# Birnbaum Model Plots
        
if(parms==2){
	for(i in 1:nitems){
		prob <- (1/(1+exp(-(data$disc[i])*(theta - (data$diff[i])))))
                plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1),
		col=lcolors[i], main="ICC's for all items")
                
		par(new=T)
	}
        if(legend==T){
          legend(-4,1,legend=c(data$inum[1:nitems]),lty=c(rep(1,nitems)),col=c(lcolors[1:nitems]))
        }
}

# 3PL Birnbaum Model Plots
        
if(parms==3){
	for(i in 1:nitems){
		prob <- (1/(1+exp(-(data$disc[i])*(theta - (data$diff[i]))))) + guess[i]
                plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1),
		col=lcolors[i], main="ICC's for all items")
                par(new=T)
	}
        if(legend==T){
          legend(-4,1,legend=c(data$inum[1:nitems]),lty=c(rep(1,nitems)),col=c(lcolors[1:nitems]))
        }
}
par(new=F, bg="white")
}


#--------------------------------------------------------------------------------
# FUNCTION 2: piicc or "Plot Individual Item Characteristic Curves"
#--------------------------------------------------------------------------------


# Syntax:
# idat = R object (dataframe) containing item level data
# pdat = file containing the raw responses - doesn't make sense if it is not binary right now
# edat = file containing the person (ability) parameters
# item = the item number you wish to plot
# idlength = the offset from the first column to the first item - often used for the person ID

piicc  <- function(idat,pdat,edat,item,parms,offset=0){
	theta <- seq(-4, 4, by=.001)
	selout <- offset + item - 1
        if(parms==1){
          prob <- (1/(1+exp(-1.7*(theta - (idat$diff[item])))))
        }
        if(parms==2){
          prob <- (1/(1+exp(-(idat$disc[item])*(theta - (idat$diff[item])))))
        }
        if(parms==3){
          prob <- (1/(1+exp(-(idat$disc[item])*(theta - (idat$diff[item]))))) + guess[item]
        }
	par(xaxp=c(-4,4,20),yaxp=c(0,1,10),xaxs="i",yaxs="i")
	plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1), main="ICC and observed data for item", col="blue" )
	par(new=T)
	responses <- read.fwf(pdat,width=c(selout,1))
	edat <- scan(edat, list(ptheta=0, junk=""))
	freqdat <- data.frame(resp=responses$V2,ptheta=signif(edat$ptheta, digits=1))
	gres <- aggregate(freqdat, list(agtheta=freqdat$ptheta), FUN=mean)
	gresn <- aggregate(freqdat, list(agtheta=freqdat$ptheta), FUN=sum)
	agtheta <- as.vector(gres$agtheta)
	symbols(agtheta, gres$resp, circles=gresn$resp, xlab="", ylab="", xlim=c(-4,4), ylim=c(0,1), col = "red")
	par(new=F)
} 

## END IPU


ioc <- function(x,ID=NULL,first.item=NULL,nitems=NULL,responses=NULL,anchor=1){

# -----------------------------------------------------------------------------
#  IOC function to compute both item and person parameters from polytomous data
#  Authors:  Patrick E. McKnight, Ph.D.
#            Dept. of Psychology
#            University of Arizona
#
#            Robert Massof, Ph.D.
#            Lions Vision Center
#            Johns Hopkins University
# -----------------------------------------------------------------------------
#    x = data.frame to analyze

#    ID = subject identifier - very important for person parameters

#    first.item = what column contains the first item - assumes the
#    items are contiguously organized with no unnecessary variables
#    between the ones to be analyzed.

#    nitems = the number of items to analyze

#    responses = the possible responses to the items.  For example, a
#    five-response category Likert scale would use the following code:
#    responses=c(1,2,3,4,5).

#    anchor = the item you wish to anchor all item parameters.  The
#    default is item 1 but it may be specified to correspond to any item
#    in the dataset.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#  STEP ONE:  Compute Item Parameters
# -----------------------------------------------------------------------------
  
# Convert dataframe into matrix and get a list of names to use for the variables

  attach(x)
  dat <- as.matrix(x)
  item.dat <- dat[,first.item:(first.item + nitems - 1)]
  inames <- dimnames(item.dat)[[2]]

# Create Table of zscore transformed Category Endorsements to be Analyzed

  maxresp <- max(responses)
  dimresp <- dim(as.array(responses))
  zitems <- matrix(nrow=maxresp,ncol=nitems)
	for (i in 1:nitems){
          for (j in 1:dimresp){
              normdat <- qnorm(cumsum(table(item.dat[,i]))/(sum(table(item.dat[,i]))))                
              zitems[j,i] <- normdat[j]
            }
        }
#  Trim off the "Inf" from the matrix

  zfin <- zitems[1:dim(zitems)[1] - 1,]
  zfin[zfin=="Inf"] <- NA
  
# Construct two matrices to fill with slope and intercept values

  slopemat <- matrix(nrow=nitems,ncol=nitems)
  intmat <- matrix(nrow=nitems,ncol=nitems)
  
# Now populate the matrices with slope and intercept values from the bivariate regressions

  	for (i in 1:nitems){
  	  for (j in 1:nitems){
  	    iv <- zfin[,i]   
	    dv <- zfin[,j]
	      lm1 <- lm(dv~iv)
	      slopemat[j,i] <- as.numeric(coef(lm1)[2])
	      intmat[j,i] <- as.numeric(coef(lm1)[1])
	  }
	}

# Construct new matrix to hold mean slopes and populate matrix with the difference values

  xslope <- matrix(nrow=nitems,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:nitems){
      xslope[j,i] <- (log(slopemat[i,j])-log(slopemat[j,i]))/2
    }
  }

# Anchor the scalar

  lnsd <- exp(xslope[,anchor])
  deltax <- matrix(nrow=nitems,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:nitems){
      deltax[j,i] <- ((intmat[j,i])*(lnsd[j]))
    }
  }

# Means minus anchor matrix

  xminanc <- matrix(nrow=nitems,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:nitems){
      xminanc[j,i] <- deltax[anchor,i]-deltax[j,i]
    }
  }

# Item Parameter Results

  parms <- matrix(nrow=nitems,ncol=3)
  for (i in 1:nitems){
      parms[i,1] <- inames[i]
      parms[i,2] <- mean(xminanc[i,]) * -1
      parms[i,3] <- sd(xminanc[i,])
  }
  
item.parms <<- data.frame(item=parms[,1],difficulty=as.numeric(parms[,2]),sd=as.numeric(parms[,3]))


# -----------------------------------------------------------------------------
# STEP TWO:  Compute Person Parameters
# -----------------------------------------------------------------------------

# Convert dataframe into matrix and get a list of names to use for the variables

  person.dat <<- t(dat[,first.item:(first.item + nitems - 1)])
  pnames <<- ID

# Create Person Table of zscore transformed Category Endorsements to be Analyzed

   maxresp <<- max(responses)
   dimresp <<- dim(as.array(responses))
   zitems <<- matrix(nrow=maxresp,ncol=nrow(item.dat))
 	for (i in 1:nrow(item.dat)){
           for (j in 1:dimresp){
               normdat <<- qnorm(cumsum(table(person.dat[,i]))/(sum(table(person.dat[,i]))))                
               zitems[j,i] <<- normdat[j]
             }
         }

# #  Trim off the "Inf" from the matrix

#   zfin <<- zitems[1:dim(zitems)[1] - 1,]
#   zfin[zfin=="Inf"] <<- NA
  
# # Construct two matrices to fill with slope and intercept values

#   slopemat <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   intmat <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
  
# # Now populate the matrices with slope and intercept values from the bivariate regressions

#   	for (i in 1:nrow(item.dat)){
#   	  for (j in 1:nrow(item.dat)){
#   	    iv <<- zfin[,i]   
# 	    dv <<- zfin[,j]
# 	      lm1 <<- lm(dv~iv)
# 	      slopemat[j,i] <<- as.numeric(coef(lm1)[2])
# 	      intmat[j,i] <<- as.numeric(coef(lm1)[1])
# 	  }
# 	}

# # Construct new matrix to hold mean slopes and populate matrix with the difference values

#   xslope <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   for (i in 1:nrow(item.dat)){
#     for (j in 1:nrow(item.dat)){
#       xslope[j,i] <<- (log(slopemat[i,j])-log(slopemat[j,i]))/2
#     }
#   }

# # Anchor the scalar

#   lnsd <<- exp(xslope[,anchor])
#   deltax <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   for (i in 1:nrow(item.dat)){
#     for (j in 1:nrow(item.dat)){
#       deltax[j,i] <<- ((intmat[j,i])*(lnsd[j]))
#     }
#   }

# # Means minus anchor matrix

#   xminanc <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   for (i in 1:nrow(item.dat)){
#     for (j in 1:nrow(item.dat)){
#       xminanc[j,i] <<- deltax[anchor,i]-deltax[j,i]
#     }
#   }

# # Item Parameter Results

#   parms <<- matrix(nrow=nrow(item.dat),ncol=3)
#   for (i in 1:nrow(item.dat)){
#       parms[i,1] <<- pnames[i]
#       parms[i,2] <<- mean(xminanc[i,]) * -1
#       parms[i,3] <<- sd(xminanc[i,])
#   }
  
# person.parms <<- data.frame(item=parms[,1],difficulty=as.numeric(parms[,2]),sd=as.numeric(parms[,3]))

  
}

# Create test data

testme <- data.frame(i1=c(1,2,3,4,3,2,1,4,5,4,2,4,5,3,5),i2=c(2,3,5,1,3,5,2,2,5,2,4,5,2,4,2))

### END IOC


## Intraclass correlations - both fixed (icc.fix) and random (icc.ran)

icc.fix <- function(rating,judge,person){
     ICC <- lm(rating~judge + person)
     edf <- df.residual(ICC)
     ems <- deviance(ICC)/edf
     bdf <- anova(ICC)[1,1]
     bms <- anova(ICC)[1,2]/bdf
     msb <- bms
     jdf <- anova(ICC)[2,1]
     jms <- anova(ICC)[2,2]/jdf
     k <- jdf + 1
     msw <- ((ems*edf)+(jms*jdf))/(edf+jdf)
     wms <- msw
     n <- bdf + 1
     theta <- (msb - msw)/(k*msw)
     winer11 <- theta / (1+theta)
     winer1k <- (k*theta)/(1+k*theta)
     sf11 <- (bms - wms)/(bms+(k-1)*wms)
     sf21 <- (bms - ems)/((bms)+((k-1)*ems)+((k*(jms-ems))/n))
     sf31 <- (bms - ems)/(bms+((k-1)*ems))
     sf1k <- (bms - wms)/bms
     sf2k <- (bms - ems)/(bms+((jms-ems)/n))
     sf3k <- (bms - ems)/bms
     cat("ICC(1,1)\tICC(2,1)\tICC(3,1)\tICC(1,k)\tICC(2,k)\tICC(3,k)\n")
     cat(sf11,"\t",sf21,"\t",sf31,"\t",sf1k,"\t",sf2k,"\t",sf3k,"\n")
   }

icc.ran <- function(resp,fx,rx){
     require(Devore5,quietly=T)
     require(ctest,quietly=T)
     fix <- factor(fx)
     ran <- factor(rx)
     res0 <- aov(resp~fix+ran)
     res1 <- anova(res0)
     print(res1)
     print(shapiro.test(res0$resid))
     cat("Fixed means\n")
     medias <- tapply(resp,fx,mean)
     print(medias)
     cat("\nTukey multiple comparison of means")
     cat("\n  95% family-wise confidence level\n")
     res2 <- Tukey(res0,conf.level=0.95)
     print(res2$fix)
     ems <- res1[[3,3]]
     bms <- res1[[1,3]]
     jms <- res1[[2,3]]
     k <- res1[[2,1]]+1
     n <- res1[[1,1]]+1
     se <- sqrt(ems/k)*qt(1-(1-0.95)/2,k-1)
     tmp <- barplot(medias,ylim=range(resp),main="Confidence Limits",col="light gray",xlab="Fixed",ylab="Response")
     arrows(tmp,medias-se,tmp,medias+se,code=3,angle=90,length=0.1)
     icc <- k*(jms-ems)/(k*jms+n*bms+(n*k-n-k)*ems)
     cat("\nICC:",icc,"\n")
   }

## END ICC

##########################################################################
## Missing Data Diagnostic Routines                                     ##
##                                                                      ##
## Purpose:  Analyze missing data patterns, plot patterns, test for     ##
##           for most stringent MCAR assumption                         ##
## Version 0.4                                                          ##
## Date Updated: 03/29/05                                               ##
## Author:  Patrick E. McKnight, Ph.D.                                  ##
##          Evaluation Group for Analysis of Data                       ##
##          Department of Psychology                                    ##
##          University of Arizona                                       ##
##          Tucson, AZ  85721                                           ##
##          pem@u.arizona.edu                                           ##
##                                                                      ##
## Distributed under the GPL                                            ##
## Copyright (C) 2005 Patrick E. McKnight                               ##
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
  library(mvnmle)
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

# Display patterns by a grouping variable in a dotchart (Cleveland, 1993)

plot.dotmiss <- function(x,gvar){
    misst <- table(x$mdp,gvar)
    dotchart(misst,main="Plot of missing data patterns for grouping variable")
  }

#####################################
# Perturb data with these functions #
#####################################

cmcar <- function(x,prop){
  outdat <- x
  for(i in 1:as.integer(prop*prop*nrow(x)*ncol(x))){
    obs.num <- sample(1:(nrow(x)),1)
    var.num <- sample(1:(ncol(x)),1)
    outdat[obs.num,var.num] <- NA
  }
  return(outdat)
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
  X11()  
}


############################
##  plottwo function      ##
##  added 10/17/05        ##
##  description in header ##
############################

plottwo <- function(pre,post,skip=3,title,items=1,ident=1,scores=0,xlab=NULL,ylab=NULL){
if (scores==0){
	if (items==1){
		pre <- read.table(pre, header=T, skip=skip)
		post <- read.table(post, header=T, skip=skip)
		pre <- subset(pre,pre$ST==1)
		post <- subset(post,post$ST==1)
		plot(pre$MEASURE,post$MEASURE,xlim=c(-5,5),ylim=c(-5,5),xlab=xlab,ylab=ylab,main=title)
		par(new=T)
		x1 <- seq(-5,5)
		y1 <- seq(-5,5)
		plot(x1,y1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l")
		xse <- mean(pre$ERROR)
		yse <- mean(post$ERROR)
		x2plus <- x1 + 3*xse
		x2minus <- x1 - 3*xse
		y2plus <- y1 + 3*yse
		y2minus <- y1 - 3*yse
		par(new=T)
		plot(x2minus,y2plus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=T)
		plot(x2plus,y2minus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=F)
		if (ident==1){
			identify(pre$MEASURE,post$MEASURE,labels=pre$NAME)
		}
	}
	if (items==0){
		pre <- read.table(pre, header=T, skip=skip)
		post <- read.table(post, header=T, skip=skip)
		pre <- subset(pre,pre$ST==1)
		post <- subset(post,post$ST==1)
		comb <- merge(pre,post,by="NAME")
		comb <- subset(comb)
		plot(comb$MEASURE.x,comb$MEASURE.y,xlim=c(-5,5),ylim=c(-5,5),xlab=xlab,ylab=ylab,main=title)
		par(new=T)
		x1 <- seq(-5,5)
		y1 <- seq(-5,5)
		plot(x1,y1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l")
		xse <- mean(pre$ERROR)
		yse <- mean(post$ERROR)
		x2plus <- x1 + 3*xse
		x2minus <- x1 - 3*xse
		y2plus <- y1 + 3*yse
		y2minus <- y1 - 3*yse
		par(new=T)
		plot(x2minus,y2plus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=T)
		plot(x2plus,y2minus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=F)
	}	
}
if (scores==1){
		pre <- read.table(pre, header=T, skip=skip)
		post <- read.table(post, header=T, skip=skip)
		pre <- subset(pre,pre$ST==1)
		post <- subset(post,post$ST==1)
		comb <- merge(pre,post,by="NAME")
		comb <- subset(comb)
		plot(comb$SCORE.x,comb$SCORE.y,xlim=c(0,40),ylim=c(0,40),xlab=xlab,ylab=ylab,main=title)
		par(new=T)
		x1 <- seq(0,40)
		y1 <- seq(0,40)
		plot(x1,y1,xlab="",ylab="",xlim=c(0,40),ylim=c(0,40),type="l")
		xse <- sd(comb$SCORE.x)*sqrt(1-.8)
		yse <- sd(comb$SCORE.y)*sqrt(1-.8)
		x2plus <- x1 + 2*xse
		x2minus <- x1 - 2*xse
		y2plus <- y1 + 2*yse
		y2minus <- y1 - 2*yse
		par(new=T)
		plot(x2minus,y2plus,xlab="",ylab="",xlim=c(0,40),ylim=c(0,40),type="l", lty="dashed")
		par(new=T)
		plot(x2plus,y2minus,xlab="",ylab="",xlim=c(0,40),ylim=c(0,40),type="l", lty="dashed")
		par(new=F)
	}
}

## pre post plot by two different versions - plots the item numbers - useful for BRS plots

plotitwobytwo <- function(x1,x2,y1,y2,skip=3,title,xlab=NULL,ylab=NULL){
  library(sfsmisc)

  # Read in data from Big/Winsteps - note, only if files
  pre1 <- read.table(x1, header=T, skip=skip)
  post1 <- read.table(x2, header=T, skip=skip)
  pre2 <- read.table(y1, header=T, skip=skip)
  post2 <- read.table(y2, header=T, skip=skip)
  pre1 <- subset(pre1,pre1$ST==1)
  post1 <- subset(post1,post1$ST==1)
  pre2 <- subset(pre2,pre2$ST==1)
  post2 <- subset(post2,post2$ST==1)
  
  # plot 1
  n.plot(pre1$MEASURE,post1$MEASURE,nam=pre1$NAME,col="red",xlim=c(-5,5),ylim=c(-5,5),xlab=xlab,ylab=ylab,main=title)
#  abline(0,1) # add identity line
  par(new=T) # prepare for additional graph

  # plot 2
  n.plot(pre2$MEASURE,post2$MEASURE,nam=pre2$NAME,col="blue",xlim=c(-5,5),ylim=c(-5,5),xlab="",ylab="",main="")
  par(new=T)

  x1 <- seq(-5,5)
  y1 <- seq(-5,5)

  # old code for identity line - see x1 and y1 below
  plot(x1,y1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l")
  par(new=T)

  # plot 1 standard error bars

  xse1 <- mean(pre1$ERROR)
  yse1 <- mean(post1$ERROR)
  x2plus1 <- x1 + 3*xse1
  x2minus1 <- x1 - 3*xse1
  y2plus1 <- y1 + 3*yse1
  y2minus1 <- y1 - 3*yse1
  par(new=T)
  plot(x2minus1,y2plus1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="red",type="l", lty="dashed")
  par(new=T)
  plot(x2plus1,y2minus1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="red",type="l", lty="dashed")
  par(new=T)

  # plot 2 standard error bars
  xse2 <- mean(pre2$ERROR)
  yse2 <- mean(post2$ERROR)
  x2plus2 <- x1 + 3*xse2
  x2minus2 <- x1 - 3*xse2
  y2plus2 <- y1 + 3*yse2
  y2minus2 <- y1 - 3*yse2
  par(new=T)
  plot(x2minus2,y2plus2,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="blue",type="l", lty="dashed")
  par(new=T)
  plot(x2plus2,y2minus2,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="blue",type="l", lty="dashed")
  par(new=F)
  
}

###############################
##  GTcoef function          ##
##  computes G-coef          ##
##  added: 6/6/06            ##
##  description in header    ##
##  latest revision: 6/6/06  ##
###############################

##############################################################################
### Compute generalizability theory coefficients with the following code:
# x is an aov output object and the ANOVA must be a fully-crossed factorial design

GTcoef <- function(x){
  Nparms <- length(all.vars(formula(x))) - 1
  ExpParms <- c(3,7,15,31)
  Effects <- attributes(terms(x))$term.labels
  DF <- summary(x)[[1]][1:length(Effects),1]
  MS <- summary(x)[[1]][1:length(Effects),3]
  Effects.bin <- data.frame(t(attributes(terms(x))$factors))[,-1]
  
  calcGT <- function(x,Effects,DF,MS,Nparms){
     # create container
     Nx <- DF + 1
     out <- data.frame(DF=DF,MS=MS,VarComp=NA,PercVar=NA)
     index <- seq(nrow(out),1)
     row.names(out) <- attributes(terms(x))$term.labels
     if(Nparms==2){
       out[3,3] <- MS[3]
       out[2,3] <- (MS[2] - MS[3])/(Nx[1])
       out[1,3] <- (MS[1] - MS[3])/(Nx[2])
       out[out[,3] < 0,3] <- 0
       out[,4] <- round(out[,3]/sum(out[,3]),digits=2)*100
     }
     if(Nparms==3){
       out[7,3] <- MS[7]
       out[6,3] <- (MS[6] - MS[7])/(Nx[1])
       out[5,3] <- (MS[5] - MS[7])/(Nx[2])
       out[4,3] <- (MS[4] - MS[7])/(Nx[3])
       out[3,3] <- (MS[3] - (Nx[1]*out[6,3]) - (Nx[2]*out[5,3]) - MS[7])/((Nx[1])*(Nx[2]))
       out[2,3] <- (MS[2] - (Nx[1]*out[6,3]) - (Nx[3]*out[4,3]) - MS[7])/((Nx[1])*(Nx[3]))
       out[1,3] <- (MS[1] - (Nx[2]*out[5,3]) - (Nx[3]*out[4,3]) - MS[7])/((Nx[2])*(Nx[3]))
       out[out[,3] < 0,3] <- 0
       out[,4] <- round(out[,3]/sum(out[,3]),digits=2)*100
     }
     if(Nparms==4){
       out[15,3] <- MS[15]
       out[14,3] <- (MS[14] - MS[15])/(Nx[1])
       out[13,3] <- (MS[13] - MS[15])/(Nx[2])
       out[12,3] <- (MS[12] - MS[15])/(Nx[3])
       out[11,3] <- (MS[11] - MS[15])/(Nx[4])
       out[10,3] <- (MS[10] - (Nx[1]*out[14,3]) - (Nx[2]*out[13,3]) - MS[15])/(Nx[1]*Nx[2])
       out[9,3] <- (MS[9] - (Nx[1]*out[14,3]) - (Nx[3]*out[12,3]) - MS[15])/(Nx[1]*Nx[3])
       out[8,3] <- (MS[8] - (Nx[2]*out[13,3]) - (Nx[3]*out[12,3]) - MS[15])/(Nx[2]*Nx[3])
       out[7,3] <- (MS[7] - (Nx[1]*out[14,3]) - (Nx[4]*out[11,3]) - MS[15])/(Nx[1]*Nx[4])
       out[6,3] <- (MS[6] - (Nx[2]*out[13,3]) - (Nx[4]*out[11,3]) - MS[15])/(Nx[2]*Nx[4])
       out[5,3] <- (MS[5] - (Nx[3]*out[12,3]) - (Nx[4]*out[11,3]) - MS[15])/(Nx[3]*Nx[4])
       out[4,3] <- (MS[4] - (Nx[2]*Nx[3]*out[8,3]) - (Nx[1]*Nx[3]*out[9,3]) - (Nx[1]*Nx[2]*out[10,3]) - (Nx[3]*out[12,3]) - (Nx[2]*out[13,3]) - (Nx[1]*out[14,3]) - MS[15])/(Nx[1]*Nx[2]*Nx[3])
       out[3,3] <- (MS[3] - (Nx[2]*Nx[4]*out[6,3]) - (Nx[1]*Nx[4]*out[7,3]) - (Nx[1]*Nx[2]*out[10,3]) - (Nx[4]*out[11,3]) - (Nx[2]*out[13,3]) - (Nx[1]*out[14,3]) - MS[15])/(Nx[1]*Nx[2]*Nx[4])
       out[2,3] <- (MS[2] - (Nx[3]*Nx[4]*out[5,3]) - (Nx[1]*Nx[4]*out[7,3]) - (Nx[1]*Nx[3]*out[9,3]) - (Nx[4]*out[11,3]) - (Nx[3]*out[12,3]) - (Nx[1]*out[14,3]) - MS[15])/(Nx[1]*Nx[3]*Nx[4])
       out[1,3] <- (MS[1] - (Nx[3]*Nx[4]*out[5,3]) - (Nx[2]*Nx[4]*out[6,3]) - (Nx[2]*Nx[3]*out[8,3]) - (Nx[4]*out[11,3]) - (Nx[3]*out[12,3]) - (Nx[2]*out[13,3]) - MS[15])/(Nx[2]*Nx[3]*Nx[4])
       out[out[,3] < 0,3] <- 0
       out[,4] <- round(out[,3]/sum(out[,3]),digits=2)*100
     }
     if(Nparms==5){
       out[31,3] <- MS[31]
       out[30,3] <- (MS[30] - MS[31])/Nx[1]
       out[29,3] <- (MS[29] - MS[31])/Nx[2]
       out[28,3] <- (MS[28] - MS[31])/Nx[3]
       out[27,3] <- (MS[27] - MS[31])/Nx[4]
       out[26,3] <- (MS[26] - MS[31])/Nx[5]
       out[25,3] <- (MS[25] - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - MS[31])/(Nx[1]*Nx[2])
       out[24,3] <- (MS[24] - (Nx[1]*out[30,3]) - (Nx[3]*out[28,3]) - MS[31])/(Nx[1]*Nx[3])
       out[23,3] <- (MS[23] - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - MS[31])/(Nx[2]*Nx[3])
       out[22,3] <- (MS[22] - (Nx[1]*out[30,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[1]*Nx[4])
       out[21,3] <- (MS[21] - (Nx[2]*out[29,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[2]*Nx[4])
       out[20,3] <- (MS[20] - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[3]*Nx[4])
       out[19,3] <- (MS[19] - (Nx[1]*out[30,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[5])
       out[18,3] <- (MS[18] - (Nx[2]*out[29,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[2]*Nx[5])
       out[17,3] <- (MS[17] - (Nx[3]*out[28,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[3]*Nx[5])
       out[16,3] <- (MS[16] - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[4]*Nx[5])
       out[15,3] <- (MS[15] - (Nx[1]*Nx[2]*out[25,3]) - (Nx[1]*Nx[3]*out[24,3]) - (Nx[2]*Nx[3]*out[23,3]) - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - MS[31])/(Nx[1]*Nx[2]*Nx[3])
       out[14,3] <- (MS[14] - (Nx[1]*Nx[2]*out[25,3]) - (Nx[1]*Nx[4]*out[22,3]) - (Nx[2]*Nx[4]*out[21,3]) - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[1]*Nx[2]*Nx[4])
       out[13,3] <- (MS[13] - (Nx[1]*Nx[3]*out[24,3]) - (Nx[1]*Nx[4]*out[22,3]) - (Nx[3]*Nx[4]*out[20,3]) - (Nx[1]*out[30,3]) - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[1]*Nx[3]*Nx[4])
       out[12,3] <- (MS[12] - (Nx[2]*Nx[3]*out[23,3]) - (Nx[2]*Nx[4]*out[21,3]) - (Nx[3]*Nx[4]*out[20,3]) - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[2]*Nx[3]*Nx[4])
       out[11,3] <- (MS[11] - (Nx[1]*Nx[2]*out[25,3]) - (Nx[1]*Nx[5]*out[19,3]) - (Nx[2]*Nx[5]*out[18,3]) - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[2]*Nx[5])
       out[10,3] <- (MS[10] - (Nx[1]*Nx[3]*out[24,3]) - (Nx[1]*Nx[5]*out[19,3]) - (Nx[3]*Nx[5]*out[17,3]) - (Nx[1]*out[30,3]) - (Nx[3]*out[28,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[3]*Nx[5])
       out[9,3] <- (MS[9] - (Nx[2]*Nx[3]*out[23,3]) - (Nx[2]*Nx[5]*out[18,3]) - (Nx[3]*Nx[5]*out[17,3]) - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[2]*Nx[3]*Nx[5])
       out[8,3] <- (MS[8] - (Nx[1]*Nx[4]*out[22,3]) - (Nx[1]*Nx[5]*out[19,3]) - (Nx[4]*Nx[5]*out[16,3]) - (Nx[1]*out[30,3]) - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[4]*Nx[5])
       out[7,3] <- (MS[7] - (Nx[2]*Nx[4]*out[21,3]) - (Nx[2]*Nx[5]*out[18,3]) - (Nx[4]*Nx[5]*out[16,3]) - (Nx[2]*out[29,3]) - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[2]*Nx[4]*Nx[5])
       out[6,3] <- (MS[6] - (Nx[3]*Nx[4]*out[20,3]) - (Nx[3]*Nx[5]*out[17,3]) - (Nx[4]*Nx[5]*out[16,3]) - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[3]*Nx[4]*Nx[5])

       out[5,3] <- (MS[5] - (Nx[1]*Nx[2]*Nx[3]*out[15,3]) - (Nx[1]*Nx[2]*Nx[4]*out[14,3]) - (Nx[1]*Nx[3]*Nx[4]*out[13,3]) - (Nx[2]*Nx[3]*Nx[4]*out[12,3]) - (Nx[1]*Nx[2]*out[25,3]) - (Nx[1]*Nx[3]*out[24,3]) - (Nx[2]*Nx[3]*out[23,3]) - (Nx[1]*Nx[4]*out[22,3]) - (Nx[2]*Nx[4]*out[21,3]) - (Nx[3]*Nx[4]*out[20,3]) - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - MS[31])/(Nx[1]*Nx[2]*Nx[3]*Nx[4])

       out[4,3] <- (MS[4] - (Nx[1]*Nx[2]*Nx[3]*out[15,3]) - (Nx[1]*Nx[2]*Nx[5]*out[11,3]) - (Nx[1]*Nx[3]*Nx[5]*out[10,3]) - (Nx[2]*Nx[3]*Nx[5]*out[9,3]) - (Nx[1]*Nx[2]*out[25,3]) - (Nx[1]*Nx[3]*out[24,3]) - (Nx[2]*Nx[3]*out[23,3]) - (Nx[1]*Nx[5]*out[19,3]) - (Nx[2]*Nx[5]*out[18,3]) - (Nx[3]*Nx[5]*out[17,3]) - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[2]*Nx[3]*Nx[5])

       out[3,3] <- (MS[3] - (Nx[1]*Nx[2]*Nx[4]*out[14,3]) - (Nx[1]*Nx[2]*Nx[5]*out[11,3]) - (Nx[1]*Nx[4]*Nx[5]*out[8,3]) - (Nx[2]*Nx[4]*Nx[5]*out[7,3]) - (Nx[1]*Nx[2]*out[25,3]) - (Nx[1]*Nx[4]*out[22,3]) - (Nx[2]*Nx[4]*out[21,3]) - (Nx[1]*Nx[5]*out[19,3]) - (Nx[2]*Nx[5]*out[18,3]) - (Nx[4]*Nx[5]*out[16,3]) - (Nx[1]*out[30,3]) - (Nx[2]*out[29,3]) - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[2]*Nx[4]*Nx[5])

       out[2,3] <- (MS[2] - (Nx[1]*Nx[3]*Nx[4]*out[13,3]) - (Nx[1]*Nx[3]*Nx[5]*out[10,3]) - (Nx[1]*Nx[4]*Nx[5]*out[8,3]) - (Nx[3]*Nx[4]*Nx[5]*out[6,3]) - (Nx[1]*Nx[3]*out[24,3]) - (Nx[1]*Nx[4]*out[22,3]) - (Nx[3]*Nx[4]*out[20,3]) - (Nx[1]*Nx[5]*out[19,3]) - (Nx[3]*Nx[5]*out[17,3]) - (Nx[4]*Nx[5]*out[16,3]) - (Nx[1]*out[30,3]) - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[1]*Nx[3]*Nx[4]*Nx[5])

       out[1,3] <- (MS[1] - (Nx[2]*Nx[3]*Nx[4]*out[12,3]) - (Nx[2]*Nx[3]*Nx[5]*out[9,3]) - (Nx[2]*Nx[4]*Nx[5]*out[7,3]) - (Nx[3]*Nx[4]*Nx[5]*out[6,3]) - (Nx[1]*Nx[3]*out[23,3]) - (Nx[2]*Nx[4]*out[21,3]) - (Nx[3]*Nx[4]*out[20,3]) - (Nx[2]*Nx[5]*out[18,3]) - (Nx[3]*Nx[5]*out[17,3]) - (Nx[4]*Nx[5]*out[16,3]) - (Nx[2]*out[29,3]) - (Nx[3]*out[28,3]) - (Nx[4]*out[27,3]) - (Nx[5]*out[26,3]) - MS[31])/(Nx[2]*Nx[3]*Nx[4]*Nx[5])
       out[out[,3] < 0,3] <- 0
       out[,4] <- round(out[,3]/sum(out[,3]),digits=2)*100
      }
   return(out)
   }



  if(length(Effects) != ExpParms[Nparms-1]){
    print("Error: Unbalanced Model")
  } else {
  res <- calcGT(x,Effects,DF,MS,Nparms)
  return(round(res,digits=2))
  }
}


Dcoef <- function(x,gt,evar=.01){
  Nparms <- length(all.vars(formula(x))) - 1
  xlevels <- x$xlevels
  if (Nparms == 2){
    xlvl1 <- as.numeric(xlevels[[1]])
    xlvl2 <- as.numeric(xlevels[[2]])
    rVarmax <- round(gt[3,3]/evar)
    aVarmax <- round((gt[2,3] + gt[3,3])/evar)

    newlvls <- c(min(xlvl2),round(mean(xlvl2)),max(xlvl2),2*max(xlvl2),rVarmax,aVarmax)
    relVar <- gt[3,3] / newlvls
    absVar <- (gt[2,3] + gt[3,3])/newlvls
    gCoef <- gt[1,3]/(gt[1,3] + relVar)
    phiCoef <- gt[1,3]/(gt[1,3] + absVar)
    #rho <- gt[1,3]/((gt[1,3] + gt[3,3])/newlvls) ### values do not look right with this formula CHECK
    
    out <- round(data.frame(cbind(newlvls,relVar,absVar,gCoef,phiCoef)),2)
    names(out)[1] <- paste(rownames(gt)[2],"levels",sep=".")
    
    }

  if (Nparms == 3){
    xlvl1 <- as.numeric(xlevels[[1]])
    xlvl2 <- as.numeric(xlevels[[2]])
    xlvl3 <- as.numeric(xlevels[[3]])
    newlvls2 <- c(min(xlvl2,na.rm=T),round(mean(xlvl2)),max(xlvl2,na.rm=T),2*max(xlvl2,na.rm=T),100,500)
    newlvls3 <- c(min(xlvl3,na.rm=T),round(mean(xlvl3,na.rm=T)),max(xlvl3,na.rm=T),2*max(xlvl3,na.rm=T),5*max(xlvl3,na.rm=T),10*max(xlvl3,na.rm=T))
    levels <- expand.grid(newlvls2,newlvls3)
    out <- data.frame(iter=1:length(newlvls2)*length(newlvls3),levels,relVar=0,absVar=0,gCoef=0,phi=0)

    for(i in 1:nrow(out)){
        out$relVar[i] <- (gt[4,3]/out[i,2]) + (gt[5,3]/out[i,3]) + (gt[7,3]/(out[i,2]*out[i,3]))
        out$absVar[i] <- (gt[2,3]/out[i,2]) + (gt[3,3]/out[i,3]) + (gt[6,3]/(out[i,2]*out[i,3])) + (gt[4,3]/out[i,2]) + (gt[5,3]/out[i,3]) + (gt[7,3]/(out[i,2]*out[i,3]))
    }
    out$gCoef <- gt[1,3]/(gt[1,3] + out$relVar)
    out$phi <- gt[1,3]/(gt[1,3] + out$absVar)
    names(out)[2:3] <- c(paste(rownames(gt)[2],"levels",sep="."),paste(rownames(gt)[3],"levels",sep="."))
    out <- round(out[,-1],2)
  }
  return(out)
}

plot.Dcoef <- function(x){}



#################################################################
### igc - compute individual growth curve parameters using lm
### added:  2/10/2007
###
### usage igc(dat,idvar="",ivar="",dvar="",parms=2 or 3

igc <- function(x,...){
  UseMethod("igc")
}

igc.default <- function(x,idvar,ivar,dvar,cvar=NULL,parms=2,method="OLS"){

  # First index the data file (x) with the names

  ID <- match(idvar,names(x))
  IV <- match(ivar,names(x))
  DV <- match(dvar,names(x))
  CV <- match(cvar,names(x))

  # subset data so we have only the three relevant variables
  x <- data.frame(ID=x[,ID],IV=x[,IV],DV=x[,DV])
  rm(ID,IV,DV) # get rid of these values because they mess up things below

 
  #### Missing Data Handling
  
  # Second get rid of NA's in the dv vector
  x <- x[!is.na(x$DV),]

  # Third, make sure we have more than 3 observations per subject and only retain those that do have more than 3
  x.sums <- data.frame(ID=row.names(table(x$ID,x$IV)),counts=rowSums(table(x$ID,x$IV)))
  x.lim <- x.sums[x.sums[,2] > 2,]
  x <- x[x$ID %in% x.lim[,1],]
  IDuniq <- unique(x$ID)


  
  # Fourth, store some useful values for later analysis
  if (min(x$IV,na.rm=T) < 0){
    xlim <- c(min(x$IV,na.rm=T),max(x$IV,na.rm=T))
  } else {
    xlim <- c(0,max(x$IV,na.rm=T))
  }

  if (min(x$DV,na.rm=T) < 0){
    ylim <- c(min(x$DV,na.rm=T),max(x$DV,na.rm=T))
  } else {
    ylim <- c(0,max(x$DV,na.rm=T))
  }
  

  ### Now get the parameter estimates

  if (parms == 2){
    if (method == "OLS"){  # Compute the IGC's using OLS via lm
      gc.out <- data.frame(id=IDuniq,Intparm=0,Lparm=0,IntSE=0,LparmSE=0,Rsq=0)
      for (i in 1:length(IDuniq)){
        dat <- subset(x,x$ID == IDuniq[i])
        lm.tmp <- lm(DV~IV,data=dat)
        gc.out[i,2] <- coef(lm.tmp)[[1]]
        gc.out[i,3] <- coef(lm.tmp)[[2]]
        gc.out[i,4] <- summary(lm.tmp)$coefficients[3]
        gc.out[i,5] <- summary(lm.tmp)$coefficients[4]
        gc.out[i,6] <- summary(lm.tmp)[[9]]
      }
    }
    if (method == "ML"){ # Compute the IGC's using ML via lmer
      library(lme4)
      lme.tmp <- lmer(DV~IV + (IV|ID),data=x,na.action=na.exclude)
      gc.out <- data.frame(id=as.numeric(as.character(rownames(coef(lme.tmp)[[1]]))),Intparm=coef(lme.tmp)[[1]][,1],Lparm=coef(lme.tmp)[[1]][,2])
      }
    # compute fixed parameters for graphiing the results
    lm.fixed <- lm(DV~IV,data=x)
    fixed.parms <- c(coef(lm.fixed)[[1]],coef(lm.fixed)[[2]],summary(lm.fixed)$coefficients[3],summary(lm.fixed)$coefficients[4])
    }

  if (parms == 3){
    if (method == "OLS"){
      gc.out <- data.frame(id=IDuniq,Intparm=0,Lparm=0,Qparm=0,IntSE=0,LparmSE=0,QparmSE=0,Rsq=0)
      for (i in 1:length(IDuniq)){
        dat <- subset(x,x$ID == IDuniq[i])
        lm.tmp <- lm(DV~IV+I(IV*IV),data=dat)
        gc.out[i,2] <- coef(lm.tmp)[[1]]
        gc.out[i,3] <- coef(lm.tmp)[[2]]
        gc.out[i,4] <- coef(lm.tmp)[[3]]
        gc.out[i,5] <- summary(lm.tmp)$coefficients[4]
        gc.out[i,6] <- summary(lm.tmp)$coefficients[5]
        gc.out[i,7] <- summary(lm.tmp)$coefficients[6]        
        gc.out[i,8] <- summary(lm.tmp)[[9]]
      }
    }
    if (method == "ML"){
      library(nlme)
      dat.grpd <- groupedData(DV~IV | ID, data=x)
      lme.tmp <- lme(DV~IV + I(IV*IV),data=dat.grpd)
      gc.out <- data.frame(id=rownames(coef(lme.tmp)),Intparm=coef(lme.tmp)[,1],Lparm=coef(lme.tmp)[,2],Qparm=coef(lme.tmp)[,3])
    }
    # compute fixed parameters for graphing the results
    lm.fixed <- lm(DV~IV+I(IV*IV),data=x)
    fixed.parms <- c(coef(lm.fixed)[[1]],coef(lm.fixed)[[2]],coef(lm.fixed)[[3]],summary(lm.fixed)$coefficients[4],summary(lm.fixed)$coefficients[5],summary(lm.fixed)$coefficients[6])
  }
  
  res <- list(params=gc.out,parms=parms,method=method,xlim=xlim,ylim=ylim,fixed.parms=fixed.parms)
  class(res) <- "igc"
  return(res)
}


plot.igc <- function(x,xlab="",ylab="",main="",...){
#  if(length(ylim)==0){
    ylim <- x$ylim
#  }
  xlim <- x$xlim
  gcdat <- x$params
  fixed <- x$fixed.parms
  if ((x$method == "OLS" & ncol(x$params) == 6) | (x$method == "ML" & ncol(x$params) == 3)){
    curve(gcdat[1,3]*x + gcdat[1,2],min(xlim),max(xlim),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
    for (i in 2:nrow(gcdat)){
      curve(gcdat[i,3]*x + gcdat[i,2],min(xlim),max(xlim),add=T)
    }
    curve(fixed[2]*x + fixed[1],min(xlim),max(xlim),lwd=2,col="red",add=T)
  }
  if ((x$method == "OLS" & ncol(x$params) == 7) | (x$method == "ML" & ncol(x$params) == 4)){
    curve(gcdat[1,3]*x + gcdat[1,4]*x^2 + gcdat[1,2],min(xlim),max(xlim),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
    for (i in 2:nrow(gcdat)){
      curve(gcdat[i,3]*x + gcdat[i,4]*x^2 + gcdat[i,2],min(xlim),max(xlim),lwd=2,col="red",add=T)
    }
  }
}


summary.igc <- function(x){
  if (length(x) == 6){
    out <- matrix(0,3,2)
    out[1,1] <- round(mean(x[[2]],na.rm=T),2)
    out[2,1] <- round(mean(x[[3]],na.rm=T),2)
    out[3,1] <- round(mean(x[[4]],na.rm=T),2)
    out[1,2] <- round(sd(x[[2]],na.rm=T),2)
    out[2,2] <- round(sd(x[[3]],na.rm=T),2)
    out[3,2] <- round(sd(x[[4]],na.rm=T),2)    
    rownames(out) <- c("Intercept","Linear.Slope","R-squared")
    colnames(out) <- c("Mean","SD")
  }

  if (length(x) == 7){
    out <- matrix(0,4,2)
    out[1,1] <- round(mean(x[[2]],na.rm=T),2)
    out[2,1] <- round(mean(x[[3]],na.rm=T),2)
    out[3,1] <- round(mean(x[[4]],na.rm=T),2)
    out[4,1] <- round(mean(x[[5]],na.rm=T),2)    
    out[1,2] <- round(sd(x[[2]],na.rm=T),2)
    out[2,2] <- round(sd(x[[3]],na.rm=T),2)
    out[3,2] <- round(sd(x[[4]],na.rm=T),2)
    out[4,2] <- round(sd(x[[5]],na.rm=T),2)
    rownames(out) <- c("Intercept","Linear.Slope","Quadratic.Slope","R-squared")
    colnames(out) <- c("Mean","SD")
  }
  return(out)
}

coef.igc <- function(x,prefix){
  dat <- x$params
  if (x$parms == 2){
    if (x$method == "OLS"){
      names(dat) <- c("id",paste(prefix,"I",sep=""),paste(prefix,"L",sep=""),paste(prefix,"Ise",sep=""),paste(prefix,"Lse",sep=""),paste(prefix,"Rsq",sep=""))
    }
    if (x$method == "ML"){
      names(dat) <- c("id",paste(prefix,"I",sep=""),paste(prefix,"L",sep=""))
    }
  }
  if (x$parms == 3){
   if (x$method == "OLS"){
      names(dat) <- c("id",paste(prefix,"I",sep=""),paste(prefix,"L",sep=""),paste(prefix,"Q",sep=""),paste(prefix,"Ise",sep=""),paste(prefix,"Lse",sep=""),paste(prefix,paste(prefix,"Qse",sep=""),"Rsq",sep=""))
    }
    if (x$method == "ML"){
      names(dat) <- c("id",paste(prefix,"I",sep=""),paste(prefix,"L",sep=""),paste(prefix,"Q",sep=""))
    }
  } 
  return(dat)
}



  
##############################
##  corSum function
##  added:  11/1/2006
##  latest edited: 11/1/2006

corSum <- function(x,...){
  # in memory of Jeff and his jeffs.cor function
  UseMethod("corSum")
}

corSum.default <- function(x){
#  x <- x[,sapply(x,is.numeric)]
  out <- round(cor(x),digits=2)
  diag(out) <- 0
  out.tab <- as.data.frame.table(out)
  out.tab <- out.tab[c(lower.tri(out)),]
  out.tab$abscor <- abs(out.tab[,3])
  out.tab <- out.tab[order(out.tab$abscor,decreasing=T),]
  names(out.tab)[3:4] <- c("Cor","AbsCor")
  res <- list(out.tab)
  class(res) <- "corSum"
  invisible(res)
}

summary.corSum <- function(object){
  z <- object[[1]]
  mincor <- min(z[,4])
  maxcor <- max(z[,4])
  meancor <- mean(z[,3])
  sdcor <- sd(z[,3])
  vars.min <- z[z[,4]==mincor,1:2]
  vars.max <- z[z[,4]==maxcor,1:2]
  cat("\n Summary of Correlation Matrix \n\n")
  cat("Minimum Correlation:")
  print(mincor)
  cat("\n Maximum Correlation:")
  print(maxcor)
  cat("\n Mean Correlation:")
  print(meancor)
  cat("\n SD of Correlations:")
  print(sdcor)
  cat("\n\n -------------------------- \n")
  cat("\n Variables that have the minimum correlation \n")
  print(vars.min)
  cat("\n Variables that have the maximum correlation \n")
  print(vars.max)
}

plot.corSum <- function(x,...){
  hist(x[[1]][,3],xlab="Observed Correlations",main="")
}

