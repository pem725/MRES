##############################################################################
### Compute generalizability theory coefficients with the following
### code: x is an aov output object and the ANOVA must be a
### fully-crossed factorial design

gtheory <- function(x,data,gcrit=.8,phicrit=.8){
  UseMethod("gtheory")
}

Gstudy <- function(x){
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
  return(res)
  }
}


Dstudy <- function(aout,gt,gcrit=gcrit,phicrit=phicrit){
  Nparms <- length(all.vars(formula(aout))) - 1
  xlevels <- aout$xlevels
  if (Nparms == 2){
    a <- gt[1,3]
    b <- gt[2,3]
    c <- gt[3,3]
    xlvl1 <- 1:length(xlevels[[1]])
    xlvl2 <- 1:length(xlevels[[2]])

    targetg <- round(gcrit*c/(a - gcrit*a))
    targetphi <- round((phicrit*(b+c))/(a - phicrit*a))

##     NO LONGER WANT TO CONDITION ON EVAR
##     rVarmax <- round(gt[3,3]/evar)  
##     aVarmax <- round((gt[2,3] + gt[3,3])/evar)

    newlvls <- c(min(xlvl2),round(mean(xlvl2)),max(xlvl2),2*max(xlvl2),targetg,targetphi)
    relVar <- c / newlvls
    absVar <- (b + c)/newlvls
    gCoef <- a/(a + relVar)
    phiCoef <- a/(a + absVar)
    #rho <- gt[1,3]/((gt[1,3] + gt[3,3])/newlvls) ### values do not look right with this formula CHECK
    
    out <- round(data.frame(cbind(newlvls,relVar,absVar,gCoef,phiCoef)),2)
    names(out)[1] <- paste(rownames(gt)[2],"levels",sep=".")
    out <- out[is.finite(out[,1]),]
  }

  if (Nparms == 3){
    a <- gt[1,3]
    b <- gt[2,3]
    c <- gt[3,3]
    d <- gt[4,3]
    e <- gt[5,3]
    f <- gt[6,3]
    g <- gt[7,3]

    xlvl1 <- 1:length(xlevels[[1]])
    xlvl2 <- 1:length(xlevels[[2]])
    xlvl3 <- 1:length(xlevels[[3]])
    bN <- max(xlvl2)
    cN <- max(xlvl3)

    ## simplified for n by hand and produced the following equations
    targetg1 <- round((gcrit*((d*cN) + g))/((a*cN) - gcrit*((a*cN) - e)))
    targetg2 <- round((gcrit*(g + (bN*e)))/((a*bN)- gcrit*((a*bN)-d)))

    targetphi1 <- round((phicrit*((b*cN)+(d*cN)+f+g))/((a*cN) - (phicrit*((a*cN) - c - e))))
    targetphi2 <- round((phicrit*((c*bN) + (e*bN) + f + g))/((a*bN) - b - d))
    
    newlvls2 <- c(min(xlvl2,na.rm=T),round(mean(xlvl2)),max(xlvl2,na.rm=T),2*max(xlvl2,na.rm=T),targetg1,targetphi1)
    newlvls3 <- c(min(xlvl3,na.rm=T),round(mean(xlvl3,na.rm=T)),max(xlvl3,na.rm=T),2*max(xlvl3,na.rm=T),targetg2,targetphi2)
    levels <- expand.grid(newlvls2,newlvls3)
    out <- data.frame(iter=1:length(newlvls2)*length(newlvls3),levels,relVar=0,absVar=0,gCoef=0,phi=0)

    for(i in 1:nrow(out)){
        out$relVar[i] <- (gt[4,3]/out[i,2]) + (gt[5,3]/out[i,3]) + (gt[7,3]/(out[i,2]*out[i,3]))
        out$absVar[i] <- (gt[2,3]/out[i,2]) + (gt[3,3]/out[i,3]) + (gt[6,3]/(out[i,2]*out[i,3])) + (gt[4,3]/out[i,2]) + (gt[5,3]/out[i,3]) + (gt[7,3]/(out[i,2]*out[i,3]))
    }
    out$gCoef <- gt[1,3]/(gt[1,3] + out$relVar)
    out$phi <- gt[1,3]/(gt[1,3] + out$absVar)
    names(out)[2:3] <- c(paste(rownames(gt)[2],"levels",sep="."),paste(rownames(gt)[3],"levels",sep="."))
    out <- out[is.finite(out[,1]),-1]

  }

  ## from Charles C. Berry via R-help posting: used to eliminate redundancies in the output
  
  count.rows <- function(x){
    order.x <- do.call(order,as.data.frame(x))
    equal.to.previous <- rowSums(x[tail(order.x,-1),] != x[head(order.x,-1),])==0
    tf.runs <- rle(equal.to.previous)
    counts <- c(1,unlist(mapply( function(x,y) if (y) x+1 else (rep(1,x)),tf.runs$length, tf.runs$value )))
    counts <- counts[ c(diff(counts) <= 0, TRUE ) ]
    unique.rows <- which( c(TRUE, !equal.to.previous ) )
    cbind( counts, x[order.x[ unique.rows ], ,drop=F] )
  }

  out <- count.rows(out)[,-1]
  row.names(out) <- 1:nrow(out)
  return(out)
}


gtheory.default <- function(formula,data=NULL,gcrit=.8,phicrit=.8){
  Nparms <- length(all.vars(formula(formula))) - 1
  aout <- aov(formula,data)
  gs <- Gstudy(aout)
  ds <- Dstudy(aout,gs,gcrit=gcrit,phicrit=phicrit)
  res <- list(aout=aout,gs=gs,gcrit=gcrit,phicrit=phicrit,Nparms=Nparms,ds=ds) ## ds=ds, ## note I pulled out the Dstudy stuff until I fix it
  class(res) <- "gtheory"
  invisible(res)
  return(res)
}

summary.gtheory <- function(x){
  aout <- x$aout
  gs <- x$gs
#  ds <- x$ds
  cat("\n Generalizability Theory Results: \n\n")
  cat("\n ANOVA Table Summary: \n\n")
  print(summary(aout),digits=2)
  cat("\n G-Study Results: \n\n")
  gmat <- gs[,3]/(gs[,3] + ((gs[nrow(gs),3])/length(aout$xlevels[[1]])))
  gsfin <- cbind(gs,gmat)
  names(gsfin) <- c("df","MS","Var Comp","Percent Var","G Coef")
  print(round(gsfin,2))
#  cat("\n D-Study Results: \n\n")
#  print(round(ds,2))
}
 

plot.gtheory <- function(x,coef="g",legend=T,ablcol=1,xlab=NULL,ltitle=NULL,legposx="topleft",legposy=NULL,...){
  aout <- x$aout
  gout <- x$gs
  dout <- x$ds
  Nparms <- x$Nparms

  if(Nparms == 2){

    a <- gout[1,3]
    b <- gout[2,3]
    c <- gout[3,3]
    
    if(coef=="g"){
      xlab <- switch(1+is.null(xlab),xlab,paste("Levels of",terms(aout)[2][[3]]))
      curve(a/(a + (c/x)),min(dout[,1]),max(dout[,1])*1.5,xlab=xlab,ylab="G-Coefficient Value",ylim=c(0,1))
      abline(x$gcrit,0,lwd=2,lty=2,col=ablcol)
    }
    if(coef=="phi"){
      xlab <- switch(1+is.null(xlab),xlab,paste("Levels of",terms(aout)[2][[3]]))
      curve(a/(a + ((b/x) + (c/x))),min(dout[,1]),max(dout[,1])*1.5,xlab=xlab,ylab="Phi-Coefficient Value (Index of Dependability)")
      abline(x$phicrit,0,lwd=2,lty=2,col=ablcol)
    }
  }

  if(Nparms==3){
    a <- gout[1,3]
    b <- gout[2,3]
    c <- gout[3,3]
    d <- gout[4,3]
    e <- gout[5,3]
    f <- gout[6,3]
    g <- gout[7,3]

    if(coef=="g"){
      d2lvls <- sort(unique(dout[,2]))
      xlab <- switch(1+is.null(xlab),xlab,paste("Levels of",terms(aout)[2][[3]]))
      curve(a/(a + (d/x) + (e/d2lvls[1]) + (g/(x*d2lvls[1]))),min(dout[,1]),max(dout[,1]),xlab=xlab,ylab="G-Coefficient Value",col=1,lwd=2,ylim=c(0,1))
      for(i in 2:length(d2lvls)){
        curve(a/(a + (d/x) + (e/d2lvls[i]) + (g/(x * d2lvls[i]))),min(dout[,1]),max(dout[,1]),xlab="",ylab="",add=T,col=i,lwd=2,ylim=c(0,1))
      }
      abline(x$gcrit,0,lwd=2,lty=2,col=ablcol)
      if(legend==T){
        ltitle <- switch(1+is.null(ltitle),ltitle,names(dout)[2])
        legposx <- switch(1+is.null(legposx), legposx, "topleft")
        legposy <- switch(1+is.null(legposy), legposy, NULL)
        legend(x=legposx,y=legposy,lty=1,lwd=2,col=1:length(d2lvls),legend=d2lvls,title=ltitle,bg="white")
      }
    }

    if(coef=="phi"){
      d2lvls <- sort(unique(dout[,2]))
      xlab <- switch(1+is.null(xlab),xlab,paste("Levels of",terms(aout)[2][[3]]))
      curve(a/(a + (b/x) + (c/d2lvls[1]) + (c/d2lvls[1]) + (d/x) + (e/d2lvls[1]) + (f/(x*d2lvls[1])) + (g/(x*d2lvls[1]))),min(dout[,1]),max(dout[,1]),xlab=xlab,ylab="Phi-Coefficient Value (Index of Dependability)",col=1,lwd=2,ylim=c(0,1))
      for(i in 2:length(d2lvls)){
        curve(a/(a + (b/x) + (c/d2lvls[i]) + (c/d2lvls[i]) + (d/x) + (e/d2lvls[i]) + (f/(x*d2lvls[i])) + (g/(x*d2lvls[i]))),min(dout[,1]),max(dout[,1]),add=T,xlab="",ylab="",col=i,ylim=c(0,1))
      }
      abline(x$phicrit,0,lwd=2,lty=2,col=ablcol)
      if(legend==T){
        ltitle <- switch(1+is.null(ltitle),ltitle,names(dout)[2])
        legpos <- switch(1+is.null(legpos),legpos,"topleft")
        legposy <- switch(1+is.null(legposy), legposy, NULL)
        legend(x=legposx,y=legposy,lty=1,lwd=2,col=1:length(d2lvls),legend=d2lvls,title=ltitle,bg="white")
      }
    }
  }
}

xtable.gtheory <- function(x,gout="g",caption=NULL,label=NULL,align=NULL,digits=NULL,display=NULL,...){
  Nparms <- x$Nparms
  if(gout=="g"){
    x <- x$gs
    digits <- c(0,0,2,2,0)
    }
  if(gout=="d"){
    x <- x$ds
    digits <- c(rep(0,Nparms),2,2,2,2)
  }
  class(x) <- c("xtable","data.frame")
  caption(x)
  label(x) <- label
  align(x) <- switch(1+is.null(align),align,c("l",rep("r",ncol(x))))
#  digits(x) <- switch(1+is.null(digits),digits,c(0,0,2,2,0))
  display(x) <- switch(1+is.null(display),display,c("s",rep("f",ncol(x))))
  return(xtable(x,caption=caption,label=label,align=align,digits=digits,display=display))
}



