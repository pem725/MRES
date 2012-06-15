
#################################################################
### igc - compute individual growth curve parameters using lm
### updated:  11/05/2008
###
### usage igc(dat,idvar="",ivar="",dvar="",parms=2 or 3

igc <- function(x,...){
  UseMethod("igc")
}

igc.default <- function(x,idvar,ivar,dvar,byvar=NULL,cvar=NULL,parms=2,method="OLS"){

  # First index the data file (x) with the names

  dat.orig <- x
  idv <- match(idvar,names(x))
  ID <- match(idvar,names(x))
  IV <- match(ivar,names(x))
  DV <- match(dvar,names(x))
  if(!is.null(cvar)){
    CV <- match(cvar,names(x))
  }
  if(!is.null(byvar)){
    BV <- match(byvar,names(x))
  }

  # subset data so we have only the three relevant variables

  if (is.null(cvar) & is.null(byvar)){
    x <- data.frame(ID=x[,ID],IV=x[,IV],DV=x[,DV])
    rm(ID,IV,DV) # get rid of these values because they mess up things below
  } ## else if (is.null(cvar)

 
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

  if(!is.null(cvar)){
    cvar <- dat.orig[,c(idv,CV)]
    cvar <- cvar[cvar[,1] %in% gc.out$id,]
    cvar <- cvar[!duplicated(cvar),2]
  }
  
  res <- list(params=gc.out,parms=parms,method=method,xlim=xlim,ylim=ylim,fixed.parms=fixed.parms,cvar=cvar,dvar=dvar)
  class(res) <- "igc"
  return(res)
}


plot.igc <- function(x,xlab="",ylab="",main="",selines=T,cplot=F...){
  ylim <- x$ylim
  xlim <- x$xlim
  gcdat <- x$params
  fixed <- x$fixed.parms
  cvar <- x$cvar
  if (x$parms == 2){
    curve(gcdat[1,3]*x + gcdat[1,2],min(xlim),max(xlim),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
    for (i in 2:nrow(gcdat)){
      curve(gcdat[i,3]*x + gcdat[i,2],min(xlim),max(xlim),add=T)
    }
    curve(fixed[2]*x + fixed[1],min(xlim),max(xlim),lwd=2,col="red",add=T)
    if (selines==T){
      curve((fixed[2]+1.96*fixed[4])*x + fixed[1]+1.96*fixed[3],min(xlim),max(xlim),lwd=2,lty=2,col="red",add=T)
      curve((fixed[2]-1.96*fixed[4])*x + fixed[1]-1.96*fixed[3],min(xlim),max(xlim),lwd=2,lty=2,col="red",add=T)
    }    
  }
  if (x$parms == 3){
    curve(gcdat[1,3]*x + gcdat[1,4]*x^2 + gcdat[1,2],min(xlim),max(xlim),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
    for (i in 2:nrow(gcdat)){
      curve(gcdat[i,3]*x + gcdat[i,4]*x^2 + gcdat[i,2],min(xlim),max(xlim),lwd=2,add=T)
    }
    curve(fixed[3]*x^2 + fixed[2]*x + fixed[1],min(xlim),max(xlim),lwd=2,col="red",add=T)
    if (selines==T){
      curve((fixed[3]+1.96*fixed[6])*x^2 + (fixed[2]+1.96*fixed[5])*x + fixed[1]+1.96*fixed[4],min(xlim),max(xlim),lwd=2,lty=2,col="red",add=T)
      curve((fixed[3]-1.96*fixed[6])*x^2 + (fixed[2]-1.96*fixed[5])*x + fixed[1]-1.96*fixed[4],min(xlim),max(xlim),lwd=2,lty=2,col="red",add=T)
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

coef.igc <- function(x,prefix=x$dvar){
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

