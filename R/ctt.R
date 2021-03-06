ctt <- function(x,...){
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
  ## out <- aggregate(scores,list(x[,1]),mean,na.rm=T)
  ## for (i in 2:ncol(x)){
  ##   out <- merge(out,aggregate(scores,list(x[,i]),mean,na.rm=T),by="Group.1",all=T)
  ## }

  for (i in 1:ncol(x)){
    assign(paste("var",i,sep=""),aggregate(scores,list(x[,i]),mean,na.rm=T))
    tmp <- get(paste("var",i,sep=""))
    names(tmp) <- c("RespOpt",paste("MeanItem",i,sep=""))
    assign(paste("var",i,sep=""),tmp)
  }

  out <- get(paste("var",1,sep=""))
  for (i in 2:ncol(x)){
    out <- merge(out,get(paste("var",i,sep="")),by="RespOpt",all=T)
  }

  # garbage cleanup
  rm(list=paste("var",1:ncol(x),sep=""))
  
  out <- t(out)[-1,]
  row.names(out) <- names(x)
  out <- as.data.frame(out)
  return(out)
}

resp.cat.usage <- function(x){
  initial <- data.frame(resp=unique(as.vector(as.matrix(x)))[order(unique(as.vector(as.matrix(x))))],dummy=NA) 
  
  out <- as.data.frame(table(x[,1],useNA="ifany"))
  names(out) <- c("resp",names(x)[1])
  out <- merge(initial,out,by="resp",all=T)[,-2]
  
  for (i in 2:ncol(x)){
    out <- merge(out,as.data.frame(table(x[,i])),by.x="resp",by.y="Var1")
  }
  names(out)[-1] <- names(x)
  out <- t(out)
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

