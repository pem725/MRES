sbl <- function(pred,med,out,time,group,id,boot=1000){
  UseMethod("sbl")
}

sbl.default <- function(pred,med,out,time=NULL,group=NULL,id,boot=1000){
  dat <- data.frame(IV=pred,MED=med,DV=out,time=time,group=group,id=id)
  library(lme4)
  library(nlme)

  if (is.null(time)){
    if (is.null(group)){
      m1.lme <- summary(lme(DV~IV,random=~1|id,data=dat))$tTable
      m2.lme <- summary(lme(DV~MED+IV,random=~1|id,data=dat))$tTable
      m3.lme <- summary(lme(MED~IV,random=~1|id,data=dat))$tTable
    } else {
      m1.lme <- summary(lme(DV~time+group+IV,random=~time+1|id,data=dat))$tTable
      m2.lme <- summary(lme(DV~time+group+MED+IV,random=~time+1|id,data=dat))$tTable
      m3.lme <- summary(lme(MED~time+group+IV,random=~time+1|id,data=dat))$tTable
    }
  } else {
    if (is.null(group)){
      m1.lme <- summary(lme(DV~time+IV,random=~time+1|id,data=dat))$tTable
      m2.lme <- summary(lme(DV~time+MED+IV,random=~time+1|id,data=dat))$tTable
      m3.lme <- summary(lme(MED~time+IV,random=~time+1|id,data=dat))$tTable
    } else {
      m1.lme <- summary(lme(DV~time+group+IV,random=~time+1|id,data=dat))$tTable
      m2.lme <- summary(lme(DV~time+group+MED+IV,random=~time+1|id,data=dat))$tTable
      m3.lme <- summary(lme(MED~time+group+IV,random=~time+1|id,data=dat))$tTable
    }
            
  indir <- m3.lme[5,1]*m2.lme[5,1]
  effvar <- (m3.lme[5,1])^2 * (m2.lme[5,2])^2 + (m2.lme[5,1])^2 * (m3.lme[5,2])^2
  serr <- sqrt(effvar)
  zvalue <- indir/serr

  c.sig <- 0
  c.sig[m1.lme[5,5] < .05] <- 1
  a.sig <- 0
  a.sig[m3.lme[5,5] < .05] <- 1
  b.sig <- 0
  b.sig[m2.lme[5,5] < .05] <- 1
  cprime.sig <- 0
  cprime.sig[m2.lme[6,5] < .05] <- 1
  suppression <- 0
  suppression[abs(m2.lme[6,1]) > abs(m1.lme[5,1])] <- 1

  mediation.test <- list(c.sig=c.sig,a.sig=a.sig,b.sig=b.sig,cprime.sig=cprime.sig,suppression=suppression)

  calcMED <- function(x){
    dat <- x
    m1 <- lmer(DV~time+group+IV+(time-1|id)+(1|id),data=dat)
    m2 <- lmer(DV~time+group+MED+IV+(time-1|id)+(1|id),data=dat)
    m3 <- lmer(MED~time+group+IV+(time-1|id)+(1|id),data=dat)
    a <- summary(m3)@coefs[5,1]
    b <- summary(m2)@coefs[5,1]
    c <- summary(m1)@coefs[5,1]
    c.prime <- summary(m2)@coefs[6,1]
    se.a <- summary(m3)@coefs[5,2]
    se.b <- summary(m2)@coefs[5,2]
    se.c <- summary(m1)@coefs[5,2]
    se.c.prime <- summary(m2)@coefs[6,2]
    ind.eff <- a*b
    se.ind.eff <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
    zvalue <- ind.eff/se.ind.eff
    out <- c(a,b,c,c.prime,se.a,se.b,se.c,se.c.prime,ind.eff,se.ind.eff,zvalue)
    return(out)
  }
  sobel.out <- calcMED(dat)
  outboot <- matrix(0,boot,11)
  outboot <- data.frame(outboot)  
  for(i in 1:boot){
    dat.new <- dat[sample(1:nrow(dat),nrow(dat),replace=TRUE),]
    outboot[i,] <- calcMED(dat.new)
  }
  names(outboot) <- c("a","b","c","c.prime","se.a","se.b","se.c","se.c.prime","ind.eff","se.ind.eff","zvalue")

  boot.se <- sd(outboot$ind.eff)
  boot.mean <- mean(outboot$ind.eff)

  UL95 <- quantile(outboot$ind.eff,.95)
  LL95 <- quantile(outboot$ind.eff,.05)
  UL99 <- quantile(outboot$ind.eff,.99)
  LL99 <- quantile(outboot$ind.eff,.01)

  conf.limits <- c(LL95,UL95,LL99,UL99)

  res <- list(m1=m1.lme,m2=m2.lme,m3=m3.lme,mediation.test=mediation.test,fixed=sobel.out,boot.results=outboot,boot.conf=conf.limits,sobel.effect=sobel.out[9],sobel.se=sobel.out[10],boot.effect=boot.mean,boot.se=boot.se)
  class(res) <- "sbl"
  return(res)
}


summary.sbl <- function(x){
  cat("Mediation Results\n")
  cat("--------------------------------\n")
  cat("Model 1\n")
  print(round(x$m1,3))
  cat("--------------------------------\n")
  cat("Model 2\n")
  print(round(x$m2,3))
  cat("--------------------------------\n")
  cat("Model 3\n")
  print(round(x$m3,3))
  cat("---------------------------------------------------\n\n")


  cat("Standard Models\n")
  
  cat("    IV -- c --> DV\n\n")

  cat("          M           \n")
  cat("         / \\         \n")
  cat("      a /   \\ b      \n")
  cat("       /     \\       \n")
  cat("      /       \\      \n")
  cat("     IV - c' -> DV    \n")

  cat("--------------------------------\n")
  
  cat("Mediation Diagnostics:\n")
  cat("Path c significant: ")
  if(x$mediation.test[[1]] == 1) cat("YES\n")  else cat("NO\n")
  cat("Path a significant: ")
  if(x$mediation.test[[2]] == 1) cat("YES\n")  else cat("NO\n")
  cat("Path b significant: ")
  if(x$mediation.test[[3]] == 1) cat("YES\n")  else cat("NO\n")
  cat("Path c' significant: ")
  if(x$mediation.test[[4]] == 1) cat("YES\n")  else cat("NO\n")
  cat("Evidence of suppression: ")
  if(x$mediation.test[[5]] == 1) cat("YES\n")  else cat("NO\n")

  cat("--------------------------------\n")
  cat("Mean Results from Bootstrap:\n")
  cat("a: ")
  print(round(mean(x$boot.results$a,na.rm=T),2))
  cat(" (")
  print(round(mean(x$boot.results$se.a,na.rm=T),2))
  cat(")\n")
  cat("b: ")
  print(round(mean(x$boot.results$b,na.rm=T),2))
  cat(" (")
  print(round(mean(x$boot.results$se.b,na.rm=T),2))
  cat(")\n")
  cat("c: ")
  print(round(mean(x$boot.results$c,na.rm=T),2))
  cat(" (")
  print(round(mean(x$boot.results$se.a,na.rm=T),2))
  cat(")\n")
  cat("c prime: ")
  print(round(mean(x$boot.results$c.prime,na.rm=T),2))
  cat(" (")
  print(round(mean(x$boot.results$se.c.prime,na.rm=T),2))
  cat(")\n")
  cat("ab: ")
  print(round(mean(x$boot.results$ind.eff,na.rm=T),2))
  cat(" (")
  print(round(mean(x$boot.results$se.ind.eff,na.rm=T),2))
  cat(")\n")
}


plot.sbl <- function(x,main="",col="lightyellow",border="gray",in.color=T,...){
  if(in.color==T){
    hist(x$boot.results$ind.eff,col=col,border=border,xlab="Estimated Indirect Effect",main=main)
    abline(v=x$boot.effect,col="red",lwd=2)
    abline(v=x$boot.conf[[1]],col="red",lty=2,lwd=2)
    abline(v=x$boot.conf[[2]],col="red",lty=2,lwd=2)
    abline(v=x$sobel.effect,col="blue",lwd=2)
    legend("topleft",c("Bootstrap Mean","95% Confidence Limits","Sobel Estimate"),col=c("red","red","blue"),lty=c(1,2,1),lwd=2)
  } else {
    hist(x$boot.results$ind.eff,col="white",border="black",xlab="Estimated Indirect Effect",main=main)
    abline(v=x$boot.effect,lwd=3)
    abline(v=x$boot.conf[[1]],lty=2,lwd=3)
    abline(v=x$boot.conf[[2]],lty=2,lwd=3)
    abline(v=x$sobel.effect,lwd=2,lty=3)
    legend("topleft",c("Bootstrap Mean","95% Confidence Limits","Sobel Estimate"),col="black",lty=c(1,2,3),lwd=c(3,3,2))
  } 
}
  
qq.sbl <- function(x,main="",...){
  qqnorm(x$boot.results$ind.eff)
  qqline(x$boot.results$ind.eff)
}

## drawfig.sbl <- function(x,iv="IV",med="MED",dv="DV"){
##   library(RGraphviz)
  
  
