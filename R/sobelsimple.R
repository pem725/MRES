sobelsimple <- function(pred,med,out,boot=1000){
  UseMethod("sobelsimple")
}

sobelsimple.default <- function(pred,med,out,boot=1000){
  dat <- data.frame(IV=pred,MED=med,DV=out)

  m1.lm <- summary(lm(DV~IV,data=dat))$coefficients
  m2.lm <- summary(lm(DV~MED+IV,data=dat))$coefficients
  m3.lm <- summary(lm(MED~IV,data=dat))$coefficients

  indir <- m3.lm[2,1]*m2.lm[2,1]
  effvar <- (m3.lm[2,1])^2 * (m2.lm[2,2])^2 + (m2.lm[2,1])^2 * (m3.lm[2,2])^2
  serr <- sqrt(effvar)
  zvalue <- indir/serr

  c.sig <- 0
  c.sig[m1.lm[2,4] < .05] <- 1
  a.sig <- 0
  a.sig[m3.lm[2,4] < .05] <- 1
  b.sig <- 0
  b.sig[m2.lm[2,4] < .05] <- 1
  cprime.sig <- 0
  cprime.sig[m2.lm[3,4] < .05] <- 1
  suppression <- 0
  suppression[abs(m2.lm[3,1]) > abs(m1.lm[2,1])] <- 1

  mediation.test <- list(c.sig=c.sig,a.sig=a.sig,b.sig=b.sig,cprime.sig=cprime.sig,suppression=suppression)

  calcMED <- function(x){
    dat <- x
    m1 <- lm(DV~IV,data=dat)
    m2 <- lm(DV~MED+IV,data=dat)
    m3 <- lm(MED~IV,data=dat)
    a <- summary(m3)$coefficients[2,1]
    b <- summary(m2)$coefficients[2,1]
    c <- summary(m1)$coefficients[2,1]
    c.prime <- summary(m2)$coefficients[3,1]
    se.a <- summary(m3)$coefficients[2,2]
    se.b <- summary(m2)$coefficients[2,2]
    se.c <- summary(m1)$coefficients[2,2]
    se.c.prime <- summary(m2)$coefficients[3,2]
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

  res <- list(m1=m1.lm,m2=m2.lm,m3=m3.lm,mediation.test=mediation.test,fixed=sobel.out,boot.results=outboot,boot.conf=conf.limits,sobel.effect=sobel.out[9],sobel.se=sobel.out[10],boot.effect=boot.mean,boot.se=boot.se)
  class(res) <- "sobelsimple"
  return(res)
}


summary.sobelsimple <- function(x){
  cat("Mediation Results\n")
  cat("--------------------------------\n")
  cat("Model 1: c\n")
  print(round(x$m1,3))
  cat("--------------------------------\n")
  cat("Model 2: b and c'\n")
  print(round(x$m2,3))
  cat("--------------------------------\n")
  cat("Model 3: a\n")
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
  cat(paste("a: ",round(mean(x$boot.results$a,na.rm=T),2)," (",round(mean(x$boot.results$se.a,na.rm=T),2),")\n",sep=""))
  cat(paste("b: ",round(mean(x$boot.results$b,na.rm=T),2)," (",round(mean(x$boot.results$se.b,na.rm=T),2),")\n",sep=""))
  cat(paste("c: ",round(mean(x$boot.results$c,na.rm=T),2)," (",round(mean(x$boot.results$se.a,na.rm=T),2),")\n",sep=""))
  cat(paste("c': ",round(mean(x$boot.results$c.prime,na.rm=T),2)," (",round(mean(x$boot.results$se.c.prime,na.rm=T),2),")\n",sep=""))
  cat(paste("ab: ",round(mean(x$boot.results$ind.eff,na.rm=T),2)," (",round(mean(x$boot.results$se.ind.eff,na.rm=T),2),")\n",sep=""))
}


plot.sobelsimple <- function(x,main="",col="lightyellow",border="gray",in.color=T,...){
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
  
qq.sobelsimple <- function(x,main="",...){
  qqnorm(x$boot.results$ind.eff)
  qqline(x$boot.results$ind.eff)
}
  
DrawfigSS <- function(x=NULL,names=NULL,main=""){
  ##require(Rgraphviz)
  require(diagram)
  # x is the sobelsimple object
  # names is a vector of names specifying the names you wish to have in your diagram (MED, IV, DV) in that order.
  openplotmat()
  elpos <- coordinates(c(1,2))
  fromto <- matrix(ncol=2,byrow=T,data=c(1,3,2,1,2,3))
  nr <- nrow(fromto)
  arrpos <- matrix(ncol=2,nrow=nr)
  for (i in 1:nr){
    arrpos[i,] <- straightarrow(to=elpos[fromto[i,2],],from = elpos[fromto[i,1],],lwd=2,arr.pos=.75,arr.length=.5)
  }
  textrect(elpos[1,], 0.1, lab=names[1])
  textrect(elpos[2,], 0.1, lab=names[2])
  textrect(elpos[3,], 0.1, lab=names[3])

  if(is.null(x)){
    plabs <- c("a","b","c")
  } else {
    parms <- c(round(mean(x$boot.results$a),2),round(mean(x$boot.results$b),2),round(mean(x$boot.results$c.prime),2))
    plabs <- parms
    plabsse <- c(round(mean(x$boot.results$se.a),2),round(mean(x$boot.results$se.b),2),round(mean(x$boot.results$se.c.prime),2))
    if (x$mediation.test$a.sig==1) plabs[1] <- paste(plabs[1],"*",sep=" ")
    if (x$mediation.test$b.sig==1) plabs[2] <- paste(plabs[2],"*",sep=" ")
    if (x$mediation.test$c.sig==1) plabs[3] <- paste(plabs[3],"*",sep=" ")
    plabs <- paste(plabs,"\n(",plabsse,")",sep="")
  }
  text(arrpos[1,1] + .05, arrpos[1,2] + .15, plabs[2])
  text(arrpos[2,1] -.15, arrpos[2,2] - .1, plabs[1])
  text(arrpos[3,1] -.12, arrpos[3,2] -.05, plabs[3])
}

## now need an xtable function and this puppy is ready for others
