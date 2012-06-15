JND <- function(dv,iv,lambda=0){
  UseMethod("JND")
}
  

JND.default <- function(dv,iv,lambda=0){
  glm.out <- glm(dv~iv,family=binomial)
  gamma <- predict(glm.out,data.frame(iv=0),type="response",se.fit=T)
  out <- list(gout=glm.out,gamma=gamma,dv=dv,iv=iv)
  class(out) <- "JND"
  invisible(out)
  return(out)
}

summary.JND <- function(x){
  gout <- x$gout
  gamma <- x$gamma
  print(summary(gout),digits=2)
  cat("-------------------------\n\n")
  cat("gamma = ",round(gamma[[1]][[1]],2),"\n")
  cat("gamma(SE) = ",round(gamma[[2]][[1]],2),"\n")
}


plot.JND <- function(x,col.cur="red",lwd.cur=4,annotate=T,atype="e",xlab="Difference",ylab="Probability",legend=T){
  gamma <- x$gamma
  x.new <- seq(0,1.5*max(x$iv), len=1000)
  y.new <- predict(x$gout,data.frame(iv=x.new),type="response",se.fit=T)
  dat <- data.frame(x.new,yhat=y.new[[1]],se=y.new[[2]])
  pse <- dat[dat$yhat>.49999999,][1,1]
  jnd <- dat[dat$yhat>.7499999,][1,1]
  pse.se <- dat[dat$yhat > (dat[dat$yhat>.4999999,][1,3] + dat[dat$yhat>.49999999,][1,2]),1] - pse
  jnd.se <- dat[dat$yhat > (dat[dat$yhat>.7499999,][1,3] + dat[dat$yhat>.74999999,][1,2]),1] - jnd


  
  plot(jitter(x$dv,amount=0)~x$iv,xlab=xlab,ylab=ylab,xlim=c(0,max(x.new)))
  xs <- c(x.new,sort(x.new,decreasing=T))
  ys <- c(y.new[[1]] + 1.96*y.new[[2]],rev(y.new[[1]] - 1.96*y.new[[2]]))
  polygon(xs,ys,col="orange",density=20,angle=100,border=NA)
  
  lines(x.new, y.new[[1]], lwd = lwd.cur, col = col.cur)
  lines(x.new, y.new[[1]] + 1.96*y.new[[2]], lwd = 2, lty=2, col=col.cur)
  lines(x.new, y.new[[1]] - 1.96*y.new[[2]], lwd = 2, lty=2, col=col.cur)
 

  if(annotate==T){
    arrows(0,x$gamma[[1]],0,.5,length=.1,angle=90,code=3)
    text(.05*max(x.new),median(c(x$gamma[[1]],.5)),expression(.5 - gamma))
    segments(0,.5,pse,.5,col="blue",lwd=2)
    segments(pse,0,pse,.5,col="blue",lwd=2)
    segments(0,.75,jnd,.75,col="blue",lwd=2)
    segments(jnd,0,jnd,.75,col="blue",lwd=2)
    text(pse+.05*max(x.new),.25,"PSE")
    text(jnd+.05*max(x.new),.75/2,"JND")
  }

  if(legend==T){
    L <- list("PSE =",
              "JND =",
              bquote(paste(gamma,"  =  ")),
              bquote(.(round(pse,2)) %+-% .(round(pse.se,2))),
              bquote(.(round(jnd,2)) %+-% .(round(jnd.se,2))),
              bquote(.(round(gamma[[1]][[1]],2)) %+-% .(round(gamma[[2]][[1]],2))))
    legend("bottomright",legend=do.call("expression",L),ncol=2,pch=c('','','',''),x.intersp=0.1,title="Critical Values")
  }
#  return(dat)  
}

