
## Intraclass correlations - both fixed (icc.fix) and random (icc.ran)

sficc <- function(rating,judge,person,data,fixed=T){
  UseMethod("sficc")
}

sficc.default <- function(rating,judge,person,data,fixed=T){
  if(fixed==T){
    rating <- data[,rating]
    judge <- data[,judge]
    person <- data[,person]
    judge <- as.factor(judge)
    person <- as.factor(person)
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
    res <- data.frame(round(c(sf11,sf21,sf31,sf1k,sf2k,sf3k),2))
    names(res) <- "icc"
    row.names(res) <- c("ICC(1,1)","ICC(2,1)","ICC(3,1)","ICC(1,k)","ICC(2,k)","ICC(3,k)")
    out <- list(res)
    class(out) <- "sficc"
    return(out)
  }


    if(fixed==F){
#      icc.ran <- function(resp,fx,rx){
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
}

 
## END ICC


##      cat("ICC(1,1)\tICC(2,1)\tICC(3,1)\tICC(1,k)\tICC(2,k)\tICC(3,k)\n")
##      cat(sf11,"\t",sf21,"\t",sf31,"\t",sf1k,"\t",sf2k,"\t",sf3k,"\n")


summary.sficc <- function(x){
  ICC <- x[[1]]$icc
  DETAILS <- c("Each subject rated by multiple, same N raters who are assumed to be randomly assigned to subjects.","All subjects rated by same raters who are assumed to be a random subset of all possible raters.","All subjects rated by the same raters who are the entire population of raters.","same assumptions as ICC(1,1) but reliability for the mean of k ratings.","same assumptions as ICC(2,1) but reliability for the mean of k ratings.","same assumptions as ICC(3,1) but reliability for the mean of k ratings and no subject by judge interaction.")
  out <- data.frame(ICC,DETAILS)
  row.names(out) <- c("ICC(1,1)","ICC(2,1)","ICC(3,1)","ICC(1,k)","ICC(2,k)","ICC(3,k)")
  return(out)
}


print.sficc <- function(x){
  return(x[[1]])
}
