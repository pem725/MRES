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

