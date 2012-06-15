## This is the beginning code to test Bob Holt's agreement indices
## Latest Revision Date:  11/04/2011

## generate example data
dat <- data.frame(id=1:10,rater1=c(0,1,2,3,4,5,6,7,8,9),rater2=c(1,2,3,4,5,6,7,8,9,9),rater3=c(9,8,7,6,5,4,3,2,1,0),rater4=c(3,2,5,6,9,8,1,2,4,0),group=gl(2,5))

dat.l <- reshape(IRRdatEx,varying=2:5,timevar="rater",v.names="rating",sep="",direction="long")


omegasq <- function(x,rater=NULL,data=NULL){
  ## sensitivity formula for ratings
  ## data must be in long format
  ##   where x is a formula as rating~group
  ##         rater is the rater variable
  ##         group is the grouping variable to compare raters across groups
  ##         data is the data.frame that contains the above information

  ## ANOVA by rater to determine sensitivity by group assignment
  raterVar <- match(rater,names(data))
  raters <- unique(data[,raterVar])
  out <- as.data.frame(matrix(NA,length(raters),2))
  for (i in 1:length(raters)){
    a <- summary(aov(x,data=subset(data,data[,raterVar]==raters[i])))
    out[i,] <- c(raters[i],a[[1]][[3]][1] / sum(a[[1]][[3]]))
  }
  names(out) <- c("rater","OmegaSqrd")
  return(out)
}

## run code above with this line:
omegasq(rating~group,rater="rater",data=dat.l)


sysdiff <- function(x,obj,rater,data){
  # systematic (holt's t) differences formula for ratings
  # consists of a t-test between the rater and group averages
  # where...
  #  x is the rating
  #  obj is the variable that contains the object being rated
  #  rater is the variable that contains the rater variable
  #  data is the data.frame in long format
  x <- data[,match(x,names(data))]
  obj <- data[,match(obj,names(data))]
  rater <- data[,match(rater,names(data))]


  out <- data.frame(rater=unique(rater),t.Test=NA)
  
  grpmeans <- aggregate(x,list(obj),mean)
  

  
  num <- mean(x$rater - x$group)
  den <- var(x$rater - x$group)
  t.out <- num/den
}


rwg <- function(x,scale=NULL,group=NULL,data=NULL){
  ## agreement formula for ratings
  ## rwg (holt's agreement) = 1 - (var(ratings)/E(uniform ratings))
  ## where x is a vector of ratings
  ##       scale is a vector that includes all possible ratings
  ##       group is an optional string denoting a conditioning variable
  ##       data is the data frame
  rwgfcn <- function(x,scale=scale){
    return(1-(var(x)/var(rep(scale,1000))))
  }
  x <- data[,match(x,names(data))]
  if(is.null(group)){
    rwg <- rwgfcn(x,scale)
  } else {
    grplvls <- unique(data[,match(group,names(data))])
    rwg <- data.frame(group=grplvls,rwg=NA)
    for (i in 1:length(grplvls)){
      tmpdat <- x[data[,match(group,names(data))]==grplvls[i]]
      rwg[i,2] <- rwgfcn(tmpdat,scale=scale)
    }
    names(rwg) <- c(group,"rwg")
  }
  return(rwg)
}


rwg("rating",0:9,data=dat.l)
rwg("rating",0:9,"group",data=dat.l)


consist <- function(x){
  # consistency formula for ratings
  
}

congruency <- function(x){
  # congruency formula for ratings
}
