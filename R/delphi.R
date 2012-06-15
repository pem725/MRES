## usage delphi(proj,iter)
## example:  delphi("mTBI",1)
## where x is the project name
## iter is the iteration number

## delphi <- function(proj,iter){
##     UseMethod("delphi")
## }

# the following function retrieves the data from the local host.  I
# need to make sure it agrees with the data structure from the
# database.

## Get.Data <- function(proj,iter){
##   library(RMySQL)
##   drv <- dbDriver("MySQL")
##   con <- dbConnect(drv, user="delphi",password="delphi-user",host="localhost",dbname="delphi")
##   proj.table <- paste(x,"responses",sep=".")
##   # need to include question number, question type, agreement tolerance in the sql statement below
##   sql <- paste("SELECT * FROM project_items_by_user WHERE iter =",iter,"AND project_id =",proj)
##   res <- dbSendQuery(con, statement=sql)
##   dat <- fetch(res)
##   return(dat)
## }

# the following function uses the data from the database to compute
# agreement.

delphi <- function(x,qtype){
  #library(irr)
  #dat <- Get.Data(x,iter)
  
  # dat should include question number, question type, agreement
  # tolerance, rater, and rater responses
  
  # quests <- unique(x$question)
  
  ## out <- rep(NA,length(quests))
  ## for (i in 1:length(quests)){
  ##   tmp <- dat[dat$question == quests[i],]
  ##   # process the agreement method based upon the question types

    if (qtype == 1){
    # out <- sd(tmp$rating)/sqrt(nrow(tmp)) # continuous data qtype=1
      Range <- max(x) - min(x)
      out <- (sd(x)/length(x))/Range # continuous data qtype=1
    }
    if (qtype == 2){
      out <- 1 - (mean(tmp$rating)*(1-mean(tmp$rating))/.5) # binary data qtype=2
    }
    if (qtype == 3){
      out <- 1 - chisq.test(tmp$rating,simulate.p.value=T)$p.value # categorical values qtype=3
    }
    # create graph here - now created using the google API
  #} # end of for loop
  out <- round(out,2)
  return(out)
}


# upload data back to the database NOW HANDLED BY PHP

## Post.Data <- function(x){
##   #TBD
## }


# create some bogus data for 3 questions and 5 raters

dat <- data.frame(proj=rep("test",15),iter=rep(1,15),question=c(rep(1,5),rep(2,5),rep(3,5)),qtype=c(rep(1,5),rep(2,5),rep(3,5)),rater=rep(1:5,3),rating=c(rnorm(5),rbinom(5,2,.2),1:5))

res <- data.frame(project=unique(dat$proj),question=unique(dat$question),iter=unique(dat$iter),agree=rep(NA,length(unique(dat$question))))

## for (i in 1:nrow(res)){
##      res$agree[i]
##    }

## check how the delphi agreement indices operate with bogus data

# continuous data

out <- c(NA,NA,NA)
for (i in 1:100){
  for (j in 1:20){
    k <- j/8
    dat <- rnorm(10,mean=5,sd=k)
    res <- delphi(dat,1)
    res <- c(i,j,res)
    out <- rbind(out,res)
  }
}

out.cont <- data.frame(out)
plot(out[,2],out[,3])
abline(lm(out[,3]~out[,2]),col="red",lwd=2)




## compare my agreement index with irr results

# continuous data

out <- c(NA,NA,NA)
for (i in 1:200){
  for (j in 1:5){
    dat <- rnorm(10,mean=5,sd=j)
    res1 <- delphi(dat,1)
    res2 <- 1 - finn(dat)$p.value
    res <- c(i*j,res1,res2)
    out <- rbind(out,res)
  }
}

out.cont <- data.frame(out)
plot(out[,2],out[,3])
abline(1,0)


# binary data

# categorical data



##
     
