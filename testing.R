## The following code are tests and examples for the MRES package

CTBS <- read.csv("./data/CTBS.dat")
CTBS.long <- reshape(CTBS,idvar="Person",varying=list(names(CTBS[-1])),v.names="response",times=1:8,timevar="item",direction="long")
CTBS.long$person.f <- as.factor(CTBS.long$Person)
CTBS.long$item.f <- as.factor(CTBS.long$item)
CTBS.gt <- gtheory(response~person.f*item.f,data=CTBS.long)
summary(CTBS.gt)

a.test <- function(formula,data){
  aout <- aov(formula,data)
  return(summary(aout))
}

a.test.out <- a.test(response~person.f*item.f,data=CTBS.long)


aout <- aov(response~person.f*item.f,data=CTBS.long)

BehObsMeas <- data.frame(child=as.factor(rep(1:13,4)),occ=as.factor(c(rep(1,26),rep(2,26))),rater=as.factor(rep(c(rep(1,13),rep(2,13)),2)),rating=c(0,3,2,1,1,4,1,2,1,1,1,1,2,1,4,2,2,2,4,1,2,1,1,2,2,1,1,1,1,0,2,3,2,0,1,1,1,1,1,2,2,0,1,1,4,1,0,2,0,1,1,1))


