
##########################################################################
# This section contains two different functions, one to plot all ICC's
# for a test and one to plot individual item ICC with observed values
#----------------------------------------------------------------
# FUNCTION 1:  paicc or "Plot all item characteristic curves"
#----------------------------------------------------------------
#
# Syntax:
#
# data = an R dataframe
# nitmes = the number of items to read
# parms = the number of parameters to read
# legend = do you want a legend on the plot - very messy with a large number of items


paicc <- function(data=NULL, nitems=NULL, parms=NULL, legend=F){
	par(bg="grey")
	theta <- seq(-4, 4, by=.001)
	lcolors <- seq(20,657, by=5)

# Rasch Model Plots
        
if(parms==1){
  
	for(i in 1:nitems){
          	prob <- (1/(1+exp(-1.7*(theta - (data$V2[i])))))
              
		plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1),
		col=lcolors[i], main="ICC's for all items")
           
		par(new=T)
	}
        
        if(legend==T){
          legend(-4,1,legend=c(data$V15[1:nitems]),lty=c(rep(1,nitems)),col=c(lcolors[1:nitems]))
        }
}
        

# Birnbaum Model Plots
        
if(parms==2){
	for(i in 1:nitems){
		prob <- (1/(1+exp(-(data$disc[i])*(theta - (data$diff[i])))))
                plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1),
		col=lcolors[i], main="ICC's for all items")
                
		par(new=T)
	}
        if(legend==T){
          legend(-4,1,legend=c(data$inum[1:nitems]),lty=c(rep(1,nitems)),col=c(lcolors[1:nitems]))
        }
}

# 3PL Birnbaum Model Plots
        
if(parms==3){
	for(i in 1:nitems){
		prob <- (1/(1+exp(-(data$disc[i])*(theta - (data$diff[i]))))) + guess[i]
                plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1),
		col=lcolors[i], main="ICC's for all items")
                par(new=T)
	}
        if(legend==T){
          legend(-4,1,legend=c(data$inum[1:nitems]),lty=c(rep(1,nitems)),col=c(lcolors[1:nitems]))
        }
}
par(new=F, bg="white")
}


#--------------------------------------------------------------------------------
# FUNCTION 2: piicc or "Plot Individual Item Characteristic Curves"
#--------------------------------------------------------------------------------


# Syntax:
# idat = R object (dataframe) containing item level data
# pdat = file containing the raw responses - doesn't make sense if it is not binary right now
# edat = file containing the person (ability) parameters
# item = the item number you wish to plot
# idlength = the offset from the first column to the first item - often used for the person ID

piicc  <- function(idat,pdat,edat,item,parms,offset=0){
	theta <- seq(-4, 4, by=.001)
	selout <- offset + item - 1
        if(parms==1){
          prob <- (1/(1+exp(-1.7*(theta - (idat$diff[item])))))
        }
        if(parms==2){
          prob <- (1/(1+exp(-(idat$disc[item])*(theta - (idat$diff[item])))))
        }
        if(parms==3){
          prob <- (1/(1+exp(-(idat$disc[item])*(theta - (idat$diff[item]))))) + guess[item]
        }
	par(xaxp=c(-4,4,20),yaxp=c(0,1,10),xaxs="i",yaxs="i")
	plot(theta, prob, type="l", xlim=c(-4,4), ylim=c(0,1), main="ICC and observed data for item", col="blue" )
	par(new=T)
	responses <- read.fwf(pdat,width=c(selout,1))
	edat <- scan(edat, list(ptheta=0, junk=""))
	freqdat <- data.frame(resp=responses$V2,ptheta=signif(edat$ptheta, digits=1))
	gres <- aggregate(freqdat, list(agtheta=freqdat$ptheta), FUN=mean)
	gresn <- aggregate(freqdat, list(agtheta=freqdat$ptheta), FUN=sum)
	agtheta <- as.vector(gres$agtheta)
	symbols(agtheta, gres$resp, circles=gresn$resp, xlab="", ylab="", xlim=c(-4,4), ylim=c(0,1), col = "red")
	par(new=F)
} 




############################
##  plottwo function      ##
##  added 10/17/05        ##
##  description in header ##
############################

plottwo <- function(pre,post,skip=3,title,items=1,ident=1,scores=0,xlab=NULL,ylab=NULL){
if (scores==0){
	if (items==1){
		pre <- read.table(pre, header=T, skip=skip)
		post <- read.table(post, header=T, skip=skip)
		pre <- subset(pre,pre$ST==1)
		post <- subset(post,post$ST==1)
		plot(pre$MEASURE,post$MEASURE,xlim=c(-5,5),ylim=c(-5,5),xlab=xlab,ylab=ylab,main=title)
		par(new=T)
		x1 <- seq(-5,5)
		y1 <- seq(-5,5)
		plot(x1,y1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l")
		xse <- mean(pre$ERROR)
		yse <- mean(post$ERROR)
		x2plus <- x1 + 3*xse
		x2minus <- x1 - 3*xse
		y2plus <- y1 + 3*yse
		y2minus <- y1 - 3*yse
		par(new=T)
		plot(x2minus,y2plus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=T)
		plot(x2plus,y2minus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=F)
		if (ident==1){
			identify(pre$MEASURE,post$MEASURE,labels=pre$NAME)
		}
	}
	if (items==0){
		pre <- read.table(pre, header=T, skip=skip)
		post <- read.table(post, header=T, skip=skip)
		pre <- subset(pre,pre$ST==1)
		post <- subset(post,post$ST==1)
		comb <- merge(pre,post,by="NAME")
		comb <- subset(comb)
		plot(comb$MEASURE.x,comb$MEASURE.y,xlim=c(-5,5),ylim=c(-5,5),xlab=xlab,ylab=ylab,main=title)
		par(new=T)
		x1 <- seq(-5,5)
		y1 <- seq(-5,5)
		plot(x1,y1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l")
		xse <- mean(pre$ERROR)
		yse <- mean(post$ERROR)
		x2plus <- x1 + 3*xse
		x2minus <- x1 - 3*xse
		y2plus <- y1 + 3*yse
		y2minus <- y1 - 3*yse
		par(new=T)
		plot(x2minus,y2plus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=T)
		plot(x2plus,y2minus,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l", lty="dashed")
		par(new=F)
	}	
}
if (scores==1){
		pre <- read.table(pre, header=T, skip=skip)
		post <- read.table(post, header=T, skip=skip)
		pre <- subset(pre,pre$ST==1)
		post <- subset(post,post$ST==1)
		comb <- merge(pre,post,by="NAME")
		comb <- subset(comb)
		plot(comb$SCORE.x,comb$SCORE.y,xlim=c(0,40),ylim=c(0,40),xlab=xlab,ylab=ylab,main=title)
		par(new=T)
		x1 <- seq(0,40)
		y1 <- seq(0,40)
		plot(x1,y1,xlab="",ylab="",xlim=c(0,40),ylim=c(0,40),type="l")
		xse <- sd(comb$SCORE.x)*sqrt(1-.8)
		yse <- sd(comb$SCORE.y)*sqrt(1-.8)
		x2plus <- x1 + 2*xse
		x2minus <- x1 - 2*xse
		y2plus <- y1 + 2*yse
		y2minus <- y1 - 2*yse
		par(new=T)
		plot(x2minus,y2plus,xlab="",ylab="",xlim=c(0,40),ylim=c(0,40),type="l", lty="dashed")
		par(new=T)
		plot(x2plus,y2minus,xlab="",ylab="",xlim=c(0,40),ylim=c(0,40),type="l", lty="dashed")
		par(new=F)
	}
}

## pre post plot by two different versions - plots the item numbers - useful for BRS plots

plotitwobytwo <- function(x1,x2,y1,y2,skip=3,title,xlab=NULL,ylab=NULL){
  library(sfsmisc)

  # Read in data from Big/Winsteps - note, only if files
  pre1 <- read.table(x1, header=T, skip=skip)
  post1 <- read.table(x2, header=T, skip=skip)
  pre2 <- read.table(y1, header=T, skip=skip)
  post2 <- read.table(y2, header=T, skip=skip)
  pre1 <- subset(pre1,pre1$ST==1)
  post1 <- subset(post1,post1$ST==1)
  pre2 <- subset(pre2,pre2$ST==1)
  post2 <- subset(post2,post2$ST==1)
  
  # plot 1
  n.plot(pre1$MEASURE,post1$MEASURE,nam=pre1$NAME,col="red",xlim=c(-5,5),ylim=c(-5,5),xlab=xlab,ylab=ylab,main=title)
#  abline(0,1) # add identity line
  par(new=T) # prepare for additional graph

  # plot 2
  n.plot(pre2$MEASURE,post2$MEASURE,nam=pre2$NAME,col="blue",xlim=c(-5,5),ylim=c(-5,5),xlab="",ylab="",main="")
  par(new=T)

  x1 <- seq(-5,5)
  y1 <- seq(-5,5)

  # old code for identity line - see x1 and y1 below
  plot(x1,y1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),type="l")
  par(new=T)

  # plot 1 standard error bars

  xse1 <- mean(pre1$ERROR)
  yse1 <- mean(post1$ERROR)
  x2plus1 <- x1 + 3*xse1
  x2minus1 <- x1 - 3*xse1
  y2plus1 <- y1 + 3*yse1
  y2minus1 <- y1 - 3*yse1
  par(new=T)
  plot(x2minus1,y2plus1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="red",type="l", lty="dashed")
  par(new=T)
  plot(x2plus1,y2minus1,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="red",type="l", lty="dashed")
  par(new=T)

  # plot 2 standard error bars
  xse2 <- mean(pre2$ERROR)
  yse2 <- mean(post2$ERROR)
  x2plus2 <- x1 + 3*xse2
  x2minus2 <- x1 - 3*xse2
  y2plus2 <- y1 + 3*yse2
  y2minus2 <- y1 - 3*yse2
  par(new=T)
  plot(x2minus2,y2plus2,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="blue",type="l", lty="dashed")
  par(new=T)
  plot(x2plus2,y2minus2,xlab="",ylab="",xlim=c(-5,5),ylim=c(-5,5),col="blue",type="l", lty="dashed")
  par(new=F)
  
}



## END IPU


