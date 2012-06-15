
ioc <- function(x,ID=NULL,first.item=NULL,nitems=NULL,responses=NULL,anchor=1){

# -----------------------------------------------------------------------------
#  IOC function to compute both item and person parameters from polytomous data
#  Authors:  Patrick E. McKnight, Ph.D.
#            Dept. of Psychology
#            University of Arizona
#
#            Robert Massof, Ph.D.
#            Lions Vision Center
#            Johns Hopkins University
# -----------------------------------------------------------------------------
#    x = data.frame to analyze

#    ID = subject identifier - very important for person parameters

#    first.item = what column contains the first item - assumes the
#    items are contiguously organized with no unnecessary variables
#    between the ones to be analyzed.

#    nitems = the number of items to analyze

#    responses = the possible responses to the items.  For example, a
#    five-response category Likert scale would use the following code:
#    responses=c(1,2,3,4,5).

#    anchor = the item you wish to anchor all item parameters.  The
#    default is item 1 but it may be specified to correspond to any item
#    in the dataset.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#  STEP ONE:  Compute Item Parameters
# -----------------------------------------------------------------------------
  
# Convert dataframe into matrix and get a list of names to use for the variables

  attach(x)
  dat <- as.matrix(x)
  item.dat <- dat[,first.item:(first.item + nitems - 1)]
  inames <- dimnames(item.dat)[[2]]

# Create Table of zscore transformed Category Endorsements to be Analyzed

  maxresp <- max(responses)
  dimresp <- dim(as.array(responses))
  zitems <- matrix(nrow=maxresp,ncol=nitems)
	for (i in 1:nitems){
          for (j in 1:dimresp){
              normdat <- qnorm(cumsum(table(item.dat[,i]))/(sum(table(item.dat[,i]))))                
              zitems[j,i] <- normdat[j]
            }
        }
#  Trim off the "Inf" from the matrix

  zfin <- zitems[1:dim(zitems)[1] - 1,]
  zfin[zfin=="Inf"] <- NA
  
# Construct two matrices to fill with slope and intercept values

  slopemat <- matrix(nrow=nitems,ncol=nitems)
  intmat <- matrix(nrow=nitems,ncol=nitems)
  
# Now populate the matrices with slope and intercept values from the bivariate regressions

  	for (i in 1:nitems){
  	  for (j in 1:nitems){
  	    iv <- zfin[,i]   
	    dv <- zfin[,j]
	      lm1 <- lm(dv~iv)
	      slopemat[j,i] <- as.numeric(coef(lm1)[2])
	      intmat[j,i] <- as.numeric(coef(lm1)[1])
	  }
	}

# Construct new matrix to hold mean slopes and populate matrix with the difference values

  xslope <- matrix(nrow=nitems,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:nitems){
      xslope[j,i] <- (log(slopemat[i,j])-log(slopemat[j,i]))/2
    }
  }

# Anchor the scalar

  lnsd <- exp(xslope[,anchor])
  deltax <- matrix(nrow=nitems,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:nitems){
      deltax[j,i] <- ((intmat[j,i])*(lnsd[j]))
    }
  }

# Means minus anchor matrix

  xminanc <- matrix(nrow=nitems,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:nitems){
      xminanc[j,i] <- deltax[anchor,i]-deltax[j,i]
    }
  }

# Item Parameter Results

  parms <- matrix(nrow=nitems,ncol=3)
  for (i in 1:nitems){
      parms[i,1] <- inames[i]
      parms[i,2] <- mean(xminanc[i,]) * -1
      parms[i,3] <- sd(xminanc[i,])
  }
  
item.parms <<- data.frame(item=parms[,1],difficulty=as.numeric(parms[,2]),sd=as.numeric(parms[,3]))


# -----------------------------------------------------------------------------
# STEP TWO:  Compute Person Parameters
# -----------------------------------------------------------------------------

# Convert dataframe into matrix and get a list of names to use for the variables

  person.dat <<- t(dat[,first.item:(first.item + nitems - 1)])
  pnames <<- ID

# Create Person Table of zscore transformed Category Endorsements to be Analyzed

   maxresp <<- max(responses)
   dimresp <<- dim(as.array(responses))
   zitems <<- matrix(nrow=maxresp,ncol=nrow(item.dat))
 	for (i in 1:nrow(item.dat)){
           for (j in 1:dimresp){
               normdat <<- qnorm(cumsum(table(person.dat[,i]))/(sum(table(person.dat[,i]))))                
               zitems[j,i] <<- normdat[j]
             }
         }

# #  Trim off the "Inf" from the matrix

#   zfin <<- zitems[1:dim(zitems)[1] - 1,]
#   zfin[zfin=="Inf"] <<- NA
  
# # Construct two matrices to fill with slope and intercept values

#   slopemat <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   intmat <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
  
# # Now populate the matrices with slope and intercept values from the bivariate regressions

#   	for (i in 1:nrow(item.dat)){
#   	  for (j in 1:nrow(item.dat)){
#   	    iv <<- zfin[,i]   
# 	    dv <<- zfin[,j]
# 	      lm1 <<- lm(dv~iv)
# 	      slopemat[j,i] <<- as.numeric(coef(lm1)[2])
# 	      intmat[j,i] <<- as.numeric(coef(lm1)[1])
# 	  }
# 	}

# # Construct new matrix to hold mean slopes and populate matrix with the difference values

#   xslope <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   for (i in 1:nrow(item.dat)){
#     for (j in 1:nrow(item.dat)){
#       xslope[j,i] <<- (log(slopemat[i,j])-log(slopemat[j,i]))/2
#     }
#   }

# # Anchor the scalar

#   lnsd <<- exp(xslope[,anchor])
#   deltax <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   for (i in 1:nrow(item.dat)){
#     for (j in 1:nrow(item.dat)){
#       deltax[j,i] <<- ((intmat[j,i])*(lnsd[j]))
#     }
#   }

# # Means minus anchor matrix

#   xminanc <<- matrix(nrow=nrow(item.dat),ncol=nrow(item.dat))
#   for (i in 1:nrow(item.dat)){
#     for (j in 1:nrow(item.dat)){
#       xminanc[j,i] <<- deltax[anchor,i]-deltax[j,i]
#     }
#   }

# # Item Parameter Results

#   parms <<- matrix(nrow=nrow(item.dat),ncol=3)
#   for (i in 1:nrow(item.dat)){
#       parms[i,1] <<- pnames[i]
#       parms[i,2] <<- mean(xminanc[i,]) * -1
#       parms[i,3] <<- sd(xminanc[i,])
#   }
  
# person.parms <<- data.frame(item=parms[,1],difficulty=as.numeric(parms[,2]),sd=as.numeric(parms[,3]))

  
}

# Create test data

testme <- data.frame(i1=c(1,2,3,4,3,2,1,4,5,4,2,4,5,3,5),i2=c(2,3,5,1,3,5,2,2,5,2,4,5,2,4,2))

### END IOC
