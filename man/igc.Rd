\name{igc}
\alias{igc}
\title{Individual Growth Curve Calculation}
\description{Computes linear growth curves paramters for longitudinal
  data using either ordinary least squares or maximum likelihood.
}
\usage{
igc(x,idvar="",ivar="",dvar="",cvar=NULL,parms=2,method="OLS")
\method{summary}{igc}(object)
\method{coef}{igc}(object,prefix=dvar)
\method{plot}{igc}(object,cvar=NULL)
}
\arguments{
  \item{x}{Data frame or matrix that provides a long format for repeated
    measures.}
  \item{idvar}{The variable name - specified in quotes - that identifies
    the unit of analysis (usually subjects) measured over time.}
  \item{ivar}{The variable name that identifies each of the repeated
    measures by idvar.}
  \item{dvar}{The variable name that identifies the repeated measure.}
  \item{cvar}{An optional variable name for a conditioning variable to
    be used in the plot statement.}
  \item{parms}{The number of growth curve parameters to compute.  Parms
    defaults to 2 indicating that only linear slopes and intercepts are
    computed.  Alternatively, three parameters (3) may be specified to
    compute linear, quadratic, and intercept parameters by idvar.}
  \item{method}{The estimation method for computing individual growth
    curves.  igc defaults to "OLS" but users may specify "ML" where
    growth curve parameters are constrained (see details) in lmer.}
}
\value{
  \code{igc} returns a list of objects of class \code{"igc"} including the
  components
  \item{params}{The individual growth curve parameters returned as a
    data frame with idvar as the first column, intercept as the second,
    linear slope as the third, and the optional quadratic term as the
    fourth.  If parameters were computed with OLS then the data frame
    includes standard errors for each parameter and the adjusted $R^2$  }
  \item{parms}{The number of parameters computed - a value passed to
    other methods.}
  \item{method}{The specific method used to compute the growth curve parameters.}
  \item{xlim}{A vector for the x-axis used for the plot method.}
  \item{ylim}{A vector for the y-axis used for the plot method.}
  \item{fixed.parms}{The fixed effects parameters for the entire sample.}
  \item{cvar}{The name of the conditioning variable used later for plot
    and summary methods.}
  \item{dvar}{The name of the dependent variable used to model }
}
\details{Individual growth curve modeling is a flexible method for
  reducing repeated measures data.  Typically, data analysts use a
  single-step procedure such as mixed-effects models to estimate both
  within-subjects and between-subjects parameters.  The advantage of
  individual growth curves over mixed-effects models lies in the
  decoupling of the steps so that analysts can diagnose each step and
  formulate alternative hypotheses from additional analyses.
  Furthermore, growth curve parameter computation outside lmer (or its
  equivalent) often allows for larger sample sizes and greater
  flexibility for missing data handling procedures.

  The igc function uses the following generic model to estimate
  individual growth curves:

  \deqn{y = b_2 * x^2 + b_1 * x + b_0}

  where y is the specified dvar, x is the specified ivar, \eqn{b_0} is
  the estimated intercept, \eqn{b_1} is the linear parameter, and
  \eqn{b_2} is the optional quadratic parmaeter.  
  
}
\references{
}
\seealso{\code{\link{lmer}}
}
\examples{
}
\author{Patrick E. McKnight <pmcknigh@gmu.edu>}
\keyword{lmer}
\keyword{longitudinal}
\keyword{data reduction}
