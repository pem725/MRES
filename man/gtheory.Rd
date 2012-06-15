\name{gtheory}
\alias{Gstudy}
\alias{Dstudy}
\title{Generalizability Theory Functions for G-studies and D-studies}
\description{Computes g-study variance components for fully crossed
  models and d-study coefficients including phi coefficients and
  generalizability coefficients for pre-specified and user-specified
  alternative designs.
}
\usage{gtheory(formula,data=NULL,gcrit=.8,phicrit=.8)
\method{summary}{gtheory}(object)
\method{print}{gtheory}(object)
\method{plot}{gtheory}(object,coef=c("g","phi"),legend=T,ablcol=1,xlab=NULL,ltitle=NULL,legposx="topleft",legposy=NULL)
\method{xtable}{gtheory}(object,gout=c("g","d"),caption=NULL,label=NULL,align=NULL,digits=NULL,display=NULL,...)
}
\arguments{
  \item{formula}{a formal two-sided expression similar to those used in
    aov and lm functions.  The left hand side is the response variable
    and the right hand side are the relevant facets expressed in a
    full-factorial model}
  \item{data}{a data frame object from a fully-crossed model whereby each subject has
    values for each of the facets of interest.  The data must be formatted in a long shape (see example below) where
    at least one facet is combined with persons.}
  \item{gcrit,phicrit}{critical values for the g-coefficient and phi-coefficient that influence
    the computation of the d-study and the plot output. Arbitrary
    values of 0.8 were set as defaults but users should select more
    appropriate values according to their research question and field expectations.}
  \item{object}{an object of class \code{"gtheory"}}
  \item{coef}{used in the plot method to differentiate between
    g-coefficients (relative decisions) by using \dQuote{g} (the default) and phi-coefficients (absolute
    decisions) by using \dQuote{phi}.}
  \item{legend}{a logical indicating whether an legend ought to be
    produced for the d-study plot.}
  \item{ablcol}{a numeric or character vector indicating a color for the
    critical value line in the d-study plot.  The default value is 1 for
    the color black.}
  \item{xlab}{an optional x-axis label to be specified as quoted text. If no text is
  provided, the default is to use paste the string \dQuote{Levels of} to
  to the name of the second predictor.}
  \item{ltitle}{an optional legend title to be specified as quoted
    text.  If no legend title is provided, the default title used is the
    name of the third predictor.  This feature is relevant only for
    2-facet or three predictor g-theory models.}
  \item{legposx,legposy}{legend position in x and y coordinates,
    respectively.  The values should be within the range of plot values
    for the second predictor in a 2-facet (i.e., three predictor)
    g-theory model.}
}
\value{
  \code{gtheory} returns an list object of class \code{"gtheory"} including
  the components
  \item{aout}{an aov object used to compute the generalizability model}
  \item{gs}{a data.frame consisting of the results from the g-study}
  \item{ds}{a data.frame consisting of the results from the d-study}
  \item{gcrit,phicrit}{the specified critical values for g and phi}
  \item{Nparms}{the number of parameters in the formula (i.e.,
    predictors)}
}

\details{Generalizability theory is a statistical method originally
  developed for psychometricians to estimate the influence of various
  facets to the measurement of abilities in people.  The theory expanded
  to a broader context where researchers were able to estimate the
  contribution of predictors in standard linear models to not just
  persons but to groups.  Thus, generalizability is general tool for
  estimating variance components (g-studies) and forecasting decisions
  (d-studies) for future research.

  The purpose of this function and its associated methods is to allow
  researchers to use generalizability theory without having to conduct
  hand calculations for any part.  Previously, data analysts needed to
  either use specialized software to compute expected variance
  components and then compute generalizability study and decision study
  coefficients.  These steps were tedious and often precluded most people
  from using generalizability theory - even where it may be appropriate.
  This gstudy function does all the computations necessary for one-facet
  and two-facet, fully-crossed designs (later versions will include more
  complicated models and designs).  Furthermore,

  All calculations follow the standard EMS method first introduced by
  Cronbach et. al. (1963) and later popularized by Shavelson and Webb
  (1991).  Specifically, the variance components are estimated using a
  fully-saturated ANOVA model (via \link{\code{aov}}) and then using
  Kirk's (1982) reverse method of computing expected mean squares and
  then variance components.  Negative variance components are set to
  zero to eliminate any problems in estimation.

  The first facet specified in each model is
  the facet of generalization.  In many cases, researchers will specify
  persons but other options may be raters, items, occasions, etc.
  
}

\references{
  Brennan, R.L. (2001).  \emph{Generalizability theory}.  New York; Springer-Verlag.
  
  Cronbach, L.J., Gleser, G.C., Nanda, H., and Rajaratnam,
  N. (1972). \emph{The dependability of behavioral measurements: Theory
  of generalizability for scores and profiles}. New York: John Wiley.

  Cronbach, L.J., Nageswari, R., and Gleser, G.C. (1963). Theory of
  generalizability: A liberation of reliability theory. \emph{The British
  Journal of Statistical Psychology}, \bold{16}, 137-163.

  Shavelson, R.J. and Webb, N.M. (1991). \emph{Generalizability theory: A
  primer}. Newbury Park, CA: Sage.
}
\seealso{
}
\examples{
data(CTBS)
CTBS.long <- reshape(CTBS,idvar="Person",varying=list(names(CTBS[-1])),v.names="response",times=1:8,timevar="item",direction="long")
CTBS.gt <- gtheory(response~as.factor(Person)*as.factor(item),data=CTBS.long)
summary(CTBS.gt)
}
\author{Patrick E. McKnight <pmcknigh@gmu.edu>}
\keyword{generalizability theory}
\keyword{variance components}
\keyword{intraclass correlation}
\keyword{aov}