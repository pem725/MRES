\name{zfa}
\alias{zfa}
\title{Compute Unit-Weighted Factor Scores}
\description{Standardizes items and computes a mean score (single factor) of available
  items.  Also provides information concerning the coherency (item total
  correlation) of the items to the computed factor score.  Useful in
  data reduction procedures where differential weighting as used in
  \code{\link{factanal}} via regression scores might not be desirable.
}
\usage{
zfa(x,use="complete.obs")
\method{summary}{zfa}(object)
\method{residual}{zfa}(object)
\method{plot}{zfa}(object,...)
}
\arguments{
  \item{x}{Matrix or data.frame of variables or items to be reduced.}
  \item{use}{Method of handling missing observations for the item-total correlations.}
}
\value{
  \code{zfa} returns objects of class \code{"zfa"} including the
  components
  \item{zdat}{Standardized scores of the original data.}
  \item{scores}{Unit-weighted scores computed by taking the mean of the
    available standardized values.}
  \item{itcor}{Item-total correlations between the unit-weighted factor
    scores and the original items/variables.}
  \item{resid}{Residual scores for each of the standardized
    items/variables.}
}
\details{Items or variables are standardized (see \code{\link{scale}})
  with a mean of zero (0) and a standard deviation of 1 by variable.
  The standardized values are then aggregated into a single score via
  \code{\link{rowMeans}} to form a composite score that is
  unit-weighted.  The function produces an item-total correlation matrix
  whereby each item/variable is correlated with the single,
  unit-weighted factor score.  Furthermore, residuals are retained to
  allow analysts to assess whether the unit-weighted factor score
  sufficiently accounts for invidual item/variable variance.

  Summary returns the item-total correlation for each item/variable in
  the analysis in a column vector.  Item-total correlations exclude the
  specific item so the values are not inflated by including the target
  item.  Some refer to this method as "corrected" item-total
  correlations.  The plot method returns a bivariate plot for each
  standardized item/variable in the analysis on the x-axis and the
  unit-weighted factor scores on the y-axis.  Each variable is colored
  according to the primary color palette and may be changed according to
  the \code{\link{col}} specification in \code{\link{par}}.

  The procedure of unit-weighting item totals in psychometric scaling
  has been around for a long time but recently has come back as a
  preferred method for instrument scoring with missing data.  According
  to several monte carlo simulations, the available item method - the
  very procedure used in this function - tends to provide robust
  estimates of construct scores.
}
\references{
  Grice, J.W. and Harris, R.J. (1998).  A comparison of regression and
  loading weights for the computation of factor scores.
  \emph{Multivariate Behavioral Research}, \bold{33}(2), 221-247.

 Morris, J.D. (1979).  A comparison of regression prediction accuracy on
 several types of factor scores.  \emph{American Educational Research
   Journal}, \bold{16}(1), 17-24.

 Rosenberg, N., Kaplan, H., and Skordahl, D.M. (1961). Validation of the
 Army Fixed-Wing Aptitude Batter against success in Army flight
 training.  \emph{U.S. Army TAG Research & Development Command, Human
 Factors Research Branch Technical Research Note}, \bold{No. 112}, 19.
}
\seealso{\code{\link{factanal}}
}
\examples{
data(fat,package="UsingR")
fat.zfa <- zfa(fat[,10:19])
summary(fat.zfa)
plot(fat.zfa)
}
\author{Patrick E. McKnight <pmcknigh@gmu.edu>}
\keyword{factor}
\keyword{scores}
\keyword{scale}
