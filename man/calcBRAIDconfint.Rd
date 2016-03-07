\name{calcBRAIDconfint}
\alias{calcBRAIDconfint}
\title{
Calculate General BRAID Confidence Intervals
}
\description{
Produces confidence intervals for any general calculable property of a BRAID surface, based on a
given BRAID surface fit (with bootstrapped coefficients).
}
\usage{
calcBRAIDconfint(bfit, parfunc, civals = NULL)
}
\arguments{
  \item{bfit}{ A BRAID surface fit of class \code{braidrm}, with bootstrapped coefficients }
  \item{parfunc}{ A function to calculate a property of a response surface.  \code{parfunc} must take a ten-element BRAID 
parameter vector as input, and return a single numeric value or a numeric vector of a fixed length as output.}
  \item{civals}{ Values specifying the lower and upper bound of the confidence interval(s) to be calculated. If \code{NULL}
(the default), the values of the \code{ciLevs} property of \code{bfit} will be used. NOTE: using a wider percentile
here than was used to estimate the bootstrapped coefficients can lead to an unstable confidence interval.  }
}
\value{
If \code{parfunc} produces a single value as output, a vector of three values: the second value is the estimate of the best
BRAID surface fit, while the first and third values represent the lower and upper bounds of the calculated confidence interval.
If \code{parfunc} produces a vector as output, an array of width three: the second column is the estimate of the best BRAID
surface fit, while the first and third columns represent the lower and upper bounds of the calculated confidence intervals.
}
\author{
Nathaniel R. Twarog
}
\examples{
data(es8olatmz)
\dontshow{
es8olatmz <- es8olatmz[(es8olatmz$plate==8&grepl("[AB]",es8olatmz$well))|es8olatmz$plate==9|
							(es8olatmz$plate==27&grepl("(01|02)",es8olatmz$well)),]
}
brd <- braidrm(act~conc1+conc2,es8olatmz,fixed=c(NA,NA,NA,NA,1,NA,NA,-4,-4,-4))
potci <- calcBRAIDconfint(brd,function(parv) invertBRAIDrsm(-1,DA=10^-7,parv=parv))
\dontrun{
pconc <- exp(seq(log(10^-6),log(10^-3),length=100))
effci <- calcBRAIDconfint(brd,function(parv) evalBRAIDrsm(10^-7,pconc,parv))
}
}

\seealso{\code{\link{braidrm}},\code{\link{getBRAIDbootstrap}}}
\keyword{ regression }