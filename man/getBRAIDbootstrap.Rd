\name{getBRAIDbootstrap}
\alias{getBRAIDbootstrap}
\title{
Bootstrap BRAID Parameter Confidence Intervals
}
\description{
Estimate bootstrapped confidence intervals of a BRAID surface fit's parameters, by sampling residuals
and refitting.
}
\usage{
getBRAIDbootstrap(bfit, ciLevs = c(0.025, 0.975),numBoot=NULL,parallel=TRUE,cl=NULL)
}
\arguments{
  \item{bfit}{ a BRAID surface fit of class 'braidrm' }
  \item{ciLevs}{ a two-element vector specifying the lower and upper percentiles of the desired confidence interval. The default
  value of \code{c(0.025, 0.975)} results in a 95\% confidence interval.}
  \item{numBoot}{ the number of bootstrapped coefficients to be calculated.  If \code{NULL} (the default), number of
  bootstrapped coefficients is determined from the width of the interval.}
  \item{parallel}{ If \code{TRUE}, the function will attempt to calculate the bootstrapped parameter fits using multiple
parallel workers, implemented by the \code{boot} package.  If this package is not installed, the parameter will be
ignored. }
  \item{cl}{ An optional reference to an existing \code{parallel} or \code{snow} cluster, to be passed to \code{\link{boot}} }
}
\details{
This function constructs bootstrapped confidence intervals for objects of class \code{braidrm}.  If the input already has
bootstrapped confidence intervals (or has tried and failed to generate them), this function will throw a warning and
return the input model.

Parallelized bootstrapping is implemented by the \code{\link{boot}} function from the \code{boot} package, a
recommended package but not part of base R.  By default, the \code{parallel} parameter of \code{link{boot}}
is set to "snow", but can be altered by changing the option \code{"boot.parallel"} to one the following: \code{"no"},
signifying no parallel processing; \code{"multicore"} (not available on Windows) which spawns relatively lightweight
worker processes; and \code{"snow"} (the default) a slower and more resource intensive parallel worker option that
works on Windows systems.  The number of CPUs used by \code{\link{boot}} is by default set to one less than the number
of cores available, determined using the function \code{\link{detectCores}} in the \code{parallel} package.  If a
different number of CPUs should be used (for example if an existing cluster is passed using the \code{cl} parameter)
the number of CPUS can be changed by setting the option \code{"boot.ncpus"}.
}
\value{
An object of class \code{braidrm}, with the new elements:
\item{ciPass}{ TRUE if a sufficient proportion of bootstrapping trials successfully converged, FALSE otherwise }
\item{ciLevs}{ Two-element vector specifying the lower and upper percentiles of the desired confidence interval. Defaults
to 0.025 and 0.975 for a 95\% confidence interval. }
\item{ciVec}{ Assuming 'ciPass' is true, a vector containing the lower and upper bounds on the confidence intervals of
all free parameters }
\item{bCoefs}{ Array containing the best-fit parameter values from all bootstrapping trials.  Useful for calculating
confidence intervals on other values. }
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
brd <- braidrm(act~conc1+conc2,es8olatmz,getCIs=FALSE,fixed=c(NA,NA,NA,NA,1,NA,NA,-4,-4,-4))
summary(getBRAIDbootstrap(brd))
\dontrun{summary(getBRAIDbootstrap(brd,ciLevs=c(0.05,0.95)))}
}

\seealso{
\code{\link{braidrm}}, \code{\link{calcBRAIDconfint}}
}
\keyword{ regression }