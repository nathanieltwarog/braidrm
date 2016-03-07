\name{getHillBootstrap}
\alias{getHillBootstrap}
\title{
Bootstrap Hill Equation Parameter Confidence Intervals
}
\description{
Uses sampling residuals with replacement to bootstrap a distribution of Hill equation parameters for the construction
of confidence intervals on those parameters.
}
\usage{
getHillBootstrap(hfitl, ciLevs = c(0.025,0.975), mi = NULL)
}
\arguments{
  \item{hfitl}{a list of Hill equation model fits produced by \code{\link{findBestHill}}}
  \item{ciLevs}{a two-element vector specifying the lower and upper percentiles of the desired confidence intervals. Default
is \code{c(0.025,0.975)}, which yields a 95\% confidence interval.}
  \item{mi}{index of the desired model to bootstrap and construct confidence intervals on.  If \code{NULL} (the default), the
value of \code{bestModIdx} in \code{hfitl} will be used.}
}
\value{
A list of the same form as produced by \code{\link{findBestHill}}, with the following additional components:
\item{ciPass}{TRUE if at least 50\% of bootstrapping trials result in a converged fit}
\item{ciLevs}{Equal to the input parameter \code{ciLevs}}
\item{ciMInd}{Index of the model used to construct the bootstrapped confidence interval}
\item{ciVec}{Vector of lower and upper bounds of confidence intervals on all freely varying parameters, in order}
\item{bCoefs}{Four column array containing complete Hill equation parameters for all resampled fits}
}
\author{
Nathaniel R. Twarog
}
\examples{
data(es8olatmz)
drv <- es8olatmz$compound2=="DMSO"
hfit <- findBestHill(act~conc1,es8olatmz[drv,],defaults=c(0,-2.7))
hfit1 <- getHillBootstrap(hfit)

# Setting 'mi' to 4 constructs a confidence interval on the four-
# parameter Hill model, overriding the 'best' model index
hfit2 <- getHillBootstrap(hfit,ciLevs=c(0.05,0.95),mi=4)
}

\seealso{\code{\link{findBestHill}}, \code{\link{evalHillEqn}}}
\keyword{regression}
