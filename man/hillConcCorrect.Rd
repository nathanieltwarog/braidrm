\name{hillConcCorrect}
\alias{hillConcCorrect}
\title{
Hill-Based Concentration Correction
}
\description{
Estimates actual underlying concentrations leading to a given set of response measurements based on the assumption that
actual concentrations are log-normally distributed around target concentrations, response errors are normally distributed, and
the actual underlying relationship between concentration and response is represented by the given Hill dose-response model.
}
\usage{
hillConcCorrect(conc, act, parv, sigr = 1)
}
\arguments{
  \item{conc}{ a vector of expected or target concentrations, around which actual concentrations are assumed to be log-normally
distributed }
  \item{act}{ a vector of response values }
  \item{parv}{ a four-parameter vector specifying a Hill model as described in \code{\link{evalHillEqn}} which is assumed to
be the actual relationship between concentration and response  }
  \item{sigr}{ the estimated ratio of the noises in response- and log (base10) concentration-space }
}
\details{
Suppose that \eqn{\hat{c}}{c'} is a given target concentration, and \eqn{c} is the actual concentration in given well, plate, or
condition.  Suppose also that \eqn{y} is the actual response that would result from the concentration in the given Hill
dose-response model, and \eqn{\hat{y}}{y'} is the measured response value.  This function assumes that
\deqn{\hat{y} \sim N(y,\sigma)}{y' ~ N(y,sigma)}
\deqn{\log_{10} c \sim N(\log_{10}\hat{c},\frac{\sigma}{r})}{log10(c) ~ N(log10(c'),sigma/r)}
for some \eqn{\sigma}{sigma}, where \eqn{N} is a normal distribution, and \eqn{r} is the ratio specified by the parameter \code{sigr}.
Based on these assumptions, the function uses Bayes' rule to calculate the maximum likelihood estimate of \eqn{c} for every given
value of \eqn{\hat{c}}{c'} and \eqn{\hat{y}}{y'}.
}
\value{
A vector of concentrations representing the maximum likelihood estimates of the actual concentrations which produced the
given responses.
}
\author{
Nathaniel R. Twarog
}
\examples{
data(es8olatmz)
drv <- es8olatmz$compound2=="DMSO"
hfit <- findBestHill(act~conc1,es8olatmz[drv,],defaults=c(0,-2.7))
drvpos <- drv & es8olatmz$conc1>0
cconc <- hillConcCorrect(es8olatmz$conc1[drvpos],es8olatmz$act[drvpos],
			coef(hfit$allfits[[hfit$bestModIdx]]),sigr=1)
}
\seealso{\code{\link{evalHillEqn}}, \code{\link{runBRAIDanalysis}}}
\keyword{regression}
