\name{evalHillEqn}
\alias{evalHillEqn}
\alias{invertHillEqn}
\title{
Calculate and Invert Hill Equation
}
\description{
Calculates the response of a four-parameter Hill (or log-logistic) dose-response model, or determines the concentrations which
yield the given response.
}
\usage{
evalHillEqn(conc, parv)
invertHillEqn(val, parv)
}
\arguments{
  \item{conc}{ a vector of concentrations to be passed as input into the model}
  \item{parv}{ a four-element vector specifying the full parameter set of a particular Hill model }
  \item{val}{ a vector of responses to be inverted in a particular Hill model }
}
\details{
There is some ambiguity in how the parameters of a Hill or log-logistic model are specified.  We have chosen to use the following
equation:
\deqn{E(D) = E_0 + \frac{E_f-E_0}{1+\left(\frac{D}{{ID}_M}\right)^{-n}}}{E(D) = E0 + (Ef-E0)/(1+(D/IDM)^(-n))}
where the Hill slope \eqn{n} is always positive.  Using this equation, the parameter vector for a Hill model is
(\eqn{n},\eqn{E_0}{E0},\eqn{E_f}{Ef},\eqn{ln{{ID}_M}}{ln(IDM)}) (note that the fourth parameter value is log-transformed in the
parameter vector).
}
\value{
For \code{evalHillEqn} a vector of responses resulting from the given concentrations.  For \code{invertHillEqn}, a vector
of concentrations that produce the given responses; responses beyond the model's maximal effect will produce a value of \code{Inf},
while responses that do not reach the model's baseline effect produce a value of 0.
}
\author{
Nathaniel R. Twarog
}
\examples{
act <- evalHillEqn(exp(seq(log(10^-8),log(10^-5),length=100)),parv=c(2,0,100,log(10^-6)))
ec90 <- invertHillEqn(90,parv=c(2,0,100,log(10^-6)))
}

\seealso{\code{\link{findBestHill}}, \code{\link{hillConcCorrect}}}
\keyword{regression}
