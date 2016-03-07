\name{invertBRAIDrsm}
\alias{invertBRAIDrsm}
\title{ Invert BRAID Response Function }
\description{
Determines which drug concentrations, when paired with the given input concentrations, produce the given response
values.  Useful for estimating potentiation of one compound by presence of the other.
}
\usage{
invertBRAIDrsm(val, DA=NULL, DB=NULL, parv)
}
\arguments{
  \item{val}{ Effect value or values to be inverted.  If a single value, it will be repeated for all values of
\code{DA} or \code{DB} which are input. }
  \item{DA}{ Concentrations of drug A to be inverted.  If NULL (the default), the concentrations in \code{DB} will
be inverted and the concentrations of drug A will be estimated. }
  \item{DB}{ Concentrations of drug B to be inverted.  If NULL (the default), the concentrations in \code{DA} will
be inverted and the concentrations of drug B will be estimated. }
  \item{parv}{ 10-element vector specifying the full set of parameters for the BRAID surface }
}
\value{
A vector of concentrations which, when paired with the given concentrations, produce the given values.  If \code{DA}
is NULL, the vector contains the concentrations of drug A corresponding to the concentrations of drug B given in
\code{DB}.  If \code{concB} is null, the vector contains the concentrations of drug B corresponding to the
concentrations of drug A give in \code{DA}.  Exactly one of these two inputs must be NULL.
}
\author{
Nathaniel R. Twarog
}
\examples{
level <- 90
doses <- c(0,10^-8,10^-7,10^-6)
ec90A <- invertBRAIDrsm(level,DB=doses,parv=c(10^-6,10^-6,1.5,1.5,1,1.8,0,100,100,100))
}

\seealso{
\code{\link{braidrm}}, \code{\link{evalBRAIDrsm}}
}
\keyword{ regression }
