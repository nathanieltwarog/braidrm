\name{findBestBRAID}
\alias{findBestBRAID}
\alias{findBestBRAID.default}
\alias{findBestBRAID.formula}
\title{ Select Best Fitting BRAID Surface Model }
\description{
Fits several BRAID surface models to the given data, and selects the most parsimonious model using the Akaike information criterion.
}
\usage{
\method{findBestBRAID}{default}(model, data, defaults, startparv=NULL, llims=NULL,
				ulims=NULL, itype=1, getCIs=TRUE, crossval=TRUE, ...)
\method{findBestBRAID}{formula}(model, data, ...)
}
\arguments{
  \item{model}{ a two-column array containing concentrations of Drug 1 and Drug 2 in each dose pair, or a symbolic formula (e.g.
\code{act ~ conc1+conc2}) specifying which variables are to be fit }
  \item{data}{ if \code{model} is an array, a vector of measurements of  response to the concentrations of Drug 1 and Drug 2;
if \code{model} is a formula, a data frame containing the columns specified in \code{formula} }
  \item{defaults}{ two-element vector specifying the default initial and maximal effects for the response surface.  These values
will be used in several of the models that are tried (see Details below).}
  \item{startparv}{ an optional parameter specifying starting parameter values for the optimization }
  \item{llims}{ a ten-element vector of lower limits on parameters being fit.   Any parameters that do not require a limit can
have a value of \code{NA}.  If \code{NULL} (the default), \code{llims} is calculated from the starting values in \code{startparv}
(or the values calculated for \code{startparv} if \code{startparv} is not specified). }
  \item{ulims}{ a vector of upper limits on parameters being fit.  Follows same behavior as \code{llims}. }
  \item{itype}{ an integer that specifies the type of interaction(s) that is assumed in the models.  The default is 1, which
assumes that the interaction parameter \eqn{\kappa}{kappa} is varying.  See details below for other possible values. }
  \item{getCIs}{ determines if confidence intervals will be calculated for all response surface parameters being fit.
Parameters are fit using a bootstrapping approach which resamples residuals and refits the new data. }
  \item{crossval}{ if \code{TRUE}, goodness of fit is determined by randomly assigning data points to four blocks, and evaluating
goodness of fit on each block by fitting the remaining three.  If \code{FALSE}, all data points are fit once, and goodness of fit
is determined using the residuals from that fit.  It is advisable not to use cross-validation when a relatively small number of
data points are available, especially on the margins (when drug A or drug B has concentration 0). }
  \item{\dots}{ Not used }
}
\details{
Because experiments do not reliably capture the full range of responses that a combination can produce, estimation of the initial
and maximal effect parameters for a BRAID surface can be highly unstable.  This function fits at least 10 distinct BRAID surface
models to the given data, and selects the model which best balances simplicity with goodness of fit.  For each interaction type
(see below), the following 10 models are fit:
\itemize{
  \item The initial effect parameter \eqn{E_0}{E0} varies freely, and the two maximal effect parameters \eqn{E_{f,A}}{EfA} and
\eqn{E_{f,B}}{EfB} vary freely independently of one another.
  \item The initial effect parameter \eqn{E_0}{E0} is fixed at the default value, and the two maximal effect parameters
\eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB} vary freely independently of one another.
  \item The initial effect parameter \eqn{E_0}{E0} the maximal effect parameter \eqn{E_{f,A}}{EfA} vary freely, and the
maximal effect parameter \eqn{E_{f,B}}{EfB} is fixed at the default value.
  \item The initial effect parameter \eqn{E_0}{E0} the maximal effect parameter \eqn{E_{f,B}}{EfB} vary freely, and the
maximal effect parameter \eqn{E_{f,A}}{EfA} is fixed at the default value.
  \item The initial effect parameter \eqn{E_0}{E0} varies freely, and the two maximal effect parameters \eqn{E_{f,A}}{EfA} and
\eqn{E_{f,B}}{EfB} are constrained to vary as a single parameter \eqn{E_f}{Ef}.
  \item The initial effect parameter \eqn{E_0}{E0} varies freely, and the two maximal effect parameters \eqn{E_{f,A}}{EfA} and
\eqn{E_{f,B}}{EfB} are fixed at the default value.
  \item The initial effect parameter \eqn{E_0}{E0} the two maximal effect parameter \eqn{E_{f,B}}{EfB} are fixed at the default
values, and the maximal effect parameter \eqn{E_{f,A}}{EfA} varies freely.
  \item The initial effect parameter \eqn{E_0}{E0} the two maximal effect parameter \eqn{E_{f,A}}{EfA} are fixed at the default
values, and the maximal effect parameter \eqn{E_{f,B}}{EfB} varies freely.
  \item The initial effect parameter \eqn{E_0}{E0} is fixed at the default value,, and the two maximal effect parameters
\eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB} are constrained to vary as a single parameter \eqn{E_f}{Ef}.
  \item The initial effect parameter \eqn{E_0}{E0} and the two maximal effect parameters \eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB}
are all fixed at the default values.
}
In all models, the potencies of the two drugs (represented by \eqn{{ID}_{M,A}}{IDMA} and \eqn{{ID}_{M,B}}{IDMB}) and the Hill
slopes of both drugs (represented by \eqn{n_a}{na} and \eqn{n_b}{nb}) vary freely.  Which of the interaction parameters
\eqn{\kappa} and \eqn{\delta} varies depends on the parameter \code{itype}, as follows:
\itemize{
\item \code{itype = 1}: \eqn{\kappa} varies freely in all models; \eqn{\delta} is fixed at 1 (10 models total).
\item \code{itype = 2}: \eqn{\delta} varies freely in all models; \eqn{\kappa} is fixed at 0 (10 models total).
\item \code{itype = 3}: Either \eqn{\kappa} or \eqn{\delta} (but not both) vary freely in all models (20 models total).
\item \code{itype = 4}: Either \eqn{\kappa} or \eqn{\delta} or both vary freely in all models (30 models total).
\item \code{itype = 5}: \eqn{\kappa} is fixed at 0 and \eqn{\delta} is fixed at 1 in all models (10 models total).
}
}
\value{
An object of the class 'braidrm', with elements as described in \code{\link{braidrm}}.
}
\author{
Nathaniel R. Twarog
}
\examples{
data(es8olatmz)
\dontrun{summary(findBestBRAID(cbind(es8olatmz$conc1,es8olatmz$conc2),
				es8olatmz$act,defaults=c(0,-2.7)))}
\dontrun{summary(findBestBRAID(act~conc1+conc2,es8olatmz,defaults=c(0,-2.7),itype=2))}
summary(findBestBRAID(act~conc1+conc2,es8olatmz,defaults=c(0,-4),getCIs=FALSE))
}

\seealso{ \code{\link{braidrm}}, \code{\link{getBRAIDbootstrap}}, \code{\link{runBRAIDanalysis}} }
\keyword{regression}