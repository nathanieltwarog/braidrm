\name{findBestHill}
\alias{findBestHill}
\alias{findBestHill.default}
\alias{findBestHill.formula}
\title{
Fit Hill Model with Model Selection
}
\description{
Fits the given concentration and response data with several Hill (or log-likelihood) models, and selects the model that best
balances simplicity and accuracy according to the Akaike Information Criterion.
}
\usage{
\method{findBestHill}{default}(conc,act,defaults,startparv=NULL,llims=NULL,ulims=NULL,...)
\method{findBestHill}{formula}(formula,data,...)
}
\arguments{
  \item{conc}{ a vector of concentrations }
  \item{act}{ a vector of measurements of response to the concentrations set down in \code{conc} }
  \item{formula}{ a symbolic formula (e.g. \code{act ~ conc}) specifying which variables are to be fit }
  \item{data}{ a data frame containing the columns specified in \code{formula} }
  \item{defaults}{ two-element vector specifying the default initial and maximal effects for the Hill model.  These values will
be used in several of the models that are tried (see Details below).}
  \item{startparv}{ an optional parameter specifying starting parameter values for the optimization.  Any parameters that do not
have a starting value can have a value of NA (and will be estimated from the concentrations and responses.}
  \item{llims}{ a four-element vector of lower limits on parameters being fit.  Any parameters that do not require a limit can
have a value of NA.  If NULL (the default), \code{llims} is calculated from the starting values in \code{startparv} (or the values
calculated for \code{startparv} if \code{startparv} is not specified). }
  \item{ulims}{ a vector of upper limits on parameters being fit.  Follows same behavior as \code{llims}. }
  \item{\dots}{ Not used }
}
\details{
As many dose response experiments fail to capture the full sigmoidal behavior of a compound, estimation of the initial and maximal
effect parameters for a four-parameter Hill equation (see \code{\link{evalHillEqn}} for a description of the equation used in this
package) can be quite unstable.  To avoid estimation of wildly divergent parameters and overfitting, this function fits a
four-parameter Hill equation to the given data four times, with the initial and maximal effect parameters either varying freely or
fixed to the provided default values.  The most parsimonious model is selected from the four fits according to the Akaike Information
Criterion.  The model also fits a constant flat function to the data; though this fifth "model" is not included in the model
selection, it is included in the output for reference, as a check that the given data shows any signal at all.
}
\value{
A list describing the results of the model fits, with the following components:
\item{ conc }{The vector of concentrations used in the fit}
\item{ act }{The vector of responses used in the fit}
\item{ bestModIdx }{The index of the model among those fit which yielded the most parsimonious fit}
\item{ bestModName }{The name of the best-fit model.  The names of the models are
\enumerate{
  \item "m2p": A model in which both the initial and maximal effect parameters are fixed to the default values
  \item "m3puc": A model in which the initial effect parameter varies freely
  \item "m3plc": A model in which the maximal effect parameter varies freely
  \item "m4p": A model in which both the initial and maximal effect parameters vary freely
}
The fifth "model", a constant linear fit, is named "Lin", but will not be selected by the model selection.
}
\item{ mlims }{A matrix of upper and lower limits applied to the four Hill parameters in the fit}
\item{ allfits }{A list of the five fits (including the trivial linear fit) applied to the data.  Each fit is described by a list
with the following components:
\itemize{
  \item coefficients: The full four-element parameter vector describing the best-fit Hill equation in this model (or the single
parameter for the constant linear fit)
  \item parv: The vector of best-fit values for only those parameters which varied freely in the model
  \item pinds: A vector of parameter indices indicating which of the four Hill parameters varied in this model
  \item AIC: The Akaike Information Criterion value for this model
}}
}
\author{
Nathaniel R. Twarog
}
\examples{
data(es8olatmz)
drv <- es8olatmz$compound2=="DMSO"
hfit <- findBestHill(es8olatmz$conc1[drv],es8olatmz$act[drv],defaults=c(0,-2.7))
coef(hfit$allfits[[hfit$bestModIdx]])
hfit <- findBestHill(act~conc1,es8olatmz[drv,],defaults=c(0,-4))
coef(hfit$allfits[[hfit$bestModIdx]])
}

\seealso{
\code{\link{evalHillEqn}}, \code{\link{getHillBootstrap}}, \code{\link{runBRAIDanalysis}}
}
\keyword{regression}
