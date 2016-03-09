\name{braidrm}
\alias{braidrm}
\alias{braidrm.default}
\alias{braidrm.formula}
\alias{print.braidrm}
\alias{summary.braidrm}
\alias{print.summary.braidrm}
\title{
BRAID Synergy Dose-Response Modeling
}
\description{
Calculates the best fit BRAID response surface given the concentrations of two drugs and the 
response of a measured target to the combination of those drugs.
}
\usage{
\method{braidrm}{default}(model,data,getCIs=TRUE,fixed="kappa2",startparv=NULL,llims=NULL,ulims=NULL,...)
\method{braidrm}{formula}(model,data,...)
\method{print}{braidrm}(x, ...)
\method{summary}{braidrm}(object, ...)
}
\arguments{
  \item{model}{ a two-column array containing concentrations of Drug 1 and Drug 2 in each dose pair, or a symbolic formula (e.g.
\code{act ~ conc1+conc2}) specifying which variables are to be fit }
  \item{data}{ if \code{model} is an array, a vector of measurements of  response to the concentrations of Drug 1 and Drug 2;
if \code{model} is a formula, a data frame containing the columns specified in \code{formula} }
  \item{getCIs}{ determines if confidence intervals will be calculated for all response surface parameters being fit.
Parameters are fit using a bootstrapping approach which resamples residuals and refits the new data. }
  \item{fixed}{ parameter specifying which parameters of the full BRAID model will be fit.  See Details for a full description
of this highly critical parameter }
  \item{startparv}{ an optional parameter specifying starting parameter values for the optimization.  The relationship between
\code{startparv} and \code{fixed} is rather subtle, and is discussed below in the Details section.}
  \item{llims}{ a vector of lower limits on parameters being fit.  May have length 10, or length equal to the number of
free parameters being fit.  Any parameters that do not require a limit can have a value of \code{NA}.  If \code{NULL} (the
default), \code{llims} is calculated from the starting values in \code{startparv} (or the values calculated for \code{startparv}
if \code{startparv} is not specified). }
  \item{ulims}{ a vector of upper limits on parameters being fit.  Follows same behavior as \code{llims}. }
  \item{x}{ the object of class "braidrm" to be printed }
  \item{object}{ the object of class "braidrm" to be summarized }
  \item{\dots}{ Not used }
}
\details{
This function is designed to give as much control as is reasonably possible over which parameters are optimized and how the
optimization behaves.  However, implementing this much control can be quite complicated, and despite our efforts to make the
function intuitive and transparent, the interface is still quite unwieldy.  A great deal can be accomplished with the function
simply using the model nicknames "kappa1" and "kappa2" (which constrain the two maximal effects to be equal or allow them to vary
independently, respectively); if these are sufficient for you, feel free to ignore the remainder of this section.

The parameter \code{fixed} is used to control which parameters are fit by the optimization, and in the case of the three maximal
effect parameters, how the parameters are constrained with respect to one another. \code{fixed} may be specified in one of three
forms: a raw vector, an index vector, or a model nickname.  There are seven model nicknames; using them is equivalent to using
the corresponding index vector as specified in the following table:
\tabular{ll}{
kappa1 \tab 1,2,3,4,6,7,10 \cr
kappa2 (default) \tab 1,2,3,4,6,7,8,9 \cr
kappa3 \tab 1,2,3,4,6,7,8,9,10 \cr
delta1 \tab 1,2,3,4,5,7,10 \cr
delta2 \tab 1,2,3,4,5,7,8,9 \cr
delta3 \tab 1,2,3,4,5,7,8,9,10 \cr
ebraid \tab 1,2,3,4,5,6,7,8,9,10 \cr
}
An index vector specifies which parameters vary in the optimization by listing their indices in the full ten parameter vector; a
raw vector specifies which parameters vary by setting the corresponding values equal to \code{NA}.  The remaining values in a raw
vector specify the values at which the fixed parameters are fixed, unless these values are overridden by \code{startparv}

\code{startparv} specifies the starting values for the optimization; one can input a vector in \code{startparv} that specifies
only the values of varying parameters, but the remaining values must be specified in a raw vector in \code{fixed}.  If a full-length
ten-element vector is provided for \code{starparv}, the values of fixed parameters are taken from that vector, regardless of the
type of input in \code{fixed}.

For parameters one through seven, the presence or absence of each parameter in the optimization is quite simple: either the parameter
is fixed at the specified (or calculated) value, or it varies between the specified (or calculated) optimization limits.  Parameters
eight through ten, however, which specify the maximal effect parameters, are more complicated.  The complication is introduced by the
fact that the parameters can be fixed \emph{or} constrained to be equal to one another, which introduce the same number of degrees of
freedom, but very different optimization behaviors.  How the different possibilities for \code{fixed} and \code{startparv} influence
the optimization behavior (for as many possibilities as we could think of) is described by the following table. In the table \eqn{W},
\eqn{X}, \eqn{Y}, and \eqn{Z} represent arbitrary but distinct valid values while ~ represents any valid (non-\code{NA}) value; in
addition, all \code{startparv} vectors are assumed to be complete ten-element vectors (incomplete vectors will simply be extended
with the corresponding values in \code{fixed}).
\tabular{lll}{
\bold{fixed} \tab \bold{startparv} \tab \bold{Behavior} \cr
Index: (\eqn{\ldots},8,9,10) or \tab Any \tab  All three maximal effect parameters vary independently,\cr
Raw: (\eqn{\ldots},\code{NA},\code{NA},\code{NA}) \tab \tab with \eqn{E_f}{Ef} constrained to be a larger effect than \eqn{E_{f,A}}{EfA} and \cr
\tab \tab \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},8,9) or \tab \code{NULL} \tab \eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB} vary independently and \eqn{E_f}{Ef} is \cr
Raw: (\eqn{\ldots},\code{NA},\code{NA},~) \tab \tab constrained to vary with the larger effect of \eqn{E_{f,A}}{EfA} and \cr
\tab \tab \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},8,9) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{X}) or \tab \eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB} vary independently and \eqn{E_f}{Ef} is \cr
Raw: (\eqn{\ldots},\code{NA},\code{NA},~) \tab (\eqn{\ldots}, \eqn{X}, \eqn{Y}, \eqn{X}) or \tab constrained to vary with the larger effect of \eqn{E_{f,A}}{EfA} and \cr
\tab (\eqn{\ldots}, \eqn{Y}, \eqn{X}, \eqn{X}) \tab \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},8,9) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{Y}) or \tab  \eqn{E_f}{Ef} is fixed at the value \eqn{Y}, while \eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB} vary \cr
Raw: (\eqn{\ldots},\code{NA},\code{NA},~) \tab (\eqn{\ldots}, \eqn{X}, \eqn{Z}, \eqn{Y}) \tab independently and are constrained to be smaller effects \cr
\tab \tab than \eqn{Y}. \cr
Index: (\eqn{\ldots},10) or \tab \code{NULL} \tab All maximal parameters are constrained to vary as a \cr
Raw: (\eqn{\ldots},\eqn{X},\eqn{X},\code{NA}) \tab \tab single parameter. \cr
Raw: (\eqn{\ldots},\eqn{X},\eqn{Y},\code{NA}) \tab \code{NULL} \tab If \eqn{X} is larger than \eqn{Y}, \eqn{E_f}{Ef} and \eqn{E_{f,A}}{EfA} are constrained to \cr
\tab \tab vary as a single parameter above \eqn{Y}, and \eqn{E_{f,B}}{EfB} is fixed\cr
\tab \tab  at \eqn{Y}. Otherwise, the roles of A and B are reversed. \cr
Index: (\eqn{\ldots},10) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{X}) \tab All maximal parameters are constrained to vary as a \cr
Raw: (\eqn{\ldots},~,~,\code{NA}) \tab \tab single parameter. \cr
Index: (\eqn{\ldots},10) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{Y}, \eqn{X}) \tab \eqn{E_f}{Ef} and \eqn{E_{f,A}}{EfA} are constrained to vary as a single \cr
Raw: (\eqn{\ldots},~,~,\code{NA}) \tab \tab parameter above \eqn{Y}, and \eqn{E_{f,B}}{EfB} is fixed at \eqn{Y}. \cr
Index: (\eqn{\ldots},10) or \tab (\eqn{\ldots}, \eqn{Y}, \eqn{X}, \eqn{X}) \tab \eqn{E_f}{Ef} and \eqn{E_{f,B}}{EfB} are constrained to vary as a single \cr
Raw: (\eqn{\ldots},~,~,\code{NA}) \tab \tab parameter above \eqn{Y}, and \eqn{E_{f,A}}{EfA} is fixed at \eqn{Y}. \cr
Index: (\eqn{\ldots},10) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{Y}, \eqn{Z}) \tab \eqn{E_f}{Ef} is constrained to vary above the larger effect of \eqn{X} \cr
Raw: (\eqn{\ldots},~,~,\code{NA}) \tab \tab and \eqn{Y}, \eqn{E_{f,A}}{EfA} is fixed at \eqn{X}, and \eqn{E_{f,B}}{EfB} is fixed at \eqn{Y}. \cr
Index: (\eqn{\ldots},8,10) or \tab \code{NULL} or \tab  \eqn{E_{f,A}}{EfA} varies freely, and \eqn{E_{f,B}}{EfB} and \eqn{E_f}{Ef} are constrained to \cr
Raw: (\eqn{\ldots},\code{NA},~,\code{NA}) \tab (\eqn{\ldots}, ~, \eqn{X}, \eqn{X}) \tab vary as a single parameter with a larger effect than \eqn{E_{f,A}}{EfA}. \cr
Index: (\eqn{\ldots},8,10) or \tab (\eqn{\ldots}, ~, \eqn{X}, \eqn{Y}) \tab \eqn{E_{f,A}}{EfA} varies freely, \eqn{E_{f,B}}{EfB} is fixed at \eqn{X} and \eqn{E_f}{Ef} varies \cr
Raw: (\eqn{\ldots},\code{NA},~,\code{NA}) \tab \tab freely constrained to be a larger effect than \eqn{E_{f,A}}{EfA} and \cr
\tab \tab \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},9,10) or \tab \code{NULL} or \tab \eqn{E_{f,B}}{EfB} varies freely, and \eqn{E_{f,A}}{EfA} and \eqn{E_f}{Ef} are constrained to \cr
Raw: (\eqn{\ldots},~,\code{NA},\code{NA}) \tab (\eqn{\ldots}, \eqn{X}, ~, \eqn{X}) \tab vary as a single parameter with a larger effect than \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},9,10) or \tab (\eqn{\ldots}, \eqn{X}, ~, \eqn{Y}) \tab \eqn{E_{f,B}}{EfB} varies freely, \eqn{E_{f,A}}{EfA} is fixed at \eqn{X} and \eqn{E_f}{Ef} varies \cr
Raw: (\eqn{\ldots},~,\code{NA},\code{NA}) \tab \tab freely constrained to be a larger effect than \eqn{E_{f,A}}{EfA} and \cr
\tab \tab \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},8) or \tab \code{NULL} \tab \eqn{E_{f,A}}{EfA} varies freely, \eqn{E_{f,B}}{EfB} is fixed at \eqn{X} (or calculated \cr
Raw: (\eqn{\ldots},\code{NA},\eqn{X},\eqn{X}) \tab \tab starting value), and \eqn{E_f}{Ef} is constrained to vary with \cr
\tab \tab the larger effect of \eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB}. \cr
Raw: (\eqn{\ldots},\code{NA},\eqn{X},\eqn{Y}) \tab \code{NULL} \tab \eqn{E_{f,B}}{EfB} is fixed at \eqn{X}, \eqn{E_f}{Ef} is fixed at \eqn{Y}, and \eqn{E_{f,A}}{EfA} varies \cr
\tab \tab freely constrained to be a smaller effect than \eqn{Y}. \cr
Index: (\eqn{\ldots},8) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{X}) or \tab \eqn{E_{f,A}}{EfA} varies freely, \eqn{E_{f,B}}{EfB} is fixed at \eqn{X}, and \eqn{E_f}{Ef} is \cr
Raw: (\eqn{\ldots},\code{NA},~,~) \tab (\eqn{\ldots}, \eqn{Y}, \eqn{X}, \eqn{X}) or \tab constrained to vary with the larger effect of \eqn{E_{f,A}}{EfA} \cr
\tab (\eqn{\ldots}, \eqn{Y}, \eqn{X}, \eqn{Y}) \tab and \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},8) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{Y}) or \tab \eqn{E_{f,B}}{EfB} is fixed at \eqn{X}, \eqn{E_f}{Ef} is fixed at \eqn{Y}, and \eqn{E_{f,A}}{EfA} varies \cr
Raw: (\eqn{\ldots},\code{NA},~,~) \tab (\eqn{\ldots}, \eqn{Z}, \eqn{X}, \eqn{Y}) \tab freely constrained to be a smaller effect than \eqn{Y}. \cr
Index: (\eqn{\ldots},9) or \tab \code{NULL} \tab \eqn{E_{f,B}}{EfB} varies freely, \eqn{E_{f,A}}{EfA} is fixed at \eqn{X} (or calculated \cr
Raw: (\eqn{\ldots},\eqn{X},\code{NA},\eqn{X}) \tab \tab starting value), and \eqn{E_f}{Ef} is constrained to vary with \cr
\tab \tab the larger effect of \eqn{E_{f,A}}{EfA} and \eqn{E_{f,B}}{EfB}. \cr
Raw: (\eqn{\ldots},\eqn{X},\code{NA},\eqn{Y}) \tab \code{NULL} \tab \eqn{E_{f,A}}{EfA} is fixed at \eqn{X}, \eqn{E_f}{Ef} is fixed at \eqn{Y}, and \eqn{E_{f,B}}{EfB} varies \cr
\tab \tab freely constrained to be a smaller effect than \eqn{Y}. \cr
Index: (\eqn{\ldots},9) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{X}) or \tab \eqn{E_{f,B}}{EfB} varies freely, \eqn{E_{f,A}}{EfA} is fixed at \eqn{X}, and \eqn{E_f}{Ef} is \cr
Raw: (\eqn{\ldots},~,\code{NA},~) \tab (\eqn{\ldots}, \eqn{X}, \eqn{Y}, \eqn{X}) or \tab constrained to vary with the larger effect of \eqn{E_{f,A}}{EfA} \cr
\tab (\eqn{\ldots}, \eqn{X}, \eqn{Y}, \eqn{Y}) \tab and \eqn{E_{f,B}}{EfB}. \cr
Index: (\eqn{\ldots},9) or \tab (\eqn{\ldots}, \eqn{X}, \eqn{X}, \eqn{Y}) or \tab \eqn{E_{f,A}}{EfA} is fixed at \eqn{X}, \eqn{E_f}{Ef} is fixed at \eqn{Y}, and \eqn{E_{f,B}}{EfB} varies \cr
Raw: (\eqn{\ldots},~,\code{NA},~) \tab (\eqn{\ldots}, \eqn{X}, \eqn{Z}, \eqn{Y}) \tab freely constrained to be a smaller effect than \eqn{Y}. \cr
Index: (\eqn{\ldots}) or \tab Any \tab All three values are fixed, and do not vary. \cr
Raw: (\eqn{\ldots},~,~,~) \tab \tab 
}
}
\value{
An object of class 'braidrm', including the elements
\item{conc1}{ Concentrations of drug 1 used in the fit }
\item{conc2}{ Concentrations of drug 2 used in the fit }
\item{act}{ Reponse to drugs 1 and 2 used in the fit }
\item{fitted.values}{ Value of best-fit response surface at concentrations of drugs 1 and 2 used in the fit }
\item{residuals}{ Difference between actual measurements and fitted values }
\item{ostart}{ The full 10-parameter starting parameter vector used in the optimization }
\item{fixed}{ Vector describing which parameters were fit }
\item{mlims}{ Array containing upper and lower parameter bounds used in the optimization }
\item{coefficients}{ Estimates of fitted parameters }
\item{fullpar}{ Full parameter vector including best-fit parameters }
\item{convergence}{ From \code{optim}, indicates if optimization successfully converged }
\item{message}{ From \code{optim}, describes results of non-linear optimization }
\item{call}{ Function call }
If \code{getCIs} is TRUE, then the following elements are also included
\item{ciPass}{ TRUE if a sufficient proportion of bootstrapping trials successfully converged, FALSE otherwise }
\item{ciLevs}{ Two-element vector specifying the lower and upper percentiles of the desired confidence interval. Defaults
to 0.025 and 0.975 for a 95\% confidence interval. }
\item{ciVec}{ Assuming 'ciPass' is true, a vector containing the lower and upper bounds on the confidence intervals of
all free parameters }
\item{bCoefs}{ A matrix containing the best-fit parameter values from all bootstrapping trials.  Useful for calculating
confidence intervals on other values, such as EC99 }
}
\author{ Nathaniel R. Twarog }
\examples{
data(es8olatmz)
summary(braidrm(act~conc1+conc2,es8olatmz,getCIs=FALSE))
\dontrun{
summary(braidrm(cbind(es8olatmz$conc1,es8olatmz$conc2),es8olatmz$act))
summary(braidrm(act~conc1+conc2,es8olatmz,fixed="delta2"))
summary(braidrm(act~conc1+conc2,es8olatmz,fixed=c(1,2,3,4,6,8,9)))
summary(braidrm(act~conc1+conc2,es8olatmz,llims=c(NA,NA,NA,NA,NA,NA,NA,-4,-4,-4)))
}
}
\seealso{ \code{\link{getBRAIDbootstrap}}, \code{\link{calcBRAIDconfint}}, \code{\link{evalBRAIDrsm}} }
\keyword{ regression }