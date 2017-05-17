\name{braidrm-package}
\alias{braidrm-package}
\docType{package}
\title{
Fitting Dose Response with the BRAID Combined Action Model
}
\description{
Contains functions necessary for evaluating, analyzing, and fitting combined action dose response surfaces with the
Bivariate Response to Additive Interacting Dose (BRAID) model of combined action.
}
\details{
\tabular{ll}{
Package: \tab braidrm\cr
Type: \tab Package\cr
Version: \tab 0.82\cr
Date: \tab 2017-05-16\cr
License: \tab GPL (>= 3)\cr
}
A mathematical description of the BRAID surface model can be found in \code{\link{evalBRAIDrsm}}. The heart of the package is the
function \code{\link{braidrm}}, which takes a set of paired concentrations and measured responses and fits a particular
instantiation of the BRAID model to the data.  Bootstrapped confidence intervals on all fit parameters can be constructed by
\code{\link{braidrm}}, or afterwards by the function \code{\link{getBRAIDbootstrap}}. Known BRAID surfaces can be evaluated or
inverted using the calculating functions \code{\link{evalBRAIDrsm}} and \code{\link{invertBRAIDrsm}}. As analysis of combined dose
response is closely linked to single agent dose response, we have provided several simple functions for analyzing such
relationships, such as \code{\link{findBestHill}}.  Finally, we also provide the more complete analytic functions
\code{\link{findBestBRAID}} and \code{\link{runBRAIDanalysis}}, which perform a series of BRAID surface fits (and in the case of
\code{\link{runBRAIDanalysis}}, single agent fits) to simplify the process of implementing a BRAID analysis pipeline.
}
\author{
Nathaniel R. Twarog <nathaniel.twarog@stjude.org>
}
\keyword{ package }
