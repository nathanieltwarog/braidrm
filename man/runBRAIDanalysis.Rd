\name{runBRAIDanalysis}
\alias{runBRAIDanalysis}
\title{
Run Full BRAID Analysis
}
\description{
Performs a complete analysis of combined action data using the BRAID model.  The behaviors of both drugs alone (if such
doses are included) are fit using a Hill equation; if desired, the resulting Hill equation fits can be used to correct
uncertain or noisy input concentrations for use in combination analysis.  All measurements of the two drugs (in isolation
or in combination) are then fit using the BRAID equation.
}
\usage{
runBRAIDanalysis(data, defaults, llims = NULL, ulims = NULL, itype = 1,
					compounds = NULL, corrconc = FALSE, corrsigr = 1)
}
\arguments{
  \item{data}{ a data frame containing all measurements to be fit.  Must conatin columns named "conc1", "conc2", and "act"; other
columns may be required (see Details below). }
  \item{defaults}{ a two-element vector containing default values for the initial and maximal effects of the combination; used by
\code{\link{findBestBRAID}} to select best BRAID model }
  \item{llims}{ a ten-element vector of lower limits on parameters being fit.   Any parameters that do not require a limit can
have a value of \code{NA}. }
  \item{ulims}{ a vector of upper limits on parameters being fit.  Follows same behavior as \code{llims}. }
  \item{itype}{ an integer that specifies the type of interaction(s) that is assumed in the models.  Used by
\code{\link{findBestBRAID}} }
  \item{compounds}{ an optional parameter specifying which compounds are to be fit; this requires that compounds be identified in
the input data frame \code{data}. }
  \item{corrconc}{ boolean specifying whether input concentrations will be corrected according to individual dose response using
\code{\link{hillConcCorrect}}  }
  \item{corrsigr}{ the ratio of noise in measurement to noise the base-10 logarithm of concentration; used by
\code{\link{hillConcCorrect}} if \code{corrconc} is \code{TRUE} }
}
\details{
This convenience function is intended as a blueprint for a complete BRAID model analysis.  Though users of this package can
develop their own approaches to applying the BRAID model to combined action data (indeed we encourage it), this function
encapsulates our overall strategy for fitting the BRAID equation.  The input parameter \code{data} must contain all measurements
to be fit, with the columns "conc1", "conc2", and "act" containing the corresponding inputs and outputs.  If no other columns
required (based on the values of other parameters), "conc1" is assumed to describe the concentration of drug A, and "conc2" the
concentration of drug B, with values of 0 for "conc1" or "conc2" corresponding to conditions with Drug B or Drug A alone,
respectively.

If \code{data} contains measurements from multiple combinations, the identity of the compounds to be analyzed must be specified
in the parameter \code{compounds}.  Further, the input data frame \code{data} must contain columns named "compound1" and "compound2"
that specify which two compounds were used in each measurement.  Measurements will be used in which "compound1" is equal to the
first element of \code{compounds} or "conc1" is equal to 0, and in which "compound2" is equal to the second element of
\code{compounds} or "conc2 is equal to 0.  Note that this means measurements in which "conc1" and "conc2" are both equal to 0
will be used in any analysis.

For concentrations to be corrected in the combination analysis, there must be a correspondence between concentrations of each drug
alone and concentrations used in combination; in our implementation, this is assumed to be the result of a shared drugging plate,
such that measurements corresponding to the same well are given the same actual concentration of drug.  As a result, if
\code{corrconc} is set to \code{TRUE}, each measurement must be associated with a well specified in a column named "well".
}
\value{
A list with the following componenents:
\item{ concs }{ A two-column array containing the concentrations of the first and second drugs used in the analysis.  This will
only inlcude measurements from the input frame \code{data} that contained drug A alone, drug B alone, or drug A and B in
combination. Even if \code{corrconc} is \code{TRUE}, this array will contain the original, uncorrected concentrations. To access
the corrected concentrations, use the \code{conc1} and \code{conc2} elements of the component \code{braidFit} (described below). }
\item{ act }{ The measured response values corresponding to the concentration pairs in \code{concs}. }
\item{ corrconc }{ The value of the input parameter \code{corrconc}. }
\item{ corrsigr }{ The value of the input parameter \code{corrsigr}. }
\item{ braidFit }{ The best BRAID model fit as determined by \code{\link{findBestBRAID}}. }
\item{ hfit1 }{ A set of Hill model fits to the behavior of drug A in isolation, as output by \code{\link{findBestHill}}. If no
data points for drug A in isolation are included, this component is \code{NULL}. }
\item{ hfit2 }{ A set of Hill model fits to the behavior of drug B in isolation, as output by \code{\link{findBestHill}}. If no
data points for drug B in isolation are included, this component is \code{NULL}. }
}
\author{
Nathaniel R. Twarog
}
\examples{
data(es8olatmz)
\dontshow{
es8olatmz <- es8olatmz[(es8olatmz$plate==8&grepl("[AB]",es8olatmz$well))|es8olatmz$plate==9|
							(es8olatmz$plate==27&grepl("(01|02)",es8olatmz$well)),]
}# Note that 'es8olatmz' contains all necessary column names, including
# 'conc1','conc2','act','compound1','compound2', and, for concentration
# correction, 'well'
\dontrun{brdAnalysis <- runBRAIDanalysis(es8olatmz,defaults=c(0,-2.7),corrconc=TRUE)}
brdAnalysis <- runBRAIDanalysis(es8olatmz,defaults=c(0,-4))
summary(brdAnalysis$braidFit)
}

\seealso{\code{\link{braidrm}}, \code{\link{findBestBRAID}}, \code{\link{findBestHill}}, \code{\link{hillConcCorrect}}}
\keyword{regression}
