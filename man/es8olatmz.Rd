\name{es8olatmz}
\docType{data}
\alias{es8olatmz}
\title{Results from Ewings-Sarcoma Experiments}
\description{
  This data frames contain the results of a single combination of the PARP inhibitor
  olaparib and the DNA-damaging agent temozolomide in a cytotoxicity experiment on the
  Ewings' sarcoma cell line ES8 run on April 4, 2014. 
}
\usage{data(es8olatmz)}
\format{A data frame containing the following columns:
  \describe{
    \item{plate}{The number of the plate this measurement was taken from}
    \item{well}{The name of the well (e.g. "A01") this measurement was taken from}
    \item{compound1}{The compound drugged in decreasing concentration across this plate (DMSO if none)}
    \item{compound2}{The compound drugged in decreasing concentration down this plate (DMSO if none)}
    \item{conc1}{The concentration (in molar) of the first compound in this well}
    \item{conc2}{The concentration (in molar) of the second compound in this well}
    \item{act}{The base-10 logarithm of CellTiter-Glo intensity in this well normalized to intensity
    in negative control wells on the same plate}
    \item{date}{The date on which the experiment using this plate was run}
  }
}
\keyword{datasets}