\name{es1data}
\docType{data}
\alias{es1data}
\alias{es1extdata}
\alias{es8data}
\alias{es8extdata}
\alias{ew8data}
\title{Results from Ewings-Sarcoma Experiments}
\description{
  These data frames contain the results of a series of cytotoxicity experiments
  on cell-lines from the pediatric cancer Ewings' sarcoma.  Data frames describe
  results from several different experiments or sets of experiments:
  \describe{
    \item{es1data}{Experiment combining three PARP inhibitors with standard-of-care DNA damaging
    agents in the ES1 cell line.}
    \item{es1extdata}{Follow up experiments combining BMN 673 and SN-38 with several additional
    mechanistic classes in the ES1 cell line.}
    \item{es8data}{Experiment combining three PARP inhibitors with standard-of-care DNA damaging
    agents in the ES8 cell line.}
    \item{es8extdata}{Follow up experiments combining BMN 673 and SN-38 with several additional
    mechanistic classes in the ES8 cell line.}
    \item{ew8data}{Experiment combining three PARP inhibitors with standard-of-care DNA damaging
    agents in the EW8 cell line.}
  }
}
\usage{data(ewings)}
\format{Data frames containing the following columns:
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