\name{es1plates}
\docType{data}
\alias{es1plates}
\alias{es1extplates}
\alias{es8plates}
\alias{es8extplates}
\alias{ew8plates}
\title{Plate Information from Ewings-Sarcoma Experiments}
\description{
  These data frames describe the compounds used in each of plates tested
  in a series of cytotoxicity experiments on cell-lines from the pediatric
  cancer Ewings' sarcoma.  Data frames describe plates from several different
  experiments or sets of experiments:
  \describe{
    \item{es1plates}{Experiment combining three PARP inhibitors with standard-of-care DNA damaging
    agents in the ES1 cell line.}
    \item{es1extplates}{Follow up experiments combining BMN 673 and SN-38 with several additional
    mechanistic classes in the ES1 cell line.}
    \item{es8plates}{Experiment combining three PARP inhibitors with standard-of-care DNA damaging
    agents in the ES8 cell line.}
    \item{es8extplates}{Follow up experiments combining BMN 673 and SN-38 with several additional
    mechanistic classes in the ES8 cell line.}
    \item{ew8plates}{Experiment combining three PARP inhibitors with standard-of-care DNA damaging
    agents in the EW8 cell line.}
  }
}
\usage{data(ewings)}
\format{Data frames containing the following columns:
  \describe{
    \item{plate}{The number of a plate run in a particular experiment}
    \item{compound1}{The compound drugged in decreasing concentration across the plate (DMSO if none)}
    \item{compound2}{The compound drugged in decreasing concentration down the plate (DMSO if none)}
    \item{posmean}{The average base-10 logarithm of CellTiter-Glo intensity in positive control wells
    on the plate, normalized to intensity in negative control wells on the same plate}
    \item{date}{The date on which the experiment using that plate was run}
  }
}
\keyword{datasets}