# data.R

#' Example Additive Surface
#'
#' A synthetically generated response surface using a additive parameter vector.
#' The surface was generated with IDMA of 1, IDMB of 1, na of 3, nb of 3,
#' kappa of 0, E0 of 0, EfA of 1, EfB of 1, and Ef of 1. Every pair of
#' concentrations is sampled once, with concentrations of 0 and a seven-point
#' two-fold dilution from 0.125 to 8. "Measurements" were sampled from a
#' normal noise distribution around ground truth values with a standard
#' deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"additiveExample"

#' Example Synergistic Surface
#'
#' A synthetically generated response surface using a synergistic parameter vector.
#' The surface was generated with IDMA of 1, IDMB of 1, na of 3, nb of 3,
#' kappa of 2, E0 of 0, EfA of 1, EfB of 1, and Ef of 1. Every pair of
#' concentrations is sampled once, with concentrations of 0 and a seven-point
#' two-fold dilution from 0.125 to 8. "Measurements" were sampled from a
#' normal noise distribution around ground truth values with a standard
#' deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"synergisticExample"

#' Example Antagonistic Surface
#'
#' A synthetically generated response surface using a antagonistic parameter vector.
#' The surface was generated with IDMA of 1, IDMB of 1, na of 3, nb of 3,
#' kappa of -1, E0 of 0, EfA of 1, EfB of 1, and Ef of 1. Every pair of
#' concentrations is sampled once, with concentrations of 0 and a seven-point
#' two-fold dilution from 0.125 to 8. "Measurements" were sampled from a
#' normal noise distribution around ground truth values with a standard
#' deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"antagonisticExample"

#' Example Partial or Incomplete Surface
#'
#' A synthetically generated response surface using parameter vector describing
#' a surface with one only barely detectable effect. The surface was generated
#' with IDMA of 1, IDMB of 100, na of 3, nb of 3, kappa of 0, E0 of 0, EfA of 1,
#' EfB of 0.1, and Ef of 1. Every pair of concentrations is sampled once, with
#' concentrations of 0 and a seven-point two-fold dilution from 0.125 to 8.
#' "Measurements" were sampled from a normal noise distribution around ground
#' truth values with a standard deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"incompleteExample"

#' Example Protective Surface
#'
#' A synthetically generated response surface using a flipped "protective"
#' parameter vector. The surface was generated with IDMA of 0.5, IDMB of 2, na
#' of 3, nb of 3, kappa of 2, E0 of 0, EfA of 1, EfB of 0, and Ef of 0; the
#' surface was flipped along the drug A axis (so `flip` was set to "A"). Every
#' pair of concentrations is sampled once, with concentrations of 0 and a
#' seven-point two-fold dilution from 0.125 to 8. "Measurements" were sampled
#' from a normal noise distribution around ground truth values with a standard
#' deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"protectiveExample"

#' Example Oppositional Surface
#'
#' A synthetically generated response surface using a flipped "oppositional"
#' parameter vector. The surface was generated with IDMA of 1, IDMB of 1, na
#' of 3, nb of 3, kappa of -0.5, E0 of 0, EfA of 1, EfB of -0.5, and Ef of 1;
#' the surface was flipped along the drug B axis (so `flip` was set to "B").
#' Every pair of concentrations is sampled once, with concentrations of 0 and a
#' seven-point two-fold dilution from 0.125 to 8. "Measurements" were sampled
#' from a normal noise distribution around ground truth values with a standard
#' deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"oppositionalExample"

#' Example Coactive Surface
#'
#' A synthetically generated response surface using a flipped "coactive"
#' parameter vector. The surface was generated with IDMA of 0.5, IDMB of 0.5, na
#' of 3, nb of 3, kappa of 0, E0 of 0, EfA of 0, EfB of 0, and Ef of 1;
#' the surface was flipped along both drug axes (so `flip` was set to "both").
#' Every pair of concentrations is sampled once, with concentrations of 0 and a
#' seven-point two-fold dilution from 0.125 to 8. "Measurements" were sampled
#' from a normal noise distribution around ground truth values with a standard
#' deviation of 0.07.
#'
#' @format
#' A data frame with 64 rows and 4 columns
#' \describe{
#'   \item{concA}{The concentration of drug A}
#'   \item{concB}{The concentration of drug B}
#'   \item{truth}{The true response surface value at the given dose pair}
#'   \item{measure}{The sampled noisy measurement of the response surface at
#'   the given dose pair}
#' }
"coactiveExample"
