

#' Braid kappa Bayesian Prior
#'
#' Generates a Bayesian prior object on the BRAID parameter kappa to stabilize
#' parameter fitting
#'
#' @param spread Rough estimate of the standard deviation of measurement noise
#' or errors expected in a given data set.  Commonly used values are standard
#' deviation of negative/positive controls or root mean squared error of a
#' preliminary surface fit.
#' @param strength String indicating the influence of the BRAID prior on the
#' resulting fit.  Must be one of "mild", "moderate" (the default), "high", or
#' "none".
#'
#' @return An object of class `kappaPrior` containing two numeric elements,
#' `spread`, and `strength`.  Used in BRAID fitting functions to stabilize
#' the parameter kappa
#' @export
#'
#' @examples
#' prior <- kappaPrior(0.05,"mild")
#'
#' bfit <- braidrm(measure ~ concA + concB, incompleteExample, prior=prior)
#' summary(bfit)
kappaPrior <- function(spread,strength="moderate") {
	if (is.character(strength)) {
		strength <- switch(strength,
			none = 0,
			mild = 2/3,
			moderate = 1,
			high = 3/2,
			stop(sprintf("Unrecognized prior strength '%s'.",strength))
		)
	}
	structure(
		list(spread=spread,strength=strength),
		class = "kappaPrior"
	)
}

defaultSurfaceSpread <- function(concs,act,weights=NULL,start=NULL) {
	if (length(act)<=9) { return(stats::sd(act)) }
	model <- 1:9
	if (is.null(weights)) {
		weights <- rep(1,length(act))
	}
	if (isBraidParameter(start)) {
		start <- fillOutBraidPar(start)
	} else {
		start <- defaultStartingVector(concs,act,model,"",start,0)
	}
	pbounds <- fillParameterBounds(NULL,NULL,model,concs,start)
	bfit <- fitBraidScenario_IV_1(concs,act,model,weights,start,0,pbounds,0)
	spread <- sqrt(sum(bfit$residuals^2)/(length(act)-9))
	spread
}
