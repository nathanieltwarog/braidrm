

#' BRAID Parameter Confidence Intervals
#'
#' Uses residuals-based bootstrapping to estimate confidence intervals on a
#' BRAID fit's response surface parameters
#'
#' @param bfit A BRAID fit object of class `braidrm`.  If this object already
#' has bootstrapped coefficients, a warning will be given, and they will be
#' overwritten
#' @param ciLevs The lower and upper quantiles at which the confidence intervals
#' should be estimated.  Default is 0.025 and 0.975, producing 95% confidence
#' intervals
#' @param numBoot The number of bootstrapped coefficient values to estimate.
#' Defaults to a  number large enough that at least 10 measurement should lie
#' outside the estimated interval; but held to a minimum value of 100 and a
#' maximum value of 1000.
#'
#' @return An object of class `braidrm` but with three additional elements:
#'
#' * `ciLevs`: The two quantiles at which the confidence intervals are set
#' * `ciCoefs`: An array of bootstrapped coefficients.  The number of rows is
#' the number of *successful* bootstrapped fits; there is one column for each
#' the nine BRAID parameters
#' * `ciMat`: An array of confidence intervals on the fitted parameters. The
#' rows correspond to the *free* parameters included in the fit, and are named
#' for those parameters; the firt column contains the lower bound of the
#' confidence intervals, the second column the upper bound.
#'
#' @export
#'
#' @examples
#' bfit <- braidrm(measure ~ concA + concB, synergisticExample, getCIs = FALSE)
#' summary(bfit)
#'
#' bfit_ci <- calcBraidBootstrap(bfit)
#' summary(bfit_ci)
calcBraidBootstrap <- function(bfit,ciLevs=c(0.025,0.975),numBoot=NULL) {
	if (!inherits(bfit,"braidrm")) {
		stop("Object 'hfit' must be of class 'braidrm'.")
	}
	if (!is.null(bfit$ciLevs)) {
		warning("This BRAID fit already has parameter confidence intervals ",
				"estimated. These will be deleted and replaced.")
		bfit$ciLevs <- NULL
		bfit$ciCoefs <- NULL
		bfit$ciMat <- NULL
	}
	if (is.null(numBoot)) { numBoot <- round(max(min(10/(1-ciLevs[2]+ciLevs[1]),1000),100)) }

	bcoefs <- array(NA,dim=c(numBoot,9))
	scenarioFunction <- getScenarioFunction(bfit$scenario)
	for (i in 1:numBoot) {
		bact <- bfit$fitted.values+sample(bfit$residuals,length(bfit$residuals),replace=TRUE)
		tfit <- try({
			do.call(scenarioFunction,
					list(bfit$concs,bact,bfit$model,bfit$weights,bfit$coefficients,bfit$direction,bfit$pbounds,bfit$kweight))
		})

		if (!inherits(tfit,"try-error")) { bcoefs[i,] <- tfit$coefficients }
	}
	bcoefs <- bcoefs[!is.na(bcoefs[,1]),]
	qmat <- t(apply(bcoefs,2,stats::quantile,probs=ciLevs))
	rownames(qmat) <- names(bfit$coefficients)
	qmat <- qmat[bfit$model,]

	structure(
		c(unclass(bfit),list(ciLevs=ciLevs,ciCoefs=bcoefs,ciMat=qmat)),
		class="braidrm"
	)
}


#' Generic BRAID confidence intervals
#'
#' Generates confidence intervals on derived BRAID response surface values
#'
#' @param bfit A BRAID fit object of class `braidrm` which contains a full set
#' of bootstrapped response surface coefficients
#' @param parfunc A function that takes a full-length BRAID parameter vector as
#' an input and gives a single numeric value or numeric vector as an output. If
#' the function produces a vector, it must produce the same length vector for
#' all inputs
#' @param civals If given, the lower and upper quantile values at which the
#' confidence intervals are set.  Defaults to the `ciLevs` parameter of the
#' bootstrapped BRAID fit
#'
#' @details
#' In come cases, it is desirable to estimate a confidence interval on a value
#' derived from or dependent on a BRAID surface model that is not a parameter of
#' the model itself. For example, one might want a confidence interval on a
#' given index of achievable efficacy value, or the predicted effect at a
#' certain set of dose pairs. This function replicates confidence interval
#' calculations on any such derive values
#'
#'
#' @return An n-by-3 array, where n is the length of the output produced by
#' `parfunc`.  The first column is the lower bound of the confidence interval;
#' the second column is the derived value for the best fit coefficients; and the
#' third column is the upper bound of the confidence intervals. Note that is
#' possible for the lower bound of the confidence interval to lie above the
#' central value, or for the upper bound to lie below it; though this is only
#' likely to occur in the case of a poorly determined fit.
#' @export
#'
#' @examples
#' bfit <- braidrm(measure ~ concA + concB, synergisticExample, getCIs=TRUE)
#'
#' calcBraidConfInt(bfit, function(p) evalBraidModel(10, 10, p))
#' calcBraidConfInt(bfit, function(p) estimateIAE(p, c(0.5, 0.9), c(10, 10)))
calcBraidConfInt <- function(bfit,parfunc,civals=NULL) {
	if (!inherits(bfit,"braidrm")) {
		stop("Object 'bfit' must be of class 'braidrm'.")
	}
	if (!inherits(parfunc,"function")) { stop("Input 'parfunc' must be a function of one variable.") }
	if (is.null(bfit$ciLevs)) { stop("Input 'bfit' must have bootstrapped coefficients.") }
	if (is.null(civals)) { civals <- bfit$ciLevs }

	outval <- parfunc(bfit$coefficients)
	outmat <- cbind(array(NA,dim=c(length(outval),nrow(bfit$ciCoefs))))
	for (b in 1:nrow(bfit$ciCoefs)) { outmat[,b] <- parfunc(bfit$ciCoefs[b,]) }

	outci <- apply(outmat,1,stats::quantile,probs=civals)
	fullout <- cbind(as.vector(outci[1,]),outval,as.vector(outci[2,]))
	colnames(fullout) <- c("Lo","Est","Hi")
	return(fullout)
}
