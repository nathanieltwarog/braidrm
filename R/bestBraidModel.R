
#' Select Best BRAID Response Fit
#'
#' Picks the most parsimonious BRAID fit from a standard set of commonly used
#' variants
#'
#' @inheritParams braidrm
#' @param defaults Default minimal and maximal effect values used to fix effect
#' parameters during model selection.
#' @param extended Should models with an additional freely varying Ef parameter
#' be included.  If `FALSE` (the default), ten models in which the maximal
#' effect parameter Ef is constrained to be equal to one or more of the two
#' individual maximal effect parameters will be tested; if `TRUE`, an additional
#' two models in which Ef varies freely will be included.
#' @param useBIC If `TRUE` (the default), the best (read: most parsimonious)
#' model will be selected from all tested models using the Bayesian information
#' criterion (Schwarz 1978).  If `FALSE` the function will follow the convention
#' of earlier versions of the `braidrm` package and use the Akaike information
#' criterion (Akaike 1974).
#'
#' @details
#' When fitting real experimental data, it is not uncommon for a measured
#' surface to contain such incomplete or noisy data that many of the parameters
#' are highly underdetermined.  Unfortunately, in such cases, non-linear
#' optimization can often resort to wildly implausible values to explain small
#' variations in the data. To address this, this function runs multiple BRAID
#' response fits, including some in which the minimal and maximal effect
#' parameters are constrained to reasonable default values, to test if
#' additional free parameters offer sufficiently improved fits to be included.
#'
#' When the parameter `extended` is set to `FALSE`, the function runs ten BRAID
#' scenarios: five in which the minimal effect parameter is allowed to vary
#' freely, and five in which it is fixed at the first default value.  The five
#' tested models in each set represent five distinct configurations of the
#' maximal effect parameters:
#'
#' * Both maximal effects are fixed the same value (the second default)
#' * Maximal effect EfA (and when it is larger, Ef) varies freely, but effect
#' EfB is fixed at the second default
#' * Maximal effect EfB (and when it is larger, Ef) varies freely, but effect
#' EfA is fixed at the second default
#' * The maximal effect Ef varies freely, and both EfA and EfB are constrained
#' to be equal to it
#' * The maximal effects EfA and EfB both vary freely, and Ef is constrained to
#' be equal to the larger of the two
#'
#' When `extended` is `TRUE`, two additional models (one with E0 fixed and one
#' in which it varies freely) are included, in which all three maximal effect
#' parameters are allowed to vary freely and independently.
#'
#'
#' @return An object of class `braidrm`.  It will contain all the fields of a
#' standard `braidrm` object, and also an additional field, `allfits`
#' containing a summary of the best fit model from each of the 10 or 12
#' candidate models tested.
#'
#' @references
#' Akaike, Hirotugu. 1974. “A New Look at the Statistical Model Identification.”
#' *IEEE Transactions on Automatic Control* **19** (6): 716–23.
#'
#' Schwarz, Gideon. 1978. “Estimating the Dimension of a Model.”
#' *The Annals of Statistics*, 461–64.
#'
#' @export
#'
#' @examples
#' bfit1 <- findBestBraid(measure ~ concA + concB, additiveExample,
#'                        defaults=c(0,1))
#' summary(bfit1)
#' length(bfit1$allfits)
#'
#' bfit2 <- findBestBraid(measure ~ concA + concB, additiveExample,
#'                        defaults=c(0,2), extended=TRUE, getCIs = FALSE)
#' summary(bfit2)
#' length(bfit2$allfits)
findBestBraid <- function(formula,data,defaults,extended=FALSE,weights=NULL,
						  start=NULL,direction=0,lower=NULL,upper=NULL,
						  prior="moderate",getCIs=TRUE,useBIC=TRUE) UseMethod("findBestBraid")

#' @export
#' @rdname findBestBraid
findBestBraid.formula <- function(formula,data,defaults,extended=FALSE,weights=NULL,
								  start=NULL,direction=0,lower=NULL,upper=NULL,
								  prior="moderate",getCIs=TRUE,useBIC=TRUE) {
	mf <- stats::model.frame(formula=formula, data=data)
	concs <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	bfit <- findBestBraid.default(concs,act,defaults,extended,weights,
								  start=start,direction=direction,lower=lower,upper=upper,
								  prior=prior,getCIs=getCIs,useBIC=useBIC)
	bfit$call <- match.call()
	return(bfit)
}

#' @export
#' @rdname findBestBraid
findBestBraid.default <- function(formula,data,defaults,extended=FALSE,weights=NULL,
							start=NULL,direction=0,lower=NULL,upper=NULL,
							prior="moderate",getCIs=TRUE,useBIC=TRUE) {
	concs <- formula
	act <- data

	# Fill weights
	if (is.null(weights)) { weights <- rep(1,length(act)) }

	# Fill starting vector
	originalStart <- start
	if (isBraidParameter(start)) {
		start <- fillOutBraidPar(start)
		start[6:9] <- c(defaults[[1]],rep(defaults[[2]],3))
	} else {
		start <- defaultStartingVector(concs,act,1:8,"F",defaults,direction)
		start[6:9] <- c(defaults[[1]],rep(defaults[[2]],3))
	}

	# Fill parameter bounds
	pbounds <- fillParameterBounds(lower,upper,1:9,concs,originalStart)
	if (isBraidParameter(originalStart) && is.null(lower)) {
		pbounds[1,1:5] <- pmin(pbounds[1,1:5],start[1:5])
	} else {
		start[1:5] <- pmax(pbounds[1,1:5],start[1:5])
	}
	if (isBraidParameter(originalStart) && is.null(upper)) {
		pbounds[2,1:5] <- pmax(pbounds[2,1:5],start[1:5])
	} else {
		start[1:5] <- pmin(pbounds[2,1:5],start[1:5])
	}
	pbounds[1,6:9] <- pmin(pbounds[1,6:9],start[6:9])
	pbounds[2,6:9] <- pmax(pbounds[2,6:9],start[6:9])

	# Fill in kappa prior
	if (is.null(prior) || is.character(prior)) {
		if (is.null(prior)) { prior <- "moderate" }
		if (prior=="none") { spread <- 1 }
		else { spread <- defaultSurfaceSpread(concs,act,weights,originalStart) }
		prior <- kappaPrior(spread,prior)
	}
	if (!inherits(prior,"kappaPrior")) {
		stop("Parameter 'prior' must be a character string specifying a prior strength",
			 "\nor an object of class 'kappaPrior' produced by the function 'kappaPrior'.")
	}
	kweight <- (prior$spread^2)*(prior$strength^2)

	if (extended) { candidateModels <- bestBraidModels_extended() }
	else { candidateModels <- bestBraidModels() }
	allfits <- list()
	allmodels <- list()
	for (miter in seq_along(candidateModels$models)) {
		model <- candidateModels$models[[miter]]
		scenario <- candidateModels$scenarios[[miter]]
		scenarioFunction <- getScenarioFunction(scenario)

		cbfit <- do.call(scenarioFunction,
						list(concs,act,model,weights,start,direction,pbounds[,model],kweight))
		allmodels[[miter]] <- cbfit
		allfits[[miter]] <- list(scenario=scenario,
								 model=model,
								 coefficients=cbfit$coefficients,
								 AIC=basicdrm::estimateAIC(weights*cbfit$residuals,length(model)),
								 BIC=basicdrm::estimateBIC(weights*cbfit$residuals,length(model)))
	}

	if (useBIC) { icv <- sapply(allfits,function(m) m$BIC) }
	else { icv <- sapply(allfits,function(m) m$AIC) }
	# parv <- sapply(allfits,function(m) length(m$model))
	# biv <- basicdrm:::modelSelect(icv,parv,ik=5)
	biv <- which.min(icv)

	bfit <- allmodels[[biv]]
	bfit$allfits <- allfits

	if (getCIs) { bfit <- calcBraidBootstrap(bfit) }
	bfit
}


bestBraidModels <- function() {
	list(
		models = list(
			c(1:5),         # E0 fixed and all maximal effects fixed
			c(1:5,6),       # E0 varying and all all maximal effects fixed
			c(1:5,9),       # E0 fixed, Ef varying, EfA and EfB equal to Ef
			c(1:5,6,9),     # E0 and Ef varying, EfA and EfB equal to Ef
			c(1:5,7),       # E0 and EfB fixed, EfA varying, Ef equal to max
			c(1:5,6,7),     # E0 and EfA varying, EfB fixed, Ef equal to max
			c(1:5,8),       # E0 and EfA fixed, EfB varying, Ef equal to max
			c(1:5,6,8),     # E0 and EfB varying, EfA fixed, Ef equal to max
			c(1:5,7,8),     # E0 fixed, EfA and EfB varying, Ef equal to max
			c(1:5,6,7,8)    # E0, EfA, and EfB varying, Ef equal to max
		),
		links = c(
			"",             # E0 fixed and all maximal effects fixed
			"",             # E0 varying and all all maximal effects fixed
			"AB",           # E0 fixed, Ef varying, EfA and EfB equal to Ef
			"AB",           # E0 and Ef varying, EfA and EfB equal to Ef
			"F",            # E0 and EfB fixed, EfA varying, Ef equal to max
			"F",            # E0 and EfA varying, EfB fixed, Ef equal to max
			"F",            # E0 and EfA fixed, EfB varying, Ef equal to max
			"F",            # E0 and EfB varying, EfA fixed, Ef equal to max
			"F",            # E0 fixed, EfA and EfB varying, Ef equal to max
			"F"             # E0, EfA, and EfB varying, Ef equal to max
		),
		scenarios = c(
			"Z_1",          # E0 fixed and all maximal effects fixed
			"I_1",          # E0 varying and all all maximal effects fixed
			"I_7",          # E0 fixed, Ef varying, EfA and EfB equal to Ef
			"II_6",         # E0 and Ef varying, EfA and EfB equal to Ef
			"I_5A",         # E0 and EfB fixed, EfA varying, Ef equal to max
			"II_4A",        # E0 and EfA varying, EfB fixed, Ef equal to max
			"I_5B",         # E0 and EfA fixed, EfB varying, Ef equal to max
			"II_4B",        # E0 and EfB varying, EfA fixed, Ef equal to max
			"II_10",        # E0 fixed, EfA and EfB varying, Ef equal to max
			"III_4"         # E0, EfA, and EfB varying, Ef equal to max
		)
	)
}

bestBraidModels_extended <- function() {
	baseModels <- bestBraidModels()
	list(
		models = c(baseModels$models,list(
			c(1:5,7,8,9),   # E0 fixed, EfA, EfB, and Ef all varying
			c(1:5,6,7,8,9)  # E0, EfA, EfB, and Ef all varying
		)),
		links = c(baseModels$links,
			"",             # E0 fixed, EfA, EfB, and Ef all varying
			""              # E0, EfA, EfB, and Ef all varying
		),
		scenarios = c(baseModels$scenarios,
			"III_5",        # E0 fixed, EfA, EfB, and Ef all varying
			"IV_1"          # E0, EfA, EfB, and Ef all varying
		)
	)
}
