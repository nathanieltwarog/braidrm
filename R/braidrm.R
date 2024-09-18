
#' BRAID Response Surface Fitting
#'
#' Fits a BRAID response surface model to data.
#'
#' @param formula Either an object of class `formula` such as would be provided
#' to a modeling function like [stats::lm()], or a width-2 numeric array vector
#' of concentration pairs (including 0 or Inf).  A formula should specify a
#' single output as a function of two inputs, eg. `activity ~ conc1 + conc2`.
#' @param data If `forumula` is a symbolic formula, a data frame containing the
#' specified values. If `formula` is a numeric array of concentrations, a
#' numeric vector of response values, the same length as the number of rows of
#' `formula`.
#' @param model,links Parameters `model` and `links` are used to specify which
#' variant of the BRAID model is fit to data.  Model may be one of the
#' following character strings: "kappa1", "kappa2", or "kappa3" (see Details),
#' or a subset of the numbers 1 through 9 specifying which of the nine BRAID
#' response surface parameters is allowed to vary when fitting.  `links` allows
#' the user to further specify constraints on the three BRAID maximal effect
#' parameters (see Details for more).  If `model` is one of the supported
#' character strings, the parameter `links` will be ignored.
#' @param weights A vector of weights (between 0 and 1) the same length as
#' the data which determines the weight with which each measurement
#' will impact the the sum of squared errors.  Weights will be multiplied by
#' errors *before* squaring.  If `NULL` (the default) all weights will be set
#' to 1. Can be a numeric vector, or the name of a column in `data` if `formula`
#' is a symbolic formula
#' @param start A BRAID parameter vector specifying the first guess where the
#' non-linear optimization should begin.  May be a length 7, 8, or 9 vector,
#' though a full length vector is always preferable.  If `NULL` (the default),
#' it will be estimated from the data.
#' @param direction Determines the possible directionality of the BRAID
#' model.  If 0 (the default) no additional constraints are placed on the
#' parameters.  If greater than 0, the fitting will require that the maximal
#' effects are all *greater* than or equal to the minimal effect.  If less
#' than 0, the fitting will require that all maximal effect is *less* than or
#' equal to the minimal effect.
#' @param lower A numeric vector of lower bounds on the fitted parameter values.
#' May be the same length as the number of fitted parameters, or a full,
#' length-9 vector. Missing or unspecified lower bounds may be included as `NA`
#' or `Inf`; if unspecified, lower bounds on the first five parameters (IDMA,
#' IDMB, na, nb, and kappa) will be automatically estimated from the data.
#' Bounds on the minimal and maximal effect parameters however (E0, EfA, EfB,
#' and Ef) will be assumed to be infinite unless specified.  A value of `NULL`,
#' the default, will be treated as all lower parameter bounds being
#' unspecified.
#' @param upper A numeric vector of upper bounds on the fitted parameter values.
#' Used in the same way as `lower`.
#' @param prior A character string specifying the desired Bayesian prior term
#' for kappa, or an object of class `kappaPrior` genererated by the function
#' [kappaPrior()].  Allowed strings are "mild", "moderate" (the default),
#' "high", or "none".  If a string is given, the kappa prior object will be
#' estimated from the data using an initial ten-parameter fit to approximate
#' measurement noise.
#' @param getCIs Should bootstrapped confidence intervals be estimated and
#' added to the BRAID fit object. Default value is `TRUE`.
#' @param object An object of class `braidrm` to be summarized
#' @param x An object of class `braidrm` or `summary.braidrm` to be printed
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class `braidrm` containing the following elements:
#'
#' * `concs`: A width-two array of the concentration pairs fit by the model
#' * `act`: A vector of responses fit by the model
#' * `weights`: The vector of weights (the same length as `act`) specifying
#' the relative weight of each measurement
#' * `coefficients`: A full length-9 named BRAID parameter vector representing
#' the best fit BRAID surface for the data
#' * `fitted.values`: A vector of responses (the same length of `act`) given by
#' the best fit response surface as a function of the concentrations in `concs`
#' * `residuals`: The fitting errors of the best fit model, equal to the
#' `fitted.values` subtracted from `act`
#' `scenario`: A character string specifying one of 32 distinct fitting
#' scenarios determined from the parameters `model`, `links`, and `start`. Used
#' in bootstrapping confidence intervals
#' * `model`: The model vector (a subset of values between 1 and 9) specifying
#' which BRAID parameters were varying freely in the fit
#' * `start`: The length-9 starting BRAID parameter vector used in non-linear
#' optimization
#' * `direction`: Like the input parameter, a value of -1, 0, or 1 specifying
#' the constraint on the directionality of the fitted surface
#' * `pbounds`: A 2-by-k array of values specifying the lower and upper bounds
#' of all varying parameters (where k is the number of free parameters).
#' * `kweight`: A numeric value summarizing the relative Bayesian influence of
#' the BRAID parameter kappa on optimized objective function.
#'
#' Fit objects with bootstrapped confidence intervals include several additional
#' elements derived from that; see [calcBraidBootstrap()] for details.
#'
#' @details
#' One of the hairiest and most confusing aspects of fitting a combined
#' response surface model is handling the relationships between maximal effects.
#' Unlike a simple dose response model such as those fit in `basicdrm` in which
#' all parameters can be treated as fairly independent, response surface models
#' are often considered with constraints that cannot be expressed as simply one
#' parameter being fixed a particular value.  Many response surface models
#' assume that the two drugs in combination (and the overall combination) share
#' a single common maximal effect; others might assume that the maximal effects
#' of the two drugs should differ but that the overall maximal effect must be
#' equal to one of these; still others may wish to fit with no constraints on
#' maximal effect at all beyond guaranteeing that they lie above a fixed
#' minimal effect.  All these approaches are valid, and creating a functional
#' interface to support them all is a challenge. The parameters given here
#' represent our best effort to balance ease-of-use with flexibility.
#'
#' The primary interface for model selection and customization is the paired
#' paramters `model` and `links`.  For the first six parameters of the BRAID
#' surface, `model` is the only relevant control, and operates much as it would
#' in any fitting function.  If a given parameter (say IDMB) is included in
#' `model` (as index 2) then it will vary freely within the provided bounds to
#' best fit the data.  If it is not, the value will be fixed at the value given
#' in `start` (or if `start` is `NULL`, estimated from the data), and will
#' remain fixed at that value in the best fit surface.
#'
#' Parameters EfA, EfB, and Ef (the maximal effect parameters) require slightly
#' more care. Relationships between these values is represented by the `links`
#' parameter, which can take on one of the following five values:
#'
#' * "AB": Indicating that the overall Ef parameter is the driver, both EfA and
#' EfB are constrained to be equal to Ef, whatever its value.  Can be used
#' when indices 7 and 8 (EfA and EfB) are absent from `model`, but index 9 (Ef)
#' is present
#' * "F": Indicating instead that the individual parameters EfA and EfB are
#' the drivers, and EfA is constrained to be equal to larger magnitude of the
#' two. Can be used when indices 7 and 8 are present `model`, but index 9 is
#' absent
#' * "A": Specifies that the overall maximal effect is equal to that of the
#' first trug (and consequentially that the effect of drug A must be of greater
#' or equal magnitude to that of drug B).  Can be used when index 7 is
#' absent in `model` and index 9 is present; index 8 may be present or absent
#' * "B": Specifies that the overall maximal effect is equal to that of the
#' second drug.  Can be used when index 8 is  absent in `model` and index 9 is
#' present; index 7 may be present or absent
#' * "" (the empty string): Indicates no equality between maximal effects.
#' Parameters that are present in `model` vary freely in fitting, those that
#' are absent are fixed at constant values.
#'
#' For example, if the maximal effects *should* be fit, but should be
#' constrained to be all equal, then it is the parameter Ef that varies freely
#' in the fitting; the fact that EfA and EfB are always equal to this value is
#' represented by setting the `links` parameter to "AB".  Contrast this with the
#' (admittedly much less common) scenario in which the `links` parameter is set
#' to the empty string "", representing no link between maximal effects. In this
#' case the parameter Ef will indeed vary freely in the fitting, but EfA and
#' EfB will instead always be held at the constant initial values in the
#' starting parameter vector.  On the other hand, if we were to assume that the
#' two individual maximum effect parameters can vary independently, but that
#' the overall maximal effect should be equal to the larger of the two, then
#' indices 7 and 8 (representing EfA and EfB) would be included in `model`;
#' index 9 (representing Ef) would be excluded, and the `links` parameter would
#' be set to "F" indicating that Ef is tied to the larger of the two.
#'
#' Note also that the default value for `links` is *not* the empty string, but
#' instead `NULL`.  By default the value of `links` will be guessed from the
#' model vector, based on the scenarios that we have encountered most often.
#' If `model` inclues Ef (index 9) but not EfA or EfB (indices 7 and 8), `links`
#' is assumed to be "F"; if EfA and EfB are present but not Ef, `links` is set
#' to "AB". In the vast majority of cases, you will not need to specify this
#' parameter yourself.  This is especially true when `model` is specified with
#' one of the model strings, in which any provided value for `links` will be
#' discarded and replaced witht the following preset models:
#'
#' * kappa1: Model vector includes (1, 2, 3, 4, 5, 6, 9) and `links` is set to
#' "AB"
#' * kappa2: Model vector includes (1, 2, 3, 4, 5, 6, 7, 8) and `links` is set
#' to "F"
#' * kappa3: Model vector includes (1, 2, 3, 4, 5, 6, 7, 8, 9) and `links` is
#' set to the empty string
#'
#'
#' @export
#'
#' @examples
#' bfit1 <- braidrm(measure ~ concA + concB, additiveExample)
#' summary(bfit1)
#'
#' bfit2 <- braidrm(measure ~ concA + concB, synergisticExample,
#'                  model = c(1,2,3,4,5,6,9),
#'                  lower = c(NA,NA,NA,NA,NA,0,0),
#'                  prior = "none",
#'                  getCIs = FALSE)
#' summary(bfit2)
braidrm <- function(formula,data,model="kappa2",links=NULL,weights=NULL,
					start=NULL,direction=0,lower=NULL,upper=NULL,
					prior="moderate",getCIs=TRUE) UseMethod("braidrm")

#' @export
#' @rdname braidrm
summary.braidrm <- function(object, ...) {
	if (is.null(object$ciCoefs)) {
		tab <- cbind(object$coefficients)
		colnames(tab) <- "Est"
	} else {
		tab <- cbind(NA_real_,object$coefficients,NA_real_)
		tab[object$model,1] <- object$ciMat[,1]
		tab[object$model,3] <- object$ciMat[,2]
		colnames(tab) <- c("Lo","Est","Hi")
	}
	if (!is.null(object$call)) {
		res <- list(call=object$call,coefficients=tab)
	} else {
		res <- list(coefficients=tab)
	}
	class(res) <- "summary.braidrm"
	return(res)
}

#' @export
#' @rdname braidrm
print.summary.braidrm <- function(x, ...) {
	if (!is.null(x$call)) {
		cat("Call:\n")
		print(x$call)
	}
	cat("\n")
	print(round(x$coefficients,digits=4))
}

#' @export
#' @rdname braidrm
print.braidrm <- function(x, ...) {
	if (!is.null(x$call)) {
		cat("Call:\n")
		print(x$call)
	}
	cat("\nCoefficients:\n")
	print(x$coefficients)
}

#' @export
#' @rdname braidrm
braidrm.formula <- function(formula,data,model="kappa2",links=NULL,weights=NULL,
							start=NULL,direction=0,lower=NULL,upper=NULL,
							prior="moderate",getCIs=TRUE) {
	mf <- stats::model.frame(formula=formula, data=data)
	concs <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	bfit <- braidrm.default(concs,act,model,links=links,weights=weights,
							start=start,direction=direction,lower=lower,upper=upper,
							prior=prior,getCIs=getCIs)
	bfit$call <- match.call()
	return(bfit)
}

#' @export
#' @rdname braidrm
braidrm.default <- function(formula,data,model="kappa2",links=NULL,weights=NULL,
							start=NULL,direction=0,lower=NULL,upper=NULL,
							prior="moderate",getCIs=TRUE) {
	concs <- formula
	act <- data

	# Rectify model and links
	originalModel <- model
	originalLinks <- links
	if (is.character(model)) {
		model <- switch(originalModel,
						kappa1=c(1,2,3,4,5,6,9),
						kappa2=c(1,2,3,4,5,6,7,8),
						kappa3=c(1,2,3,4,5,6,7,8,9),
						stop(sprintf("Unknown model name '%s'.",model)))
		links <- switch(originalModel,
						kappa1="AB",
						kappa2="F",
						kappa3="",
						stop(sprintf("Unknown model name '%s'.",model)))
	} else if (is.null(links)) {
		if (all((7:9)%in%model)||!any((7:9)%in%model)) { links <- "" }
		else if (!(9%in%model)) { links <- "F" }
		else if (8%in%model) { links <- "A" }
		else if (7%in%model) { links <- "B" }
		else { links <- "AB" }
	}
	if (!is.null(originalLinks) && originalLinks!=links) {
		warning("Warning: The named model implies a value for 'links' that does not match the specified value.",
				"\nIn such cases, the value implied by the value of 'model' will be used.")
	}

	# Fill weights
	if (is.null(weights)) { weights <- rep(1,length(act)) }

	# Fill starting vector
	originalStart <- start
	if (isBraidParameter(start)) {
		start <- fillOutBraidPar(start)
		if (!checkLinks(start,links)) {
			warning("Warning: given starting parameter vector does not satisfy the equality constraints implied",
					"\nby parameters 'model' and 'links'. Starting parameter vector will be coerced to satisfy",
					"\nthese constraints")
			start <- coerceLinks(start,links)
		}
	} else {
		start <- defaultStartingVector(concs,act,model,links,start,direction)
	}

	# Fill parameter bounds
	pbounds <- fillParameterBounds(lower,upper,model,concs,originalStart)

	if (isBraidParameter(originalStart) && is.null(lower)) {
		pbounds[1,] <- pmin(pbounds[1,],start[model])
	} else {
		start[model] <- pmax(pbounds[1,],start[model])
	}
	if (isBraidParameter(originalStart) && is.null(upper)) {
		pbounds[2,] <- pmax(pbounds[2,],start[model])
	} else {
		start[model] <- pmin(pbounds[2,],start[model])
	}

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

	# Determine scenario
	scenarioFunction <- getScenario(model,links)

	bfit <- do.call(scenarioFunction,
			list(concs,act,model,weights,start,direction,pbounds,kweight))

	if (getCIs) { bfit <- calcBraidBootstrap(bfit) }
	bfit
}

defaultStartingVector <- function(concs,act,model,links="",start=NULL,direction=0) {
	if (isBraidParameter(start)) {
		return(fillOutBraidPar(start))
	} else if (is.null(start)) {
		e0 <- NULL
		ef <- NULL
	} else if (is.numeric(start) && length(start)==2) {
		e0 <- start[[1]]
		ef <- start[[2]]
	}

	conc1 <- concs[,1]
	conc2 <- concs[,2]
	prel <- conc1>0 & is.finite(conc1) & conc2>0 & is.finite(conc2)
	crng1 <- range(conc1[prel])
	crng2 <- range(conc2[prel])
	dcenter <- c(exp(mean(log(crng1))),exp(mean(log(crng2))),1,1,0)
	if (is.null(e0)) {
		erange <- stats::quantile(act[prel],c(0.05,0.95),names=FALSE)
		if (direction==0) {
			logconc <- log(conc1[prel])+log(conc2[prel])
			lact <- act[prel]
			dlm <- stats::lm(lact~logconc)
			if (stats::coef(dlm)[[2]]>=0) {
				e0 <- erange[[1]]
				ef <- erange[[2]]
			} else {
				e0 <- erange[[2]]
				ef <- erange[[1]]
			}
		} else if (direction>0) {
			e0 <- erange[[1]]
			ef <- erange[[2]]
		} else {
			e0 <- erange[[2]]
			ef <- erange[[1]]
		}
	}

	c(dcenter,e0,ef,ef,ef)
}

checkScenario <- function(scenario,model,links=NULL) {
	scheck = switch(scenario,
		   Z_1  = checkModelValues(model,c(FALSE, FALSE, FALSE, FALSE)),
		   I_1  = checkModelValues(model,c(TRUE,  FALSE, FALSE, FALSE)),
		   I_2  = checkModelValues(model,c(FALSE, FALSE, FALSE, TRUE )),
		   I_3A = checkModelValues(model,c(FALSE, TRUE,  FALSE, FALSE)),
		   I_3B = checkModelValues(model,c(FALSE, FALSE, TRUE,  FALSE)),
		   I_4A = checkModelValues(model,c(FALSE, FALSE, FALSE, TRUE )),
		   I_4B = checkModelValues(model,c(FALSE, FALSE, FALSE, TRUE )),
		   I_5A = checkModelValues(model,c(FALSE, TRUE,  FALSE, FALSE)),
		   I_5B = checkModelValues(model,c(FALSE, FALSE, TRUE,  FALSE)),
		   I_7  = checkModelValues(model,c(FALSE, FALSE, FALSE, TRUE )),
		   II_1  = checkModelValues(model,c(TRUE,  FALSE, FALSE, TRUE )),
		   II_2A = checkModelValues(model,c(TRUE,  TRUE,  FALSE, FALSE)),
		   II_2B = checkModelValues(model,c(TRUE,  FALSE, TRUE,  FALSE)),
		   II_3A = checkModelValues(model,c(TRUE,  FALSE, FALSE, TRUE )),
		   II_3B = checkModelValues(model,c(TRUE,  FALSE, FALSE, TRUE )),
		   II_4A = checkModelValues(model,c(TRUE,  TRUE,  FALSE, FALSE)),
		   II_4B = checkModelValues(model,c(TRUE,  FALSE, TRUE,  FALSE)),
		   II_6  = checkModelValues(model,c(TRUE,  FALSE, FALSE, TRUE )),
		   II_7A = checkModelValues(model,c(FALSE, TRUE,  FALSE, TRUE )),
		   II_7B = checkModelValues(model,c(FALSE, FALSE, TRUE,  TRUE )),
		   II_8  = checkModelValues(model,c(FALSE, TRUE,  TRUE,  FALSE)),
		   II_9A = checkModelValues(model,c(FALSE, TRUE,  FALSE, TRUE )),
		   II_9B = checkModelValues(model,c(FALSE, FALSE, TRUE,  TRUE )),
		   II_10 = checkModelValues(model,c(FALSE, TRUE,  TRUE,  FALSE)),
		   III_1A = checkModelValues(model,c(TRUE,  TRUE,  FALSE, TRUE )),
		   III_1B = checkModelValues(model,c(TRUE,  FALSE, TRUE,  TRUE )),
		   III_2  = checkModelValues(model,c(TRUE,  TRUE,  TRUE,  FALSE)),
		   III_3A = checkModelValues(model,c(TRUE,  TRUE,  FALSE, TRUE )),
		   III_3B = checkModelValues(model,c(TRUE,  FALSE, TRUE,  TRUE )),
		   III_4  = checkModelValues(model,c(TRUE,  TRUE,  TRUE,  FALSE)),
		   III_5  = checkModelValues(model,c(FALSE, TRUE,  TRUE,  TRUE )),
		   IV_1 = checkModelValues(model,c(TRUE,  TRUE,  TRUE,  TRUE )),
		   stop(sprintf("Unrecognized fitting scenario '%s'.",scenario))
	)
	if (!is.null(links)) { scheck <- scheck && links==getLinks(scenario) }
	scheck
}
getLinks <- function(scenario) {
	switch(scenario,
		   Z_1  = "",
		   I_1  = "",
		   I_2  = "",
		   I_3A = "",
		   I_3B = "",
		   I_4A = "A",
		   I_4B = "B",
		   I_5A = "F",
		   I_5B = "F",
		   I_7  = "AB",
		   II_1  = "",
		   II_2A = "",
		   II_2B = "",
		   II_3A = "A",
		   II_3B = "B",
		   II_4A = "F",
		   II_4B = "F",
		   II_6  = "AB",
		   II_7A = "",
		   II_7B = "",
		   II_8  = "",
		   II_9A = "B",
		   II_9B = "A",
		   II_10 = "F",
		   III_1A = "",
		   III_1B = "",
		   III_2  = "",
		   III_3A = "B",
		   III_3B = "A",
		   III_4  = "F",
		   III_5  = "",
		   IV_1 = "",
		   stop(sprintf("Unrecognized fitting scenario '%s'.",scenario))
	)
}
getScenario <- function(model,links) {
	longlinks <- paste0("link",links)
	if (6 %in% model) {
		# E0 is varying
		if (!any((7:9) %in% model)) {
			# None of the maximal effects are varying freely
			scenario <- "I_1"
		} else if (all((7:9) %in% model)) {
			# All of the maximal effects are varying freely
			scenario <- "IV_1"
		} else if (!(9 %in% model)) {
			# Ef is not varying freely, at least one of EfA or EfB is
			if (all((7:8) %in% model)) {
				# EfA and EfB are varying freely
				scenario <- switch(longlinks, linkF="III_4", link="III_2",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			} else if (7 %in% model) {
				# EfA varies freely
				scenario <- switch(longlinks, linkF="II_4A", link="II_2A",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			} else {
				# EfB varies freely
				scenario <- switch(longlinks, linkF="II_4B", link="II_2B",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			}
		} else if (!any((7:8) %in% model)) {
			# Ef is varying, EfA and EfB are not free
			scenario <- switch(longlinks, linkAB="II_6", linkA="II_3A", linkB="II_3B", link="II_1",
							   linkF=stop("This 'links' value is not compatible with this model vector"),
							   stop(sprintf("Unsupported 'links' value '%s'.",links)))
		} else {
			# Ef and one of EfA or EfB is varying, the other is not
			if (7 %in% model) {
				# EfA varies freely
				scenario <- switch(longlinks, linkB="III_3A", link="III_1A",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   linkF=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			} else {
				# EfB varies freely
				scenario <- switch(longlinks, linkA="III_3B", link="III_1B",
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   linkF=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			}
		}
	} else {
		# E0 is fixed
		if (!any((7:9) %in% model)) {
			# None of the maximal effects are varying freely
			scenario <- "Z_1"
		} else if (all((7:9) %in% model)) {
			# All of the maximal effects are varying freely
			scenario <- "III_5"
		} else if (!(9 %in% model)) {
			# Ef is not varying freely, at least one of EfA or EfB is
			if (all((7:8) %in% model)) {
				# EfA and EfB are varying freely
				scenario <- switch(longlinks, linkF="II_10", link="II_8",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			} else if (7 %in% model) {
				# EfA varies freely
				scenario <- switch(longlinks, linkF="I_5A", link="I_3A",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			} else {
				# EfB varies freely
				scenario <- switch(longlinks, linkF="I_5B", link="I_3B",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			}
		} else if (!any((7:8) %in% model)) {
			# Ef is varying, EfA and EfB are not free
			scenario <- switch(longlinks, linkAB="I_7", linkA="I_4A", linkB="I_4B", link="I_2",
							   linkF=stop("This 'links' value is not compatible with this model vector"),
							   stop(sprintf("Unsupported 'links' value '%s'.",links)))
		} else {
			# Ef and one of EfA or EfB is varying, the other is not
			if (7 %in% model) {
				# EfA varies freely
				scenario <- switch(longlinks, linkB="II_9A", link="II_7A",
								   linkA=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   linkF=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			} else {
				# EfB varies freely
				scenario <- switch(longlinks, linkA="II_9B", link="II_7B",
								   linkB=stop("This 'links' value is not compatible with this model vector"),
								   linkAB=stop("This 'links' value is not compatible with this model vector"),
								   linkF=stop("This 'links' value is not compatible with this model vector"),
								   stop(sprintf("Unsupported 'links' value '%s'.",links)))
			}
		}
	}

	getScenarioFunction(scenario)
}
getScenarioFunction <- function(scenario) {
	switch(scenario,
		   Z_1  = fitBraidScenario_Z_1,
		   I_1  = fitBraidScenario_I_1,
		   I_2  = fitBraidScenario_I_2,
		   I_3A = fitBraidScenario_I_3A,
		   I_3B = fitBraidScenario_I_3B,
		   I_4A = fitBraidScenario_I_4A,
		   I_4B = fitBraidScenario_I_4B,
		   I_5A = fitBraidScenario_I_5A,
		   I_5B = fitBraidScenario_I_5B,
		   I_7  = fitBraidScenario_I_7,
		   II_1  = fitBraidScenario_II_1,
		   II_2A = fitBraidScenario_II_2A,
		   II_2B = fitBraidScenario_II_2B,
		   II_3A = fitBraidScenario_II_3A,
		   II_3B = fitBraidScenario_II_3B,
		   II_4A = fitBraidScenario_II_4A,
		   II_4B = fitBraidScenario_II_4B,
		   II_6  = fitBraidScenario_II_6,
		   II_7A = fitBraidScenario_II_7A,
		   II_7B = fitBraidScenario_II_7B,
		   II_8  = fitBraidScenario_II_8,
		   II_9A = fitBraidScenario_II_9A,
		   II_9B = fitBraidScenario_II_9B,
		   II_10 = fitBraidScenario_II_10,
		   III_1A = fitBraidScenario_III_1A,
		   III_1B = fitBraidScenario_III_1B,
		   III_2  = fitBraidScenario_III_2,
		   III_3A = fitBraidScenario_III_3A,
		   III_3B = fitBraidScenario_III_3B,
		   III_4  = fitBraidScenario_III_4,
		   III_5  = fitBraidScenario_III_5,
		   IV_1 = fitBraidScenario_IV_1,
		   stop(sprintf("Unrecognized fitting scenario '%s'.",scenario))
	)
}

fillParameterBounds <- function(llims,ulims,model,concs,start=NULL) {
	if (isBraidParameter(start)) { start <- fillOutBraidPar(start) }
	else { start <- NULL }

	conc1 <- concs[,1]
	crng1 <- range(conc1[conc1>0 & is.finite(conc1)])
	conc2 <- concs[,2]
	crng2 <- range(conc2[conc2>0 & is.finite(conc2)])
	dllims <- c(exp((1.5*log(crng1[1])-0.5*log(crng1[2]))),exp((1.5*log(crng2[1])-0.5*log(crng2[2]))),0.1,0.1,-1.96)
	dulims <- c(exp((1.5*log(crng1[2])-0.5*log(crng1[1]))),exp((1.5*log(crng2[2])-0.5*log(crng2[1]))),10,10,100)

	if (is.null(llims)) { llims <- rep(NA,length(model)) }
	else if (length(model)<9 & length(llims)==9) { llims <- llims[model] }
	if (is.null(ulims)) { ulims <- rep(NA,length(model)) }
	else if (length(model)<9 & length(ulims)==9) { ulims <- ulims[model] }
	if (any(c(length(llims),length(ulims))!=length(model))) {
		stop("Upper and lower bounding vectors must be full length (length 9) or the same length as the number of free parameters.")
	}

	for (i in 1:4) {
		if (!(i%in%model)) { next }
		ci <- which(model==i)
		if (is.null(start)) {
			if (!is.finite(log(llims[[ci]]))) { llims[[ci]] <- dllims[[i]] }
			if (!is.finite(log(ulims[[ci]]))) { ulims[[ci]] <- dulims[[i]] }
		} else {
			if (!is.finite(log(llims[[ci]]))) { llims[[ci]] <- min(start[[i]],dllims[[i]]) }
			if (!is.finite(log(ulims[[ci]]))) { ulims[[ci]] <- max(start[[i]],dulims[[i]]) }
		}
	}
	if (5 %in% model) {
		ci <- which(model==5)
		if (is.null(start)) {
			if (!is.finite(log(llims[[ci]]+2))) { llims[[ci]] <- dllims[[5]] }
			if (!is.finite(log(ulims[[ci]]+2))) { ulims[[ci]] <- dulims[[5]] }
		} else {
			if (!is.finite(log(llims[[ci]]+2))) { llims[[ci]] <- min(start[[5]],dllims[[5]]) }
			if (!is.finite(log(ulims[[ci]]+2))) { ulims[[ci]] <- max(start[[5]],dulims[[5]]) }
		}
	}
	for (i in 6:9) {
		if (!i%in%model) { next }
		ci <- which(model==i)
		if (is.na(llims[[ci]])) { llims[[ci]] <- -Inf }
		if (is.na(ulims[[ci]])) { ulims[[ci]] <- Inf }
	}
	return(rbind(llims,ulims))
}
coerceLinks <- function(par,links) {
	if (links=="") { return(par) }
	if (links=="F") {
		incr <- sign((par[[7]]+par[[8]])/2-par[[6]])
		par[[9]] <- incr*max(incr*par[7:8])
	} else if (links=="A") { par[[7]] <- par[[9]] }
	else if (links=="B") { par[[8]] <- par[[9]] }
	else { par[7:8] <- par[[9]] }
	return(par)
}
checkLinks <- function(par,links) {
	if (links=="") { return(TRUE) }
	if (links=="F") {
		incr <- sign((par[[7]]+par[[8]])/2-par[[6]])
		return(par[[9]]==incr*max(incr*par[7:8]))
	}
	return(switch(links,
				  A=(par[[7]]==par[[9]]),
				  B=(par[[8]]==par[[9]]),
				  AB=(par[[7]]==par[[9]]&&par[[8]]==par[[9]]),
				  TRUE))
}

isBraidParameter <- function(par) {
	is.numeric(par) && (length(par)>=7) && (length(par)<=9)
}
