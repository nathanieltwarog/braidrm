
#' Evaluate Flipped BRAID Surfaces
#'
#' @inheritParams evalBraidModel
#' @param bpar Flipped-BRAID parameter of the flipped response surface. See
#' details for more information on specifying atypical surfaces
#' @param flip String specifying the direction or directions of the surface's
#' flip.  Must be one of "A", "B", or "both".
#'
#' @details
#' While the BRAID model generates a fairly versatile range of combined
#' behaviors, the traditional model is still strictly constrained in certain
#' respects.  Response surfaces must exhibit a change in response to both drugs
#' and this change must be in the same direction.  Furthermore, the model is
#' only suited to surfaces in which the change resulting from single drugs is
#' larger than the additional effect of the combination.  However, modifying
#' the model equation by inverting one or both of the slope parameters can
#' produce a set of qualitatively distinct "flipped" surfaces, allowing the
#' model to produce a much wider range of behaviors.  For example, flipping the
#' model along the axis representing drug A can produce a surface in which drug
#' B has no effect in isolation, but attenuates or eliminates the effect of
#' drug A, a pattern we call a "protective" surface.  See [fitBraidFlipped()]
#' for the possible range of surfaces.
#'
#' An important note: a flipped BRAID surface, like a traditional BRAID surface,
#' is represented by a parameter vector of up to 9 values.  While these
#' functions will attempt to fill a 7- or 8- value parameter factor to a full
#' 9-element vector, it is strongly recommended that you specify the response
#' surface with the full 9-element vector, as the precise ordering with which
#' implicit values are arranged can be extremely confusing.  Note also that
#' flipped parameter vectors, regardless of the underlying mathematical
#' representation, should always be specified in the same order as in a
#' traditional vector.  So in a full 9-element vector:
#'
#' * Parameter 6 (E0) should always specify the expected effect when both drugs
#' are absent
#' * Parameter 7 (EfA) should always specify the expected effect at high
#' concentrations of drug A when drug B is absent
#' * Parameter 8 (EfB) should always specify the expected effect at high
#' concentrations of drug B when drug A is absent
#' * Parameter 9 (Ef) should always specify the expected effect when both drugs
#' are present at high concentrations
#'
#'
#' @return A numeric vector the same length as `DA`  and/or `DB` with the
#' predicted flipped BRAID response surface values.
#' @export
#'
#' @examples
#' concentrations <- c(0, 2^(-3:3))
#' surface <- data.frame(
#'     concA = rep(concentrations,each=length(concentrations)),
#'     concB = rep(concentrations,times=length(concentrations))
#' )
#'
#' surface$protective <- evalFlippedBraidModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, 0, 0, 100, 0, 10),
#'     flip="A"
#' )
#' surface$coactive <- evalFlippedBraidModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, 0, 0, 0, 0, 100),
#'     flip="both"
#' )
#'
#' head(surface)
evalFlippedBraidModel <- function(DA,DB,bpar,flip) {
	anchors <- bpar[1:2]^2
	flipbpar <- flippedFillOutBraidPar(bpar,flip,anchors)
	if (flip %in% c("A","both")) {
		flipDA <- anchors[[1]]/DA
	} else { flipDA <- DA }
	if (flip %in% c("B","both")) {
		flipDB <- anchors[[2]]/DB
	} else { flipDB <- DB }

	evalBraidModel(flipDA,flipDB,flipbpar)
}

#' Invert Flipped BRAID Surfaces
#'
#' Given a particular effect and one of the doses in a flipped combined action
#' response surface, this function calculates the other dose that will produce
#' the desired effect.`invertFlippedBraidModelA` and
#' `invertFlippedBraidModelB` are convenience wrapper functions that set `DA`
#' or `DB` to `NULL` to estimate the necessary concentrations of drug A and
#' drug B respectively.
#'
#' @inheritParams invertBraidModel
#' @param bpar Flipped-BRAID parameter of the flipped response surface. See
#' [evalFlippedBraidModel()] for more information on specifying atypical
#' surfaces
#' @param flip String specifying the direction or directions of the surface's
#' flip.  Must be one of "A", "B", or "both".
#'
#' @return A vector of concentrations the same length as either `DA` or `DB`
#' (whichever is not `NULL`) and/or `effect`, representing the concentration of
#' the other drug producing the specified effect in combination with the given
#' dose of the provided drug
#' @export
#'
#' @examples
#' fbfit <- fitProtectiveBraid_A(measure ~ concA + concB,
#'                               protectiveExample, getCIs=FALSE)
#'
#' invertFlippedBraidModel_A(DB=0, effect=0.5, coef(fbfit), fbfit$flip)
#' invertFlippedBraidModel_A(DB=0.75, effect=0.5, coef(fbfit), fbfit$flip)
invertFlippedBraidModel <- function(DA=NULL,DB=NULL,effect,bpar,flip,invalidNA=FALSE,lowerBound=FALSE) {
	anchors <- bpar[1:2]^2
	flipbpar <- flippedFillOutBraidPar(bpar,flip,anchors)
	if (flip %in% c("A","both") && !is.null(DA)) {
		flipDA <- anchors[[1]]/DA
	} else { flipDA <- DA }
	if (flip %in% c("B","both") && !is.null(DB)) {
		flipDB <- anchors[[2]]/DB
	} else { flipDB <- DB }

	if (flip %in% c("A","both") && !is.null(DB)) {
		output <- anchors[[1]]/invertBraidModel(flipDA,flipDB,effect,flipbpar,invalidNA,lowerBound)
	} else if (flip %in% c("B","both") && !is.null(DA)) {
		output <- anchors[[2]]/invertBraidModel(flipDA,flipDB,effect,flipbpar,invalidNA,lowerBound)
	} else {
		output <- invertBraidModel(flipDA,flipDB,effect,flipbpar,invalidNA,lowerBound)
	}

	output
}

#' @export
#' @rdname invertFlippedBraidModel
invertFlippedBraidModel_A <- function(DB,effect,bpar,flip,invalidNA=FALSE,lowerBound=FALSE) {
	invertFlippedBraidModel(DA=NULL,DB=DB,effect,bpar,flip,invalidNA=invalidNA,lowerBound=lowerBound)
}

#' @export
#' @rdname invertFlippedBraidModel
invertFlippedBraidModel_B <- function(DA,effect,bpar,flip,invalidNA=FALSE,lowerBound=FALSE) {
	invertFlippedBraidModel(DA=DA,DB=NULL,effect,bpar,flip,invalidNA=invalidNA,lowerBound=lowerBound)
}

#' Fit Flipped BRAID Surfaces
#'
#' Functions to fit protective, oppoistional and coactive BRAID surfaces, and
#' more specific flipped surfaces if necessary.
#'
#' @inheritParams braidrm
#' @inheritParams evalFlippedBraidModel
#' @param ... Additional parameters to be passed to `braidrm`
#'
#' @details
#' Though `fitBraidFlipped` offers the option of fitting any flipped BRAID
#' surface model specified by `flip`, this is not recommended, as the interplay
#' between flipping paramteers and parameter constraints becomes very confusing
#' very quickly.  In nearly all cases, it is preferable to use one of the
#' pre-defined flipped fitting functions.
#'
#' `fitProtectiveBraid_A` and `fitProtectiveBraid_B` fit "protective" surfaces
#' in which one drug has no effect in isolation, but attenuates or eliminates
#' the effect of the other.  `fitProtectiveBraid_A` generates a surface in which
#' drug A is active and is attenuated by drug B; `fitProtectiveBraid_B`
#' generates the reverse.
#'
#' `fitOppositinalBraid_A` and `fitOppositionalBraid_B` produce "oppositional"
#' surfaces in which a second drug produces an effect that is in the opposite
#' direction to the first drug, but which is then overwhelmed by the effect of
#' the first drug at higher concentrations.  `fitProtectiveBraid_A` generates a
#' surface in which the maximal effect of drug A dominates a high concentrations,
#' `fitProtectiveBraid_B` generates the reverse.  Note that the `A` and `B` in
#' the function names specify which compound's effect is dominant, not the
#' direction of the underlying flip; in actuality the surfaces generated by
#' `fitProtectiveBraid_A` are produced by flipping along the `B` axis.
#'
#' `fitCoactiveBraid_pure` and `fitCoactiveBraid_partial` produce "coactive"
#' surfaces, in which both drugs have no or minimal effect in isolation, but
#' produce a pronounced effect when both are present.  `fitCoactiveBraid_pure`
#' generates surfaces in which both drugs have no effect at all in isolation;
#' `fitCoactiveBraid_partial` generates surfaces in which either drug may have
#' a smaller partial effect in isolation.
#'
#' @return A fit object of class `braidrmflip`.  This structure contains the
#' exact same elements as an object of class `braidrm` (see [braidrm()] for
#' details) along with one additional element: `flip`, a character value
#' specifying the direction that the surface is flipped.  The object's
#' `coefficients` and `flip` fields can be used to evaluate and invert the best
#' fit response surface using [evalFlippedBraidModel()] and
#' [invertFlippedBraidModel()].
#' @export
#'
#' @examples
#' fbfit1 <- fitProtectiveBraid_A(measure ~ concA + concB,
#'                                protectiveExample, getCIs=FALSE)
#' coef(fbfit1)
#'
#' fbfit2 <- fitOppositionalBraid_A(measure ~ concA + concB,
#'                                  oppositionalExample, getCIs=FALSE)
#' coef(fbfit2)
#'
#' fbfit3 <- fitCoactiveBraid_pure(measure ~ concA + concB,
#'                                 coactiveExample, getCIs=FALSE)
#' coef(fbfit3)
fitBraidFlipped <- function(formula,data,flip,model,links=NULL,...) {
	fitBRAIDflipped_internal(formula,data,flip,dirflip=TRUE, model,links,...)
}


#' @export
#' @rdname fitBraidFlipped
fitProtectiveBraid_A <- function(formula,data,...) {
	fitBRAIDflipped_internal(formula,data,flip="A",dirflip=TRUE,model=c(1:5,7,8,9),links="A",...)
}
#' @export
#' @rdname fitBraidFlipped
fitProtectiveBraid_B <- function(formula,data,...) {
	fitBRAIDflipped_internal(formula,data,flip="B",dirflip=TRUE,model=c(1:5,7,8,9),links="B",...)
}


#' @export
#' @rdname fitBraidFlipped
fitOppositionalBraid_A <- function(formula,data,...) {
	fitBRAIDflipped_internal(formula,data,flip="B",dirflip=FALSE,model=c(1:5,6,7,8),links="A",...)
}
#' @export
#' @rdname fitBraidFlipped
fitOppositionalBraid_B <- function(formula,data,...) {
	fitBRAIDflipped_internal(formula,data,flip="A",dirflip=FALSE,model=c(1:5,6,7,8),links="B",...)
}


#' @export
#' @rdname fitBraidFlipped
fitCoactiveBraid_pure <- function(formula,data,...) {
	fitBRAIDflipped_internal(formula,data,flip="both",dirflip=TRUE,model=c(1:5,6,9),links="AB",...)
}
#' @export
#' @rdname fitBraidFlipped
fitCoactiveBraid_partial <- function(formula,data,...) {
	fitBRAIDflipped_internal(formula,data,flip="both",dirflip=TRUE,model=c(1:5,6,7,8,9),links="",...)
}


fitBRAIDflipped_internal <- function(formula,data,flip,dirflip=TRUE,
									 model,links=NULL,weights=NULL,...) {
	if (inherits(formula,"formula")) {
		fitBRAIDflipped_internal.formula(formula,data,flip,dirflip,model,links,weights,...)
	} else {
		fitBRAIDflipped_internal.default(formula,data,flip,dirflip,model,links,weights,...)
	}
}

fitBRAIDflipped_internal.formula <- function(formula,data,flip,dirflip=TRUE,
											 model,links=NULL,weights=NULL,...) {
	mf <- stats::model.frame(formula=formula, data=data)
	concs <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	bfit <- fitBRAIDflipped_internal.default(concs,act,flip,dirflip,model,links,weights,...)
	return(bfit)
}

fitBRAIDflipped_internal.default <- function(formula,data,flip,dirflip=TRUE,
									 model,links=NULL,weights=NULL,
									 start=NULL,direction=0,lower=NULL,upper=NULL,
									 prior="moderate",getCIs=TRUE) {
	if (dirflip) { direction <- -direction }

	concs <- formula
	act <- data

	c1v <- concs[,1]
	c2v <- concs[,2]
	c1range <- range(c1v[is.finite(log(c1v))])
	c2range <- range(c2v[is.finite(log(c2v))])
	anchors <- c(prod(c1range),prod(c2range))
	if (flip %in% c("A","both")) { c1v <- anchors[[1]]/c1v }
	if (flip %in% c("B","both")) { c2v <- anchors[[2]]/c2v }
	flipconcs <- cbind(c1v,c2v)

	originalModel <- model
	if (is.character(model)) {
		model <- switch(originalModel,
						kappa1=c(1,2,3,4,5,6,9),
						kappa2=c(1,2,3,4,5,6,7,8),
						kappa3=c(1,2,3,4,5,6,7,8,9),
						stop(sprintf("Unknown model name '%s'.",model)))
	}
	flippedModel <- flipModelVector(model,flip)

	originalStart <- start
	if (!is.null(start)) {
		if (length(start)==9) { start <- flipParameterVector(start,flip,anchors) }
		else {
			stop("For simplicity, start vectors passed to flipped fitting",
				 " functions must be full 9-length BRAID vectors.")
		}
	}

	originalLower <- lower
	originalUpper <- upper
	lower <- flipBoundVector(lower,model,flip,anchors)
	upper <- flipBoundVector(upper,model,flip,anchors)
	if ((flip %in% c("A","both")) && 1 %in% model) {
		if (!is.null(lower)) {
			if (length(lower)==9) { lowind <- 1 }
			else { lowind <- which(model==1) }
			lowval <- lower[lowind]
		} else { lowval <- NA }
		if (!is.null(upper)) {
			if (length(upper)==9) { highind <- 1 }
			else { highind <- which(model==1) }
			highval <- upper[highind]
		} else { highval <- NA }
		if (is.null(upper)) {
			if (is.finite(log(lowval))) {
				upper <- rep(NA,9)
				upper[[1]] <- lowval
			}
		} else { upper[[1]] <- lowval }
		if (is.null(lower)) {
			if (is.finite(log(highval))) {
				lower <- rep(NA,9)
				lower[[1]] <- highval
			}
		} else { lower[[1]] <- highval }
	}
	if ((flip %in% c("B","both")) && 2 %in% model) {
		if (!is.null(lower)) {
			if (length(lower)==9) { lowind <- 2 }
			else { lowind <- which(model==2) }
			lowval <- lower[lowind]
		} else { lowval <- NA }
		if (!is.null(upper)) {
			if (length(upper)==9) { highind <- 2 }
			else { highind <- which(model==2) }
			highval <- upper[highind]
		} else { highval <- NA }
		if (is.null(upper)) {
			if (is.finite(log(lowval))) {
				upper <- rep(NA,9)
				upper[[2]] <- lowval
			}
		} else { upper[[2]] <- lowval }
		if (is.null(lower)) {
			if (is.finite(log(highval))) {
				lower <- rep(NA,9)
				lower[[2]] <- highval
			}
		} else { lower[[2]] <- highval }
	}

	bfit <- braidrm.default(flipconcs,act,model=flippedModel,links=links,weights=weights,
							start=start,direction=direction,lower=lower,upper=upper,
							prior=prior,getCIs=getCIs)

	flipBraidFitObject(bfit,flip,anchors)
}

flipBraidFitObject <- function(bfit,flip,anchors) {
	concs <- bfit$concs
	if (flip %in% c("A","both")) { concs[,1] <- anchors[[1]]/concs[,1] }
	if (flip %in% c("B","both")) { concs[,2] <- anchors[[2]]/concs[,2] }
	model <- flipModelVector(bfit$model,flip)
	coefficients = flipParameterVector(unname(bfit$coefficients),flip,anchors)
	names(coefficients) <- names(bfit$coefficients)

	fbfit <- list(
		concs = concs,
		act = bfit$act,
		flip = flip,
		weights = bfit$weights,
		coefficients = coefficients,
		par = bfit$par,
		fitted.values = bfit$fitted,
		residuals = bfit$residuals,
		scenario = bfit$scenario,
		model = model,
		start = flipParameterVector(bfit$start,flip,anchors),
		direction = bfit$direction,
		pbounds = flipBoundMat(bfit$pbounds,bfit$model,flip,anchors),
		kweight = bfit$kweight
	)

	if (!is.null(bfit$ciLevs)) {
		ciCoefs <- bfit$ciCoefs[,flipIndices(flip)]
		if (flip %in% c("A","both")) { ciCoefs[,1] <- anchors[[1]]/ciCoefs[,1] }
		if (flip %in% c("B","both")) { ciCoefs[,2] <- anchors[[2]]/ciCoefs[,2] }
		ciMat <- t(flipBoundMat(t(unname(bfit$ciMat)),bfit$model,flip,anchors))
		rownames(ciMat) <- names(coefficients)[model]
		fbfit <- c(fbfit,list(ciLevs=bfit$ciLevs,ciCoefs=ciCoefs,ciMat=ciMat))
	}
	structure(fbfit,class="braidrmflip")
}

flipBoundMat <- function(boundmat,model,flip,anchors) {
	fboundmat <- array(NA,dim=c(2,9))
	fboundmat[,model] <- boundmat
	fboundmat <- fboundmat[,flipIndices(flip)]
	if (flip %in% c("A","both")) {
		fboundmat[,1] <- anchors[[1]]/fboundmat[c(2,1),1]
	}
	if (flip %in% c("B","both")) {
		fboundmat[,2] <- anchors[[2]]/fboundmat[c(2,1),2]
	}
	newboundmat <- fboundmat[,flipModelVector(model,flip)]
	newboundmat
}
flipBoundVector <- function(bound,model,flip,anchors) {
	if (is.null(bound)) { return(NULL) }
	else if (length(bound)==9) { return(flipParameterVector(bound,flip,anchors)) }
	else if (length(bound)==length(model)) {
		fbound <- rep(NA,9)
		fbound[model] <- bound
		fbound <- flipParameterVector(fbound,flip,anchors)
		newbound <- fbound[flipModelVector(model,flip)]
		return(newbound)
	} else {
		stop("Upper and lower bound paramters, if not null, must be full, 9-",
			 "length vectors or the same length as the model selected.")
	}
}
flippedFillOutBraidPar <- function(bpar,flip,anchors=NULL) {
	if (is.null(anchors)) { anchors <- bpar[1:2]^2 }
	if (length(bpar)==9) { newbpar <- bpar }
	else if (length(bpar)==7) {
		if (flip=="A") { newpbar <- c(bpar[1:5],bpar[[6]],bpar[[7]],bpar[[6]],bpar[[6]]) }
		else if (flip=="B") { newbpar <- c(bpar[1:5],bpar[[6]],bpar[[6]],bpar[[7]],bpar[[6]]) }
		else { newbpar <- c(bpar[1:5],bpar[[6]],bpar[[6]],bpar[[6]],bpar[[7]]) }
	} else if (length(bpar)==8) {
		if (flip=="A") {
			if (bpar[[6]]!=bpar[[7]]) {
				if ((bpar[[6]]<bpar[[7]] && bpar[[8]]>bpar[[7]]) ||
					(bpar[[6]]>bpar[[7]] && bpar[[8]]<bpar[[7]])) {
					stop("Not a valid flipped BRAID parameter vector.")
				}
				bsign <- sign(bpar[[7]]-bpar[[6]])
				EfB <- bsign*min(bsign*bpar[c(6,8)])
			} else {
				EfB <- bpar[[7]]
			}
			newbpar <- c(bpar[1:7],EfB,bpar[8])
		} else if (flip=="B") {
			if (bpar[[6]]!=bpar[[7]]) {
				if ((bpar[[6]]<bpar[[7]] && bpar[[8]]>bpar[[7]]) ||
					(bpar[[6]]>bpar[[7]] && bpar[[8]]<bpar[[7]])) {
					stop("Not a valid flipped BRAID parameter vector.")
				}
				bsign <- sign(bpar[[7]]-bpar[[6]])
				EfA <- bsign*min(bsign*bpar[c(6,8)])
			} else {
				EfA <- bpar[[7]]
			}
			newbpar <- c(bpar[1:6],EfA,bpar[7:8])
		} else {
			if (bpar[[6]]!=bpar[[8]]) {
				if ((bpar[[6]]<bpar[[8]] && bpar[[7]]>bpar[[8]]) ||
					(bpar[[6]]>bpar[[8]] && bpar[[7]]<bpar[[8]])) {
					stop("Not a valid flipped BRAID parameter vector.")
				}
				bsign <- sign(bpar[[8]]-bpar[[6]])
				E0 <- bsign*min(bsign*bpar[6:7])
			} else {
				E0 <- bpar[[7]]
			}
			newbpar <- c(bpar[1:5],E0,bpar[6:8])
		}
	} else {
		stop("Not a valid flipped BRAID parameter vector.")
	}

	flipParameterVector(newbpar,flip,anchors)
}
flipParameterVector <- function(bpar,flip,anchors) {
	newbpar <- pmax(bpar,c(0,0,0,0,-Inf,-Inf,-Inf,-Inf,-Inf))
	if (flip %in% c("A","both")) {
		newbpar[[1]] <- anchors[[1]]/newbpar[[1]]
	}
	if (flip %in% c("B","both")) {
		newbpar[[2]] <- anchors[[2]]/newbpar[[2]]
	}
	newbpar[flipIndices(flip)]
}
flipModelVector <- function(model,flip) {
	indicator <- rep(FALSE,9)
	indicator[model] <- TRUE
	indicator <- indicator[flipIndices(flip)]
	which(indicator)
}
flipIndices <- function(flip) {
	switch(
		flip,
		A = c(1:5,7,6,9,8),
		B = c(1:5,8,9,6,7),
		both = c(1:5,9,8,7,6),
		stop(sprintf("Unrecognized flip paramter '%s'.",flip))
	)
}
