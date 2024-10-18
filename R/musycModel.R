
#' MuSyC Response Surface Fitting
#'
#' Fits the Multidimensional Synergy of Combinations (MuSyC) model of combined
#' action to the given data (Wooten *et al.* 2021).
#'
#' @inheritParams braidrm
#' @param variant String specifying which variant of the MuSyC model is to be
#' fit to the data.  If "standard" (the default), all MuSyC parameters except
#' `gamma12` and `gamma21` will be fit (these will be fixed at 1).  If
#' "independent", the four individual dose-response parameters (`IDMA`, `IDMB`,
#' `na`, and `nb`) and the four maximal effect parameters (`E0`, `EfA`, `EfB`
#' and `Ef`) will be fit, while the four interaction parameters (`alpha12`,
#' `alpha21`, `gamma12`, and `gamma12`) will all be fixed at 1.  If "full", the
#' full twelve-parameter MuSyC vector will be fit.
#' @param weights An optional vector of weights the same length as `act`.  If
#' `NULL` (the default), will be set to 1 for all measurements
#' @param lower An optional set of lower bounds on the fitted MuSyC response
#' parameters.  Any values set to NA will be filled with default calculated
#' bounds. May be length 4 (will be treated as a set of lower bounds on the
#' minimal and maximal effect parameters only), length 8 (will be treated as
#' lower bounds on the four individual dose response parameters and the four
#' minimal and maximla effect parameters), the same length as the space of
#' parameters being optimized (8 for "independent", 10 for "standard", or 12
#' for "full"), or length 12.
#' @param upper An optional set of lower bounds on the fitted MuSyC response
#' parameters.  Behaves the same as `lower`.
#'
#' @return An object of class `braidAltFit` with the following values:
#'
#' * `concs`: The array of concentrations passed to the functions
#' * `act`: The vector of measurements associated with the given dose pairs
#' * `weights`: The vector of weights for the given measurements, set to 1 for
#' all measurements by default
#' * `method`: Specifying the alternate surface model being used (in this case
#' "MuSyC")
#' * `variant`: A string specifying which MuSyC variant was fit: "independent",
#' "standard", or "full"
#' * `coefficients`: A parameter vector of the appropriate length for `variant`
#' specifying the best fit response surface
#' * `fitted.values`: The predicted response surface value for the given dose
#' pairs and best-fit response surface
#' * `residuals`: The difference between the predicted and measured values for
#' the given dose pairs, always equal to "measured minus predicted"
#' * `direction`: The direction value passed to the function
#' * `pbounds`': A 2-by-k array of bounds on the MuSyC parameters used in the
#' parameter optimization
#'
#' @references
#' Wooten, David J, Christian T Meyer, Alexander LR Lubbock, Vito Quaranta,
#' and Carlos F Lopez. 2021. “MuSyC Is a Consensus Framework That Unifies
#' Multi-Drug Synergy Metrics for Combinatorial Drug Discovery.”
#' *Nature Communications* **12** (1): 4607.
#'
#' @export
#'
#' @examples
#' mfit1 <- fitMusycModel(measure ~ concA + concB, synergisticExample)
#' coef(mfit1)
#'
#' mfit2 <- fitMusycModel(measure ~ concA + concB, oppositionalExample,
#'                        variant = "independent")
#' coef(mfit2)
fitMusycModel <- function(formula,data,variant="standard",weights=NULL,
								  direction=0,lower=NULL,upper=NULL)
	UseMethod("fitMusycModel")

#' @export
#' @rdname fitMusycModel
fitMusycModel.formula <- function(formula,data,variant="standard",weights=NULL,
								  direction=0,lower=NULL,upper=NULL) {
	mf <- stats::model.frame(formula=formula, data=data)
	concs <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	mfit <- fitMusycModel.default(concs,act,variant=variant,weights=weights,
								 direction=direction,lower=lower,upper=upper)
	mfit$call <- match.call()
	return(mfit)
}

#' @export
#' @rdname fitMusycModel
fitMusycModel.default <- function(formula,data,variant="standard",weights=NULL,
						  direction=0,lower=NULL,upper=NULL) {
	concs <- formula
	act <- data

	variantLength = switch(
		variant,
		independent = 8,
		standard = 10,
		full = 12,
		stop(sprintf("Unrecognized MuSyC variant '%s'.",variant))
	)

	# Fill weights
	if (is.null(weights)) { weights <- rep(1,length(act)) }

	# Fill parameter bounds
	if (!is.null(lower)) {
		if (length(lower)==4) { tlower <- c(NA,NA,NA,NA,NA,lower) }
		else if (length(lower)==8) { tlower <- c(lower[1:4],NA,lower[5:8]) }
		else if (length(lower)==12) { tlower <- c(lower[1:4],NA,lower[9:12]) }
		else if (length(lower)==variantLength) { tlower <- c(lower[1:4],NA,lower[7:10]) }
		else {
			stop("Parameter 'lower' must be of length 4 (bounding only the four",
				 " effect values), 12 (bounding all twelve full MuSyC parameters)",
				 " 8 (bounding all non-interaction parameters), or, if 'variant'",
				 " is 'standard', 10 (bounding all ten standard MuSyC parameters).")
		}
	} else { tlower <- NULL }
	if (!is.null(upper)) {
		if (length(upper)==4) { tupper <- c(NA,NA,NA,NA,NA,upper) }
		else if (length(upper)==8) { tupper <- c(upper[1:4],NA,upper[5:8]) }
		else if (length(upper)==12) { tupper <- c(upper[1:4],NA,upper[9:12]) }
		else if (length(upper)==variantLength) { tupper <- c(upper[1:4],NA,upper[7:10]) }
		else {
			stop("Parameter 'upper' must be of length 4 (bounding only the four",
				 " effect values), 12 (bounding all twelve full MuSyC parameters)",
				 " 8 (bounding all non-interaction parameters), or, if 'variant'",
				 " is 'standard', 10 (bounding all ten standard MuSyC parameters).")
		}
	} else { tupper <- NULL }
	pbounds <- fillParameterBounds(tlower,tupper,c(1:5,6,7,8,9),concs,NULL)
	pbounds <- pbounds[,c(1:4,6:9)]
	if (any(pbounds[2,]<pbounds[1,])) { stop("Unsatisfiable bounds,") }
	ebounds <- pbounds[,5:8]

	# Rectify direction
	if (any(pbounds[1,6:8]>pbounds[2,5])) { direction <- setDirection(1,direction) }
	if (any(pbounds[2,6:8]<pbounds[1,5])) { direction <- setDirection(-1,direction) }

	# Fill starting vector
	start <- defaultStartingVector(concs,act,c(1:5,6,9),"AB",NULL,direction)
	start <- start[c(1:4,6:9)]

	# Initialize non-linear starting vector and bounds
	if (variant=="independent") {
		nstart <- log(start[1:4])
		nbounds <- log(pbounds[,1:4])
	} else if (variant=="standard") {
		nstart <- log(c(start[1:4],start[1:2]))
		nbounds <- log(cbind(pbounds[,1:4],pbounds[,1:2]))
	} else if (variant=="full") {
		nstart <- log(c(start[1:4],start[1:4]))
		nbounds <- log(cbind(pbounds[,1:4],pbounds[,1:4]))
	}

	# Generate linear bounds on effects
	Wcoef <- array(0,dim=c(4,0))
	zval <- c()
	for (eind in 1:4) {
		wvec <- rep(0,4)
		wvec[[eind]] <- 1
		for (dind in 1:2) {
			sn <- sign(1.5-dind)
			if (is.finite(ebounds[dind,eind])) {
				Wcoef <- cbind(Wcoef,sn*wvec)
				zval <- c(zval,sn*ebounds[dind,eind])
			}
		}
	}
	if (direction!=0) {
		sn <- sign(direction)
		for (eind in 2:4) {
			wvec <- c(-sn,0,0,0)
			wvec[[eind]] <- sn
			Wcoef <- cbind(Wcoef,wvec)
			zval <- c(zval,0)
		}
	}
	if (length(zval)==0) {
		Wcoef <- NULL
		zval <- NULL
	}

	imat <- function(n) { diag(rep(1,n)) }
	variantMat <- switch(
		variant,
		independent = rbind(imat(4),imat(4)),
		standard = rbind(imat(6),cbind(array(0,dim=c(2,2)),imat(2),array(0,dim=c(2,2)))),
		full = imat(8),
		stop(sprintf("Unrecognized MuSyC variant '%s'.",variant))
	)
	variantMat <- variantMat[c(1,3,2,4,5,7,6,8),]

	dmats <- NULL
	parCQfunc <- function(nlpar) {
		mupar_bal <- exp(as.numeric(variantMat%*%cbind(nlpar)))

		mures <- evalMuSyC_balanced(concs[,1],concs[,2],mupar_bal,calcderivs=TRUE)

		wtmat <- diag(weights^2)
		Amat <- t(mures$value) %*% ( wtmat %*% mures$value )
		bvec <- t(mures$value) %*% ( wtmat %*% cbind(act) )
		cscalar <- sum((weights*act)^2)

		val <- constrainedQuadratic(Amat,bvec,cscalar,Wcoef,zval)
		derivs <- list()
		for (dindex in seq_along(nlpar)) {
			fAmat <- array(0,c(4,4))
			fbvec<- array(0,c(4,1))
			for (pindex in which(variantMat[,dindex]>0)) {
				Amat <- t(mures$derivatives[[pindex]]) %*% ( wtmat %*% mures$value )
				Amat <- Amat+t(Amat)
				fAmat <- fAmat+Amat
				bvec <- t(mures$derivatives[[pindex]]) %*% ( wtmat %*% cbind(act) )
				fbvec <- fbvec+bvec
			}
			derivs[[dindex]] <- constrainedQuadratic(fAmat*exp(nlpar[[dindex]]),
													 fbvec*exp(nlpar[[dindex]]),
													 0,NULL,NULL)
		}
		dmats <<- derivs
		return(val)
	}
	dCQfunc <- function(nlpar) {
		return(dmats)
	}
	nlns <- fitLNLfunction(nstart,parCQfunc,dCQfunc,
						   method="L-BFGS-B",lower=nbounds[1,],upper=nbounds[2,])


	mupar_bal <- exp(as.numeric(variantMat%*%cbind(nlns$par[seq_len(variantLength-4)])))
	mupar_sf <- c(mupar_bal[c(1,3,2,4)],
				  mupar_bal[c(3,1)]/mupar_bal[c(7,5)],
				  mupar_bal[c(8,6)]/mupar_bal[c(4,2)])
	mupar_full <- c(mupar_sf,nlns$par[(variantLength-3):variantLength])
	names(mupar_full) <- c("IDMA","IDMB","na","nb","alpha12","alpha21",
						   "gamma12","gamma21","E0","EfA","EfB","Ef")
	fullpar <- switch(
		variant,
		independent = mupar_full[c(1:4,9:12)],
		standard = mupar_full[c(1:6,9:12)],
		full = mupar_full,
		stop(sprintf("Unrecognized MuSyC variant '%s'.",variant))
	)

	fit <- evalMusycModel(concs[,1],concs[,2],mupar_full)

	structure(
		list(concs=concs,act=act,weights=weights,method="MuSyC",
			 variant=variant,coefficients=fullpar,fullvector=mupar_full,
			 fitted.values=fit,residuals=act-fit,
			 direction=direction,pbounds=pbounds),
		class="braidAltFit"
	)
}

#' Evaluate MuSyC Response Surfaces
#'
#' Evaluates the Mulitdimensional Synergy of Combinations (MuSyC) model of
#' combined action for the given values and parameters (Wooten *et al.* 2021).
#'
#' @inheritParams evalBraidModel
#' @param mupar A MuSyC response surface parameter vector; may be length 8, 10,
#' or 12 (see details for specifics of MuSyC parameters)
#'
#' @details
#' The multi-dimensional synergy of combinatoins, or MySyC, model is a
#' parametric response surface model introduced by Wooten et al. in 2021.  The
#' method models the effect of combination by simulating occupancy in four
#' compartments in which compounds are affected or unaffected by either drug.
#' The full MuSyC model can be specified by a total of twelve parameters:
#'
#' * `IDMA`: dose of median effect of first drug
#' * `IDMB`: dose of median effect of second drug
#' * `na`: Hill slope of first drug
#' * `nb`: Hill slope of second drug
#' * `alpha12`: factor by which first drug potentiates the second
#' * `alpha21`: factor by which second drug potentiates the first
#' * `gamma12`: factor by which first drug increases second drug's Hill slope
#' * `gamma21`: factor by which second drug increases first drug's Hill slope
#' * `E0` - the observed effect when unaffected by either drug
#' * `EfA` - the observed effect when affected by drug 1 but not drug 2
#' * `EfB` - the observed effect when affected by drug 2 but not drug 1
#' * `Ef` - the observed effect when affected by both drugs
#'
#' In practice, `gamma12` and `gamma21` are rarely used, so a ten-element
#' parameter vector specifies the other 10 values and assumes that `gamma12`
#' and `gamma21` are both equal to 1.  In some cases it is even useful to
#' specify a MuSyC surface with no interaction at all with an eight-element
#' vector, in which case `alpha12`, `alpha21`, `gamma12`, and `gamma21` are all
#' set equal to 1.
#'
#' @return If `calcderivs` is `FALSE`, a numeric vector the same length as `DA`
#' and/or `DB` with the predicted MuSyC response surface values.  If
#' `calcderivs` is `TRUE`, a list with two elements: `value`, containing the
#' response surface values, and `derivatives`, a matrix with as many rows as
#' `value` has elements, and all columns containing the partial derivatives of
#' the response surface with respect to the fitted MuSyC response surface
#' parameters
#'
#' @references
#' Wooten, David J, Christian T Meyer, Alexander LR Lubbock, Vito Quaranta,
#' and Carlos F Lopez. 2021. “MuSyC Is a Consensus Framework That Unifies
#' Multi-Drug Synergy Metrics for Combinatorial Drug Discovery.”
#' *Nature Communications* **12** (1): 4607.
#'
#' @export
#'
#' @examples
#' efficacyPar <- c(
#'     1, 1, 3, 3,
#'                       # Omitted shape synergy parameters assume to be 1
#'     0, 100, 100, 125  # Elevated Ef indicates efficacy synergy
#' )
#' potencyPar <- c(
#'     1, 1, 3, 3,
#'     10, 15,           # alphas above 1 indicate potency synergy
#'     0, 100, 100, 100  # No efficacy synergy
#' )
#'
#' concentrations <- c(0, 2^(-3:3))
#' surface <- data.frame(
#'     concA = rep(concentrations,each=length(concentrations)),
#'     concB = rep(concentrations,times=length(concentrations))
#' )
#' surface$efficacy <- evalMusycModel(surface$concA, surface$concB, efficacyPar)
#' surface$potency  <- evalMusycModel(surface$concA, surface$concB, potencyPar)
#'
#' head(surface)
evalMusycModel <- function(DA,DB,mupar,calcderivs=FALSE) {
	# mupar is a full 12-value MuSyC parameter vector. The first
	# eight values are a scale-free MuSyC parameter set described
	# above, and the remaining four represent the observed effects
	# in each of the four compartments:
	#
	# E0 - the observed effect when unaffected by either drug
	# E1 - the observed effect when affected by drug 1 but not drug 2
	# E2 - the observed effect when affected by drug 2 but not drug 1
	# E3 - the observed effect when affected by both drugs
	conc1 <- DA
	conc2 <- DB

	originalMupar <- mupar
	mupar <- fillOutMuSyCpar(mupar)

	mures <- piecesOfMuSyC(conc1,conc2,mupar[1:8],calcderivs)
	if (calcderivs) {
		muvalue <- mures$value
		muderivs <- mures$derivatives
	} else { muvalue <- mures }

	Evec <- cbind(mupar[9:12])
	vals <- as.numeric(muvalue%*%Evec)
	if (!calcderivs) { return(vals) }

	derivs <- array(0,dim=c(length(conc1),8))
	for (i in 1:8) {
		derivs[,i] <- muderivs[[i]]%*%Evec
	}
	derivs <- cbind(derivs,muvalue)

	if (length(originalMupar)==8) { derivs <- derivs[,c(1:4,9:12)] }
	else if (length(originalMupar)==10) { derivs <- derivs[,c(1:6,9:12)] }
	else if (length(originalMupar)==12) { derivs <- derivs }

	return(list(value=vals,derivatives=derivs))
}

piecesOfMuSyC <- function(conc1,conc2,mupar_sf,calcderivs=FALSE) {
	# mupar_sf is a scale free eight-parameter vector describing
	# the potency, sigmoidicity and potency interacion of the
	# two drugs. The values are:
	# IDMA - potency (or EC50) of first drug
	# IDMB - potency (or EC50) of second drug
	# na - Hill slope of first drug
	# nb - Hill slope of second drug
	# alpha1 - factor by which first drug potentiates the second
	# alpha2 - factor by which second drug potentiates the first
	# gamma1 - factor by which first drug increases second drug's Hill slope
	# gamma2 - factor by which second drug increases first drug's Hill slope
	#
	# There is an additional possible parameter:
	# zeta - ratio of second reverse rate to first
	# But its effects are small and not described in the MuSyC
	# paper, so it is assumed to be 1

	mupar_bal <- c(mupar_sf[c(1,3)],
				   mupar_sf[c(2,4)],
				   mupar_sf[[1]]/mupar_sf[[6]],mupar_sf[[3]]*mupar_sf[[8]],
				   mupar_sf[[2]]/mupar_sf[[5]],mupar_sf[[4]]*mupar_sf[[7]])
	balres <- evalMuSyC_balanced(conc1,conc2,mupar_bal,calcderivs)

	if (!calcderivs) { return(balres) }

	balderivs <- balres$derivatives

	sfderivs <- c(
		balderivs[c(1,3,2,4)],
		list (
			-balderivs[[7]]*mupar_bal[[7]]/mupar_sf[[5]],
			-balderivs[[5]]*mupar_bal[[5]]/mupar_sf[[6]],
			balderivs[[8]]*mupar_sf[[4]],
			balderivs[[6]]*mupar_sf[[3]]
		)
	)

	return(list(value=balres$value,derivatives=sfderivs))
}

evalMuSyC_balanced <- function(conc1,conc2,mupar_bal,calcderivs=FALSE) {
	X1 <- (conc1/mupar_bal[1])^(mupar_bal[2])
	Y1 <- (conc2/mupar_bal[3])^(mupar_bal[4])
	X2 <- (conc1/mupar_bal[5])^(mupar_bal[6])
	Y2 <- (conc2/mupar_bal[7])^(mupar_bal[8])

	flat <- array(1,dim=c(4,4))
	vmat <- cbind(
		2+X2+Y2,
		2*X1+(X1+Y1)*X2,
		2*Y1+(X1+Y1)*Y2,
		X1*Y2+Y1*X2+(X1+Y1)*X2*Y2
	)

	vmat_A <- cbind(0,1,0,Y2)
	vmat_B <- cbind(0,0,1,X2)
	vmat_AB <- cbind(0,0,0,rep(1,length(conc1)))
	relA <- is.infinite(X1+X2) & is.finite(Y1+Y2)
	relB <- is.finite(X1+X2) & is.infinite(Y1+Y2)
	relAB<- is.infinite(X1+X2) & is.infinite(Y1+Y2)
	vmat[relA,] <- vmat_A[relA,]
	vmat[relB,] <- vmat_B[relB,]
	vmat[relAB,] <- vmat_AB[relAB,]

	vtmat <- vmat%*%flat
	fmat <- vmat/vtmat
	if (!calcderivs) { return(fmat) }

	dvdxy <- list()
	dvdxy[[1]] <- cbind( 0, 2+X2, Y2, Y2+X2*Y2 )
	dvdxy[[1]][relA|relB|relAB,] <- 0
	dvdxy[[2]] <- cbind( 0, X2, 2+Y2, X2+X2*Y2 )
	dvdxy[[2]][relA|relB|relAB,] <- 0
	dvdxy[[3]] <- cbind( 1, X1+Y1, 0, Y1+(X1+Y1)*Y2 )
	dvdxy[[3]][relA|relB|relAB,] <- 0
	dvdxy[[3]][relB,4] <- 1
	dvdxy[[4]] <- cbind( 1, 0, X1+Y1, X1+(X1+Y1)*X2 )
	dvdxy[[4]][relA|relB|relAB,] <- 0
	dvdxy[[4]][relA,4] <- 1

	dfdxy <- list()
	for (i in 1:4) {
		dvtdxy <- dvdxy[[i]]%*%flat
		dfdxy[[i]] <- (vtmat*dvdxy[[i]]-vmat*dvtdxy)/(vtmat^2)
	}

	rep4 <- function(v) clip_finite(cbind(v,v,v,v))

	derivs <- list(
		dfdxy[[1]]*rep4(-mupar_bal[[2]]*X1/mupar_bal[[1]]),
		dfdxy[[1]]*rep4(clip_finite(log(conc1/mupar_bal[[1]]))*X1),
		dfdxy[[2]]*rep4(-mupar_bal[[4]]*Y1/mupar_bal[[3]]),
		dfdxy[[2]]*rep4(clip_finite(log(conc2/mupar_bal[[3]]))*Y1),
		dfdxy[[3]]*rep4(-mupar_bal[[6]]*X2/mupar_bal[[5]]),
		dfdxy[[3]]*rep4(clip_finite(log(conc1/mupar_bal[[5]]))*X2),
		dfdxy[[4]]*rep4(-mupar_bal[[8]]*Y2/mupar_bal[[7]]),
		dfdxy[[4]]*rep4(clip_finite(log(conc2/mupar_bal[[7]]))*Y2)
	)

	return(list(value=fmat,derivatives=derivs))
}

fillOutMuSyCpar <- function(mupar) {
	if (length(mupar)==12) { return(mupar) }
	else if (length(mupar)==10) { return(c(mupar[1:6],1,1,mupar[7:10])) }
	else if (length(mupar)==8) { return(c(mupar[1:4],1,1,1,1,mupar[5:8])) }
	else {
		stop("MuSyC parameter vectors must be of length 12 (specifying",
				" all surface parameters), length 10 (specifying both potency",
				" synergy parameters but not slope interaction parameters), or",
				" 8 (specifying neither potency synergy nor slope interaction",
				" parameters).")
	}
}

# start = a k-dimensional starting vector of non-linear parameters
# parCQfunc = a function which takes a length k non-linear parameter vector, and
# outputs an d-dimensional constrained quadratic
# dCQvfunc = an optional parameter that outputs a list of k d-dimensional
# constrained quadratics, which represent the derivative of the output CQ with
# respect to the k input parameters
fitLNLfunction <- function(start,parCQfunc,dCQfunc=NULL,...) {
	derivs <- NULL

	if (is.null(dCQfunc)) {
		parfunc <- function(nlpar) {
			minCQ <- fullyMinimizeConstrainedQuadratic(parCQfunc(nlpar))
			return(minCQ$Oval)
		}

		nls <- stats::optim(start,parfunc,gr=NULL,...)
	} else {
		derivs <- NULL
		parfunc <- function(nlpar) {
			parCQ <- parCQfunc(nlpar)
			dCQ <- dCQfunc(nlpar)
			minCQ <- fullyMinimizeConstrainedQuadratic(parCQ,dCQ)
			derivs <<- minCQ$Oderiv
			return(minCQ$Oval)
		}
		derivfunc <- function(nlpar) { return(derivs) }

		nls <- stats::optim(start,parfunc,gr=derivfunc,...)
	}

	minCQ <- fullyMinimizeConstrainedQuadratic(parCQfunc(nls$par))
	return(list(par=c(nls$par,as.numeric(minCQ$Evec)),
				fit=nls))
}

constrainedQuadratic <- function(Amat,bvec,cscalar,Wcoef=NULL,zval=NULL) {
	if (!is.numeric(Amat) || is.null(dim(Amat)) ||
		length(dim(Amat))!=2 || dim(Amat)[1]!=dim(Amat)[2] ||
		any(t(Amat)!=Amat)) {
		stop("Parameter 'Amat' must be a symmetric matrix.")
	}
	spaceDim <- dim(Amat)[1]

	if (!is.numeric(bvec) || (length(bvec)!=1 && length(bvec)!=spaceDim)) {
		stop("Parameter 'bvec' must be a vector of the same length as either",
			 " dimension of 'Amat'.")
	}
	if (length(bvec)==1) { bvec <- rep(bvec,spaceDim) }
	bvec <- cbind(as.numeric(bvec))

	if (!is.numeric(cscalar) || length(cscalar)!=1) {
		stop("Parameter 'cvec' must be a scalar number.")
	}

	if (!is.null(Wcoef)) {
		if (is.null(zval)) {
			stop("If parameter 'Wcoef' is specified, 'zval' must be also, and",
				 " vice versa.")
		}
		if (!is.numeric(Wcoef)) {
			stop("Parameter 'Wcoef' should be a vector of the same",
				 " length as either dimension of 'Amat' or a matrix/array",
				 " with as many rows as 'Amat'.")
		}
		if (is.null(dim(Wcoef))) {
			if (length(Wcoef)!=spaceDim) {
				stop("Parameter 'Wcoef' should be a vector of the same",
					 " length as either dimension of 'Amat' or a matrix/array",
					 " with as many rows as 'Amat'.")
			}
			Wcoef <- cbind(Wcoef)
		} else {
			if (length(dim(Wcoef))!=2 || dim(Wcoef)[1]!=spaceDim) {
				stop("Parameter 'Wcoef' should be a vector of the same",
					 " length as either dimension of 'Amat' or a matrix/array",
					 " with as many rows as 'Amat'.")
			}
		}
		constraintNum <- dim(Wcoef)[2]

		if (!is.numeric(zval) || length(zval)!=constraintNum) {
			stop("Parameter 'zval' must be a vector of the same length as the",
				 " number of columns (constraints) in 'Wcoef'.")
		}
		zval <- cbind(as.numeric(zval))

	} else if (!is.null(zval)) {
		stop("If parameter 'Wcoef' is specified, 'zval' must be also, and",
			 " vice versa.")
	}
	new_constrainedQuadratic(Amat,bvec,cscalar,Wcoef,zval)
}
new_constrainedQuadratic <- function(Amat,bvec,cscalar,Wcoef,zval) {
	structure(list(Amat=Amat,bvec=bvec,cscalar=cscalar,
				   Wcoef=Wcoef,zval=zval),
			  class="constrainedQuadratic")
}


fullyMinimizeConstrainedQuadratic <- function(cq,cqderivs=NULL) {
	if (!is.null(cqderivs) &&
		!all(vapply(cqderivs,function(d) checkConstrainedQuadraticDims(cq,d),FALSE))) {
		stop("Derivative quadratics must match the dimensions of the optimized
			 quadratic.")
	}

	spaceDim <- dim(cq$Amat)[1]
	if (is.null(cq$Wcoef)) {
		cq$Wcoef <- array(0,dim=c(spaceDim,0))
		cq$zval <- array(0,dim=c(0,1))
	}

	subModels <- list(list(set=c(),value=Inf))
	curInd <- 1
	while(curInd<=length(subModels)) {
		curset <- subModels[[curInd]]$set
		diff <- setdiff(seq_len(ncol(cq$Wcoef)),curset)
		curCQ <- subsetConstrainedQuadratic(cq,curset)

		minCQ <- minimizeConstrainedQuadratic(curCQ)
		if (length(diff)>0) {
			invalid <- t(cq$Wcoef[,diff,drop=FALSE])%*%cbind(minCQ$Evec) < cq$zval[diff]
			if (any(invalid)) {
				for (index in diff[invalid]) {
					nextset <- sort(unique(c(curset,index)))
					if (!any(vapply(subModels,function(m) checkSet(m$set,nextset),TRUE))) {
						subModels <- append(subModels,list(list(set=nextset,value=Inf)))
					}
				}
			} else {
				subModels[[curInd]]$value <- minCQ$Oval
			}
		} else {
			subModels[[curInd]]$value <- minCQ$Oval
		}
		curInd <- curInd+1
	}
	subModelScores <- vapply(subModels,function(m) m$value,Inf)
	bestInd <- which.min(subModelScores)
	bestSet <- subModels[[bestInd]]$set
	bestCQ <- subsetConstrainedQuadratic(cq,bestSet)

	if (is.null(cqderivs)) {
		return(minimizeConstrainedQuadratic(bestCQ))
	} else {
		bestCQderivs <- lapply(cqderivs, function(cqd) subsetConstrainedQuadratic(cqd,bestSet))
		return(minimizeConstrainedQuadraticDerivatives(bestCQ,bestCQderivs))
	}
}
checkConstrainedQuadraticDims <- function(cq1,cq2) {
	if (!all(dim(cq1$Amat)==dim(cq2$Amat))) { return(FALSE) }
	if (!all(dim(cq1$bvec)==dim(cq2$bvec))) { return(FALSE) }
	if (length(cq1$cscalar)!=length(cq2$cscalar)) { return(FALSE) }
	if (is.null(cq1$Wcoef)) {
		if (!is.null(cq2$Wcoef)) { return(FALSE) }
	} else {
		if (!is.null(cq2$Wcoef)) {
			if (!all(dim(cq1$Wcoef)==dim(cq2$Wcoef))) { return(FALSE)}
			if (!all(dim(cq1$zval)==dim(cq2$zval))) { return(FALSE)}
		}
	}
	return(TRUE)
}
subsetConstrainedQuadratic <- function(cq,set) {
	new_cq <- cq
	if (!is.null(cq$Wcoef)) {
		new_cq$Wcoef <- cq$Wcoef[,set,drop=FALSE]
		new_cq$zval <- cq$zval[set,,drop=FALSE]
	}
	new_cq
}
checkSet <- function(set1,set2) {
	if (length(set1)!=length(set2)) { FALSE }
	else { all(set1==set2) }
}

minimizeQuadratic <- function(cq) {
	Evec <- solve(cq$Amat,cq$bvec,tol=0)
	Oval <- cq$cscalar-as.numeric(t(cq$bvec)%*%Evec)
	return(list(Evec=Evec,Oval=Oval))
}
minimizeQuadraticDerivatives <- function(cq,cqderivs) {
	Evec <- solve(cq$Amat,cq$bvec,tol=0)
	Oval <- cq$cscalar-as.numeric(t(cq$bvec)%*%Evec)

	Oderiv <- rep(0,length(cqderivs))
	for (i in seq_along(cqderivs)) {
		cqd <- cqderivs[[i]]
		Oderiv[i] <-
			t(Evec)%*%(cqd$Amat%*%Evec) -
			2*t(Evec)%*%cqd$bvec +
			cqd$cscalar
	}
	return(list(Evec=Evec,Oval=Oval,Oderiv=Oderiv))
}
minimizeConstrainedQuadratic <- function(cq) {
	mq <- minimizeQuadratic(cq)
	if (is.null(cq$Wcoef)||length(cq$Wcoef)==0) { return(mq) }
	Aw <- solve(cq$Amat,tol=0)%*%cq$Wcoef
	altZ <- -(t(cq$Wcoef)%*%mq$Evec-cq$zval)
	wAwaltZ <- solve(t(cq$Wcoef)%*%Aw,altZ,tol=0)
	list(Evec=mq$Evec+Aw%*%wAwaltZ,
		 Oval=mq$Oval+t(altZ)%*%wAwaltZ)
}
minimizeConstrainedQuadraticDerivatives <- function(cq,cqderivs) {
	if (is.null(cq$Wcoef)||length(cq$Wcoef)==0) {
		return(minimizeQuadraticDerivatives(cq,cqderivs))
	} else { mq <- minimizeQuadratic(cq) }
	Aw <- solve(cq$Amat,tol=0)%*%cq$Wcoef
	altZ <- -(t(cq$Wcoef)%*%mq$Evec-cq$zval)
	wAwaltZ <- solve(t(cq$Wcoef)%*%Aw,altZ,tol=0)
	Evec <- as.numeric(mq$Evec+Aw%*%wAwaltZ)

	Oderiv <- rep(0,length(cqderivs))
	for (i in seq_along(cqderivs)) {
		cqd <- cqderivs[[i]]
		Oderiv[i] <-
			rbind(Evec)%*%(cqd$Amat%*%cbind(Evec)) -
			2*rbind(Evec)%*%cqd$bvec +
			cqd$cscalar
		if (is.null(cqd$Wcoef)||length(cqd$Wcoef)==0) { next }
		Oderiv[i] <- Oderiv[i] -
			2*rbind(Evec)%*%(cqd$Wcoef%*%wAwaltZ) +
			2*t(wAwaltZ)%*%cqd$zval
	}
	return(list(Evec=Evec,
				Oval=mq$Oval+t(altZ)%*%wAwaltZ,
				Oderiv=Oderiv))
}
