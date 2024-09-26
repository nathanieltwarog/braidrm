
#' URSA Response Surface Fitting
#'
#' Fits the universal response surface approach (URSA) model to the given
#' data (Greco, Park, and Rustum 1990)
#'
#' @inheritParams braidrm
#' @param lower An optional set of lower bounds on the seven URSA response
#' parameters.  Any values set to NA will be filled with default calculated
#' bounds.
#' @param upper An optional set of upper bounds on the seven URSA response
#' parameters.  Any values set to NA will be filled with default calculated
#' bounds.
#'
#' @return An object of class `braidAltFit` with the following values:
#'
#' * `concs`: The array of concentrations passed to the functions
#' * `act`: The vector of measurements associated with the given dose pairs
#' * `weights`: The vector of weights for the given measurements, set to 1 for
#' all measurements by default
#' * `method`: Specifying the alternate surface model being used (in this case
#' "URSA")
#' * `coefficients`: A length-seven parameter vector specifying the URSA
#' response surface
#' * `fitted.values`: The predicted response surface value for the given dose
#' pairs and best-fit response surface
#' * `residuals`: The difference between the predicted and measured values for
#' the given dose pairs, always equal to "measured minus predicted"
#' * `direction`: The direction value passed to the function
#' * `pbounds`': A 2-by-7 array of bounds on the URSA parameters used in the
#' parameter optimization
#'
#' @references
#' Greco, William R, Hyoung Sook Park, and Youcef M Rustum. 1990.
#' “Application of a New Approach for the Quantitation of Drug Synergism to
#' the Combination of Cis-Diamminedichloroplatinum and
#' 1-b-d-Arabinofuranosylcytosine.” *Cancer Research* **50** (17): 5318–27.
#'
#' @export
#'
#' @examples
#' ufit1 <- fitUrsaModel(measure ~ concA + concB, additiveExample)
#' coef(ufit1)
#'
#' ufit2 <- fitUrsaModel(measure ~ concA + concB, synergisticExample,
#'                       direction = 1, lower=c(NA, NA, NA, NA, NA, 0, 0))
#' coef(ufit2)
fitUrsaModel<- function(formula,data,weights=NULL,direction=0,lower=NULL,upper=NULL)
	UseMethod("fitUrsaModel")

#' @export
#' @rdname fitUrsaModel
fitUrsaModel.formula <- function(formula,data,weights=NULL,direction=0,lower=NULL,upper=NULL) {
	mf <- stats::model.frame(formula=formula, data=data)
	concs <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	ufit <- fitUrsaModel.default(concs,act,weights=weights,
							direction=direction,lower=lower,upper=upper)
	ufit$call <- match.call()
	return(ufit)

}

#' @export
#' @rdname fitUrsaModel
fitUrsaModel.default <- function(formula,data,weights=NULL,direction=0,lower=NULL,upper=NULL) {
	concs <- formula
	act <- data

	# Fill weights
	if (is.null(weights)) { weights <- rep(1,length(act)) }

	# Fill starting vector
	start <- defaultStartingVector(concs,act,c(1:5,6,9),"AB",NULL,direction)
	start <- start[c(1:5,6,9)]
	start[[5]] <- 0

	# Fill parameter bounds
	pbounds <- fillParameterBounds(lower,upper,c(1:5,6,9),concs,NULL)
	pbounds[,5] <- c(-0.9975,400)
	start <- pmax(pbounds[1,],start)
	start <- pmin(pbounds[2,],start)

	nstart <- start[1:5]
	nstart[1:4] <- log(nstart[1:4])
	nstart[[5]] <- log(nstart[[5]]+1)
	nbounds <- pbounds[,1:5]
	nbounds[,1:4] <- log(nbounds[,1:4])
	nbounds[,5] <- log(nbounds[,5]+1)

	obounds <- getOuterBounds(direction,pbounds[,6:7])

	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,obounds) {
		sfpar <- c(exp(parv[1:4]),exp(parv[[5]])-1,0,1)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalUrsaModel(concs[,1],concs[,2],sfpar)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,obounds)
		fpar <- c(sfpar[1:5],ebnds[1],ebnds[2])
		return(fpar)
	}
	# Define valderivfunc
	valfunc <- function(parv,concs,act,weights,start,obounds) {
		sfpar <- c(exp(parv[1:4]),exp(parv[[5]])-1,0,1)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalUrsaModel(concs[,1],concs[,2],sfpar)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,obounds)
		sfact <- (ebnds[2]-ebnds[1])*sfact+ebnds[1]-act
		ovalue <- sum((weights*sfact)^2)
		return(ovalue)
	}

	vfunc <- function(parv) valfunc(parv,concs,act,weights,start,obounds)
	fpfunc <- function(parv) par2fullpar(parv,concs,act,weights,start,obounds)

	nls <- stats::optim(nstart,vfunc,gr=NULL,method="L-BFGS-B",lower=nbounds[1,],
						upper=nbounds[2,],control=list(maxit=1000))

	fullpar <- fpfunc(nls$par)
	names(fullpar) <- c("IDMA","IDMB","na","nb","alpha","E0","Ef")
	fit <- evalUrsaModel(concs[,1],concs[,2],fullpar)

	structure(
		list(concs=concs,act=act,weights=weights,method="URSA",
			 coefficients=fullpar,fitted.values=fit,residuals=act-fit,
			 direction=direction,pbounds=pbounds),
		class="braidAltFit"
	)
}


#' Evaluate URSA response surface model
#'
#' Numerically estimates the universal response surface approach (URSA) model
#' for the given data and parameters (Greco, Park, and Rustum 1990).
#'
#' @inheritParams evalBraidModel
#' @param upar A length seven URSA response surface parameter vector (see
#' Details)
#'
#' @return A numeric vector the same length as `DA`
#' and/or `DB` with the predicted URSA response surface values.
#'
#' @details
#' The URSA model is described by the following seven values
#'
#' * IDMA: The dose of median effect of drug A, also called the EC50
#' * IDMB: The dose of median effect of drug B
#' * na: The Hill slope, or sigmoidicity, of drug A
#' * nb: The Hill slope of drug B
#' * alpha: The URSA interaction parameter, indicating additivity (alpha = 0),
#' antagonism (alpha < 0), or synergy (alpha > 0)
#' * E0: The minimal effect, the effect observed when neither drug is present
#' * Ef: The maximal effect of the drugs, theoretically observed when
#' either drug is present at infinite concentration
#'
#' @references
#' Greco, William R, Hyoung Sook Park, and Youcef M Rustum. 1990.
#' “Application of a New Approach for the Quantitation of Drug Synergism to
#' the Combination of Cis-Diamminedichloroplatinum and
#' 1-b-d-Arabinofuranosylcytosine.” *Cancer Research* **50** (17): 5318–27.
#'
#' @export
#'
#' @examples
#' concentrations <- c(0, 2^(-3:3))
#' surface <- data.frame(
#'     concA = rep(concentrations,each=length(concentrations)),
#'     concB = rep(concentrations,times=length(concentrations))
#' )
#'
#' surface$uadditive <- evalUrsaModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, 0, 0, 100)
#' )
#'
#' surface$usynergy <- evalUrsaModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, 5, 0, 80)
#' )
#'
#' surface$uantagonism <- evalUrsaModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, -0.5, 0, 100)
#' )
#'
#' head(surface)
evalUrsaModel <- function(DA,DB,upar) {
	if (length(DA)>1 && length(DB)>1 && length(DA)!=length(DB)) {
		stop("Paramters 'DA' and 'DB' must have length 1 or equal length.")
	}
	if (length(DA)==1 && length(DB)>1) {
		DA <- rep(DA,length(DB))
	} else if (length(DB)==1 && length(DA)>1) {
		DB <- rep(DB,length(DA))
	}

	sact <- rep(0.5*(upar[6]+upar[7]),length(DA))
	erng <- upar[7]-upar[6]
	sact[DA==0&DB==0] <- upar[6]
	sact[is.infinite(DA)|is.infinite(DB)] <- upar[7]
	rel <- DA==0 & is.finite(log(DB))
	sact[rel] <- upar[6]+(upar[7]-upar[6])/(1+(DB[rel]/upar[2])^(-upar[4]))
	rel <- is.finite(log(DA)) & DB==0
	sact[rel] <- upar[6]+(upar[7]-upar[6])/(1+(DA[rel]/upar[1])^(-upar[3]))
	rel <- is.finite(log(DA))  & is.finite(log(DB))
	for (iter in 1:25) {
		Ef <- (sact[rel]-upar[6])/(upar[7]-sact[rel])
		D1 <- Ef^(1/upar[3])
		D2 <- Ef^(1/upar[4])
		G <- DA[rel]/(upar[1]*D1) + DB[rel]/(upar[2]*D2) +
			upar[5]*DA[rel]*DB[rel]/(upar[1]*upar[2]*sqrt(D1*D2))
		sn <- rep(1,length(Ef))
		sn[G<1] <- -1
		sact[rel] <- sact[rel]+sn*erng*(2^(-iter-1))
	}
	return(sact)
}
