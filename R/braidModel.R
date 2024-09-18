

#' Evaluate the BRAID response surface model
#'
#' @param DA A vector of concentrations of drug A in a combination (values 0
#' and `Inf` are permitted). Must be length 1 or the same length as `DB`.
#' @param DB A vector of concentrations of drug B in a combination (values 0
#' and `Inf` are permitted). Must be length 1 or the same length as `DA`.
#' @param bpar A BRAID response surface parameter vector (see Details)
#' @param calcderivs Primarily used by fitting functions for non-linear
#' optimization.  If `FALSE` (the default), the function returns a vector of
#' response values; if `TRUE`, it returns a list including the partial
#' derivatives of the BRAID parameters.
#'
#' @return If `calcderivs` is `FALSE`, a numeric vector the same length as `DA`
#' and/or `DB` with the predicted BRAID response surface values.  If
#' `calcderivs` is `TRUE`, a list with two elements: `value`, containing the
#' response surface values, and `derivatives`, a matrix with as many rows as
#' `value` has elements, and nine columns containing the partial derivatives of
#' the response surface with respect to the nine BRAID response surface
#' parameters
#' @export
#'
#' @details
#' The BRAID response model is, in total, described by nine response surface
#' parameters.  A BRAID parameter vector should uniquely determine all these
#' values. They are
#'
#' * IDMA: The dose of median effect of drug A, also called the EC50
#' * IDMB: The dose of median effect of drug B
#' * na: The Hill slope, or sigmoidicity, of drug A
#' * nb: The Hill slope of drug B
#' * kappa: The BRAID interaction parameter, indicating additivity (kappa = 0),
#' antagonism (2 < kappa < 0), or synergy (kappa > 0)
#' * E0: The minimal effect, the effect observed when neither drug is present
#' * EfA: The maximal effect of drug A, the effect theoretically observed when
#' drug B is absent and drug A is present at infinite concentration
#' * EfB: The maximal effect of drug B,
#' * Ef: The maximal effect of the combination, theoretically observed when
#' both drugs are present at infinite concentration. It may be (but often is
#' not) further from E0 than either EfA or EfB.
#'
#' In many cases, however, it is easier to specify only some of the final three
#' parameters.  [braidrm] functions therefore support BRAID parameter vectors
#' of length 7 (in which the sixth and seventh values are assumed to be E0 and
#' Ef, and EfA and EfB are assumed to be equal to Ef), length 8 (in which the
#' seventh and eighth values are EfA and EfB, and Ef is assumed to be equal to
#' whichever of these two values is further from E0), or the full length 9
#' parameter vector.
#'
#' @examples
#' concentrations <- c(0, 2^(-3:3))
#' surface <- data.frame(
#'     concA = rep(concentrations,each=length(concentrations)),
#'     concB = rep(concentrations,times=length(concentrations))
#' )
#'
#' surface$additive <- evalBraidModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, 0, 0, 100, 100, 100)
#' )
#'
#' surface$synergy <- evalBraidModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, 2, 0, 80, 90)
#' )
#'
#' surface$antagonism <- evalBraidModel(
#'     surface$concA,
#'     surface$concB,
#'     c(1, 1, 3, 3, -1, 0, 100)
#' )
#'
#' head(surface)
evalBraidModel <- function(DA,DB,bpar,calcderivs=FALSE) {
	if (length(DA)>1 && length(DB)>1 && length(DA)!=length(DB)) {
		stop("Paramters 'DA' and 'DB' must have length 1 or equal length.")
	}
	if (length(DA)==1 && length(DB)>1) {
		DA <- rep(DA,length(DB))
	} else if (length(DB)==1 && length(DA)>1) {
		DB <- rep(DB,length(DA))
	}
	fbpar <- fillOutBraidPar(bpar)
	parsf <- braidParToSF(fbpar)

	E0 <- fbpar[[6]]
	Ed <- fbpar[[9]]-E0
	result <- evalBraidModel_sf(DA,DB,parsf,calcderivs)
	if (!calcderivs) {
		return(E0+Ed*result)
	} else {
		value <- E0+Ed*result$value
		derivs <- cbind(Ed*result$derivs[,1:5,drop=FALSE],
						-result$value,
						result$derivs[,6:7,drop=FALSE],
						result$value)
		return(list(value=value,derivatives=derivs))
	}
}

#' Invert a BRAID Response Surface Model
#'
#' Given a particular effect and one of the doses in a combined action response
#' surface, this function calculates the other dose that will produce the
#' desired effect.  Used in the estimation of the IAE (see [estimateIAE()]) but
#' also useful for calculating something like the IC50 of one drug in the
#' presence of various doses of the other. `invertBraidModelA` and
#' `invertBraidModelB` are convenience wrapper functions that set `DA` or `DB`
#' to `NULL` to estimate the necessary concentrations of drug A and drug B
#' respectively.
#'
#' @param DA If not `NULL`, a vector of doses of drug A. Must be length 1 or
#' the same length as `effect`.  Only one of `DA` and `DB` may be not null.
#' @param DB If not `NULL`, a vector of doses of drug B. Must be length 1 or
#' the same length as `effect`.  Only one of `DA` and `DB` may be not null.
#' @param effect A vector of desired effect values to be reached.  Must be
#' length 1 or the same length as whichever of `DA` or `DB` is not null.
#' @param bpar A BRAID response surface parameter vector (see
#' [evalBraidModel()] for details)
#' @param invalidNA Specifies what to do with values that are outside the range
#' of the given BRAID model or doses.  If `FALSE` (the default), values "below"
#' the given range will be set to zero, and values "above" the given range will
#' be set to Inf.  If `TRUE`, all invalid values will be set to `NA`.
#' @param lowerBound Primarily used by [estimateIAE()].  If set to TRUE, will
#' return the lowest non-negative dose that produces an effect no greater than
#' the specified effect, rather than the highest
#'
#' @return A vector of concentrations the same length as either `DA` or `DB`
#' (whichever is not `NULL`) and/or `effect`, representing the concentration of
#' the other drug producing the specified effect in combination with the given
#' dose of the provided drug
#' @export
#'
#' @examples
#' baseIC <- invertBraidModel_A(
#'     DB=0,
#'     effect=seq(10,90,by=10),
#'     bpar=c(1, 1, 3, 3, 2, 0, 100, 100, 100)
#' )
#'
#' potentiatedIC <- invertBraidModel_A(
#'     DB=1,
#'     effect=seq(10,90,by=10),
#'     bpar=c(1, 1, 3, 3, 2, 0, 100, 100, 100)
#' )
invertBraidModel <- function(DA=NULL,DB=NULL,effect,bpar,invalidNA=FALSE,lowerBound=FALSE) {
	if (!is.null(DA)) {
		if (!is.null(DB)) { stop("Exactly one of 'DA' and 'DB must be null.") }
		if (length(DA)>1 && length(effect)>1 && length(DA)!=length(effect)) {
			stop("If 'DA' is not null, paramters 'DA' and 'effect' must have length 1 or equal length.")
		}
		if (length(DA)==1 && length(effect)>1) {
			DA <- rep(DA,length(effect))
		} else if (length(DA)>1 && length(effect)==1) {
			effect <- rep(effect,length(DA))
		}
	} else if (!is.null(DB)) {
		if (length(DB)>1 && length(effect)>1 && length(DB)!=length(effect)) {
			stop("If 'DB' is not null, paramters 'DB' and 'effect' must have length 1 or equal length.")
		}
		if (length(DB)==1 && length(effect)>1) {
			DB <- rep(DB,length(effect))
		} else if (length(DB)>1 && length(effect)==1) {
			effect <- rep(effect,length(DB))
		}
	} else { stop("Exactly one of 'DA' and 'DB must be null.") }

	fbpar <- fillOutBraidPar(bpar)
	parsf <- braidParToSF(fbpar)

	E0 <- fbpar[[6]]
	Ef <- fbpar[[9]]
	effect_sf <- (effect-E0)/(Ef-E0)
	lob <- effect_sf<0
	uob <- effect_sf>1
	effect_sf[lob] <- 0
	effect_sf[uob] <- 1

	if (invalidNA) { values <- c(NA,NA) }
	else { values <- c(0,Inf) }

	result <- invertBraidModel_sf(DA,DB,effect_sf,parsf,values,lowerBound)
	if (lowerBound) {
		result[lob] <- values[[2]]
		result[uob] <- values[[1]]
	} else {
		result[lob] <- values[[1]]
		result[uob] <- values[[2]]
	}
	result
}

#' @export
#' @rdname invertBraidModel
invertBraidModel_A <- function(DB,effect,bpar,invalidNA=FALSE,lowerBound=FALSE) {
	invertBraidModel(DA=NULL,DB=DB,effect,bpar,invalidNA=invalidNA,lowerBound=lowerBound)
}

#' @export
#' @rdname invertBraidModel
invertBraidModel_B <- function(DA,effect,bpar,invalidNA=FALSE,lowerBound=FALSE) {
	invertBraidModel(DA=DA,DB=NULL,effect,bpar,invalidNA=invalidNA,lowerBound=lowerBound)
}

fillOutBraidPar <- function(bpar) {
	if (length(bpar)>9) { stop("Not a valid BRAID parameter vector.") }
	else if (length(bpar)==9) { return(bpar) }
	else {
		if (length(bpar)<7) { stop("Not a valid BRAID parameter vector.") }
		if (length(bpar)==7) {
			return(c(bpar[1:6],rep(bpar[[7]],3)))
		} else {
			if (bpar[[6]]!=bpar[[7]]) {
				if ((bpar[[6]]<bpar[[7]] && bpar[[6]]>bpar[[8]]) ||
					(bpar[[6]]>bpar[[7]] && bpar[[6]]<bpar[[8]])) {
					stop("Not a valid BRAID parameter vector.")
				}
				bsign <- sign(bpar[[7]]-bpar[[6]])
				Ef <- bsign*max(bsign*bpar[7:8])
			} else {
				Ef <- bpar[[8]]
			}
			return(c(bpar,Ef))
		}
	}
}
braidParToSF <- function(bpar) {
	if (length(bpar)!=9) { stop("Must be a full BRAID parameter vector") }
	if (bpar[[9]]==bpar[[6]]) {
		if (bpar[[7]]!=bpar[[9]] || bpar[[8]]!=bpar[[9]]) {
			stop("Not a valid BRAID parameter vector.")
		}
		return(c(bpar[1:5],1,1))
	}
	FA <- (bpar[[7]]-bpar[[6]])/(bpar[[9]]-bpar[[6]])
	FB <- (bpar[[8]]-bpar[[6]])/(bpar[[9]]-bpar[[6]])
	if (FA<0 || FA>1 || FB<0 || FB>1) {
		stop("Not a valid BRAID parameter vector.")
	}
	return(c(bpar[1:5],FA,FB))
}

clip_lo <- function(n) {
	n[n < -.Machine$double.xmax] <- -.Machine$double.xmax
	n
}
clip_hi <- function(n) {
	n[n > .Machine$double.xmax] <- .Machine$double.xmax
	n
}
clip_finite <- function(n) {
	n[n < -.Machine$double.xmax] <- -.Machine$double.xmax
	n[n > .Machine$double.xmax] <- .Machine$double.xmax
	n
}
clip_positive <- function(n) {
	n[n < .Machine$double.xmin] <- .Machine$double.xmin
	n[n > .Machine$double.xmax] <- .Machine$double.xmax
	n
}
invertBraidModel_sf <- function(DA=NULL,DB=NULL,Rf,parsf,values=c(0,Inf),lowerBound=FALSE) {
	CA <- clip_finite(log(parsf[1]))  # (-Inf, Inf)
	CB <- clip_finite(log(parsf[2]))  # (-Inf, Inf)
	na <- clip_positive(parsf[3])  # (0, Inf)
	nb <- clip_positive(parsf[4])  # (0, Inf)
	kappa <- clip_hi(parsf[5])  # [-2, Inf)
	if (parsf[[6]] < 0) { FA <- 0 }
	else if (parsf[[6]] > 1 ) { FA <- 1 }
	else { FA <- parsf[[6]] } # [0, 1]
	if (parsf[[7]] < 0) { FB <- 0 }
	else if (parsf[[7]] > 1 ) { FB <- 1 }
	else { FB <- parsf[[7]] } # [0, 1]

	n <- clip_positive(sqrt(na)*sqrt(nb))  # (0, Inf)

	lob <- rep(FALSE,length(Rf))
	uob <- rep(FALSE,length(Rf))

	DTAB <- (Rf/(1-Rf))^(1/n)
	if (is.null(DA)) {
		lnDB <- CB-log(DB)  # [-Inf, Inf]
		pDB <- exp(nb*lnDB)  # [0, Inf]
		DTB <- FB/(1+pDB-FB)  # [0, Inf]
		CTB <- DTB^(1/n) # [0, Inf]

		RB <- DTAB/clip_positive(CTB) # [0, Inf]
		determinant <- kappa^2 - 4 + 4*RB
		determinant[determinant<0] <- 0
		if (kappa>0) {
			if (lowerBound) {
				uob <- uob | (RB<1)
				kappaRatio <- 0*RB
			} else {
				lob <- lob | (RB<1)
				kappaRatio <- ( -kappa + sqrt(determinant) )/2
			}
		} else {
			thresh <- 1-(kappa^2)/4
			if (lowerBound) {
				uob <- uob | (RB < thresh)
				lob <- lob | (RB > 1)
				kappaRatio <- ( -kappa - sqrt(determinant) )/2
			} else {
				lob <- lob | (RB < thresh)
				kappaRatio <- ( -kappa + sqrt(determinant) )/2
			}
		}

		kappaRatio[kappaRatio<0] <- 0
		kappaRatio <- clip_positive(kappaRatio^2) # (0, Inf)
		DTA <- ifelse(DB==0,DTAB*(!lowerBound),(CTB*kappaRatio)^n) # [0, Inf]
		if (lowerBound) { lob <- lob | (DTA > FA/(1-FA)) }
		else { uob <- uob | (DTA > FA/(1-FA)) }
		if (FA==1) { CCA <- DTA }
		else {
			CCA <- (FA/(FA-(1-FA)*DTA) - 1)/(1-FA)
		}
		DA <- clip_positive(exp(CA))*(CCA^(1/na))
		result <- DA
	} else {
		lnDA <- CA-log(DA)  # [-Inf, Inf]
		pDA <- exp(na*lnDA)  # [0, Inf]
		DTA <- FA/(1+pDA-FA)  # [0, Inf]
		CTA <- DTA^(1/n) # [0, Inf]

		RA <- DTAB/clip_positive(CTA) # [0, Inf]
		determinant <- kappa^2 - 4 + 4*RA
		determinant[determinant<0] <- 0
		if (kappa>0) {
			if (lowerBound) {
				uob <- uob | (RA < 1)
				kappaRatio <- 0*RA
			} else {
				lob <- lob | (RA < 1)
				kappaRatio <- ( -kappa + sqrt(determinant) )/2
			}
		} else {
			thresh <- 1-(kappa^2)/4
			if (lowerBound) {
				uob <- uob | (RA < thresh)
				# lob <- lob | (RA > 1)
				kappaRatio <- ( -kappa - sqrt(determinant) )/2
			} else {
				lob <- lob | (RA < thresh)
				kappaRatio <- ( -kappa + sqrt(determinant) )/2
			}
		}

		kappaRatio[kappaRatio<0] <- 0
		kappaRatio <- clip_positive(kappaRatio^2) # (0, Inf)
		DTB <- ifelse(DA==0,DTAB*(!lowerBound),(CTA*kappaRatio)^n) # [0, Inf]
		if (lowerBound) { lob <- lob | (DTB > FB/(1-FB)) }
		else { uob <- uob | (DTB > FB/(1-FB)) }
		if (FB==1) { CCB <- DTB }
		else {
			CCB <- (FB/(FB-(1-FB)*DTB) - 1)/(1-FB)
		}
		DB <- clip_positive(exp(CB))*(CCB^(1/nb))
		result <- DB
	}

	result[lob] <- values[[1]]
	result[uob] <- values[[2]]
	result
}

evalBraidModel_sf <- function(DA,DB,parsf,calcderivs=FALSE) {
	CA <- clip_finite(log(parsf[1]))  # (-Inf, Inf)
	CB <- clip_finite(log(parsf[2]))  # (-Inf, Inf)
	na <- clip_positive(parsf[3])  # (0, Inf)
	nb <- clip_positive(parsf[4])  # (0, Inf)
	kappa <- clip_hi(parsf[5])  # [-2, Inf)
	FA <- parsf[6]
	FB <- parsf[7]
	if (FA < 0 ) { FA <- 0 } else if (FA > 1) { FA <- 1} # [0 ,1]
	if (FB < 0 ) { FB <- 0 } else if (FB > 1) { FB <- 1} # [0 ,1]

	n <- clip_positive(sqrt(na)*sqrt(nb))  # (0, Inf)

	lnDA <- CA-log(DA)  # [-Inf, Inf]
	lnDB <- CB-log(DB)  # [-Inf, Inf]
	pDA <- exp(na*lnDA)  # [0, Inf]
	pDB <- exp(nb*lnDB)  # [0, Inf]

	DTA <- FA/(1+pDA-FA)  # [0, Inf]
	DTB <- FB/(1+pDB-FB)  # [0, Inf]

	pDTA <- DTA^(1/n)  # [0, Inf]
	pDTB <- DTB^(1/n)  # [0, Inf]
	rDAB <- clip_hi(sqrt(clip_hi(pDTA))*sqrt(clip_hi(pDTB)))  # [0, Inf)

	hDTA <- pDTA+clip_lo((kappa/2)*rDAB)  # (-Inf, Inf]
	hDTB <- pDTB+clip_lo((kappa/2)*rDAB)  # (-Inf, Inf]
	DTAB <- hDTA+hDTB  # [0, Inf]

	DT <- DTAB^n  # [0, Inf]
	R0 <- 1/(1+DT)  # [0, 1]
	Rf <- 1-R0  # [0, 1]

	if (!calcderivs) { return(Rf) }

	IDMA <- clip_positive(parsf[1])  # (0, Inf)
	IDMB <- clip_positive(parsf[2])  # (0, Inf)
	Rp <- R0*Rf  # [0, 1]
	dScA <- Rp*(clip_hi(hDTA)/clip_positive(DTAB))  # (-Inf, Inf)
	dScB <- Rp*(clip_hi(hDTB)/clip_positive(DTAB))  # (-Inf, Inf)
	dRfdn <- -clip_finite(dScA*clip_finite(log(pDTA/clip_positive(DTAB)))+
						  	dScB*clip_finite(log(pDTB/clip_positive(DTAB))))  # (-Inf, Inf)

	dRfdIDMA <- -clip_finite(dScA*clip_hi(na*(clip_hi(DTA)*clip_hi(pDA))/clip_positive(IDMA*FA))) # (-Inf, Inf)
	dRfdIDMB <- -clip_finite(dScB*clip_hi(nb*(clip_hi(DTB)*clip_hi(pDB))/clip_positive(IDMB*FB))) # (-Inf, Inf)
	dRfdna <- clip_finite(-dScA*clip_finite(clip_finite(clip_finite(lnDA)*clip_hi(DTA))*clip_hi(pDA)/clip_positive(FA))+
						  	clip_finite(n/(2*na))*dRfdn) # (-Inf, Inf)
	dRfdnb <- clip_finite(-dScB*clip_finite(clip_finite(clip_finite(lnDB)*clip_hi(DTB))*clip_hi(pDB)/clip_positive(FB))+
						  	clip_finite(n/(2*nb))*dRfdn) # (-Inf, Inf)
	dRfdkappa <- clip_hi(n*(Rp*(rDAB/clip_positive(DTAB))))  # [0, Inf)
	dRfdFA <- clip_finite(dScA*clip_hi((1+DTA)/clip_positive(FA))) # (-Inf, Inf)
	dRfdFB <- clip_finite(dScB*clip_hi((1+DTB)/clip_positive(FB))) # (-Inf, Inf)

	# derivs <- cbind(dRfdIDMA, dRfdIDMB,
	# 				dRfdna,
	# 				dRfdnb,
	# 				dRfdkappa, dRfdFA, dRfdFB)
	derivs <- cbind(dRfdIDMA, dRfdIDMB, dRfdna, dRfdnb,
					dRfdkappa, dRfdFA, dRfdFB)
	return(list(value=Rf,derivatives=derivs))
}
