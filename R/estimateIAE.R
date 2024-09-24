
#' BRAID Response Surface Combined Potency
#'
#' Estimates the BRAID index of achievable efficacy (or IAE) for a given
#' response surface
#'
#' @param bpar The response surface to be evaluated.  Can be a numeric vector
#' (which will be treated as a standard BRAID parameter vector), a BRAID fit
#' object of class `braidrm`, or a flipped BRAID fit object of class
#' `braidrmflip`
#' @param levels The effect level or levels at which the index is to be
#' estimated
#' @param limits The upper concentration limits beneath which the IAE is to be
#' estimated.  Could be known toxicity thresholds, limits on pharmacokinetic
#' availability, or simply a convenient and consistent reference concentration
#' @param lowerLimits By default, the IAE is calculated by comparing the area
#' of dose space below which a given effect level is not reached with the total
#' achievable dose space specified by `limits`.  However, in some cases, it is
#' not desirable to allow the sub-threshold area to become infinitesimally
#' small. If lowerLimits is included, doses that lie below both lower limits
#' will always be included in the sub-threshold area, placing an effective upper
#' bound on possible IAE values
#'
#' @details
#' The index of achievable efficacy is an aggregate measure of combined potency,
#' and a useful first pass for quantifying the efficacy of a given response
#' surface.  Formally, it is equal to the area of achievable dose pairs divided
#' by the area of achievable doses below which a desired effect level is not
#' reached (then passed through a square root to give a more dimensionless
#' result).  If the surface is more potent, the area of sub-threshold achievable
#' doses is smaller, and the IAE is larger.  If the surface is less efficacious,
#' the doses at which a desired effect is reach will be larger, the
#' sub-threshold area will increase, and the IAE will decrease.  By default, the
#' IAE can range from 1 to infinity, but upper bounds can be placed by setting
#' the values in `lowerLimits` to concentrations above 0.  For convenience, the
#' function takes BRAID parameter vectors, `braidrm` fit objects, and
#' `braidrmflip` flipped BRAID fit objects.  However, flipped BRAID response
#' surface parameters cannot be passed to the function as is, so the function
#' `estimateFlippedIAE` is also included specifically for flipped parameter
#' vectors.
#'
#' @return A numeric vector, the same length as `levels` with the estimated
#' index of achievable efficacy (IAE) values for each of those levels.
#' @export
#'
#' @examples
#'
#' estimateIAE(c(1,1,3,3,0,0,100),c(50,90),c(5,5))
#'
#' bfit <- braidrm(measure ~ concA + concB, synergisticExample, getCIs = FALSE)
#' estimateIAE(bfit, c(0,0.25,0.5,0.75,1), c(10,1), lowerLimits=c(0.01,0.01))
estimateIAE <- function(bpar,levels,limits,lowerLimits=c(0,0))
	UseMethod("estimateIAE")

#' @export
#' @rdname estimateIAE
estimateIAE.braidrm <- function(bpar,levels,limits,lowerLimits=c(0,0)) {
	estimateIAE_internal(bpar$coefficients,levels,limits,lowerLimits)
}

#' @export
#' @rdname estimateIAE
estimateIAE.braidrmflip <- function(bpar,levels,limits,lowerLimits=c(0,0)) {
	estimateFlippedIAE_internal(bpar$coefficients,bpar$flip,levels,limits,lowerLimits)
}

#' @export
#' @rdname estimateIAE
estimateIAE.default <- function(bpar,levels,limits,lowerLimits=c(0,0)) {
	estimateIAE_internal(bpar,levels,limits,lowerLimits)
}

#' @param flip String specifying the direction or directions of the surface's
#' flip.  Must be one of "A", "B", or "both".
#' @export
#' @rdname estimateIAE
estimateFlippedIAE <- function(bpar,flip,levels,limits,lowerLimits=c(0,0)) {
	estimateFlippedIAE_internal(bpar,flip,levels,limits,lowerLimits)
}

estimateIAE_internal <- function(bpar,levels,limits,lowerLimits=c(0,0)) {
	bpar <- fillOutBraidPar(bpar)

	numpts <- 100

	if (length(limits)==1) {
		limits <- c(limits,limits)
	}
	if (length(lowerLimits)==1) {
		lowerLimits <- c(lowerLimits,lowerLimits)
	}

	iae <- rep(NA,length(levels))
	for (lind in seq_along(levels)) {
		level <- levels[[lind]]
		xybd_hi <- calculateXYBounds(bpar,level)

		# xbd_hi1 <- invertBraidModel(DB=0,effect=level,bpar=bpar,invalidNA=FALSE)
		xbd_hi1 <- xybd_hi[[1]]
		if (xbd_hi1==0) { area1 <- prod(lowerLimits) }
		else {
			xbd_hi1 <- min(max(xbd_hi1,lowerLimits[[1]]),limits[[1]])
			c1vec <- (seq_len(numpts)-0.5)*xbd_hi1/numpts
			c2out <- invertBraidModel(DA=c1vec,effect=level,bpar=bpar,invalidNA=FALSE)
			c2out[c2out>limits[[2]]] <- limits[[2]]
			c2out[c1vec<lowerLimits[[1]] & c2out<lowerLimits[[2]]] <- lowerLimits[[2]]
			if (bpar[[5]]<0) {
				c2out_lo <- invertBraidModel(DA=c1vec,effect=level,bpar=bpar,invalidNA=FALSE,lowerBound=TRUE)
				c2out_lo[c2out_lo>c2out] <- c2out[c2out_lo>c2out]
				c2out_lo[c1vec<lowerLimits[[1]]] <- 0
				c2out <- c2out-c2out_lo
			}
			area1 <- mean(c2out)*xbd_hi1
		}

		# ybd_hi2 <- invertBraidModel(DA=0,effect=level,bpar=bpar,invalidNA=FALSE)
		ybd_hi2 <- xybd_hi[[2]]
		if (ybd_hi2==0) { area2 <- prod(lowerLimits) }
		else {
			ybd_hi2 <- min(max(ybd_hi2,lowerLimits[[2]]),limits[[2]])
			c2vec <- (seq_len(numpts)-0.5)*ybd_hi2/numpts
			c1out <- invertBraidModel(DB=c2vec,effect=level,bpar=bpar,invalidNA=FALSE)
			c1out[c1out>limits[[1]]] <- limits[[1]]
			c1out[c2vec<lowerLimits[[2]] & c1out<lowerLimits[[1]]] <- lowerLimits[[1]]
			if (bpar[[5]]<0) {
				c1out_lo <- invertBraidModel(DB=c2vec,effect=level,bpar=bpar,invalidNA=FALSE,lowerBound=TRUE)
				c1out_lo[c1out_lo>c1out] <- c1out[c1out_lo>c1out]
				c1out_lo[c2vec<lowerLimits[[2]]] <- 0
				c1out <- c1out-c1out_lo
			}
			area2 <- mean(c1out)*ybd_hi2
		}

		iae[[lind]] <- sqrt((2*prod(limits))/(area1+area2))
	}
	iae
}

estimateFlippedIAE_internal <- function(bpar,flip,levels,limits,lowerLimits=c(0,0)) {
	bpar <- flippedFillOutBraidPar(bpar,flip)

	numpts <- 100

	if (length(limits)==1) {
		limits <- c(limits,limits)
	}
	if (length(lowerLimits)==1) {
		lowerLimits <- c(lowerLimits,lowerLimits)
	}

	iae <- rep(NA,length(levels))
	for (lind in seq_along(levels)) {
		level <- levels[[lind]]

		# xbd_hi1 <- invertFlippedBraidModel(DB=0,effect=level,bpar=bpar,flip=flip,invalidNA=FALSE)
		xbd_hi1 <- limits[[1]]
		if (xbd_hi1==0) { area1 <- prod(lowerLimits) }
		else {
			xbd_hi1 <- min(max(xbd_hi1,lowerLimits[[1]]),limits[[1]])
			c1vec <- (seq_len(numpts)-0.5)*xbd_hi1/numpts
			c2out <- invertFlippedBraidModel(DA=c1vec,effect=level,bpar=bpar,flip=flip,invalidNA=FALSE)
			c2out[c2out>limits[[2]]] <- limits[[2]]
			c2out[c1vec<lowerLimits[[1]] & c2out<lowerLimits[[2]]] <- lowerLimits[[2]]
			if (bpar[[5]]<0) {
				c2out_lo <- invertFlippedBraidModel(DA=c1vec,effect=level,bpar=bpar,flip=flip,invalidNA=FALSE,lowerBound=TRUE)
				c2out_lo[c2out_lo>c2out] <- c2out[c2out_lo>c2out]
				c2out_lo[c1vec<lowerLimits[[1]]] <- 0
				c2out <- c2out-c2out_lo
			}
			area1 <- mean(c2out)*xbd_hi1
		}

		# ybd_hi2 <- invertFlippedBraidModel(DA=0,effect=level,bpar=bpar,flip=flip,invalidNA=FALSE)
		ybd_hi2 <- limits[[2]]
		if (ybd_hi2==0) { area2 <- prod(lowerLimits) }
		else {
			ybd_hi2 <- min(max(ybd_hi2,lowerLimits[[2]]),limits[[2]])
			c2vec <- (seq_len(numpts)-0.5)*ybd_hi2/numpts
			c1out <- invertFlippedBraidModel(DB=c2vec,effect=level,bpar=bpar,flip=flip,invalidNA=FALSE)
			c1out[c1out>limits[[1]]] <- limits[[1]]
			c1out[c2vec<lowerLimits[[2]] & c1out<lowerLimits[[1]]] <- lowerLimits[[1]]
			if (bpar[[5]]<0) {
				c1out_lo <- invertFlippedBraidModel(DB=c2vec,effect=level,bpar=bpar,flip=flip,invalidNA=FALSE,lowerBound=TRUE)
				c1out_lo[c1out_lo>c1out] <- c1out[c1out_lo>c1out]
				c1out_lo[c2vec<lowerLimits[[2]]] <- 0
				c1out <- c1out-c1out_lo
			}
			area2 <- mean(c1out)*ybd_hi2
		}

		iae[[lind]] <- sqrt((2*prod(limits))/(area1+area2))
	}
	iae
}

calculateXYBounds <- function(bpar,level) {
	if (bpar[[5]]>=0) {
		return(c(
			invertBraidModel(DB=0,effect=level,bpar=bpar,invalidNA=FALSE),
			invertBraidModel(DA=0,effect=level,bpar=bpar,invalidNA=FALSE))
		)
	}
	na <- clip_positive(bpar[[3]])  # (0, Inf)
	nb <- clip_positive(bpar[[4]])  # (0, Inf)
	E0 <- bpar[[6]]
	FA <- (bpar[[7]]-bpar[[6]])/(bpar[[9]]-bpar[[6]])
	FB <- (bpar[[8]]-bpar[[6]])/(bpar[[9]]-bpar[[6]])
	Ef <- bpar[[9]]

	n <- clip_positive(sqrt(na)*sqrt(nb))  # (0, Inf)

	if (E0==Ef) { return(c(Inf,Inf)) }
	Rf <- (level-E0)/(Ef-E0)
	if (Rf>=1) { return(c(Inf,Inf)) }
	else if (Rf<=0) { return(c(0,0)) }

	kappa <- bpar[[5]]
	DAB <- Rf/(1-Rf)
	RT <- DAB/((1-(kappa^2)/4)^n)

	out <- c(Inf,Inf)
	if (FA>=1 || RT< (FA/(1-FA))) {
		RTA <- RT/(FA-(1-FA)*RT)
		out[[1]] <- bpar[[1]]*(RTA^(1/na))
	}
	if (FB>=1 || RT< (FB/(1-FB))) {
		RTB <- RT/(FB-(1-FB)*RT)
		out[[2]] <- bpar[[2]]*(RTB^(1/nb))
	}
	return(out)
}
