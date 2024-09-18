
#' Non-interacting Reference Surfaces
#'
#' Estimate the best fitting non-interacting reference surface according to
#' multiple methods, including Loewe, Bliss, and HSA
#'
#' @param concs A width-two array of concentrations representing all measured
#' dose pairs
#' @param act A vector of measured activity or effect values
#' @param method A string specifying which model of non-interaction should be
#' used; possible values are "Bliss" (the default), "HSA", "Loewe", and "ZIP"
#' @param ... Additional parameters to be passed to the method-specific
#' deviation or reference surface functions
#' @param range For Bliss calculations, the range of effects assumed by Bliss
#' independence; a  two-element numeric vector containing the minimal effect
#' and the maximal effect, in that order. For ZIP calculations, the initial
#' estimate of the minimal and maximal effects used in fitting the individual
#' dose response curves.
#' @param clip Clipping method for Bliss reference and deviation calculations.
#' Possible values are "pre", "post", and "none".  See details for specifics.
#' @param increasing For HSA calculations, is the effect increasing (TRUE)
#' meaning the "highest" single agent activity is numerically greater; or
#' decreasing (FALSE), meaning the "highest" single agent activity is
#' numerically lower. The latter may be appropriate when the modeled response is
#' target growth or survival.
#' @param limits For Loewe and ZIP calculations, the fixed values of the
#' minimal and maximal effects of the drugs and the combination. By default,
#' both values are set to `NA`; any value set to `NA` will fit from the data.
#'
#' @details
#' This collection of functions can be used to implement a family of combination
#' analysis methods known as "deviation" methods.  The details of the methods
#' differ, but the core strategy is common to all of them: estimate what a given
#' response measurement *would be* based on individual behaviors and some model
#' of non-interaction, and use the deviation of the measured response from that
#' expected response as a measure of the degree of synergy or antagonism.
#'
#' Bliss independence is the most widely used of these, and can be described as
#' the assumption that any response represents a fraction of the target
#' population being unaffected, and that a combined response corresponds to the
#' product of these fractions for both drugs (Bliss, 1939). It is extremely simple to
#' calculate, and relies on the intuitive model of probabilisticalyly independent
#' events. Because it treats responses as a scaled representation of
#' probabilities, it requires that all values be expressed relative to two
#' limiting values: the response seen when **all** targets are unaffected (the
#' minimal effect) and the response seen when **none** of the targets remain
#' unaffected (the maximal effect).  For a Bliss independent surface to be
#' estimated, these two values must be provided, using the parameter `range`.
#' Further, because values outside of this range would represent proportions
#' above 1 or below 0, most Bliss calculations involve some adjustment of the
#' data to ensure they always lie within the specified range.  The Bliss
#' functions support two ways of doing this to generate a reference surface:
#' "pre" will clip all values to the range immediately; "post" will clip all
#' calculated responses to the given range after they have been combined. A
#' third option, "none", performs no clipping at all, and allows for proportions
#' outside of the 0 to 1 range.
#'
#' The highest-single-agent, or HSA model, is even simpler than Bliss. The
#' effect of a combined pair of doses is simply the "greater" of the individual
#' effects produced by the two drugs at those levels.  The word "greater" is
#' placed in quotes here as taking the larger response value is only appropriate
#' when a numerically larger measurement corresponds to greater activity. For
#' survival or growth inhibition studies, the reverse may be true; the parameter
#' `increasing` allows the user to specify this directionality.
#'
#' Loewe additivity is one of the oldest models of non-interaction, and the
#' inspiration for BRAID additivity (Loewe and Muischnek, 1926).  According to
#' Loewe additivity, the combined response to a pair of doses the the effect
#' such that the two  individual doses represent complementary fractions of the
#' individual doses of both drugs that produce the same effect in isolation. It
#' is considered the gold standard of non-interaction by many researchers, but =
#' has several significant pragmatic disadvantages.  It requires that the full
#' dose response behavior of both drugs is known, and that they produce an
#' identical range of effects.  The `loeweReference` and `loeweDeviation`
#' functions perform basic dose response fitting with the additional constraint
#' of matching minimal and maximal effects, either fixed by the `limits`
#' parameter or estimated directly from the data.
#'
#' The zero-interaction potency, or ZIP, model is a variant of Bliss
#' independence that uses smoothing to give more robust values (Wooten *et al.*
#' 2015). The reference surface is calculated by fitting the dose response of
#' the individual drugs and then combining them using Bliss independence; the
#' method then adds additional robustness and smooths the measured surface
#' itself by fitting each constant-dose set of data points as its own dose
#' response curve, constrained to match the other drug's effect when the first
#' drug is zero.  Fitting these partial dose-response curves in either
#' direction produces a smoothed version of the original measurements (which
#' can be accessed directly using the function `zipSmoothed`), from which the
#' reference surface is subtracted to get the deviation (or "delta") surface.
#'
#' @return For the deviation functions (`deviationSurface`, `blissDeviation`,
#' `hsaDeviation`, `loeweDeviation`, and `zipDeviation`), a vector of values
#' the same length as `act` and/or `concs` representing the deviation of the
#' measurement from the specified reference surface.  For the reference
#' functions (`referenceSurface`, `blissReference`, `hsaReference`,
#' `loeweReference`, and `zipReferences`), a vector of values the same length
#' as `act` and/or `concs` containing the appropriate non-interacting reference
#' surface itself.  For `zipSmoothed`, the smoothed measurement surface given
#' by ZIP's dose-response-based smoothing method (see Details).
#'
#' @references
#' Bliss, Chester I. 1939. “The Toxicity of Poisons Applied Jointly 1.”
#' *Annals of Applied Biology* **26** (3): 585–615.
#'
#' Loewe, S, and H Muischnek. 1926. “Uber Kombinationswirkungen.”
#' Naunyn. Schmiedebergs. Arch. Pharmacol. 114: 313–26.
#'
#' Yadav, Bhagwan, Krister Wennerberg, Tero Aittokallio, and Jing Tang. 2015.
#' “Searching for Drug Synergy in Complex Dose–Response Landscapes Using an
#' Interaction Potency Model.”
#' *Computational and Structural Biotechnology Journal* **13**: 504–13.
#'
#' @export
#'
#' @examples
#' surface <- additiveExample
#' concs1 <- cbind(surface$concA,surface$concB)
#' act1 <- surface$measure
#'
#' sum(deviationSurface(concs1,act1,"Bliss",range=c(0,1)))
#' sum(deviationSurface(concs1,act1,"Loewe"))
#' surface$hsa <- hsaReference(concs1,act1,increasing=TRUE)
#'
#' surface <- synergisticExample
#' concs2 <- cbind(surface$concA,surface$concB)
#' act2 <- surface$measure
#'
#' sum(deviationSurface(concs2,act2,"ZIP",range=c(0,1)))
#' sum(deviationSurface(concs2,act2,"Loewe"))
#' surface$smooth <- zipSmoothed(concs2,act2,range=c(0,1))
deviationSurface <- function(concs,act,method="Bliss",...) {
	switch(
		method,
		Bliss = blissDeviation(concs,act,...),
		HSA   = hsaDeviation(concs,act,...),
		Loewe = loeweDeviation(concs,act,...),
		ZIP   = zipDeviation(concs,act,...),
		stop(sprintf("Unrecognized null reference model '%s'.",method))
	)
}

#' @export
#' @rdname deviationSurface
referenceSurface <- function(concs,act,method="Bliss",...) {
	switch(
		method,
		Bliss = blissReference(concs,act,...),
		HSA   = hsaReference(concs,act,...),
		Loewe = loeweReference(concs,act,...),
		ZIP   = zipReference(concs,act,...),
		stop(sprintf("Unrecognized null reference model '%s'.",method))
	)
}

#' @export
#' @rdname deviationSurface
blissDeviation <- function(concs,act,range,clip="none") {
	(act-blissReference(concs,act,range,clip))
}

#' @export
#' @rdname deviationSurface
blissReference <- function(concs,act,range,clip="none") {
	conc1 <- concs[,1]
	conc2 <- concs[,2]

	rawbconc1 <- conc1[conc1>0 & conc2==0]
	if (length(rawbconc1)==0) {
		stop("Surface must contain measurements of both drugs in the absence of the other")
	}
	rawact1 <- act[conc1>0 & conc2==0]
	bvalues1 <- fillValues(unique(conc1[conc1>0]),rawbconc1,rawact1)

	rawbconc2 <- conc2[conc2>0 & conc1==0]
	if (length(rawbconc2)==0) {
		stop("Surface must contain measurements of both drugs in the absence of the other")
	}
	rawact2 <- act[conc2>0 & conc1==0]
	bvalues2 <- fillValues(unique(conc2[conc2>0]),rawbconc2,rawact2)

	bvalues1$values <- (bvalues1$values-range[[1]])/(diff(range))
	bvalues2$values <- (bvalues2$values-range[[1]])/(diff(range))

	bact <- (act-range[[1]])/diff(range)

	if (clip=="pre") {
		bvalues1$values <- pmin(pmax(bvalues1$values,0),1)
		bvalues2$values <- pmin(pmax(bvalues2$values,0),1)
		bact <- pmin(pmax(bact,0),1)
	}
	for (i in seq_along(act)) {
		if (conc1[[i]]==0 || conc2[[i]]==0) { next }
		bact1 <- bvalues1$values[which(bvalues1$conc==conc1[[i]])]
		bact2 <- bvalues2$values[which(bvalues2$conc==conc2[[i]])]
		bact[[i]] <- bact1+bact2-(bact1*bact2)
	}
	if (clip=="post") {
		bact <- pmin(pmax(bact,0),1)
	}
	bact <- range[[1]] + diff(range)*bact
	bact
}

#' @export
#' @rdname deviationSurface
hsaDeviation <- function(concs,act,increasing) {
	(act-hsaReference(concs,act,increasing))
}

#' @export
#' @rdname deviationSurface
hsaReference <- function(concs,act,increasing) {
	if (increasing) { incr <- 1 } else { incr <- -1 }
	conc1 <- concs[,1]
	conc2 <- concs[,2]

	rawbconc1 <- conc1[conc1>0 & conc2==0]
	if (length(rawbconc1)==0) {
		stop("Surface must contain measurements of both drugs in the absence of the other")
	}
	rawact1 <- act[conc1>0 & conc2==0]
	bvalues1 <- fillValues(unique(conc1[conc1>0]),rawbconc1,rawact1)

	rawbconc2 <- conc2[conc2>0 & conc1==0]
	if (length(rawbconc2)==0) {
		stop("Surface must contain measurements of both drugs in the absence of the other")
	}
	rawact2 <- act[conc2>0 & conc1==0]
	bvalues2 <- fillValues(unique(conc2[conc2>0]),rawbconc2,rawact2)

	bact <- act
	for (i in seq_along(act)) {
		if (conc1[[i]]==0 || conc2[[i]]==0) { next }
		bact1 <- bvalues1$values[which(bvalues1$conc==conc1[[i]])]
		bact2 <- bvalues2$values[which(bvalues2$conc==conc2[[i]])]
		bact[[i]] <- incr*max(incr*bact1,incr*bact2)
	}
	bact
}

#' @inheritParams fitUrsaModel
#' @export
#' @rdname deviationSurface
loeweDeviation <- function(concs,act,weights=NULL,limits=c(NA,NA)) {
	(act-loeweReference(concs,act,weights,limits))
}

#' @export
#' @rdname deviationSurface
loeweReference <- function(concs,act,weights=NULL,limits=c(NA,NA)) {
	if (is.null(weights)) { weights <- rep(1,length(act)) }
	rel <- concs[,1]==0 | concs[,2]==0
	dhpar <- dualHillFit(concs[rel,1],concs[rel,2],act[rel],weights[rel],limits)

	upar <- c(dhpar[1:4],0,dhpar[5:6])
	evalUrsaModel(concs[,1],concs[,2],upar)
}

dualHillFit <- function(conc1,conc2,act,weights=NULL,limits=c(NA,NA)) {
	if (is.null(weights)) { weights <- rep(1,length(act)) }
	wzrel <- weights[conc1==0 & conc2==0]
	azrel <- act[conc1==0 & conc2==0]
	c1rel <- conc1[conc1>0 & conc2==0]
	w1rel <- weights[conc1>0 & conc2==0]
	a1rel <- act[conc1>0 & conc2==0]
	c2rel <- conc2[conc1==0 & conc2>0]
	w2rel <- weights[conc1==0 & conc2>0]
	a2rel <- act[conc1==0 & conc2>0]

	trwts <- c(wzrel,w1rel,w2rel)
	tract <- c(azrel,a1rel,a2rel)

	concs <- cbind(conc1,conc2)
	crng1 <- range(conc1[is.finite(log(conc1))])
	crng2 <- range(conc2[is.finite(log(conc2))])
	start <- c(exp(mean(log(crng1))),exp(mean(log(crng2))),1,1)
	pbounds <- fillParameterBounds(NULL,NULL,c(1:5,6,9),concs,NULL)

	nstart <- log(start[1:4])
	nbounds <- log(pbounds[,1:4])

	if (is.na(limits[[1]])) {
		if (is.na(limits[[2]])) {
			obounds <- basicdrm:::getHillOuterBounds(0,pbounds[,6:7])
			fpfunc <- function(parv) {
				sfact1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=FALSE)
				sfact2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=FALSE)
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				wt2 <- (trwts^2)/mean(trwts^2)

				mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*tract),mean(wt2*sfact*tract))
				ebnds <- basicdrm:::boundedOpt2d(mnv,obounds)
				fpar <- c(exp(parv),ebnds[1],ebnds[2])
			}
			# Define valderivfunc
			vdfunc <- function(parv) {
				sfres1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=TRUE)
				sfact1 <- sfres1$value
				sfres2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=TRUE)
				sfact2 <- sfres2$value
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				wt2 <- (trwts^2)/mean(trwts^2)

				mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*tract),mean(wt2*sfact*tract))
				ebnds <- basicdrm:::boundedOpt2d(mnv,obounds)
				ebnds[4:5] <- ebnds[4:5]*mean(trwts^2)
				sfact <- (ebnds[2]-ebnds[1])*sfact+ebnds[1]-tract
				ovalue <- sum((trwts*sfact)^2)
				sfact1 <- (ebnds[2]-ebnds[1])*sfact1+ebnds[1]-a1rel
				sfact2 <- (ebnds[2]-ebnds[1])*sfact2+ebnds[1]-a2rel

				derivs1 <- (2*(ebnds[2]-ebnds[1]))*as.vector(rbind(w1rel*w1rel*sfact1)%*%sfres1$derivatives)
				derivs2 <- (2*(ebnds[2]-ebnds[1]))*as.vector(rbind(w2rel*w2rel*sfact2)%*%sfres2$derivatives)
				derivs <- c(derivs1[[1]],derivs2[[1]],derivs1[[2]],derivs2[[2]])*exp(parv)

				return(list(value=ovalue,derivatives=derivs))
			}
		} else {
			obounds <- as.numeric(pbounds[,6])
			# Define par2fullpar
			fpfunc <- function(parv) {
				sfact1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=FALSE)
				sfact2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=FALSE)
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				wt2 <- (trwts^2)/mean(trwts^2)

				sfact <- 1-sfact
				mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*tract),mean(wt2*sfact*tract))
				ebnd <- basicdrm:::boundedOpt1d(mnv,limits[[2]],obounds)
				fpar <- c(exp(parv),ebnd[[1]],limits[[2]])
				return(fpar)
			}
			# Define valderivfunc
			vdfunc <- function(parv) {
				sfres1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=TRUE)
				sfact1 <- sfres1$value
				sfres2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=TRUE)
				sfact2 <- sfres2$value
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				wt2 <- (trwts^2)/mean(trwts^2)

				sfact <- 1-sfact
				mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*tract),mean(wt2*sfact*tract))
				sfact <- 1-sfact
				ebnd <- basicdrm:::boundedOpt1d(mnv,limits[[2]],obounds)
				ebnd[[3]] <- ebnd[[3]]*mean(trwts^2)
				sfact <- (limits[[2]]-ebnd[[1]])*sfact+ebnd[[1]]-act
				ovalue <- sum((trwts*sfact)^2)
				sfact1 <- (limits[[2]]-ebnd[[1]])*sfact1+ebnd[[1]]-a1rel
				sfact2 <- (limits[[2]]-ebnd[[1]])*sfact2+ebnd[[1]]-a2rel

				derivs1 <- (2*(limits[[2]]-ebnd[[1]]))*as.vector(rbind(w1rel*w1rel*sfact1)%*%sfres1$derivatives)
				derivs2 <- (2*(limits[[2]]-ebnd[[1]]))*as.vector(rbind(w2rel*w2rel*sfact2)%*%sfres2$derivatives)
				derivs <- c(derivs1[[1]],derivs2[[1]],derivs1[[2]],derivs2[[2]])*exp(parv)

				return(list(value=ovalue,derivatives=derivs))
			}
		}
	} else {
		if (is.na(limits[[2]])) {
			obounds <- as.numeric(pbounds[,7])
			# Define par2fullpar
			fpfunc <- function(parv) {
				sfact1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=FALSE)
				sfact2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=FALSE)
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				wt2 <- (trwts^2)/mean(trwts^2)

				mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*tract),mean(wt2*sfact*tract))
				ebnd <- basicdrm:::boundedOpt1d(mnv,limits[[1]],obounds)
				fpar <- c(exp(parv),limits[[1]],ebnd[[1]])
				return(fpar)
			}
			# Define valderivfunc
			vdfunc <- function(parv) {
				sfres1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=TRUE)
				sfact1 <- sfres1$value
				sfres2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=TRUE)
				sfact2 <- sfres2$value
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				wt2 <- (trwts^2)/mean(trwts^2)

				mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*tract),mean(wt2*sfact*tract))
				ebnd <- basicdrm:::boundedOpt1d(mnv,limits[[1]],obounds)
				ebnd[[3]] <- ebnd[[3]]*mean(trwts^2)
				sfact <- (ebnd[[1]]-limits[[1]])*sfact+limits[[1]]-tract
				ovalue <- sum((trwts*sfact)^2)
				sfact1 <- (ebnd[[1]]-limits[[1]])*sfact1+limits[[1]]-a1rel
				sfact2 <- (ebnd[[1]]-limits[[1]])*sfact2+limits[[1]]-a2rel

				derivs1 <- (2*(ebnd[[1]]-limits[[1]]))*as.vector(rbind(w1rel*w1rel*sfact1)%*%sfres1$derivatives)
				derivs2 <- (2*(ebnd[[1]]-limits[[1]]))*as.vector(rbind(w2rel*w2rel*sfact2)%*%sfres2$derivatives)
				derivs <- c(derivs1[[1]],derivs2[[1]],derivs1[[2]],derivs2[[2]])*exp(parv)

				return(list(value=ovalue,derivatives=derivs))
			}
		} else  {
			fpfunc <- function(parv) { fpar <- c(exp(parv),limits) }
			vdfunc <- function(parv) {
				sfres1 <- basicdrm:::evalHillModel_sf(c1rel,exp(parv[c(1,3)]),calcderivs=TRUE)
				sfact1 <- sfres1$value
				sfres2 <- basicdrm:::evalHillModel_sf(c2rel,exp(parv[c(2,4)]),calcderivs=TRUE)
				sfact2 <- sfres2$value
				sfact <- c(rep(0,length(azrel)),sfact1,sfact2)
				sfact <- (limits[2]-limits[1])*sfact+limits[1]-tract
				ovalue <- sum((trwts*sfact)^2)
				sfact1 <- (limits[2]-limits[1])*sfact1+limits[1]-a1rel
				sfact2 <- (limits[2]-limits[1])*sfact2+limits[1]-a2rel

				derivs1 <- (2*(limits[2]-limits[1]))*as.vector(rbind(w1rel*w1rel*sfact1)%*%sfres1$derivatives)
				derivs2 <- (2*(limits[2]-limits[1]))*as.vector(rbind(w2rel*w2rel*sfact2)%*%sfres2$derivatives)
				derivs <- c(derivs1[[1]],derivs2[[1]],derivs1[[2]],derivs2[[2]])*exp(parv)
				return(list(value=ovalue,derivatives=derivs))
			}
		}
	}
	# Run optim
	nls <- basicdrm:::runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	nls$fullpar
}

#' @export
#' @rdname deviationSurface
zipDeviation <- function(concs,act,range,weights=NULL,limits=c(NA,NA)) {
	zipfit <- zipFit_internal(concs,act,range,weights,limits)
	reference <- zipReference_internal(concs,act,range,weights,limits,zipfit)
	smoothed <- zipSmoothed_internal(concs,act,range,weights,limits,zipfit)
	(smoothed-reference)
}

#' @export
#' @rdname deviationSurface
zipReference <- function(concs,act,range,weights=NULL,limits=c(NA,NA)) {
	zipReference_internal(concs,act,range,weights,limits,NULL)
}

#' @export
#' @rdname deviationSurface
zipSmoothed <- function(concs,act,range,weights=NULL,limits=c(NA,NA)) {
	zipSmoothed_internal(concs,act,range,weights,limits,NULL)
}

zipReference_internal <- function(concs,act,range,weights=NULL,limits=c(NA,NA),zipfit=NULL) {
	if (is.null(zipfit)) {
		zipfit <- zipFit_internal(concs,act,range,weights,limits)
	}
	hpar1 <- zipfit$hpar1
	hpar2 <- zipfit$hpar2
	conc1 <- concs[,1]
	conc2 <- concs[,2]
	act1 <- basicdrm::evalHillModel(conc1,hpar1)
	act2 <- basicdrm::evalHillModel(conc2,hpar2)

	bact1 <- pmin(pmax((act1-range[[1]])/diff(range),0),1)
	bact2 <- pmin(pmax((act2-range[[1]])/diff(range),0),1)
	bact <- bact1+bact2-(bact1*bact2)
	act <- range[[1]]+diff(range)*bact
	act
}

zipSmoothed_internal <- function(concs,act,range,weights=NULL,limits=c(NA,NA),zipfit=NULL) {
	if (is.null(zipfit)) {
		zipfit <- zipFit_internal(concs,act,range,weights,limits)
	}
	hpar1 <- zipfit$hpar1
	hpar2 <- zipfit$hpar2
	conc1 <- concs[,1]
	conc2 <- concs[,2]
	if (is.null(weights)) { weights <- rep(1,length(act)) }

	sact1 <- basicdrm::evalHillModel(conc1,hpar1)
	for (cval2 in unique(conc2)) {
		if (cval2==0) { next }
		crel <- conc2==cval2
		e0 <- basicdrm::evalHillModel(cval2,hpar2)
		if (is.na(limits[[2]])) {
			start <- c(hpar1[1:2],e0,range[[2]])
			model <- c(1,2,4)
		} else {
			start <- c(hpar1[1:2],e0,limits[[2]])
			model <- c(1,2)
		}
		chfit1 <- basicdrm::fitHillModel(conc1[crel],act[crel],model,weights[crel],start)
		chpar1 <- stats::coef(chfit1)
		sact1[crel] <- basicdrm::evalHillModel(conc1[crel],chpar1)
	}
	sact2 <- basicdrm::evalHillModel(conc2,hpar2)
	for (cval1 in unique(conc1)) {
		if (cval1==0) { next }
		crel <- conc1==cval1
		e0 <- basicdrm::evalHillModel(cval1,hpar1)
		if (is.na(limits[[2]])) {
			start <- c(hpar2[1:2],e0,range[[2]])
			model <- c(1,2,4)
		} else {
			start <- c(hpar2[1:2],e0,limits[[2]])
			model <- c(1,2)
		}
		chfit2 <- basicdrm::fitHillModel(conc2[crel],act[crel],model,weights[crel],start)
		chpar2 <- stats::coef(chfit2)
		sact2[crel] <- basicdrm::evalHillModel(conc2[crel],chpar2)
	}

	(sact1+sact2)/2
}

zipFit_internal <- function(concs,act,range,weights,limits) {
	conc1 <- concs[,1]
	conc2 <- concs[,2]
	if (is.null(weights)) { weights <- rep(1,length(act)) }
	prel <- is.finite(log(conc1)) & is.finite(log(conc2))

	start1 <- c(exp(mean(range(log(conc1[prel])))),1,range)
	start2 <- c(exp(mean(range(log(conc2[prel])))),1,range)
	if (is.na(limits[[1]])) {
		if (is.na(limits[[2]])) { model <- c(1,2,3,4) }
		else {
			start1[[4]] <- limits[[2]]
			start2[[4]] <- limits[[2]]
			model <- c(1,2,3)
		}
	} else {
		if (is.na(limits[[2]])) {
			start1[[3]] <- limits[[1]]
			start2[[3]] <- limits[[1]]
			model <- c(1,2,4)
		} else {
			start1[3:4] <- limits
			start2[3:4] <- limits
			model <- c(1,2)
		}
	}
	hfit1 <- basicdrm::fitHillModel(conc1[conc2==0],act[conc2==0],model,weights[conc2==0],start1)
	hpar1 <- stats::coef(hfit1)
	hfit2 <- basicdrm::fitHillModel(conc2[conc1==0],act[conc1==0],model,weights[conc1==0],start2)
	hpar2 <- stats::coef(hfit2)
	return(list(hpar1=hpar1,hpar2=hpar2))
}

fillValues <- function(conc,rawconc,rawact) {
	rawact <- rawact[order(rawconc)]
	rawconc <- rawconc[order(rawconc)]

	rconc <- sort(unique(rawconc[rawconc>0]))
	ract <- rep(NA,length(rconc))
	for (i in seq_along(rconc)) {
		ract[[i]] <- mean(rawact[rawconc==rconc[[i]]])
	}

	values  <- rep(NA,length(conc))
	for (i in seq_along(conc)) {
		if (is.na(conc[[i]]) || conc[[i]]<=0) { next }
		if (is.infinite(conc[[i]])) {
			if (any(is.infinite(rconc))) {
				values[[i]] <- ract[[which(is.infinite(rconc))[[1]]]]
			} else {
				values[[i]] <- ract[[which.max(rconc)]]
			}
		} else {
			if (any(rconc==conc[[i]])) {
				values[[i]] <- ract[[which(rconc==conc[[i]])[[1]]]]
			} else if (conc[[i]]>max(rconc)) {
				values[[i]] <- ract[[which.max(rconc)]]
			} else if ((conc[[i]]<min(rconc))) {
				values[[i]] <- ract[[which.min(rconc)]]
			} else {
				lower <- which(rconc<=conc[[i]])
				i1 <- lower[[which.max(rconc[lower])]]
				upper <- which(rconc>=conc[[i]])
				i2 <- upper[[which.min(rconc[upper])]]
				if (rconc[[i1]]==rconc[[i2]]) {
					values[[i]] <- ract[[i1]]
				} else {
					conc1 <- log(rconc[[i1]])
					act1 <- ract[[i1]]
					conc2 <- log(rconc[[i2]])
					act2 <- ract[[i2]]
					values[[i]] <- (act1*(conc2-conc[[i]])+act2*(conc[[i]]-conc1))/(conc2-conc1)
				}
			}
		}
	}

	list(conc=conc,values=values)
}
