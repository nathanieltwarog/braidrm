
#' Combination Index
#'
#' Estimates the combination index using the median effet method of Chou and
#' Talalay (1984) or the more robust method of non-linear optimization.
#'
#' @inheritParams deviationSurface
#' @param level A numeric vector of one or more effect levels at which to
#' estimate the combination index
#' @param limits The fixed values of the minimal and maximal effects of the
#' drugs and the combination. By default, both values are set to `NA`; any
#' value set to `NA` will fit from the data.
#' @param dr1 A data frame with two columns, `conc` and `act` reflecting the
#' dose response behavior of the first drug alone
#' @param dr2 A data frame with two columns, `conc` and `act` reflecting the
#' dose response behavior of the second drug alone
#' @param drc A data frame with two columns, `conc` and `act` reflecting the
#' dose response behavior of a constant ratio combination; `conc` should be the
#' combined concentrations of the two drugs
#' @param ratio The ratio of the two drugs in the constant ratio combination
#' (dose A to dose B)
#' @param range The range of effects assumed by the median effect model; a
#' two-element numeric vector containing the minimal effect and the maximal
#' effect, in that order.
#' @param excess For `estimateChouIndices` and `estimateChouIndex`, what should
#' be done with values outside the expected range.  If "clip" (the default),
#' values will be clipped to 0.1% and 99.9% of the expected range; if "drop",
#' values outside the range will be dropped
#'
#' @details
#' The combination index is a peculiar method, as it does not produce values
#' corresponding to each measurement (as the deviation methods do), nor does it
#' produce a value shared by the entire surface (as response surface methods
#' do).  It instead, produces a value associated with a particular effect level
#' and a particular *dose ratio*.  This makes implementing the method
#' consistently for a wide range of possible data sources quite tricky;
#' nevertheless, we have attempted to do so here.  In brief, the combination
#' index method involves fitting the dose response of both individual drugs and
#' and a constant ratio combination of the drugs (treated as a virtual third
#' drug).  It the compares the potency of the constant-ratio combination
#' (estimated a particular effect level) with the expected potency according to
#' Loewe additivity, and returns the degree to which the combination is *more*
#' potent (with a combination index less than 1) or *less* potent (with a
#' combination greater than one) than expected by additivity.  This method is
#' also in turns known as the sum of fractional inhibitory coefficients (FICs),
#' observed-over-expected, or originally as Berenbaum's interaction index.
#'
#' Because the method was originally built for three distinct sets of dose
#' response measurements, we have included the `estimateCombinationIndex` and
#' `estimateChouIndex` functions which operate on three separate data frames.
#' However, in most cases, it will be easier to use the
#' `estimateCombinationIndices` and `estimateChouIndices` functions, which
#' operate on an array of concentrations and a vector of responses, just like
#' the numerous other functions in this package, and generate a set of
#' combination index values with level and ratio information  included.
#'
#' The only difference between the `estimateCombination*` and `estimateChou*`
#' functions is the way in which they perform dose response fitting.
#' The `estimateCombination*` functions use non-linear least squares
#' optimization (based on the package `basicdrm`) to estimate dose-response
#' parameters.  The `estimateChou*` functions use the median-effect method
#' described by Chou and Talalay in their 1984 paper, which linearized all
#' measurements and performs linear regression.  We *do not* recommend using
#' these methods, as they are much less reliable than the non-linear
#' optimization approach and extremely susceptible to noise at extreme values.
#'
#'
#' @return For `estimateCombinationIndices` and `estimateChouIndices`, a data
#' frame with the following columns:
#'
#' * `ratio`: The ratio of doses (dose A to dose B) along which the index
#' was estimated
#' * `level`: The effect level at which the index was estimated
#' * `ci`: The estimated combination index at that dose ratio and effect level
#'
#' Combination index estimates will be included for all provided effect levels
#' and all dose ratios present. For `estimateCombinationIndex` and
#' `estimateChouIndex`, a vector of estimated combinatoin indices the same
#' length as `level`.
#'
#' @references
#' Berenbaum, MC. 1978. “A Method for Testing for Synergy with Any Number of
#' Agents.” *Journal of Infectious Diseases* **137** (2): 122–30.
#'
#' Berenbaum, MC. 1989. “What Is Synergy.” *Pharmacol Rev* **41**: 93–141.
#'
#' Chou, Ting-Chao, and Paul Talalay. 1984. “Quantitative Analysis of
#' Dose-Effect Relationships: The Combined Effects of Multiple Drugs or Enzyme
#' Inhibitors.” *Advances in Enzyme Regulation* **22**: 27–55.
#'
#' @export
#'
#' @examples
#' surface <- synergisticExample
#' concs1 <- cbind(surface$concA, surface$concB)
#' act1 <- surface$measure
#'
#' estimateCombinationIndices(concs1,act1,c(0.5))
#'
#' dr1 <- surface[surface$concB==0, c("concA","measure")]
#' names(dr1) <- c("conc","act")
#' dr2 <- surface[surface$concA==0, c("concB","measure")]
#' names(dr2) <- c("conc","act")
#' drc <- surface[surface$concA==surface$concB,]
#' drc$conc <- drc$concA+drc$concB
#' drc <- drc[,c("conc","measure")]
#' names(drc) <- c("conc","act")
#'
#' estimateChouIndex(dr1,dr2,drc,ratio=1,
#'                   level=c(0.5,0.9,0.99),
#'                   range=c(0,1))
estimateCombinationIndices <- function(concs,act,level,weights=NULL,limits=c(NA,NA)) {
	c1rel <- concs[,2]==0
	c2rel <- concs[,1]==0
	cprel <- is.finite(log(concs[,1])) & is.finite(log(concs[,2]))
	dr1 <- data.frame(conc=concs[c1rel,1],act=act[c1rel])
	dr2 <- data.frame(conc=concs[c2rel,2],act=act[c2rel])

	if (length(which(c1rel & concs[,1]>0))<=1 ||
		length(which(c2rel & concs[,2]>0))<=1) {
		stop("Data must include at least two dose response measurments",
			 " for each drug in isolation (where the concentration of",
			 " the other drug is zero).")
	}
	if (is.null(weights)) { weights <- rep(1,length(act)) }

	c1vec <- concs[c1rel|c2rel,1]
	c2vec <- concs[c1rel|c2rel,2]
	actvec <- act[c1rel|c2rel]
	wtvec <- weights[c1rel|c2rel]
	dhpar <- dualHillFit(c1vec,c2vec,actvec,wtvec,limits=limits)

	hpar1 <- dhpar[c(1,3,5,6)]
	hpar2 <- dhpar[c(2,4,5,6)]
	range <- dhpar[c(5,6)]

	drcf <- data.frame(conc=concs[cprel,1]+concs[cprel,2],
					   act=act[cprel],
					   weight=weights[cprel],
					   ratio=concs[cprel,2]/concs[cprel,1])
	drcf$ratioNumber <- round(100*log(drcf$ratio))
	rnvalues <- sort(unique(drcf$ratioNumber))
	outdf <- data.frame()
	for (rind in seq_along(rnvalues)) {
		crel <- drcf$ratioNumber==rnvalues[[rind]]
		if (length(which(crel))<=1) { next }
		ratio <- exp(mean(log(drcf$ratio[crel])))

		civec <- estimateCombinationIndex_internal(hpar1,hpar2,drcf[crel,c("conc","act")],ratio,level,range)

		outdf <- rbind(outdf,data.frame(ratio=ratio,level=level,ci=civec))
	}
	if (nrow(outdf)==0) {
		stop("No dose ratio was measured at at least two concentrations, so",
			 " no combination index values could be estimated.")
	}

	outdf
}

#' @export
#' @rdname estimateCombinationIndices
estimateCombinationIndex <- function(dr1,dr2,drc,ratio,level,limits=c(NA,NA)) {
	if (is.null(dr1$weight)) { dr1$weight <- 1 }
	if (is.null(dr2$weight)) { dr2$weight <- 1 }
	if (is.null(drc$weight)) { drc$weight <- 1 }

	if (nrow(dr1)<=1 || nrow(dr2)<=1) {
		stop("Data must include at least two dose response measurments",
			 " for each drug in isolation (where the concentration of",
			 " the other drug is zero).")
	}

	c1vec <- c(dr1$conc,rep(0,nrow(dr2)))
	c2vec <- c(rep(0,nrow(dr1),dr2$conc))
	actvec <- c(dr1$act,dr2$act)
	wtvec <- c(dr1$weight,dr2$weight)
	dhpar <- dualHillFit(c1vec,c2vec,actvec,wtvec,limits=limits)

	hpar1 <- dhpar[c(1,3,5,6)]
	hpar2 <- dhpar[c(2,4,5,6)]
	range <- dhpar[c(5,6)]
	estimateCombinationIndex_internal(hpar1,hpar2,drc,ratio,level,range)
}

estimateCombinationIndex_internal <- function(hpar1,hpar2,drc,ratio,level,range) {
	if (any(hpar1[3:4]!=range) || any(hpar2[3:4]!=range)) {
		stop("Hill parameters must match range.")
	}

	civalues <- rep(NA,length(level))
	if (nrow(drc)<=1) { return(civalues) }

	start <- c(mean(c(hpar1[[1]],hpar2[[1]])),
			   exp(mean(c(log(hpar1[[2]]),log(hpar2[[2]])))),
			   range)
	hfitc <- basicdrm::fitHillModel(drc$conc,drc$act,weights=drc$weight,model=1:2,start=start)
	hparc <- stats::coef(hfitc)

	id1 <- basicdrm::invertHillModel(level,hpar1)
	id2 <- basicdrm::invertHillModel(level,hpar2)
	idc <- basicdrm::invertHillModel(level,hparc)
	for (li in seq_along(level)) {
		if (!all(is.finite(log(c(id1[[li]],id2[[li]],idc[[li]]))))) { next }
		idc1 <- idc[[li]]/(1+ratio)
		idc2 <- idc[[li]]*ratio/(1+ratio)

		civalues[[li]] <- (idc1/id1[[li]]) + (idc2/id2[[li]])
	}
	return(civalues)
}

#' @export
#' @rdname estimateCombinationIndices
estimateChouIndices <- function(concs,act,level,range,excess="clip") {
	c1rel <- is.finite(log(concs[,1])) & concs[,2]==0
	c2rel <- concs[,1]==0 & is.finite(log(concs[,2]>0))
	cprel <- is.finite(log(concs[,1])) & is.finite(log(concs[,2]))
	dr1 <- data.frame(conc=concs[c1rel,1],act=act[c1rel])
	dr2 <- data.frame(conc=concs[c2rel,2],act=act[c2rel])

	if (nrow(dr1)<=1 || nrow(dr2)<=1) {
		stop("Data must include at least two dose response measurments",
			 " for each drug in isolation (where the concentration of",
			 " the other drug is zero).")
	}

	drcf <- data.frame(conc=concs[cprel,1]+concs[cprel,2],act=act[cprel],
					   ratio=concs[cprel,2]/concs[cprel,1])
	drcf$ratioNumber <- round(100*log(drcf$ratio))
	rnvalues <- sort(unique(drcf$ratioNumber))
	outdf <- data.frame()
	for (rind in seq_along(rnvalues)) {
		crel <- drcf$ratioNumber==rnvalues[[rind]]
		if (length(which(crel))<=1) { next }
		ratio <- exp(mean(log(drcf$ratio[crel])))
		civec <- estimateChouIndex(dr1,dr2,drcf[crel,c("conc","act")],ratio,level,range,excess)
		outdf <- rbind(outdf,data.frame(ratio=ratio,level=level,ci=civec))
	}
	if (nrow(outdf)==0) {
		stop("No dose ratio was measured at at least two concentrations, so",
			 " no combination index values could be estimated.")
	}

	outdf
}

#' @export
#' @rdname estimateCombinationIndices
estimateChouIndex <- function(dr1,dr2,drc,ratio,level,range,excess="clip") {
	prepdr <- function(dr) {
		dr <- dr[is.finite(log(dr$conc)),]
		dr$norm <- (dr$act-range[[1]])/diff(range)
		if (excess=="clip") { dr$norm <- pmin(pmax(dr$norm,0.001),0.999) }
		else { dr1 <- dr[dr$norm>0 & dr$norm<1,] }
		dr$logconc <- log(dr$conc)
		dr$medeff <- log(dr$norm/(1-dr$norm))
		dr
	}

	dr1 <- prepdr(dr1)
	dr2 <- prepdr(dr2)
	drc <- prepdr(drc)

	civalues <- rep(NA,length(level))
	if (nrow(dr1)<=1 || nrow(dr2)<=1 || nrow(drc)<=1) { return(civalues) }

	medlevs <- (level-range[[1]])/(range[[2]]-level)
	medlevs[medlevs < 0] <- NA
	medlevs <- log(medlevs)

	lm1 <- stats::coef(stats::lm(medeff~logconc,dr1))
	lm2 <- stats::coef(stats::lm(medeff~logconc,dr2))
	lmc <- stats::coef(stats::lm(medeff~logconc,drc))

	for (mi in seq_along(medlevs)) {
		if (is.na(medlevs[[mi]])) { next }

		medlev <- medlevs[[mi]]
		id1 <- exp((medlev-lm1[[1]])/lm1[[2]])
		id2 <- exp((medlev-lm2[[1]])/lm2[[2]])
		idc <- exp((medlev-lmc[[1]])/lmc[[2]])
		idc1 <- idc/(1+ratio)
		idc2 <- idc*ratio/(1+ratio)

		civalues[[mi]] <- (idc1/id1) + (idc2/id2)
	}
	return(civalues)
}
