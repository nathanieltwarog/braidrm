
# E0, EfA, EfB, and Ef all fixed
fitBraidScenario_Z_1 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,FALSE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	if (start[[9]]>start[[6]]) { direction <- setDirection(1,direction) }
	else if (start[[9]]<start[[6]]) { direction <- setDirection(-1,direction) }
	else { stop("Parameter 'start` specifies a constant surface.") }

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	nstart <- bANDs$start
	nbounds <- bANDs$bounds


	# Define par2fullpar
	par2fullpar <- function(parv,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		return(c(sfpar,start[6:9]))
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- (start[7:8]-start[[6]])/(start[[9]]-start[[6]])
		sfpar <- c(sfpar,fpar)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (start[[9]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(start[[9]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5])]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,start,model)

	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"Z_1",concs,act,model,weights,start,direction,pbounds,kweight)
}

# EfA, EfB and Ef fixed, E0 varies below
fitBraidScenario_I_1 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,FALSE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	diffA <- start[[9]]-start[[7]]
	diffB <- start[[9]]-start[[8]]
	if (diffA==0 && diffB==0) {
		return(fitBraidScenario_I_1s(concs,act,model,weights,start,direction,pbounds,kweight))
	}
	if (diffA>0) { direction <- setDirection(1,direction) }
	else if (diffA<0) { direction <- setDirection(-1,direction) }
	if (diffB>0) { direction <- setDirection(1,direction) }
	else if (diffB<0) { direction <- setDirection(-1,direction) }
	anchorA <- abs(diffA)>=abs(diffB)
	if (anchorA) { anchorRatio <- diffB/diffA }
	else { anchorRatio <- diffA/diffB }

	# Rectify bounds
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[7]],start[[8]])
		if (anchorA) { ibounds <- (start[[9]]-start[[7]])/(start[[9]]-ebounds) }
		else { ibounds <- (start[[9]]-start[[8]])/(start[[9]]-ebounds) }
	} else if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[7]],start[[8]])
		if (anchorA) { ibounds <- (start[[9]]-start[[7]])/(start[[9]]-ebounds[c(2,1),,drop=FALSE]) }
		else { ibounds <- (start[[9]]-start[[8]])/(start[[9]]-ebounds[c(2,1),,drop=FALSE]) }
	}
	ibounds[1,1] <- max(ibounds[1,1],0.001)
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	if (anchorA) { istart <- (start[[9]]-start[[7]])/(start[[9]]-start[[6]]) }
	else { istart <- (start[[9]]-start[[8]])/(start[[9]]-start[[6]]) }
	istart <- min(max(istart,ibounds[1,1]),ibounds[2,1])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)

	# Define par2fullpar
	par2fullpar <- function(parv,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- parv[which(model>5)]
		if (anchorA) { E0 <- start[[9]]-(start[[9]]-start[[7]])/fpar }
		else { E0 <- start[[9]]-(start[[9]]-start[[8]])/fpar }
		sfpar <- c(sfpar,E0,start[7:9])
		return(sfpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar1 <- parv[which(model>5)]
		if (anchorA) {
			E0 <- start[[9]]-(start[[9]]-start[[7]])/fpar1
			fpar <- c(1-fpar1,1-(fpar1*anchorRatio))
		} else {
			E0 <- start[[9]]-(start[[9]]-start[[8]])/fpar1
			fpar <- c(1-(fpar1*anchorRatio),1-fpar1)
		}
		sfpar <- c(sfpar,fpar)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (start[[9]]-E0)*sfact+E0-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		Ederiv <- 2*sum(weights*weights*sfact*(1-sfres$value))/(fpar1^2)
		derivs <- (2*(start[[9]]-E0))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		if (anchorA) {
			derivs[6] <- -derivs[6]+Ederiv*(start[[9]]-start[[7]])
			derivs <- derivs[c(model[model<=5],6)]
		} else {
			derivs[7] <- -derivs[7]+Ederiv*(start[[9]]-start[[7]])
			derivs <- derivs[c(model[model<=5],7)]
		}
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,start,model)

	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"I_1",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0 varies freely, EfA, EfB, and Ef fixed at the *same* constant value
fitBraidScenario_I_1s <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,FALSE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}
	if (start[[7]]!=start[[9]] || start[[8]]!=start[[9]]) {
		stop("Scenario I_1s is only valid when all maximal effects are equal.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[1,1]>start[[9]]) { direction <- setDirection(-1,direction) }
	if (ebounds[2,1]<start[[9]]) { direction <- setDirection(1,direction) }
	if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[9]])
	} else if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[9]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Get (1) outer bound on E0
	base_obounds <- ebounds

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	nstart <- bANDs$start
	nbounds <- bANDs$bounds

	# Define parOuterBounds1
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,base_obounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(1,1)
		bounds <- base_obounds
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)
		mnv <- c(mean(wt2*(1-sfact)),mean(wt2*(1-sfact)^2),mean(wt2*act),mean(wt2*(1-sfact)*act))
		ebnd <- boundedOpt1d(mnv,start[[9]],bsf$bounds)
		fpar <- c(bsf$sfpar[1:5],ebnd[[1]],rep(start[[9]],3))
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		mnv <- c(mean(wt2*(1-sfact)),mean(wt2*(1-sfact)^2),mean(wt2*act),mean(wt2*(1-sfact)*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)
		sfact <- (start[[9]]-ebnd[[1]])*sfact+ebnd[[1]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(start[[9]]-ebnd[[1]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[model[model<=5]]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"I_1",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0, EfA, and EfB are fixed, Ef varies above
fitBraidScenario_I_2 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	diffA <- start[[7]]-start[[6]]
	diffB <- start[[8]]-start[[6]]
	if (diffA==0 && diffB==0) {
		stop("'model' parameter does not match scenario.")
	}
	if (diffA>0) { direction <- setDirection(1,direction) }
	else if (diffA<0) { direction <- setDirection(-1,direction) }
	if (diffB>0) { direction <- setDirection(1,direction) }
	else if (diffB<0) { direction <- setDirection(-1,direction) }
	anchorA <- abs(diffA)>=abs(diffB)
	if (anchorA) { anchorRatio <- diffB/diffA }
	else { anchorRatio <- diffA/diffB }

	# Rectify bounds
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[7]],start[[8]])
		if (anchorA) { ibounds <- (start[[7]]-start[[6]])/(ebounds[c(2,1),,drop=FALSE]-start[[6]]) }
		else { ibounds <- (start[[8]]-start[[6]])/(ebounds[c(2,1),,drop=FALSE]-start[[6]]) }
	} else if (direction<0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[7]],start[[8]])
		if (anchorA) { ibounds <- (start[[7]]-start[[6]])/(ebounds-start[[6]]) }
		else { ibounds <- (start[[8]]-start[[6]])/(ebounds-start[[6]]) }
	}
	ibounds[1,1] <- max(ibounds[1,1],0.001)
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	if (anchorA) { istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]]) }
	else { istart <- (start[[8]]-start[[6]])/(start[[9]]-start[[6]]) }
	istart <- min(max(istart,ibounds[1,1]),ibounds[2,1])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)

	# Define par2fullpar
	par2fullpar <- function(parv,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- parv[which(model>5)]
		if (anchorA) { Ef <- start[[6]]+(start[[7]]-start[[6]])/fpar }
		else { Ef <- start[[6]]+(start[[8]]-start[[6]])/fpar }
		sfpar <- c(sfpar,start[6:8],Ef)
		return(sfpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar1 <- parv[which(model>5)]
		if (anchorA) {
			Ef <- start[[6]]+(start[[7]]-start[[6]])/fpar1
			fpar <- c(fpar1,fpar1*anchorRatio)
		} else {
			Ef <- start[[6]]+(start[[8]]-start[[6]])/fpar1
			fpar <- c(fpar1*anchorRatio,fpar1)
		}
		sfpar <- c(sfpar,fpar)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (Ef-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		Ederiv <- 2*sum(weights*weights*sfact*sfres$value)/(fpar1^2)
		derivs <- (2*(start[[9]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		if (anchorA) {
			derivs[6] <- derivs[6]-Ederiv*(start[[7]]-start[[6]])
			derivs <- derivs[c(model[model<=5],6)]
		} else {
			derivs[7] <- derivs[7]-Ederiv*(start[[8]]-start[[6]])
			derivs <- derivs[c(model[model<=5],7)]
		}
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,start,model)

	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"I_2",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0, EfB, and Ef fixed, EfA varies between E0 and Ef
fitBraidScenario_I_3A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction and inner bound
	if (start[[9]]>start[[6]]) { direction <- setDirection(1,direction) }
	else if (start[[9]]<start[[6]]) { direction <- setDirection(-1,direction) }
	else { stop("Parameter 'start` specifies a constant surface.") }
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[1,1] <- pmax(ebounds[1,1],start[[6]])
		ebounds[2,1] <- pmin(ebounds[2,1],start[[9]])
		ibounds <- (ebounds-start[[6]])/(start[[9]]-start[[6]])
	} else if (direction<0) {
		ebounds[2,1] <- pmin(ebounds[2,1],start[[6]])
		ebounds[1,1] <- pmax(ebounds[1,1],start[[9]])
		ibounds <- (ebounds[c(2,1),,drop=FALSE]-start[[6]])/(start[[9]]-start[[6]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[1,1]),ibounds[2,1])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)

	# Define par2fullpar
	par2fullpar <- function(parv,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		sfpar <- c(sfpar,start[6:9])
		sfpar[7] <- parv[which(model>5)]*(start[[9]]-start[[6]])+start[[6]]
		return(sfpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[which(model>5)],(start[[8]]-start[[6]])/(start[[9]]-start[[6]]))
		sfpar <- c(sfpar,fpar)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (start[[9]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(start[[9]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,start,model)

	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"I_3A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0, EfA, and Ef fixed, EfB varies between E0 and Ef
fitBraidScenario_I_3B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_I_3A(concs[,c(2,1)],
									 act,
									 swapModelVector(model),
									 weights,
									 swapParameterVector(start),
									 direction,
									 swapBoundMat(pbounds,model),
									 kweight)
	swapBraidFitObject(flipfit,"I_3B")
}

# E0 and EfB fixed, Ef varies above, EfA fixed equal to Ef
fitBraidScenario_I_4A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	if (start[[8]]>start[[6]]) { direction <- setDirection(1,direction) }
	else if (start[[8]]<start[[6]]) { direction <- setDirection(-1,direction) }
	else { stop("Scenario allows effect of only one drug.") }

	# Rectify bounds
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[8]])
		ibounds <- (start[[8]]-start[[6]])/(ebounds[c(2,1),,drop=FALSE]-start[[6]])
	} else if (direction<0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[8]])
		ibounds <- (start[[8]]-start[[6]])/(ebounds-start[[6]])
	}
	ibounds[1,1] <- max(ibounds[1,1],0.001)
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[8]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[1,1]),ibounds[2,1])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)

	# Define par2fullpar
	par2fullpar <- function(parv,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- parv[which(model>5)]
		Ef <- start[[6]]+(start[[8]]-start[[6]])/fpar
		sfpar <- c(sfpar,start[[6]],Ef,start[[8]],Ef)
		return(sfpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar1 <- parv[which(model>5)]
		Ef <- start[[6]]+(start[[8]]-start[[6]])/fpar1
		fpar <- c(1,fpar1)
		sfpar <- c(sfpar,fpar)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (Ef-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		Ederiv <- 2*sum(weights*weights*sfact*sfres$value)/(fpar1^2)
		derivs <- (2*(start[[9]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs[7] <- derivs[7]-Ederiv*(start[[8]]-start[[6]])
		derivs <- derivs[c(model[model<=5],7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,start,model)

	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"I_4A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0 and EfA fixed, Ef varies above, EfB fixed equal to Ef
fitBraidScenario_I_4B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_I_4A(concs[,c(2,1)],
									 act,
									 swapModelVector(model),
									 weights,
									 swapParameterVector(start),
									 direction,
									 swapBoundMat(pbounds,model),
									 kweight)
	swapBraidFitObject(flipfit,"I_4B")
}

# E0 fixed, EfB fixed at constant, EfA varies freely, Ef fixed equal to maximum
fitBraidScenario_I_5A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	if (start[[8]]>start[[6]]) { direction <- setDirection(1,direction) }
	else if (start[[8]]<start[[6]]) { direction <- setDirection(-1,direction) }
	else { stop("Scenario allows effect of only one drug.") }
	ebounds <- pbounds[,which(model>5),drop=FALSE]

	checkModel1 <- checkModel2 <- TRUE
	if (direction>0) {
		if (ebounds[1,1]>=start[[8]]) { checkModel1 <- FALSE }
		else if (ebounds[2,1]<=start[[8]]) { checkModel2 <- FALSE }
	} else {
		if (ebounds[1,1]>=start[[8]]) { checkModel2 <- FALSE }
		else if (ebounds[2,1]<=start[[8]]) { checkModel1 <- FALSE }
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }

	if (checkModel1) {
		model1 <- c(model[which(model<=5)],7)
		start1 <- start
		if ((direction>0 && start[[7]]>start[[8]]) ||
			(direction<0 && start[[7]]<start[[8]])) {
			start1[[7]] <- start[[8]]
		}
		start1[[9]] <- start[[8]]
		bfit1 <- fitBraidScenario_I_3A(concs,act,model1,weights,start1,direction,pbounds,kweight=0)
	}
	if (checkModel2) {
		model2 <- c(model[which(model<=5)],9)
		start2 <- start
		if ((direction>0 && start[[7]]<start[[8]]) ||
			(direction<0 && start[[7]]>start[[8]])) {
			start2[[7]] <- start[[8]]
		}
		start2[[9]] <- start2[[7]]
		bfit2 <- fitBraidScenario_I_4A(concs,act,model2,weights,start2,direction,pbounds,kweight=0)
	}

	if (checkModel1 && checkModel2) {
		sse1 <- sum((bfit1$residuals*bfit1$weights)^2)
		sse2 <- sum((bfit2$residuals*bfit2$weights)^2)
		if (sse1<=sse2) { bfit <- bfit1 }
		else { bfit <- bfit2 }
	} else if (checkModel1) {
		bfit <- bfit1
	} else if (checkModel2) {
		bfit <- bfit2
	} else {
		stop("Unsatisfiable bounds.")
	}

	# bfit$scenario <- "I_5A"
	# bfit$model <- model
	# bfit$start <- start
	# bfit$pbounds <- pbounds
	bfit
}
# E0 fixed, EfA fixed at constant, EfB varies freely, Ef fixed equal to maximum
fitBraidScenario_I_5B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_I_5A(concs[,c(2,1)],
									 act,
									 swapModelVector(model),
									 weights,
									 swapParameterVector(start),
									 direction,
									 swapBoundMat(pbounds,model),
									 kweight)
	# swapBraidFitObject(flipfit,"I_5B")
	if (flipfit$scenario=="I_3A") {
		swapBraidFitObject(flipfit,"I_3B")
	} else if (flipfit$scenario=="I_4A") {
		swapBraidFitObject(flipfit,"I_4B")
	} else { stop("Invalid sub-scenario.")}
}

# E0 fixed, Ef varies freely, EfA and EfB fixed equal to Ef
fitBraidScenario_I_7 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[1,1]>start[[6]]) { direction <- setDirection(1,direction) }
	if (ebounds[2,1]<start[[6]]) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[6]])
	} else if (direction<0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[6]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Get (1) outer bound on Ef
	base_obounds <- ebounds

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	nstart <- bANDs$start
	nbounds <- bANDs$bounds

	# Define parOuterBounds1
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,base_obounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(1,1)
		bounds <- base_obounds
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)
		fpar <- c(bsf$sfpar[1:5],start[[6]],rep(ebnd[[1]],3))
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)
		sfact <- (ebnd[[1]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnd[[1]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[model[model<=5]]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"I_7",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0 and Ef vary, EfA and EfB fixed at  constant
fitBraidScenario_II_1 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	if (start[[7]]==start[[8]]) {
		return(fitBraidScenario_II_1s(concs,act,model,weights,start,direction,pbounds,kweight))
	}
	# Rectify direction
	erng <- start[7:8]
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[1,2]<erng[[2]] || ebounds[2,1]>erng[[1]]) {
		direction <- setDirection(1,direction)
	}
	if (ebounds[1,1]>erng[[1]] || ebounds[2,2]<erng[[2]]) {
		direction <- setDirection(-1,direction)
	}

	if (direction>=0) {
		start1 <- start
		start1[[6]] <- min(start1[[6]],erng[[1]])
		start1[[9]] <- max(start1[[9]],erng[[2]])
		bfit1 <- fitBraidScenario_II_1d(concs,act,model,weights,start1,direction=1,pbounds,kweight=0)
	}
	if (direction<=0) {
		start2 <- start
		start2[[6]] <- max(start2[[6]],erng[[2]])
		start2[[9]] <- min(start2[[9]],erng[[1]])
		bfit2 <- fitBraidScenario_II_1d(concs,act,model,weights,start2,direction=-1,pbounds,kweight=0)
	}

	if (direction==0) {
		sse1 <- sum((bfit1$residuals*bfit1$weights)^2)
		sse2 <- sum((bfit2$residuals*bfit2$weights)^2)
		if (sse1<=sse2) { bfit <- bfit1 }
		else { bfit <- bfit2 }
	} else if (direction>0) {
		bfit <- bfit1
	} else if (direction<0) {
		bfit <- bfit2
	}

	bfit$scenario <- "II_1"
	bfit$model <- model
	bfit$start <- start
	bfit$direction <- direction
	bfit$pbounds <- pbounds
	bfit
}
fitBraidScenario_II_1d <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[7]],start[[8]])
		ebounds[1,2] <- max(ebounds[1,2],start[[7]],start[[8]])
	} else if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[7]],start[[8]])
		ebounds[2,2] <- min(ebounds[2,2],start[[7]],start[[8]])
	} else {
		stop("Scenario II_1d is only valid when a direction is set.")
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }

	# Ensure finite effect bounds
	eminmax <- start[c(6,9)]
	eparscale <- abs(diff(eminmax))
	if (any(is.infinite(ebounds[1,]))) {
		loshot <- min(eminmax)-2*eparscale
		newEmin <- pmin(ebounds[2,],loshot)
		ebounds[1,is.infinite(ebounds[1,])] <- newEmin[is.infinite(ebounds[1,])]
	}
	if (any(is.infinite(ebounds[2,]))) {
		hishot <- max(eminmax)+2*eparscale
		newEmax <- pmax(ebounds[1,],hishot)
		ebounds[2,is.infinite(ebounds[2,])] <- newEmax[is.infinite(ebounds[2,])]
	}
	pbounds[,which(model>5)] <- ebounds

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	nstart <- c(bANDs$start,start[[6]],start[[9]])
	nbounds <- cbind(bANDs$bounds,ebounds)

	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		epar <- parv[(length(parv)-1):length(parv)]
		fpar <- (start[7:8]-epar[[1]])/(epar[[2]]-epar[[1]])
		fpar <- c(sfpar,epar[[1]],start[7:8],epar[[2]])
		return(fpar)
	}
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		epar <- parv[(length(parv)-1):length(parv)]
		fpar <- (start[7:8]-epar[[1]])/(epar[[2]]-epar[[1]])
		sfpar <- c(sfpar,fpar)

		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (epar[[2]]-epar[[1]])*sfact+epar[[1]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		delE <- epar[[2]]-epar[[1]]
		derivs <- (2*(epar[[2]]-epar[[1]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)

		ederivs <- c(2*sum(weights*weights*sfact*(1-sfres$value)),
					 2*sum(weights*weights*sfact*sfres$value))
		ederivs[[1]] <- ederivs[[1]]+(derivs[[6]]*(start[[7]]-epar[[2]])+
									  	derivs[[7]]*(start[[8]]-epar[[2]]))/(delE^2)
		ederivs[[2]] <- ederivs[[2]]-(derivs[[6]]*(start[[7]]-epar[[1]])+
									  	derivs[[7]]*(start[[8]]-epar[[1]]))/(delE^2)
		derivs[6:7] <- ederivs

		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model)

	# Run optim
	parscale <- rep(1,length(nstart))
	parscale[(length(nstart)-1):length(nstart)] <- eparscale
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds,parscale=parscale)

	braidFitObject(nls,"II_1",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0 and Ef vary, EfA and EfB fixed at the *same* constant
fitBraidScenario_II_1s <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}
	if (start[[7]]!=start[[8]]) {
		stop("Scenario I_1s is only valid when both agents' maximal effects are equal.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[2,1]<=start[[7]] || ebounds[1,2]>=start[[7]]) {
		direction <- setDirection(1,direction)
	}
	if (ebounds[1,1]>=start[[7]] || ebounds[2,2]<=start[[7]]) {
		direction <- setDirection(-1,direction)
	}
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[7]])
		ebounds[1,2] <- max(ebounds[1,2],start[[7]])
	} else if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[7]])
		ebounds[2,2] <- min(ebounds[2,2],start[[7]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,cbind(ebounds[,1],c(start[[7]],start[[7]]),
											  c(start[[8]],start[[8]]),ebounds[,2]))
	ibounds <- ibounds[,1]
	ibounds[[1]] <- max(ibounds[[1]],0.001)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[[1]]),ibounds[[2]])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],parv[length(parv)])
		fA <- fpar[[1]]
		if (fA>.Machine$double.eps) {
			bounds <- initialBound1d(start[[7]]*(1+1/fA)-rev(abbounds[,1])/fA,c("E0max","E0min"))
			if (fA<1) {
				bounds <- addBound1d(bounds,c(start[[7]]*(1-1/(1-fA))+abbounds[1,2]/(1-fA),1),"Efmin")
				bounds <- addBound1d(bounds,c(start[[7]]*(1-1/(1-fA))+abbounds[2,2]/(1-fA),-1),"Efmax")
			}
		} else {
			bounds <- initialBound1d(start[[7]]*(1-1/(1-fA))+abbounds[,2]/(1-fA),c("Efmin","Efmax"))
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,abbounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		fA <- bsf$sfpar[[6]]
		sfact <- sfact-fA
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- sfact+fA
		ebnd <- boundedOpt1d(mnv,start[[7]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		e0 <- (1+fA)*start[[7]]-fA*ebnd[[1]]
		ef <- (1-fA)*ebnd[[1]]+fA*start[[7]]
		fpar <- c(bsf$sfpar[1:5],e0,start[7:8],ef)
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,abbounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		fA <- bsf$sfpar[[6]]
		sfact <- sfact-fA
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- sfact+fA
		ebnd <- boundedOpt1d(mnv,start[[7]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)


		e0 <- (1+fA)*start[[7]]-fA*ebnd[[1]]
		ef <- (1-fA)*ebnd[[1]]+fA*start[[7]]

		sfact <- (ef-e0)*sfact+e0-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ef-e0))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		derivs[[6]] <- derivs[[6]]+derivs[[7]]
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }
			if (bname %in% c("E0min","E0max")) {
				bi <- which(c("E0min","E0max")==bname)
				derivs[[6]] <- derivs[[6]] - length(act)*ebnd[[3]]*(start[[7]]-abbounds[bi,1])/(fA^2)
			} else if (bname %in% c("Efmin","Efmax")) {
				bi <- which(c("Efmin","Efmax")==bname)
				derivs[[6]] <- derivs[[6]] + length(act)*ebnd[[3]]*(abbounds[bi,2]-start[[7]])/((1-fA)^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,abbounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"II_1s",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0 varies, EfB and Ef fixed, EfA varies
fitBraidScenario_II_2A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	if (start[[8]]<start[[9]]) { direction <- setDirection(1,direction) }
	else if (start[[8]]>start[[9]]) { direction <- setDirection(-1,direction) }
	else { return(fitBraidScenario_II_2As(concs,act,model,weights,start,direction,pbounds,kweight)) }
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction<0) {
		ebounds[1,2] <- max(ebounds[1,2],start[[9]])
		ebounds[1,1] <- max(c(ebounds[1,],start[[8]]))
		ebounds[2,2] <- min(ebounds[2,])
	} else if (direction>0) {
		ebounds[2,2] <- min(ebounds[2,2],start[[9]])
		ebounds[2,1] <- min(c(ebounds[2,],start[[8]]))
		ebounds[1,2] <- max(ebounds[1,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }

	# Ensure finite effect bounds
	eminmax <- start[c(6,9)]
	eparscale <- abs(diff(eminmax))
	if (any(is.infinite(ebounds[1,]))) {
		loshot <- min(eminmax)-2*eparscale
		newEmin <- pmin(ebounds[2,],loshot)
		ebounds[1,is.infinite(ebounds[1,])] <- newEmin[is.infinite(ebounds[1,])]
	}
	if (any(is.infinite(ebounds[2,]))) {
		hishot <- max(eminmax)+2*eparscale
		newEmax <- pmax(ebounds[1,],hishot)
		ebounds[2,is.infinite(ebounds[2,])] <- newEmax[is.infinite(ebounds[2,])]
	}
	pbounds[,which(model>5)] <- ebounds

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	if (direction>0) {
		if (start[[6]]>ebounds[1,2]) {
			fpar <- (start[[7]]-start[[6]])/(ebounds[2,2]-start[[6]])
		} else {
			fpar <- (start[[7]]-ebounds[1,2])/(ebounds[2,2]-ebounds[1,2])
		}
	} else {
		if (start[[6]]<ebounds[2,2]) {
			fpar <- (start[[7]]-ebounds[1,2])/(start[[6]]-ebounds[1,2])
		} else {
			fpar <- (start[[7]]-ebounds[1,2])/(ebounds[2,2]-ebounds[1,2])
		}
	}
	nstart <- c(bANDs$start,fpar,start[[6]])
	nbounds <- cbind(bANDs$bounds,c(0,1),ebounds[,1])
	abbounds <- ebounds[,2]

	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		e0 <- parv[[length(parv)]]
		fpar <- parv[[length(parv)-1]]
		if (direction>0) {
			if (e0>abbounds[[1]]) {
				efA <- e0+(abbounds[[2]]-e0)*fpar
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
			}
		} else {
			if (e0<abbounds[[2]]) {
				efA <- abbounds[[1]]+(e0-abbounds[[1]])*fpar
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
			}
		}
		fpar <- c(sfpar,e0,efA,start[[8]],start[[9]])
		return(fpar)
	}
	valderivfunc <- function(parv,concs,act,weights,start,model,abbounds,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		e0 <- parv[[length(parv)]]
		fpar <- parv[[length(parv)-1]]
		if (direction>0) {
			if (e0>abbounds[[1]]) {
				efA <- e0+(abbounds[[2]]-e0)*fpar
				escale <- (abbounds[[2]]-e0)
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
				escale <- (abbounds[[2]]-abbounds[[1]])
			}
		} else {
			if (e0<abbounds[[2]]) {
				efA <- abbounds[[1]]+(e0-abbounds[[1]])*fpar
				escale <- (e0-abbounds[[1]])
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
				escale <- (abbounds[[2]]-abbounds[[1]])
			}
		}
		fpar <- (c(efA,start[[8]])-e0)/(start[[9]]-e0)
		sfpar <- c(sfpar,fpar)

		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (start[[9]]-e0)*sfact+e0-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		delE <- start[[9]]-e0
		derivs <- (2*(start[[9]]-e0))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)

		derivs[[7]] <- (start[[8]]-start[[9]])*derivs[[7]]/(delE^2)
		derivs[[7]] <- derivs[[7]] + (efA-start[[9]])*derivs[[6]]/(delE^2)
		derivs[[7]] <- derivs[[7]] + 2*sum(weights*weights*sfact*(1-sfres$value))
		derivs[[6]] <- escale*derivs[[6]]/delE

		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,abbounds)

	# Run optim
	parscale <- rep(1,length(nstart))
	parscale[[length(nstart)]] <- eparscale
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds,parscale=parscale)

	braidFitObject(nls,"II_2A",concs,act,model,weights,start,direction,pbounds,kweight)

}
# E0 varies, EfB and Ef fixed at the *same* constant value, EfA varies
fitBraidScenario_II_2As <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}
	if (start[[8]]!=start[[9]]) {
		stop("Scenario II_2s is only valid when the fixed maximal effects are equal.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (any(ebounds[1,]>start[[9]])) { direction <- setDirection(-1,direction) }
	if (any(ebounds[2,]<start[[9]])) { direction <- setDirection(1,direction) }
	if (direction<0) {
		ebounds[1,2] <- max(ebounds[1,2],start[[9]])
		ebounds[1,1] <- max(ebounds[1,])
		ebounds[2,2] <- min(ebounds[2,])
	} else if (direction>0) {
		ebounds[2,2] <- min(ebounds[2,2],start[[9]])
		ebounds[2,1] <- min(ebounds[2,])
		ebounds[1,2] <- max(ebounds[1,])
	} else {
		ebounds[1,2] <- max(ebounds[1,])
		ebounds[2,2] <- min(ebounds[2,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (1) inner bound on fA
	ibounds <- getInnerBounds(direction,cbind(ebounds,c(start[[8]],start[[8]]),c(start[[9]],start[[9]])))
	ibounds <- ibounds[,1]
	if (ebounds[1,2]>start[[9]] || ebounds[2,2]<start[[9]]) {
		ibounds[[2]] <- min(ibounds[[2]],0.999)
	}

	# Get (1) initial outer bound on E0
	base_obounds <- initialBound1d(ebounds[,1],c("E0min","E0max"))

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[[1]]),ibounds[[2]])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,2,drop=FALSE]


	# Define parOuterBounds1
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],1)
		bounds <- base_obounds
		bnms <- c("EfA")
		btps <- c("min","max")
		for (j in 1:2) {
			i <- 1
			if (is.finite(abbounds[j,i]) && fpar[i]<1) {
				thisName <- paste0(bnms[i],btps[j])
				thresh <- (abbounds[j,i] - fpar[[i]]*start[[9]])/(1-fpar[i])
				nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)

				if (!is.null(nbounds)) { bounds <- nbounds }
				else {
					if (j==1) {
						bounds$lower <- bounds$upper
						bounds$lower$name <- thisName
						fpar[[i]] <- (abbounds[j,i]-bounds$lower$value)/(start[[9]]-bounds$lower$value)
					} else {
						bounds$upper <- bounds$lower
						bounds$upper$name <- thisName
						fpar[[i]] <- (abbounds[j,i]-bounds$upper$value)/(start[[9]]-bounds$upper$value)
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		sfact <- 1-sfact
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[9]],bsf$bounds)

		fpar <- c(bsf$sfpar[1:5],ebnd[[1]],start[[9]],start[[9]],start[[9]])
		fpar[7:8] <- ebnd[[1]]+(start[[9]]-ebnd[[1]])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		sfact <- 1-sfact
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- 1-sfact
		ebnd <- boundedOpt1d(mnv,start[[9]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		sfact <- (start[[9]]-ebnd[[1]])*sfact+ebnd[[1]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(start[[9]]-ebnd[[1]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }

			if (bname %in% c("EfAmin","EfAmax")) {
				bi <- which(c("EfAmin","EfAmax")==bname)
				derivs[[6]] <- derivs[[6]] + length(act)*ebnd[[3]]*(abbounds[bi,1]-start[[9]])/((1-bsf$sfpar[[6]])^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"II_2A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0 varies, EfA and Ef fixed, EfB varies
fitBraidScenario_II_2B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_II_2A(concs[,c(2,1)],
									  act,
									  swapModelVector(model),
									  weights,
									  swapParameterVector(start),
									  direction,
									  swapBoundMat(pbounds,model),
									  kweight)
	swapBraidFitObject(flipfit,"II_2B")
}

# E0 and Ef vary, EfB fixed at constant, EfA fixed equal to Ef
fitBraidScenario_II_3A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_II_3B(concs[,c(2,1)],
									  act,
									  swapModelVector(model),
									  weights,
									  swapParameterVector(start),
									  direction,
									  swapBoundMat(pbounds,model),
									  kweight)
	swapBraidFitObject(flipfit,"II_3A")
}
# E0 and Ef vary, EfA fixed at constant, EfB fixed equal to Ef
fitBraidScenario_II_3B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[2,1]<=start[[7]] || ebounds[1,2]>=start[[7]]) {
		direction <- setDirection(1,direction)
	}
	if (ebounds[1,1]>=start[[7]] || ebounds[2,2]<=start[[7]]) {
		direction <- setDirection(-1,direction)
	}
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[7]])
		ebounds[1,2] <- max(ebounds[1,2],start[[7]])
	} else if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[7]])
		ebounds[2,2] <- min(ebounds[2,2],start[[7]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,cbind(ebounds[,1],c(start[[7]],start[[7]]),
											  c(start[[8]],start[[8]]),ebounds[,2]))
	ibounds <- ibounds[,1]

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[[1]]),ibounds[[2]])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],1)
		fA <- fpar[[1]]
		if (fA>.Machine$double.eps) {
			bounds <- initialBound1d(start[[7]]*(1+1/fA)-rev(abbounds[,1])/fA,c("E0max","E0min"))
			if (fA<1) {
				for (j in 1:2) {
					thisName <- c("Efmin","EfMax")[[j]]
					thresh <- start[[7]]*(1-1/(1-fpar[[1]]))+abbounds[j,2]/(1-fpar[[1]])
					nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)

					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						if (j==1) {
							bounds$lower <- bounds$upper
							bounds$lower$name <- thisName
							fpar[[1]] <- (bounds$lower$value-abbounds[j,2])/(bounds$lower$value-start[[7]])
						} else {
							bounds$upper <- bounds$lower
							bounds$upper$name <- thisName
							fpar[[1]] <- (bounds$lower$value-abbounds[j,2])/(bounds$upper$value-start[[7]])
						}
					}
				}
				# bounds <- addBound1d(bounds,c(start[[7]]*(1-1/(1-fA))+abbounds[1,2]/(1-fA),1),"Efmin")
				# bounds <- addBound1d(bounds,c(start[[7]]*(1-1/(1-fA))+abbounds[2,2]/(1-fA),-1),"Efmax")
			}
		} else {
			bounds <- initialBound1d(start[[7]]*(1-1/(1-fA))+abbounds[,2]/(1-fA),c("Efmin","Efmax"))
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,abbounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		fA <- bsf$sfpar[[6]]
		sfact <- sfact-fA
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- sfact+fA
		ebnd <- boundedOpt1d(mnv,start[[7]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		e0 <- (1+fA)*start[[7]]-fA*ebnd[[1]]
		ef <- (1-fA)*ebnd[[1]]+fA*start[[7]]
		fpar <- c(bsf$sfpar[1:5],e0,start[[7]],ef,ef)
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,abbounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		fA <- bsf$sfpar[[6]]
		sfact <- sfact-fA
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- sfact+fA
		ebnd <- boundedOpt1d(mnv,start[[7]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)


		e0 <- (1+fA)*start[[7]]-fA*ebnd[[1]]
		ef <- (1-fA)*ebnd[[1]]+fA*start[[7]]

		sfact <- (ef-e0)*sfact+e0-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ef-e0))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }
			if (bname %in% c("E0min","E0max")) {
				bi <- which(c("E0min","E0max")==bname)
				derivs[[6]] <- derivs[[6]] - length(act)*ebnd[[3]]*(start[[7]]-abbounds[bi,1])/(fA^2)
			} else if (bname %in% c("Efmin","Efmax")) {
				bi <- which(c("Efmin","Efmax")==bname)
				derivs[[6]] <- derivs[[6]] + length(act)*ebnd[[3]]*(abbounds[bi,2]-start[[7]])/((1-fA)^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,abbounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"II_3B",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0 and EfA vary, EfB fixed at constant, Ef fixed equal to maximum
fitBraidScenario_II_4A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,FALSE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[2,1]<=start[[8]]) {
		direction <- setDirection(1,direction)
	}
	if (ebounds[1,1]>=start[[8]]) {
		direction <- setDirection(-1,direction)
	}
	checkModel1 <- checkModel2 <- TRUE
	if (direction>0) {
		if (ebounds[1,2]>=start[[8]]) { checkModel1 <- FALSE }
		if (ebounds[2,2]<=start[[8]]) { checkModel2 <- FALSE }
	} else if (direction<0) {
		if (ebounds[2,2]<=start[[8]]) { checkModel1 <- FALSE }
		if (ebounds[1,2]>=start[[8]]) { checkModel2 <- FALSE }
	}

	if (checkModel1) {
		model1 <- c(model[which(model<=5)],6,7)
		start1 <- start
		if ((start[[8]]>=start[[6]] && start[[7]]>start[[8]]) ||
			(start[[8]]<=start[[6]] && start[[7]]<start[[8]])) {
			start1[[7]] <- start[[8]]
		}
		start1[[9]] <- start1[[8]]
		bfit1 <- fitBraidScenario_II_2As(concs,act,model1,weights,start1,direction,pbounds,kweight=0)
	}
	if (checkModel2) {
		model2 <- c(model[which(model<=5)],6,9)
		start2 <- start
		if ((start[[8]]>=start[[6]] && start[[7]]<start[[8]]) ||
			(start[[8]]<=start[[6]] && start[[7]]>start[[8]])) {
			start2[[7]] <- start[[8]]
		}
		start2[[9]] <- start2[[7]]
		bfit2 <- fitBraidScenario_II_3B(concs,act,model2,weights,start2,direction,pbounds,kweight=0)
	}

	if (checkModel1 && checkModel2) {
		sse1 <- sum((bfit1$residuals*bfit1$weights)^2)
		sse2 <- sum((bfit2$residuals*bfit2$weights)^2)
		if (sse1<=sse2) { bfit <- bfit1 }
		else { bfit <- bfit2 }
	} else if (checkModel1) {
		bfit <- bfit1
	} else if (checkModel2) {
		bfit <- bfit2
	} else {
		stop("Unsatisfiable bounds.")
	}

	# bfit$scenario <- "II_4A"
	# bfit$model <- model
	# bfit$start <- start
	# bfit$direction <- direction
	# bfit$pbounds <- pbounds
	bfit
}
# E0 and EfB vary, EfA fixed at constant, Ef fixed equal to maximum
fitBraidScenario_II_4B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_II_4A(concs[,c(2,1)],
									  act,
									  swapModelVector(model),
									  weights,
									  swapParameterVector(start),
									  direction,
									  swapBoundMat(pbounds,model),
									  kweight)
	# swapBraidFitObject(flipfit,"II_4B")
	if (flipfit$scenario=="II_2A") {
		swapBraidFitObject(flipfit,"II_2B")
	} else if (flipfit$scenario=="II_3B") {
		swapBraidFitObject(flipfit,"II_3A")
	} else { stop("Invalid sub-scenario.") }
}

# Scenario II 6: E0 and Ef vary freely, EfA and EfB fixed equal to Ef
fitBraidScenario_II_6 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,FALSE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (ebounds[1,2]>ebounds[2,1]) { direction <- setDirection(1,direction) }
	if (ebounds[2,2]<ebounds[1,1]) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,])
		ebounds[1,2] <- max(ebounds[1,])
	} else if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,])
		ebounds[2,2] <- min(ebounds[2,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Get (2) outer bounds on E0 and Ef
	base_obounds <- getOuterBounds(direction,ebounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	# Add starting elision of fA and fB
	nstart <- bANDs$start
	nbounds <- bANDs$bounds

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(1,1)
		bounds <- base_obounds
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,bsf$bounds)
		fpar <- c(bsf$sfpar[1:5],ebnds[1],ebnds[2],ebnds[2],ebnds[2])
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		obounds <- bsf$bounds
		ebnds <- boundedOpt2d(mnv,obounds)
		ebnds[4:5] <- ebnds[4:5]*mean(weights^2)
		sfact <- (ebnds[2]-ebnds[1])*sfact+ebnds[1]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnds[2]-ebnds[1]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[model[model<=5]]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"II_6",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0 and EfB fixed, EfA and Ef vary
fitBraidScenario_II_7A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	if (start[[8]]>start[[6]]) { direction <- setDirection(1,direction) }
	else if (start[[8]]<start[[6]]) { direction <- setDirection(-1,direction) }
	else { stop("Scenario allows effect of only one drug.") }
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[6]])
		ebounds[1,2] <- max(c(ebounds[1,],start[[8]]))
		ebounds[2,1] <- min(ebounds[2,])
	} else if (direction<0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[6]])
		ebounds[2,2] <- min(c(ebounds[2,],start[[8]]))
		ebounds[1,1] <- max(ebounds[1,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }

	# Ensure finite effect bounds
	eminmax <- start[c(6,9)]
	eparscale <- abs(diff(eminmax))
	if (any(is.infinite(ebounds[1,]))) {
		loshot <- min(eminmax)-2*eparscale
		newEmin <- pmin(ebounds[2,],loshot)
		ebounds[1,is.infinite(ebounds[1,])] <- newEmin[is.infinite(ebounds[1,])]
	}
	if (any(is.infinite(ebounds[2,]))) {
		hishot <- max(eminmax)+2*eparscale
		newEmax <- pmax(ebounds[1,],hishot)
		ebounds[2,is.infinite(ebounds[2,])] <- newEmax[is.infinite(ebounds[2,])]
	}
	pbounds[,which(model>5)] <- ebounds

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	if (direction>0) {
		if (start[[9]]<ebounds[2,1]) {
			fpar <- (start[[7]]-ebounds[1,1])/(start[[9]]-ebounds[1,1])
		} else {
			fpar <- (start[[7]]-ebounds[1,1])/(ebounds[2,1]-ebounds[1,1])
		}
	} else {
		if (start[[9]]>ebounds[1,1]) {
			fpar <- (start[[7]]-start[[9]])/(ebounds[2,1]-start[[9]])
		} else {
			fpar <- (start[[7]]-ebounds[1,1])/(ebounds[2,1]-ebounds[1,1])
		}
	}
	nstart <- c(bANDs$start,fpar,start[[9]])
	nbounds <- cbind(bANDs$bounds,c(0,1),ebounds[,2])
	abbounds <- ebounds[,1]

	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		ef <- parv[[length(parv)]]
		fpar <- parv[[length(parv)-1]]
		if (direction>0) {
			if (ef<abbounds[[2]]) {
				efA <- abbounds[[1]]+(ef-abbounds[[1]])*fpar
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
			}
		} else {
			if (ef>abbounds[[1]]) {
				efA <- ef+(abbounds[[2]]-ef)*fpar
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
			}
		}
		fpar <- c(sfpar,start[[6]],efA,start[[8]],ef)
		return(fpar)
	}
	valderivfunc <- function(parv,concs,act,weights,start,model,abbounds,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		ef <- parv[[length(parv)]]
		fpar <- parv[[length(parv)-1]]
		if (direction>0) {
			if (ef<abbounds[[2]]) {
				efA <- abbounds[[1]]+(ef-abbounds[[1]])*fpar
				escale <- (ef-abbounds[[1]])
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
				escale <- (abbounds[[2]]-abbounds[[1]])
			}
		} else {
			if (ef>abbounds[[1]]) {
				efA <- ef+(abbounds[[2]]-ef)*fpar
				escale <- (abbounds[[2]]-ef)
			} else {
				efA <- abbounds[[1]]+(abbounds[[2]]-abbounds[[1]])*fpar
				escale <- (abbounds[[2]]-abbounds[[1]])
			}
		}
		fpar <- (c(efA,start[[8]])-start[[6]])/(ef-start[[6]])
		sfpar <- c(sfpar,fpar)

		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (ef-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		delE <- ef-start[[6]]
		derivs <- (2*(ef-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)

		derivs[[7]] <- -(start[[8]]-start[[6]])*derivs[[7]]/(delE^2)
		derivs[[7]] <- derivs[[7]] - (efA-start[[6]])*derivs[[6]]/(delE^2)
		derivs[[7]] <- derivs[[7]] + 2*sum(weights*weights*sfact*sfres$value)
		derivs[[6]] <- escale*derivs[[6]]/delE

		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,abbounds)

	# Run optim
	parscale <- rep(1,length(nstart))
	parscale[[length(nstart)]] <- eparscale
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds,parscale=parscale)

	braidFitObject(nls,"II_7A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0 and EfA fixed, EfB and Ef vary
fitBraidScenario_II_7B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_II_7A(concs[,c(2,1)],
									  act,
									  swapModelVector(model),
									  weights,
									  swapParameterVector(start),
									  direction,
									  swapBoundMat(pbounds,model),
									  kweight)
	swapBraidFitObject(flipfit,"II_7B")
}

# E0 and Ef fixed, EfA and EfB vary freely between them
fitBraidScenario_II_8 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,TRUE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction and inner bounds
	if (start[[9]]>start[[6]]) { direction <- setDirection(1,direction) }
	else if (start[[9]]<start[[6]]) { direction <- setDirection(-1,direction) }
	else { stop("Parameter 'start` specifies a constant surface.") }
	ebounds <- pbounds[,which(model>5)]
	if (direction>0) {
		ebounds[1,1:2] <- pmax(ebounds[1,1:2],start[[6]])
		ebounds[2,1:2] <- pmin(ebounds[2,1:2],start[[9]])
		ibounds <- (ebounds-start[[6]])/(start[[9]]-start[[6]])
	} else if (direction<0) {
		ebounds[2,1:2] <- pmin(ebounds[2,1:2],start[[6]])
		ebounds[1,1:2] <- pmax(ebounds[1,1:2],start[[9]])
		ibounds <- (ebounds[c(2,1),]-start[[6]])/(start[[9]]-start[[6]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[7:8]-start[[6]])/(start[[9]]-start[[6]])
	istart <- pmin(pmax(istart,ibounds[1,]),ibounds[2,])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)

	# Define par2fullpar
	par2fullpar <- function(parv,start,model) {
		sfpar <- preFillScaleFree(parv,start,model)
		sfpar <- c(sfpar,start[[6]],rep(start[[9]],3))
		sfpar[7:8] <- parv[which(model>5)]*(start[[9]]-start[[6]])+start[[6]]
		return(sfpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,kweight) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- parv[which(model>5)]
		sfpar <- c(sfpar,fpar)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (start[[9]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(start[[9]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(sfpar[1:4],sfpar[5]+2,1,1)
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,kweight)
	fpfunc <- function(p) par2fullpar(p,start,model)

	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"II_8",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0 fixed, EfA and Ef vary, EfB fixed equal to Ef
fitBraidScenario_II_9A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}
	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (any(ebounds[1,]>start[[6]])) { direction <- setDirection(1,direction) }
	if (any(ebounds[2,]<start[[6]])) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[6]])
		ebounds[1,2] <- max(ebounds[1,])
		ebounds[2,1] <- min(ebounds[2,])
	} else if (direction<0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[6]])
		ebounds[2,2] <- min(ebounds[2,])
		ebounds[1,1] <- max(ebounds[1,])
	} else {
		ebounds[1,1] <- max(ebounds[1,])
		ebounds[2,1] <- min(ebounds[2,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (1) inner bound on fA
	ibounds <- getInnerBounds(direction,cbind(c(start[[6]],start[[6]]),ebounds[,1],
											  c(start[[8]],start[[8]]),ebounds[,2]))
	ibounds <- ibounds[,1]
	if (ebounds[1,2]>start[[6]] || ebounds[2,2]<start[[6]]) {
		ibounds[[1]] <- max(ibounds[[1]],0.001)
	}

	# Get (1) initial outer bound on Ef
	base_obounds <- initialBound1d(ebounds[,2],c("Efmin","Efmax"))

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[[1]]),ibounds[[2]])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,1,drop=FALSE]


	# Define parOuterBounds1
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],1)
		bounds <- base_obounds
		bnms <- c("EfA")
		btps <- c("min","max")
		for (j in 1:2) {
			i <- 1
			if (is.finite(abbounds[j,i]) && fpar[i]>.Machine$double.eps) {
				thisName <- paste0(bnms[i],btps[j])
				thresh <- (abbounds[j,i] - (1-fpar[[i]])*start[[6]])/fpar[i]
				nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)

				if (!is.null(nbounds)) { bounds <- nbounds }
				else {
					if (j==1) {
						bounds$lower <- bounds$upper
						bounds$lower$name <- thisName
						fpar[[i]] <- (abbounds[j,i]-start[[6]])/(bounds$lower$value-start[[6]])
					} else {
						bounds$upper <- bounds$lower
						bounds$upper$name <- thisName
						fpar[[i]] <- (abbounds[j,i]-start[[6]])/(bounds$upper$value-start[[6]])
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)

		fpar <- c(bsf$sfpar[1:5],start[[6]],ebnd[[1]],ebnd[[1]],ebnd[[1]])
		fpar[7:8] <- start[[6]]+(ebnd[[1]]-start[[6]])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		sfact <- (ebnd[[1]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnd[[1]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }

			if (bname %in% c("EfAmin","EfAmax")) {
				bi <- which(c("EfAmin","EfAmax")==bname)
				derivs[[6]] <- derivs[[6]] - length(act)*ebnd[[3]]*(abbounds[bi,1]-start[[6]])/(bsf$sfpar[[6]]^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	braidFitObject(nls,"II_9A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0 fixed, EfB and Ef vary, EfA fixed equal to Ef
fitBraidScenario_II_9B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_II_9A(concs[,c(2,1)],
									  act,
									  swapModelVector(model),
									  weights,
									  swapParameterVector(start),
									  direction,
									  swapBoundMat(pbounds,model),
									  kweight)
	swapBraidFitObject(flipfit,"II_9B")
}

# E0 fixed, EfA and EfB vary freely Ef fixed equal to maximum
fitBraidScenario_II_10 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,TRUE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (any(ebounds[1,]>=start[[6]])) { direction <- setDirection(1,direction) }
	if (any(ebounds[2,]<=start[[6]])) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[1,] <- pmax(ebounds[1,],start[[6]])
		efbound <- c(max(ebounds[1,]),max(ebounds[2,]))
	} else if (direction<0) {
		ebounds[2,] <- pmin(ebounds[2,],start[[6]])
		efbound <- c(min(ebounds[1,]),min(ebounds[2,]))
	} else { efbound <- c(min(ebounds[1,]),max(ebounds[2,])) }
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (1) inner bound on the elision of fA and fB
	ibounds <- getInnerBounds(direction,cbind(c(start[[6]],start[[6]]),ebounds,efbound))
	if (all(ibounds[2,]<1)) { stop("Unsatisfiable bounds.") }
	else if (ibounds[2,1]<1) { ibounds <- as.numeric(ibounds[,1]) }
	else if (ibounds[2,2]<1) { ibounds <- rev(2-as.numeric(ibounds[,2])) }
	else { ibounds <- c(ibounds[1,1],2-ibounds[1,2]) }

	# Get (1) initial outer bound on Ef
	base_obounds <- initialBound1d(efbound,c("Efmin","Efmax"))

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	# Add starting elision of fA and fB
	nstart <- c(bANDs$start,pmin(pmax(1,ibounds[1]),ibounds[2]))
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],2-parv[length(parv)])
		fpar[fpar>1] <- 1
		bounds <- base_obounds
		bnms <- c("EfA","EfB")
		btps <- c("min","max")
		if (fpar[1]==1) { iv <- 1:2 } else { iv <- 2:1 }
		for (j in 1:2) {
			for (i in iv) {
				if (is.finite(abbounds[j,i]) && fpar[i]>.Machine$double.eps) {
					thisName <- paste0(bnms[i],btps[j])
					thresh <- (abbounds[j,i] - (1-fpar[i])*start[[6]])/fpar[[i]]
					nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)

					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						if (j==1) {
							bounds$lower <- bounds$upper
							bounds$lower$name <- thisName
							fpar[[i]] <- (abbounds[j,i]-start[[6]])/(bounds$lower$value-start[[6]])
						} else {
							bounds$upper <- bounds$lower
							bounds$upper$name <- thisName
							fpar[[i]] <- (abbounds[j,i]-start[[6]])/(bounds$upper$value-start[[6]])
						}
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)

		fpar <- c(bsf$sfpar[1:5],start[[6]],ebnd[[1]],ebnd[[1]],ebnd[[1]])
		fpar[7:8] <- start[[6]]+(ebnd[[1]]-start[[6]])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		sfact <- (ebnd[[1]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnd[[1]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }

			if (bname %in% c("EfAmin","EfAmax")) {
				bi <- which(c("EfAmin","EfAmax")==bname)
				derivs[[6]] <- derivs[[6]] - length(act)*ebnd[[3]]*(abbounds[bi,1]-start[[6]])/(bsf$sfpar[[6]]^2)
			} else if (bname %in% c("EfBmin","EfBmax")) {
				bi <- which(c("EfBmin","EfBmax")==bname)
				derivs[[7]] <- derivs[[7]] - length(act)*ebnd[[3]]*(abbounds[bi,2]-start[[6]])/(bsf$sfpar[[7]]^2)
			}
		}
		if (parv[length(parv)]>1) { derivs[6] <- -derivs[7] }
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	npar <- length(nstart)
	if (nbounds[1,npar]<0.5 || nbounds[2,npar]>1.5) {
		tnbounds <- nbounds
		tnbounds[,npar] <- c(max(nbounds[1,npar],0.5),min(nbounds[2,npar],1.5))
		tnls <- runBoundedOptim(vdfunc,fpfunc,nstart,tnbounds)
		nls <- runBoundedOptim(vdfunc,fpfunc,tnls$par,nbounds)
	} else {
		nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	}

	braidFitObject(nls,"II_10",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0, EfA and Ef vary, EfB fixed at constant
fitBraidScenario_III_1A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (ebounds[2,1]<=start[[8]] || ebounds[1,3]>=start[[8]]) {
		direction <- setDirection(1,direction)
	}
	if (ebounds[1,1]>=start[[8]] || ebounds[2,3]<=start[[8]]) {
		direction <- setDirection(-1,direction)
	}
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,1],start[[8]])
		ebounds[2,2] <- min(ebounds[2,2:3])
		ebounds[1,2] <- max(ebounds[1,1:2])
		ebounds[1,3] <- max(ebounds[1,3],start[[8]])
	} else if (direction>0) {
		ebounds[1,1] <- max(ebounds[1,1],start[[8]])
		ebounds[1,2] <- max(ebounds[1,2:3])
		ebounds[2,2] <- min(ebounds[2,1:2])
		ebounds[2,3] <- min(ebounds[2,3],start[[8]])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	if (ebounds[1,2]>start[[8]] || ebounds[2,2]<start[[8]]) {
		return(fitBraidScenario_III_1As(concs,act,model,weights,start,direction,pbounds,kweight=0))
	}

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,cbind(ebounds[,1:2],c(start[[8]],start[[8]]),ebounds[,3]))
	ibounds[1,2] <- max(ibounds[1,2],0.001)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[7:8]-start[[6]])/(start[[9]]-start[[6]])
	istart <- pmin(pmax(istart,ibounds[1,]),ibounds[2,])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- parv[(length(parv)-1):length(parv)]

		bounds <- emptyBound1d()
		btps <- c("min","max")
		if (fpar[[2]]>.Machine$double.eps) {
			for (j in 1:2) {
				thisName <- paste0("E0",btps[j])
				if (is.finite(abbounds[j,1])) {
					thresh <- start[[8]] + (start[[8]]-abbounds[j,1])/fpar[[2]]
					bounds <- addBound1d(bounds,c(thresh,-sign(1.5-j)),thisName)
				}
			}
		}
		if (fpar[[2]]<1) {
			for (j in 1:2) {
				dprecision <- abs(start[[9]]-start[[6]])/1e8
				thisName <- paste0("Ef",btps[j])
				if (is.finite(abbounds[j,3])) {
					thresh <- start[[8]] + (abbounds[j,3]-start[[8]])/(1-fpar[[2]])
					if (j==1 && abs(thresh-bounds$lower$value)<dprecision) {
						bounds$lower$name <- "lower"
					} else if (j==2 && abs(thresh-bounds$upper$value)<dprecision) {
						bounds$upper$name <- "upper"
					} else {
						bounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)
					}
				}
			}
		}
		if (fpar[[1]]!=fpar[[2]]) {
			for (j in 1:2) {
				thisName <- paste0("EfA",btps[j])
				if (is.finite(abbounds[j,2])) {
					thresh <- start[[8]] + (abbounds[j,2]-start[[8]])/(fpar[[1]]-fpar[[2]])
					nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)*sign(fpar[[1]]-fpar[[2]])),thisName)
					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						if ((j==1 && fpar[[1]]>fpar[[2]]) || (j==2 && fpar[[1]]<fpar[[2]])) {
							bounds$lower <- bounds$upper
							bounds$lower$name <- thisName
							fpar[[1]] <- (abbounds[j,2]-start[[8]])/(bounds$lower$value-start[[8]])
						} else {
							bounds$upper <- bounds$lower
							bounds$upper$name <- thisName
							fpar[[1]] <- (abbounds[j,2]-start[[8]])/(bounds$upper$value-start[[8]])
						}
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,abbounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		fB <- bsf$sfpar[[7]]
		sfact <- sfact-fB
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[8]],bsf$bounds)

		e0 <- (1+fB)*start[[8]]-fB*ebnd[[1]]
		ef <- (1-fB)*ebnd[[1]]+fB*start[[8]]
		efA <- e0+(ef-e0)*bsf$sfpar[[6]]
		fpar <- c(bsf$sfpar[1:5],e0,efA,start[[8]],ef)
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,abbounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		fA <- bsf$sfpar[[6]]
		fB <- bsf$sfpar[[7]]
		sfact <- sfact-fB
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- sfact+fB
		ebnd <- boundedOpt1d(mnv,start[[8]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		e0 <- (1+fB)*start[[8]]-fB*ebnd[[1]]
		ef <- (1-fB)*ebnd[[1]]+fB*start[[8]]
		efA <- e0+(ef-e0)*bsf$sfpar[[6]]

		sfact <- (ef-e0)*sfact+e0-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ef-e0))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }
			if (bname %in% c("E0min","E0max")) {
				bi <- which(c("E0min","E0max")==bname)
				derivs[[7]] <- derivs[[7]] - length(act)*ebnd[[3]]*(start[[8]]-abbounds[bi,1])/(fB^2)
			} else if (bname %in% c("Efmin","Efmax")) {
				bi <- which(c("Efmin","Efmax")==bname)
				derivs[[7]] <- derivs[[7]] + length(act)*ebnd[[3]]*(abbounds[bi,3]-start[[8]])/((1-fB)^2)
			} else if (bname %in% c("EfAmin","EfAmx")) {
				bi <- which(c("EfAmin","EfAmax")==bname)
				derivs[[6]] <- derivs[[6]] - length(act)*ebnd[[3]]*(abbounds[bi,2]-start[[8]])/((fA-fB)^2)
				derivs[[7]] <- derivs[[7]] + length(act)*ebnd[[3]]*(abbounds[bi,3]-start[[8]])/((fA-fB)^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,abbounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"III_1A",concs,act,model,weights,start,direction,pbounds,kweight)
}
fitBraidScenario_III_1As <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	ebounds <- pbounds[,which(model==7)]
	if (ebounds[[2]]<start[[8]]) { erng <- c(ebounds[[2]],start[[8]]) }
	else if (ebounds[[1]]>start[[8]]) { erng <- c(start[[8]],ebounds[[1]])}
	else {
		stop("Scenario III_1As is only used when the bounds on EfA do not include the fixed value for EfB.")
	}

	if (direction>=0) {
		start1 <- start
		start1[[6]] <- min(start1[[6]],erng[[1]])
		start1[[9]] <- max(start1[[9]],erng[[2]])
		bfit1 <- fitBraidScenario_III_1Ad(concs,act,model,weights,start1,direction=1,pbounds,kweight=0)
	}
	if (direction<=0) {
		start2 <- start
		start2[[6]] <- max(start2[[6]],erng[[2]])
		start2[[9]] <- min(start2[[9]],erng[[1]])
		bfit2 <- fitBraidScenario_III_1Ad(concs,act,model,weights,start2,direction=-1,pbounds,kweight=0)
	}

	if (direction==0) {
		sse1 <- sum((bfit1$residuals*bfit1$weights)^2)
		sse2 <- sum((bfit2$residuals*bfit2$weights)^2)
		if (sse1<=sse2) { bfit <- bfit1 }
		else { bfit <- bfit2 }
	} else if (direction>0) {
		bfit <- bfit1
	} else if (direction<0) {
		bfit <- bfit2
	}

	bfit$scenario <- "II_1"
	bfit$model <- model
	bfit$start <- start
	bfit$direction <- direction
	bfit$pbounds <- pbounds
	bfit
}
fitBraidScenario_III_1Ad <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Rectify direction
	ebounds <- pbounds[,which(model>5),drop=FALSE]
	if (direction>0) {
		ebounds[2,1] <- min(c(ebounds[2,1:2],start[[8]]))
		ebounds[2,2] <- min(ebounds[2,2:3])
		ebounds[1,2] <- max(ebounds[1,1:2])
		ebounds[1,3] <- max(c(ebounds[1,2:3],start[[8]]))
	} else if (direction<0) {
		ebounds[1,1] <- max(c(ebounds[1,1:2],start[[8]]))
		ebounds[1,2] <- max(ebounds[1,2:3])
		ebounds[2,2] <- min(ebounds[2,1:2])
		ebounds[2,3] <- min(c(ebounds[2,2:3],start[[8]]))
	} else {
		stop("Scenario III_1Ad is only valid when a direction is set.")
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }

	# Ensure finite effect bounds
	eminmax <- start[c(6,9)]
	eparscale <- abs(diff(eminmax))
	if (any(is.infinite(ebounds[1,]))) {
		loshot <- min(eminmax)-2*eparscale
		newEmin <- pmin(ebounds[2,],loshot)
		ebounds[1,is.infinite(ebounds[1,])] <- newEmin[is.infinite(ebounds[1,])]
	}
	if (any(is.infinite(ebounds[2,]))) {
		hishot <- max(eminmax)+2*eparscale
		newEmax <- pmax(ebounds[1,],hishot)
		ebounds[2,is.infinite(ebounds[2,])] <- newEmax[is.infinite(ebounds[2,])]
	}
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	ibounds <- getInnerBounds(direction,cbind(c(start[[6]],start[[6]]),ebounds[,2],
											  c(start[[8]],start[[8]]),c(start[[9]],start[[9]])))
	ibounds <- ibounds[,1]
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	if (ibounds[[1]]==ibounds[[2]]) { istart <- 1 }
	else { istart <- min(max((istart-ibounds[[1]])/(ibounds[[2]]-ibounds[[1]]),0),1) }
	nstart <- c(bANDs$start,istart,start[[6]],start[[9]])
	nbounds <- cbind(bANDs$bounds,c(0,1),ebounds[,c(1,3)])
	abbounds <- ebounds[,2,drop=FALSE]


	# Define parOuterBounds2
	# Define par2sfpar
	getBounds1AndScaleFree <- function(parv,start,model,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		epar <- parv[(length(parv)-1):length(parv)]
		fpar <- parv[[(length(parv)-2)]]
		ibounds <- getInnerBounds(direction,cbind(c(epar[[1]],epar[[1]]),abbounds,
												  c(start[[8]],start[[8]]),c(epar[[2]],epar[[2]])))
		ibounds <- ibounds[,1]
		fpar <- c(ibounds[[1]]+(ibounds[[2]]-ibounds[[1]])*fpar,
				  (start[[8]]-epar[[1]])/(epar[[2]]-epar[[1]]))
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=ibounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,abbounds) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		epar <- parv[(length(parv)-1):length(parv)]
		fpar <- (start[7:8]-epar[[1]])/(epar[[2]]-epar[[1]])
		fpar <- c(bsf$sfpar[1:5],epar[[1]],epar[[2]],start[[8]],epar[[2]])
		fpar[[7]] <- epar[[1]]+(epar[[2]]-epar[[1]])*bsf$sfpar[[6]]
		return(fpar)
	}
	valderivfunc <- function(parv,concs,act,weights,start,model,abbounds,kweight) {
		bsf <- getBounds1AndScaleFree(parv,start,model,abbounds)
		epar <- parv[(length(parv)-1):length(parv)]

		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		sfact <- (epar[[2]]-epar[[1]])*sfact+epar[[1]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		delE <- epar[[2]]-epar[[1]]
		derivs <- (2*(epar[[2]]-epar[[1]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)

		if (delE != 0) {
			ederivs <- c(2*sum(weights*weights*sfact*(1-sfres$value)),
						 2*sum(weights*weights*sfact*sfres$value))
			ederivs[[1]] <- ederivs[[1]]+derivs[[7]]*(start[[8]]-epar[[2]])/(delE^2)
			ederivs[[2]] <- ederivs[[2]]-derivs[[7]]*(start[[8]]-epar[[1]])/(delE^2)
		} else { ederivs <- c(0,0) }
		derivs[[6]] <- derivs[[6]]*(bsf$bounds[[2]]-bsf$bounds[[1]])

		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- c(derivs[c(model[model<=5],6)],ederivs)
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,abbounds)

	# Run optim
	parscale <- rep(1,length(nstart))
	parscale[(length(nstart)-1):length(nstart)] <- eparscale
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds,parscale=parscale)

	braidFitObject(nls,"III_1A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0, EfB and Ef vary, EfA fixed at constant
fitBraidScenario_III_1B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_III_1A(concs[,c(2,1)],
									   act,
									   swapModelVector(model),
									   weights,
									   swapParameterVector(start),
									   direction,
									   swapBoundMat(pbounds,model),
									   kweight)
	swapBraidFitObject(flipfit,"III_1B")
}

# Ef fixed, E0, EfA, and EfB all vary freely
fitBraidScenario_III_2 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,TRUE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}
	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (any(ebounds[1,]>=start[[9]])) { direction <- setDirection(-1,direction) }
	if (any(ebounds[2,]<=start[[9]])) { direction <- setDirection(1,direction) }
	if (direction>0) {
		ebounds[2,] <- pmin(ebounds[2,],start[[9]])
		ebounds[1,2:3] <- pmax(ebounds[1,2:3],ebounds[1,1])
		ebounds[2,1] <- min(ebounds[2,])
	} else if (direction<0) {
		ebounds[1,] <- pmax(ebounds[1,],start[[9]])
		ebounds[2,2:3] <- pmin(ebounds[2,2:3],ebounds[2,1])
		ebounds[1,1] <- max(ebounds[1,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,cbind(ebounds,c(start[[9]],start[[9]])))

	# Get (1) inital outer bound on Ef
	base_obounds <- initialBound1d(ebounds[,1],c("Efmin","Efmax"))

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[7:8]-start[[6]])/(start[[9]]-start[[6]])
	istart <- pmin(pmax(istart,ibounds[1,]),ibounds[2,])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,2:3]

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)-1],parv[length(parv)])
		bounds <- base_obounds
		bnms <- c("EfA","EfB")
		btps <- c("min","max")
		if (fpar[1]==1) { iv <- 1:2 } else { iv <- 2:1 }
		for (j in 1:2) {
			for (i in iv) {
				if (is.finite(abbounds[j,i]) && fpar[i]<1) {
					thisName <- paste0(bnms[i],btps[j])
					thresh <- (abbounds[j,i] - fpar[[i]]*start[[9]])/(1-fpar[[i]])
					nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)

					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						if (j==1) {
							bounds$lower <- bounds$upper
							bounds$lower$name <- thisName
							fpar[[i]] <- (abbounds[j,i]-bounds$lower$value)/(start[[9]]-bounds$lower$value)
						} else {
							bounds$upper <- bounds$lower
							bounds$upper$name <- thisName
							fpar[[i]] <- (abbounds[j,i]-bounds$upper$value)/(start[[9]]-bounds$upper$value)
						}
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		sfact <- 1-sfact
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[9]],bsf$bounds)

		fpar <- c(bsf$sfpar[1:5],ebnd[[1]],start[[9]],start[[9]],start[[9]])
		fpar[7:8] <- ebnd[[1]]+(start[[9]]-ebnd[[1]])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		sfact <- 1-sfact
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		sfact <- 1-sfact
		ebnd <- boundedOpt1d(mnv,start[[9]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		sfact <- (start[[9]]-ebnd[[1]])*sfact+ebnd[[1]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(start[[9]]-ebnd[[1]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }

			if (bname %in% c("EfAmin","EfAmax")) {
				bi <- which(c("EfAmin","EfAmax")==bname)
				derivs[[6]] <- derivs[[6]] + length(act)*ebnd[[3]]*(abbounds[bi,1]-start[[9]])/((1-bsf$sfpar[[6]])^2)
			} else if (bname %in% c("EfBmin","EfBmax")) {
				bi <- which(c("EfBmin","EfBmax")==bname)
				derivs[[7]] <- derivs[[7]] + length(act)*ebnd[[3]]*(abbounds[bi,2]-start[[9]])/((1-bsf$sfpar[[7]])^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	npar <- length(nstart)
	if (nbounds[1,npar-1]<0.5 || nbounds[1,npar]<0.5) {
		tnbounds <- nbounds
		tnbounds[,(npar-1):npar] <- pmax(nbounds[,(npar-1):npar],0.5)
		tnls <- runBoundedOptim(vdfunc,fpfunc,nstart,tnbounds)
		nls <- runBoundedOptim(vdfunc,fpfunc,tnls$par,nbounds)
	} else {
		nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	}

	braidFitObject(nls,"III_2",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0, EfA, and Ef all vary, EfB fixed equal to Ef
fitBraidScenario_III_3A <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,FALSE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (any(ebounds[1,2:3]>ebounds[2,1])) { direction <- setDirection(1,direction) }
	if (any(ebounds[2,2:3]<ebounds[1,1])) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,])
		ebounds[1,2] <- max(ebounds[1,1:2])
		ebounds[2,2] <- min(ebounds[2,2:3])
		ebounds[1,3] <- max(ebounds[1,])
	} else if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,])
		ebounds[2,2] <- min(ebounds[2,1:2])
		ebounds[1,2] <- max(ebounds[1,2:3])
		ebounds[2,3] <- min(ebounds[3,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,cbind(ebounds[,1:2],ebounds[,3],ebounds[,3]))
	ibounds <- ibounds[,1]

	# Get (2) outer bounds on E0 and Ef
	base_obounds <- getOuterBounds(direction,ebounds[,c(1,3)])

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[[7]]-start[[6]])/(start[[9]]-start[[6]])
	istart <- min(max(istart,ibounds[[1]]),ibounds[[2]])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,2,drop=FALSE]

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],1)
		bounds <- base_obounds
		bnms <- c("EfA")
		btps <- c("min","max")
		if (fpar[1]==1) { iv <- 1:2 } else { iv <- 2:1 }
		for (j in 1:2) {
			i <- 1
			if (is.finite(abbounds[j,i])) {
				thisName <- paste0(bnms[i],btps[j])
				nbounds <- addBound2d(bounds,c(1-fpar[i],fpar[i],abbounds[j,i])*sign(1.5-j),thisName)
				if (!is.null(nbounds)) { bounds <- nbounds }
				else {
					fV <- (abbounds[j,i]-bounds$verts[,1])/(bounds$verts[,2]-bounds$verts[,1])
					if (fV[1]>=fpar[i]) { fpar[i] <- min(fV) } else { fpar[i] <- max(fV) }
					bounds <- addBound2d(bounds,c(1-fpar[i],fpar[i],abbounds[j,i])*sign(1.5-j),thisName)
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,bsf$bounds)
		fpar <- c(bsf$sfpar[1:5],ebnds[1],ebnds[2],ebnds[2],ebnds[2])
		fpar[7:8] <- ebnds[1]+(ebnds[2]-ebnds[1])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		obounds <- bsf$bounds
		ebnds <- boundedOpt2d(mnv,obounds)
		ebnds[4:5] <- ebnds[4:5]*mean(weights^2)
		sfact <- (ebnds[2]-ebnds[1])*sfact+ebnds[1]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnds[2]-ebnds[1]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnds[3]>0 && ebnds[1]!=ebnds[2]) {
			if ((2*ebnds[3])%%2==0) { eind <- ebnds[3] }
			else {
				eind <- c(floor(ebnds[3]),ceiling(ebnds[3]))
				if (eind[1]==0) { eind[1] <- nrow(obounds$edges) }
			}
			bnames <- row.names(obounds$edges)[eind]
			if (any(c("EfAmin","EfAmax")%in%bnames)) {
				bi <- which(c("EfAmin","EfAmax")%in%bnames)[1]
				ediv <- c(abbounds[bi,1]-ebnds[2],ebnds[1]-abbounds[bi,1])/((ebnds[1]-ebnds[2])^2)
				ediv <- ediv/sum(ediv^2)
				derivs[6] <- derivs[6]+(ediv[1]*ebnds[4]+ediv[2]*ebnds[5])*length(act)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)

	braidFitObject(nls,"III_3A",concs,act,model,weights,start,direction,pbounds,kweight)
}
# E0, EfA, and Ef all vary, EfB fixed equal to Ef
fitBraidScenario_III_3B <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	flipfit <- fitBraidScenario_III_3A(concs[,c(2,1)],
									   act,
									   swapModelVector(model),
									   weights,
									   swapParameterVector(start),
									   direction,
									   swapBoundMat(pbounds,model),
									   kweight)
	swapBraidFitObject(flipfit,"III_3B")
}

# Scenario III 4: E0, EfA, and EfB vary freely, Ef fixed equal to maximum
fitBraidScenario_III_4 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,TRUE,FALSE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (any(ebounds[1,2:3]>ebounds[2,1])) { direction <- setDirection(1,direction) }
	if (any(ebounds[2,2:3]<ebounds[1,1])) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,])
		ebounds[1,2:3] <- pmax(ebounds[1,2:3],ebounds[1,1])
		efbound <- c(max(ebounds[1,2:3]),max(ebounds[2,2:3]))
	} else if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,])
		ebounds[2,2:3] <- pmin(ebounds[2,2:3],ebounds[2,1])
		efbound <- c(min(ebounds[1,2:3]),min(ebounds[2,2:3]))
	} else { efbound <- c(min(ebounds[1,2:3]),max(ebounds[2,2:3])) }
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (1) inner bound on the elision of fA and fB
	ibounds <- getInnerBounds(direction,cbind(ebounds,efbound))
	if (all(ibounds[2,]<1)) { stop("Unsatisfiable bounds.") }
	else if (ibounds[2,1]<1) { ibounds <- as.numeric(ibounds[,1]) }
	else if (ibounds[2,2]<1) { ibounds <- rev(2-as.numeric(ibounds[,2])) }
	else { ibounds <- c(ibounds[1,1],2-ibounds[1,2]) }

	# Get (2) outer bounds on E0 and Ef
	base_obounds <- getOuterBounds(direction,cbind(ebounds[,1],efbound))

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	# Add starting elision of fA and fB
	nstart <- c(bANDs$start,pmin(pmax(1,ibounds[1]),ibounds[2]))
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,2:3]

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)],2-parv[length(parv)])
		fpar[fpar>1] <- 1
		bounds <- base_obounds
		bnms <- c("EfA","EfB")
		btps <- c("min","max")
		if (fpar[1]==1) { iv <- 1:2 } else { iv <- 2:1 }
		for (j in 1:2) {
			for (i in iv) {
				if (is.finite(abbounds[j,i])) {
					thisName <- paste0(bnms[i],btps[j])
					nbounds <- addBound2d(bounds,c(1-fpar[i],fpar[i],abbounds[j,i])*sign(1.5-j),thisName)
					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						fV <- (abbounds[j,i]-bounds$verts[,1])/(bounds$verts[,2]-bounds$verts[,1])
						if (fV[1]>=fpar[i]) { fpar[i] <- min(fV) } else { fpar[i] <- max(fV) }
						bounds <- addBound2d(bounds,c(1-fpar[i],fpar[i],abbounds[j,i])*sign(1.5-j),thisName)
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,bsf$bounds)
		fpar <- c(bsf$sfpar[1:5],ebnds[1],ebnds[2],ebnds[2],ebnds[2])
		fpar[7:8] <- ebnds[1]+(ebnds[2]-ebnds[1])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		obounds <- bsf$bounds
		ebnds <- boundedOpt2d(mnv,obounds)
		ebnds[4:5] <- ebnds[4:5]*mean(weights^2)
		sfact <- (ebnds[2]-ebnds[1])*sfact+ebnds[1]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnds[2]-ebnds[1]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnds[3]>0 && ebnds[1]!=ebnds[2]) {
			if ((2*ebnds[3])%%2==0) { eind <- ebnds[3] }
			else {
				eind <- c(floor(ebnds[3]),ceiling(ebnds[3]))
				if (eind[1]==0) { eind[1] <- nrow(obounds$edges) }
			}
			bnames <- row.names(obounds$edges)[eind]
			if (any(c("EfAmin","EfAmax")%in%bnames)) {
				bi <- which(c("EfAmin","EfAmax")%in%bnames)[1]
				ediv <- c(abbounds[bi,1]-ebnds[2],ebnds[1]-abbounds[bi,1])/((ebnds[1]-ebnds[2])^2)
				ediv <- ediv/sum(ediv^2)
				derivs[6] <- derivs[6]+(ediv[1]*ebnds[4]+ediv[2]*ebnds[5])*length(act)
			}
			if (any(c("EfBmin","EfBmax")%in%bnames)) {
				bi <- which(c("EfBmin","EfBmax")%in%bnames)[1]
				ediv <- c(abbounds[bi,2]-ebnds[2],ebnds[1]-abbounds[bi,2])/((ebnds[1]-ebnds[2])^2)
				ediv <- ediv/sum(ediv^2)
				derivs[7] <- derivs[7]+(ediv[1]*ebnds[4]+ediv[2]*ebnds[5])*length(act)
			}
		}
		if (parv[length(parv)]>1) { derivs[6] <- -derivs[7] }
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	npar <- length(nstart)
	if (nbounds[1,npar]<0.5 || nbounds[2,npar]>1.5) {
		tnbounds <- nbounds
		tnbounds[,npar] <- c(max(nbounds[1,npar],0.5),min(nbounds[2,npar],1.5))
		tnls <- runBoundedOptim(vdfunc,fpfunc,nstart,tnbounds)
		nls <- runBoundedOptim(vdfunc,fpfunc,tnls$par,nbounds)
	} else {
		nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	}

	braidFitObject(nls,"III_4",concs,act,model,weights,start,direction,pbounds,kweight)
}

# E0 fixed, EfA, EfB, and Ef all vary freely
fitBraidScenario_III_5 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(FALSE,TRUE,TRUE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (any(ebounds[1,]>=start[[6]])) { direction <- setDirection(1,direction) }
	if (any(ebounds[2,]<=start[[6]])) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[1,] <- pmax(ebounds[1,],start[[6]])
		ebounds[2,1:2] <- pmin(ebounds[2,1:2],ebounds[2,3])
		ebounds[1,3] <- max(ebounds[1,])
	} else if (direction<0) {
		ebounds[2,] <- pmin(ebounds[2,],start[[6]])
		ebounds[1,1:2] <- pmax(ebounds[1,1:2],ebounds[1,3])
		ebounds[2,3] <- min(ebounds[2,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,cbind(c(start[[6]],start[[6]]),ebounds))

	# Get (1) inital outer bound on Ef
	base_obounds <- initialBound1d(ebounds[,3],c("Efmin","Efmax"))

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[7:8]-start[[6]])/(start[[9]]-start[[6]])
	istart <- pmin(pmax(istart,ibounds[1,]),ibounds[2,])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,1:2]

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- c(parv[length(parv)-1],parv[length(parv)])
		bounds <- base_obounds
		bnms <- c("EfA","EfB")
		btps <- c("min","max")
		if (fpar[1]==1) { iv <- 1:2 } else { iv <- 2:1 }
		for (j in 1:2) {
			for (i in iv) {
				if (is.finite(abbounds[j,i]) && fpar[i]>.Machine$double.eps) {
					thisName <- paste0(bnms[i],btps[j])
					thresh <- (abbounds[j,i] - (1-fpar[i])*start[[6]])/fpar[[i]]
					nbounds <- addBound1d(bounds,c(thresh,sign(1.5-j)),thisName)

					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						if (j==1) {
							bounds$lower <- bounds$upper
							bounds$lower$name <- thisName
							fpar[[i]] <- (abbounds[j,i]-start[[6]])/(bounds$lower$value-start[[6]])
						} else {
							bounds$upper <- bounds$lower
							bounds$upper$name <- thisName
							fpar[[i]] <- (abbounds[j,i]-start[[6]])/(bounds$upper$value-start[[6]])
						}
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)

		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)

		fpar <- c(bsf$sfpar[1:5],start[[6]],ebnd[[1]],ebnd[[1]],ebnd[[1]])
		fpar[7:8] <- start[[6]]+(ebnd[[1]]-start[[6]])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[6]],bsf$bounds)
		ebnd[[3]] <- ebnd[[3]]*mean(weights^2)

		sfact <- (ebnd[[1]]-start[[6]])*sfact+start[[6]]-act
		ovalue <- sum((weights*sfact)^2)


		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnd[[1]]-start[[6]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnd[[2]]!=0) {
			if (ebnd[[2]]>0) { bname <- bsf$bounds$lower$name }
			else { bname <- bsf$bounds$upper$name }

			if (bname %in% c("EfAmin","EfAmax")) {
				bi <- which(c("EfAmin","EfAmax")==bname)
				derivs[[6]] <- derivs[[6]] - length(act)*ebnd[[3]]*(abbounds[bi,1]-start[[6]])/(bsf$sfpar[[6]]^2)
			} else if (bname %in% c("EfBmin","EfBmax")) {
				bi <- which(c("EfBmin","EfBmax")==bname)
				derivs[[7]] <- derivs[[7]] - length(act)*ebnd[[3]]*(abbounds[bi,2]-start[[6]])/(bsf$sfpar[[7]]^2)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	npar <- length(nstart)
	if (nbounds[1,npar-1]<0.5 || nbounds[1,npar]<0.5) {
		tnbounds <- nbounds
		tnbounds[,(npar-1):npar] <- pmax(nbounds[,(npar-1):npar],0.5)
		tnls <- runBoundedOptim(vdfunc,fpfunc,nstart,tnbounds)
		nls <- runBoundedOptim(vdfunc,fpfunc,tnls$par,nbounds)
	} else {
		nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	}

	braidFitObject(nls,"III_5",concs,act,model,weights,start,direction,pbounds,kweight)
}

# Scenario IV 1: All values vary freely
fitBraidScenario_IV_1 <- function(concs,act,model,weights,start,direction,pbounds,kweight=0) {
	# Check model values
	if (!checkModelValues(model,c(TRUE,TRUE,TRUE,TRUE))) {
		stop("'model' parameter does not match scenario.")
	}

	# Rectify direction
	ebounds <- pbounds[,which(model>5)]
	if (any(ebounds[1,2:4]>ebounds[2,1])) { direction <- setDirection(1,direction) }
	if (any(ebounds[2,2:4]<ebounds[1,1])) { direction <- setDirection(-1,direction) }
	if (direction>0) {
		ebounds[2,1] <- min(ebounds[2,])
		ebounds[1,2:3] <- pmax(ebounds[1,2:3],ebounds[1,1])
		ebounds[2,2:3] <- pmin(ebounds[2,2:3],ebounds[2,4])
		ebounds[1,4] <- max(ebounds[1,])
	} else if (direction<0) {
		ebounds[1,1] <- max(ebounds[1,])
		ebounds[2,2:3] <- pmin(ebounds[2,2:3],ebounds[2,1])
		ebounds[1,2:3] <- pmax(ebounds[1,2:3],ebounds[1,4])
		ebounds[2,4] <- min(ebounds[2,])
	}
	if (any(ebounds[1,]>ebounds[2,])) { stop("Unsatisfiable bounds.") }
	pbounds[,which(model>5)] <- ebounds

	# Rectify start
	start <- rectifyStart(model,start,pbounds)

	# Rectify (2) inner bounds on fA and fB
	ibounds <- getInnerBounds(direction,ebounds)

	# Get (2) outer bounds on E0 and Ef
	base_obounds <- getOuterBounds(direction,ebounds[,c(1,4)])

	# Specify nbounds
	# Specify nstart
	bANDs <- basicPotencyValues(model,start,pbounds)
	istart <- (start[7:8]-start[[6]])/(start[[9]]-start[[6]])
	istart <- pmin(pmax(istart,ibounds[1,]),ibounds[2,])
	nstart <- c(bANDs$start,istart)
	nbounds <- cbind(bANDs$bounds,ibounds)
	abbounds <- ebounds[,2:3]

	# Define parOuterBounds2
	# Define par2sfpar
	getBounds2AndScaleFree <- function(parv,start,model,base_obounds,abbounds) {
		sfpar <- preFillScaleFree(parv,start,model)
		fpar <- parv[(length(parv)-1):length(parv)]
		bounds <- base_obounds
		bnms <- c("EfA","EfB")
		btps <- c("min","max")
		if (fpar[1]==1) { iv <- 1:2 } else { iv <- 2:1 }
		for (j in 1:2) {
			for (i in iv) {
				if (is.finite(abbounds[j,i])) {
					thisName <- paste0(bnms[i],btps[j])
					nbounds <- addBound2d(bounds,c(1-fpar[i],fpar[i],abbounds[j,i])*sign(1.5-j),thisName)
					if (!is.null(nbounds)) { bounds <- nbounds }
					else {
						fV <- (abbounds[j,i]-bounds$verts[,1])/(bounds$verts[,2]-bounds$verts[,1])
						if (fV[1]>=fpar[i]) { fpar[i] <- min(fV) } else { fpar[i] <- max(fV) }
						bounds <- addBound2d(bounds,c(1-fpar[i],fpar[i],abbounds[j,i])*sign(1.5-j),thisName)
					}
				}
			}
		}
		sfpar <- c(sfpar,fpar)
		return(list(sfpar=sfpar,bounds=bounds))
	}
	# Define par2fullpar
	par2fullpar <- function(parv,concs,act,weights,start,model,base_obounds,abbounds) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfact <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=FALSE)
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,bsf$bounds)
		fpar <- c(bsf$sfpar[1:5],ebnds[1],ebnds[2],ebnds[2],ebnds[2])
		fpar[7:8] <- ebnds[1]+(ebnds[2]-ebnds[1])*bsf$sfpar[6:7]
		return(fpar)
	}
	# Define valderivfunc
	valderivfunc <- function(parv,concs,act,weights,start,model,base_obounds,abbounds,kweight) {
		bsf <- getBounds2AndScaleFree(parv,start,model,base_obounds,abbounds)
		wt2 <- (weights^2)/mean(weights^2)
		sfres <- evalBraidModel_sf(concs[,1],concs[,2],bsf$sfpar,calcderivs=TRUE)
		sfact <- sfres$value
		mnv <- c(mean(wt2*sfact),mean(wt2*sfact^2),mean(wt2*act),mean(wt2*sfact*act))
		obounds <- bsf$bounds
		ebnds <- boundedOpt2d(mnv,obounds)
		ebnds[4:5] <- ebnds[4:5]*mean(weights^2)
		sfact <- (ebnds[2]-ebnds[1])*sfact+ebnds[1]-act
		ovalue <- sum((weights*sfact)^2)

		if (5 %in% model) { epsilon <- parv[which(model==5)] }
		else { epsilon <- 0 }
		ovalue <- ovalue+kweight*(epsilon^2)

		derivs <- (2*(ebnds[2]-ebnds[1]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)
		derivs <- derivs*c(bsf$sfpar[1:4],bsf$sfpar[5]+2,1,1)
		if (ebnds[3]>0 && ebnds[1]!=ebnds[2]) {
			if ((2*ebnds[3])%%2==0) { eind <- ebnds[3] }
			else {
				eind <- c(floor(ebnds[3]),ceiling(ebnds[3]))
				if (eind[1]==0) { eind[1] <- nrow(obounds$edges) }
			}
			bnames <- row.names(obounds$edges)[eind]
			if (any(c("EfAmin","EfAmax")%in%bnames)) {
				bi <- which(c("EfAmin","EfAmax")%in%bnames)[1]
				ediv <- c(abbounds[bi,1]-ebnds[2],ebnds[1]-abbounds[bi,1])/((ebnds[1]-ebnds[2])^2)
				ediv <- ediv/sum(ediv^2)
				derivs[6] <- derivs[6]+(ediv[1]*ebnds[4]+ediv[2]*ebnds[5])*length(act)
			}
			if (any(c("EfBmin","EfBmax")%in%bnames)) {
				bi <- which(c("EfBmin","EfBmax")%in%bnames)[1]
				ediv <- c(abbounds[bi,2]-ebnds[2],ebnds[1]-abbounds[bi,2])/((ebnds[1]-ebnds[2])^2)
				ediv <- ediv/sum(ediv^2)
				derivs[7] <- derivs[7]+(ediv[1]*ebnds[4]+ediv[2]*ebnds[5])*length(act)
			}
		}
		derivs[5] <- derivs[5]+2*kweight*epsilon
		derivs <- derivs[c(model[model<=5],6,7)]
		return(list(value=ovalue,derivatives=derivs))
	}

	vdfunc <- function(p) valderivfunc(p,concs,act,weights,start,model,base_obounds,abbounds,kweight)
	fpfunc <- function(p) par2fullpar(p,concs,act,weights,start,model,base_obounds,abbounds)

	# Run optim
	npar <- length(nstart)
	if (nbounds[1,npar-1]<0.5 || nbounds[1,npar]<0.5) {
		tnbounds <- nbounds
		tnbounds[,(npar-1):npar] <- pmax(nbounds[,(npar-1):npar],0.5)
		tnls <- runBoundedOptim(vdfunc,fpfunc,nstart,tnbounds)
		nls <- runBoundedOptim(vdfunc,fpfunc,tnls$par,nbounds)
	} else {
		nls <- runBoundedOptim(vdfunc,fpfunc,nstart,nbounds)
	}

	braidFitObject(nls,"IV_1",concs,act,model,weights,start,direction,pbounds,kweight)
}

swapBraidFitObject <- function(bfit,newScenario=NULL) {
	coefficients = swapParameterVector(unname(bfit$coefficients))
	names(coefficients) <- names(bfit$coefficients)
	if (is.null(newScenario)) { newScenario <- bfit$scenario }
	structure(
		list(
			concs = bfit$concs[,c(2,1)],
			act = bfit$act,
			weights = bfit$weights,
			coefficients = coefficients,
			par = bfit$par,
			fitted.values = bfit$fitted,
			residuals = bfit$residuals,
			scenario = newScenario,
			model = swapModelVector(bfit$model),
			start = swapParameterVector(bfit$start),
			direction = bfit$direction,
			pbounds = swapBoundMat(bfit$pbounds,bfit$model),
			kweight = bfit$kweight
		),
		class="braidrm"
	)
}
swapParameterVector <- function(bpar) {
	bpar[braidSwapInds()]
}
swapModelVector <- function(model) {
	modelvec <- rep(FALSE,9)
	modelvec[model] <- TRUE
	modelvec <- modelvec[braidSwapInds()]
	which(modelvec)
}
swapBoundMat <- function(pbounds,model) {
	boundmat <- array(NA,dim=c(2,9))
	boundmat[,model] <- pbounds
	boundmat <- boundmat[,braidSwapInds()]
	boundmat[,swapModelVector(model)]
}
braidSwapInds <- function() { c(2,1,4,3,5,6,8,7,9) }

braidFitObject <- function(nls,scenario,concs,act,model,weights,start,direction,pbounds,kweight) {
	fullpar <- nls$fullpar
	names(fullpar) <- c("IDMA","IDMB","na","nb","kappa","E0","EfA","EfB","Ef")
	fitted <- evalBraidModel(concs[,1],concs[,2],fullpar)
	structure(
		list(
			concs = concs,
			act = act,
			weights = weights,
			coefficients = fullpar,
			par = nls$par,
			fitted.values = fitted,
			residuals = act-fitted,
			scenario = scenario,
			model = model,
			start = start,
			direction = direction,
			pbounds = pbounds,
			kweight = kweight
		),
		class="braidrm"
	)
}

checkModelValues <- function(model,checks) {
	pinds <- 6:9
	return(!((any(checks) && !all(pinds[checks]%in%model)) ||
			 	(any(!checks) && any(pinds[!checks]%in%model))))
}
setDirection <- function(newdirection,direction) {
	if (newdirection*direction<0) { stop("Invalid effect signs.") }
	return(newdirection)
}
rectifyStart <- function(model,start,pbounds) {
	newstart <- start
	newstart[model] <- pmin(pmax(newstart[model],pbounds[1,]),pbounds[2,])
	return(newstart)
}
getOuterBounds <- function(direction,spbounds) {
	base_obounds <- emptyBound2d()
	if (is.finite(spbounds[1,1])) { base_obounds <- addBound2d(base_obounds,c(1,0,spbounds[1,1]),"E0min") }
	if (is.finite(spbounds[2,1])) { base_obounds <- addBound2d(base_obounds,c(-1,0,-spbounds[2,1]),"E0max") }
	if (is.finite(spbounds[1,2])) { base_obounds <- addBound2d(base_obounds,c(0,1,spbounds[1,2]),"Efmin") }
	if (is.finite(spbounds[2,2])) { base_obounds <- addBound2d(base_obounds,c(0,-1,-spbounds[2,2]),"Efmax") }
	if (direction!=0) { base_obounds <- addBound2d(base_obounds,c(-direction,direction,0),"Direction") }
	if (is.null(base_obounds)) { stop("Bounds specified by 'bounds' and 'direction' parameters cannot be satisified.") }
	return(base_obounds)
}
getInnerBounds <- function(direction,spb) {
	bib1 <- bib2 <- cbind(c(0,1),c(0,1))
	if (direction>=0) {
		if (is.finite(spb[1,1])&&is.finite(spb[1,4])&&is.finite(spb[2,2])&&(spb[2,2]<spb[1,4])) {
			bib1[2,1] <- min(bib1[2,1],(spb[2,2]-spb[1,1])/(spb[1,4]-spb[1,1])) }
		if (is.finite(spb[1,1])&&is.finite(spb[1,4])&&is.finite(spb[2,3])&&(spb[2,3]<spb[1,4])) {
			bib1[2,2] <- min(bib1[2,2],(spb[2,3]-spb[1,1])/(spb[1,4]-spb[1,1])) }
		if (is.finite(spb[2,1])&&is.finite(spb[2,4])&&is.finite(spb[1,2])&&(spb[1,2]>spb[2,1])) {
			bib1[1,1] <- max(bib1[1,1],(spb[1,2]-spb[2,1])/(spb[2,4]-spb[2,1])) }
		if (is.finite(spb[2,1])&&is.finite(spb[2,4])&&is.finite(spb[1,3])&&(spb[1,3]>spb[2,1])) {
			bib1[1,2] <- max(bib1[1,2],(spb[1,3]-spb[2,1])/(spb[2,4]-spb[2,1])) }
	}
	if (direction<=0) {
		if (is.finite(spb[1,1])&&is.finite(spb[1,4])&&is.finite(spb[2,2])&&(spb[2,2]<spb[1,1])) {
			bib2[1,1] <- max(bib2[1,1],(spb[2,2]-spb[1,1])/(spb[1,4]-spb[1,1])) }
		if (is.finite(spb[1,1])&&is.finite(spb[1,4])&&is.finite(spb[2,3])&&(spb[2,3]<spb[1,1])) {
			bib2[1,2] <- max(bib2[1,2],(spb[2,3]-spb[1,1])/(spb[1,4]-spb[1,1])) }
		if (is.finite(spb[2,1])&&is.finite(spb[2,4])&&is.finite(spb[1,2])&&(spb[1,2]>spb[2,4])) {
			bib2[2,1] <- min(bib2[2,1],(spb[1,2]-spb[2,1])/(spb[2,4]-spb[2,1])) }
		if (is.finite(spb[2,1])&&is.finite(spb[2,4])&&is.finite(spb[1,3])&&(spb[1,3]>spb[2,4])) {
			bib2[2,2] <- min(bib2[2,2],(spb[1,3]-spb[2,1])/(spb[2,4]-spb[2,1])) }
	}
	if (direction==0) {
		base_ibounds <- rbind(pmin(bib1[1,],bib2[1,]),pmax(bib1[2,],bib2[2,]))
	} else if (direction>0) { base_ibounds <- bib1 }
	else { base_ibounds <- bib2 }
	if (any(base_ibounds[1,]>1)||any(base_ibounds[2,]<0)) { stop("Effect bounds cannot be satisfied.") }
	return(base_ibounds)
}
basicPotencyValues <- function(model,start,pbounds) {
	nbounds <- pbounds[,which(model<=5)]
	nstart <- start[which(model<=5)]
	if (any(model<5)) {
		nstart[which(model<5)] <- log(nstart[which(model<5)])
		nbounds[,which(model<5)] <- log(nbounds[,which(model<5)])
	}
	if (5%in%model) {
		nstart[which(model==5)] <- log((nstart[which(model==5)]+2)/2)
		nbounds[,which(model==5)] <- log((nbounds[,which(model==5)]+2)/2)
	}
	return(list(start=nstart,bounds=nbounds))
}
preFillScaleFree <- function(parv,start,model) {
	sfpar <- c(log(start[1:4]),log((start[5]+2)/2))
	sfpar[model[model<=5]] <- parv[which(model<=5)]
	sfpar <- c(exp(sfpar[1:4]),2*exp(sfpar[5])-2)
	return(sfpar)
}

initialBound1d <- function(values,names) {
	stopifnot(values[[2]]>=values[[1]])
	structure(
		list(lower=list(name=names[[1]],value=values[[1]]),
			 upper=list(name=names[[2]],value=values[[2]])),
		class="boundobj1d"
	)
}
