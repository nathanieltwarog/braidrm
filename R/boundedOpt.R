# boundedOpt.R
# Included for compatibility with optimization function from `basicdrm`

runBoundedOptim <- function(vdfunc,fpfunc,nstart,nbounds,parscale=NULL) {
	derivs <- NULL
	parfunc <- function(p) {
		p <- pmin(pmax(p,nbounds[1,]),nbounds[2,])
		vd <- vdfunc(p)
		derivs <<- vd$derivatives
		return(vd$value)
	}
	derivfunc <- function(p) { return(derivs) }

	if (is.null(parscale)) { control <- list(maxit=1000) }
	else { control <- list(maxit=1000,parscale=parscale) }

	nls <- stats::optim(nstart,parfunc,derivfunc,method="L-BFGS-B",lower=nbounds[1,],
						upper=nbounds[2,],control=control)
	nls$par <- pmin(pmax(nls$par,nbounds[1,]),nbounds[2,])
	nls$fullpar <- fpfunc(nls$par)
	return(nls)
}

boundedOpt2d <- function(mnvec,bounds) {
	if (is.null(bounds)) { return(NULL) }
	if (!inherits(bounds,"boundobj2d")) {
		stop("Parameter 'bounds' must be of class 'boundobj2d'.")
	}

	SA <- mnvec[2]
	SB <- mnvec[1]-SA
	SC <- 1-2*SB-SA
	SD <- mnvec[4]
	SE <- mnvec[3]-SD
	OE0 <- (SB*SD-SA*SE)/(SB^2-SA*SC)
	OEf <- (SB*SE-SC*SD)/(SB^2-SA*SC)

	# If minimum lies within bounds, return it
	edges <- bounds$edges
	if ((SB^2-SA*SC)!=0 && !any((OE0*edges[,1]+OEf*edges[,2])<edges[,3])) { return(c(OE0,OEf,0,0,0)) }

	# Otherwise, find minimum on violated edges
	verts <- bounds$verts
	if ((SB^2-SA*SC)!=0) { viol <- which((OE0*edges[,1]+OEf*edges[,2])<edges[,3]) }
	else {
		den <- mnvec[1]*edges[,1]-(1-mnvec[1])*edges[,2]
		if ((mnvec[1]>0 && any(den==0 & (mnvec[1]*edges[,3]>mnvec[3]*edges[,2])))||
			(mnvec[1]==0 && any(den==0 & (edges[,3]>mnvec[3]*edges[,1])))) {
			if (mnvec[1]>0) { rel <- which(den==0 & mnvec[1]*edges[,3]>mnvec[3]*edges[,2]) }
			else { rel <- which(den==0 & edges[,3]>mnvec[3]*edges[,1]) }
			if (is.infinite(verts[rel,1]) && is.infinite(verts[(if (rel==nrow(verts)) 1 else (rel+1)),1])) {
				if (edges[rel,2]!=0) { return(c(0,edges[rel,3]/edges[rel,1],0,0,0)) }
				else { return(c(edges[rel,3]/edges[rel,2],0,0,0,0)) }
			} else {
				if (is.finite(verts[rel,1])) { outv <- c(verts[rel,],rel-0.5) }
				else { outv <- c(verts[(if (rel==nrow(verts)) 1 else (rel+1)),],(if (rel==nrow(verts)) 1 else (rel+1))-0.5) }
				return(c(outv,coefObjective(outv[1:2],mnvec)[2:3]))
			}
		}
		if (all(den==0)) {
			if (mnvec[1]>0) { return(c(0,mnvec[3]/mnvec[1],0,0,0)) }
			else { return(c(mnvec[3],0,0,0,0)) }
		}
		viol <- which(den!=0)
	}
	cfp <- cbind(edges[viol,1]^2,edges[viol,1]*edges[viol,2],edges[viol,2]^2,
				 edges[viol,1]*edges[viol,3],edges[viol,2]*edges[viol,3])
	num1 <- -cfp[,5]*mnvec[1]+(cfp[,4]+cfp[,5])*mnvec[2]+cfp[,3]*mnvec[3]-(cfp[,2]+cfp[,3])*mnvec[4]
	num2 <- -(cfp[,4]+2*cfp[,5])*mnvec[1]+(cfp[,4]+cfp[,5])*mnvec[2]-cfp[,2]*mnvec[3]+(cfp[,1]+cfp[,2])*mnvec[4]+cfp[,5]
	den <- -2*(cfp[,2]+cfp[,3])*mnvec[1]+(cfp[,1]+2*cfp[,2]+cfp[,3])*mnvec[2]+cfp[,3]
	ccf <- cbind(num1/den,num2/den,NA,NA,NA)
	if (nrow(edges)==1) { return(c(ccf[1,1:2],1,coefObjective(ccf[1,1:2],mnvec)[2:3])) }
	for (i in 1:length(viol)) {
		nbrs <- c((if (viol[i]==1) nrow(edges) else (viol[i]-1)),(if (viol[i]==nrow(edges)) 1 else (viol[i]+1)))
		if (all((ccf[i,1]*edges[nbrs,1]+ccf[i,2]*edges[nbrs,2])>=edges[nbrs,3])) {
			ccf[i,3:5] <- coefObjective(ccf[i,1:2],mnvec)
		}
	}
	if (any(!is.na(ccf[,3]))) {
		bnd <- which.min(ccf[,3])
		return(c(ccf[bnd,1:2],viol[bnd],ccf[bnd,4:5]))
	}

	bndmat <- cbind(verts,NA,NA,NA)
	for (i in 1:nrow(bndmat)) {
		if (is.finite(bndmat[i,1])) { bndmat[i,3:5] <- coefObjective(bndmat[i,1:2],mnvec) }
	}
	bnd <- which.min(bndmat[,3])
	return(c(bndmat[bnd,1:2],bnd-0.5,bndmat[bnd,4:5]))
}

coefObjective <- function(cf,mnv,deriv=TRUE) {
	obj <- (2*(mnv[4]-mnv[3])*cf[1]-2*mnv[4]*cf[2]+(1-2*mnv[1]+mnv[2])*cf[1]^2
			+2*(mnv[1]-mnv[2])*cf[1]*cf[2]+mnv[2]*cf[2]^2)
	if (deriv) { obj <- c(obj,2*((mnv[4]-mnv[3])+(1-2*mnv[1]+mnv[2])*cf[1]+(mnv[1]-mnv[2])*cf[2]),
						  2*(-mnv[4]+(mnv[1]-mnv[2])*cf[1]+mnv[2]*cf[2])) }
	return(obj)
}

emptyBound2d <- function() {
	structure(list(),class="boundobj2d")
}

addBound2d <- function(bounds,newbound,name="") {
	if (is.null(bounds)) { return(NULL) }
	if (!inherits(bounds,"boundobj2d")) {
		stop("Parameter 'bounds' must be of class 'boundobj2d'.")
	}

	if (name!="") {
		newbound <- rbind(newbound)
		row.names(newbound) <- name
	}
	nboundob <- addBound2dInt(unclass(bounds),newbound)
	if (is.null(nboundob)) { return(nboundob) }
	nverts <- nboundob$verts
	if (nrow(nverts)>1) {
		rem <- which(is.finite(nverts[,1]) & nverts[,1]==nverts[c(2:nrow(nverts),1),1]
					 & nverts[,2]==nverts[c(2:nrow(nverts),1),2])
		if (length(rem)>0 & length(rem)<nrow(nverts)) {
			nboundob$edges <- nboundob$edges[-rem,]
			nboundob$verts <- nboundob$verts[-rem,]
		}
	}
	return(structure(nboundob,class="boundobj2d"))
}

addBound2dInt <- function(boundob,newbound) {
	if (is.null(boundob$edges)) { return(list(edges=rbind(newbound),verts=rbind(c(Inf,Inf)))) }
	if (newbound[1]==0 && newbound[2]==0) { return(boundob) }
	bounds <- boundob$edges
	nrb <- nrow(bounds)
	nbv <- 1:nrb
	verts <- boundob$verts
	vsign <- verts[,1]*newbound[1]+verts[,2]*newbound[2]-newbound[3]
	if (all(!is.finite(vsign))||all(vsign[is.finite(vsign)]<0)) {
		# All vertices negative or nonfinite
		if (all(is.finite(vsign))) { return(NULL) }
		if (nrb==1) {
			nvert <- interpoints(bounds,newbound)
			if (!is.finite(nvert[1,1])) {
				if (is.nan(nvert[1,1]) && sign(bounds[1,1]*newbound[1]+bounds[1,2]*newbound[2])>0) { return(boundob) }
				if (!is.nan(nvert[1,1]) && nvert[1,1]<0 && !is.nan(nvert[1,2]) && nvert[1,2]<0) { return(NULL) }
				else if (!is.nan(nvert[1,1]) && nvert[1,1]<0) { return(list(edges=rbind(newbound),verts=rbind(c(Inf,Inf)))) }
				else if (!is.nan(nvert[1,2]) && nvert[1,2]<0) { return(boundob) }
				else { nvert[1,] <- c(Inf,Inf) }
			}
			return(list(edges=rbind(bounds,newbound),verts=rbind(verts,nvert)))
		}
		nvert <- interpoints(bounds[c(1,nrb),],newbound)
		if (!any(is.finite(nvert[,1]))) {
			nvert[is.nan(nvert)] <- 0
			if (all(nvert[,1]>=0)) { return(boundob) }
			else if (all(nvert[,1]<0)) { return(NULL) }
			subind <- (c(1,nrb))[which(nvert[,1]<0)]
			bounds[subind,] <- newbound
			if (subind==1) {
				newvert <- interpoints(rbind(bounds[2,]),newbound)
				newvert[!is.finite(newvert)] <- Inf
				verts[2,] <- newvert
			} else {
				newvert <- interpoints(rbind(bounds[nrb-1,]),newbound)
				newvert[!is.finite(newvert)] <- Inf
				verts[nrb-1,] <- newvert
			}
			return(list(edges=bounds,verts=verts))
		} else if (!is.finite(nvert[1,1])) {
			if (!is.nan(nvert[1,2])&& nvert[1,2]<0) { return(NULL) }
			nvert[1,] <- Inf
			return(list(edges=rbind(newbound,bounds[nrb,,drop=FALSE]),verts=nvert))
		} else if (!is.finite(nvert[2,1])) {
			if (!is.nan(nvert[2,2])&& nvert[2,2]<0) { return(NULL) }
			nvert[2,] <- Inf
			return(list(edges=rbind(bounds[1,,drop=FALSE],newbound),verts=nvert))
		}
		if (is.finite(verts[2,1]) && (nvert[1,1]*bounds[2,1]+nvert[1,2]*bounds[2,2]<=bounds[2,3])) { return(NULL) }
		return(list(edges=rbind(bounds[nbv==1,,drop=FALSE],newbound,bounds[nbv==nrb,,drop=FALSE]),verts=rbind(c(Inf,Inf),nvert)))
	} else if (all(vsign[is.finite(vsign)]>=0)) {
		# All vertices positive or nonfinite, at least one positive
		if (all(is.finite(vsign))) { return(boundob) }
		nvert <- interpoints(bounds[c(1,nrb),],newbound)
		if (!any(is.finite(nvert[,1]))) { return(boundob) }
		else if (!is.finite(nvert[1,1])) {
			if (!is.nan(nvert[1,2])&& nvert[1,2]<0) { return(boundob) }
			return(list(edges=rbind(bounds,newbound),verts=rbind(verts,nvert[2,])))
		} else if (!is.finite(nvert[2,1])) {
			if (!is.nan(nvert[2,2])&& nvert[2,2]<0) { return(boundob) }
			nvert[2,] <- Inf
			return(list(edges=rbind(newbound,bounds),verts=rbind(nvert[c(2,1),],verts[-1,])))
		}
		if (is.finite(verts[2,1]) && (nvert[1,1]*bounds[2,1]+nvert[1,2]*bounds[2,2]<=bounds[2,3])) { return(boundob) }
		return(list(edges=rbind(newbound,bounds),verts=rbind(nvert[c(2,1),],verts[-1,])))
	} else if (all(is.finite(vsign)) || (vsign[2]>=0 && vsign[nrb]>=0) || (vsign[2]<0 && vsign[nrb]<0)) {
		# Bound intersects two finite edges
		if ((is.finite(vsign[1])&&vsign[1]<0)||(!is.finite(vsign[1])&&vsign[2]<0)) {
			pos <- which(vsign>=0)
			nvert <- interpoints(bounds[c(max(pos),min(pos)-1),],newbound)
			return(list(edges=rbind(newbound,bounds[c(min(pos)-1,pos),,drop=FALSE]),verts=rbind(nvert,verts[pos,])))
		} else {
			neg <- which(vsign<0)
			nvert <- interpoints(bounds[c(min(neg)-1,max(neg)),],newbound)
			newverts <- rbind(verts[1:(min(neg)-1),],nvert)
			if (max(neg)<nrow(verts)) { newverts <- rbind(newverts,verts[(max(neg)+1):nrow(verts),]) }
			return(list(edges=rbind(bounds[1:(min(neg)-1),,drop=FALSE],newbound,bounds[max(neg):nrb,,drop=FALSE]),verts=newverts))
		}
	} else {
		# New bound splits existing verticies in two
		if (vsign[2]>=0) {
			mxv <- max(which(vsign>=0))
			nvert <- interpoints(bounds[c(1,mxv,nrb),],newbound)
			nvert[!is.finite(nvert)] <- Inf
			if (is.finite(nvert[1,1]) && (nvert[1,1]*bounds[2,1]+nvert[1,2]*bounds[2,2]>=bounds[2,3])) {
				return(list(edges=rbind(newbound,bounds[1:mxv,,drop=FALSE]),verts=rbind(nvert[c(2,1),],verts[2:mxv,])))
			} else if (is.finite(nvert[3,1]) && (nvert[3,1]*bounds[nrb-1,1]+nvert[3,2]*bounds[nrb-1,2]>=bounds[nrb-1,3])) {
				return(list(edges=rbind(bounds[1:mxv,,drop=FALSE],newbound,bounds[nrb,,drop=FALSE]),
							verts=rbind(verts[1:mxv,],nvert[2:3,])))
			} else { return(list(edges=rbind(bounds[1:mxv,,drop=FALSE],newbound),verts=rbind(verts[1:mxv,],nvert[2,]))) }
		} else {
			mnv <- min(which(vsign>=0))
			nvert <- interpoints(bounds[c(1,mnv-1,nrb),],newbound)
			if (is.finite(nvert[2,1]) && (nvert[2,1]*bounds[nrb-1,1]+nvert[2,2]*bounds[nrb-1,2]>=bounds[nrb-1,3])) {
				return(list(edges=rbind(newbound,bounds[(mnv-1):nrb,,drop=FALSE]),
							verts=rbind(nvert[c(3,2),],verts[mnv:nrow(verts),])))
			} else if (is.finite(nvert[1,1]) && (nvert[1,1]*bounds[2,1]+nvert[1,2]*bounds[2,2]>=bounds[2,3])) {
				return(list(edges=rbind(bounds[1,,drop=FALSE],newbound,bounds[(mnv-1):nrb,,drop=FALSE]),
							verts=rbind(verts[1,],nvert[1:2,],verts[mnv:nrow(verts),])))
			} else {
				return(list(edges=rbind(newbound,bounds[(mnv-1):nrb,,drop=FALSE]),
							verts=rbind(verts[1,],nvert[2,],verts[mnv:nrow(verts),])))
			}
		}
	}
}

interpoints <- function(bounds,newbound) {
	den <- newbound[1]*bounds[,2]-bounds[,1]*newbound[2]
	nvert <- cbind((newbound[3]*bounds[,2]-bounds[,3]*newbound[2])/den,
				   (newbound[1]*bounds[,3]-bounds[,1]*newbound[3])/den)
	if (any(den==0)) {
		rel <- den==0&bounds[,1]!=0
		nvert[rel,1] <- sign(bounds[rel,3]*newbound[1]/bounds[rel,1]-newbound[3])*Inf
		rel <- den==0&bounds[,1]==0
		nvert[rel,1] <- sign(bounds[rel,3]*newbound[2]/bounds[rel,2]-newbound[3])*Inf
		if (newbound[1]!=0) { nvert[den==0,2] <- sign(newbound[3]*bounds[den==0,1]/newbound[1]-bounds[den==0,3])*Inf }
		else { nvert[den==0,2] <- sign(newbound[3]*bounds[den==0,2]/newbound[2]-bounds[den==0,3])*Inf }
	}
	return(nvert)
}

boundedOpt1d <- function(mnvec,other,bounds) {
	if (is.null(bounds)) { return(NULL) }
	if (is.numeric(bounds)) {
		basicbounds <- bounds
	} else if (inherits(bounds,"boundobj1d")) {
		basicbounds <- c(bounds$lower$value,bounds$upper$value)
	} else {
		stop("Parameter 'bounds' must be of class 'boundobj1d'.")
	}

	if (mnvec[2]==0) {
		if (is.finite(basicbounds[[1]])) { return(c(basicbounds[1],1,0)) }
		else if (is.finite(basicbounds[[2]])) { return(c(basicbounds[2],-1,0)) }
		else { return(c(other,0,0)) }
	}
	E <- (1-mnvec[1]/mnvec[2])*other+mnvec[4]/mnvec[2]
	if (E<basicbounds[1]) {
		return(c(basicbounds[1],1,2*(-mnvec[4]+basicbounds[1]*mnvec[2]+other*(mnvec[1]-mnvec[2]))))
	} else if (E>basicbounds[2]) {
		return(c(basicbounds[2],-1,2*(-mnvec[4]+basicbounds[2]*mnvec[2]+other*(mnvec[1]-mnvec[2]))))
	} else { return(c(E,0,0)) }
}

emptyBound1d <- function() {
	structure(
		list(lower=list(name="lower",value=-Inf),
			 upper=list(name="upper",value=Inf)),
		class="boundobj1d"
	)
}

addBound1d <- function(bounds,newbound,name="") {
	if (is.null(bounds)) { return(NULL) }
	if (!inherits(bounds,"boundobj1d")) {
		stop("Parameter 'bounds' must be of class 'boundobj1d'.")
	}

	if (name!="") { names(newbound) <- c(name,"Dir") }

	basicbounds <- c(bounds$lower$value,bounds$upper$value)
	if (newbound[[2]]<0) {
		if (newbound[[1]] < basicbounds[[1]]) { return(NULL) }
		else if (newbound[[1]] >= basicbounds[[2]]) { return(bounds) }
		else {
			if (name=="") { name <- "upper"}
			newbounds <- structure(
				list (lower=bounds$lower,
					  upper=list(name=name,value=newbound[[1]])),
				class="boundobj1d"
			)
		}
	} else if (newbound[[2]]>0) {
		if (newbound[[1]] > basicbounds[[2]]) { return(NULL) }
		else if (newbound[[1]] <= basicbounds[[1]]) { return(bounds) }
		else {
			if (name=="") { name <- "lower"}
			newbounds <- structure(
				list (lower=list(name=name,value=newbound[[1]]),
					  upper=bounds$upper),
				class="boundobj1d"
			)
		}
	} else {
		stop("The second value of the parameter 'newbound' should be -1",
			 " (indicating an upper bound) or 1 (indicating a lower bound).")
	}
	newbounds
}
