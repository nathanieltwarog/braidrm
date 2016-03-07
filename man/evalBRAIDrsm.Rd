\name{evalBRAIDrsm}
\alias{evalBRAIDrsm}
\title{ Calculate BRAID Surface Values }
\description{
Calculates the value of the Bivariate Response to Additive Interacting Doses (BRAID) surface model for the given concentration
pairs.
}
\usage{
evalBRAIDrsm(DA, DB, parv)
}
\arguments{
  \item{DA}{ vector of concentrations of drug A }
  \item{DB}{ vector of concentrations of drug B }
  \item{parv}{ ten-element vector specifying the full set of parameters for the BRAID surface (see Details below) }
}
\details{
The full ten-parameter BRAID model, which we refer to as the \emph{extended BRAID} or \emph{eBRAID} model is defined as:
\deqn{E(D_A,D_B) = E_0 + \frac{E_f-E_0}{1+{{\tilde{D}}_{AB}}^{-\delta\sqrt{n_a n_b}}}}{E(DA,DB) = E0 + (Ef-E0)/(1+<DAB>^(-\delta*\sqrt(na*nb)))}
\deqn{{\tilde{D}}_{AB} = {{\tilde{D}}_A}^\frac{1}{\delta\sqrt{n_a n_b}}+{{\tilde{D}}_B}^\frac{1}{\delta\sqrt{n_a n_b}}
		+\kappa\sqrt{{{\tilde{D}}_A}^\frac{1}{\delta\sqrt{n_a n_b}}
		{{\tilde{D}}_B}^\frac{1}{\delta\sqrt{n_a n_b}}}}{<DAB> = <DA>^(1/(\delta*\sqrt(na*nb)))
		+<DB>^(1/(\delta*\sqrt(na*nb)))+\kappa*\sqrt((<DA>^(1/(\delta*\sqrt(na*nb))))*(<DB>^(1/(\delta*\sqrt(na*nb)))))}
\deqn{{\tilde{D}}_A = \frac{\left(\frac{E_{f,A}-E_0}{E_f-E_0}\right)\left(\frac{D_A}{{ID}_{M,A}}\right)^{n_a}}
		{1+\left(1-\frac{E_{f,A}-E_0}{E_f-E_0}\right)\left(\frac{D_A}{{ID}_{M,A}}\right)^{n_a}}}{<DA> = 
		(((EfA-E0)/(Ef-E0))*((DA/IDMA)^na))/(1+(1-(EfA-E0)/(Ef-E0))*((DA/IDMA)^na))}
\deqn{{\tilde{D}}_B = \frac{\left(\frac{E_{f,B}-E_0}{E_f-E_0}\right)\left(\frac{D_B}{{ID}_{M,B}}\right)^{n_b}}
		{1+\left(1-\frac{E_{f,B}-E_0}{E_f-E_0}\right)\left(\frac{D_B}{{ID}_{M,B}}\right)^{n_b}}}{<DB> 
		(((EfB-E0)/(Ef-E0))*((DB/IDMB)^nb))/(1+(1-(EfB-E0)/(Ef-E0))*((DB/IDMB)^nb))}
The parameters of this equation must satisfy \eqn{n_a>0}{na>0}, \eqn{n_b>0}{nb>0}, \eqn{\delta>0}, \eqn{\kappa> -2},
\eqn{sign(E_f-E_0)=sign(E_{f,A}-E_0)=sign(E_{f,B}-E_0)}{sign(Ef-E0)=sign(EfA-E0)=sign(EfB-E0)},
\eqn{|E_f-E_0|\geq|E_{f,A}-E_0|}{|Ef-E0|\ge|EfA-E0|}, and \eqn{|E_f-E_0|\geq|E_{f,B}-E_0|}{|Ef-E0|\ge|EfB-E0|}. With this
definition, the ten-element parameter vector is [\eqn{{ID}_{M,A}}{IDMA}, \eqn{{ID}_{M,B}}{IDMB}, \eqn{n_a}{na}, \eqn{n_b}{nb},
\eqn{\delta}, \eqn{\kappa}, \eqn{E_0}{E0}, \eqn{E_{f,A}}{EfA}, \eqn{E_{f,B}}{EfB}, \eqn{E_f}{Ef}]. The simpler standard BRAID
model, as described in Twarog \emph{et al.} is obtained by setting \eqn{\delta} equal to 1 and setting \eqn{E_f}{Ef} such that
\eqn{|E_f-E_0|}{|Ef-E0|} is equal to the maximum of \eqn{|E_{f,A}-E_0|}{|EfA-E0|} and \eqn{|E_{f,B}-E_0|}{|EfB-E0|}.  Assuming that
this sets \eqn{E_f}{Ef} equal to \eqn{E_{f,A}}{EfA}, this causes the equation for \eqn{{\tilde{D}}_A}{<DA>} to simplify to
\deqn{{\tilde{D}}_A = \left(\frac{D_A}{{ID}_{M,A}}\right)^{n_a}}{<DA> = (DA/IDMA)^na}
}
\value{ A vector of response values corresponding to the pairs of concentrations in \code{DA} and \code{DB} }
\author{
Nathaniel R. Twarog
}
\examples{
conc1 <- rep(seq(0,3*10^-6,length=50),each=50)
conc2 <- rep(seq(0,3*10^-6,length=50),times=50)

# Additive surface
act <- evalBRAIDrsm(conc1,conc2,parv=c(10^-6,10^-6,1.5,1.5,1,0,0,100,100,100))
# A BRAID additive surface is not a Loewe additive surface
act <- evalBRAIDrsm(conc1,conc2,parv=c(10^-6,10^-6,1,3,1,0,0,100,100,100))
# BRAID antagonism
act <- evalBRAIDrsm(conc1,conc2,parv=c(10^-6,10^-6,1.5,1.5,1,-1,0,100,100,100))
# delta-BRAID synergy
act <- evalBRAIDrsm(conc1,conc2,parv=c(10^-6,10^-6,1.5,1.5,1.75,0,0,100,100,100))
# Differing final effects
act <- evalBRAIDrsm(conc1,conc2,parv=c(10^-6,10^-6,1.5,1.5,1,0,0,75,100,100))
}

\seealso{
\code{\link{braidrm}}, \code{\link{invertBRAIDrsm}}
}
\keyword{ regression }
