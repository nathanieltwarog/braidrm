library(tidyverse)
library(braidrm)

set.seed(20240801)

additivePar <- c(
	1, 1, 3, 3,
	0, # kappa = 0 indicates additive
	0, 1, 1, 1
)
synergisticPar <- c(
	1, 1, 3, 3,
	2, # kappa > 0 indicates synergy
	0, 1, 1, 1
)
antagonisticPar <- c(
	1, 1, 3, 3,
	-1, # kappa < 0 indicates antagonism
	0, 1, 1, 1
)
incompletePar <- c(
	1, 100, 3, 3, 0,
	0, 1, 0.1, 1
)
protectivePar <- c(
	0.5, 2, 3, 3, 2,
	0, 1, 0, 0
)
oppositionalPar <- c(
	1, 1, 3, 3, -0.5,
	0, 1, -0.5, 1
)
coactivePar <- c(
	0.3, 0.3, 3, 3, 0,
	0, 0, 0, 1
)

measured <- c(0, 2^(-3:3))
basicExample <- expand_grid(concA=measured,concB=measured)

additiveExample <- basicExample %>%
	mutate(truth=evalBraidModel(concA,concB,additivePar),
		   measure=truth+rnorm(n(),sd=0.07))

synergisticExample <- basicExample %>%
	mutate(truth=evalBraidModel(concA,concB,synergisticPar),
		   measure=truth+rnorm(n(),sd=0.07))

antagonisticExample <- basicExample %>%
	mutate(truth=evalBraidModel(concA,concB,antagonisticPar),
		   measure=truth+rnorm(n(),sd=0.07))

incompleteExample <- basicExample %>%
	mutate(truth=evalBraidModel(concA,concB,incompletePar),
		   measure=truth+rnorm(n(),sd=0.07))

protectiveExample <- basicExample %>%
	mutate(truth=evalFlippedBraidModel(concA,concB,protectivePar,flip="A"),
		   measure=truth+rnorm(n(),sd=0.07))

oppositionalExample <- basicExample %>%
	mutate(truth=evalFlippedBraidModel(concA,concB,oppositionalPar,flip="B"),
		   measure=truth+rnorm(n(),sd=0.07))

coactiveExample <- basicExample %>%
	mutate(truth=evalFlippedBraidModel(concA,concB,coactivePar,flip="both"),
		   measure=truth+rnorm(n(),sd=0.07))

usethis::use_data(additiveExample, overwrite =TRUE)
usethis::use_data(synergisticExample, overwrite =TRUE)
usethis::use_data(antagonisticExample, overwrite =TRUE)
usethis::use_data(incompleteExample, overwrite =TRUE)
usethis::use_data(protectiveExample, overwrite =TRUE)
usethis::use_data(oppositionalExample, overwrite =TRUE)
usethis::use_data(coactiveExample, overwrite =TRUE)
