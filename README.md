
<!-- README.md is generated from README.Rmd. Please edit that file -->

# braidrm <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The goal of braidrm is to to make the best combination analysis
available more robust, more accessible, and easier to use, so that drug
combinations can be understood more completely and new therapies can be
discovered more quickly.

## Example

This example shows how to fit a BRAID response surface to data, and
print a summary of the resulting fit.

``` r
library(braidrm)

# Fit a basic braid surface
braidFit <- braidrm(measure ~ concA + concB, synergisticExample,
                    model = "kappa2", getCIs=TRUE)

summary(braidFit)
#> Call:
#> braidrm.formula(formula = measure ~ concA + concB, data = synergisticExample, 
#>     model = "kappa2", getCIs = TRUE)
#> 
#>            Lo     Est     Hi
#> IDMA   0.9129  1.0398 1.1477
#> IDMB   0.8690  1.0259 1.1829
#> na     2.3488  2.9116 3.6568
#> nb     2.0612  2.4990 3.1371
#> kappa  1.6104  2.1258 2.6469
#> E0    -0.0762 -0.0300 0.0227
#> EfA    0.9394  1.0080 1.0281
#> EfB    0.9053  0.9848 1.0264
#> Ef         NA  1.0080     NA
```
