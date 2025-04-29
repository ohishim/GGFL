
# GGFL (v1.0.2)

<!-- badges: start -->
<!-- badges: end -->

This package gives GGFL estimates via coordinate descent algorithm.

**cite this package of the latest version**:  
Ohishi, M. (2025).
GGFL: Coordinate optimization for GGFL.
R package version 1.0.2.
https://github.com/ohishim/GGFL  

**information**:  

- (2024/08/14) <span style="color: red; ">attention !</span> There is a mistake in the v0.3.0 and v0.3.1. The mistake has been corrected in v0.4.0.


## Installation

You can install the current version of GGFL like so:

``` r
devtools::install_github("ohishim/GGFL")
```

## Example 1

This package has an example dataset as the two objects `exdata` and `adj`.
The `exdata` has 500 individuals and 22 items of which labels are "y", "X1", ..., "X20", and "group".
The "y" is a response variable, "X1", ..., "X20" are explanatory variables, and "group" is a variable of spatial information which has 20 areas.  

The GGFL procedure can be executed as follows:
``` r
library(GGFL)

y <- exdata$y
X <- dplyr::select(exdata, X1:X20)
group <- exdata$group

res <- GGFL(y, X, group, adj)
```
In default, `y` and all column vectors of `X` are standardized in the sense of norm.
If you do not want to do it, you can use the option `standardize = FALSE`.
If you want to have the estimates and fitted values, you can respectively obtain as `res$coefficients` and `res$fitted.values`.
Here, these results are for the tuning parameter optimized by EGCV criterion. 
To obtain results for all candidates of tuning parameter, you can use the option `out.all = TRUE`.

## Example 2

You can partially apply GGFL like this:
``` r
y <- exdata$y
X <- dplyr::select(exdata, X1:X15)   # for local estimation
Z <- dplyr::select(exdata, X16:X20)  # for global estimation
group <- exdata$group

res <- pGGFL(y, X, Z, group, adj)
```

## Reference

1. Ohishi, M., Okamura, K., Itoh, Y., Wakaki, H. & Yanagihara, H. (2025).
Coordinate descent algorithm for generalized group fused Lasso.
*Behaviormetrika*, **52** (1), 105--137.
doi: [10.1007/s41237-024-00233-6](https://doi.org/10.1007/s41237-024-00233-6)
