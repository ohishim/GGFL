
# GGFL (v0.2.0; in development)

<!-- badges: start -->
<!-- badges: end -->

This package gives GGFL estimates via coordinate descent algorithm.

**cite this package of the latest version**:  
Ohishi, M. (2024).
GGFL: Coordinate optimization for GGFL.
R package version 0.2.0.
https://github.com/ohishim/GGFL


## Installation

You can install the development version of GGFL like so:

``` r
devtools::install_github("ohishim/GGFL")
```

## Example 1

This package has an example dataset as the two objects `exdata` and `adj`.
The `exdata` has 500 individuals and 22 items of which labels are "y", "X1", ..., "X20", and "area".
The "y" is a response variable, "X1", ..., "X20" are explanatory variables, and "area" is a variable of spatial information which has 20 areas.  

This example requires the following package:
```r
library(GGFL)
library(dplyr)
library(magrittr)
```
First, `exdata` is sorted by "area", like this:
```r
exdata1 <- arrange(exdata, area)
```
Moreover, column names of `adj` are redefined as follows:
```r
adj1 <- select(adj, area=areaNO, adj=adjNO)
```
Next, variables required for GGFL are extracted as follwos:
``` r
y <- exdata1$y
X <- select(exdata1, X1:X20) %>% data.matrix
area <- exdata1$area
```
Then, GGFL procedure can be executed as follows:
``` r
res <- GGFL(y, X, area, adj1)
```
In default, `y` and all column vectors of `X` are standardized in the sense of norm.
If you do not want to do it, you can use the option `standardize = FALSE`.
If you want to have the estimates and fitted values, you can respectively obtain as `res$coefGGFL` and `res$fitGGFL`.
Here, these results are for the tuning parameter optimized by EGCV criterion. 
To obtain results for all candidates of tuning parameter, you can use the option `out.all = TRUE`.

## Example 2

In the above example, all variables are penalized.
When only partial variables are penalized, you can use `pGGFL` instead of `GGFL`, like this:
``` r
y <- exdata$y
X <- select(exdata, X1:X15) %>% cbind(X0=1, .)  #variables for penalized estimation
Z <- select(exdata, X16:X20)                    #variables for non-penalized estimation
area <- exdata$area

yli <- split(y, area)
Xli <- lapply(split(X, area), data.matrix)
Zli <- lapply(split(Z, area), data.matrix)
D <- split(adj$adjNO, adj$areaNO)

res <- pGGFL(yli, Xli, Zli, D)
```

## Reference

1. Ohishi, M., Okamura, K., Itoh, Y. & Yanagihara, H. (2021).
Coordinate descent algorithm for generalized group fused Lasso.
*Hiroshima Statistical Research Group Technical Report*, TR-No. 21-02, Hiroshima University.
[[PDF](http://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR21/TR21-02.pdf)]
