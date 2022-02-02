
# GGFL

<!-- badges: start -->
<!-- badges: end -->

This package gives GGFL estimates via coordinate descent algorithm.

## Installation

You can install the development version of GGFL like so:

``` r
devtools::install_github("ohishim/GGFL")
```

## Example 1

This package has an example dataset as the two objects `exdata` and `adj`.
The `exdata` has 500 individuals and 22 items of which labels are "y", "X1", ..., "X20", and "area".
The "y" is a response variable, "X1", ..., "X20" are explanatory variables, and "area" is a variable of spatial information which has 20 areas.  

First, packages are read, and each variable is extracted, like this:
``` r
library(GGFL)
library(dplyr)
library(magrittr)

y <- exdata$y
X <- select(exdata, X1:X20) %>% cbind(X0=1, .)
area <- exdata$area
```
Then, GGFL procedure can be executed as follows:
``` r
yli <- split(y, area)
Xli <- lapply(split(X, area), data.matrix)
D <- split(adj$adjNO, adj$areaNO)

res <- GGFL.cda(yli, Xli, D)
```
If you want to have the estimates and predictive values, you can respectively obtain as `res$coef.opt` and `res$pred`.
Here, these results are for the tuning parameter optimized by EGCV criterion. 
To obtain results for all candidates of tuning parameter, you can use the option `out.all = TRUE`.

## Example 2

In the above example, all variables are penalized.
When only partial variables are penalized, you can use `GGFLa.cda` instead of `GGFL.cda`, like this:
``` r
y <- exdata$y
X <- select(exdata, X1:X15) %>% cbind(X0=1, .)  #variables for penalized estimation
Z <- select(exdata, X16:X20)                    #variables for non-penalized estimation
area <- exdata$area

yli <- split(y, area)
Xli <- lapply(split(X, area), data.matrix)
Zli <- lapply(split(Z, area), data.matrix)
D <- split(adj$adjNO, adj$areaNO)

res <- GGFLa.cda(yli, Xli, Zli, D)
```


