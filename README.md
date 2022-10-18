
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggsdt

<!-- badges: start -->
<!-- badges: end -->

The `ggsdt` package implements generalized gaussian signal detection
theory analysis.

## Installation

You can install the development version of `ggsdt` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kiyomiyoshi/ggsdt")
```

## Example

The `fit_ggsdt` function requires arguments as specified in Maniscalco &
Law’s webpage: <http://www.columbia.edu/~bsm2105/type2sdt/>

The same arguments are also employed in `metaSDT` package:
<https://github.com/craddm/metaSDT>

nR_S1 and nR_S2 are response frequency vectors, ordered from highest
confidence S2 to highest confidence S1.

add_constant = T adds a small value to the response frequency vectors
primarily for avoiding zero-cell-related issues (see above pages).

``` r
library(ggsdt)

nR_S1 <- c(40, 45, 25, 50, 120, 170) # response frequency for S1 stimuli
nR_S2 <- c(240, 70, 20, 30, 50, 40)  # response frequency for S2 stimuli

f1 <- fit_ggsdt(nR_S1, nR_S2, add_constant = F)
f1
#>        mu2   alpha2     beta   LogLike    sigma1    sigma2  kurtosis         X1
#> 1 1.090098 1.306408 1.782542 -1338.253 0.7556205 0.9871485 0.2562262 -0.2240238
#>          X2        X3        X4       X5
#> 1 0.2713651 0.5001792 0.6354313 1.011963
```

``` r
ggdistr(f1[1, 1], f1[1, 2], f1[1, 3])
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="50%" />

``` r
ggroc1(f1[1, 1], f1[1, 2], f1[1, 3])
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="50%" />

``` r
ggzroc1(f1[1, 1], f1[1, 2], f1[1, 3])
```

<img src="man/figures/README-unnamed-chunk-2-3.png" width="50%" />
