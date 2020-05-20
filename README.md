
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MCID

<!-- badges: start -->

<!-- badges: end -->

‘MCID’ is an R package used to provide the point and interval estimation
on the minimal clinically important difference (MCID) at the population
and individual level. For population level, it produces a constant value
for the estimation of MCID. For individual level, the MCID is defined as
a linear function of patient’s clinical characteristics, and the
estimated linear coefficients can be obtained.

## Installation

You can install the released version of MCID from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("zzhou0721/MCID")
```

## Example

### Generate data

``` r
rm(list = ls())
library(MCID)
#> Welcome to MCID package
n <- 500
lambdaseq <- 10 ^ seq(-3, 3, 0.1)
deltaseq <- seq(0.1, 0.3, 0.1)
a <- 0.1
b <- 0.55
c <- -0.1
d <- 0.45

set.seed(115)
p <- 0.5
y <- 2 * rbinom(n, 1, p) - 1
z <- rnorm(n, 1, 0.1)
y_1 <- which(y == 1)
y_0 <- which(y == -1)
x <- c()
x[y_1] <- a + z[y_1] * b + rnorm(length(y_1), 0, 0.1)
x[y_0] <- c + z[y_0] * d + rnorm(length(y_0), 0, 0.1)
```

### Determine MCID at the population level

To determine MCID at the populaton level, we first need to select an
optimal value for the tuning parameter , which is used to control the
difference between 0-1 loss and surrogate loss.

``` r
sel <- cv.pmcid(x = x, y = y, delseq = deltaseq, 
                k = 5, maxit = 100, tol = 1e-02)
delsel <- sel$'Selected delta'
delsel
#> [1] 0.1
```

Then with selected optimal value of , we can determine the point and
interval estimation of MCID at the population level. The confidence
interval is constructed based on the asymptotic normality. In our
simulated data, the true population MCID is 0.5.

``` r
result <- pmcid(x = x, y = y, n = n, delta = delsel, 
                maxit = 100, tol = 1e-02, alpha = 0.05)
result$'Point estimate'
#> [1] 0.4961217
result$'Standard error'
#> [1] 0.0101217
result$'Confidence interval'
#> [1] 0.4762836 0.5159599
```

### Determine MCID at the individual level

To determine MCID at the individual level, we first need to select a
combination of the optimal values for the tuning parameters and . is
used to control the difference between 0-1 loss and surrogate loss. is
the coefficient of the penalty term used to avoid the overfitting issue.

``` r
sel <- cv.imcid(x = x, y = y, z = z, lamseq = lambdaseq, 
                delseq = deltaseq, k = 5, maxit = 100, tol = 1e-02)
lamsel <- sel$'Selected lambda'
delsel <- sel$'Selected delta'
lamsel
#> [1] 0.004496472
delsel
#> [1] 0.2
```

Then with selected and , we can determine the point and interval
estimation for the linear coefficients of the individualized MCID
function. The confidence intervals are constructed based on the
asymptotic normalities. In our simulated data, the true linear
coefficients of the individualized MCID are \_0=0 and \_1=0.5

``` r
result <- imcid(x = x, y = y, z = z, n = n, lambda = lamsel, 
                delta = delsel, maxit = 100, tol = 1e-02, alpha = 0.05)
result$'Point estimates'
#>              [,1]
#> beta0 -0.01551779
#> beta1  0.51193914
result$'Standard errors'
#>            [,1]
#> beta0 0.1292652
#> beta1 0.1227669
result$'Confidence intervals'
#>       Lower bound Upper bound
#> beta0  -0.2688728   0.2378373
#> beta1   0.2713205   0.7525577
```
