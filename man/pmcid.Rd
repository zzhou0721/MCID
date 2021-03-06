% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmcid.R
\name{pmcid}
\alias{pmcid}
\title{Point and interval estimation for the MCID at the population level}
\usage{
pmcid(x, y, n, delta, maxit = 100, tol = 0.01, alpha = 0.05)
}
\arguments{
\item{x}{a continuous variable denoting the outcome change of interest}

\item{y}{a binary variable indicating the patient-reported outcome derived from the anchor question}

\item{n}{the sample size}

\item{delta}{the selected tuning parameter \eqn{\delta}, can be returned by \code{cv.pmcid}}

\item{maxit}{the maximum number of iterations. Defaults to 100}

\item{tol}{the convergence tolerance. Defaults to 0.01}

\item{alpha}{nominal level of the confidence interval. Defaults to 0.05}
}
\value{
a list including the point estimate of the population MCID and its standard error, and the confidence interval based on the asymptotic normality
}
\description{
\code{pmcid} returns the point estimate for the MCID at the population level
}
\examples{
rm(list = ls())
n <- 500
deltaseq <- seq(0.1, 1, 0.1)
a <- 0.2
b <- -0.1
p <- 0.5
### True MCID is 0.5 ###

set.seed(115)
y <- 2 * rbinom(n, 1, p) - 1
y_1 <- which(y == 1)
y_0 <- which(y == -1)
x <- c()
x[y_1] <- rnorm(length(y_1), a, 0.1)
x[y_0] <- rnorm(length(y_0), b, 0.1)

sel <- cv.pmcid(x = x, y = y, delseq = deltaseq, k = 5,
         maxit = 100, tol = 1e-02)
delsel <- sel$'Selected delta'

result <- pmcid(x = x, y = y, n = n, delta = delsel,
            maxit = 100, tol = 1e-02, alpha = 0.05)
result$'Point estimate'
result$'Standard error'
result$'Confidence interval'


}
