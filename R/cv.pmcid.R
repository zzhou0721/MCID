#' Selection of the tuning parameter for determining the MCID at the population level
#'
#'
#' \code{cv.pmcid} returns the optimal tuning parameter \eqn{\delta} selected from a given grid by using k-fold cross-validation.
#' The tuning parameter is selected for determining the MCID at the population level
#'
#'
#' @param x a continuous variable denoting the outcome change of interest
#' @param y a binary variable indicating the patient-reported outcome derived from the anchor question
#' @param delseq a vector containing the candidate values for the tuning parameter \eqn{\delta}, where \eqn{\delta} is used to control the difference between the 0-1 loss and the surrogate loss.
#' We recommend selecting the possible values from the neighborhood of the standard deviation of x
#' @param k the number of groups into which the data should be split to select the tuning parameter \eqn{\delta} by cross-validation. Defaults to 5
#' @param maxit the maximum number of iterations. Defaults to 100
#' @param tol the convergence tolerance. Defaults to 0.01
#'
#' @return a list including the selected tuning parameter and the value of the corresponding target function
#' @export
#'
#' @importFrom stats median optim
#' @examples
#' rm(list = ls())
#' n <- 500
#' deltaseq <- seq(0.1, 1, 0.1)
#' a <- 0.2
#' b <- -0.1
#' p <- 0.5
#'
#' set.seed(115)
#' y <- 2 * rbinom(n, 1, p) - 1
#' y_1 <- which(y == 1)
#' y_0 <- which(y == -1)
#' x <- c()
#' x[y_1] <- rnorm(length(y_1), a, 0.1)
#' x[y_0] <- rnorm(length(y_0), b, 0.1)
#'
#' sel <- cv.pmcid(x = x, y = y, delseq = deltaseq, k = 5,
#'          maxit = 100, tol = 1e-02)
#' sel$'Selected delta'
#' sel$'Function value'
#'
cv.pmcid <- function(x, y, delseq, k = 5, maxit = 100, tol = 1e-02) {
  funval_sum <- c()
  for (i in 1:length(delseq)) {
    delta_i <- delseq[i]
    funval_sum[i] <- pmcidcv.smooth(x = x, y = y, delta = delta_i, fold = k, maxit = maxit, tol = tol)
  }
  index <- which(funval_sum == min(funval_sum))
  delsel <- median(delseq[index])
  funval_min <- min(funval_sum)
  return(list("Selected delta" = delsel, "Function value" = funval_min))
}

