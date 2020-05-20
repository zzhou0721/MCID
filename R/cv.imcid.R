#' Selection of the tuning parameters for determining the MCID at the individual level
#'
#'
#' \code{cv.imcid} returns the optimal tuning parameter \eqn{\delta} and \eqn{\lambda} selected from a given grid by using k-fold cross-validation.
#' The tuning parameters are selected for determining the MCID at the individual level
#'
#'
#' @param x a continuous variable denoting the outcome change of interest
#' @param y a binary variable denoting the patient-reported outcome derived from the anchor question
#' @param z a vector or matrix denoting the patient's clinical profiles
#' @param lamseq a vector containing the candidate values for the tuning parameter \eqn{\lambda}, where \eqn{\lambda} is the coefficient of the penalty term, used for avoiding the issue of model overfitting
#' @param delseq a vector containing the candidate values for the tuning parameter \eqn{\delta}, where \eqn{\delta} is used to control the difference between the 0-1 loss and the surrogate loss.
#' We recommend selecting the possible values from the neighborhood of the standard deviation of x
#' @param k the number of groups into which the data should be split to select the tuning parameter \eqn{\delta} by cross-validation. Defaults to 5
#' @param maxit the maximum number of iterations. Defaults to 100
#' @param tol the convergence tolerance. Defaults to 0.01
#'
#' @return a list including the combinations of the selected tuning parameters and the value of the corresponding target function
#' @export
#'
#' @importFrom stats median optim
#' @examples
#' rm(list = ls())
#' n <- 500
#' lambdaseq <- 10 ^ seq(-3, 3, 0.1)
#' deltaseq <- seq(0.1, 0.3, 0.1)
#' a <- 0.1
#' b <- 0.55
#' c <- -0.1
#' d <- 0.45
#'
#' set.seed(721)
#' p <- 0.5
#' y <- 2 * rbinom(n, 1, p) - 1
#' z <- rnorm(n, 1, 0.1)
#' y_1 <- which(y == 1)
#' y_0 <- which(y == -1)
#' x <- c()
#' x[y_1] <- a + z[y_1] * b + rnorm(length(y_1), 0, 0.1)
#' x[y_0] <- c + z[y_0] * d + rnorm(length(y_0), 0, 0.1)
#'
#' sel <- cv.imcid(x = x, y = y, z = z, lamseq = lambdaseq,
#'          delseq = deltaseq, k = 5, maxit = 100, tol = 1e-02)
#' sel$'Selected lambda'
#' sel$'Selected delta'
#' sel$'Function value'
#'
cv.imcid <- function(x, y, z, lamseq, delseq, k = 5, maxit = 100, tol = 1e-02) {
  funval_sum <- matrix(NA, nrow = length(lamseq), ncol = length(delseq))
  for (i in 1:length(lamseq)) {
    for (j in 1:length(delseq)) {
      lambda_i <- lamseq[i]
      delta_j <- delseq[j]
      funval_sum[i, j] <- imcidcv.smooth(x = x, y = y, z = z, lambda = lambda_i, delta = delta_j, fold = k, maxit = maxit, tol = tol)
    }
  }
  index <- which(funval_sum == min(funval_sum), arr.ind = TRUE)
  lamsel <- median(lamseq[index[,1]])
  delsel <- median(delseq[index[,2]])
  funval_min <- min(funval_sum)
  return(list("Selected lambda" = lamsel, "Selected delta" = delsel, "Function value" = funval_min))
}

