#' Point and interval estimation for the MCID at the population level
#'
#' \code{pmcid} returns the point estimate for the MCID at the population level
#'
#' @param x a continuous variable denoting the outcome change of interest
#' @param y a binary variable indicating the patient-reported outcome derived from the anchor question
#' @param n the sample size
#' @param delta the selected tuning parameter \eqn{\delta}, can be returned by \code{cv.pmcid}
#' @param maxit the maximum number of iterations. Defaults to 100
#' @param tol the convergence tolerance. Defaults to 0.01
#' @param alpha nominal level of the confidence interval. Defaults to 0.05
#'
#' @return a list including the point estimate of the population MCID and its standard error, and the confidence interval based on the asymptotic normality
#'
#' @importFrom stats optim qnorm
#'
#' @export
#'
#' @examples
#' rm(list = ls())
#' n <- 500
#' deltaseq <- seq(0.1, 1, 0.1)
#' a <- 0.2
#' b <- -0.1
#' p <- 0.5
#' ### True MCID is 0.5 ###
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
#' delsel <- sel$'Selected delta'
#'
#' result <- pmcid(x = x, y = y, n = n, delta = delsel,
#'             maxit = 100, tol = 1e-02, alpha = 0.05)
#' result$'Point estimate'
#' result$'Standard error'
#' result$'Confidence interval'
#'
#'
pmcid <- function(x, y, n, delta, maxit = 100, tol = 1e-02, alpha = 0.05) {
  w <- ifelse(y == 1, n / sum(y == 1), n / sum(y == -1))
  opt_init <- optim(par = 0, fn = pmcid.hinge.smooth, x = x, y = y, w = w, n = n, delta = delta, method = "L-BFGS-B")
  tau <- opt_init$par
  iter <- 0
  oldvalfun <- 2000
  newvalfun <- 1000
  while ((oldvalfun - newvalfun) >= tol & iter <= maxit) {
    iter<- iter + 1
    oldvalfun <- newvalfun
    u_new <- y * (x - tau)
    t_new <- -1 / n * ((-2 / delta) * ((1 - 2 / delta * u_new) * ifelse(u_new < delta / 2 & u_new >= 0, 1, 0) + ifelse(u_new < 0, 1, 0)))
    opt <- optim(par = 0, fn = pmcid.ramp.smooth, x = x, y = y, w = w, n = n, t = t_new, delta = delta, method = "L-BFGS-B")
    tau <- opt$par
    newvalfun <- opt$value
  }
  h_f <- c()
  g_f <- c()
  for (i in 1:n) {
    val <- 1 / delta * y[i] * (x[i] - tau)
    ind1 <- ifelse(val > 0.5 & val <= 1, 1, 0)
    ind2 <- ifelse(val > 0 & val <= 0.5, 1, 0)
    h_f[i] <- matrix(h.pfun(ind1, ind2, w[i]))
    g_f[i] <- matrix(g.pfun(val, ind1, ind2, w[i]))
  }
  hmean <- mean(h_f)
  gmean <- mean(g_f)
  var_try <- try(delta ^ 2 * solve(hmean) %*% gmean %*% solve(hmean))
  if (!is.matrix(var_try)) {
    hmean <- hmean + 10 ^ (-6)
    var_est <- delta ^ 2 * solve(hmean) %*% gmean %*% solve(hmean)
  } else {
    var_est <- var_try
  }
  pmcid_ci <- c()
  pmcid_ci <- c(tau - qnorm(1 - alpha / 2) * sqrt(var_est / n), tau + qnorm(1 - alpha / 2) * sqrt(var_est / n))
  result <- list("Point estimate" = tau, "Standard error" = drop(sqrt(var_est / n)), "Confidence interval" = pmcid_ci)
  return(result)
}

