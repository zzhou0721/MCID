#' Point and interval estimation for the MCID at the individual level
#'
#' We formulate the individualized MCID as a linear function of the patients' clinical profiles. \code{imcid} returns the point estimate for the linear coefficients of the MCID at the individual level
#'
#' @param x a continuous variable denoting the outcome change of interest
#' @param y a binary variable indicating the patient-reported outcome derived from the anchor question
#' @param z a vector or matrix denoting the patient's clinical profiles
#' @param n the sample size
#' @param lambda the selected tuning parameter \eqn{\lambda}, can be returned by \code{cv.imcid}
#' @param delta the selected tuning parameter \eqn{\delta}, can be returned by \code{cv.imcid}
#' @param maxit the maximum number of iterations. Defaults to 100
#' @param tol the convergence tolerance. Defaults to 0.01
#' @param alpha nominal level of the confidence interval. Defaults to 0.05
#'
#' @return a list including the point estimates for the linear coefficients of the individualized MCID and their standard errors, and the corresponding confidence intervals based on the asymptotic normality
#'
#' @importFrom stats optim qnorm
#'
#' @export
#'
#' @examples
#' \donttest{
#' rm(list = ls())
#' n <- 500
#' lambdaseq <- 10 ^ seq(-3, 3, 0.1)
#' deltaseq <- seq(0.1, 0.3, 0.1)
#' a <- 0.1
#' b <- 0.55
#' c <- -0.1
#' d <- 0.45
#' ### True linear coefficients of the individualized MCID: ###
#' ### beta0=0, beta1=0.5 ###
#'
#' set.seed(115)
#' p <- 0.5
#' y <- 2 * rbinom(n, 1, p) - 1
#' z <- rnorm(n, 1, 0.1)
#' y_1 <- which(y == 1)
#' y_0 <- which(y == -1)
#' x <- c()
#' x[y_1] <- a + z[y_1] * b + rnorm(length(y_1), 0, 0.1)
#' x[y_0] <- c + z[y_0] * d + rnorm(length(y_0), 0, 0.1)
#' sel <- cv.imcid(x = x, y = y, z = z, lamseq = lambdaseq,
#'          delseq = deltaseq, k = 5, maxit = 100, tol = 1e-02)
#' lamsel <- sel$'Selected lambda'
#' delsel <- sel$'Selected delta'
#' result <- imcid(x = x, y = y, z = z, n = n, lambda = lamsel,
#'          delta = delsel, maxit = 100, tol = 1e-02, alpha = 0.05)
#' result$'Point estimates'
#' result$'Standard errors'
#' result$'Confidence intervals'
#' }
imcid <- function(x, y, z, n, lambda, delta, maxit = 100, tol = 1e-02, alpha = 0.05) {
  w <- ifelse(y == 1, n / sum(y == 1), n / sum(y == -1))
  opt_init <- optim(par = rep(0, ncol(cbind(x, y, z)) - 1), fn = imcid.hinge.smooth, x = x, y = y, z = z, w = w, n = n, lambda = lambda, delta = delta, method = "L-BFGS-B")
  par_new <- opt_init$par
  iter <- 0
  oldvalfun <- 2000
  newvalfun <- 1000
  while ((oldvalfun - newvalfun) >= tol & iter <= maxit) {
    iter <- iter + 1
    oldvalfun <- newvalfun
    u_new <- y * (x - par_new[1] - z %*% matrix(par_new[-1]))
    t_new <- -1 / n * ((-2 / delta) * ((1 - 2 / delta * u_new) * ifelse(u_new < delta / 2 & u_new >= 0, 1, 0) + ifelse(u_new < 0, 1, 0)))
    opt <- optim(par = rep(0, ncol(cbind(x, y, z)) - 1), fn = imcid.ramp.smooth, x = x, y = y, z = z, w = w, n = n, lambda = lambda, delta = delta, t = t_new, method = "L-BFGS-B")
    par_new <- opt$par
    newvalfun <- opt$value
  }
  h_f <- list()
  g_f <- list()
  for (i in 1:n) {
    if (is.null(dim(z))) {
      z_tilde <- c(1, z[i])
    } else {
      z_tilde <- c(1, z[i, ])
    }
    val <- 1 / delta * y[i] * (x[i] - z_tilde %*% par_new)
    ind1 <- ifelse(val > 0.5 & val <= 1, 1, 0)
    ind2 <- ifelse(val > 0 & val <= 0.5, 1, 0)
    h_f[[i]] <- h.ifun(z_tilde, ind1, ind2, w[i])
    g_f[[i]] <- g.ifun(z_tilde, val, ind1, ind2, w[i])
  }
  hmean <- apply(simplify2array(h_f), 1:2, mean)
  gmean <- apply(simplify2array(g_f), 1:2, mean)
  var_try <- try(delta ^ 2 * solve(hmean) %*% gmean %*% solve(hmean))
  if (!is.matrix(var_try)) {
    hmean <- hmean + 10 ^ (-6)
    var_est <- delta ^ 2 * solve(hmean) %*% gmean %*% solve(hmean)
  } else {
    var_est <- var_try
  }
  coef_name <- paste0("beta",0:(ncol(cbind(x, y, z))-2))
  coef_est <- matrix(par_new)
  row.names(coef_est) <- coef_name
  coef_se <- matrix(sqrt(diag(var_est) / n))
  row.names(coef_se) <- coef_name
  imcid_ci <- matrix(NA, nrow = nrow(var_est), ncol = 2)
  imcid_ci <- t(apply(cbind(par_new, diag(var_est)), 1, function(x){c(x[1] - qnorm(1 - alpha / 2) * sqrt(x[2] / n), x[1] + qnorm(1 - alpha / 2) * sqrt(x[2] / n))}))
  row.names(imcid_ci) <- coef_name
  colnames(imcid_ci) <- c("Lower bound", "Upper bound")
  result <- list("Point estimates" = coef_est, "Standard errors" = coef_se, "Confidence intervals" = imcid_ci)
  return(result)
}
