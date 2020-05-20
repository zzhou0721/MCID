#' The value function used to search the initial value for determining the individualized MCID
#' @keywords internal
imcid.hinge.smooth <- function(x, y, z, w, n, parm, lambda, delta) {
  beta0 <- parm[1]
  beta <- matrix(parm[-1])
  u <- y * (x - beta0 - z %*% beta)
  val1 <- w * 2 * (1 - 1 / delta * u) ^ 2 * ifelse(u < delta & u >= delta / 2, 1, 0)
  val2 <- w * (1.5 - 2 / delta * u) * ifelse(u < delta / 2, 1, 0)
  val <- n * lambda / 2 * sum(beta ^ 2) + sum(val1 + val2)
  return(val)
}

#' The value function needed to be optimized for determining the individualized MCID
#' @keywords internal
imcid.ramp.smooth <- function(x, y, z, w, n, parm, lambda, delta, t) {
  beta0 <- parm[1]
  beta <- matrix(parm[-1])
  u <- y * (x - beta0 - z %*% beta)
  val1 <- w * 2 * (1 - 1 / delta * u) ^ 2 * ifelse(u < delta & u >= delta / 2, 1, 0)
  val2 <- w * (1.5 - 2 / delta * u) * ifelse(u < delta / 2, 1, 0)
  val3 <- w * u * t
  val <- n * lambda / 2 * sum(beta ^ 2) + sum(val1 + val2) + n * sum(val3)
  return(val)
}


#' The hassen matrix function for determining the individualized MCID
#' @keywords internal
h.ifun <- function(z_tilde, ind1, ind2, w) {
  h_fun <- drop(4 * (ind1 - ind2)) * w * (z_tilde %*% t(z_tilde))
  return(h_fun)
}

#' The square of the score function for determining the individualized MCID
#' @keywords internal
g.ifun <- function(z_tilde, val, ind1, ind2, w) {
  g_fun <- drop(4 * ((2 - 2 * val) ^ 2 * ind1 + (2 * val) ^ 2 * ind2)) * w ^ 2 * (z_tilde %*% t(z_tilde))
  return(g_fun)
}
