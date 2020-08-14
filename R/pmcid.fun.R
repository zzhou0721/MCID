#' The value function used to search the initial value for determining the population MCID
#' @keywords internal
pmcid.hinge.smooth <- function(x, y, w, n, tau, delta) {
  u <- y * (x - tau)
  val1 <- w * 2 * (1 - 1 / delta * u) ^ 2 * ifelse(u < delta & u >= delta / 2, 1, 0)
  val2 <- w * (1.5 - 2 / delta * u) * ifelse(u < delta / 2, 1, 0)
  val <- sum(val1 + val2)
  return(val)
}


#' The value function needed to be optimized for determining the population MCID
#' @keywords internal
pmcid.ramp.smooth <- function(x, y, w, n, tau, t, delta) {
  u <- y * (x - tau)
  val1 <- w * 2 * (1 - u / delta) ^ 2 * ifelse(u < delta & u >= delta / 2, 1, 0)
  val2 <- w * (1.5 - 2 * u / delta) * ifelse(u < delta / 2, 1, 0)
  val3 <- w * u * t
  val <- sum(val1 + val2) + n * sum(val3)
  return(val)
}


#' The hassen matrix function for determining the population MCID
#' @keywords internal
h.pfun <- function(ind1, ind2, w) {
  h_fun <- drop(4 * (ind1 - ind2) * w)
  return(h_fun)
}


#' The square of the score function for determining the population MCID
#' @keywords internal
g.pfun <- function(val, ind1, ind2, w) {
  g_fun <- drop(4 * ((2 - 2 * val) ^ 2 * ind1 + (2 * val) ^ 2 * ind2) * w ^ 2)
  return(g_fun)
}
