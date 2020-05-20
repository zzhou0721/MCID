#' The function used to determine the value of target function for different value of delta in the individualized MCID setting
#' @keywords internal
imcidcv.smooth <- function(x, y, z, lambda, delta, fold, maxit, tol) {
  if (fold == 1) {stop("Fold must be >1")}
  if (maxit < 1) {stop("Maxit must be >=1")}
  if (tol < 0) {stop("Tol must be >0")}
  funval <- c()
  #data <- cbind(x, y, z)[sample(nrow(cbind(x, y, z))), ]
  data <- cbind(x, y, z)
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = FALSE)
  for (jj in 1:fold) {
    testid <- which(folds == jj, arr.ind = TRUE)
    datat <- data[testid, ]
    datatr <- data[-testid, ]
    ntr <- nrow(datatr)
    xtr <- datatr[ ,1]
    ytr <- datatr[ ,2]
    ztr <- datatr[ ,-(1:2)]
    wtr <- ifelse(ytr == 1, ntr / sum(ytr == 1), ntr / sum(ytr == -1))
    xtt <- datat[ ,1]
    ytt <- datat[ ,2]
    ztt <- datat[ ,-(1:2)]
    wtt <- ifelse(ytt == 1, length(ytt) / sum(ytt==1), length(ytt) / sum(ytt == -1))
    opt_init <- optim(par = rep(0, ncol(data) - 1), fn = imcid.hinge.smooth, x = xtr, y = ytr, z = ztr, w = wtr, n = ntr, lambda = lambda, delta = delta, method = "L-BFGS-B")
    par_new <- opt_init$par
    iter <- 0
    oldvalfun <- 2000
    newvalfun <- 1000
    while ((oldvalfun - newvalfun) >= tol & iter <= maxit) {
      iter <- iter + 1
      oldvalfun <- newvalfun
      u_new <- ytr * (xtr - par_new[1] - ztr %*% matrix(par_new[-1]))
      t_new <- -1 / ntr * ((-2 / delta) * ((1 - 2 / delta * u_new) * ifelse(u_new < delta / 2 & u_new >= 0, 1, 0) + ifelse(u_new < 0, 1, 0)))
      opt <- optim(par = rep(0, ncol(data) - 1), fn = imcid.ramp.smooth, x = xtr, y = ytr, z = ztr, w = wtr, n = ntr, lambda = lambda, delta = delta, t = t_new, method = "L-BFGS-B")
      par_new <- opt$par
      newvalfun <- opt$value
    }
    mcid <- par_new[1] + ztt %*% matrix(par_new[-1])
    funval[jj] <- 0.5 * mean(wtt * (1 - ytt * ifelse(xtt - mcid >= 0, 1, -1)))
  }
  return(mean(funval))
}
