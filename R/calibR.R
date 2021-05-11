#' @title  Calibration by using raking ratio
#'
#' @description 
#' 
#' This function is inspired by the function \code{\link{}}
#'
#' @param Xs 
#' @param d 
#' @param total 
#' @param q 
#' @param max_iter 
#' 
#' @importFrom MASS ginv
#' 
#' @return Return the g-weights.
#'
#' @examples
calibR <- function (Xs, d, total, q = rep(1, length(d)), max_iter = 500) {
  EPS = .Machine$double.eps
  EPS1 = 1e-09
  n = length(d)
  lambda = as.matrix(rep(0, n))
  lambda1 = MASS::ginv(t(Xs * d * q) %*% Xs, tol = EPS) %*% (total - as.vector(t(d) %*% Xs))
  lambda = as.matrix(rep(0, ncol(Xs)))
  w1 = as.vector(d * exp(Xs %*% lambda * q))
  for (l in 1:max_iter) {
    phi = t(Xs) %*% w1 - total
    T1 = t(Xs * w1)
    phiprim = T1 %*% Xs
    lambda = lambda - MASS::ginv(phiprim, tol = EPS) %*% phi
    w1 = as.vector(d * exp(Xs %*% lambda * q))
    if (any(is.na(w1)) | any(is.infinite(w1))) {
      warning("No convergence")
      g = NULL
      break
    }
    tr = crossprod(Xs, w1)
    expression = max(abs(tr - total)/total)
    if (any(total == 0)) {
      expression = max(abs(tr - total))
    }
    if (expression < EPS1) {
      break
    }
  }
  if (l == max_iter) {
    warning("No convergence")
    g = NULL
  }else{
    g = w1/d
  } 
}
