# Introducing two helper functions to remove dependency from extraDistr (orphaned since 01-2026)

#' Log PMF of multivariate hypergeometric
#'
#' @param x integer vector of counts drawn in each category
#' @param n integer vector of population sizes in each category
#' @param k integer total draws (should equal sum(x))
#' @param log logical; return log-probability?
#' @keywords internal
dmvhyper_base <- function(x, n, k = sum(x), log = FALSE) {
  x <- as.integer(x)
  n <- as.integer(n)
  k <- as.integer(k)
  
  if (length(x) != length(n)) stop("x and n must have the same length.")
  if (anyNA(x) || anyNA(n) || is.na(k)) return(if (log) NA_real_ else NA_real_)
  if (any(n < 0L)) stop("n must be non-negative.")
  if (k < 0L) return(if (log) -Inf else 0)
  
  # invalid support
  if (any(x < 0L) || any(x > n) || sum(x) != k) {
    return(if (log) -Inf else 0)
  }
  
  N <- sum(n)
  if (k > N) return(if (log) -Inf else 0)
  
  lp <- sum(lchoose(n, x)) - lchoose(N, k)
  if (log) lp else exp(lp)
}


#' Random samples from multivariate hypergeometric
#'
#' @param nn number of samples
#' @param n integer vector of population sizes
#' @param k integer total draws
#' @return integer matrix with nn rows and length(n) columns
#' @keywords internal
rmvhyper_base <- function(nn, n, k) {
  n <- as.integer(n)
  k <- as.integer(k)
  nn <- as.integer(nn)
  
  if (nn < 0L) stop("nn must be non-negative.")
  if (any(n < 0L)) stop("n must be non-negative.")
  N <- sum(n)
  if (k < 0L || k > N) stop("k must be between 0 and sum(n).")
  
  p <- length(n)
  out <- matrix(0L, nrow = nn, ncol = p)
  
  if (nn == 0L) return(out)
  if (p == 0L) return(out)
  if (k == 0L) return(out)
  
  k_rem <- rep.int(k, nn)
  N_rem <- N
  
  # sequential draws for categories 1..(p-1)
  if (p > 1L) {
    for (i in seq_len(p - 1L)) {
      ni <- n[i]
      # rhyper(nn, m=ni successes, n=N_rem-ni failures, k=k_rem draws)
      xi <- stats::rhyper(nn, m = ni, n = N_rem - ni, k = k_rem)
      out[, i] <- xi
      k_rem <- k_rem - xi
      N_rem <- N_rem - ni
    }
  }
  
  # last category gets remainder
  out[, p] <- k_rem
  out
}