#' Find the expected loss in a conversion experiment
#'
#' This is a wrapper around \code{\link{expected_loss}}. It applies
#' a prior and calculates in terms of successes and totals of A and B,
#' instead of the beta distribution. Note that the "error" is treatment
#' A being better than B.
#'
#' @param success_A number of successes (e.g. clicks) in condition A
#' @param total_A total number of exposures in A
#' @param success_B number of successes (e.g. clicks) in condition B
#' @param total_B total number of exposures in B
#' @param prior_alpha prior parameter to use for alpha (default uniform)
#' @param prior_beta prior parameter to use for beta (default uniform)
#' @param approx whether to use a normal approximation for h
#' @param ... Extra arguments, passed on to \code{\link{expected_loss}}
#'
#' @export
conversion_expected_loss <- function(sA, nA,
                                     sB, nB,
                                     prior_alpha = 1, prior_beta = 1,
                                     approx = FALSE, ...) {
  if (!approx[1] & length(sA) > 1) {
    # cannot vectorize; must sapply
    a <- sA + prior_alpha
    b <- nA - sA + prior_beta
    c <- sB + prior_alpha
    d <- nB - sB + prior_beta
    ret <- sapply(seq_along(a), function(i) {
      expected_loss(a[i], b[i], c[i], d[i], ...)
    })
    return(ret)
  }

  expected_loss(sA + prior_alpha, nA - sA + prior_beta,
                sB + prior_alpha, nB - sB + prior_beta,
                approx = approx, ...)
}
