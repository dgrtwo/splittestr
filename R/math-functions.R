# here go mathematical functions (as opposed to programming/utility ones)

#' Calculate the probability one random beta variable is greater than another
#'
#' Find the probability that \code{Beta(a, b) > Beta(c, d)}. This comes from
#' https://www.chrisstucchio.com/blog/2014/bayesian_ab_decision_rule.html. Note
#' that the \code{approx = TRUE} version is vectorized, but the exact version
#' is not.
#'
#' @param a alpha parameter for first Beta
#' @param b beta parameter for second Beta
#' @param c alpha parameter for first Beta
#' @param d beta parameter for second Beta
#' @param approx whether to use a normal approximation to the beta
#' @param log_h whether to return \code{log(h(a, b, c, d))} rather than
#' \code{h(a, b, c, d)}
#'
#' @export
h <- function(a, b, c, d, approx = FALSE, log_h = FALSE) {
  if (approx[1]) {
    # use normal approximation to the beta
    u1 <- a / (a + b)
    u2 <- c / (c + d)
    var1 <- a * b / ((a + b) ^ 2 * (a + b + 1))
    var2 <- c * d / ((c + d) ^ 2 * (c + d + 1))
    return(pnorm(0, u2 - u1, sqrt(var1 + var2), log.p = log_h))
  }

  j <- seq(0, c - 1)
  log_vals <- lbeta(a + j, b + d) - log(d + j) - lbeta(1 + j, d) - lbeta(a, b)

  # due to floating point error it is possible to be *very* slightly below 0
  ret <- max(1 - sum(exp(log_vals)), 0)
  if (log_h) {
    return(log(ret))
  }
  ret
}


#' Calculate the expected loss from assuming the wrong beta
#'
#' Find the expected loss from two Beta posteriors if one assumes (wrongly)
#' that a / b > c / d. That is, \code{E[max(Beta(c, d) - Beta(a, b), 0)]}. This comes from
#' https://www.chrisstucchio.com/blog/2014/bayesian_ab_decision_rule.html. Note
#' that the \code{approx = TRUE} version is vectorized, but the exact
#' version is not.
#'
#' @param a alpha parameter for first Beta
#' @param b beta parameter for second Beta
#' @param c alpha parameter for first Beta
#' @param d beta parameter for second Beta
#' @param approx whether to use a normal approximation when calculating h
#' @param ... Extra arguments, ignored
#'
#' @export
expected_loss <- function(a, b, c, d, approx = FALSE, ...) {
  v1 <- lbeta(a + 1, b) - lbeta(a, b) + h(a + 1, b, c, d, approx = approx, log_h = TRUE)
  v2 <- lbeta(c + 1, d) - lbeta(c, d) + h(a, b, c + 1, d, approx = approx, log_h = TRUE)
  exp(v1) - exp(v2)
}
