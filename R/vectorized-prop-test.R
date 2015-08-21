#' Perform vectorized hypothesis tests for the difference of two proportions
#'
#' Given vectors a, b, c and d, return a vector of
#' p-values that test whether a / (a + b) is the same proportion
#' as c / (c + d). This is used in the simulations of p-value
#' behavior.
#'
#' @param a Number of successes in first condition
#' @param b Number of failures in first condition
#' @param c Number of successes in second condition
#' @param d Number of failures in second condition
#'
#' @export
vectorized_prop_test <- function(a, b, c, d) {
  # if any values are < 30, use exact
  exact <- (a < 30 | b < 30 | c < 30 | d < 30)

  result <- rep(NA, length(a))
  result[exact] <- vectorized_prop_test_exact(a[exact], b[exact], c[exact], d[exact])
  result[!exact] <- vectorized_prop_test_approx(a[!exact], b[!exact], c[!exact], d[!exact])

  result
}


#' A vectorized version of the prop.test function
vectorized_prop_test_approx <- function(a, b, c, d) {
  n1 <- a + b
  n2 <- c + d
  n <- n1 + n2
  p <- (a + c) / n
  E <- cbind(p * n1, (1 - p) * n1, p * n2, (1 - p) * n2)

  x <- cbind(a, b, c, d)

  DELTA <- a / n1 - c / n2
  YATES <- pmin(.5, abs(DELTA) / sum(1 / n1 + 1 / n2))

  STATISTIC <- rowSums((abs(x - E) - YATES)^2 / E)
  PVAL <- pchisq(STATISTIC, 1, lower.tail = FALSE)
  PVAL
}


#' Perform a Fisher exact test for each set of a, b, c, d
vectorized_prop_test_exact <- function(a, b, c, d) {
  sapply(seq_along(a), function(i) {
    fisher.test(cbind(c(a[i], c[i]), c(b[i], d[i])))$p.value
  })
}
