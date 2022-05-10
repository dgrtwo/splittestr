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
vectorized_prop_test <- function(x1, n1, x2, n2, conf.level = .95) {
  a <- x1
  b <- n1 - x1
  c <- x2
  d <- n2 - x2

  # if any values are < 20, use Fisher's exact test
  exact <- (a < 20 | b < 20 | c < 20 | d < 20)

  pvalue <- rep(NA, length(a))

  if (any(exact)) {
    pvalue[exact] <- vectorized_prop_test_exact(a[exact], b[exact], c[exact], d[exact])
  }
  if (any(!exact)) {
    pvalue[!exact] <- vectorized_prop_test_approx(a[!exact], b[!exact], c[!exact], d[!exact])
  }

  mu1 <- a / (a + b)
  mu2 <- c / (c + d)

  ## confidence interval
  alpha2 <- (1 - conf.level) / 2
  DELTA <- mu2 - mu1
  WIDTH <- qnorm(alpha2)
  alpha <- (a + .5) / (a + b + 1)
  beta <- (c + .5) / (c + d + 1)

  n <- n1 + n2
  YATES <- pmin(.5, abs(DELTA) / sum(1 / n1 + 1 / n2))

  # TODO: add yates correction
  z <- qnorm((1 + conf.level) / 2)
  WIDTH <- z * sqrt(mu1 * (1 - mu1) / n1 + mu2 * (1 - mu2) / n2)

  tibble::tibble(estimate = DELTA,
                 conf.low = pmax(DELTA - WIDTH, -1),
                 conf.high = pmin(DELTA + WIDTH, 1),
                 p.value = pvalue)
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
