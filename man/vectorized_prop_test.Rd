% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vectorized-prop-test.R
\name{vectorized_prop_test}
\alias{vectorized_prop_test}
\title{Perform vectorized hypothesis tests for the difference of two proportions}
\usage{
vectorized_prop_test(x1, n1, x2, n2, conf.level = 0.95)
}
\arguments{
\item{a}{Number of successes in first condition}

\item{b}{Number of failures in first condition}

\item{c}{Number of successes in second condition}

\item{d}{Number of failures in second condition}
}
\description{
Given vectors a, b, c and d, return a vector of
p-values that test whether a / (a + b) is the same proportion
as c / (c + d). This is used in the simulations of p-value
behavior.
}
