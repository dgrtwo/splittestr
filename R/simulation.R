#' Perform a simulation of A/B testing based on a data frame of parameters
#'
#' @param params A data frame of parameters, each of which will be passed on to
#' \code{\link{conversion_expected_loss}}
#' @param days Number of days to run
#' @param per_day Number of "clicks" per day; overridden if present in params
#' @param proportion_A "Clickthrough rate" in treatment A
#' @param ... Extra arguments, passed on to \code{\link{conversion_expected_loss}}
#'
#' @import dplyr
#'
#' @examples
#'
#' library(dplyr)
#'
#' sim <- tibble(replicate = seq_len(nreps)) %>%
#'   mutate(proportion_A = .001, effect = 0, per_day = 10000) %>%
#'   perform_simulation()
#'
#' library(ggplot2)
#' ggplot(sim, aes(day, expected_loss, group = replicate)) +
#'   geom_line(alpha = .5) +
#'   scale_y_log10()
#'
#' @export
perform_simulation <- function(params, days = 20, per_day = 100,
                               proportion_A = .1, ...) {
  if (is.null(params$per_day)) {
    params$per_day <- per_day
  }
  if (is.null(params$proportion_A)) {
    params$proportion_A <- proportion_A
  }

  params <- cbind(params, ...)

  # simulate an A/B test for each day
  ret <- params %>%
    tidyr::crossing(day = seq_len(days)) %>%
    tidyr::crossing(type = c("A", "B")) %>%   # two types, A and B
    mutate(type = factor(type)) %>%
    mutate(total = rbinom(n(), .$per_day, .5)) %>%  # half A, half B
    mutate(success = rbinom(n(), total, .$proportion_A + effect * (type == "B"))) %>%            # this is the simulation
    group_by(effect, replicate, type, per_day) %>%
    mutate(n = cumsum(total),                      # count up cumulative
           s = cumsum(success)) %>%
    ungroup()

  # spread into columns for sA, nA, sB, nB
  ret %>%
    select(-success, -total) %>%
    tidyr::gather(metric, value, n:s) %>%
    tidyr::unite(metric2, metric, type, sep = "") %>%
    tidyr::spread(metric2, value) %>%
    mutate(expected_loss = do.call(conversion_expected_loss, .))
}
