
#' Log likelihood of the stochastic SIR model
#'
#'
#' allows non-Markovian dynamics through Weibull-distribution infectious periods (Streftaris and Gibson, 2004)
#' allows generalized SIR process (Severo, 1969)
#' Severo, N. C. (1969). Generalizations of some stochastic epidemic models. Mathematical Biosciences, 4(3-4), 395-402.
#' Streftaris, G., & Gibson, G. J. (2004). Bayesian inference for stochastic epidemics in closed populations. Statistical Modelling, 4(1), 63-75.
#'
#'
#' @inheritParams rprop_x
#'
#' @param SS sufficient statistics of the current configuration of the latent data
#'
#' @return log likelihood
#' @export
#'
f_log <- function(theta, SS, gener, b, iota_dist = "exponential") {

      # Setup
      beta   <- theta[["beta"  ]]
      gamma  <- theta[["gamma" ]]
      shape  <- theta[["shape" ]]
      lambda <- theta[["lambda"]]

      n_I             <- SS[["n_I"            ]]
      I_tau_I         <- SS[["I_tau_I"        ]]
      S_tau_I         <- SS[["S_tau_I"        ]]
      integral_SI     <- SS[["integral_SI"    ]]
      iota_removed    <- SS[["iota_removed"   ]]
      iota_infectious <- SS[["iota_infectious"]]

      # # Loglik
      # loglik <- # agent-based likelihood
      #   n_I   * log(beta) + sum(log(I_tau_I)) +
      #   n_R   * log(gamma) -
      #   beta  * integral_SI -
      #   gamma * integral_I

      # contribution of infections
      loglik_infec <- if(gener){
            n_I * log(beta) + sum(log(I_tau_I) - b * log(S_tau_I)) - beta  * integral_SI
      } else {
            n_I * log(beta) + sum(log(I_tau_I)                   ) - beta  * integral_SI
      }

      # contribution of removals
      loglik_remov <-
            if(iota_dist == "exponential") {
                  sum(stats::dexp(iota_removed   , gamma        , log   = TRUE                    )) + # removals before t_end (observed)
                        sum(stats::pexp(iota_infectious, gamma        , log.p = TRUE, lower.tail = FALSE))   # removals after  t_end (not observed)
            } else if(iota_dist == "weibull") {
                  sum(dweibull2  (iota_removed   , shape, lambda, log   = TRUE                    )) +
                        sum(pweibull2  (iota_infectious, shape, lambda, log.p = TRUE, lower.tail = FALSE))
            }

      loglik <- loglik_infec + loglik_remov

      return(loglik)

}
