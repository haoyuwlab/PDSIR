

#' Probability of being removed before next observation time
#'
#' Computes p_i, the probability of being removed before t_end given an infection at time tau_T.
#'
#' @param theta parameters
#' @param tau_T observed event times
#' @param iota_dist distribution of infection periods
#' @param t_end time before which individuals may be removed
#'
#' @return a vector of the probabilities that individual infected at time tau_T are removed before t_end
#' @export
#'
compute_pi <- function(theta, tau_T, iota_dist, t_end) {

      # Setup
      gamma  <- theta[["gamma" ]]
      lambda <- theta[["lambda"]]
      shape  <- theta[["shape" ]]

      # Compute p_i
      p_i <- if(iota_dist == "exponential") {
            stats::pexp(t_end - tau_T, gamma)
      } else if(iota_dist == "weibull") {
            pweibull2(t_end - tau_T, shape, lambda)
      }

      # Output
      return(p_i)

}


#' Generates infection times
#'
#' @param n number of infection times to generate
#' @param mu individual rate of infection
#' @param lower lower bound
#' @param upper upper bound
#' @param approx poisson approximation
#'
#' @return infection times
#' @export
#'
propose_tau_T <- function(n, mu, lower, upper, approx) {

      stopifnot(approx %in% c("poisson", "ldp"))

      tau_T <- if(approx == "poisson") {
            stats::runif(n,     lower, upper)
      } else if(approx == "ldp") {
            rexp_trunc  (n, mu, lower, upper)
      }

      return(tau_T)

}



#' Generates removal times
#'
#' @inheritParams compute_pi
#'
#' @return removal times
#' @export
#'
propose_tau_J <- function(theta, tau_T, iota_dist, t_end) {

      # Setup
      gamma  <- theta[["gamma" ]]
      lambda <- theta[["lambda"]]
      shape  <- theta[["shape" ]]

      # Compute p_i
      iota_obs <- if(iota_dist == "exponential") {
            rexp_trunc    (length(tau_T), gamma         , 0, t_end - tau_T)
      } else      if(iota_dist == "weibull") {
            rweibull2_trunc(length(tau_T), shape, lambda, 0, t_end - tau_T)
      }

      tau_J <- tau_T + iota_obs

      # Output
      return(tau_J)

}



#' Generates PD-SIR process conditionally on the observed data
#'
#' @param theta current values of the parameters for the SIR
#' @param Y     observed data
#' @param gener     generalized
#' @param b     parameter of generalized
#' @param iota_dist     distribution of infection periods
#' @param approx     poisson approximation
#' @param x     current configuration of the latent data
#' @param rho     proportion
#'
#' @return latent data for the SIR
#' @export
#'
rprop_x <- function(
            theta, Y,
            gener, b, iota_dist, approx,
            x = NULL, rho = 1
) {

      # Setup

      I_k    <- Y[["I_k"  ]]
      I0     <- Y[["I0"   ]]
      S0     <- Y[["S0"   ]]
      ts     <- Y[["ts"   ]]
      t_end  <- Y[["t_end"]]

      theta <- complete_theta(theta, iota_dist, S0)

      beta   <- theta[["beta"  ]]
      gamma  <- theta[["gamma" ]]
      lambda <- theta[["lambda"]]
      shape  <- theta[["shape" ]]

      K        <- length(I_k)
      I_cumsum <- cumsum(I_k)
      I_sum    <- I_cumsum[K] + I0
      I_t_k    <- J_k <- mu_k <- rep(NA, K)
      U_k      <- c(0, I_cumsum)
      S_k      <- (S0 - U_k)[1 : K]
      p_i      <- rep(NA, I_sum)

      if(is.null(x)) { # for initial iteration of Markov chain
            tau_T  <- rep(Inf, I_sum)
            tau_J  <- rep(Inf, I_sum)
            rho    <- 1
      } else {
            tau_T <- x[["tau_T"]]
            tau_J <- x[["tau_J"]]
      }

      # Indices of Updates
      i_update        <- which(as.logical(stats::rbinom(I_sum, 1, rho)))
      tau_J[i_update] <- Inf # tau_J that are updated are set to Inf by default

      # Initially Infectious
      I_t_k[1]   <- I0
      i_k        <- 1 : I0
      i_k_update <- intersect(i_k, i_update)
      n_k_update <- length(i_k_update)

      if(n_k_update > 0) {

            tau_T[i_k_update] <- 0 # should not be 0 for non-Markovian process

            # Compute p_i
            p_i[i_k_update] <- compute_pi(theta, tau_T[i_k_update], iota_dist, t_end)

            # Sample particles recovering before t_end
            is_obs    <- as.logical(stats::rbinom(n_k_update, 1, p_i[i_k_update]))
            tau_J_obs <- i_k_update[is_obs] # indices of particles recovering before t_end

            # Propose tau_J (the tau_J not updated are Inf by default)
            tau_J[tau_J_obs] <- propose_tau_J(theta, tau_T[tau_J_obs], iota_dist, t_end)

      }


      # Update trajectories of particles infected in interval k
      for(k in 1 : K) {

            if(I_k[k] > 0) {

                  # Verify that I(t_k) > 0
                  if(I_t_k[k] == 0)  return(list(compatible = FALSE))

                  i_k        <- I0 + (U_k[k] + 1) : U_k[k + 1] # indices of particles infected during interval k
                  i_k_update <- intersect(i_k, i_update)
                  n_k_update <- length(i_k_update)

                  if(n_k_update > 0) {

                        # Propose tau_T
                        mu_k[k]  <- if(gener)  beta * I_t_k[k] * S_k[k]^(-b)  else  beta * I_t_k[k]

                        tau_T[i_k_update] <- propose_tau_T(n_k_update, mu_k[k], ts[k], ts[k+1], approx)

                        # Propose tau_J
                        p_i[i_k_update]  <- compute_pi(theta, tau_T[i_k_update], iota_dist, t_end)
                        is_obs           <- as.logical(stats::rbinom(n_k_update, 1, p_i[i_k_update])) # particles recovering before t_end
                        tau_J_obs        <- i_k_update[is_obs]
                        tau_J[tau_J_obs] <- propose_tau_J(theta, tau_T[tau_J_obs], iota_dist, t_end = t_end)

                  } # end-if(n_k_update)

            } # end-if(I_k[k])

            # Update I_t_k
            if(k < K) {
                  J_k[k    ] <- sum(dplyr::between(tau_J, ts[k], ts[k + 1]))
                  I_t_k[k + 1] <- I_t_k[k] + I_k[k] - J_k[k]
            }

      } # end-for

      # Output
      x_new <- list(
            compatible = TRUE, tau_T = tau_T, tau_J = tau_J, i_update = i_update, I_t_k = I_t_k, S_k = S_k
      )
      return(x_new)

}



#' Log proposal density of observed removals
#'
#' @param theta parameters of SIR model
#' @param tau_T infection times
#' @param tau_J removal times
#' @param i_update_obs indices of observations to update
#' @param iota_dist distribution of infection periods
#' @param t_end end of observation window
#'
#' @return log density
#' @export
#'
contribution_observed_removal <- function(
            theta, tau_T, tau_J, i_update_obs, iota_dist, t_end
) {

      # Setup
      gamma  <- theta[["gamma" ]]
      lambda <- theta[["lambda"]]
      shape  <- theta[["shape" ]]

      # loglik
      iota_obs <- tau_J[i_update_obs] - tau_T[i_update_obs]
      loglik <- if(iota_dist == "exponential") {
            dexp_trunc_log    (iota_obs, gamma         , 0, t_end - tau_T[i_update_obs])
      } else if(iota_dist == "weibull") {
            dweibull2_trunc_log(iota_obs, shape, lambda, 0, t_end - tau_T[i_update_obs])
      }

      # Output
      return(loglik)

}

#
#' Log proposal density for the latent data
#'
#' @inheritParams rprop_x
#'
#' @param i_update index set of particles whose values are updated
#'
#' @return log density
#' @export
#'
dprop_x <- function(
            theta, Y, x, i_update,
            gener, b, iota_dist, approx
) {

      # Setup

      I_k    <- Y[["I_k"  ]]
      I0     <- Y[["I0"   ]]
      ts     <- Y[["ts"   ]]
      t_end  <- Y[["t_end"]]
      K      <- length(I_k)
      I_sum  <- sum(I_k) + I0

      beta   <- theta[["beta"  ]]
      gamma  <- theta[["gamma" ]]
      lambda <- theta[["lambda"]]
      shape  <- theta[["shape" ]]

      tau_T  <- x    [["tau_T"]]
      tau_J  <- x    [["tau_J"]]
      I_t_k  <- x    [["I_t_k"]]
      S_k    <- x    [["S_k"  ]]


      # Contribution of infections
      contribution_infection <- if(approx == "poisson") {
            0
      } else if(approx == "ldp") { # linear death process

            # exclude initially infectious particles
            i_update_tau_T <- setdiff(i_update, 1 : I0)

            mu_k  <- if(gener)  beta * I_t_k * S_k^(-b)  else  beta * I_t_k
            low   <- ts[1 : K      ] # TODO: compute only once and attach to Y
            upp   <- ts[2 : (K + 1)]

            mu_i  <- rep(mu_k, I_k)
            low_i <- rep(low , I_k)
            upp_i <- rep(upp , I_k)

            sum(dexp_trunc_log(
                  tau_T[i_update_tau_T     ], mu_i [i_update_tau_T - I0],
                  low_i[i_update_tau_T - I0], upp_i[i_update_tau_T - I0]
            ))
      } # end-if(approximation)

      # Contribution of removal
      p_i <- compute_pi(theta, tau_T, iota_dist, t_end)

      i_update_obs     <- setdiff(i_update, which(is.infinite(tau_J)))
      i_update_not_obs <- setdiff(i_update, which(is.finite  (tau_J)))

      contribution_removal <-
            sum(
                  log(1 - p_i[i_update_not_obs])
            ) +
            sum(
                  log(p_i[i_update_obs]) +
                        contribution_observed_removal(
                              theta, tau_T, tau_J, i_update_obs, iota_dist, t_end
                        )
            )

      # Output
      loglik <- contribution_infection + contribution_removal
      return(loglik)

}
