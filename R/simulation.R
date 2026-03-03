#' Simulation of a Poisson process
#' 
#' Simulation of a (possibly inhomogeneous) Poisson process using the thinning 
#' algorithm in  Laub, Taimre, and Pollett (2021) Chapter 4 (Algorithm 1). It 
#' works by generating a homogeneous Poisson process with rate `M` on 
#' \eqn{[0,T]}, and then, in the inhomogeneous case, thinning this process using
#' the provided intensity function `fn`
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param fun A non-negative function (or numeric if simulating a homogeneous
#'   Poisson process) representing the intensity function, \eqn{\lambda (t)} of 
#'   the Poisson process.
#' @param M A non-negative numeric value - upper bound on `fun`. Ignored for 
#'   a homogeneous Poisson process.
#'
#' @return A vector of simulated arrival times from the process.
#'
#' @export
#'
#' @examples
#' # simulate an inhomogeneous Poisson process with Exp(2) intensity
#' t <- simulate_pp(T = 5, fn = function(y\t) {2*exp(-2*t)}, M = 2)
#' # Simulate a homogeneous Poisson process with λ = 2
#' x <- simulate_pp(T = 5, fn = 2)
sim_pp <- function(T, fn, M = NULL) {
  
  if (is.numeric(fn)) {
    # number of events in interval of length T is Poisson(λT) distributed
    N_T <- rpois(1, fn*T)
    # given N(T) = n, the times of arrivals tn,...,t_n are U(0, T) distributed
    arrival_times <- sort(runif(N_T, 0, T))
    return(arrival_times)
  }
  
  # some guardrails
  if (is.null(M)) {
    stop("For inhomogeneous case, 'sim_pp' requires an upper bound of the function over [0, T].")
  }
  
  if (!is.function(fn)) {
    stop("'fn' needs to be a non-negative function or numeric.")
  }
  
  # simulate homogeneous poisson process with rate M on [0, T]
  N_T <- rpois(1, M*T)
  unthinned_arrival_times <- sort(runif(N_T, 0, T))
  # keep arrivals with some probability
  u <- runif(N_T, 0, M)
  keep <- ifelse(u <= fn(unthinned_arrival_times), TRUE, FALSE)
  
  sort(unthinned_arrival_times[keep])
  
}


#' Simulation of a Hawkes process by thinning
#' 
#' Simulation of a Hawkes process using Ogata's modified thinning algorithm, 
#' following Chapter 4 (Algorithm 2) of Laub, Taimre, and Pollett (2021). This 
#' is similar to the thinning method for an inhomogeneous Poisson process, but
#' in this case we have no a.s. asymptotic bound M for the conditional intensity 
#' \eqn{\lambda^*(t)}. Instead, we restrict the function space to be only 
#' non-increasing functions of $t$, which allows us to simply use the value of
#' `fn(t_i)` where $t_i$ is the time just after an arrival.
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param lambda A non-negative numeric value - the background arrival 
#'   component \eqn{\lambda} of the Hawkes conditional intensity
#' @param mu A non-negative function - the excitation function/kernel 
#'   \eqn{\mu(\cdot)} of the Hawkes conditional intensity
#'
#' @return A vector of simulated arrival times from the process.
#'
#' @export
#'
#' @examples
#' # simulate a Hawkes process with background rate 1 and Exp(2) kernel
#' t <- simulate_hp(T=5, lambda = 1, mu = function(t) {2*exp(-2*t)})
sim_hp <- function(T, lambda, mu) {
  
  # some guardrails
  if (!is.function(mu)) {
    stop("'mu' needs to be a non-negative function.")
  }
  
  if (!is.numeric(lambda) & lambda < 0) {
    stop("'lambda' needs to be a non-negative numeric.")
  }
  
  # initialise simulation objects
  arrival_times <- c()
  t <- 0
  while (t < T) {
    # find new upper bound
    M <- lambda + sum(mu(t - arrival_times))
    # generate next candidate point
    t <- t + rexp(1, M)
    # keep it with some probability
    u <- runif(1, 0, M)
    
    if (t < T & u <= (lambda + sum(mu(t - arrival_times)))) {
      arrival_times <- c(arrival_times, t)
    }
  }
  
  arrival_times
  
}